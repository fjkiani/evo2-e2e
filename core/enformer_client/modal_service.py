"""
modal_enformer.py — CrisPRO.ai Enformer Service (v2)
Fixes: torch>=2.6 required by transformers CVE-2025-32434
Deploy: modal deploy modal_enformer.py
"""
import modal, os

model_cache = modal.Volume.from_name("crispro-model-cache", create_if_missing=True)

enformer_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install([
        # torch 2.6+ required — CVE-2025-32434 fix
        "torch==2.6.0",
        "enformer-pytorch>=0.8.10",
        "einops", "pandas", "numpy", "requests",
        "scikit-learn", "huggingface_hub>=0.23",
    ])
    .env({"HF_HOME": "/model-cache"})
)

app = modal.App("crispro-enformer")

@app.cls(
    image=enformer_image,
    gpu="T4",
    timeout=120,
    volumes={"/model-cache": model_cache},
    min_containers=0,
    scaledown_window=300,
)
class EnformerService:
    SEQ_LEN      = 196_608
    DNASE_TRACKS = list(range(121, 245))
    ATAC_TRACKS  = list(range(684, 724))

    @modal.enter()
    def load(self):
        import torch
        from enformer_pytorch import Enformer
        print("[enformer] Loading...", flush=True)
        self.model = Enformer.from_pretrained(
            "EleutherAI/enformer-official-rough", target_length=-1).eval()
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.model = self.model.to(self.device)
        print(f"[enformer] Ready on {self.device}", flush=True)

    def _one_hot(self, seq):
        import torch
        m = {"A":0,"C":1,"G":2,"T":3}
        oh = torch.zeros(len(seq), 4)
        for i, c in enumerate(seq.upper()):
            idx = m.get(c)
            if idx is not None: oh[i,idx] = 1.0
        return oh

    def _pad_trim(self, seq):
        if len(seq) < self.SEQ_LEN: return seq + "N"*(self.SEQ_LEN-len(seq))
        s = (len(seq)-self.SEQ_LEN)//2
        return seq[s:s+self.SEQ_LEN]

    def _predict(self, seq):
        import torch
        seq = self._pad_trim(seq)
        enc = self._one_hot(seq).unsqueeze(0).to(self.device)
        with torch.no_grad(): preds = self.model(enc)
        return preds["human"][0].cpu().numpy()

    @modal.method()
    def score_accessibility(self, sequence: str, track_type: str = "dnase") -> dict:
        human  = self._predict(sequence)
        tracks = self.DNASE_TRACKS if track_type == "dnase" else self.ATAC_TRACKS
        score  = float(human[:, tracks].mean())
        return {"accessibility_score": score, "track_type": track_type,
                "n_tracks": len(tracks), "sequence_length": len(sequence)}

    @modal.method()
    def score_vs_target(self, sequence: str, target_peaks: list = None) -> dict:
        import numpy as np
        from sklearn.metrics import roc_auc_score
        human   = self._predict(sequence)
        preds   = human[:, self.DNASE_TRACKS].mean(axis=1)
        n_bins  = len(preds)
        target  = np.zeros(n_bins)
        for p in (target_peaks or []):
            pb = p["position_bp"]//128
            wb = max(1, p["width_bp"]//128)
            target[max(0,pb-wb//2):min(n_bins,pb+wb//2)] = 1.0
        if target.sum() in (0,n_bins): return {"auroc":0.5,"accessibility_score":float(preds.mean())}
        try: auroc = float(roc_auc_score(target, preds))
        except: auroc = 0.5
        return {"auroc":auroc,"accessibility_score":float(preds.mean()),"n_bins":n_bins}

    @modal.method()
    def health(self) -> dict:
        return {"status": "ok", "model": "enformer-pytorch", "torch_version": __import__("torch").__version__}


@app.local_entrypoint()
def test():
    import random; random.seed(99)
    seq = "".join(random.choices("ACGT", k=1000))
    print("\n=== CrisPRO Enformer Smoke Test ===")
    svc = EnformerService()
    r = svc.score_accessibility.remote(seq, "dnase")
    print(f"accessibility_score: {r['accessibility_score']:.6f}")
    print(f"torch_version via health: {svc.health.remote()['torch_version']}")
    print("✓ crispro-enformer live")
