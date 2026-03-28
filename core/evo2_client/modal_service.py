"""
modal_evo2.py — CrisPRO.ai Evo2 Inference Service (v7 — field engineer architecture)

Key fixes from evo2_deployment_patch.mdc:
  1. Base image: nvidia/cuda:12.4.0-devel-ubuntu22.04 (has nvcc for flash-attn build)
  2. torch==2.3.0 (field-tested version)
  3. Global torch.load monkey-patch BEFORE evo2 import (patches.py apply_torch_patches)
  4. _FlashAttnGPUProxy to normalize fwd() return signature
  5. timeout=1800, volume at /root/.cache/huggingface

Deploy: modal deploy modal_evo2.py
"""
import modal
from pathlib import Path

# Model cache volume — 13.7GB evo2_7b.pt persists across cold starts
model_cache  = modal.Volume.from_name("crispro-model-cache", create_if_missing=True)
EVO_MODEL_DIR = "/root/.cache/huggingface"

# Field-engineer proven image — CUDA devel base so flash-attn can compile
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11"
    )
    .apt_install([
        "build-essential", "cmake", "ninja-build",
        "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"
    ])
    .env({
        "CC":  "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
        "CRISPRO_VERSION": "v9b",  # bump to force cold container
    })
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar",
    )
    .run_commands("pip install wheel setuptools")
    .run_commands(
        "pip install --no-cache-dir torch==2.3.0 torchvision torchaudio "
        "--index-url https://download.pytorch.org/whl/cu121"
    )
    .run_commands(
        "git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && "
        "cd evo2 && pip install ."
    )
    # TE 1.13 = last version supporting cuDNN 8 (our base image has cuDNN 8, not 9)
    # vortex imports transformer_engine at module level — cannot skip
    .run_commands("pip uninstall -y transformer-engine transformer_engine || true")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .run_commands(
        "bash -c \"pip install 'flash-attn-cu124==2.6.3' --no-build-isolation 2>/dev/null || "
        "pip install 'flash-attn-cu121==2.6.3' --no-build-isolation 2>/dev/null || "
        "pip install 'flash-attn==2.6.3' --no-build-isolation\""
    )
    .pip_install(
        "fastapi", "uvicorn[standard]", "loguru",
        "pydantic", "numpy", "httpx", "biopython", "requests", "pandas"
    )
    # Copy patches.py into image
    .add_local_file("/home/user/workspace/patches.py", "/patches.py", copy=True)
    # Ensure no stale pyc of modal_evo2.py in container
    .run_commands("find /root -name '*.pyc' -delete 2>/dev/null; find /root -name '__pycache__' -delete 2>/dev/null; true")
)

app = modal.App("crispro-evo2-v9")


@app.cls(
    image=evo2_image,
    gpu="A100",             # A100 = SM8.0, supports BF16; T4=SM7.5 does NOT support BF16
                            # Field engineer confirmed: H100:2 for FP8; A100 for BF16+no-FP8
    timeout=1800,           # 30min — field engineer proven value
    volumes={EVO_MODEL_DIR: model_cache},
    min_containers=0,
    scaledown_window=300,
)
class Evo2Service:

    # evo2_1b_base on A100 (SM8.0, BF16 ok, FP8 disabled via patches.py)
    # For production 7B: switch gpu="H100:2", model_name="evo2_7b"
    model_name: str = "evo2_1b_base"

    @modal.enter()
    def load(self):
        import sys
        sys.path.insert(0, "/")

        # CRITICAL: apply torch patches BEFORE importing evo2
        from patches import apply_torch_patches, apply_evo2_runtime_patches
        apply_torch_patches()
        print("[evo2] torch.load patched (weights_only=False)", flush=True)

        import torch
        from evo2 import Evo2
        apply_evo2_runtime_patches()

        print(f"[evo2] Loading {self.model_name}...", flush=True)
        self.model = Evo2(self.model_name)
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self._rc_map = str.maketrans("ACGTacgt", "TGCAtgca")
        print(f"[evo2] ✓ Loaded on {self.device}", flush=True)

    def _rc(self, seq):
        return seq.translate(self._rc_map)[::-1]

    def _ll(self, seq):
        # score_sequences returns mean LL per token — for full-sequence comparison
        scores = self.model.score_sequences(
            [seq], average_reverse_complement=True
        )
        return float(scores[0])

    def _ll_sum(self, seq):
        """Sum (not mean) LL — more sensitive for short variant windows."""
        import torch
        import torch.nn.functional as F
        # Tokenize
        tok = self.model.tokenizer
        ids = tok.tokenize(seq)
        if hasattr(ids, 'input_ids'):
            ids = ids.input_ids
        if not isinstance(ids, torch.Tensor):
            ids = torch.tensor(ids)
        ids = ids.to(self.device).unsqueeze(0)  # [1, L]
        with torch.inference_mode():
            out    = self.model.model(ids)
            logits = out[0] if isinstance(out, (tuple, list)) else out
        # logits: [1, L, V] — predict next token
        log_probs = F.log_softmax(logits[:, :-1, :], dim=-1)  # [1, L-1, V]
        target    = ids[:, 1:]                                  # [1, L-1]
        ll = log_probs[0, torch.arange(target.shape[1]), target[0]].sum().item()
        return ll

    def _tokenize(self, seq: str):
        """Tokenize a DNA sequence, returning a Long tensor on self.device.
        Evo2's byte-level tokenizer may return ByteTensor which must be cast."""
        import torch
        tok = self.model.tokenizer

        # Try the standard tokenizer interface
        try:
            ids = tok.tokenize(seq)
        except AttributeError:
            ids = tok(seq)

        if hasattr(ids, 'input_ids'):
            ids = ids.input_ids

        if not isinstance(ids, torch.Tensor):
            ids = torch.tensor(ids)

        # ByteTensor (uint8) cannot be used as embedding indices — cast to Long
        ids = ids.long().to(self.device)

        # Ensure 2D: [1, L]
        if ids.dim() == 1:
            ids = ids.unsqueeze(0)
        return ids

    def _conditional_ll_at_pos(self, seq: str, var_idx: int):
        """Compute log P(base at var_idx | preceding context) using autoregressive prediction.
        Uses Evo2's own prepare_batch + forward pass so the tensor dtype is correct.
        var_idx: 0-based index of the variant position in seq."""
        import torch
        import torch.nn.functional as F
        from evo2.scoring import prepare_batch

        # prepare_batch uses dtype=torch.long internally (same as Evo2 paper scoring)
        input_ids, _ = prepare_batch([seq], self.model.tokenizer, device=self.device)
        # input_ids: [1, L], Long

        with torch.inference_mode():
            # Evo2.forward() returns (logits, embeddings_dict_or_None)
            # StripedHyena.forward() may return just tensor or (tensor, ...)
            raw_out = self.model.forward(input_ids)

        # Unwrap: handle both (logits, extra) and plain logits
        if isinstance(raw_out, (tuple, list)):
            logits = raw_out[0]
        else:
            logits = raw_out

        # logits may still be a tuple if StripedHyena returns (tensor, state)
        if isinstance(logits, (tuple, list)):
            logits = logits[0]

        # logits: [1, L, V], dtype=bfloat16 or float32
        # logits[i] predicts token at position i+1 (causal LM convention)
        pred_pos = var_idx - 1
        if pred_pos < 0 or pred_pos >= logits.shape[1]:
            return None

        log_probs = F.log_softmax(logits[0, pred_pos, :].float(), dim=-1)  # [V]
        target_id = input_ids[0, var_idx].item()
        return log_probs[target_id].item()

    @modal.method()
    def score_variant(self, wt_seq: str, mut_seq: str) -> dict:
        """Mean LL scoring — fast, but diluted over full sequence length."""
        ll_wt  = self._ll(wt_seq)
        ll_mut = self._ll(mut_seq)
        delta  = ll_mut - ll_wt
        return {
            "delta_ll": delta,
            "ll_wt":    ll_wt,
            "ll_mut":   ll_mut,
            "interpretation": (
                "deleterious"        if delta < -5.0 else
                "likely_deleterious" if delta < -2.0 else
                "uncertain"          if abs(delta) < 2.0 else
                "likely_neutral"
            ),
            "model":   self.model_name,
            "method": "mean_ll",
        }

    @modal.method()
    def score_variant_conditional(self, wt_seq: str, mut_seq: str, var_idx: int) -> dict:
        """Position-specific conditional LL scoring — the correct Evo2 paper method.
        Computes log P(ref_base | context) - log P(alt_base | context) at var_idx.
        var_idx: 0-based position of the variant within the sequence."""
        ll_wt_cond  = self._conditional_ll_at_pos(wt_seq,  var_idx)
        ll_mut_cond = self._conditional_ll_at_pos(mut_seq, var_idx)

        if ll_wt_cond is None or ll_mut_cond is None:
            return {"error": "var_idx out of tokenized sequence range"}

        delta = ll_mut_cond - ll_wt_cond
        return {
            "delta_ll":          delta,
            "ll_wt_conditional": ll_wt_cond,
            "ll_mut_conditional":ll_mut_cond,
            "var_idx":           var_idx,
            "interpretation": (
                "deleterious"        if delta < -2.0 else
                "likely_deleterious" if delta < -1.0 else
                "uncertain"          if abs(delta) < 1.0 else
                "likely_neutral"
            ),
            "model":  self.model_name,
            "method": "conditional_ll",
        }

    @modal.method()
    def score_variants_batch(self, variants: list) -> list:
        results = []
        for v in variants:
            s = self.score_variant.local(v["wt_seq"], v["mut_seq"])
            s.update({k: v[k] for k in v if k not in ("wt_seq", "mut_seq")})
            results.append(s)
        return results

    @modal.method()
    def embed_sequence(self, sequence: str, layer: int = 26,
                       pooling: str = "mean") -> dict:
        import torch
        activations = {}
        def hook(m, i, o): activations["emb"] = o.detach().cpu()
        layers = list(self.model.model.layers)
        if layer >= len(layers):
            return {"error": f"Layer {layer} out of range"}
        h = layers[layer].register_forward_hook(hook)
        try:
            with torch.no_grad():
                _ = self.model.score_sequences([sequence], average_reverse_complement=False)
        finally: h.remove()
        emb = activations.get("emb")
        if emb is None: return {"error": "hook failed"}
        pooled = emb.mean(dim=1).squeeze().tolist() if pooling == "mean" else emb.squeeze().tolist()
        return {"embeddings": pooled, "shape": list(emb.shape),
                "layer": layer, "pooling": pooling, "model": self.model_name}

    @modal.method()
    def generate(self, prompt: str, max_length: int = 512,
                 temperature: float = 1.0, top_k: int = 4,
                 n_candidates: int = 1) -> dict:
        out = self.model.generate(
            prompt_seqs=[prompt] * n_candidates,
            n_tokens=max_length, temperature=temperature, top_k=top_k)
        seqs = out.sequences if hasattr(out, "sequences") else list(out)
        return {"sequences": [str(s) for s in seqs],
                "n_candidates": n_candidates, "model": self.model_name}

    @modal.method()
    def health(self) -> dict:
        import torch
        return {
            "status": "ok",
            "model":  self.model_name,
            "device": self.device,
            "cuda":   torch.cuda.is_available(),
            "torch":  torch.__version__,
        }


@app.local_entrypoint()
def test():
    import random
    random.seed(42)
    seq = "".join(random.choices("ACGT", k=500))
    mut = seq[:250] + ("T" if seq[250] != "T" else "A") + seq[251:]
    print("\n=== CrisPRO Evo2 v7 Test ===")
    svc = Evo2Service()
    h = svc.health.remote()
    print(f"Health: {h}")
    r = svc.score_variant.remote(seq, mut)
    print(f"delta_ll: {r['delta_ll']:.4f} [{r['interpretation']}]")
    print("✓ LIVE")
# v7b — score_sequences API fix Sat Mar 28 02:46:30 UTC 2026
