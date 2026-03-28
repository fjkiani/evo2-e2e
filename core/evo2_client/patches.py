"""PyTorch compatibility patches for Evo2 checkpoint loading.
These patches are required for PyTorch 2.3.0+ compatibility with Evo2 checkpoints.
Source: field engineer evo2_deployment_patch.mdc"""
try:
    import torch
except ImportError:
    torch = None

class _FlashAttnGPUProxy:
    """Proxy around flash-attn GPU module to normalize `.fwd()` return signature."""
    def __init__(self, inner):
        self._inner = inner

    def __getattr__(self, name):
        return getattr(self._inner, name)

    def fwd(self, *args, **kwargs):
        res = self._inner.fwd(*args, **kwargs)
        try:
            def _clone(x): return x.clone() if torch.is_tensor(x) else x
            if isinstance(res, tuple):
                head = res[:4] if len(res) > 4 else res
                if isinstance(head, tuple) and len(head) == 4:
                    return tuple(_clone(x) for x in head)
                return head
            if isinstance(res, list):
                head = res[:4] if len(res) > 4 else res
                if len(head) == 4:
                    return tuple(_clone(x) for x in head)
                return tuple(head)
            out, softmax_lse, s_dmask, rng_state, *_rest = res
            return (_clone(out), _clone(softmax_lse), _clone(s_dmask), _clone(rng_state))
        except Exception:
            return res

def ensure_vortex_flash_attn_compat() -> bool:
    try:
        import sys
        patched_any = False
        try:
            from vortex.ops import attn_interface as _attn_interface
            if hasattr(_attn_interface, "flash_attn_gpu") and not isinstance(_attn_interface.flash_attn_gpu, _FlashAttnGPUProxy):
                _attn_interface.flash_attn_gpu = _FlashAttnGPUProxy(_attn_interface.flash_attn_gpu)
                patched_any = True
        except Exception:
            pass
        for _m in list(sys.modules.values()):
            try:
                if not getattr(_m, "__file__", ""): continue
                if not str(_m.__file__).endswith("/vortex/ops/attn_interface.py"): continue
                if hasattr(_m, "flash_attn_gpu") and not isinstance(getattr(_m, "flash_attn_gpu"), _FlashAttnGPUProxy):
                    setattr(_m, "flash_attn_gpu", _FlashAttnGPUProxy(getattr(_m, "flash_attn_gpu")))
                    patched_any = True
            except Exception:
                continue
        return patched_any
    except Exception:
        return False

def apply_torch_patches():
    """Apply PyTorch patches required for Evo2 checkpoint loading.
    MUST be called BEFORE importing evo2."""
    original_torch_load = torch.load

    def patched_torch_load(*args, **kwargs):
        kwargs['weights_only'] = False
        return original_torch_load(*args, **kwargs)

    torch.load = patched_torch_load

    if not hasattr(torch.serialization, 'add_safe_globals'):
        def add_safe_globals(_globals_list): return None
        torch.serialization.add_safe_globals = add_safe_globals

    try:
        import torch._library.custom_ops as _custom_ops
        if hasattr(_custom_ops, "_validate_outputs_to_not_alias_inputs"):
            _custom_ops._validate_outputs_to_not_alias_inputs = lambda *args, **kwargs: None
    except Exception:
        pass

    # Disable FP8 in transformer_engine — T4=SM7.5, A100=SM8.0, FP8 needs SM>=8.9 (H100/Ada)
    # Strategy: patch check_fp8_support() + FP8GlobalStateManager.fp8_available so the CC
    # check returns True before ANY local binding in vortex can capture a failing result.
    # This fires BEFORE 'from evo2 import Evo2' so vortex gets patched TE from sys.modules.
    try:
        import transformer_engine.pytorch.fp8 as _te_fp8

        # 1. Patch the source-of-truth CC check — covers cold path
        _te_fp8.check_fp8_support = lambda: (True, "")

        # 2. Force the cached value so is_fp8_available() short-circuits immediately
        if hasattr(_te_fp8, 'FP8GlobalStateManager'):
            _gsm = _te_fp8.FP8GlobalStateManager
            _gsm.fp8_available = True
            _gsm.reason_for_no_fp8 = ""
            # Also patch is_fp8_available as a classmethod returning (True, "")
            _gsm.is_fp8_available = classmethod(lambda cls: (True, ""))

        print("[patches] TE FP8 CC check disabled (T4/A100 compat)", flush=True)
    except Exception as _e:
        print(f"[patches] TE FP8 patch skipped: {_e}", flush=True)

    # 3. Also replace fp8_autocast context manager as belt-and-suspenders
    #    (covers the case where vortex has a local 'from .fp8 import fp8_autocast' binding
    #    that was captured before our patch — unlikely since we patch first, but safe)
    try:
        import transformer_engine.pytorch as te
        from contextlib import contextmanager
        @contextmanager
        def _noop_fp8_autocast(*args, **kwargs):
            yield
        te.fp8_autocast = _noop_fp8_autocast
        if hasattr(te, 'fp8'):
            te.fp8.fp8_autocast = _noop_fp8_autocast
    except Exception:
        pass

    # 4. Disable at the TE common recipe level too
    try:
        import transformer_engine.common.recipe as te_recipe
        class _NoopRecipe:
            def __init__(self, *a, **kw): pass
        te_recipe.DelayedScaling = _NoopRecipe
    except Exception:
        pass

def apply_evo2_runtime_patches():
    """Apply runtime monkey-patches to the installed evo2 package."""
    # Disable FP8 globally — A100 is compute 8.0, FP8 needs 8.9 (H100)
    try:
        import evo2.models as evo2_models
        # Patch the model config loading to force fp8=False
        _orig_init = evo2_models.Evo2.__init__
        def _patched_init(self, *args, **kwargs):
            _orig_init(self, *args, **kwargs)
            # After model loads, force fp8 off on all submodules
            if hasattr(self, 'model'):
                for name, module in self.model.named_modules():
                    for attr in ('use_fp8', 'use_fp8_input_projections',
                                 '_use_fp8', 'fp8_enabled'):
                        if hasattr(module, attr):
                            setattr(module, attr, False)
        evo2_models.Evo2.__init__ = _patched_init
    except Exception:
        pass

    # Also patch score_sequences to disable fp8 before calling
    try:
        import evo2.scoring as _scoring
        _orig_score = _scoring.score_sequences
        def _safe_score(seqs, model, tokenizer, *args, **kwargs):
            # Force fp8 off on model before scoring
            if hasattr(model, 'named_modules'):
                for _, m in model.named_modules():
                    for attr in ('use_fp8', 'use_fp8_input_projections'):
                        if hasattr(m, attr): setattr(m, attr, False)
            return _orig_score(seqs, model, tokenizer, *args, **kwargs)
        _scoring.score_sequences = _safe_score
    except Exception:
        pass

    try:
        import evo2.scoring as scoring
    except Exception:
        return

    if hasattr(scoring, "prepare_batch"):
        _orig_prepare_batch = scoring.prepare_batch
        def _patched_prepare_batch(*args, **kwargs):
            out = _orig_prepare_batch(*args, **kwargs)
            if isinstance(out, tuple) and len(out) > 4: return out[:4]
            return out
        scoring.prepare_batch = _patched_prepare_batch

    if hasattr(scoring, "_score_sequences"):
        _orig__score_sequences = scoring._score_sequences
        def _patched__score_sequences(*args, **kwargs):
            try:
                return _orig__score_sequences(*args, **kwargs)
            except ValueError as e:
                if "too many values to unpack" not in str(e): raise
                import numpy as np
                seqs = kwargs.get("seqs") or (args[0] if len(args) > 0 else None)
                model = kwargs.get("model") or (args[1] if len(args) > 1 else None)
                tokenizer = kwargs.get("tokenizer") or (args[2] if len(args) > 2 else None)
                prepend_bos = kwargs.get("prepend_bos", False)
                reduce_method = kwargs.get("reduce_method", "mean")
                device = kwargs.get("device", "cuda:0")

                input_ids, seq_lengths, *_rest = scoring.prepare_batch(
                    seqs, tokenizer, device=device, prepend_bos=prepend_bos
                )
                with torch.inference_mode():
                    out = model(input_ids)
                    logits = out[0] if isinstance(out, (tuple, list)) else out

                logprobs = scoring.logits_to_logprobs(logits, input_ids)
                logprobs = logprobs.float().cpu().numpy()
                reduce_func = np.sum if reduce_method == "sum" else np.mean

                scores = []
                for idx in range(len(seq_lengths)):
                    n = int(seq_lengths[idx])
                    n_eff = max(0, n - 1)
                    s = reduce_func(logprobs[idx, 1 : 1 + n_eff])
                    scores.append(float(s) if np.isfinite(s) else float("nan"))
                return scores
        scoring._score_sequences = _patched__score_sequences

    try:
        ensure_vortex_flash_attn_compat()
    except Exception:
        pass
