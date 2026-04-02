import argparse
import json
import os
import sys


def _load_mono(path: str):
    import soundfile as sf
    y, sr = sf.read(path, always_2d=False)
    if hasattr(y, "ndim") and y.ndim > 1:
        y = y.mean(axis=1)
    return y, int(sr)


def _ensure_same_length(a, b):
    n = min(len(a), len(b))
    return a[:n], b[:n]


def _pesq_mos(ref, deg, sr):
    # PESQ supports only 8k (nb) or 16k (wb). We'll resample to 16k and use wideband.
    import numpy as np
    import librosa
    from pesq import pesq

    target_sr = 16000
    if sr != target_sr:
        ref = librosa.resample(ref.astype(np.float32), orig_sr=sr, target_sr=target_sr)
        deg = librosa.resample(deg.astype(np.float32), orig_sr=sr, target_sr=target_sr)
        sr = target_sr
    ref, deg = _ensure_same_length(ref, deg)
    if len(ref) < sr // 10:  # too short
        return None
    # pesq expects float32 in [-1..1]
    ref = ref.astype(np.float32)
    deg = deg.astype(np.float32)
    return float(pesq(sr, ref, deg, "wb"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True, help="reference clean wav")
    p.add_argument("--deg", required=True, help="degraded/noisy/processed wav")
    p.add_argument("--out", required=True, help="output json path")
    args = p.parse_args()

    payload = {
        "snr_db": None,
        "sdr_db": None,
        "si_sdr_db": None,
        "pesq_mos": None,
        "error": None,
    }

    try:
        import numpy as np
        import torch
        from torchmetrics.audio import (
            SignalNoiseRatio,
            SignalDistortionRatio,
            ScaleInvariantSignalDistortionRatio,
        )

        ref, sr_ref = _load_mono(args.ref)
        deg, sr_deg = _load_mono(args.deg)
        if sr_ref != sr_deg:
            # Keep TorchMetrics on the same sample rate: resample deg to ref sr.
            import librosa
            deg = librosa.resample(deg.astype(np.float32), orig_sr=sr_deg, target_sr=sr_ref)
            sr_deg = sr_ref

        ref, deg = _ensure_same_length(ref, deg)
        if len(ref) == 0:
            raise RuntimeError("Empty audio after alignment")

        x = torch.from_numpy(deg.astype(np.float32)).unsqueeze(0)
        y = torch.from_numpy(ref.astype(np.float32)).unsqueeze(0)

        # TorchMetrics return tensors
        payload["snr_db"] = float(SignalNoiseRatio()(x, y).item())
        payload["sdr_db"] = float(SignalDistortionRatio()(x, y).item())
        payload["si_sdr_db"] = float(ScaleInvariantSignalDistortionRatio()(x, y).item())

        try:
            payload["pesq_mos"] = _pesq_mos(ref, deg, sr_ref)
        except Exception:
            payload["pesq_mos"] = None
    except Exception as e:
        payload["error"] = str(e)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False)

    if payload["error"]:
        print(payload["error"], file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

