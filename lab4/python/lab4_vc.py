"""Lab 4 — Voice Conversion.

Two modes (selected with --mode):

  * knn  — official kNN-VC (Voice Conversion With Just k-Nearest Neighbors)
           loaded via torch.hub from bshall/knn-vc:master. WavLM features +
           HiFi-GAN vocoder. Optional --vocoder-ckpt to use a fine-tuned
           HiFi-GAN produced by python\lab4_finetune_colab.ipynb.

  * cosy — CosyVoice-300M zero-shot synthesis used as a voice conversion
           baseline (per the textbook's chapter 2.4). Requires source text
           and target text to be supplied.

Run inside the dedicated venv created by setup_cosyvoice.bat:

    python\.cosyvenv\Scripts\python.exe python\lab4_vc.py --mode knn \
        --src source.wav --tgt target.wav --out out.wav

The source's spoken content is preserved (kNN-VC operates on acoustic
features, not phonemes), and the timbre/voice colour is taken from the
target speaker reference.
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path


# ---------------------------------------------------------------------------
# kNN-VC backend
# ---------------------------------------------------------------------------

def _load_knn_vc(device: str, vocoder_ckpt: str | None):
    import torch

    print(f"Loading kNN-VC (bshall/knn-vc:master) on {device} ...", flush=True)
    knn_vc = torch.hub.load(
        "bshall/knn-vc:master", "knn_vc",
        prematched=True, trust_repo=True, pretrained=True, device=device,
    )

    if vocoder_ckpt:
        ckpt_path = Path(vocoder_ckpt).resolve()
        if not ckpt_path.exists():
            raise FileNotFoundError(f"Fine-tuned vocoder not found: {ckpt_path}")
        print(f"Loading fine-tuned HiFi-GAN weights: {ckpt_path}", flush=True)
        state = torch.load(str(ckpt_path), map_location=device)
        if isinstance(state, dict) and "generator" in state:
            state = state["generator"]
        knn_vc.hifigan.load_state_dict(state)
    return knn_vc


def _to_mono_16k(wav_path: str, device: str):
    import torch
    import torchaudio

    wav, sr = torchaudio.load(wav_path)
    if sr != 16000:
        wav = torchaudio.functional.resample(wav, sr, 16000)
    if wav.shape[0] > 1:
        wav = wav.mean(dim=0, keepdim=True)
    return wav.to(device)


def _align_TD(query, matching):
    """Make sure both tensors are shape (T, D) before kNN matching."""
    q = query.squeeze()
    m = matching.squeeze()
    if q.dim() == 3:
        q = q.reshape(-1, q.shape[-1])
    if m.dim() == 3:
        m = m.reshape(-1, m.shape[-1])

    if q.shape[-1] == m.shape[-1]:
        return q, m
    if q.shape[-1] == m.shape[0]:
        return q, m.T
    if q.shape[0] == m.shape[-1]:
        return q.T, m
    if q.shape[0] == m.shape[0]:
        return q.T, m.T
    if q.shape[-1] > 2000:
        q = q.T
    if m.shape[-1] > 2000:
        m = m.T
    return q, m


def run_knn(args) -> None:
    import torch
    import torchaudio

    device = "cuda" if torch.cuda.is_available() else "cpu"
    knn_vc = _load_knn_vc(device, args.vocoder_ckpt)

    src_wav = _to_mono_16k(args.src, device)
    tgt_wav = _to_mono_16k(args.tgt, device)

    if args.len_sec and args.len_sec > 0:
        n = max(1, int(args.len_sec * 16000))
        tgt_wav = tgt_wav[:, :n]

    tmp_ref = None
    try:
        if args.len_sec and args.len_sec > 0:
            tmp_ref = tempfile.NamedTemporaryFile(suffix=".wav", delete=False)
            tmp_ref.close()
            torchaudio.save(tmp_ref.name, tgt_wav.cpu(), 16000)
            ref_paths = [tmp_ref.name]
        else:
            ref_paths = [args.tgt]

        query_seq = knn_vc.get_features(src_wav)
        matching_set = knn_vc.get_matching_set(ref_paths)
        query_seq, matching_set = _align_TD(query_seq, matching_set)

        out_audio = knn_vc.match(query_seq, matching_set, topk=args.topk)
    finally:
        if tmp_ref is not None:
            try:
                Path(tmp_ref.name).unlink(missing_ok=True)
            except Exception:
                pass

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    torchaudio.save(str(out_path), out_audio.unsqueeze(0).cpu(), 16000)
    print(f"OK kNN-VC: {out_path}", flush=True)


# ---------------------------------------------------------------------------
# CosyVoice zero-shot backend (used for the comparison experiment, p.2.4)
# ---------------------------------------------------------------------------

def _resolve_cosy(explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit).resolve()
        if not path.exists():
            raise FileNotFoundError(f"--cosy-dir not found: {path}")
        return path
    candidate = Path(__file__).resolve().parent / "CosyVoice"
    if candidate.exists():
        return candidate
    raise FileNotFoundError(
        "CosyVoice repository not found. Run python\\setup_cosyvoice.bat first."
    )


def run_cosy_zero_shot(args) -> None:
    if not args.src_text or not args.tgt_text:
        raise ValueError("--src-text and --tgt-text are required for cosy mode")

    cosy = _resolve_cosy(args.cosy_dir)
    sys.path.insert(0, str(cosy))
    matcha = cosy / "third_party" / "Matcha-TTS"
    if matcha.exists():
        sys.path.insert(0, str(matcha))

    import torchaudio
    from cosyvoice.cli.cosyvoice import AutoModel

    model_dir = cosy / "pretrained_models" / "CosyVoice-300M"
    if not model_dir.exists():
        raise FileNotFoundError(
            f"Model directory missing: {model_dir}. Run setup_cosyvoice.bat to download."
        )

    print(f"Loading CosyVoice from {model_dir} ...", flush=True)
    cosyvoice = AutoModel(model_dir=str(model_dir))

    print("Synthesising in zero-shot mode (source words spoken in target voice) ...", flush=True)
    gen = cosyvoice.inference_zero_shot(
        args.src_text, args.tgt_text, args.tgt, stream=False,
    )
    first = next(iter(gen))
    audio = first["tts_speech"]

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    torchaudio.save(str(out_path), audio, cosyvoice.sample_rate)
    print(f"OK CosyVoice zero-shot VC: {out_path}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["knn", "cosy"], default="knn")
    parser.add_argument("--src", required=True, help="Source WAV path")
    parser.add_argument("--tgt", required=True, help="Target speaker reference WAV path")
    parser.add_argument("--out", required=True, help="Output WAV path")
    parser.add_argument("--topk", type=int, default=4, help="(knn) k for kNN matching")
    parser.add_argument("--vocoder-ckpt", default=None,
                        help="(knn) optional path to a fine-tuned HiFi-GAN .pt")
    parser.add_argument("--len-sec", type=float, default=None,
                        help="Trim target reference to this many seconds before VC")
    parser.add_argument("--cosy-dir", default=None,
                        help="(cosy) path to cloned CosyVoice repo")
    parser.add_argument("--src-text", default=None,
                        help="(cosy) text spoken in the source recording")
    parser.add_argument("--tgt-text", default=None,
                        help="(cosy) text spoken in the target/prompt recording")
    args = parser.parse_args()

    if args.mode == "knn":
        run_knn(args)
    else:
        run_cosy_zero_shot(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
