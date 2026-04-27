"""Lab 4 — Text-to-Speech via the official CosyVoice project.

Two modes:

  * sft       — CosyVoice-300M-SFT with a fixed roster of preset voices
                (e.g. '中文女', '英文男').
  * zero_shot — CosyVoice-300M zero-shot voice cloning. Requires a short
                reference WAV plus the transcript of what is said in it.

Run with the dedicated venv created by setup_cosyvoice.bat:

    python\.cosyvenv\Scripts\python.exe python\lab4_tts.py --mode sft \
        --text "Привет, мир" --speaker "英文女" --out out.wav

The script auto-locates the CosyVoice repo at python/CosyVoice/ unless
--cosy-dir is provided.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _force_utf8_stdio() -> None:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")


def _resolve_cosy_dir(explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit).resolve()
        if not path.exists():
            raise FileNotFoundError(f"--cosy-dir not found: {path}")
        return path
    here = Path(__file__).resolve().parent
    candidate = here / "CosyVoice"
    if candidate.exists():
        return candidate
    raise FileNotFoundError(
        "CosyVoice repository not found. Run python\\setup_cosyvoice.bat first."
    )


def _enable_imports(cosy: Path) -> None:
    sys.path.insert(0, str(cosy))
    matcha = cosy / "third_party" / "Matcha-TTS"
    if matcha.exists():
        sys.path.insert(0, str(matcha))


def _load_model(model_dir: Path):
    from cosyvoice.cli.cosyvoice import AutoModel

    print(f"Loading CosyVoice from {model_dir} ...", flush=True)
    try:
        return AutoModel(model_dir=str(model_dir))
    except TypeError:
        return AutoModel(model_dir=str(model_dir))


def _save_first(generator, out_path: Path, sample_rate: int) -> None:
    import torchaudio

    out_path.parent.mkdir(parents=True, exist_ok=True)
    first = next(iter(generator))
    audio = first["tts_speech"]
    torchaudio.save(str(out_path), audio, sample_rate)


def main() -> int:
    _force_utf8_stdio()

    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["sft", "zero_shot"], default="sft")
    parser.add_argument("--text", required=True, help="Text to synthesise")
    parser.add_argument("--speaker", default="中文女",
                        help="(sft) preset speaker, e.g. 中文女/中文男/英文女/英文男/日语男/粤语女/韩语女")
    parser.add_argument("--prompt-audio", default=None, help="(zero_shot) reference WAV")
    parser.add_argument("--prompt-text", default=None, help="(zero_shot) transcript of the reference WAV")
    parser.add_argument("--out", required=True, help="Output WAV path")
    parser.add_argument("--cosy-dir", default=None, help="Path to cloned CosyVoice repo")
    args = parser.parse_args()

    cosy = _resolve_cosy_dir(args.cosy_dir)
    _enable_imports(cosy)

    if args.mode == "sft":
        model_dir = cosy / "pretrained_models" / "CosyVoice-300M-SFT"
    else:
        model_dir = cosy / "pretrained_models" / "CosyVoice-300M"

    if not model_dir.exists():
        raise FileNotFoundError(
            f"Model directory missing: {model_dir}. Run setup_cosyvoice.bat to download."
        )

    cosyvoice = _load_model(model_dir)
    out_path = Path(args.out).resolve()

    if args.mode == "sft":
        print(f"Mode: SFT, speaker={args.speaker}", flush=True)
        gen = cosyvoice.inference_sft(args.text, args.speaker, stream=False)
    else:
        if not args.prompt_audio or not args.prompt_text:
            raise ValueError("--prompt-audio and --prompt-text are required for zero_shot mode")
        prompt_path = Path(args.prompt_audio).resolve()
        if not prompt_path.exists():
            raise FileNotFoundError(f"--prompt-audio not found: {prompt_path}")
        print(f"Mode: zero-shot, prompt={prompt_path}", flush=True)
        gen = cosyvoice.inference_zero_shot(
            args.text, args.prompt_text, str(prompt_path), stream=False
        )

    _save_first(gen, out_path, cosyvoice.sample_rate)
    print(f"OK CosyVoice TTS: {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
