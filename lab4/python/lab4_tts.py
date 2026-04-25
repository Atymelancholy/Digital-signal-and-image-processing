import argparse
import subprocess
import sys
from pathlib import Path


def _ensure_package(import_name: str, pip_name: str | None = None):
    try:
        __import__(import_name)
        return
    except Exception:
        pass
    name = pip_name or import_name
    subprocess.check_call([sys.executable, "-m", "pip", "install", name])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--text", required=True, help="Input text")
    parser.add_argument("--out", required=True, help="Output wav path")
    args = parser.parse_args()

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    _ensure_package("pyttsx3")
    import pyttsx3

    engine = pyttsx3.init()
    engine.setProperty("rate", 175)
    engine.save_to_file(args.text, str(out_path))
    engine.runAndWait()

    if not out_path.exists() or out_path.stat().st_size == 0:
        raise RuntimeError(f"TTS did not produce output file: {out_path}")

    print(f"OK: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

