"""One-shot installer for the Lab 4 environment.

Runs inside the dedicated `.cosyvenv` virtualenv (created by setup_cosyvoice.bat).
Performs the heavy lifting:

  1. Installs PyTorch + torchaudio with CUDA 12.1 wheels (or falls back to CPU).
  2. Clones the FunAudioLLM/CosyVoice repository (with submodules).
  3. Installs CosyVoice runtime dependencies that are Windows-friendly.
     (deepspeed / vllm / pynini are skipped: not available on Windows pip.)
  4. Downloads CosyVoice-300M and CosyVoice-300M-SFT model checkpoints
     from ModelScope into CosyVoice/pretrained_models/.

After this script finishes:

  - python/.cosyvenv/Scripts/python.exe  -> dedicated CUDA Python
  - python/CosyVoice/                    -> upstream repo
  - python/CosyVoice/pretrained_models/  -> downloaded models

Re-running the script is safe: every step is idempotent.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

PYDIR = Path(__file__).resolve().parent
COSY = PYDIR / "CosyVoice"
COSY_REPO = "https://github.com/FunAudioLLM/CosyVoice.git"
TORCH_CUDA_INDEX = "https://download.pytorch.org/whl/cu121"

# Direct CDN URLs for torch + torchaudio CUDA 12.1 wheels (Python 3.11, Windows x64).
# We download these via curl with resume support to dodge `pip` getting stuck on
# 2 GB streams (Windows TCP error 10013 on flaky links).
TORCH_VERSION = "2.5.1+cu121"
TORCHAUDIO_VERSION = "2.5.1+cu121"
TORCH_WHEELS = {
    f"torch-{TORCH_VERSION}-cp311-cp311-win_amd64.whl":
        f"https://download.pytorch.org/whl/cu121/torch-2.5.1%2Bcu121-cp311-cp311-win_amd64.whl",
    f"torchaudio-{TORCHAUDIO_VERSION}-cp311-cp311-win_amd64.whl":
        f"https://download.pytorch.org/whl/cu121/torchaudio-2.5.1%2Bcu121-cp311-cp311-win_amd64.whl",
}
WHEEL_CACHE = PYDIR / "wheel_cache"

# Windows-friendly subset of CosyVoice dependencies. Anything that requires
# pynini / WeTextProcessing is intentionally skipped; lab4_tts.py runs the
# model with text_frontend=False to bypass that path.
PIP_PACKAGES = [
    "numpy<2.0",
    "scipy",
    "librosa",
    "soundfile",
    "modelscope",
    "huggingface_hub",
    "transformers",
    "accelerate",
    "diffusers>=0.27",
    "omegaconf",
    "hydra-core",
    "hyperpyyaml",
    "conformer",
    "einops",
    "matplotlib",
    "tensorboard",
    "tqdm",
    "wget",
    "inflect",
    "protobuf",
    "pyarrow",
    "pyworld",
    "onnx",
    "onnxruntime",
    "lightning",
    "rich",
    "gdown",
    "openai-whisper",
]


def run(cmd: list[str], check: bool = True, env: dict | None = None) -> int:
    print("$", " ".join(str(part) for part in cmd), flush=True)
    result = subprocess.run(cmd, env=env)
    if check and result.returncode != 0:
        raise SystemExit(result.returncode)
    return result.returncode


def _curl_download(url: str, dest: Path) -> bool:
    """Download `url` to `dest` using curl with infinite retries and resume.

    Returns True on success, False on failure. The destination directory must
    already exist. Uses `-C -` (continue-at), `--retry`, and `--retry-all-errors`
    so transient WinError 10013 / 10054 / EOF errors don't kill the run.
    """
    if not shutil.which("curl"):
        print("[curl] curl is not in PATH, cannot do resumable download.", flush=True)
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "curl",
        "-L",                       # follow redirects
        "-C", "-",                  # resume
        "--retry", "999",
        "--retry-delay", "5",
        "--retry-max-time", "0",
        "--retry-all-errors",
        "--connect-timeout", "30",
        "--continue-at", "-",       # resume any partial bytes
        "-o", str(dest),
        url,
    ]
    print("$", " ".join(cmd), flush=True)
    # We allow many retries internally; if curl returns non-zero we fall through.
    rc = subprocess.call(cmd)
    # torchaudio wheels are only ~4 MB; the lower bound is just a sanity check.
    return rc == 0 and dest.exists() and dest.stat().st_size > 1 * 1024 * 1024


def _ensure_torch_wheels() -> list[Path]:
    WHEEL_CACHE.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for filename, url in TORCH_WHEELS.items():
        dest = WHEEL_CACHE / filename
        if dest.exists() and dest.stat().st_size > 1 * 1024 * 1024:
            print(f"[torch] Already have {dest.name} ({dest.stat().st_size/1024/1024:.1f} MB), skipping download.", flush=True)
            paths.append(dest)
            continue
        print(f"\n[torch] Downloading {filename} via curl (resumable)...", flush=True)
        ok = _curl_download(url, dest)
        if not ok:
            print(f"[torch] curl download failed for {filename}", flush=True)
            return []
        size_mb = dest.stat().st_size / 1024 / 1024
        print(f"[torch] {filename} downloaded ({size_mb:.0f} MB).", flush=True)
        paths.append(dest)
    return paths


def step_torch_cuda() -> None:
    print("\n[torch] Installing torch + torchaudio (CUDA 12.1 wheels)...", flush=True)
    wheels = _ensure_torch_wheels()
    if wheels:
        cmd = [sys.executable, "-m", "pip", "install", "--upgrade", *(str(p) for p in wheels)]
        rc = run(cmd, check=False)
        if rc != 0:
            print("[torch] Local-wheel install failed, trying pip with --index-url cu121...", flush=True)
            rc = run(
                [sys.executable, "-m", "pip", "install", "--upgrade",
                 "torch", "torchaudio", "--index-url", TORCH_CUDA_INDEX],
                check=False,
            )
    else:
        print("[torch] Falling back to direct pip install --index-url cu121...", flush=True)
        rc = run(
            [sys.executable, "-m", "pip", "install", "--upgrade",
             "torch", "torchaudio", "--index-url", TORCH_CUDA_INDEX],
            check=False,
        )

    if rc != 0:
        print("[torch] CUDA wheel install failed, falling back to CPU torch...", flush=True)
        run([sys.executable, "-m", "pip", "install", "--upgrade", "torch", "torchaudio"])

    try:
        import torch  # noqa: F401
        print(f"[torch] version={torch.__version__}, cuda_available={torch.cuda.is_available()}", flush=True)
        if torch.cuda.is_available():
            print(f"[torch] GPU: {torch.cuda.get_device_name(0)}", flush=True)
    except Exception as exc:
        print(f"[torch] WARNING: could not import torch after install: {exc}", flush=True)


def step_clone_cosy() -> None:
    print("\n[clone] Ensuring CosyVoice repository is present...", flush=True)
    if not shutil.which("git"):
        print("ERROR: git is not in PATH. Install Git for Windows and retry.", flush=True)
        raise SystemExit(2)

    if not COSY.exists():
        run(["git", "clone", "--recursive", COSY_REPO, str(COSY)])
    else:
        print(f"[clone] {COSY} exists, pulling latest...", flush=True)
        run(["git", "-C", str(COSY), "pull", "--ff-only"], check=False)
        run(["git", "-C", str(COSY), "submodule", "update", "--init", "--recursive"], check=False)

    matcha = COSY / "third_party" / "Matcha-TTS"
    if not matcha.exists():
        print("WARNING: third_party/Matcha-TTS missing; cloning manually...", flush=True)
        run(["git", "clone", "https://github.com/shivammehta25/Matcha-TTS.git", str(matcha)], check=False)


def step_install_deps() -> None:
    print("\n[deps] Installing Windows-safe runtime dependencies...", flush=True)
    run([sys.executable, "-m", "pip", "install", "--upgrade", *PIP_PACKAGES])


def step_download_models() -> None:
    print("\n[models] Downloading CosyVoice checkpoints (ModelScope)...", flush=True)
    pretrained = COSY / "pretrained_models"
    pretrained.mkdir(parents=True, exist_ok=True)

    targets = [
        ("iic/CosyVoice-300M", pretrained / "CosyVoice-300M"),
        ("iic/CosyVoice-300M-SFT", pretrained / "CosyVoice-300M-SFT"),
    ]
    for repo_id, local_dir in targets:
        marker = local_dir / "cosyvoice.yaml"
        if marker.exists():
            print(f"[models] {repo_id} already downloaded -> {local_dir}", flush=True)
            continue
        print(f"[models] Downloading {repo_id} -> {local_dir}", flush=True)
        run([
            sys.executable,
            "-c",
            (
                "from modelscope import snapshot_download; "
                f"snapshot_download(model_id='{repo_id}', local_dir=r'{local_dir}')"
            ),
        ])


def main() -> int:
    if sys.version_info[:2] < (3, 10):
        print("ERROR: Python 3.10 or newer is required (running %s)." % sys.version, flush=True)
        return 1

    in_venv = (
        getattr(sys, "base_prefix", sys.prefix) != sys.prefix
        or os.environ.get("VIRTUAL_ENV")
    )
    if not in_venv:
        print("WARNING: this script is meant to run inside .cosyvenv (use setup_cosyvoice.bat).", flush=True)

    step_torch_cuda()
    step_clone_cosy()
    step_install_deps()
    step_download_models()
    print("\nAll steps finished.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
