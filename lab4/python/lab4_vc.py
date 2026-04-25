import argparse
import importlib.metadata
import shutil
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


def _ensure_distribution(distribution_name: str, pip_name: str | None = None):
    try:
        importlib.metadata.version(distribution_name)
        return
    except importlib.metadata.PackageNotFoundError:
        pass
    subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name or distribution_name])


def _ensure_speechbrain_legacy():
    try:
        version = importlib.metadata.version("speechbrain")
        parts = tuple(int(part) for part in version.split(".")[:2])
        if parts < (1, 0):
            return
    except Exception:
        pass
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "speechbrain==0.5.16"])


def _ensure_legacy_datasets():
    try:
        import datasets
        version = tuple(int(part) for part in datasets.__version__.split(".")[:2])
        if version < (4, 0):
            return
    except Exception:
        pass
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "datasets==3.6.0"])


def _ensure_deps():
    _ensure_package("numpy")
    _ensure_package("librosa")
    _ensure_package("soundfile")


def _ensure_neural_deps():
    _ensure_deps()
    _ensure_package("torch")
    _ensure_package("transformers")
    # SpeechBrain 1.x lazy-loads optional k2 modules on this Windows/Python setup.
    # The x-vector recipe works reliably with the 0.5.x API.
    _ensure_speechbrain_legacy()
    _ensure_package("sentencepiece")


def _try_pyworld():
    try:
        import pyworld  # noqa: F401
        return True
    except Exception:
        pass
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyworld"])
        import pyworld  # noqa: F401
        return True
    except Exception as exc:
        print(f"pyworld unavailable, will use librosa fallback ({exc})", file=sys.stderr)
        return False


def _create_speaker_embedding(waveform, sr, device, mode):
    import torch
    import huggingface_hub
    from speechbrain.pretrained import EncoderClassifier

    # SpeechBrain 0.5.x passes the old use_auth_token argument, while newer
    # huggingface_hub expects token. Keep this compatibility local to this run.
    original_download = huggingface_hub.hf_hub_download

    def compat_hf_hub_download(*args, use_auth_token=None, **kwargs):
        if use_auth_token is not None and "token" not in kwargs:
            kwargs["token"] = use_auth_token
        try:
            return original_download(*args, **kwargs)
        except Exception:
            filename = kwargs.get("filename")
            if filename is None and len(args) >= 2:
                filename = args[1]
            if filename == "custom.py":
                empty_custom = Path.home() / ".cache" / "cosi_lab4" / "empty_custom.py"
                empty_custom.parent.mkdir(parents=True, exist_ok=True)
                empty_custom.write_text("# optional SpeechBrain custom hooks are not used here\n", encoding="utf-8")
                return str(empty_custom)
            raise

    huggingface_hub.hf_hub_download = compat_hf_hub_download

    original_symlink_to = Path.symlink_to

    def copy_instead_of_symlink(self, target, target_is_directory=False):
        src = Path(target)
        self.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            if self.exists():
                return None
            shutil.copytree(src, self)
        else:
            shutil.copy2(src, self)
        return None

    Path.symlink_to = copy_instead_of_symlink

    try:
        classifier = EncoderClassifier.from_hparams(
            source="speechbrain/spkrec-xvect-voxceleb",
            savedir=str(Path.home() / ".cache" / "cosi_lab4" / "spkrec-xvect-voxceleb"),
            run_opts={"device": device},
        )
    finally:
        Path.symlink_to = original_symlink_to

    with torch.no_grad():
        signal = torch.tensor(waveform, dtype=torch.float32, device=device).unsqueeze(0)
        embedding = classifier.encode_batch(signal)
        embedding = torch.nn.functional.normalize(embedding, dim=2)
        embedding = embedding.squeeze(0).cpu()
    print("Speaker embedding: extracted from uploaded target WAV", flush=True)
    return embedding.to(device)


def _speecht5_vc(src_path, tgt_path, out_path, mode, len_sec=None, source_anchor=0.18):
    _ensure_neural_deps()

    import numpy as np
    import torch
    import librosa
    import soundfile as sf
    from transformers import SpeechT5ForSpeechToSpeech, SpeechT5HifiGan, SpeechT5Processor

    sr = 16000
    source, _ = librosa.load(src_path, sr=sr, mono=True)
    target, _ = librosa.load(tgt_path, sr=sr, mono=True)
    if len_sec is not None and len_sec > 0:
        target = target[: max(1, int(len_sec * sr))]

    if source.size == 0:
        raise RuntimeError("Source audio is empty")
    if target.size < sr // 4:
        raise RuntimeError("Target audio is too short. Use at least 0.25 second.")

    print(
        f"Source WAV: {src_path} ({len(source) / sr:.2f}s), "
        f"Target WAV: {tgt_path} ({len(target) / sr:.2f}s)",
        flush=True,
    )

    # SpeechT5 is much more stable on phrase-length chunks than on very long recordings.
    max_chunk = sr * 12
    chunks = [source[i:i + max_chunk] for i in range(0, len(source), max_chunk)]

    device = "cuda" if torch.cuda.is_available() else "cpu"
    processor = SpeechT5Processor.from_pretrained("microsoft/speecht5_vc")
    model = SpeechT5ForSpeechToSpeech.from_pretrained("microsoft/speecht5_vc").to(device)
    vocoder = SpeechT5HifiGan.from_pretrained("microsoft/speecht5_hifigan").to(device)
    model.eval()
    vocoder.eval()

    speaker_embeddings = _create_speaker_embedding(target, sr, device, mode)
    generated = []
    with torch.no_grad():
        for chunk in chunks:
            if len(chunk) < sr // 4:
                continue
            inputs = processor(audio=chunk, sampling_rate=sr, return_tensors="pt")
            input_values = inputs["input_values"].to(device)
            speech = model.generate_speech(input_values, speaker_embeddings, vocoder=vocoder)
            generated.append(speech.detach().cpu().numpy())

    if not generated:
        raise RuntimeError("SpeechT5 did not generate audio")

    out = np.concatenate(generated).astype(np.float32)
    out = _anchor_to_source(out, source, source_anchor)
    peak = np.max(np.abs(out)) + 1e-9
    if peak > 1.0:
        out = out / peak

    Path(out_path).resolve().parent.mkdir(parents=True, exist_ok=True)
    sf.write(str(out_path), out, sr)


def _anchor_to_source(generated, source, amount):
    import numpy as np
    import librosa

    amount = float(np.clip(amount, 0.0, 0.5))
    if amount <= 0.0 or generated.size == 0 or source.size == 0:
        return generated

    if len(source) != len(generated):
        rate = len(source) / max(1, len(generated))
        try:
            anchored = librosa.effects.time_stretch(source.astype(np.float32), rate=rate)
        except Exception:
            anchored = source
        if len(anchored) < len(generated):
            anchored = np.pad(anchored, (0, len(generated) - len(anchored)))
        anchored = anchored[:len(generated)]
    else:
        anchored = source

    out = (1.0 - amount) * generated + amount * anchored.astype(np.float32)
    peak = np.max(np.abs(out)) + 1e-9
    if peak > 1.0:
        out = out / peak
    print(f"Source anchor mix: {amount:.2f}", flush=True)
    return out.astype(np.float32)


def _robust_median_pitch(y, sr, fmin=70.0, fmax=350.0):
    import numpy as np
    import librosa

    f0 = librosa.yin(y, fmin=fmin, fmax=fmax, sr=sr, frame_length=2048, hop_length=256)
    valid = np.isfinite(f0) & (f0 > 0.0)
    if not np.any(valid):
        return None
    return float(np.median(f0[valid]))


def _transfer_pitch(src, tgt, sr, max_abs_semitones=7.0):
    import librosa
    import numpy as np

    src_f0 = _robust_median_pitch(src, sr)
    tgt_f0 = _robust_median_pitch(tgt, sr)
    if src_f0 is None or tgt_f0 is None:
        return src
    semitones = 12.0 * np.log2((tgt_f0 + 1e-9) / (src_f0 + 1e-9))
    semitones = float(np.clip(semitones, -max_abs_semitones, max_abs_semitones))
    if abs(semitones) < 0.2:
        return src
    return librosa.effects.pitch_shift(src, sr=sr, n_steps=semitones)


def _spectral_knn_transfer(src, tgt, sr, style_strength=0.9):
    import numpy as np
    import librosa

    n_fft = 1024
    hop = 256
    win = "hann"

    src_stft = librosa.stft(src, n_fft=n_fft, hop_length=hop, window=win).astype(np.complex64)
    tgt_stft = librosa.stft(tgt, n_fft=n_fft, hop_length=hop, window=win).astype(np.complex64)

    src_mag = np.abs(src_stft).astype(np.float32) + 1e-7
    src_phase = np.angle(src_stft).astype(np.float32)
    tgt_mag = np.abs(tgt_stft).astype(np.float32) + 1e-7

    src_mfcc = librosa.feature.mfcc(S=librosa.power_to_db(src_mag ** 2), sr=sr, n_mfcc=20)
    tgt_mfcc = librosa.feature.mfcc(S=librosa.power_to_db(tgt_mag ** 2), sr=sr, n_mfcc=20)

    src_norm = src_mfcc / (np.linalg.norm(src_mfcc, axis=0, keepdims=True) + 1e-8)
    tgt_norm = tgt_mfcc / (np.linalg.norm(tgt_mfcc, axis=0, keepdims=True) + 1e-8)
    sim = tgt_norm.T @ src_norm
    nn_idx = np.argmax(sim, axis=0)
    tgt_aligned = tgt_mag[:, nn_idx]

    ratio = tgt_aligned / src_mag
    ratio = np.clip(ratio, 0.6, 1.9)
    ratio = np.power(ratio, style_strength)

    # Smoothing keeps intelligibility by preventing fast frame-to-frame envelope jumps.
    ratio = _smooth_2d(ratio, f_win=7, t_win=9)

    # Preserve source articulation in very low and very high bands.
    bins = ratio.shape[0]
    low_keep = int(0.06 * bins)
    high_keep = int(0.90 * bins)
    ratio[:low_keep, :] = 1.0
    ratio[high_keep:, :] = 1.0

    converted_mag = src_mag * ratio
    converted_stft = converted_mag * np.exp(1j * src_phase)
    out = librosa.istft(converted_stft, hop_length=hop, window=win, length=len(src))

    peak = np.max(np.abs(out)) + 1e-9
    if peak > 1.0:
        out = out / peak
    return out


def _target_average_envelope(target_mag, smooth_t=15):
    import numpy as np
    avg = np.mean(target_mag, axis=1, keepdims=True)
    avg = avg / (np.max(avg) + 1e-9)
    return avg


def _keepwords_librosa(src, tgt, sr, mode):
    """Source-as-carrier voice conversion using librosa STFT + cepstral envelope swap.

    The source waveform is preserved as the audio carrier so the spoken words and
    language are kept intact. Only the spectral envelope (timbre) is replaced with
    the target speaker's average envelope, and the average pitch is shifted toward
    the target's median f0. This is intentionally language-agnostic.
    """
    import numpy as np
    import librosa

    n_fft = 1024
    hop = 256
    win = "hann"

    # 1. Pitch shift source toward target average pitch (preserves intonation contour).
    src_shifted = _transfer_pitch(src, tgt, sr, max_abs_semitones=6.0)

    src_stft = librosa.stft(src_shifted, n_fft=n_fft, hop_length=hop, window=win)
    tgt_stft = librosa.stft(tgt, n_fft=n_fft, hop_length=hop, window=win)
    src_mag = np.abs(src_stft).astype(np.float32) + 1e-7
    src_phase = np.angle(src_stft).astype(np.float32)
    tgt_mag = np.abs(tgt_stft).astype(np.float32) + 1e-7

    # 2. Extract smooth spectral envelope per frame via cepstral lifter.
    log_src = np.log(src_mag)
    log_tgt = np.log(tgt_mag)
    n_keep = 30  # low-quefrency liftering -> envelope only

    def cepstral_envelope(log_mag):
        ceps = np.fft.rfft(log_mag, axis=0)
        lifter = np.zeros(ceps.shape[0], dtype=np.float32)
        lifter[:n_keep] = 1.0
        ceps = ceps * lifter[:, None]
        env = np.fft.irfft(ceps, n=log_mag.shape[0], axis=0)
        return env

    src_env = cepstral_envelope(log_src)
    tgt_env = cepstral_envelope(log_tgt)
    excitation = log_src - src_env  # source's fine spectrum (carries phonemes)

    # 3. Build average target envelope (the "voice color" we want to apply).
    tgt_env_avg = np.mean(tgt_env, axis=1, keepdims=True)

    # 4. Strength of color transfer.
    strength = 0.85 if mode == "cosy" else 0.95

    new_env = (1.0 - strength) * src_env + strength * tgt_env_avg
    new_log_mag = new_env + excitation
    new_mag = np.exp(new_log_mag)

    # 5. Re-synthesize with source phase => words and rhythm follow the source exactly.
    new_stft = new_mag * np.exp(1j * src_phase)
    out = librosa.istft(new_stft, hop_length=hop, window=win, length=len(src_shifted))

    peak = np.max(np.abs(out)) + 1e-9
    if peak > 1.0:
        out = out / peak
    return out.astype(np.float32)


def _keepwords_pyworld(src, tgt, sr, mode):
    """High-quality source-preserving VC via WORLD vocoder.

    Decomposes both signals into f0, spectral envelope (sp), and aperiodicity (ap).
    Re-synthesizes using SOURCE f0 contour and aperiodicity, but TARGET spectral
    envelope. Result keeps source words and language; voice color follows target.
    """
    import numpy as np
    import pyworld as pw

    src64 = src.astype(np.float64)
    tgt64 = tgt.astype(np.float64)

    f0_src, t_src = pw.dio(src64, sr)
    f0_src = pw.stonemask(src64, f0_src, t_src, sr)
    sp_src = pw.cheaptrick(src64, f0_src, t_src, sr)
    ap_src = pw.d4c(src64, f0_src, t_src, sr)

    f0_tgt, t_tgt = pw.dio(tgt64, sr)
    f0_tgt = pw.stonemask(tgt64, f0_tgt, t_tgt, sr)
    sp_tgt = pw.cheaptrick(tgt64, f0_tgt, t_tgt, sr)

    # 1. Build a robust average target spectral envelope (only voiced frames).
    voiced = f0_tgt > 0
    if np.any(voiced):
        sp_tgt_avg = np.mean(sp_tgt[voiced], axis=0, keepdims=True)
    else:
        sp_tgt_avg = np.mean(sp_tgt, axis=0, keepdims=True)
    sp_src_avg = np.mean(sp_src[f0_src > 0], axis=0, keepdims=True) if np.any(f0_src > 0) else np.mean(sp_src, axis=0, keepdims=True)

    # 2. Per-frame envelope morph in log domain.
    eps = 1e-12
    log_sp_src = np.log(sp_src + eps)
    log_sp_tgt_avg = np.log(sp_tgt_avg + eps)
    log_sp_src_avg = np.log(sp_src_avg + eps)
    delta = log_sp_tgt_avg - log_sp_src_avg  # average voice color shift

    strength = 0.85 if mode == "cosy" else 0.95
    log_sp_new = log_sp_src + strength * delta
    sp_new = np.exp(log_sp_new)

    # 3. Shift source f0 contour to target median (keeps prosody, changes register).
    voiced_src = f0_src > 0
    voiced_tgt = f0_tgt > 0
    if np.any(voiced_src) and np.any(voiced_tgt):
        med_src = float(np.median(f0_src[voiced_src]))
        med_tgt = float(np.median(f0_tgt[voiced_tgt]))
        if med_src > 0:
            ratio = med_tgt / med_src
            ratio = float(np.clip(ratio, 0.6, 1.7))
            f0_new = f0_src.copy()
            f0_new[voiced_src] = f0_src[voiced_src] * ratio
        else:
            f0_new = f0_src
    else:
        f0_new = f0_src

    # 4. Re-synthesize with source aperiodicity (preserves consonants, breathiness, voicing pattern).
    out = pw.synthesize(f0_new, sp_new, ap_src, sr)

    peak = np.max(np.abs(out)) + 1e-9
    if peak > 1.0:
        out = out / peak
    return out.astype(np.float32)


def _keepwords_vc(src_path, tgt_path, out_path, mode, len_sec=None):
    import numpy as np
    import librosa
    import soundfile as sf

    sr = 22050
    src, _ = librosa.load(src_path, sr=sr, mono=True)
    tgt, _ = librosa.load(tgt_path, sr=sr, mono=True)
    if len_sec is not None and len_sec > 0:
        tgt = tgt[: max(1, int(len_sec * sr))]

    if src.size == 0:
        raise RuntimeError("Source audio is empty")
    if tgt.size < sr // 4:
        raise RuntimeError("Target audio is too short. Use at least 0.25 second.")

    print(
        f"[keepwords] Source WAV: {src_path} ({len(src) / sr:.2f}s), "
        f"Target WAV: {tgt_path} ({len(tgt) / sr:.2f}s)",
        flush=True,
    )

    use_pyworld = _try_pyworld()
    if use_pyworld:
        print("[keepwords] backend: WORLD vocoder (source words and language preserved)", flush=True)
        out = _keepwords_pyworld(src, tgt, sr, mode)
    else:
        print("[keepwords] backend: librosa cepstral envelope swap (source words preserved)", flush=True)
        out = _keepwords_librosa(src, tgt, sr, mode)

    Path(out_path).resolve().parent.mkdir(parents=True, exist_ok=True)
    sf.write(str(out_path), np.asarray(out, dtype=np.float32), sr)


def _post_process(content, converted, mix=0.08):
    import numpy as np
    out = (1.0 - mix) * converted + mix * content
    peak = np.max(np.abs(out)) + 1e-9
    if peak > 1.0:
        out = out / peak
    return out


def _smooth_2d(x, f_win=5, t_win=7):
    import numpy as np
    if f_win < 2 and t_win < 2:
        return x
    y = x.copy()
    if f_win >= 2:
        kf = np.ones(f_win, dtype=np.float32) / float(f_win)
        y = np.apply_along_axis(lambda v: np.convolve(v, kf, mode="same"), axis=0, arr=y)
    if t_win >= 2:
        kt = np.ones(t_win, dtype=np.float32) / float(t_win)
        y = np.apply_along_axis(lambda v: np.convolve(v, kt, mode="same"), axis=1, arr=y)
    return y


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", required=True, help="Source wav")
    parser.add_argument("--tgt", required=True, help="Target wav")
    parser.add_argument("--out", required=True, help="Output wav")
    parser.add_argument("--mode", choices=["cosy", "knn"], default="knn")
    parser.add_argument("--backend",
                        choices=["auto", "keepwords", "speecht5", "dsp"],
                        default="keepwords",
                        help="keepwords (default) keeps source words/language and only changes timbre.")
    parser.add_argument("--len-sec", type=float, default=None, help="Optional truncate target reference")
    parser.add_argument("--source-anchor", type=float, default=0.18,
                        help="(speecht5 only) mix source waveform back after neural VC")
    args = parser.parse_args()

    if args.backend in ("auto", "keepwords"):
        try:
            _keepwords_vc(args.src, args.tgt, args.out, args.mode, args.len_sec)
            print(f"OK keepwords VC: {Path(args.out).resolve()}")
            return 0
        except Exception as exc:
            if args.backend == "keepwords":
                raise
            print(f"keepwords backend failed, falling back: {exc}", file=sys.stderr)

    if args.backend == "speecht5":
        try:
            _speecht5_vc(args.src, args.tgt, args.out, args.mode, args.len_sec, args.source_anchor)
            print(f"OK SpeechT5 VC: {Path(args.out).resolve()}")
            return 0
        except Exception as exc:
            raise

    _ensure_deps()
    import numpy as np
    import librosa
    import soundfile as sf

    src, sr = librosa.load(args.src, sr=22050, mono=True)
    tgt, _ = librosa.load(args.tgt, sr=sr, mono=True)

    if args.len_sec is not None and args.len_sec > 0:
        lim = int(args.len_sec * sr)
        tgt = tgt[: max(1, lim)]

    src_pitch_shifted = _transfer_pitch(src, tgt, sr, max_abs_semitones=8.0 if args.mode == "cosy" else 6.0)
    strength = 0.58 if args.mode == "cosy" else 0.46
    converted = _spectral_knn_transfer(src_pitch_shifted, tgt, sr, style_strength=strength)
    # Higher dry-content mix keeps words understandable.
    mix = 0.32 if args.mode == "cosy" else 0.40
    out = _post_process(src, converted, mix=mix)
    out = np.asarray(out, dtype=np.float32)

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sf.write(str(out_path), out, sr)
    print(f"OK: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

