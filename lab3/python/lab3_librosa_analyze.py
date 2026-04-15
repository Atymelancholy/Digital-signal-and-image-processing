import argparse
import json
import os
import sys


def _downsample_time(mel_db_t, max_frames: int = 260):
    # mel_db_t: shape [frames][mels]
    frames = len(mel_db_t)
    if frames <= max_frames or frames <= 1:
        return mel_db_t
    step = max(1, frames // max_frames)
    return mel_db_t[::step]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", dest="inp", required=True)
    parser.add_argument("--out", dest="out", required=True)
    args = parser.parse_args()

    payload = {
        "sample_rate": 0,
        "hop_length": 0,
        "mel_db": None,
        "centroid_hz": 0.0,
        "rolloff_hz": 0.0,
        "bandwidth_hz": 0.0,
        "zcr": 0.0,
        "centroid_series_hz": None,
        "rolloff_series_hz": None,
        "bandwidth_series_hz": None,
        "zcr_series": None,
        "error": None,
    }

    try:
        import numpy as np
        import librosa

        path = args.inp
        y, sr = librosa.load(path, sr=None, mono=True)
        # normalize to avoid NaNs on silence-heavy recordings
        if y.size == 0:
            raise RuntimeError("Empty audio")

        n_fft = 2048
        hop_length = 512
        n_mels = 64

        S = np.abs(librosa.stft(y, n_fft=n_fft, hop_length=hop_length, window="hann")) + 1e-10

        mel = librosa.feature.melspectrogram(
            y=y,
            sr=sr,
            n_fft=n_fft,
            hop_length=hop_length,
            n_mels=n_mels,
            power=2.0,
        )
        mel_db = librosa.power_to_db(mel, ref=np.max)

        # Spectral features (по пособию ЛР3)
        centroid = librosa.feature.spectral_centroid(S=S, sr=sr)
        rolloff = librosa.feature.spectral_rolloff(y=y, sr=sr, n_fft=n_fft, hop_length=hop_length, roll_percent=0.85)
        bandwidth = librosa.feature.spectral_bandwidth(S=S, sr=sr)
        zcr = librosa.feature.zero_crossing_rate(y, hop_length=hop_length)

        payload["sample_rate"] = int(sr)
        payload["hop_length"] = int(hop_length)
        payload["mel_db"] = _downsample_time(mel_db.T.tolist(), max_frames=260)  # [frames][mels]
        payload["centroid_hz"] = float(np.mean(centroid))
        payload["rolloff_hz"] = float(np.mean(rolloff))
        payload["bandwidth_hz"] = float(np.mean(bandwidth))
        payload["zcr"] = float(np.mean(zcr))
        payload["centroid_series_hz"] = centroid.squeeze(0).tolist()
        payload["rolloff_series_hz"] = rolloff.squeeze(0).tolist()
        payload["bandwidth_series_hz"] = bandwidth.squeeze(0).tolist()
        payload["zcr_series"] = zcr.squeeze(0).tolist()
    except Exception as e:
        payload["error"] = str(e)

    out_path = args.out
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False)

    # if there was an error, non-zero exit so Java can show stderr
    if payload["error"]:
        print(payload["error"], file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

