package org.example;

import org.jtransforms.fft.DoubleFFT_1D;

import java.util.Arrays;

public class Lab3AudioProcessor {

    /** Средние по кадрам спектральный центроид и ширина (для п. 2–3, 4–6). */
    public record AvgSpectralFeatures(double centroidHz, double bandwidthHz) {}

    public static AvgSpectralFeatures averageSpectralFeatures(double[] signal, int sampleRate, int fftSize, int hopSize) {
        double[][] spec = spectrogramLibrary(signal, fftSize, hopSize);
        double centroidSum = 0.0;
        double bandwidthSum = 0.0;
        int frames = spec.length;
        for (double[] frame : spec) {
            double c = spectralCentroid(frame, sampleRate, fftSize);
            double b = spectralBandwidth(frame, sampleRate, fftSize, c);
            centroidSum += c;
            bandwidthSum += b;
        }
        if (frames <= 0) {
            return new AvgSpectralFeatures(0.0, 0.0);
        }
        return new AvgSpectralFeatures(centroidSum / frames, bandwidthSum / frames);
    }

    public static double[] hannWindow(int size) {
        double[] w = new double[size];
        for (int n = 0; n < size; n++) {
            w[n] = 0.5 - 0.5 * Math.cos((2.0 * Math.PI * n) / (size - 1));
        }
        return w;
    }

    public static double[][] spectrogramManual(double[] signal, int fftSize, int hopSize) {
        return computeSpectrogram(signal, fftSize, hopSize, true);
    }

    public static double[][] spectrogramLibrary(double[] signal, int fftSize, int hopSize) {
        return computeSpectrogram(signal, fftSize, hopSize, false);
    }

    private static double[][] computeSpectrogram(double[] signal, int fftSize, int hopSize, boolean manualDft) {
        int frames = 1 + Math.max(0, (signal.length - fftSize) / hopSize);
        if (signal.length < fftSize) {
            frames = 1;
        }
        int bins = fftSize / 2 + 1;
        double[][] spec = new double[frames][bins];
        double[] window = hannWindow(fftSize);
        DoubleFFT_1D fft = manualDft ? null : new DoubleFFT_1D(fftSize);

        for (int frame = 0; frame < frames; frame++) {
            int start = frame * hopSize;
            double[] frameData = new double[fftSize];
            for (int n = 0; n < fftSize; n++) {
                int idx = start + n;
                double sample = idx < signal.length ? signal[idx] : 0.0;
                frameData[n] = sample * window[n];
            }

            if (manualDft) {
                for (int k = 0; k < bins; k++) {
                    double re = 0.0;
                    double im = 0.0;
                    for (int n = 0; n < fftSize; n++) {
                        double angle = -2.0 * Math.PI * k * n / fftSize;
                        re += frameData[n] * Math.cos(angle);
                        im += frameData[n] * Math.sin(angle);
                    }
                    spec[frame][k] = Math.sqrt(re * re + im * im);
                }
            } else {
                double[] complex = new double[2 * fftSize];
                for (int n = 0; n < fftSize; n++) {
                    complex[2 * n] = frameData[n];
                    complex[2 * n + 1] = 0.0;
                }
                fft.complexForward(complex);
                for (int k = 0; k < bins; k++) {
                    double re = complex[2 * k];
                    double im = complex[2 * k + 1];
                    spec[frame][k] = Math.sqrt(re * re + im * im);
                }
            }
        }
        return spec;
    }

    public static double spectralCentroid(double[] frameSpectrum, int sampleRate, int fftSize) {
        double weightedSum = 0.0;
        double magSum = 0.0;
        for (int k = 0; k < frameSpectrum.length; k++) {
            double freq = k * (sampleRate / 2.0) / (frameSpectrum.length - 1);
            double mag = frameSpectrum[k];
            weightedSum += freq * mag;
            magSum += mag;
        }
        if (magSum < 1e-12) return 0.0;
        return weightedSum / magSum;
    }

    public static double spectralBandwidth(double[] frameSpectrum, int sampleRate, int fftSize, double centroid) {
        double weighted = 0.0;
        double magSum = 0.0;
        for (int k = 0; k < frameSpectrum.length; k++) {
            double freq = k * (sampleRate / 2.0) / (frameSpectrum.length - 1);
            double mag = frameSpectrum[k];
            double diff = freq - centroid;
            weighted += mag * diff * diff;
            magSum += mag;
        }
        if (magSum < 1e-12) return 0.0;
        return Math.sqrt(weighted / magSum);
    }

    public static double[] mixBySnrDb(double[] clean, double[] noise, double snrDb) {
        int len = clean.length;
        double[] alignedNoise = repeatOrTrim(noise, len);
        double signalPower = power(clean);
        double noisePower = power(alignedNoise);
        if (noisePower < 1e-12) return Arrays.copyOf(clean, clean.length);
        double targetNoisePower = signalPower / Math.pow(10.0, snrDb / 10.0);
        double gain = Math.sqrt(targetNoisePower / noisePower);
        double[] mixed = new double[len];
        for (int i = 0; i < len; i++) {
            mixed[i] = clean[i] + gain * alignedNoise[i];
        }
        return mixed;
    }

    public static double snrDb(double[] clean, double[] estimate) {
        int len = Math.min(clean.length, estimate.length);
        double[] noise = new double[len];
        for (int i = 0; i < len; i++) noise[i] = clean[i] - estimate[i];
        double pSignal = power(Arrays.copyOf(clean, len));
        double pNoise = power(noise);
        return 10.0 * Math.log10((pSignal + 1e-12) / (pNoise + 1e-12));
    }

    public static double siSdrDb(double[] reference, double[] estimate) {
        int len = Math.min(reference.length, estimate.length);
        double[] s = Arrays.copyOf(reference, len);
        double[] x = Arrays.copyOf(estimate, len);

        double dot = dot(x, s);
        double sEnergy = dot(s, s) + 1e-12;
        double alpha = dot / sEnergy;

        double[] sTarget = new double[len];
        double[] eNoise = new double[len];
        for (int i = 0; i < len; i++) {
            sTarget[i] = alpha * s[i];
            eNoise[i] = x[i] - sTarget[i];
        }

        double targetEnergy = dot(sTarget, sTarget);
        double noiseEnergy = dot(eNoise, eNoise) + 1e-12;
        return 10.0 * Math.log10((targetEnergy + 1e-12) / noiseEnergy);
    }

    private static double[] repeatOrTrim(double[] noise, int len) {
        if (noise.length == len) return Arrays.copyOf(noise, len);
        double[] out = new double[len];
        for (int i = 0; i < len; i++) out[i] = noise[i % noise.length];
        return out;
    }

    private static double power(double[] data) {
        double sum = 0.0;
        for (double v : data) sum += v * v;
        return data.length == 0 ? 0.0 : sum / data.length;
    }

    private static double dot(double[] a, double[] b) {
        int len = Math.min(a.length, b.length);
        double sum = 0.0;
        for (int i = 0; i < len; i++) sum += a[i] * b[i];
        return sum;
    }
}
