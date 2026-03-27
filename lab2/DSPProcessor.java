package org.example;

import org.jtransforms.fft.DoubleFFT_1D;
import java.util.Arrays;

public class DSPProcessor {

    public static Complex[] dft(double[] signal) {
        int N = signal.length;
        Complex[] result = new Complex[N];
        for (int k = 0; k < N; k++) {
            double real = 0, imag = 0;
            for (int n = 0; n < N; n++) {
                double angle = -2 * Math.PI * k * n / N;
                real += signal[n] * Math.cos(angle);
                imag += signal[n] * Math.sin(angle);
            }
            result[k] = new Complex(real, imag);
        }
        return result;
    }

    public static double[] idft(Complex[] spectrum) {
        int N = spectrum.length;
        double[] result = new double[N];
        for (int n = 0; n < N; n++) {
            double real = 0;
            for (int k = 0; k < N; k++) {
                double angle = 2 * Math.PI * k * n / N;
                real += spectrum[k].re() * Math.cos(angle) - spectrum[k].im() * Math.sin(angle);
            }
            result[n] = real / N;
        }
        return result;
    }

    public static Complex[] fft(double[] signal) {
        int N = signal.length;
        if ((N & (N - 1)) != 0)
            throw new IllegalArgumentException("N должно быть степенью двойки: " + N);
        Complex[] x = new Complex[N];
        for (int i = 0; i < N; i++) x[i] = new Complex(signal[i], 0);
        return fftDit(x);
    }

    private static Complex[] fftDit(Complex[] x) {
        int N = x.length;
        if (N == 1) return new Complex[] { x[0] };
        Complex[] even = new Complex[N / 2];
        Complex[] odd = new Complex[N / 2];
        for (int k = 0; k < N / 2; k++) {
            even[k] = x[2 * k];
            odd[k] = x[2 * k + 1];
        }
        Complex[] evenFFT = fftDit(even);
        Complex[] oddFFT = fftDit(odd);
        Complex[] result = new Complex[N];
        Complex wN = new Complex(Math.cos(-2 * Math.PI / N), Math.sin(-2 * Math.PI / N));
        Complex w = new Complex(1, 0);
        for (int k = 0; k < N / 2; k++) {
            Complex t = w.multiply(oddFFT[k]);
            result[k] = evenFFT[k].add(t);
            result[k + N / 2] = evenFFT[k].subtract(t);
            w = w.multiply(wN);
        }
        return result;
    }

    public static double[] ifft(Complex[] spectrum) {
        int N = spectrum.length;
        Complex[] conj = new Complex[N];
        for (int i = 0; i < N; i++) conj[i] = spectrum[i].conjugate();
        Complex[] result = fftDit(conj);
        double[] output = new double[N];
        for (int i = 0; i < N; i++) output[i] = result[i].conjugate().re() / N;
        return output;
    }

    public static Complex[] fftLibrary(double[] signal) {
        int N = signal.length;
        double[] data = new double[2 * N];
        for (int i = 0; i < N; i++) {
            data[2 * i] = signal[i];
            data[2 * i + 1] = 0.0;
        }
        DoubleFFT_1D fft = new DoubleFFT_1D(N);
        fft.complexForward(data);
        Complex[] result = new Complex[N];
        for (int i = 0; i < N; i++) result[i] = new Complex(data[2 * i], data[2 * i + 1]);
        return result;
    }

    public static double[] ifftLibrary(Complex[] spectrum) {
        int N = spectrum.length;
        double[] data = new double[2 * N];
        for (int i = 0; i < N; i++) {
            data[2 * i] = spectrum[i].re();
            data[2 * i + 1] = spectrum[i].im();
        }
        DoubleFFT_1D fft = new DoubleFFT_1D(N);
        fft.complexInverse(data, true);
        double[] result = new double[N];
        for (int i = 0; i < N; i++) result[i] = data[2 * i];
        return result;
    }

    public static double[] convolution(double[] a, double[] b) {
        int M = a.length, N = b.length;
        int resultSize = M + N - 1;
        double[] result = new double[resultSize];
        for (int n = 0; n < resultSize; n++) {
            double sum = 0;
            for (int k = 0; k < M; k++) {
                if (n - k >= 0 && n - k < N) sum += a[k] * b[n - k];
            }
            result[n] = sum;
        }
        return result;
    }

    public static double[] correlation(double[] a, double[] b) {
        int M = a.length, N = b.length;
        int size = M + N - 1;
        double[] r = new double[size];
        for (int n = 0; n < size; n++) {
            int lag = n - (N - 1);
            double sum = 0.0;
            for (int k = 0; k < M; k++) {
                int j = k + lag;
                if (j >= 0 && j < N) sum += a[k] * b[j];
            }
            r[n] = sum;
        }
        return r;
    }

    public static double[] convolutionFFT(double[] a, double[] b) {
        int size = 1;
        while (size < a.length + b.length - 1) size *= 2;
        double[] aPadded = Arrays.copyOf(a, size);
        double[] bPadded = Arrays.copyOf(b, size);
        Complex[] fftA = fft(aPadded);
        Complex[] fftB = fft(bPadded);
        Complex[] product = new Complex[size];
        for (int i = 0; i < size; i++) product[i] = fftA[i].multiply(fftB[i]);
        double[] result = ifft(product);
        return Arrays.copyOf(result, a.length + b.length - 1);
    }

    public static double[] correlationFFT(double[] a, double[] b) {
        int M = a.length, N = b.length;
        int full = M + N - 1;
        int size = 1;
        while (size < full) size <<= 1;
        double[] aPadded = Arrays.copyOf(a, size);
        double[] bPadded = Arrays.copyOf(b, size);
        Complex[] fftA = fft(aPadded);
        Complex[] fftB = fft(bPadded);
        Complex[] product = new Complex[size];
        for (int i = 0; i < size; i++) product[i] = fftA[i].conjugate().multiply(fftB[i]);
        double[] circ = ifft(product);
        double[] lin = Arrays.copyOf(circ, full);
        double[] r = new double[full];
        int neg = N - 1;
        System.arraycopy(lin, full - neg, r, 0, neg);
        System.arraycopy(lin, 0, r, neg, M);
        return r;
    }

    public static double[] convolutionFFTLibrary(double[] a, double[] b) {
        int full = a.length + b.length - 1;
        int size = 1;
        while (size < full) size <<= 1;
        double[] aPad = Arrays.copyOf(a, size);
        double[] bPad = Arrays.copyOf(b, size);
        Complex[] A = fftLibrary(aPad);
        Complex[] B = fftLibrary(bPad);
        Complex[] Z = new Complex[size];
        for (int i = 0; i < size; i++) Z[i] = A[i].multiply(B[i]);
        return Arrays.copyOf(ifftLibrary(Z), full);
    }

    public static double[] correlationFFTLibrary(double[] a, double[] b) {
        int M = a.length, N = b.length;
        int full = M + N - 1;
        int size = 1;
        while (size < full) size <<= 1;
        double[] aPad = Arrays.copyOf(a, size);
        double[] bPad = Arrays.copyOf(b, size);
        Complex[] A = fftLibrary(aPad);
        Complex[] B = fftLibrary(bPad);
        Complex[] Z = new Complex[size];
        for (int i = 0; i < size; i++) Z[i] = A[i].conjugate().multiply(B[i]);
        double[] circ = ifftLibrary(Z);
        double[] lin = Arrays.copyOf(circ, full);
        double[] r = new double[full];
        int neg = N - 1;
        System.arraycopy(lin, full - neg, r, 0, neg);
        System.arraycopy(lin, 0, r, neg, M);
        return r;
    }

    public static double[] amplitudeSpectrum(Complex[] spectrum) {
        double[] amps = new double[spectrum.length];
        for (int i = 0; i < spectrum.length; i++) amps[i] = spectrum[i].abs();
        return amps;
    }

    public static double[] phaseForPlot(Complex[] spectrum, double relThreshold) {
        int n = spectrum.length;
        double max = 0.0;
        double[] amp = new double[n];
        for (int i = 0; i < n; i++) {
            amp[i] = spectrum[i].abs();
            if (amp[i] > max) max = amp[i];
        }
        double th = max * relThreshold;
        double[] phase = new double[n];
        for (int i = 0; i < n; i++) {
            phase[i] = (amp[i] < th) ? Double.NaN : spectrum[i].phase();
        }
        return phase;
    }

    public static double[] applyMovingAverageFilter(double[] input, int M) {
        int len = input.length;
        double[] output = new double[len];
        double invM = 1.0 / M;

        for (int n = 0; n < len; n++) {
            double sum = 0.0;
            for (int k = 0; k < M; k++) {
                int idx = n - k;
                if (idx >= 0) sum += input[idx];
            }
            output[n] = sum * invM;
        }
        return output;
    }

    public static double[] designHighpassFIR(double fc, int M, int sampleRate) {
        if (M % 2 == 0) throw new IllegalArgumentException("M должно быть нечётным");
        int N = M - 1;
        int mid = N / 2;

        double wc = 2.0 * Math.PI * fc / sampleRate;

        double[] hlp = new double[M];
        for (int n = 0; n < M; n++) {
            int m = n - mid;

            double hd = (m == 0)
                    ? (wc / Math.PI)
                    : (Math.sin(wc * m) / (Math.PI * m));

            double w = 0.54 - 0.46 * Math.cos(2.0 * Math.PI * n / N);
            hlp[n] = hd * w;
        }

        double sum = 0.0;
        for (double v : hlp) sum += v;
        for (int i = 0; i < M; i++) hlp[i] /= sum;

        double[] hhp = new double[M];
        for (int i = 0; i < M; i++) hhp[i] = -hlp[i];
        hhp[mid] += 1.0;

        return hhp;
    }

    public static double[] applyFIRFilter(double[] input, double[] h) {
        double[] full = convolution(input, h);
        int delay = (h.length - 1) / 2;

        int start = delay;
        int end = Math.min(start + input.length, full.length);

        double[] y = Arrays.copyOfRange(full, start, end);
        if (y.length < input.length) y = Arrays.copyOf(y, input.length);
        return y;
    }

    public static double[] designNotchIIR(double f0, double BW, int sampleRate) {
        double theta = 2.0 * Math.PI * f0 / sampleRate;
        double bwNorm = BW / sampleRate;
        double R = 1.0 - 3.0 * bwNorm;

        R = Math.max(1e-6, Math.min(R, 0.999999));

        double K = (1.0 - 2.0 * R * Math.cos(theta) + R * R)
                / (2.0 - 2.0 * Math.cos(theta));

        double[] coeffs = new double[5];
        coeffs[0] = K;
        coeffs[1] = -2.0 * K * Math.cos(theta);
        coeffs[2] = K;

        coeffs[3] = -2.0 * R * Math.cos(theta);
        coeffs[4] = R * R;

        return coeffs;
    }

    public static double[] applyIIRFilter(double[] input, double[] coeffs) {
        int len = input.length;
        double[] output = new double[len];
        double a0 = coeffs[0];
        double a1 = coeffs[1];
        double a2 = coeffs[2];
        double b1 = coeffs[3];
        double b2 = coeffs[4];
        for (int n = 0; n < len; n++) {
            double x0 = input[n];
            double x1 = (n > 0) ? input[n - 1] : 0;
            double x2 = (n > 1) ? input[n - 2] : 0;
            double y1 = (n > 0) ? output[n - 1] : 0;
            double y2 = (n > 1) ? output[n - 2] : 0;
            output[n] = a0 * x0 + a1 * x1 + a2 * x2 - b1 * y1 - b2 * y2;
        }
        return output;
    }

    public static double[] calculateFrequencyResponse(double[] h, boolean isIIR, double[] aCoeffs, int sampleRate, int fftSize) {
        double[] response = new double[fftSize / 2];
        if (!isIIR) {
            double[] hPadded = Arrays.copyOf(h, fftSize);
            Complex[] H = fftLibrary(hPadded);
            for (int i = 0; i < fftSize / 2; i++) {
                response[i] = H[i].abs();
            }
        } else {
            double dw = Math.PI / (fftSize / 2);
            for (int idx = 0; idx < fftSize / 2; idx++) {
                double w = idx * dw;
                Complex z = new Complex(Math.cos(w), -Math.sin(w));
                Complex z2 = new Complex(Math.cos(2 * w), -Math.sin(2 * w));
                Complex numerator = new Complex(aCoeffs[0], 0)
                        .add(new Complex(aCoeffs[1], 0).multiply(z))
                        .add(new Complex(aCoeffs[2], 0).multiply(z2));
                Complex denominator = new Complex(1, 0)
                        .add(new Complex(aCoeffs[3], 0).multiply(z))
                        .add(new Complex(aCoeffs[4], 0).multiply(z2));
                double denomAbs = denominator.abs();
                response[idx] = numerator.abs() / (denomAbs > 1e-10 ? denomAbs : 1e-10);
            }
        }
        return response;
    }

    public static double[] addLowFrequencyNoise(double[] signal, double noiseFreq, int sampleRate) {
        double[] noisy = signal.clone();
        double dt = 1.0 / sampleRate;
        for (int i = 0; i < noisy.length; i++) {
            noisy[i] += 0.5 * Math.sin(2 * Math.PI * noiseFreq * i * dt);
        }
        return noisy;
    }

    public static double[] addSingleFrequencyNoise(double[] signal, double noiseFreq, int sampleRate) {
        double[] noisy = signal.clone();
        double dt = 1.0 / sampleRate;
        for (int i = 0; i < noisy.length; i++) {
            noisy[i] += 0.7 * Math.sin(2 * Math.PI * noiseFreq * i * dt);
        }
        return noisy;
    }

    public static double[] addHighFrequencyNoise(double[] signal, double noiseFreq, double amplitude, int sampleRate) {
        double[] noisy = signal.clone();
        double dt = 1.0 / sampleRate;
        for (int i = 0; i < noisy.length; i++) {
            noisy[i] += amplitude * Math.sin(2 * Math.PI * noiseFreq * i * dt);
        }
        return noisy;
    }
}