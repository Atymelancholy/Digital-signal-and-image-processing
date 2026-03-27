package org.example;

import java.util.Arrays;

public class SignalData {
    public static final double[] A_x = {1.0, 0.8, 0.5};
    public static final double f0_x = 65.0;
    public static final int[] h_x = {1, 2, 3};
    public static final double phi_x = 0.0;

    public static final double[] A_y = {1.0, 0.6, 0.3};
    public static final double f0_y = 73.0;
    public static final int[] h_y = {1, 2, 3};
    public static final double phi_y = Math.PI / 2;

    public static final double DURATION = 0.5;
    public static final int SAMPLE_RATE = 44100;
    public static final int N = (int) (DURATION * SAMPLE_RATE);
    public static final int FFT_SIZE = nextPowerOfTwo(N);

    public static final double[] x;
    public static final double[] y;

    static {
        double[] xN = generateSignal(A_x, f0_x, h_x, phi_x, N);
        double[] yN = generateSignal(A_y, f0_y, h_y, phi_y, N);
        x = Arrays.copyOf(xN, FFT_SIZE);
        y = Arrays.copyOf(yN, FFT_SIZE);
    }

    private static int nextPowerOfTwo(int n) {
        int power = 1;
        while (power < n) power <<= 1;
        return power;
    }

    private static double[] generateSignal(double[] A, double f0, int[] h, double phi, int size) {
        double[] signal = new double[size];
        double dt = 1.0 / SAMPLE_RATE;
        for (int i = 0; i < size; i++) {
            double t = i * dt;
            double value = 0;
            for (int j = 0; j < A.length; j++) {
                value += A[j] * Math.sin(2 * Math.PI * h[j] * f0 * t + phi);
            }
            signal[i] = value;
        }
        return signal;
    }
}