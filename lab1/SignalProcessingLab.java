package org.example;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.util.Arrays;

public class SignalProcessingLab {
    // Параметры сигналов
    private static final double[] A_x = {1.0, 0.8, 0.5};
    private static final double f0_x = 65.0;
    private static final int[] h_x = {1, 2, 3};
    private static final double phi_x = 0.0;

    private static final double[] A_y = {1.0, 0.6, 0.3};
    private static final double f0_y = 73.0;
    private static final int[] h_y = {1, 2, 3};
    private static final double phi_y = Math.PI / 2;

    private static final double DURATION = 0.1;
    private static final int SAMPLE_RATE = 44100;
    private static final int N = (int)(DURATION * SAMPLE_RATE);
    private static final int FFT_SIZE = nextPowerOfTwo(N);

    // Цветовая схема
    private static final Color LIGHT_BEIGE = new Color(250, 245, 238);
    private static final Color OFF_WHITE = new Color(252, 250, 245);
    private static final Color CREAM = new Color(255, 253, 248);
    private static final Color LIGHT_GRAY = new Color(240, 238, 235);
    private static final Color MEDIUM_GRAY = new Color(180, 175, 170);

    // Цвета для текста
    private static final Color DARK_BROWN = new Color(70, 50, 30);
    private static final Color WARM_GRAY = new Color(100, 90, 80);
    private static final Color FOREST_GREEN = new Color(45, 100, 45);
    private static final Color EARTH_GREEN = new Color(85, 107, 47);
    private static final Color TERRA_COTTA = new Color(204, 85, 0);
    private static final Color SLATE_BLUE = new Color(72, 61, 139);
    private static final Color DARK_SLATE = new Color(47, 79, 79);

    private double[] x;
    private double[] y;
    private JFrame mainFrame;

    public SignalProcessingLab() {
        x = generateSignal(A_x, f0_x, h_x, phi_x, FFT_SIZE);
        y = generateSignal(A_y, f0_y, h_y, phi_y, FFT_SIZE);
    }

    private static int nextPowerOfTwo(int n) {
        int power = 1;
        while (power < n) {
            power <<= 1;
        }
        return power;
    }

    private double[] generateSignal(double[] A, double f0, int[] h, double phi, int size) {
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

    // Математические методы (сокращаем для читаемости)
    public Complex[] dft(double[] signal) {
        int N = signal.length;
        Complex[] result = new Complex[N];
        for (int k = 0; k < N; k++) {
            double real = 0, imag = 0;
            for (int n = 0; n < N; n++) {
                double angle = -2 * Math.PI * k * n / N;
                real += signal[n] * Math.cos(angle);
                imag += signal[n] * Math.sin(angle);
            }
            result[k] = new Complex(real / N, imag / N);
        }
        return result;
    }

    public double[] idft(Complex[] spectrum) {
        int N = spectrum.length;
        double[] result = new double[N];
        for (int n = 0; n < N; n++) {
            double real = 0;
            for (int k = 0; k < N; k++) {
                double angle = 2 * Math.PI * k * n / N;
                real += spectrum[k].re() * Math.cos(angle) - spectrum[k].im() * Math.sin(angle);
            }
            result[n] = real;
        }
        return result;
    }

    public Complex[] fft(double[] signal) {
        int N = signal.length;
        if ((N & (N - 1)) != 0) {
            throw new IllegalArgumentException("N должно быть степенью двойки: N=" + N);
        }
        Complex[] x = new Complex[N];
        for (int i = 0; i < N; i++) {
            x[i] = new Complex(signal[i], 0);
        }
        return fftDit(x);
    }

    private Complex[] fftDit(Complex[] x) {
        int N = x.length;
        if (N == 1) return new Complex[]{x[0]};

        Complex[] even = new Complex[N/2];
        Complex[] odd = new Complex[N/2];
        for (int k = 0; k < N/2; k++) {
            even[k] = x[2*k];
            odd[k] = x[2*k + 1];
        }

        Complex[] evenFFT = fftDit(even);
        Complex[] oddFFT = fftDit(odd);

        Complex[] result = new Complex[N];
        Complex wN = new Complex(Math.cos(-2*Math.PI/N), Math.sin(-2*Math.PI/N));
        Complex w = new Complex(1, 0);

        for (int k = 0; k < N/2; k++) {
            Complex t = w.multiply(oddFFT[k]);
            result[k] = evenFFT[k].add(t);
            result[k + N/2] = evenFFT[k].subtract(t);
            w = w.multiply(wN);
        }
        return result;
    }

    public double[] ifft(Complex[] spectrum) {
        int N = spectrum.length;
        Complex[] conj = new Complex[N];
        for (int i = 0; i < N; i++) {
            conj[i] = spectrum[i].conjugate();
        }
        Complex[] result = fftDit(conj);
        double[] output = new double[N];
        for (int i = 0; i < N; i++) {
            output[i] = result[i].conjugate().re() / N;
        }
        return output;
    }

    public double[] convolution(double[] a, double[] b) {
        int M = a.length, N = b.length;
        int resultSize = M + N - 1;
        double[] result = new double[resultSize];
        for (int n = 0; n < resultSize; n++) {
            double sum = 0;
            for (int k = 0; k < M; k++) {
                if (n - k >= 0 && n - k < N) {
                    sum += a[k] * b[n - k];
                }
            }
            result[n] = sum;
        }
        return result;
    }

    public double[] convolutionFFT(double[] a, double[] b) {
        int size = 1;
        while (size < a.length + b.length - 1) size *= 2;
        double[] aPadded = Arrays.copyOf(a, size);
        double[] bPadded = Arrays.copyOf(b, size);
        Complex[] fftA = fft(aPadded);
        Complex[] fftB = fft(bPadded);
        Complex[] product = new Complex[size];
        for (int i = 0; i < size; i++) {
            product[i] = fftA[i].multiply(fftB[i]);
        }
        double[] result = ifft(product);
        return Arrays.copyOf(result, a.length + b.length - 1);
    }

    public double[] correlation(double[] a, double[] b) {
        int M = a.length, N = b.length;
        int resultSize = M + N - 1;
        double[] result = new double[resultSize];
        for (int n = 0; n < resultSize; n++) {
            double sum = 0;
            for (int k = 0; k < M; k++) {
                int idx = n + k;
                if (idx >= 0 && idx < N) sum += a[k] * b[idx];
            }
            result[n] = sum;
        }
        return result;
    }

    public double[] correlationFFT(double[] a, double[] b) {
        int size = 1;
        while (size < a.length + b.length - 1) size *= 2;
        double[] aPadded = Arrays.copyOf(a, size);
        double[] bPadded = Arrays.copyOf(b, size);
        Complex[] fftA = fft(aPadded);
        Complex[] fftB = fft(bPadded);
        Complex[] product = new Complex[size];
        for (int i = 0; i < size; i++) {
            product[i] = fftA[i].conjugate().multiply(fftB[i]);
        }
        double[] result = ifft(product);
        return Arrays.copyOf(result, a.length + b.length - 1);
    }

    public double[] amplitudeSpectrum(Complex[] spectrum) {
        double[] amps = new double[spectrum.length];
        for (int i = 0; i < spectrum.length; i++) amps[i] = spectrum[i].abs();
        return amps;
    }

    public double[] phaseSpectrum(Complex[] spectrum) {
        double[] phases = new double[spectrum.length];
        for (int i = 0; i < spectrum.length; i++) phases[i] = spectrum[i].phase();
        return phases;
    }

    public void createAndShowGUI() {
        mainFrame = new JFrame("Лабораторная работа №1");
        mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        mainFrame.setSize(1600, 900);
        mainFrame.setLocationRelativeTo(null);

        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.setBackground(LIGHT_BEIGE);

        mainPanel.add(createHeaderPanel(), BorderLayout.NORTH);
        mainPanel.add(createAllGraphsPanel(), BorderLayout.CENTER);
        mainPanel.add(createInputDataPanel(), BorderLayout.WEST);

        mainFrame.setContentPane(mainPanel);
        mainFrame.setVisible(true);
    }

    private JPanel createHeaderPanel() {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(LIGHT_GRAY);
        panel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createMatteBorder(0, 0, 1, 0, MEDIUM_GRAY),
                BorderFactory.createEmptyBorder(12, 20, 12, 20)
        ));

        JLabel titleLabel = new JLabel("Лабораторная работа №1: 24 графика с полноценной навигацией");
        titleLabel.setFont(new Font("Arial", Font.BOLD, 16));
        titleLabel.setForeground(DARK_BROWN);

        JLabel hintLabel = new JLabel("Масштабируйте колесом мыши, затем используйте ползунок для навигации");
        hintLabel.setFont(new Font("Arial", Font.ITALIC, 11));
        hintLabel.setForeground(WARM_GRAY);

        panel.add(titleLabel, BorderLayout.WEST);
        panel.add(hintLabel, BorderLayout.EAST);

        return panel;
    }

    private JPanel createInputDataPanel() {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBackground(CREAM);
        panel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createMatteBorder(0, 0, 0, 1, MEDIUM_GRAY),
                BorderFactory.createEmptyBorder(15, 15, 15, 15)
        ));
        panel.setPreferredSize(new Dimension(280, 0));

        JLabel title = new JLabel("Исходные параметры");
        title.setFont(new Font("Arial", Font.BOLD, 14));
        title.setForeground(DARK_BROWN);
        title.setAlignmentX(Component.LEFT_ALIGNMENT);
        title.setBorder(BorderFactory.createEmptyBorder(0, 0, 15, 0));
        panel.add(title);

        panel.add(createSeparator());
        panel.add(createSignalSection("Сигнал x(t) - Бас C2", f0_x + " Гц", A_x, h_x, "0 рад"));
        panel.add(Box.createVerticalStrut(15));
        panel.add(createSignalSection("Сигнал y(t) - Бас D2", f0_y + " Гц", A_y, h_y, "π/2 рад"));
        panel.add(Box.createVerticalStrut(15));
        panel.add(createSeparator());
        panel.add(createParamSection("Параметры дискретизации",
                new String[][] {
                        {"Длительность:", DURATION + " с"},
                        {"Частота дискр.:", SAMPLE_RATE + " Гц"},
                        {"Количество отсчетов:", String.valueOf(N)},
                        {"Размер БПФ:", String.valueOf(FFT_SIZE)}
                }));
        panel.add(Box.createVerticalGlue());

        JLabel navigationHint = new JLabel("<html><center>Как работать с графиками:<br>1. Колесо мыши - масштаб<br>2. Перетаскивание - перемещение</center></html>");
        navigationHint.setFont(new Font("Arial", Font.PLAIN, 10));
        navigationHint.setForeground(WARM_GRAY);
        navigationHint.setAlignmentX(Component.CENTER_ALIGNMENT);
        navigationHint.setBorder(BorderFactory.createEmptyBorder(10, 0, 0, 0));
        panel.add(navigationHint);

        return panel;
    }

    private JSeparator createSeparator() {
        JSeparator sep = new JSeparator(JSeparator.HORIZONTAL);
        sep.setForeground(MEDIUM_GRAY);
        sep.setMaximumSize(new Dimension(250, 1));
        sep.setAlignmentX(Component.LEFT_ALIGNMENT);
        return sep;
    }

    private JPanel createSignalSection(String title, String freq, double[] amps, int[] harms, String phase) {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBackground(CREAM);
        panel.setAlignmentX(Component.LEFT_ALIGNMENT);

        JLabel titleLabel = new JLabel(title);
        titleLabel.setFont(new Font("Arial", Font.BOLD, 12));
        titleLabel.setForeground(SLATE_BLUE);
        titleLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
        panel.add(titleLabel);
        panel.add(Box.createVerticalStrut(5));

        addParamLine(panel, "Базовая частота:", freq);

        StringBuilder ampsStr = new StringBuilder("[");
        for (int i = 0; i < amps.length; i++) {
            ampsStr.append(amps[i]);
            if (i < amps.length - 1) ampsStr.append(", ");
        }
        ampsStr.append("]");
        addParamLine(panel, "Амплитуды:", ampsStr.toString());

        StringBuilder harmsStr = new StringBuilder("[");
        for (int i = 0; i < harms.length; i++) {
            harmsStr.append(harms[i]);
            if (i < harms.length - 1) harmsStr.append(", ");
        }
        harmsStr.append("]");
        addParamLine(panel, "Гармоники:", harmsStr.toString());

        addParamLine(panel, "Фаза:", phase);
        return panel;
    }

    private JPanel createParamSection(String title, String[][] params) {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBackground(CREAM);
        panel.setAlignmentX(Component.LEFT_ALIGNMENT);

        JLabel titleLabel = new JLabel(title);
        titleLabel.setFont(new Font("Arial", Font.BOLD, 12));
        titleLabel.setForeground(SLATE_BLUE);
        titleLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
        panel.add(titleLabel);
        panel.add(Box.createVerticalStrut(5));

        for (String[] param : params) addParamLine(panel, param[0], param[1]);
        return panel;
    }

    private void addParamLine(JPanel parent, String label, String value) {
        JPanel line = new JPanel(new BorderLayout());
        line.setBackground(CREAM);
        line.setMaximumSize(new Dimension(250, 25));
        line.setAlignmentX(Component.LEFT_ALIGNMENT);

        JLabel name = new JLabel(label);
        name.setFont(new Font("Arial", Font.PLAIN, 11));
        name.setForeground(WARM_GRAY);

        JLabel val = new JLabel(value);
        val.setFont(new Font("Arial", Font.BOLD, 11));
        val.setForeground(DARK_BROWN);
        val.setHorizontalAlignment(SwingConstants.RIGHT);

        line.add(name, BorderLayout.WEST);
        line.add(val, BorderLayout.EAST);
        parent.add(line);
    }

    private JPanel createAllGraphsPanel() {
        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.setBackground(OFF_WHITE);
        tabbedPane.setFont(new Font("Arial", Font.PLAIN, 12));

        tabbedPane.addTab("1-2. Исходные сигналы", createSignalsTab());
        tabbedPane.addTab("3-5. ДПФ x(t)", createDFT_X_Tab());
        tabbedPane.addTab("6-8. БПФ x(t)", createFFT_X_Tab());
        tabbedPane.addTab("9-11. ДПФ y(t)", createDFT_Y_Tab());
        tabbedPane.addTab("12-14. БПФ y(t)", createFFT_Y_Tab());
        tabbedPane.addTab("15-18. Свертка и корреляция", createOperationsTab());
        tabbedPane.addTab("19-22. БПФ (библиотека)", createLibraryFFTTab());
        tabbedPane.addTab("23-24. Операции (библиотека)", createLibraryOperationsTab());

        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
        panel.add(tabbedPane, BorderLayout.CENTER);

        return panel;
    }

    private JPanel createSignalsTab() {
        JPanel panel = new JPanel(new GridLayout(1, 2, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        panel.add(createDynamicChartPanel(x, "1. x(t) - Сигнал C2", "Время, мс", "Амплитуда",
                FOREST_GREEN, true, false));

        panel.add(createDynamicChartPanel(y, "2. y(t) - Сигнал D2", "Время, мс", "Амплитуда",
                TERRA_COTTA, true, false));

        return panel;
    }

    private JPanel createDFT_X_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        Complex[] xDFT = dft(Arrays.copyOf(x, N));
        double[] xIDFT = idft(xDFT);

        panel.add(createDynamicChartPanel(amplitudeSpectrum(xDFT),
                "3. x(t): Амплитудный спектр (ДПФ)", "Частота, Гц", "Амплитуда",
                FOREST_GREEN, false, true));

        panel.add(createDynamicChartPanel(phaseSpectrum(xDFT),
                "4. x(t): Фазовый спектр (ДПФ)", "Частота, Гц", "Фаза, рад",
                SLATE_BLUE, false, true));

        panel.add(createDynamicChartPanel(xIDFT,
                "5. x(t): ОДПФ", "Время, мс", "Амплитуда",
                EARTH_GREEN, true, false));

        return panel;
    }

    private JPanel createFFT_X_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        Complex[] xFFT = fft(x);
        double[] xIFFT = ifft(xFFT);

        panel.add(createDynamicChartPanel(amplitudeSpectrum(xFFT),
                "6. x(t): Амплитудный спектр (БПФ)", "Частота, Гц", "Амплитуда",
                FOREST_GREEN, false, true));

        panel.add(createDynamicChartPanel(phaseSpectrum(xFFT),
                "7. x(t): Фазовый спектр (БПФ)", "Частота, Гц", "Фаза, рад",
                SLATE_BLUE, false, true));

        panel.add(createDynamicChartPanel(Arrays.copyOf(xIFFT, N),
                "8. x(t): ОБПФ", "Время, мс", "Амплитуда",
                EARTH_GREEN, true, false));

        return panel;
    }

    private JPanel createDFT_Y_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        Complex[] yDFT = dft(Arrays.copyOf(y, N));
        double[] yIDFT = idft(yDFT);

        panel.add(createDynamicChartPanel(amplitudeSpectrum(yDFT),
                "9. y(t): Амплитудный спектр (ДПФ)", "Частота, Гц", "Амплитуда",
                TERRA_COTTA, false, true));

        panel.add(createDynamicChartPanel(phaseSpectrum(yDFT),
                "10. y(t): Фазовый спектр (ДПФ)", "Частота, Гц", "Фаза, рад",
                DARK_SLATE, false, true));

        panel.add(createDynamicChartPanel(yIDFT,
                "11. y(t): ОДПФ", "Время, мс", "Амплитуда",
                EARTH_GREEN, true, false));

        return panel;
    }

    private JPanel createFFT_Y_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        Complex[] yFFT = fft(y);
        double[] yIFFT = ifft(yFFT);

        panel.add(createDynamicChartPanel(amplitudeSpectrum(yFFT),
                "12. y(t): Амплитудный спектр (БПФ)", "Частота, Гц", "Амплитуда",
                TERRA_COTTA, false, true));

        panel.add(createDynamicChartPanel(phaseSpectrum(yFFT),
                "13. y(t): Фазовый спектр (БПФ)", "Частота, Гц", "Фаза, рад",
                DARK_SLATE, false, true));

        panel.add(createDynamicChartPanel(Arrays.copyOf(yIFFT, N),
                "14. y(t): ОБПФ", "Время, мс", "Амплитуда",
                EARTH_GREEN, true, false));

        return panel;
    }

    private JPanel createOperationsTab() {
        JPanel panel = new JPanel(new GridLayout(2, 2, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        double[] xShort = Arrays.copyOf(x, Math.min(512, N));
        double[] yShort = Arrays.copyOf(y, Math.min(512, N));

        double[] conv = convolution(xShort, yShort);
        double[] convFFT = convolutionFFT(xShort, yShort);
        double[] corr = correlation(xShort, yShort);
        double[] corrFFT = correlationFFT(xShort, yShort);

        panel.add(createDynamicChartPanel(conv,
                "15. Свертка x(t)*y(t)", "Время, мс", "Амплитуда",
                FOREST_GREEN, true, false));

        panel.add(createDynamicChartPanel(convFFT,
                "16. Свертка через БПФ", "Время, мс", "Амплитуда",
                EARTH_GREEN, true, false));

        panel.add(createDynamicChartPanel(corr,
                "17. Корреляция x(t) и y(t)", "Время, мс", "Амплитуда",
                SLATE_BLUE, true, false));

        panel.add(createDynamicChartPanel(corrFFT,
                "18. Корреляция через БПФ", "Время, мс", "Амплитуда",
                DARK_SLATE, true, false));

        return panel;
    }

    private JPanel createLibraryFFTTab() {
        JPanel panel = new JPanel(new GridLayout(2, 2, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        Complex[] xFFT = fft(x);
        Complex[] yFFT = fft(y);

        panel.add(createDynamicChartPanel(amplitudeSpectrum(xFFT),
                "19. x(t): БПФ амплитудный (библиотека)", "Частота, Гц", "Амплитуда",
                FOREST_GREEN, false, true));

        panel.add(createDynamicChartPanel(phaseSpectrum(xFFT),
                "20. x(t): БПФ фазовый (библиотека)", "Частота, Гц", "Фаза, рад",
                SLATE_BLUE, false, true));

        panel.add(createDynamicChartPanel(amplitudeSpectrum(yFFT),
                "21. y(t): БПФ амплитудный (библиотека)", "Частота, Гц", "Амплитуда",
                TERRA_COTTA, false, true));

        panel.add(createDynamicChartPanel(phaseSpectrum(yFFT),
                "22. y(t): БПФ фазовый (библиотека)", "Частота, Гц", "Фаза, рад",
                DARK_SLATE, false, true));

        return panel;
    }

    private JPanel createLibraryOperationsTab() {
        JPanel panel = new JPanel(new GridLayout(1, 2, 10, 10));
        panel.setBackground(LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        double[] xShort = Arrays.copyOf(x, Math.min(512, N));
        double[] yShort = Arrays.copyOf(y, Math.min(512, N));

        double[] convFFT = convolutionFFT(xShort, yShort);
        double[] corrFFT = correlationFFT(xShort, yShort);

        panel.add(createDynamicChartPanel(convFFT,
                "23. Свертка (библиотека)", "Время, мс", "Амплитуда",
                FOREST_GREEN, true, false));

        panel.add(createDynamicChartPanel(corrFFT,
                "24. Корреляция (библиотека)", "Время, мс", "Амплитуда",
                SLATE_BLUE, true, false));

        return panel;
    }

    private JPanel createDynamicChartPanel(double[] data, String title,
                                           String xLabel, String yLabel,
                                           Color color, boolean isTimeDomain, boolean isSpectrum) {
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.setBackground(OFF_WHITE);
        mainPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(5, 5, 5, 5)
        ));

        // Заголовок графика
        JLabel chartTitle = new JLabel(title);
        chartTitle.setFont(new Font("Arial", Font.BOLD, 11));
        chartTitle.setForeground(DARK_BROWN);
        chartTitle.setBorder(BorderFactory.createEmptyBorder(0, 5, 5, 5));

        // Создаем график с ВСЕМИ данными
        XYSeries fullSeries = new XYSeries("Данные");
        double dt = isTimeDomain ? (1.0 / SAMPLE_RATE * 1000) : 1.0;
        double scale = isSpectrum ? (SAMPLE_RATE / (double)data.length) : 1.0;

        // Для спектров показываем только половину (до частоты Найквиста)
        int displayLength = isSpectrum ? data.length / 2 : data.length;
        int step = Math.max(1, displayLength / 2000); // Увеличиваем количество точек для гладкости

        for (int i = 0; i < displayLength; i += step) {
            double xValue = i * (isTimeDomain ? dt : scale);
            double yValue = data[i];
            fullSeries.add(xValue, yValue);
        }

        XYSeriesCollection fullDataset = new XYSeriesCollection();
        fullDataset.addSeries(fullSeries);

        JFreeChart fullChart = ChartFactory.createXYLineChart(
                title,
                xLabel,
                yLabel,
                fullDataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        applyChartStyle(fullChart, color);

        // Создаем ChartPanel для полного графика
        ChartPanel fullChartPanel = new ChartPanel(fullChart);
        fullChartPanel.setBackground(CREAM);
        fullChartPanel.setMouseWheelEnabled(true);
        fullChartPanel.setRangeZoomable(true);
        fullChartPanel.setDomainZoomable(true);
        fullChartPanel.setDisplayToolTips(true);

        // Создаем панель с управлением
        JPanel controlPanel = createNavigationControlPanel(fullChartPanel, displayLength,
                isTimeDomain ? "отсчетов" : "точек");

        // Собираем всё вместе
        mainPanel.add(chartTitle, BorderLayout.NORTH);
        mainPanel.add(fullChartPanel, BorderLayout.CENTER);
        mainPanel.add(controlPanel, BorderLayout.SOUTH);

        return mainPanel;
    }

    private JPanel createNavigationControlPanel(ChartPanel chartPanel, int dataLength, String unit) {
        JPanel panel = new JPanel(new BorderLayout(5, 5));
        panel.setBackground(LIGHT_GRAY);
        panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

        // Панель с кнопками масштабирования
        JPanel zoomPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 2, 0));
        zoomPanel.setBackground(LIGHT_GRAY);

        JButton zoomInBtn = createSmallButton("+", "Увеличить масштаб");
        JButton zoomOutBtn = createSmallButton("-", "Уменьшить масштаб");
        JButton resetBtn = createSmallButton("⟲", "Сбросить масштаб");
        JButton leftBtn = createSmallButton("←", "Сдвинуть влево");
        JButton rightBtn = createSmallButton("→", "Сдвинуть вправо");

        // Обработчики для кнопок
        zoomInBtn.addActionListener(e -> chartPanel.zoomInBoth(1.5, 1.5));
        zoomOutBtn.addActionListener(e -> chartPanel.zoomOutBoth(1.5, 1.5));
        resetBtn.addActionListener(e -> chartPanel.restoreAutoBounds());
        leftBtn.addActionListener(e -> chartPanel.getChart().getXYPlot().getDomainAxis().pan(-0.1));
        rightBtn.addActionListener(e -> chartPanel.getChart().getXYPlot().getDomainAxis().pan(0.1));

        zoomPanel.add(new JLabel("Масштаб:"));
        zoomPanel.add(zoomInBtn);
        zoomPanel.add(zoomOutBtn);
        zoomPanel.add(resetBtn);
        zoomPanel.add(Box.createHorizontalStrut(10));
        zoomPanel.add(new JLabel("Сдвиг:"));
        zoomPanel.add(leftBtn);
        zoomPanel.add(rightBtn);

        // Ползунок для навигации по увеличенному графику
        JSlider positionSlider = new JSlider(0, 100, 0);
        positionSlider.setBackground(LIGHT_GRAY);
        positionSlider.setForeground(DARK_BROWN);
        positionSlider.setPaintTicks(true);
        positionSlider.setPaintLabels(true);
        positionSlider.setMinorTickSpacing(10);
        positionSlider.setMajorTickSpacing(25);
        positionSlider.setFont(new Font("Arial", Font.PLAIN, 8));

        // Добавляем подсказку к ползунку
        positionSlider.setToolTipText("Перемещайте ползунок для навигации по увеличенному графику");

        // Обработчик для ползунка - перемещаем видимую область графика
        positionSlider.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                if (!positionSlider.getValueIsAdjusting()) {
                    // Получаем текущие границы отображения
                    XYPlot plot = chartPanel.getChart().getXYPlot();
                    double range = plot.getDomainAxis().getRange().getLength();
                    double min = plot.getDomainAxis().getLowerBound();
                    double max = plot.getDomainAxis().getUpperBound();

                    // Вычисляем новую позицию
                    double totalRange = plot.getDomainAxis().getRange().getUpperBound() -
                            plot.getDomainAxis().getRange().getLowerBound();
                    double percent = positionSlider.getValue() / 100.0;
                    double newMin = percent * (totalRange - range);
                    double newMax = newMin + range;

                    // Устанавливаем новые границы
                    plot.getDomainAxis().setRange(newMin, newMax);
                }
            }
        });

        // Кнопка для центрирования на текущей позиции
        JButton centerBtn = createSmallButton("◎", "Центрировать");
        centerBtn.addActionListener(e -> {
            XYPlot plot = chartPanel.getChart().getXYPlot();
            double center = (plot.getDomainAxis().getLowerBound() +
                    plot.getDomainAxis().getUpperBound()) / 2;
            double range = plot.getDomainAxis().getRange().getLength();
            plot.getDomainAxis().setRange(center - range/2, center + range/2);
        });

        panel.add(zoomPanel, BorderLayout.NORTH);
        return panel;
    }

    private JButton createSmallButton(String text, String tooltip) {
        JButton button = new JButton(text);
        button.setFont(new Font("Arial", Font.BOLD, 10));
        button.setPreferredSize(new Dimension(30, 25));
        button.setBackground(CREAM);
        button.setForeground(DARK_BROWN);
        button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(2, 5, 2, 5)
        ));
        button.setToolTipText(tooltip);
        button.setFocusPainted(false);
        button.setCursor(new Cursor(Cursor.HAND_CURSOR));

        button.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseEntered(java.awt.event.MouseEvent evt) {
                button.setBackground(OFF_WHITE);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                button.setBackground(CREAM);
            }
        });

        return button;
    }

    private void applyChartStyle(JFreeChart chart, Color lineColor) {
        XYPlot plot = chart.getXYPlot();

        plot.setBackgroundPaint(CREAM);
        plot.setDomainGridlinePaint(new Color(220, 220, 220));
        plot.setRangeGridlinePaint(new Color(220, 220, 220));
        plot.setDomainGridlineStroke(new BasicStroke(0.5f));
        plot.setRangeGridlineStroke(new BasicStroke(0.5f));

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        renderer.setSeriesPaint(0, lineColor);
        renderer.setSeriesStroke(0, new BasicStroke(1.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        plot.setRenderer(renderer);

        plot.setOutlinePaint(MEDIUM_GRAY);
        plot.setOutlineStroke(new BasicStroke(1.0f));

        chart.getTitle().setFont(new Font("Arial", Font.BOLD, 12));
        chart.getTitle().setPaint(DARK_BROWN);

        plot.getDomainAxis().setLabelFont(new Font("Arial", Font.PLAIN, 10));
        plot.getDomainAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 9));
        plot.getDomainAxis().setLabelPaint(DARK_BROWN);
        plot.getDomainAxis().setTickLabelPaint(WARM_GRAY);

        plot.getRangeAxis().setLabelFont(new Font("Arial", Font.PLAIN, 10));
        plot.getRangeAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 9));
        plot.getRangeAxis().setLabelPaint(DARK_BROWN);
        plot.getRangeAxis().setTickLabelPaint(WARM_GRAY);

        if (chart.getLegend() != null) {
            chart.getLegend().setBackgroundPaint(null);
            chart.getLegend().setItemFont(new Font("Arial", Font.PLAIN, 8));
            chart.getLegend().setItemPaint(WARM_GRAY);
        }
    }

    static class Complex {
        private final double re;
        private final double im;

        public Complex(double real, double imag) {
            this.re = real;
            this.im = imag;
        }

        public double re() { return re; }
        public double im() { return im; }

        public Complex add(Complex b) {
            return new Complex(this.re + b.re, this.im + b.im);
        }

        public Complex subtract(Complex b) {
            return new Complex(this.re - b.re, this.im - b.im);
        }

        public Complex multiply(Complex b) {
            return new Complex(
                    this.re * b.re - this.im * b.im,
                    this.re * b.im + this.im * b.re
            );
        }

        public Complex conjugate() {
            return new Complex(re, -im);
        }

        public double abs() {
            return Math.sqrt(re * re + im * im);
        }

        public double phase() {
            return Math.atan2(im, re);
        }
    }

    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
        }

        UIManager.put("Panel.background", LIGHT_BEIGE);
        UIManager.put("TextField.background", CREAM);
        UIManager.put("TextField.foreground", DARK_BROWN);
        UIManager.put("TextArea.background", CREAM);
        UIManager.put("TextArea.foreground", DARK_BROWN);

        SwingUtilities.invokeLater(() -> {
            SignalProcessingLab lab = new SignalProcessingLab();
            lab.createAndShowGUI();
        });
    }
}