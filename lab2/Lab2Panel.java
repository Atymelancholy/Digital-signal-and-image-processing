package org.example;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.sound.sampled.*;
import javax.swing.*;
import java.awt.*;
import java.util.Arrays;

public class Lab2Panel extends JPanel {

    // Параметры фильтров согласно варианту
    private static final int MOVING_AVG_ORDER = 31;
    private static final double FIR_HIGHPASS_CUTOFF = 150.0;
    private static final int FIR_TAPS = 101;
    private static final double IIR_NOTCH_CENTER_FREQ = 130.0;
    private static final double IIR_NOTCH_BANDWIDTH = 8.0;

    public Lab2Panel() {
        setLayout(new BorderLayout());
        setBackground(ChartUtils.LIGHT_BEIGE);

        JTabbedPane filterTabs = new JTabbedPane();
        filterTabs.setBackground(ChartUtils.OFF_WHITE);
        filterTabs.setFont(new Font("Arial", Font.PLAIN, 12));

        filterTabs.addTab("Однородный фильтр (M=31)", createMovingAvgFilterPanel());
        filterTabs.addTab("КИХ-фильтр (ВЧ, Хэмминг, fc=150 Гц)", createFIRFilterPanel());
        filterTabs.addTab("БИХ-фильтр (Режекторный, f0=130 Гц)", createIIRFilterPanel());

        add(filterTabs, BorderLayout.CENTER);
    }

    // ================== Однородный фильтр ==================
    private JPanel createMovingAvgFilterPanel() {
        JPanel panel = new JPanel(new GridLayout(2, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        double[] filteredX = DSPProcessor.applyMovingAverageFilter(SignalData.x, MOVING_AVG_ORDER);

        int displayLen = Math.min(2000, SignalData.x.length);
        double[] xShort = Arrays.copyOf(SignalData.x, displayLen);
        double[] filteredXShort = Arrays.copyOf(filteredX, displayLen);

        // Графики с кнопками воспроизведения
        panel.add(createPlayableChartPanel(xShort, "Исходный сигнал x(t)", true, SignalData.x));
        panel.add(createPlayableChartPanel(filteredXShort, "После однородного фильтра", true, filteredX));

        // АЧХ однородного фильтра
        double[] h = new double[MOVING_AVG_ORDER];
        Arrays.fill(h, 1.0 / MOVING_AVG_ORDER);
        double[] freqResponse = DSPProcessor.calculateFrequencyResponse(h, false, null,
                SignalData.SAMPLE_RATE, SignalData.FFT_SIZE);

        XYSeries series = new XYSeries("АЧХ");
        for (int i = 1; i < freqResponse.length; i += 10) {
            double freq = i * (SignalData.SAMPLE_RATE / 2.0) / (freqResponse.length - 1);
            double magnitudeDB = 20 * Math.log10(freqResponse[i] + 1e-10);
            series.add(freq, magnitudeDB);
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "АЧХ однородного фильтра", "Частота, Гц", "Амплитуда, дБ",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartUtils.applyChartStyle(chart, ChartUtils.SLATE_BLUE);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setMouseWheelEnabled(true);
        panel.add(chartPanel);

        // Информационная панель
        JPanel infoPanel = new JPanel();
        infoPanel.setBackground(ChartUtils.CREAM);
        infoPanel.setBorder(BorderFactory.createTitledBorder("Параметры"));
        infoPanel.setLayout(new BoxLayout(infoPanel, BoxLayout.Y_AXIS));
        infoPanel.add(new JLabel("Тип: Однородный нерекурсивный"));
        infoPanel.add(new JLabel("M = " + MOVING_AVG_ORDER));
        panel.add(infoPanel);

        return panel;
    }

    // ================== КИХ ВЧ фильтр ==================
    private JPanel createFIRFilterPanel() {
        JPanel panel = new JPanel(new GridLayout(2, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        double[] firCoeffs = DSPProcessor.designHighpassFIR(FIR_HIGHPASS_CUTOFF, FIR_TAPS, SignalData.SAMPLE_RATE);
        System.out.println("sum(h) = " + Arrays.stream(firCoeffs).sum());
        double[] noisyX = DSPProcessor.addLowFrequencyNoise(SignalData.x, 30, SignalData.SAMPLE_RATE);
        double[] filteredX = DSPProcessor.applyFIRFilter(noisyX, firCoeffs);

        int displayLen = Math.min(2000, SignalData.x.length);
        double[] noisyXShort = Arrays.copyOf(noisyX, displayLen);
        double[] filteredXShort = Arrays.copyOf(filteredX, displayLen);

        panel.add(createPlayableChartPanel(noisyXShort, "Сигнал с НЧ-помехой (30 Гц)", true, noisyX));
        panel.add(createPlayableChartPanel(filteredXShort, "После КИХ ВЧ фильтра", true, filteredX));

        // АЧХ КИХ-фильтра
        double[] freqResponse = DSPProcessor.calculateFrequencyResponse(firCoeffs, false, null,
                SignalData.SAMPLE_RATE, SignalData.FFT_SIZE);
        XYSeries series = new XYSeries("АЧХ");
        for (int i = 1; i < freqResponse.length; i += 1) {
            double freq = i * (SignalData.SAMPLE_RATE / 2.0) / (freqResponse.length - 1);
            if (freq > 1000) break;
            series.add(freq, 20 * Math.log10(freqResponse[i] + 1e-12));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "АЧХ КИХ ВЧ фильтра", "Частота, Гц", "Амплитуда, дБ",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartUtils.applyChartStyle(chart, ChartUtils.SLATE_BLUE);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setMouseWheelEnabled(true);
        panel.add(chartPanel);

        // Информационная панель
        JPanel infoPanel = new JPanel();
        infoPanel.setBackground(ChartUtils.CREAM);
        infoPanel.setBorder(BorderFactory.createTitledBorder("Параметры"));
        infoPanel.setLayout(new BoxLayout(infoPanel, BoxLayout.Y_AXIS));
        infoPanel.add(new JLabel("Тип: КИХ ВЧ"));
        infoPanel.add(new JLabel("Окно: Хэмминга"));
        infoPanel.add(new JLabel("fc = " + FIR_HIGHPASS_CUTOFF + " Гц"));
        infoPanel.add(new JLabel("Порядок N = " + FIR_TAPS));
        panel.add(infoPanel);

        return panel;
    }

    // ================== БИХ режекторный фильтр ==================
    private JPanel createIIRFilterPanel() {
        JPanel panel = new JPanel(new GridLayout(2, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        double[] iirCoeffs = DSPProcessor.designNotchIIR(IIR_NOTCH_CENTER_FREQ, IIR_NOTCH_BANDWIDTH, SignalData.SAMPLE_RATE);
        double[] noisyX = DSPProcessor.addSingleFrequencyNoise(SignalData.x, IIR_NOTCH_CENTER_FREQ, SignalData.SAMPLE_RATE);
        double[] filteredX = DSPProcessor.applyIIRFilter(noisyX, iirCoeffs);

        int displayLen = Math.min(2000, SignalData.x.length);
        double[] noisyXShort = Arrays.copyOf(noisyX, displayLen);
        double[] filteredXShort = Arrays.copyOf(filteredX, displayLen);

        panel.add(createPlayableChartPanel(noisyXShort, "Сигнал с помехой " + (int) IIR_NOTCH_CENTER_FREQ + " Гц", true, noisyX));
        panel.add(createPlayableChartPanel(filteredXShort, "После БИХ режекторного", true, filteredX));

        // АЧХ БИХ-фильтра
        int RESP_FFT = 1 << 18; // 262144
        double[] freqResponse = DSPProcessor.calculateFrequencyResponse(null, true, iirCoeffs,
                SignalData.SAMPLE_RATE, RESP_FFT);
        XYSeries series = new XYSeries("АЧХ");
        for (int i = 1; i < freqResponse.length; i += 1) {
            double freq = i * (SignalData.SAMPLE_RATE / 2.0) / (freqResponse.length - 1);
            series.add(freq, 20 * Math.log10(freqResponse[i] + 1e-10));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "АЧХ БИХ режекторного фильтра", "Частота, Гц", "Амплитуда, дБ",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartUtils.applyChartStyle(chart, ChartUtils.SLATE_BLUE);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setMouseWheelEnabled(true);
        panel.add(chartPanel);

        // Информационная панель
        JPanel infoPanel = new JPanel();
        infoPanel.setBackground(ChartUtils.CREAM);
        infoPanel.setBorder(BorderFactory.createTitledBorder("Параметры"));
        infoPanel.setLayout(new BoxLayout(infoPanel, BoxLayout.Y_AXIS));
        infoPanel.add(new JLabel("Тип: БИХ режекторный"));
        infoPanel.add(new JLabel("f0 = " + IIR_NOTCH_CENTER_FREQ + " Гц"));
        infoPanel.add(new JLabel("BW = " + IIR_NOTCH_BANDWIDTH + " Гц"));
        infoPanel.add(new JLabel("a0=" + String.format("%.3f", iirCoeffs[0]) + ", a1=" + String.format("%.3f", iirCoeffs[1])));
        infoPanel.add(new JLabel("a2=" + String.format("%.3f", iirCoeffs[2]) + ", b1=" + String.format("%.3f", iirCoeffs[3])));
        infoPanel.add(new JLabel("b2=" + String.format("%.3f", iirCoeffs[4])));
        panel.add(infoPanel);

        return panel;
    }

    // ================== Вспомогательный метод для создания графика с кнопкой воспроизведения ==================
    private JPanel createPlayableChartPanel(double[] displayData, String title, boolean isTimeDomain, double[] fullSignal) {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(ChartUtils.OFF_WHITE);
        panel.setBorder(BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY, 1));

        // График (метод createDynamicChartPanel уже добавляет заголовок)
        JPanel chartPanel = ChartUtils.createDynamicChartPanel(displayData, title,
                "Время, мс", "Амплитуда", ChartUtils.FOREST_GREEN, isTimeDomain, false);

        // Кнопка воспроизведения
        JButton playButton = new JButton("▶ Воспроизвести");
        playButton.setFont(new Font("Arial", Font.PLAIN, 10));
        playButton.setBackground(ChartUtils.CREAM);
        playButton.setForeground(ChartUtils.DARK_BROWN);
        playButton.addActionListener(e -> playSignal(fullSignal));

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        buttonPanel.setBackground(ChartUtils.OFF_WHITE);
        buttonPanel.add(playButton);

        panel.add(chartPanel, BorderLayout.CENTER);
        panel.add(buttonPanel, BorderLayout.SOUTH);
        return panel;
    }

    // ================== Воспроизведение сигнала ==================
    private void playSignal(double[] signal) {
        final int sampleRate = SignalData.SAMPLE_RATE;
        final int bits = 16;
        final int channels = 1;
        final boolean signed = true;
        final boolean bigEndian = false;

        // Нормализация сигнала в диапазон [-1, 1] для 16-bit PCM
        double max = 0;
        for (double v : signal) {
            if (Math.abs(v) > max) max = Math.abs(v);
        }
        if (max == 0) return;
        final double amplification = 0.9; // избегаем клиппинга
        byte[] audioData = new byte[signal.length * 2];
        for (int i = 0; i < signal.length; i++) {
            short sample = (short) ((signal[i] / max) * amplification * Short.MAX_VALUE);
            audioData[2 * i] = (byte) (sample & 0xff);
            audioData[2 * i + 1] = (byte) ((sample >> 8) & 0xff);
        }

        try {
            AudioFormat format = new AudioFormat(sampleRate, bits, channels, signed, bigEndian);
            DataLine.Info info = new DataLine.Info(SourceDataLine.class, format);
            SourceDataLine line = (SourceDataLine) AudioSystem.getLine(info);
            line.open(format);
            line.start();
            line.write(audioData, 0, audioData.length);
            line.drain();
            line.close();
        } catch (LineUnavailableException ex) {
            JOptionPane.showMessageDialog(this, "Не удалось воспроизвести звук: " + ex.getMessage(),
                    "Ошибка", JOptionPane.ERROR_MESSAGE);
        }
    }
}