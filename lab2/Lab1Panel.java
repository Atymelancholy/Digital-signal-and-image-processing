// файл: org/example/Lab1Panel.java
package org.example;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;

public class Lab1Panel extends JPanel {
    public Lab1Panel() {
        setLayout(new BorderLayout());
        setBackground(ChartUtils.LIGHT_BEIGE);

        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.setBackground(ChartUtils.OFF_WHITE);
        tabbedPane.setFont(new Font("Arial", Font.PLAIN, 12));

        tabbedPane.addTab("1-2. Исходные сигналы", createSignalsTab());
        tabbedPane.addTab("3-5. ДПФ x(t)", createDFT_X_Tab());
        tabbedPane.addTab("6-8. БПФ x(t)", createFFT_X_Tab());
        tabbedPane.addTab("9-11. ДПФ y(t)", createDFT_Y_Tab());
        tabbedPane.addTab("12-14. БПФ y(t)", createFFT_Y_Tab());
        tabbedPane.addTab("15-18. Свертка и корреляция", createOperationsTab());
        tabbedPane.addTab("19-22. БПФ (библиотека)", createLibraryFFTTab());
        tabbedPane.addTab("23-24. Операции (библиотека)", createLibraryOperationsTab());
        tabbedPane.addTab("25. Все графики (обзор)", createAllInOneTab());

        add(tabbedPane, BorderLayout.CENTER);
    }

    private JPanel createSignalsTab() {
        JPanel panel = new JPanel(new GridLayout(1, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        int displayLen = Math.min(2000, SignalData.N);
        double[] xShort = Arrays.copyOf(SignalData.x, displayLen);
        double[] yShort = Arrays.copyOf(SignalData.y, displayLen);

        panel.add(ChartUtils.createDynamicChartPanel(xShort, "1. x(t) - Сигнал C2",
                "Время, мс", "Амплитуда", ChartUtils.FOREST_GREEN, true, false));
        panel.add(ChartUtils.createDynamicChartPanel(yShort, "2. y(t) - Сигнал D2",
                "Время, мс", "Амплитуда", ChartUtils.TERRA_COTTA, true, false));
        return panel;
    }

    private JPanel createDFT_X_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

//        Complex[] xDFT = DSPProcessor.dft(SignalData.x);
//        double[] xIDFT = DSPProcessor.idft(xDFT);
//
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.amplitudeSpectrum(xDFT),
//                "3. x(t): Амплитудный спектр (ДПФ)", "Частота, Гц", "Амплитуда",
//                ChartUtils.FOREST_GREEN, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.phaseForPlot(xDFT, 0.1),
//                "4. x(t): Фазовый спектр (ДПФ)", "Частота, Гц", "Фаза, рад",
//                ChartUtils.SLATE_BLUE, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(Arrays.copyOf(xIDFT, SignalData.N),
//                "5. x(t): ОДПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
        return panel;
    }

    private JPanel createFFT_X_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

//        Complex[] xFFT = DSPProcessor.fft(SignalData.x);
//        double[] xIFFT = DSPProcessor.ifft(xFFT);
//
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.amplitudeSpectrum(xFFT),
//                "6. x(t): Амплитудный спектр (БПФ)", "Частота, Гц", "Амплитуда",
//                ChartUtils.FOREST_GREEN, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.phaseForPlot(xFFT, 0.1),
//                "7. x(t): Фазовый спектр (БПФ)", "Частота, Гц", "Фаза, рад",
//                ChartUtils.SLATE_BLUE, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(Arrays.copyOf(xIFFT, SignalData.N),
//                "8. x(t): ОБПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
        return panel;
    }

    private JPanel createDFT_Y_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

//        Complex[] yDFT = DSPProcessor.dft(SignalData.y);
//        double[] yIDFT = DSPProcessor.idft(yDFT);
//
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.amplitudeSpectrum(yDFT),
//                "9. y(t): Амплитудный спектр (ДПФ)", "Частота, Гц", "Амплитуда",
//                ChartUtils.TERRA_COTTA, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.phaseForPlot(yDFT, 0.1),
//                "10. y(t): Фазовый спектр (ДПФ)", "Частота, Гц", "Фаза, рад",
//                ChartUtils.DARK_SLATE, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(Arrays.copyOf(yIDFT, SignalData.N),
//                "11. y(t): ОДПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
        return panel;
    }

    private JPanel createFFT_Y_Tab() {
        JPanel panel = new JPanel(new GridLayout(1, 3, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

//        Complex[] yFFT = DSPProcessor.fft(SignalData.y);
//        double[] yIFFT = DSPProcessor.ifft(yFFT);
//
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.amplitudeSpectrum(yFFT),
//                "12. y(t): Амплитудный спектр (БПФ)", "Частота, Гц", "Амплитуда",
//                ChartUtils.TERRA_COTTA, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.phaseForPlot(yFFT, 0.1),
//                "13. y(t): Фазовый спектр (БПФ)", "Частота, Гц", "Фаза, рад",
//                ChartUtils.DARK_SLATE, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(Arrays.copyOf(yIFFT, SignalData.N),
//                "14. y(t): ОБПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
        return panel;
    }

    private JPanel createOperationsTab() {
        JPanel panel = new JPanel(new GridLayout(2, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

//        double[] xShort = Arrays.copyOf(SignalData.x, Math.min(512, SignalData.N));
//        double[] yShort = Arrays.copyOf(SignalData.y, Math.min(512, SignalData.N));
//
//        double[] conv = DSPProcessor.convolution(xShort, yShort);
//        double[] convFFT = DSPProcessor.convolutionFFT(xShort, yShort);
//        double[] corr = DSPProcessor.correlation(xShort, yShort);
//        double[] corrFFT = DSPProcessor.correlationFFT(xShort, yShort);
//
//        panel.add(ChartUtils.createDynamicChartPanel(conv, "15. Свертка x(t)*y(t)",
//                "Время, мс", "Амплитуда", ChartUtils.FOREST_GREEN, true, false));
//        panel.add(ChartUtils.createDynamicChartPanel(convFFT, "16. Свертка через БПФ",
//                "Время, мс", "Амплитуда", ChartUtils.EARTH_GREEN, true, false));
//        panel.add(ChartUtils.createDynamicChartPanel(corr, "17. Корреляция x(t) и y(t)",
//                "Время, мс", "Амплитуда", ChartUtils.SLATE_BLUE, true, false));
//        panel.add(ChartUtils.createDynamicChartPanel(corrFFT, "18. Корреляция через БПФ",
//                "Время, мс", "Амплитуда", ChartUtils.DARK_SLATE, true, false));
        return panel;
    }

    private JPanel createLibraryFFTTab() {
        JPanel panel = new JPanel(new GridLayout(2, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

//        Complex[] xFFT = DSPProcessor.fftLibrary(SignalData.x);
//        Complex[] yFFT = DSPProcessor.fftLibrary(SignalData.y);
//
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.amplitudeSpectrum(xFFT),
//                "19. x(t): БПФ амплитудный (библиотека)", "Частота, Гц", "Амплитуда",
//                ChartUtils.FOREST_GREEN, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.phaseForPlot(xFFT, 0.1),
//                "20. x(t): БПФ фазовый (библиотека)", "Частота, Гц", "Фаза, рад",
//                ChartUtils.SLATE_BLUE, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.amplitudeSpectrum(yFFT),
//                "21. y(t): БПФ амплитудный (библиотека)", "Частота, Гц", "Амплитуда",
//                ChartUtils.TERRA_COTTA, false, true));
//        panel.add(ChartUtils.createDynamicChartPanel(DSPProcessor.phaseForPlot(yFFT, 0.1),
//                "22. y(t): БПФ фазовый (библиотека)", "Частота, Гц", "Фаза, рад",
//                ChartUtils.DARK_SLATE, false, true));
        return panel;
    }

    private JPanel createLibraryOperationsTab() {
        JPanel panel = new JPanel(new GridLayout(1, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

//        double[] xShort = Arrays.copyOf(SignalData.x, Math.min(512, SignalData.N));
//        double[] yShort = Arrays.copyOf(SignalData.y, Math.min(512, SignalData.N));
//
//        double[] convFFT = DSPProcessor.convolutionFFTLibrary(xShort, yShort);
//        double[] corrFFT = DSPProcessor.correlationFFTLibrary(xShort, yShort);
//
//        panel.add(ChartUtils.createDynamicChartPanel(convFFT, "23. Свертка (библиотека)",
//                "Время, мс", "Амплитуда", ChartUtils.FOREST_GREEN, true, false));
//        panel.add(ChartUtils.createDynamicChartPanel(corrFFT, "24. Корреляция (библиотека)",
//                "Время, мс", "Амплитуда", ChartUtils.SLATE_BLUE, true, false));
        return panel;
    }

    private JPanel createAllInOneTab() {
        JPanel mainPanel = new JPanel(new BorderLayout(5, 5));
        mainPanel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel graphsPanel = new JPanel(new GridLayout(6, 4, 5, 5));
        graphsPanel.setBackground(ChartUtils.LIGHT_BEIGE);
        graphsPanel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

//        Complex[] xDFT = DSPProcessor.dft(Arrays.copyOf(SignalData.x, SignalData.N));
//        double[] xIDFT = DSPProcessor.idft(xDFT);
//        Complex[] xFFT = DSPProcessor.fft(SignalData.x);
//        double[] xIFFT = DSPProcessor.ifft(xFFT);
//        Complex[] yDFT = DSPProcessor.dft(Arrays.copyOf(SignalData.y, SignalData.N));
//        double[] yIDFT = DSPProcessor.idft(yDFT);
//        Complex[] yFFT = DSPProcessor.fft(SignalData.y);
//        double[] yIFFT = DSPProcessor.ifft(yFFT);
//
//        Complex[] xFFTLib = DSPProcessor.fftLibrary(SignalData.x);
//        Complex[] yFFTLib = DSPProcessor.fftLibrary(SignalData.y);
//
//        double[] xShort = Arrays.copyOf(SignalData.x, Math.min(512, SignalData.N));
//        double[] yShort = Arrays.copyOf(SignalData.y, Math.min(512, SignalData.N));
//        double[] conv = DSPProcessor.convolution(xShort, yShort);
//        double[] convFFT = DSPProcessor.convolutionFFT(xShort, yShort);
//        double[] corr = DSPProcessor.correlation(xShort, yShort);
//        double[] corrFFT = DSPProcessor.correlationFFT(xShort, yShort);
//        double[] convFFTLib = DSPProcessor.convolutionFFTLibrary(xShort, yShort);
//        double[] corrFFTLib = DSPProcessor.correlationFFTLibrary(xShort, yShort);
//
//        graphsPanel.add(ChartUtils.createMiniChartPanel(SignalData.x, "1. x(t)", "Время, мс", "Амплитуда",
//                ChartUtils.FOREST_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(SignalData.y, "2. y(t)", "Время, мс", "Амплитуда",
//                ChartUtils.TERRA_COTTA, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.amplitudeSpectrum(xDFT), "3. x ДПФ Амп",
//                "Частота, Гц", "Амплитуда", ChartUtils.FOREST_GREEN, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.phaseForPlot(xDFT, 0.1), "4. x ДПФ Фаза",
//                "Частота, Гц", "Фаза, рад", ChartUtils.SLATE_BLUE, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(xIDFT, "5. x ОДПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.amplitudeSpectrum(xFFT), "6. x БПФ Амп",
//                "Частота, Гц", "Амплитуда", ChartUtils.FOREST_GREEN, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.phaseForPlot(xFFT, 0.1), "7. x БПФ Фаза",
//                "Частота, Гц", "Фаза, рад", ChartUtils.SLATE_BLUE, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(xIFFT, "8. x ОБПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.amplitudeSpectrum(yDFT), "9. y ДПФ Амп",
//                "Частота, Гц", "Амплитуда", ChartUtils.TERRA_COTTA, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.phaseForPlot(yDFT, 0.1), "10. y ДПФ Фаза",
//                "Частота, Гц", "Фаза, рад", ChartUtils.DARK_SLATE, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(yIDFT, "11. y ОДПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.amplitudeSpectrum(yFFT), "12. y БПФ Амп",
//                "Частота, Гц", "Амплитуда", ChartUtils.TERRA_COTTA, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.phaseForPlot(yFFT, 0.1), "13. y БПФ Фаза",
//                "Частота, Гц", "Фаза, рад", ChartUtils.DARK_SLATE, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(yIFFT, "14. y ОБПФ", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(conv, "15. Свертка", "Время, мс", "Амплитуда",
//                ChartUtils.FOREST_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(convFFT, "16. Свертка FFT", "Время, мс", "Амплитуда",
//                ChartUtils.EARTH_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(corr, "17. Корреляция", "Время, мс", "Амплитуда",
//                ChartUtils.SLATE_BLUE, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(corrFFT, "18. Корреляция FFT", "Время, мс", "Амплитуда",
//                ChartUtils.DARK_SLATE, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.amplitudeSpectrum(xFFTLib), "19. x БПФ библ Амп",
//                "Частота, Гц", "Амплитуда", ChartUtils.FOREST_GREEN, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.phaseForPlot(xFFTLib, 0.1), "20. x БПФ библ Фаза",
//                "Частота, Гц", "Фаза, рад", ChartUtils.SLATE_BLUE, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.amplitudeSpectrum(yFFTLib), "21. y БПФ библ Амп",
//                "Частота, Гц", "Амплитуда", ChartUtils.TERRA_COTTA, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(DSPProcessor.phaseForPlot(yFFTLib, 0.1), "22. y БПФ библ Фаза",
//                "Частота, Гц", "Фаза, рад", ChartUtils.DARK_SLATE, false, true));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(convFFTLib, "23. Свертка библ", "Время, мс", "Амплитуда",
//                ChartUtils.FOREST_GREEN, true, false));
//        graphsPanel.add(ChartUtils.createMiniChartPanel(corrFFTLib, "24. Корреляция библ", "Время, мс", "Амплитуда",
//                ChartUtils.SLATE_BLUE, true, false));

        mainPanel.add(graphsPanel, BorderLayout.CENTER);

        JPanel infoPanel = new JPanel(new BorderLayout());
        infoPanel.setBackground(ChartUtils.LIGHT_GRAY);
        infoPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createMatteBorder(1, 0, 0, 0, ChartUtils.MEDIUM_GRAY),
                BorderFactory.createEmptyBorder(5, 10, 5, 10)));

        JLabel infoLabel = new JLabel(
                "Обзор всех 24 графиков • Колесо мыши для масштабирования • Двойной щелчок для сброса");
        infoLabel.setFont(new Font("Arial", Font.PLAIN, 10));
        infoLabel.setForeground(ChartUtils.DARK_BROWN);
        infoPanel.add(infoLabel, BorderLayout.CENTER);

        mainPanel.add(infoPanel, BorderLayout.SOUTH);
        return mainPanel;
    }
}