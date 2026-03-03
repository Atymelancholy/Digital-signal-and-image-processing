package org.example;

import javax.sound.sampled.*;
import javax.swing.*;
import java.awt.*;
import java.util.Arrays;

public class StartPanel extends JPanel {
    private final MainFrame mainFrame;

    public StartPanel(MainFrame mainFrame) {
        this.mainFrame = mainFrame;
        setLayout(new BorderLayout(10, 10));
        setBackground(ChartUtils.LIGHT_BEIGE);
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        JLabel titleLabel = new JLabel("Главное меню — исходные сигналы", SwingConstants.CENTER);
        titleLabel.setFont(new Font("Arial", Font.BOLD, 18));
        titleLabel.setForeground(ChartUtils.DARK_BROWN);
        titleLabel.setBorder(BorderFactory.createEmptyBorder(0, 0, 10, 0));
        add(titleLabel, BorderLayout.NORTH);

        JPanel signalsPanel = createSignalsPanel();
        add(signalsPanel, BorderLayout.CENTER);

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 30, 10));
        buttonPanel.setBackground(ChartUtils.LIGHT_BEIGE);

        JButton lab1Button = new JButton("Лабораторная работа №1");
        JButton lab2Button = new JButton("Лабораторная работа №2");

        Font buttonFont = new Font("Arial", Font.BOLD, 14);
        lab1Button.setFont(buttonFont);
        lab2Button.setFont(buttonFont);
        lab1Button.setBackground(ChartUtils.CREAM);
        lab2Button.setBackground(ChartUtils.CREAM);
        lab1Button.setForeground(ChartUtils.DARK_BROWN);
        lab2Button.setForeground(ChartUtils.DARK_BROWN);
        lab1Button.setFocusPainted(false);
        lab2Button.setFocusPainted(false);
        lab1Button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(8, 20, 8, 20)));
        lab2Button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(8, 20, 8, 20)));

        lab1Button.addActionListener(e -> mainFrame.showLab1());
        lab2Button.addActionListener(e -> mainFrame.showLab2());

        buttonPanel.add(lab1Button);
        buttonPanel.add(lab2Button);

        add(buttonPanel, BorderLayout.SOUTH);
    }

    private JPanel createSignalsPanel() {
        JPanel panel = new JPanel(new GridLayout(1, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        int displayLen = Math.min(2000, SignalData.N);
        double[] xShort = Arrays.copyOf(SignalData.x, displayLen);
        double[] yShort = Arrays.copyOf(SignalData.y, displayLen);

        panel.add(createPlayableChartPanel(
                xShort,
                "x(t) - исходный сигнал",
                ChartUtils.FOREST_GREEN,
                SignalData.x
        ));

        panel.add(createPlayableChartPanel(
                yShort,
                "y(t) - исходный сигнал",
                ChartUtils.TERRA_COTTA,
                SignalData.y
        ));

        return panel;
    }

    private JPanel createPlayableChartPanel(double[] displayData, String title, Color color, double[] fullSignal) {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(ChartUtils.OFF_WHITE);
        panel.setBorder(BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY, 1));

        JPanel chartPanel = ChartUtils.createDynamicChartPanel(
                displayData,
                title,
                "Время, мс",
                "Амплитуда",
                color,
                true,
                false
        );

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

    private void playSignal(double[] signal) {
        final int sampleRate = SignalData.SAMPLE_RATE;
        final int bits = 16;
        final int channels = 1;
        final boolean signed = true;
        final boolean bigEndian = false;

        double max = 0.0;
        for (double v : signal) {
            if (Math.abs(v) > max) {
                max = Math.abs(v);
            }
        }

        if (max == 0.0) {
            return;
        }

        final double amplification = 0.9;
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
            JOptionPane.showMessageDialog(
                    this,
                    "Не удалось воспроизвести звук: " + ex.getMessage(),
                    "Ошибка",
                    JOptionPane.ERROR_MESSAGE
            );
        }
    }
}