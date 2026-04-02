package org.example;

import javax.swing.*;
import java.awt.*;

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
        JButton lab3Button = new JButton("Лабораторная работа №3");

        Font buttonFont = new Font("Arial", Font.BOLD, 14);
        lab1Button.setFont(buttonFont);
        lab2Button.setFont(buttonFont);
        lab3Button.setFont(buttonFont);
        lab1Button.setBackground(ChartUtils.CREAM);
        lab2Button.setBackground(ChartUtils.CREAM);
        lab3Button.setBackground(ChartUtils.CREAM);
        lab1Button.setForeground(ChartUtils.DARK_BROWN);
        lab2Button.setForeground(ChartUtils.DARK_BROWN);
        lab3Button.setForeground(ChartUtils.DARK_BROWN);
        lab1Button.setFocusPainted(false);
        lab2Button.setFocusPainted(false);
        lab3Button.setFocusPainted(false);
        lab1Button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(8, 20, 8, 20)));
        lab2Button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(8, 20, 8, 20)));
        lab3Button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(8, 20, 8, 20)));

        lab1Button.addActionListener(e -> mainFrame.showLab1());
        lab2Button.addActionListener(e -> mainFrame.showLab2());
        lab3Button.addActionListener(e -> mainFrame.showLab3());

        buttonPanel.add(lab1Button);
        buttonPanel.add(lab2Button);
        buttonPanel.add(lab3Button);

        add(buttonPanel, BorderLayout.SOUTH);
    }

    private JPanel createSignalsPanel() {
        JPanel panel = new JPanel(new GridLayout(1, 2, 10, 10));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        int displayLen = Math.min(2000, SignalData.N);
        double[] xShort = java.util.Arrays.copyOf(SignalData.x, displayLen);
        double[] yShort = java.util.Arrays.copyOf(SignalData.y, displayLen);

        panel.add(ChartUtils.createDynamicChartPanel(xShort, "x(t) - исходный сигнал",
                "Время, мс", "Амплитуда", ChartUtils.FOREST_GREEN, true, false));
        panel.add(ChartUtils.createDynamicChartPanel(yShort, "y(t) - исходный сигнал",
                "Время, мс", "Амплитуда", ChartUtils.TERRA_COTTA, true, false));

        return panel;
    }
}