package org.example;

import javax.swing.*;
import java.awt.*;

public class MainFrame extends JFrame {
    private final StartPanel startPanel;
    private final Lab1Panel lab1Panel;
    private final Lab2Panel lab2Panel;

    public MainFrame() {
        super("Лабораторные работы по ЦОС");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(1600, 900);
        setLocationRelativeTo(null);

        startPanel = new StartPanel(this);
        lab1Panel = new Lab1Panel();
        lab2Panel = new Lab2Panel();

        JMenuBar menuBar = new JMenuBar();
        menuBar.setBackground(ChartUtils.LIGHT_GRAY);

        JMenu mainMenu = new JMenu("Главная");
        mainMenu.setFont(new Font("Arial", Font.BOLD, 12));
        mainMenu.setForeground(ChartUtils.DARK_BROWN);
        JMenuItem startItem = new JMenuItem("Стартовый экран");
        startItem.addActionListener(e -> showStartMenu());
        mainMenu.add(startItem);
        menuBar.add(mainMenu);

        JMenu labMenu = new JMenu("Лабораторные работы");
        labMenu.setFont(new Font("Arial", Font.BOLD, 12));
        labMenu.setForeground(ChartUtils.DARK_BROWN);

        JMenuItem lab1Item = new JMenuItem("ЛР №1: Алгоритмы ЦОС");
        lab1Item.addActionListener(e -> showLab1());
        labMenu.add(lab1Item);

        JMenuItem lab2Item = new JMenuItem("ЛР №2: Проектирование фильтров");
        lab2Item.addActionListener(e -> showLab2());
        labMenu.add(lab2Item);

        menuBar.add(labMenu);
        setJMenuBar(menuBar);

        showStartMenu();
    }

    public void showStartMenu() {
        setContentPane(startPanel);
        revalidate();
        repaint();
    }

    public void showLab1() {
        setContentPane(lab1Panel);
        revalidate();
        repaint();
    }

    public void showLab2() {
        setContentPane(lab2Panel);
        revalidate();
        repaint();
    }

    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
        }
        UIManager.put("Panel.background", ChartUtils.LIGHT_BEIGE);
        UIManager.put("TextField.background", ChartUtils.CREAM);
        UIManager.put("TextField.foreground", ChartUtils.DARK_BROWN);
        UIManager.put("TextArea.background", ChartUtils.CREAM);
        UIManager.put("TextArea.foreground", ChartUtils.DARK_BROWN);

        SwingUtilities.invokeLater(() -> new MainFrame().setVisible(true));
    }
}