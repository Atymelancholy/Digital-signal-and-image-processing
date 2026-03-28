package org.example;

import javax.swing.*;
import java.awt.*;

public class Lab3Panel extends JPanel {

    public Lab3Panel() {
        setLayout(new BorderLayout());
        setBackground(ChartUtils.LIGHT_BEIGE);
        setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));

        JTabbedPane tabs = new JTabbedPane();
        tabs.setBackground(ChartUtils.OFF_WHITE);
        tabs.setFont(new Font("Arial", Font.PLAIN, 12));

        tabs.addTab("П. 2–3: музыка (4 жанра, всё рядом)", new Lab3MusicComparePanel());
        tabs.addTab("П. 4–6: голос, шум, признаки, DFN", new Lab3VoicePanel());

        add(tabs, BorderLayout.CENTER);
    }
}

//https://www.boxentriq.com/steganography/audio-spectrogram