package org.example;

import javax.sound.sampled.UnsupportedAudioFileException;
import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

public class Lab3MusicComparePanel extends JPanel {
    private static final int FFT_SIZE = 512;
    private static final int HOP_SIZE = 256;
    private static final Dimension SPEC_SIZE = new Dimension(340, 260);
    private static final String LAB_ROOT_DIR_NAME = "Lab3_COSI_audio";

    private final JTextField[] genreFields = new JTextField[4];
    private final JTextField[] pathFields = new JTextField[4];
    private final JPanel comparisonGrid = new JPanel();

    public Lab3MusicComparePanel() {
        setLayout(new BorderLayout(8, 8));
        setBackground(ChartUtils.LIGHT_BEIGE);
        setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));

        add(buildNorth(), BorderLayout.NORTH);

        comparisonGrid.setBackground(ChartUtils.OFF_WHITE);
        comparisonGrid.setLayout(new BoxLayout(comparisonGrid, BoxLayout.Y_AXIS));
        comparisonGrid.add(new JLabel(
                "<html><body style='width:720px'>Выберите до четырёх WAV (разные жанры), укажите жанр в поле — "
                        + "нажмите «Показать сравнение»: в каждом столбце под названием жанра — числа и обе спектрограммы.</body></html>"));

        JScrollPane scroll = new JScrollPane(comparisonGrid);
        scroll.getVerticalScrollBar().setUnitIncrement(16);
        add(scroll, BorderLayout.CENTER);
    }

    private JPanel buildNorth() {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBackground(ChartUtils.CREAM);
        panel.setBorder(BorderFactory.createTitledBorder("П. 2–3: сравнение музыкальных фрагментов"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(3, 4, 3, 4);
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;

        String[] hints = {"напр. классика", "напр. рок", "напр. электроника", "напр. джаз"};
        for (int i = 0; i < 4; i++) {
            int row = i;
            gbc.gridy = row;
            gbc.gridx = 0;
            gbc.weightx = 0;
            panel.add(new JLabel("Файл " + (i + 1) + ", жанр:"), gbc);
            genreFields[i] = new JTextField(hints[i], 12);
            gbc.gridx = 1;
            gbc.weightx = 0.15;
            panel.add(genreFields[i], gbc);
            pathFields[i] = new JTextField();
            pathFields[i].setEditable(false);
            gbc.gridx = 2;
            gbc.weightx = 1.0;
            panel.add(pathFields[i], gbc);
            int idx = i;
            JButton pick = new JButton("WAV…");
            pick.addActionListener(e -> pickWav(idx));
            gbc.gridx = 3;
            gbc.weightx = 0;
            panel.add(pick, gbc);

            JButton chartsBtn = new JButton("Графики…");
            chartsBtn.addActionListener(e -> showCharacteristicsWindow(idx));
            gbc.gridx = 4;
            panel.add(chartsBtn, gbc);
        }

        gbc.gridy = 4;
        gbc.gridx = 0;
        gbc.gridwidth = 5;
        JButton run = new JButton("Показать сравнение");
        run.addActionListener(e -> showComparison());
        panel.add(run, gbc);

        return panel;
    }

    private void pickWav(int index) {
        JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle("Музыкальный фрагмент " + (index + 1) + " (WAV)");
        if (chooser.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        File f = chooser.getSelectedFile();
        pathFields[index].setText(f.getAbsolutePath());
    }

    private void showComparison() {
        comparisonGrid.removeAll();
        comparisonGrid.setLayout(new BoxLayout(comparisonGrid, BoxLayout.Y_AXIS));

        AudioIOUtils.AudioData[] loaded = new AudioIOUtils.AudioData[4];
        Lab3AudioProcessor.AvgSpectralFeatures[] features = new Lab3AudioProcessor.AvgSpectralFeatures[4];
        LibrosaBridge.Result[] librosa = new LibrosaBridge.Result[4];
        boolean[] readError = new boolean[4];

        for (int i = 0; i < 4; i++) {
            String path = pathFields[i].getText().trim();
            if (path.isEmpty()) {
                continue;
            }
            try {
                loaded[i] = AudioIOUtils.readMonoWav(new File(path));
                features[i] = Lab3AudioProcessor.averageSpectralFeatures(
                        loaded[i].samples(), loaded[i].sampleRate(), FFT_SIZE, HOP_SIZE);
                librosa[i] = analyzeWithLibrosa(path);
            } catch (UnsupportedAudioFileException | IOException ex) {
                readError[i] = true;
            }
        }

        JPanel columnsRow = new JPanel(new GridLayout(1, 4, 10, 0));
        columnsRow.setBackground(ChartUtils.OFF_WHITE);
        columnsRow.setAlignmentX(Component.LEFT_ALIGNMENT);
        columnsRow.setMaximumSize(new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE));

        for (int i = 0; i < 4; i++) {
            columnsRow.add(buildColumnPanel(i, loaded[i], features[i], librosa[i], readError[i]));
        }

        comparisonGrid.add(columnsRow);
        comparisonGrid.revalidate();
        comparisonGrid.repaint();
    }

    private JPanel buildColumnPanel(int index, AudioIOUtils.AudioData data,
                                    Lab3AudioProcessor.AvgSpectralFeatures features,
                                    LibrosaBridge.Result librosa,
                                    boolean readError) {
        JPanel col = new JPanel();
        col.setLayout(new BoxLayout(col, BoxLayout.Y_AXIS));
        col.setBackground(ChartUtils.CREAM);
        col.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY),
                BorderFactory.createEmptyBorder(8, 8, 8, 8)));

        String genre = genreFields[index].getText().trim();
        if (genre.isEmpty()) {
            genre = "файл " + (index + 1);
        }
        String path = pathFields[index].getText().trim();

        JPanel header = new JPanel(new BorderLayout());
        header.setOpaque(false);
        if (path.isEmpty()) {
            header.add(new JLabel("<html><div style='text-align:center'><b>" + genre + "</b><br/>"
                    + "<span style='color:#888'>не выбран</span></div></html>", SwingConstants.CENTER),
                    BorderLayout.CENTER);
        } else if (readError) {
            header.add(new JLabel("<html><div style='text-align:center'><b>" + genre + "</b><br/>"
                    + "<span style='color:#b00'>ошибка чтения</span></div></html>", SwingConstants.CENTER),
                    BorderLayout.CENTER);
        } else {
            String shortName = new File(path).getName();
            header.add(new JLabel("<html><div style='text-align:center'><b>" + genre + "</b><br/>"
                    + "<small>" + shortName + "</small></div></html>", SwingConstants.CENTER),
                    BorderLayout.CENTER);
        }
        header.setAlignmentX(Component.CENTER_ALIGNMENT);
        col.add(header);

        col.add(Box.createVerticalStrut(10));

        JLabel cTitle = subsectionTitle("Средний спектральный центроид");
        JLabel cVal = metricValue(data == null || readError
                ? "—"
                : String.format(Locale.US, "%.1f Гц", features.centroidHz()));

        JLabel bTitle = subsectionTitle("Средняя спектральная ширина");
        JLabel bVal = metricValue(data == null || readError
                ? "—"
                : String.format(Locale.US, "%.1f Гц", features.bandwidthHz()));

        col.add(cTitle);
        col.add(cVal);
        col.add(Box.createVerticalStrut(6));
        col.add(bTitle);
        col.add(bVal);
        col.add(Box.createVerticalStrut(8));

        col.add(subsectionTitle("Спектральные признаки (librosa)"));
        col.add(Box.createVerticalStrut(4));
        col.add(metricValue("SC (центроид): " + librosaValue(librosa, "centroid")));
        col.add(metricValue("SR (спад 85%): " + librosaValue(librosa, "rolloff")));
        col.add(metricValue("SB (ширина): " + librosaValue(librosa, "bandwidth")));
        col.add(metricValue("ZCR: " + librosaValue(librosa, "zcr")));
        col.add(Box.createVerticalStrut(8));

        col.add(Box.createVerticalStrut(2));

        col.add(subsectionTitle("Самописная спектрограмма (DFT + Ханн)"));
        col.add(Box.createVerticalStrut(4));
        col.add(wrapChartCentered(spectrogramChart(data, true)));
        col.add(Box.createVerticalStrut(8));
        col.add(subsectionTitle("Мел-спектрограмма (librosa)"));
        col.add(Box.createVerticalStrut(4));
        col.add(wrapChartCentered(melSpectrogramChart(librosa)));

        return col;
    }

    private void showCharacteristicsWindow(int index) {
        String path = pathFields[index].getText().trim();
        if (path.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Сначала выберите WAV для файла " + (index + 1) + ".");
            return;
        }

        try {
            AudioIOUtils.AudioData audio = AudioIOUtils.readMonoWav(new File(path));
            LibrosaBridge.Result librosa = analyzeWithLibrosa(path);
            if (librosa.error() != null && !librosa.error().isBlank()) {
                JOptionPane.showMessageDialog(this, "Ошибка librosa: " + librosa.error());
                return;
            }

            double[] centroidManual = manualCentroidSeries(audio.samples(), audio.sampleRate());
            double[] bandwidthManual = manualBandwidthSeries(audio.samples(), audio.sampleRate());

            double[] centroidLib = librosa.centroidSeriesHz();
            double[] bandwidthLib = librosa.bandwidthSeriesHz();
            double[] rolloffLib = librosa.rolloffSeriesHz();
            double[] zcrLib = librosa.zcrSeries();
            double[] tManual = frameTimeAxis(centroidManual.length, audio.sampleRate(), HOP_SIZE);
            int libHop = librosa.hopLength() > 0 ? librosa.hopLength() : HOP_SIZE;
            double[] tLibCentroid = frameTimeAxis(nonNullSeries(centroidLib).length, audio.sampleRate(), libHop);
            double[] tLibBandwidth = frameTimeAxis(nonNullSeries(bandwidthLib).length, audio.sampleRate(), libHop);
            double[] tLibRolloff = frameTimeAxis(nonNullSeries(rolloffLib).length, audio.sampleRate(), libHop);
            double[] tLibZcr = frameTimeAxis(nonNullSeries(zcrLib).length, audio.sampleRate(), libHop);

            String genre = genreFields[index].getText().trim();
            if (genre.isEmpty()) genre = "файл " + (index + 1);
            String shortName = new File(path).getName();

            JFrame frame = new JFrame("Характеристики: " + genre + " (" + shortName + ")");
            frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
            frame.setSize(1500, 900);
            frame.setLocationRelativeTo(this);

            JPanel grid = new JPanel(new GridLayout(3, 2, 8, 8));
            grid.setBackground(ChartUtils.LIGHT_BEIGE);
            grid.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));

            grid.add(ChartUtils.createXYSeriesChartPanel(
                    tManual, centroidManual, "Центроид (самописный)", "Время, с", "Гц",
                    ChartUtils.FOREST_GREEN));
            grid.add(ChartUtils.createXYSeriesChartPanel(
                    tLibCentroid, nonNullSeries(centroidLib), "Центроид (librosa)", "Время, с", "Гц",
                    ChartUtils.SLATE_BLUE));

            grid.add(ChartUtils.createXYSeriesChartPanel(
                    tManual, bandwidthManual, "Ширина (самописная)", "Время, с", "Гц",
                    ChartUtils.EARTH_GREEN));
            grid.add(ChartUtils.createXYSeriesChartPanel(
                    tLibBandwidth, nonNullSeries(bandwidthLib), "Ширина (librosa)", "Время, с", "Гц",
                    ChartUtils.TERRA_COTTA));

            grid.add(ChartUtils.createXYSeriesChartPanel(
                    tLibRolloff, nonNullSeries(rolloffLib), "Спад / Rolloff (librosa)", "Время, с", "Гц",
                    ChartUtils.DARK_SLATE));
            grid.add(ChartUtils.createXYSeriesChartPanel(
                    tLibZcr, nonNullSeries(zcrLib), "Пересечение нуля / ZCR (librosa)", "Время, с", "Доля",
                    ChartUtils.WARM_GRAY));

            frame.setContentPane(grid);
            frame.setVisible(true);
        } catch (UnsupportedAudioFileException | IOException ex) {
            JOptionPane.showMessageDialog(this, "Ошибка чтения WAV: " + ex.getMessage(),
                    "Ошибка", JOptionPane.ERROR_MESSAGE);
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Не удалось построить графики: " + ex.getMessage(),
                    "Ошибка", JOptionPane.ERROR_MESSAGE);
        }
    }

    private JLabel subsectionTitle(String text) {
        JLabel l = new JLabel("<html><div style='text-align:center'>" + text + "</div></html>");
        l.setFont(new Font("Arial", Font.BOLD, 11));
        l.setForeground(ChartUtils.DARK_BROWN);
        l.setAlignmentX(Component.CENTER_ALIGNMENT);
        return l;
    }

    private JLabel metricValue(String text) {
        JLabel l = new JLabel(text);
        l.setFont(new Font("Arial", Font.PLAIN, 13));
        l.setForeground(ChartUtils.DARK_BROWN);
        l.setAlignmentX(Component.CENTER_ALIGNMENT);
        l.setHorizontalAlignment(SwingConstants.CENTER);
        return l;
    }

    private JPanel wrapChartCentered(JPanel chartOrPlaceholder) {
        JPanel wrap = new JPanel(new FlowLayout(FlowLayout.CENTER, 0, 0));
        wrap.setOpaque(false);
        wrap.setAlignmentX(Component.CENTER_ALIGNMENT);
        wrap.add(chartOrPlaceholder);
        return wrap;
    }

    private JPanel spectrogramChart(AudioIOUtils.AudioData data, boolean manual) {
        if (data == null) {
            JPanel empty = new JPanel(new BorderLayout());
            empty.setBackground(ChartUtils.OFF_WHITE);
            empty.setBorder(BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY));
            empty.setPreferredSize(SPEC_SIZE);
            empty.setMinimumSize(SPEC_SIZE);
            empty.add(new JLabel("—", SwingConstants.CENTER), BorderLayout.CENTER);
            return empty;
        }
        double[] signal = data.samples();
        int sr = data.sampleRate();
        double[][] spec = manual
                ? Lab3AudioProcessor.spectrogramManual(signal, FFT_SIZE, HOP_SIZE)
                : Lab3AudioProcessor.spectrogramLibrary(signal, FFT_SIZE, HOP_SIZE);
        String title = manual ? "Самописная" : "Библиотечная";
        return ChartUtils.createSpectrogramChartPanel(spec, title, sr, HOP_SIZE, SPEC_SIZE);
    }

    private JPanel melSpectrogramChart(LibrosaBridge.Result librosa) {
        if (librosa == null || librosa.melDb() == null) {
            JPanel empty = new JPanel(new BorderLayout());
            empty.setBackground(ChartUtils.OFF_WHITE);
            empty.setBorder(BorderFactory.createLineBorder(ChartUtils.MEDIUM_GRAY));
            empty.setPreferredSize(SPEC_SIZE);
            empty.setMinimumSize(SPEC_SIZE);
            String text = (librosa != null && librosa.error() != null && !librosa.error().isBlank())
                    ? "<html><div style='text-align:center'>librosa недоступна<br/><small>" + escapeHtml(librosa.error()) + "</small></div></html>"
                    : "—";
            empty.add(new JLabel(text, SwingConstants.CENTER), BorderLayout.CENTER);
            return empty;
        }
        return ChartUtils.createMelSpectrogramChartPanel(
                librosa.melDb(),
                "Mel (librosa)",
                librosa.sampleRate(),
                librosa.hopLength() > 0 ? librosa.hopLength() : 512,
                SPEC_SIZE
        );
    }

    private LibrosaBridge.Result analyzeWithLibrosa(String wavPath) {
        try {
            Path root = Paths.get(System.getProperty("user.home", ".")).resolve(LAB_ROOT_DIR_NAME);
            Path cacheDir = root.resolve("music_librosa_cache");
            Files.createDirectories(cacheDir);
            String safeName = new File(wavPath).getName().replaceAll("[^a-zA-Z0-9._-]", "_");
            Path json = cacheDir.resolve(safeName + ".json");
            return LibrosaBridge.analyzeWavWithLibrosa(wavPath, json);
        } catch (Exception ex) {
            return new LibrosaBridge.Result(0, 0, null, 0, 0, 0, 0, null, null, null, null, ex.getMessage());
        }
    }

    private double[] manualCentroidSeries(double[] signal, int sampleRate) {
        double[][] spec = Lab3AudioProcessor.spectrogramManual(signal, FFT_SIZE, HOP_SIZE);
        double[] out = new double[spec.length];
        for (int i = 0; i < spec.length; i++) {
            out[i] = Lab3AudioProcessor.spectralCentroid(spec[i], sampleRate, FFT_SIZE);
        }
        return out;
    }

    private double[] manualBandwidthSeries(double[] signal, int sampleRate) {
        double[][] spec = Lab3AudioProcessor.spectrogramManual(signal, FFT_SIZE, HOP_SIZE);
        double[] out = new double[spec.length];
        for (int i = 0; i < spec.length; i++) {
            double c = Lab3AudioProcessor.spectralCentroid(spec[i], sampleRate, FFT_SIZE);
            out[i] = Lab3AudioProcessor.spectralBandwidth(spec[i], sampleRate, FFT_SIZE, c);
        }
        return out;
    }

    private double[] nonNullSeries(double[] series) {
        if (series == null || series.length == 0) return new double[]{0.0};
        return series;
    }

    private double[] frameTimeAxis(int frames, int sampleRate, int hopSize) {
        if (frames <= 0 || sampleRate <= 0 || hopSize <= 0) return new double[]{0.0};
        double[] t = new double[frames];
        double dt = hopSize / (double) sampleRate;
        for (int i = 0; i < frames; i++) t[i] = i * dt;
        return t;
    }

    private String librosaValue(LibrosaBridge.Result r, String key) {
        if (r == null) return "—";
        if (r.error() != null && !r.error().isBlank()) return "ошибка";
        return switch (key) {
            case "centroid" -> String.format(Locale.US, "%.1f Гц", r.centroidHz());
            case "rolloff" -> String.format(Locale.US, "%.1f Гц", r.rolloffHz());
            case "bandwidth" -> String.format(Locale.US, "%.1f Гц", r.bandwidthHz());
            case "zcr" -> String.format(Locale.US, "%.4f", r.zcr());
            default -> "—";
        };
    }

    private static String escapeHtml(String s) {
        return s == null ? "" : s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;");
    }
}
