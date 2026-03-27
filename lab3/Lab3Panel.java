package org.example;

import javax.sound.sampled.UnsupportedAudioFileException;
import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class Lab3Panel extends JPanel {
    private static final int FFT_SIZE = 512;
    private static final int HOP_SIZE = 256;

    private AudioIOUtils.AudioData cleanAudio;
    private AudioIOUtils.AudioData noiseAudio;
    private double[] currentMixed;

    private final JTextField cleanPathField = new JTextField();
    private final JTextField noisePathField = new JTextField();
    private final JTextField deepFilterCommandField = new JTextField("py -3.11 -m df.enhance \"{in}\" -o \"{outDir}\" --no-suffix");
    private final JTextArea infoArea = new JTextArea();

    private final JPanel manualSpecContainer = new JPanel(new BorderLayout());
    private final JPanel librarySpecContainer = new JPanel(new BorderLayout());
    private final DefaultTableModel metricsModel = new DefaultTableModel(
            new Object[]{"SNR вход, дБ", "SNR, дБ", "SI-SDR, дБ", "SNR (после DFN2)", "SI-SDR (после DFN2)"}, 0);

    private File lastMixedFile;
    private final List<MixedCase> mixedCases = new ArrayList<>();

    private record MixedCase(int snrInputDb, int rowIndex, File mixedFile) {}

    public Lab3Panel() {
        setLayout(new BorderLayout(8, 8));
        setBackground(ChartUtils.LIGHT_BEIGE);
        setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));

        add(createControlsPanel(), BorderLayout.NORTH);
        add(createTabsPanel(), BorderLayout.CENTER);
    }

    private JPanel createControlsPanel() {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBackground(ChartUtils.CREAM);
        panel.setBorder(BorderFactory.createTitledBorder("ЛР №3: входные данные и запуск"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weightx = 0;

        cleanPathField.setEditable(false);
        noisePathField.setEditable(false);

        JButton cleanBtn = new JButton("Выбрать чистый WAV");
        JButton noiseBtn = new JButton("Выбрать шум WAV");
        JButton runBtn = new JButton("Выполнить вариант (0..18, шаг 3)");
        JButton deepFilterBtn = new JButton("Применить DeepFilterNet2 ко всем SNR (0..18)");

        cleanBtn.addActionListener(e -> pickWav(cleanPathField, true));
        noiseBtn.addActionListener(e -> pickWav(noisePathField, false));
        runBtn.addActionListener(e -> runVariant());
        deepFilterBtn.addActionListener(e -> runDeepFilterOnAllMixes());

        int row = 0;
        gbc.gridx = 0; gbc.gridy = row; panel.add(new JLabel("Чистый голос:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1.0; panel.add(cleanPathField, gbc);
        gbc.gridx = 2; gbc.weightx = 0; panel.add(cleanBtn, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row; panel.add(new JLabel("Шум:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1.0; panel.add(noisePathField, gbc);
        gbc.gridx = 2; gbc.weightx = 0; panel.add(noiseBtn, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row; panel.add(new JLabel("Команда DeepFilterNet2:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1.0; panel.add(deepFilterCommandField, gbc);
        gbc.gridx = 2; gbc.weightx = 0; panel.add(deepFilterBtn, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 3;
        panel.add(runBtn, gbc);

        return panel;
    }

    private JTabbedPane createTabsPanel() {
        JTabbedPane tabs = new JTabbedPane();
        tabs.setBackground(ChartUtils.OFF_WHITE);
        tabs.setFont(new Font("Arial", Font.PLAIN, 12));

        manualSpecContainer.setBackground(ChartUtils.OFF_WHITE);
        librarySpecContainer.setBackground(ChartUtils.OFF_WHITE);
        manualSpecContainer.add(new JLabel("После запуска здесь будет самописная спектрограмма", SwingConstants.CENTER),
                BorderLayout.CENTER);
        librarySpecContainer.add(new JLabel("После запуска здесь будет библиотечная спектрограмма", SwingConstants.CENTER),
                BorderLayout.CENTER);

        JPanel specPanel = new JPanel(new GridLayout(1, 2, 8, 8));
        specPanel.setBackground(ChartUtils.LIGHT_BEIGE);
        specPanel.add(manualSpecContainer);
        specPanel.add(librarySpecContainer);

        JTable metricsTable = new JTable(metricsModel);
        JScrollPane tableScroll = new JScrollPane(metricsTable);
        tableScroll.setBorder(BorderFactory.createTitledBorder("Метрики по диапазону SNR"));

        infoArea.setEditable(false);
        infoArea.setLineWrap(true);
        infoArea.setWrapStyleWord(true);
        infoArea.setFont(new Font("Consolas", Font.PLAIN, 12));
        infoArea.setBackground(ChartUtils.CREAM);
        infoArea.setBorder(BorderFactory.createTitledBorder("Признаки (самописные) + комментарий"));

        JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, tableScroll, new JScrollPane(infoArea));
        splitPane.setResizeWeight(0.55);

        tabs.addTab("Спектрограмма", specPanel);
        tabs.addTab("Признаки и метрики", splitPane);
        return tabs;
    }

    private void pickWav(JTextField targetField, boolean clean) {
        JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle(clean ? "Выберите чистый голос (WAV)" : "Выберите шум (WAV)");
        int result = chooser.showOpenDialog(this);
        if (result != JFileChooser.APPROVE_OPTION) return;
        File file = chooser.getSelectedFile();
        targetField.setText(file.getAbsolutePath());
        try {
            AudioIOUtils.AudioData data = AudioIOUtils.readMonoWav(file);
            if (clean) cleanAudio = data; else noiseAudio = data;
            appendInfo(String.format(Locale.US, "%s загружен: %s (%d Гц, %d сэмплов)",
                    clean ? "Чистый сигнал" : "Шум", file.getName(), data.sampleRate(), data.samples().length));
        } catch (UnsupportedAudioFileException | IOException ex) {
            JOptionPane.showMessageDialog(this, "Ошибка чтения WAV: " + ex.getMessage(),
                    "Ошибка", JOptionPane.ERROR_MESSAGE);
        }
    }

    private void runVariant() {
        if (cleanAudio == null || noiseAudio == null) {
            JOptionPane.showMessageDialog(this, "Сначала выберите чистый и шумовой WAV файлы.");
            return;
        }
        if (cleanAudio.sampleRate() != noiseAudio.sampleRate()) {
            JOptionPane.showMessageDialog(this, "Частоты дискретизации не совпадают. Выберите файлы с одинаковым sample rate.");
            return;
        }

        metricsModel.setRowCount(0);
        mixedCases.clear();
        double[] clean = cleanAudio.samples();
        double[] noise = noiseAudio.samples();

        for (int snr = 0; snr <= 18; snr += 3) {
            double[] mixed = Lab3AudioProcessor.mixBySnrDb(clean, noise, snr);
            double snrValue = Lab3AudioProcessor.snrDb(clean, mixed);
            double siSdrValue = Lab3AudioProcessor.siSdrDb(clean, mixed);
            int rowIndex = metricsModel.getRowCount();
            metricsModel.addRow(new Object[]{
                    snr,
                    formatDouble(snrValue),
                    formatDouble(siSdrValue),
                    "—",
                    "—"
            });
            try {
                Path tempMixed = Files.createTempFile("lab3_mixed_snr" + snr + "_", ".wav");
                AudioIOUtils.writeMonoWav(tempMixed.toFile(), mixed, cleanAudio.sampleRate());
                mixedCases.add(new MixedCase(snr, rowIndex, tempMixed.toFile()));
                if (snr == 6) {
                    currentMixed = mixed;
                    lastMixedFile = tempMixed.toFile();
                }
            } catch (IOException ex) {
                appendInfo("Не удалось сохранить микс для SNR=" + snr + ": " + ex.getMessage());
            }
            if (snr == 6) currentMixed = mixed;
        }

        if (currentMixed == null) {
            currentMixed = Lab3AudioProcessor.mixBySnrDb(clean, noise, 6);
        }

        updateSpectrograms(currentMixed, cleanAudio.sampleRate());
        updateFeaturesInfo(currentMixed, cleanAudio.sampleRate());
        if (lastMixedFile != null) {
            appendInfo("Последний микс (SNR=6) сохранен во временный файл: " + lastMixedFile.getAbsolutePath());
        }
    }

    private void updateSpectrograms(double[] signal, int sampleRate) {
        double[][] manual = Lab3AudioProcessor.spectrogramManual(signal, FFT_SIZE, HOP_SIZE);
        double[][] library = Lab3AudioProcessor.spectrogramLibrary(signal, FFT_SIZE, HOP_SIZE);

        manualSpecContainer.removeAll();
        librarySpecContainer.removeAll();
        manualSpecContainer.add(ChartUtils.createSpectrogramChartPanel(manual,
                "Самописная спектрограмма (DFT + Hann)", sampleRate, HOP_SIZE), BorderLayout.CENTER);
        librarySpecContainer.add(ChartUtils.createSpectrogramChartPanel(library,
                "Библиотечная спектрограмма (JTransforms)", sampleRate, HOP_SIZE), BorderLayout.CENTER);
        manualSpecContainer.revalidate();
        manualSpecContainer.repaint();
        librarySpecContainer.revalidate();
        librarySpecContainer.repaint();
    }

    private void updateFeaturesInfo(double[] signal, int sampleRate) {
        double[][] spec = Lab3AudioProcessor.spectrogramLibrary(signal, FFT_SIZE, HOP_SIZE);
        double centroidSum = 0.0;
        double bandwidthSum = 0.0;
        int frames = spec.length;
        for (double[] frame : spec) {
            double c = Lab3AudioProcessor.spectralCentroid(frame, sampleRate, FFT_SIZE);
            double b = Lab3AudioProcessor.spectralBandwidth(frame, sampleRate, FFT_SIZE, c);
            centroidSum += c;
            bandwidthSum += b;
        }

        double centroidAvg = frames > 0 ? centroidSum / frames : 0.0;
        double bandwidthAvg = frames > 0 ? bandwidthSum / frames : 0.0;

        infoArea.setText("""
                Вариант ЛР №3:
                - Тип представления: Спектрограмма
                - Самописные признаки: спектральный центроид, спектральная ширина
                - Самописные метрики: SNR, SI-SDR
                - Модель шумоподавления: DeepFilterNet2 (внешний CLI)
                - Диапазон SNR: 0..18 дБ, шаг 3 дБ

                Оценка признаков для текущего зашумленного сигнала:
                - Средний спектральный центроид: %s Гц
                - Средняя спектральная ширина: %s Гц

                Интерпретация:
                - Рост центроида обычно указывает на усиление высокочастотных компонентов (более "яркий" звук).
                - Рост ширины отражает более широкое распределение энергии по частотам, что характерно для шумового сигнала.
                """.formatted(formatDouble(centroidAvg), formatDouble(bandwidthAvg)));
    }

    private void runDeepFilterOnAllMixes() {
        if (cleanAudio == null || mixedCases.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Сначала выполните расчеты по диапазону SNR.");
            return;
        }
        String command = deepFilterCommandField.getText().trim();
        if (command.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Укажите команду DeepFilterNet2 CLI.");
            return;
        }

        appendInfo("Запуск DeepFilterNet2 для всех SNR. Это может занять некоторое время...");

        SwingWorker<Void, String> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() {
                double[] clean = cleanAudio.samples();

                for (MixedCase c : mixedCases) {
                    try {
                        Path outDir = Files.createTempDirectory("lab3_deepfilter_out_snr" + c.snrInputDb() + "_");
                        Path outWavPath = outDir.resolve(c.mixedFile().getName());
                        String rendered = command
                                .replace("{in}", c.mixedFile().getAbsolutePath())
                                .replace("{outDir}", outDir.toString());

                        publish("SNR вход=" + c.snrInputDb() + " дБ • запуск: " + rendered);

                        ProcessBuilder pb = new ProcessBuilder("cmd", "/c", rendered);
                        pb.environment().putIfAbsent("TORCHAUDIO_BACKEND", "soundfile");
                        pb.redirectErrorStream(true);
                        Process process = pb.start();
                        String processOutput = readProcessOutput(process.getInputStream());
                        int code = process.waitFor();
                        if (code != 0) {
                            publish("SNR вход=" + c.snrInputDb() + " • ошибка, код " + code + "\n" + processOutput);
                            continue;
                        }
                        if (!Files.exists(outWavPath)) {
                            publish("SNR вход=" + c.snrInputDb() + " • нет выходного файла: " + outWavPath);
                            if (!processOutput.isBlank()) {
                                publish(processOutput);
                            }
                            continue;
                        }

                        AudioIOUtils.AudioData denoised = AudioIOUtils.readMonoWav(outWavPath.toFile());
                        double snr = Lab3AudioProcessor.snrDb(clean, denoised.samples());
                        double siSdr = Lab3AudioProcessor.siSdrDb(clean, denoised.samples());

                        int row = c.rowIndex();
                        SwingUtilities.invokeLater(() -> {
                            metricsModel.setValueAt(formatDouble(snr), row, 3);
                            metricsModel.setValueAt(formatDouble(siSdr), row, 4);
                        });

                        publish("SNR вход=" + c.snrInputDb() + " • готово: SNR=" + formatDouble(snr) + ", SI-SDR=" + formatDouble(siSdr));
                    } catch (Exception ex) {
                        publish("SNR вход=" + c.snrInputDb() + " • исключение: " + ex.getMessage());
                    }
                }
                return null;
            }

            @Override
            protected void process(List<String> chunks) {
                for (String s : chunks) {
                    appendInfo(s);
                }
            }

            @Override
            protected void done() {
                appendInfo("DeepFilterNet2: обработка всех SNR завершена.");
            }
        };

        worker.execute();
    }

    private String readProcessOutput(InputStream stream) throws IOException {
        byte[] all = stream.readAllBytes();
        return new String(all, StandardCharsets.UTF_8).trim();
    }

    private String formatDouble(double value) {
        return String.format(Locale.US, "%.3f", value);
    }

    private void appendInfo(String text) {
        if (!infoArea.getText().isEmpty()) {
            infoArea.append("\n");
        }
        infoArea.append(text);
    }
}
