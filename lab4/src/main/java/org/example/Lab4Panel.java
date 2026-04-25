package org.example;

import org.jtransforms.fft.DoubleFFT_1D;

import javax.sound.sampled.UnsupportedAudioFileException;
import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class Lab4Panel extends JPanel {
    private static final String LAB_ROOT_DIR_NAME = "Lab4_COSI_tts_vc";
    private static final DateTimeFormatter TS_FORMAT = DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss");

    private final JTextArea ttsTextArea = new JTextArea(
            "Здравствуйте! Это демонстрация синтеза речи для лабораторной работы номер четыре.");
    private final JTextField ttsCommandField = new JTextField(
            "py -3.11 python\\lab4_tts.py --text \"{text}\" --out \"{out}\"");
    private final JTextField ttsOutField = new JTextField();

    private final JTextField srcVoiceField = new JTextField();
    private final JTextField tgtVoiceField = new JTextField();
    private final JTextField vcOutField = new JTextField();
    private final JTextField cosyVcCommandField = new JTextField(
            "py -3.11 python\\lab4_vc.py --backend keepwords --mode cosy --src \"{src}\" --tgt \"{tgt}\" --out \"{out}\"");
    private final JTextField knnVcCommandField = new JTextField(
            "py -3.11 python\\lab4_vc.py --backend keepwords --mode knn --src \"{src}\" --tgt \"{tgt}\" --out \"{out}\"");
    private final JComboBox<String> modelCombo = new JComboBox<>(new String[]{"CosyVoice3", "kNN-VC"});
    private final JTextField minLenField = new JTextField("1,2,3,5");

    private final JTextArea reportArea = new JTextArea();
    private final JTextArea logArea = new JTextArea();

    private final DefaultTableModel experimentsModel = new DefaultTableModel(
            new Object[]{
                    "Эксперимент",
                    "Модель",
                    "Длина ref, c",
                    "Время, c",
                    "Схожесть с target",
                    "Схожесть с source",
                    "Файл результата"
            }, 0);

    public Lab4Panel() {
        setLayout(new BorderLayout(8, 8));
        setBackground(ChartUtils.LIGHT_BEIGE);
        setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));

        JTabbedPane tabs = new JTabbedPane();
        tabs.setBackground(ChartUtils.OFF_WHITE);
        tabs.setFont(new Font("Arial", Font.PLAIN, 12));
        tabs.addTab("П.1 TTS", buildTtsTab());
        tabs.addTab("П.2 VC", buildVcTab());
        tabs.addTab("П.2.2–2.4 Эксперименты", buildExperimentsTab());
        tabs.addTab("Отчет", buildReportTab());

        add(tabs, BorderLayout.CENTER);
    }

    private JPanel buildTtsTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel controls = new JPanel(new GridBagLayout());
        controls.setBackground(ChartUtils.CREAM);
        controls.setBorder(BorderFactory.createTitledBorder("П.1 Озвучивание текста (CosyVoice3 или аналог)"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        ttsTextArea.setLineWrap(true);
        ttsTextArea.setWrapStyleWord(true);
        JScrollPane textScroll = new JScrollPane(ttsTextArea);
        textScroll.setPreferredSize(new Dimension(300, 120));

        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 0;
        controls.add(new JLabel("Текст:"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        controls.add(textScroll, gbc);

        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 0;
        controls.add(new JLabel("Команда TTS:"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        controls.add(ttsCommandField, gbc);

        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.weightx = 0;
        controls.add(new JLabel("Выходной WAV:"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        controls.add(ttsOutField, gbc);
        JButton outBtn = new JButton("Сохранить как…");
        outBtn.addActionListener(e -> chooseSavePath(ttsOutField));
        gbc.gridx = 2;
        gbc.weightx = 0;
        controls.add(outBtn, gbc);

        JButton runTtsBtn = new JButton("Озвучить текст");
        runTtsBtn.addActionListener(e -> runTts());
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.gridwidth = 3;
        controls.add(runTtsBtn, gbc);

        panel.add(controls, BorderLayout.NORTH);
        panel.add(new JScrollPane(createDescriptionArea(
                """
                        Что делает вкладка:
                        - Запускает внешний пайплайн TTS командой из поля "Команда TTS".
                        - Поддерживает шаблоны {text} и {out}.
                        - Результат можно использовать для оценки естественности (MOS) и корректности содержания.
                        """)), BorderLayout.CENTER);
        return panel;
    }

    private JPanel buildVcTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel controls = new JPanel(new GridBagLayout());
        controls.setBackground(ChartUtils.CREAM);
        controls.setBorder(BorderFactory.createTitledBorder("П.2 Преобразование голоса (CosyVoice3 / kNN-VC)"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        srcVoiceField.setEditable(false);
        tgtVoiceField.setEditable(false);

        int row = 0;
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.weightx = 0;
        controls.add(new JLabel("Исходный голос (source):"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        controls.add(srcVoiceField, gbc);
        JButton srcBtn = new JButton("WAV…");
        srcBtn.addActionListener(e -> chooseOpenPath(srcVoiceField));
        gbc.gridx = 2;
        gbc.weightx = 0;
        controls.add(srcBtn, gbc);

        row++;
        gbc.gridx = 0;
        gbc.gridy = row;
        controls.add(new JLabel("Целевой голос (target):"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        controls.add(tgtVoiceField, gbc);
        JButton tgtBtn = new JButton("WAV…");
        tgtBtn.addActionListener(e -> chooseOpenPath(tgtVoiceField));
        gbc.gridx = 2;
        gbc.weightx = 0;
        controls.add(tgtBtn, gbc);

        row++;
        gbc.gridx = 0;
        gbc.gridy = row;
        controls.add(new JLabel("Модель:"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        controls.add(modelCombo, gbc);

        row++;
        gbc.gridx = 0;
        gbc.gridy = row;
        controls.add(new JLabel("Команда CosyVoice VC:"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        gbc.gridwidth = 2;
        controls.add(cosyVcCommandField, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0;
        gbc.gridy = row;
        controls.add(new JLabel("Команда kNN-VC:"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        gbc.gridwidth = 2;
        controls.add(knnVcCommandField, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0;
        gbc.gridy = row;
        controls.add(new JLabel("Выходной WAV:"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        controls.add(vcOutField, gbc);
        JButton outBtn = new JButton("Сохранить как…");
        outBtn.addActionListener(e -> chooseSavePath(vcOutField));
        gbc.gridx = 2;
        gbc.weightx = 0;
        controls.add(outBtn, gbc);

        row++;
        JButton runVcBtn = new JButton("Сконвертировать голос");
        runVcBtn.addActionListener(e -> runSingleVoiceConversion());
        gbc.gridx = 0;
        gbc.gridy = row;
        gbc.gridwidth = 3;
        controls.add(runVcBtn, gbc);

        panel.add(controls, BorderLayout.NORTH);
        panel.add(new JScrollPane(createDescriptionArea(
                """
                        Что делает вкладка:
                        - Запускает внешний пайплайн VC для выбранной модели.
                        - Измеряет время выполнения и добавляет результат в таблицу экспериментов.
                        - Считает прокси-метрики:
                          * схожесть с target (приближение качества переноса тембра),
                          * схожесть с source (приближение сохранения содержания/ритма).
                        Шаблоны команд: {src}, {tgt}, {out}.
                        """)), BorderLayout.CENTER);
        return panel;
    }

    private JPanel buildExperimentsTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel top = new JPanel(new GridBagLayout());
        top.setBackground(ChartUtils.CREAM);
        top.setBorder(BorderFactory.createTitledBorder("П.2.2–2.4: ограниченные ресурсы, минимальная длина, сравнение моделей"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 0;
        top.add(new JLabel("Длины ref (сек, через запятую):"), gbc);
        gbc.gridx = 1;
        gbc.weightx = 1;
        top.add(minLenField, gbc);

        JButton minLenBtn = new JButton("Эксперимент min-length для выбранной модели");
        minLenBtn.addActionListener(e -> runMinLengthExperiment());
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        top.add(minLenBtn, gbc);

        JButton compareBtn = new JButton("Сравнить CosyVoice3 vs kNN-VC");
        compareBtn.addActionListener(e -> runCompareModelsExperiment());
        gbc.gridy = 2;
        top.add(compareBtn, gbc);

        JTable table = new JTable(experimentsModel);
        JScrollPane tableScroll = new JScrollPane(table);
        tableScroll.setBorder(BorderFactory.createTitledBorder("Результаты экспериментов"));

        logArea.setEditable(false);
        logArea.setLineWrap(true);
        logArea.setWrapStyleWord(true);
        logArea.setFont(new Font("Consolas", Font.PLAIN, 12));
        logArea.setBackground(ChartUtils.CREAM);
        JScrollPane logScroll = new JScrollPane(logArea);
        logScroll.setBorder(BorderFactory.createTitledBorder("Лог выполнения"));

        JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT, tableScroll, logScroll);
        split.setResizeWeight(0.55);
        split.setBorder(null);

        panel.add(top, BorderLayout.NORTH);
        panel.add(split, BorderLayout.CENTER);
        return panel;
    }

    private JPanel buildReportTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        reportArea.setEditable(true);
        reportArea.setLineWrap(true);
        reportArea.setWrapStyleWord(true);
        reportArea.setFont(new Font("Consolas", Font.PLAIN, 12));
        reportArea.setBackground(ChartUtils.CREAM);
        reportArea.setText("""
                Шаблон отчета ЛР №4

                1) Озвучивание текста (CosyVoice3)
                   - Входной текст:
                   - Команда/конфиг:
                   - Качество результата (MOS субъективно):

                2) Преобразование голоса (kNN-VC / CosyVoice)
                   - Source:
                   - Target:
                   - Модель:
                   - Наблюдения по качеству:

                2.1) Дообучение вокодера
                   - Среда (например, Colab):
                   - Параметры обучения:
                   - Изменение качества:

                2.2) Работа при ограниченных ресурсах
                   - Время выполнения:
                   - Ограничения:
                   - Вывод:

                2.3) Минимальная длина reference
                   - Тестовые длины:
                   - Минимально приемлемая длина:

                2.4) Сравнение CosyVoice3 и kNN-VC
                   - Сложность запуска:
                   - Качество голоса:
                   - Сохранение содержания:
                   - Финальный вывод:
                """);

        panel.add(new JScrollPane(reportArea), BorderLayout.CENTER);
        return panel;
    }

    private JTextArea createDescriptionArea(String text) {
        JTextArea area = new JTextArea(text);
        area.setEditable(false);
        area.setLineWrap(true);
        area.setWrapStyleWord(true);
        area.setFont(new Font("Arial", Font.PLAIN, 13));
        area.setBackground(ChartUtils.CREAM);
        return area;
    }

    private Path labRoot() {
        return Paths.get(System.getProperty("user.home", ".")).resolve(LAB_ROOT_DIR_NAME);
    }

    private Path outputsDir() {
        return labRoot().resolve("outputs");
    }

    private Path tempDir() {
        return labRoot().resolve("temp");
    }

    private void ensureLabDirs() throws IOException {
        Files.createDirectories(outputsDir());
        Files.createDirectories(tempDir());
    }

    private void chooseOpenPath(JTextField field) {
        JFileChooser chooser = new JFileChooser();
        if (chooser.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) return;
        field.setText(chooser.getSelectedFile().getAbsolutePath());
    }

    private void chooseSavePath(JTextField field) {
        JFileChooser chooser = new JFileChooser();
        chooser.setSelectedFile(new File("lab4_" + timestampNow() + ".wav"));
        if (chooser.showSaveDialog(this) != JFileChooser.APPROVE_OPTION) return;
        field.setText(chooser.getSelectedFile().getAbsolutePath());
    }

    private String timestampNow() {
        return LocalDateTime.now().format(TS_FORMAT);
    }

    private void runTts() {
        String text = ttsTextArea.getText().trim();
        if (text.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Введите текст для синтеза.");
            return;
        }
        String template = ttsCommandField.getText().trim();
        if (template.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Укажите команду TTS.");
            return;
        }

        String outPath = ttsOutField.getText().trim();
        if (outPath.isEmpty()) {
            try {
                ensureLabDirs();
                outPath = outputsDir().resolve("tts_" + timestampNow() + ".wav").toString();
                ttsOutField.setText(outPath);
            } catch (IOException ex) {
                showError("Не удалось создать директории лабораторной: " + ex.getMessage());
                return;
            }
        }

        String command = template
                .replace("{text}", escapeForQuotedArg(text))
                .replace("{out}", outPath);
        final String finalOutPath = outPath;

        appendLog("TTS запуск: " + command);
        runCommandAsync("TTS", command, finalOutPath, () -> {
            appendLog("TTS готово: " + finalOutPath);
            appendLog("Добавьте субъективную оценку MOS в вкладке Отчет.");
        });
    }

    private void runSingleVoiceConversion() {
        String src = srcVoiceField.getText().trim();
        String tgt = tgtVoiceField.getText().trim();
        if (src.isEmpty() || tgt.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Выберите source и target WAV.");
            return;
        }

        String model = selectedModel();
        String template = model.equals("CosyVoice3") ? cosyVcCommandField.getText().trim() : knnVcCommandField.getText().trim();
        if (template.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Не заполнена команда для выбранной модели.");
            return;
        }

        String out = vcOutField.getText().trim();
        if (out.isEmpty()) {
            try {
                ensureLabDirs();
                out = outputsDir().resolve(modelTag(model) + "_vc_" + timestampNow() + ".wav").toString();
                vcOutField.setText(out);
            } catch (IOException ex) {
                showError("Не удалось создать директории лабораторной: " + ex.getMessage());
                return;
            }
        }

        String command = renderVcCommand(template, src, tgt, out, null);
        long started = System.nanoTime();
        appendLog(model + " VC запуск: " + command);
        final String modelFinal = model;
        final String outFinal = out;
        runCommandAsync(model + " VC", command, out, () -> {
            double secs = (System.nanoTime() - started) / 1_000_000_000.0;
            addEvaluationRow("single", modelFinal, null, secs, src, tgt, outFinal);
        });
    }

    private void runMinLengthExperiment() {
        String src = srcVoiceField.getText().trim();
        String tgt = tgtVoiceField.getText().trim();
        if (src.isEmpty() || tgt.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Выберите source и target WAV.");
            return;
        }
        List<Double> lengths = parseLengthList(minLenField.getText());
        if (lengths.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Введите корректные длины, например: 1,2,3,5");
            return;
        }
        String model = selectedModel();
        String template = model.equals("CosyVoice3") ? cosyVcCommandField.getText().trim() : knnVcCommandField.getText().trim();
        if (template.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Не заполнена команда для выбранной модели.");
            return;
        }

        SwingWorker<Void, String> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() {
                try {
                    ensureLabDirs();
                    AudioIOUtils.AudioData tgtData = AudioIOUtils.readMonoWav(new File(tgt));
                    for (double lenSec : lengths) {
                        if (lenSec <= 0.0) continue;
                        String tag = modelTag(model) + "_len" + lenSec;
                        Path shortRefPath = tempDir().resolve(tag + "_ref.wav");
                        Path outPath = outputsDir().resolve(tag + "_" + timestampNow() + ".wav");

                        double[] trimmed = trimToSeconds(tgtData.samples(), tgtData.sampleRate(), lenSec);
                        AudioIOUtils.writeMonoWav(shortRefPath.toFile(), trimmed, tgtData.sampleRate());

                        String command = renderVcCommand(template, src, shortRefPath.toString(), outPath.toString(), lenSec);
                        long started = System.nanoTime();
                        publish("min-length: запуск " + command);
                        ExecResult result = runCommand(command);
                        if (result.exitCode != 0) {
                            publish("min-length: ошибка код " + result.exitCode + "\n" + result.output);
                            continue;
                        }
                        if (!Files.exists(outPath)) {
                            publish("min-length: выходной WAV не найден " + outPath);
                            continue;
                        }
                        double secs = (System.nanoTime() - started) / 1_000_000_000.0;
                        SwingUtilities.invokeLater(() ->
                                addEvaluationRow("min-length", model, lenSec, secs, src, tgt, outPath.toString()));
                    }
                } catch (Exception ex) {
                    publish("min-length: исключение: " + ex.getMessage());
                }
                return null;
            }

            @Override
            protected void process(List<String> chunks) {
                for (String s : chunks) appendLog(s);
            }
        };
        worker.execute();
    }

    private void runCompareModelsExperiment() {
        String src = srcVoiceField.getText().trim();
        String tgt = tgtVoiceField.getText().trim();
        if (src.isEmpty() || tgt.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Выберите source и target WAV.");
            return;
        }
        String cosyTemplate = cosyVcCommandField.getText().trim();
        String knnTemplate = knnVcCommandField.getText().trim();
        if (cosyTemplate.isEmpty() || knnTemplate.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Заполните команды для обеих моделей.");
            return;
        }

        SwingWorker<Void, String> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() {
                try {
                    ensureLabDirs();
                    runModelCompare("CosyVoice3", cosyTemplate, src, tgt, publishSink());
                    runModelCompare("kNN-VC", knnTemplate, src, tgt, publishSink());
                } catch (Exception ex) {
                    publish("compare: исключение: " + ex.getMessage());
                }
                return null;
            }

            private java.util.function.Consumer<String> publishSink() {
                return this::publish;
            }

            @Override
            protected void process(List<String> chunks) {
                for (String s : chunks) appendLog(s);
            }
        };
        worker.execute();
    }

    private void runModelCompare(String model, String template, String src, String tgt, java.util.function.Consumer<String> logSink) {
        Path outPath = outputsDir().resolve(modelTag(model) + "_compare_" + timestampNow() + ".wav");
        String command = renderVcCommand(template, src, tgt, outPath.toString(), null);
        long started = System.nanoTime();
        logSink.accept("compare " + model + ": " + command);
        ExecResult result = runCommand(command);
        if (result.exitCode != 0) {
            logSink.accept("compare " + model + ": ошибка код " + result.exitCode + "\n" + result.output);
            return;
        }
        if (!Files.exists(outPath)) {
            logSink.accept("compare " + model + ": нет результата " + outPath);
            return;
        }
        double secs = (System.nanoTime() - started) / 1_000_000_000.0;
        SwingUtilities.invokeLater(() ->
                addEvaluationRow("compare", model, null, secs, src, tgt, outPath.toString()));
    }

    private String renderVcCommand(String template, String src, String tgt, String out, Double lenSec) {
        String command = template
                .replace("{src}", src)
                .replace("{tgt}", tgt)
                .replace("{out}", out)
                .replace("{model}", selectedModel());
        if (lenSec != null && lenSec > 0) {
            command = command + " --len-sec " + String.format(Locale.US, "%.3f", lenSec);
        }
        return command;
    }

    private String selectedModel() {
        return (String) modelCombo.getSelectedItem();
    }

    private void addEvaluationRow(String experiment, String model, Double refLenSec, double runtimeSec,
                                  String srcPath, String tgtPath, String outPath) {
        try {
            AudioIOUtils.AudioData src = AudioIOUtils.readMonoWav(new File(srcPath));
            AudioIOUtils.AudioData tgt = AudioIOUtils.readMonoWav(new File(tgtPath));
            AudioIOUtils.AudioData out = AudioIOUtils.readMonoWav(new File(outPath));

            double simToTarget = similarityProxy(out.samples(), out.sampleRate(), tgt.samples(), tgt.sampleRate());
            double simToSource = similarityProxy(out.samples(), out.sampleRate(), src.samples(), src.sampleRate());

            experimentsModel.addRow(new Object[]{
                    experiment,
                    model,
                    refLenSec == null ? "—" : formatDouble(refLenSec),
                    formatDouble(runtimeSec),
                    formatDouble(simToTarget),
                    formatDouble(simToSource),
                    outPath
            });
            appendLog(String.format(Locale.US,
                    "%s/%s: готово. target=%.3f source=%.3f time=%.3fs",
                    experiment, model, simToTarget, simToSource, runtimeSec));
        } catch (UnsupportedAudioFileException | IOException ex) {
            appendLog("Не удалось оценить результат " + outPath + ": " + ex.getMessage());
        }
    }

    private double[] trimToSeconds(double[] data, int sampleRate, double sec) {
        int n = Math.max(1, (int) Math.round(sec * sampleRate));
        int len = Math.min(n, data.length);
        double[] out = new double[len];
        System.arraycopy(data, 0, out, 0, len);
        return out;
    }

    private List<Double> parseLengthList(String raw) {
        List<Double> out = new ArrayList<>();
        if (raw == null || raw.isBlank()) return out;
        String[] parts = raw.split(",");
        for (String p : parts) {
            try {
                out.add(Double.parseDouble(p.trim().replace(" ", "")));
            } catch (NumberFormatException ignored) {
            }
        }
        return out;
    }

    private double similarityProxy(double[] a, int srA, double[] b, int srB) {
        double[] fa = fingerprint(a, srA);
        double[] fb = fingerprint(b, srB);
        return cosineSimilarity(fa, fb);
    }

    private double[] fingerprint(double[] signal, int sampleRate) {
        int n = 1;
        int target = Math.min(signal.length, Math.max(2048, sampleRate * 2));
        while (n < target) n <<= 1;
        double[] fftIn = new double[2 * n];
        int copy = Math.min(signal.length, n);
        for (int i = 0; i < copy; i++) fftIn[i] = signal[i];
        DoubleFFT_1D fft = new DoubleFFT_1D(n);
        fft.realForwardFull(fftIn);

        int bins = 128;
        double[] fp = new double[bins];
        int usable = n / 2;
        for (int k = 1; k < usable; k++) {
            double re = fftIn[2 * k];
            double im = fftIn[2 * k + 1];
            double mag = Math.log1p(Math.sqrt(re * re + im * im));
            int idx = (int) ((k / (double) usable) * bins);
            if (idx >= bins) idx = bins - 1;
            fp[idx] += mag;
        }
        normalizeVector(fp);
        return fp;
    }

    private void normalizeVector(double[] v) {
        double norm = 0.0;
        for (double x : v) norm += x * x;
        norm = Math.sqrt(norm);
        if (norm < 1e-12) return;
        for (int i = 0; i < v.length; i++) v[i] /= norm;
    }

    private double cosineSimilarity(double[] a, double[] b) {
        int len = Math.min(a.length, b.length);
        double dot = 0.0;
        double na = 0.0;
        double nb = 0.0;
        for (int i = 0; i < len; i++) {
            dot += a[i] * b[i];
            na += a[i] * a[i];
            nb += b[i] * b[i];
        }
        if (na < 1e-12 || nb < 1e-12) return 0.0;
        return dot / Math.sqrt(na * nb);
    }

    private void runCommandAsync(String name, String command, String expectedWavPath, Runnable onSuccess) {
        SwingWorker<Void, String> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() {
                ExecResult r = runCommand(command);
                if (r.exitCode != 0) {
                    publish(name + " ошибка, код " + r.exitCode + "\n" + r.output);
                    return null;
                }
                if (expectedWavPath != null && !expectedWavPath.isBlank() && !Files.exists(Path.of(expectedWavPath))) {
                    publish(name + " завершен без выходного файла: " + expectedWavPath + "\n" + r.output);
                    return null;
                }
                SwingUtilities.invokeLater(onSuccess);
                return null;
            }

            @Override
            protected void process(List<String> chunks) {
                for (String s : chunks) appendLog(s);
            }
        };
        worker.execute();
    }

    private ExecResult runCommand(String command) {
        try {
            ProcessBuilder pb = new ProcessBuilder("cmd", "/c", command);
            pb.redirectErrorStream(true);
            Process process = pb.start();
            String out = readProcessOutput(process.getInputStream());
            int code = process.waitFor();
            return new ExecResult(code, out);
        } catch (Exception ex) {
            return new ExecResult(-1, ex.getMessage());
        }
    }

    private String readProcessOutput(InputStream stream) throws IOException {
        byte[] all = stream.readAllBytes();
        return new String(all, StandardCharsets.UTF_8).trim();
    }

    private String escapeForQuotedArg(String s) {
        return s.replace("\"", "\\\"").replace("\r", " ").replace("\n", " ");
    }

    private String formatDouble(double value) {
        return String.format(Locale.US, "%.3f", value);
    }

    private String modelTag(String model) {
        return model.toLowerCase(Locale.ROOT).contains("knn") ? "knn" : "cosy";
    }

    private void appendLog(String text) {
        if (!logArea.getText().isEmpty()) logArea.append("\n");
        logArea.append(text);
    }

    private void showError(String message) {
        JOptionPane.showMessageDialog(this, message, "Ошибка", JOptionPane.ERROR_MESSAGE);
    }

    private record ExecResult(int exitCode, String output) {}
}
