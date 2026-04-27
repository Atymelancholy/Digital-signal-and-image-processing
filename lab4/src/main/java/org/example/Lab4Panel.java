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

    /** Path to the venv python created by python\setup_cosyvoice.bat (relative to project root). */
    private static final String VENV_PYTHON = "python\\.cosyvenv\\Scripts\\python.exe";

    private static final String[] SFT_SPEAKERS = {
            "中文女", "中文男", "日语男", "粤语女", "英文女", "英文男", "韩语女"
    };

    // ---- TTS-SFT controls ----
    private final JTextArea sftTextArea = new JTextArea(
            "Hello, this is the CosyVoice synthesizer running on a fine-tuned 300M model.");
    private final JComboBox<String> sftSpeakerCombo = new JComboBox<>(SFT_SPEAKERS);
    private final JTextField sftOutField = new JTextField();
    private final JTextField sftCommandField = new JTextField(
            VENV_PYTHON + " python\\lab4_tts.py --mode sft --speaker \"{speaker}\" --text \"{text}\" --out \"{out}\"");

    // ---- TTS Zero-Shot controls ----
    private final JTextArea zsTextArea = new JTextArea("Текст, который надо озвучить голосом из эталонного аудио.");
    private final JTextField zsPromptAudioField = new JTextField();
    private final JTextArea zsPromptTextArea = new JTextArea("Текст, произнесённый в эталонном аудио.");
    private final JTextField zsOutField = new JTextField();
    private final JTextField zsCommandField = new JTextField(
            VENV_PYTHON + " python\\lab4_tts.py --mode zero_shot --text \"{text}\""
                    + " --prompt-audio \"{prompt_audio}\" --prompt-text \"{prompt_text}\" --out \"{out}\"");

    // ---- VC controls (kNN-VC) ----
    private final JTextField srcVoiceField = new JTextField();
    private final JTextField tgtVoiceField = new JTextField();
    private final JTextField vocoderCkptField = new JTextField();
    private final JTextField vcOutField = new JTextField();
    private final JTextField knnVcCommandField = new JTextField(
            VENV_PYTHON + " python\\lab4_vc.py --mode knn --src \"{src}\" --tgt \"{tgt}\" --out \"{out}\"");

    // ---- Experiments controls ----
    private final JTextField minLenField = new JTextField("1,2,3,5");
    private final JTextArea srcTextForCompare = new JTextArea("Что говорится в исходной записи (для CosyVoice zero-shot).");
    private final JTextArea tgtTextForCompare = new JTextArea("Что говорится в целевой записи (для CosyVoice zero-shot).");
    private final JTextField cosyVcCommandField = new JTextField(
            VENV_PYTHON + " python\\lab4_vc.py --mode cosy --src \"{src}\" --tgt \"{tgt}\""
                    + " --src-text \"{src_text}\" --tgt-text \"{tgt_text}\" --out \"{out}\"");

    private final JTextArea experimentFilesArea = new JTextArea();
    private final JTextArea logArea = new JTextArea();

    private final DefaultTableModel experimentsModel = new DefaultTableModel(
            new Object[]{
                    "Эксперимент",
                    "Модель",
                    "Source WAV",
                    "Ref/target для обработки",
                    "Target для метрик",
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
        tabs.addTab("П.1 TTS (CosyVoice)", buildTtsTab());
        tabs.addTab("П.2 VC (kNN-VC)", buildVcTab());
        tabs.addTab("П.2.2–2.4 Эксперименты", buildExperimentsTab());

        add(tabs, BorderLayout.CENTER);
    }

    // ===================================================================
    //  TTS tab — two sub-tabs: SFT and Zero-Shot
    // ===================================================================

    private JPanel buildTtsTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JTabbedPane sub = new JTabbedPane();
        sub.setBackground(ChartUtils.OFF_WHITE);
        sub.setFont(new Font("Arial", Font.PLAIN, 12));
        sub.addTab("SFT (фиксированные голоса)", buildSftSubTab());
        sub.addTab("Zero-Shot (клонирование)", buildZeroShotSubTab());

        panel.add(sub, BorderLayout.CENTER);
        return panel;
    }

    private JPanel buildSftSubTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel controls = new JPanel(new GridBagLayout());
        controls.setBackground(ChartUtils.CREAM);
        controls.setBorder(BorderFactory.createTitledBorder(
                "CosyVoice-300M-SFT — синтез текста готовым голосом"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        sftTextArea.setLineWrap(true);
        sftTextArea.setWrapStyleWord(true);
        JScrollPane textScroll = new JScrollPane(sftTextArea);
        textScroll.setPreferredSize(new Dimension(300, 100));

        int row = 0;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Текст:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        controls.add(textScroll, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Голос (speaker):"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(sftSpeakerCombo, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Команда TTS:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        controls.add(sftCommandField, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Выходной WAV:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(sftOutField, gbc);
        JButton outBtn = new JButton("Сохранить как…");
        outBtn.addActionListener(e -> chooseSavePath(sftOutField));
        gbc.gridx = 2; gbc.weightx = 0;
        controls.add(outBtn, gbc);

        row++;
        JButton runBtn = new JButton("Озвучить (SFT)");
        runBtn.addActionListener(e -> runTtsSft());
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 3;
        controls.add(runBtn, gbc);

        panel.add(controls, BorderLayout.NORTH);
        panel.add(new JScrollPane(createDescriptionArea("""
                Эта вкладка использует FunAudioLLM/CosyVoice-300M-SFT.
                Голоса берутся из встроенного списка диктаторов модели.
                Перед первым запуском выполните python\\setup_cosyvoice.bat
                — он скачивает модель и устанавливает CUDA-Python окружение.
                Шаблоны команды: {text}, {speaker}, {out}.
                """)), BorderLayout.CENTER);
        return panel;
    }

    private JPanel buildZeroShotSubTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel controls = new JPanel(new GridBagLayout());
        controls.setBackground(ChartUtils.CREAM);
        controls.setBorder(BorderFactory.createTitledBorder(
                "CosyVoice-300M Zero-Shot — клонирование голоса по короткому эталону"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        zsTextArea.setLineWrap(true);
        zsTextArea.setWrapStyleWord(true);
        JScrollPane synthTextScroll = new JScrollPane(zsTextArea);
        synthTextScroll.setPreferredSize(new Dimension(300, 80));

        zsPromptTextArea.setLineWrap(true);
        zsPromptTextArea.setWrapStyleWord(true);
        JScrollPane promptTextScroll = new JScrollPane(zsPromptTextArea);
        promptTextScroll.setPreferredSize(new Dimension(300, 60));

        int row = 0;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Текст для синтеза:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        controls.add(synthTextScroll, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Эталонное аудио:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(zsPromptAudioField, gbc);
        JButton promptBtn = new JButton("WAV…");
        promptBtn.addActionListener(e -> chooseOpenPath(zsPromptAudioField));
        gbc.gridx = 2; gbc.weightx = 0;
        controls.add(promptBtn, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Текст эталона:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        controls.add(promptTextScroll, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Команда TTS:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        controls.add(zsCommandField, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Выходной WAV:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(zsOutField, gbc);
        JButton outBtn = new JButton("Сохранить как…");
        outBtn.addActionListener(e -> chooseSavePath(zsOutField));
        gbc.gridx = 2; gbc.weightx = 0;
        controls.add(outBtn, gbc);

        row++;
        JButton runBtn = new JButton("Клонировать голос (Zero-Shot)");
        runBtn.addActionListener(e -> runTtsZeroShot());
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 3;
        controls.add(runBtn, gbc);

        panel.add(controls, BorderLayout.NORTH);
        panel.add(new JScrollPane(createDescriptionArea("""
                Zero-Shot режим CosyVoice-300M: модель клонирует тембр диктатора
                из эталонного аудио (5–15 секунд достаточно) и произносит
                заданный текст этим голосом. Лучше всего работает, если в поле
                "Текст эталона" указано именно то, что произнесено в WAV.
                Шаблоны команды: {text}, {prompt_audio}, {prompt_text}, {out}.
                """)), BorderLayout.CENTER);
        return panel;
    }

    // ===================================================================
    //  VC tab — kNN-VC only (with optional fine-tuned vocoder)
    // ===================================================================

    private JPanel buildVcTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel controls = new JPanel(new GridBagLayout());
        controls.setBackground(ChartUtils.CREAM);
        controls.setBorder(BorderFactory.createTitledBorder(
                "kNN-VC (bshall/knn-vc) — преобразование голоса"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        srcVoiceField.setEditable(false);
        tgtVoiceField.setEditable(false);

        int row = 0;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        controls.add(new JLabel("Source (ваш голос):"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(srcVoiceField, gbc);
        JButton srcBtn = new JButton("WAV…");
        srcBtn.addActionListener(e -> chooseOpenPath(srcVoiceField));
        gbc.gridx = 2; gbc.weightx = 0;
        controls.add(srcBtn, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row;
        controls.add(new JLabel("Target (целевой голос):"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(tgtVoiceField, gbc);
        JButton tgtBtn = new JButton("WAV…");
        tgtBtn.addActionListener(e -> chooseOpenPath(tgtVoiceField));
        gbc.gridx = 2; gbc.weightx = 0;
        controls.add(tgtBtn, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row;
        controls.add(new JLabel("Дообученный вокодер (.pt):"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(vocoderCkptField, gbc);
        JButton ckptBtn = new JButton("Файл…");
        ckptBtn.addActionListener(e -> chooseOpenPath(vocoderCkptField));
        gbc.gridx = 2; gbc.weightx = 0;
        controls.add(ckptBtn, gbc);

        row++;
        gbc.gridx = 0; gbc.gridy = row;
        controls.add(new JLabel("Команда kNN-VC:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        controls.add(knnVcCommandField, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row;
        controls.add(new JLabel("Выходной WAV:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1;
        controls.add(vcOutField, gbc);
        JButton outBtn = new JButton("Сохранить как…");
        outBtn.addActionListener(e -> chooseSavePath(vcOutField));
        gbc.gridx = 2; gbc.weightx = 0;
        controls.add(outBtn, gbc);

        row++;
        JButton runBtn = new JButton("Сконвертировать голос (kNN-VC)");
        runBtn.addActionListener(e -> runKnnVc());
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 3;
        controls.add(runBtn, gbc);

        panel.add(controls, BorderLayout.NORTH);
        panel.add(new JScrollPane(createDescriptionArea("""
                Использует официальный проект bshall/knn-vc через torch.hub:
                WavLM-Large (slой 6) кодирует source и target, kNN-сопоставление
                по топ-4 ближайших target-векторов, синтез HiFi-GAN-вокодером.
                Слова и язык source сохраняются — kNN-VC работает на акустических
                фичах, а не на фонемах.
                Поле "Дообученный вокодер" опциональное: указывайте путь к
                finetuned_model.pt из python\\lab4_finetune_colab.ipynb.
                Шаблоны команды: {src}, {tgt}, {out}, {vocoder_ckpt}.
                """)), BorderLayout.CENTER);
        return panel;
    }

    // ===================================================================
    //  Experiments tab
    // ===================================================================

    private JPanel buildExperimentsTab() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBackground(ChartUtils.LIGHT_BEIGE);

        JPanel top = new JPanel(new GridBagLayout());
        top.setBackground(ChartUtils.CREAM);
        top.setBorder(BorderFactory.createTitledBorder("П.2.2–2.4: ограниченные ресурсы, минимальная длина, сравнение моделей"));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        srcTextForCompare.setLineWrap(true);
        srcTextForCompare.setWrapStyleWord(true);
        tgtTextForCompare.setLineWrap(true);
        tgtTextForCompare.setWrapStyleWord(true);

        int row = 0;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        top.add(new JLabel("Длины ref (сек, через запятую):"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        top.add(minLenField, gbc);
        gbc.gridwidth = 1;

        row++;
        JButton minLenBtn = new JButton("Эксперимент min-length (kNN-VC)");
        minLenBtn.addActionListener(e -> runMinLengthExperiment());
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 3;
        top.add(minLenBtn, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        top.add(new JLabel("Текст source (для CosyVoice):"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        JScrollPane srcTextScroll = new JScrollPane(srcTextForCompare);
        srcTextScroll.setPreferredSize(new Dimension(300, 50));
        top.add(srcTextScroll, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row; gbc.weightx = 0;
        top.add(new JLabel("Текст target (для CosyVoice):"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        JScrollPane tgtTextScroll = new JScrollPane(tgtTextForCompare);
        tgtTextScroll.setPreferredSize(new Dimension(300, 50));
        top.add(tgtTextScroll, gbc);
        gbc.gridwidth = 1;

        row++;
        gbc.gridx = 0; gbc.gridy = row;
        top.add(new JLabel("Команда CosyVoice zero-shot VC:"), gbc);
        gbc.gridx = 1; gbc.weightx = 1; gbc.gridwidth = 2;
        top.add(cosyVcCommandField, gbc);
        gbc.gridwidth = 1;

        row++;
        JButton compareBtn = new JButton("Сравнить CosyVoice zero-shot vs kNN-VC");
        compareBtn.addActionListener(e -> runCompareModelsExperiment());
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 3;
        top.add(compareBtn, gbc);

        experimentFilesArea.setEditable(false);
        experimentFilesArea.setLineWrap(true);
        experimentFilesArea.setWrapStyleWord(true);
        experimentFilesArea.setFont(new Font("Consolas", Font.PLAIN, 12));
        experimentFilesArea.setBackground(ChartUtils.CREAM);
        experimentFilesArea.setBorder(BorderFactory.createTitledBorder("Файлы, которые участвуют в обработке и расчете метрик"));
        refreshExperimentFilesInfo();

        JTable table = new JTable(experimentsModel);
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.getColumnModel().getColumn(0).setPreferredWidth(110);
        table.getColumnModel().getColumn(1).setPreferredWidth(150);
        table.getColumnModel().getColumn(2).setPreferredWidth(260);
        table.getColumnModel().getColumn(3).setPreferredWidth(260);
        table.getColumnModel().getColumn(4).setPreferredWidth(260);
        table.getColumnModel().getColumn(5).setPreferredWidth(90);
        table.getColumnModel().getColumn(6).setPreferredWidth(80);
        table.getColumnModel().getColumn(7).setPreferredWidth(120);
        table.getColumnModel().getColumn(8).setPreferredWidth(120);
        table.getColumnModel().getColumn(9).setPreferredWidth(300);
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
        panel.add(experimentFilesArea, BorderLayout.SOUTH);
        panel.add(split, BorderLayout.CENTER);
        return panel;
    }

    // ===================================================================
    //  Helpers
    // ===================================================================

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

    private Path outputsDir() { return labRoot().resolve("outputs"); }
    private Path tempDir() { return labRoot().resolve("temp"); }

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

    private void refreshExperimentFilesInfo() {
        String src = srcVoiceField.getText().trim();
        String tgt = tgtVoiceField.getText().trim();
        experimentFilesArea.setText("""
                Source WAV для обработки и метрики "схожесть с source":
                %s

                Target WAV по умолчанию для обработки и метрики "схожесть с target":
                %s

                В эксперименте min-length для обработки создается отдельный короткий ref-файл,
                а метрики считаются относительно полного target WAV.
                """.formatted(displayPath(src), displayPath(tgt)));
    }

    private String displayPath(String path) {
        return path == null || path.isBlank() ? "не выбран" : path;
    }

    // ===================================================================
    //  TTS handlers
    // ===================================================================

    private void runTtsSft() {
        String text = sftTextArea.getText().trim();
        if (text.isEmpty()) { showError("Введите текст для синтеза."); return; }
        String speaker = (String) sftSpeakerCombo.getSelectedItem();
        if (speaker == null || speaker.isBlank()) { showError("Выберите голос."); return; }

        String outPath = sftOutField.getText().trim();
        if (outPath.isEmpty()) {
            try {
                ensureLabDirs();
                outPath = outputsDir().resolve("tts_sft_" + timestampNow() + ".wav").toString();
                sftOutField.setText(outPath);
            } catch (IOException ex) {
                showError("Не удалось создать директорию: " + ex.getMessage());
                return;
            }
        }

        String command = sftCommandField.getText().trim()
                .replace("{text}", escapeForQuotedArg(text))
                .replace("{speaker}", speaker)
                .replace("{out}", outPath);

        final String finalOutPath = outPath;
        appendLog("TTS-SFT запуск: " + command);
        runCommandAsync("TTS-SFT", command, finalOutPath, () -> {
            appendLog("TTS-SFT готово: " + finalOutPath);
            appendLog("Субъективную оценку MOS можно занести в свой отчет по результату.");
        });
    }

    private void runTtsZeroShot() {
        String text = zsTextArea.getText().trim();
        String promptAudio = zsPromptAudioField.getText().trim();
        String promptText = zsPromptTextArea.getText().trim();
        if (text.isEmpty() || promptAudio.isEmpty() || promptText.isEmpty()) {
            showError("Заполните все три поля: текст для синтеза, эталонное аудио, текст эталона.");
            return;
        }

        String outPath = zsOutField.getText().trim();
        if (outPath.isEmpty()) {
            try {
                ensureLabDirs();
                outPath = outputsDir().resolve("tts_zs_" + timestampNow() + ".wav").toString();
                zsOutField.setText(outPath);
            } catch (IOException ex) {
                showError("Не удалось создать директорию: " + ex.getMessage());
                return;
            }
        }

        String command = zsCommandField.getText().trim()
                .replace("{text}", escapeForQuotedArg(text))
                .replace("{prompt_audio}", promptAudio)
                .replace("{prompt_text}", escapeForQuotedArg(promptText))
                .replace("{out}", outPath);

        final String finalOutPath = outPath;
        appendLog("TTS Zero-Shot запуск: " + command);
        runCommandAsync("TTS Zero-Shot", command, finalOutPath, () -> {
            appendLog("TTS Zero-Shot готово: " + finalOutPath);
            appendLog("Субъективную оценку MOS можно занести в свой отчет по результату.");
        });
    }

    // ===================================================================
    //  VC handlers
    // ===================================================================

    private void runKnnVc() {
        String src = srcVoiceField.getText().trim();
        String tgt = tgtVoiceField.getText().trim();
        if (src.isEmpty() || tgt.isEmpty()) {
            showError("Выберите source и target WAV.");
            return;
        }

        String out = vcOutField.getText().trim();
        if (out.isEmpty()) {
            try {
                ensureLabDirs();
                out = outputsDir().resolve("knn_vc_" + timestampNow() + ".wav").toString();
                vcOutField.setText(out);
            } catch (IOException ex) {
                showError("Не удалось создать директорию: " + ex.getMessage());
                return;
            }
        }

        String command = renderKnnCommand(knnVcCommandField.getText().trim(), src, tgt, out, null);
        long started = System.nanoTime();
        appendLog("kNN-VC запуск: " + command);
        final String outFinal = out;
        runCommandAsync("kNN-VC", command, out, () -> {
            double secs = (System.nanoTime() - started) / 1_000_000_000.0;
            addEvaluationRow("single", "kNN-VC", src, tgt, tgt, null, secs, outFinal);
        });
    }

    private String renderKnnCommand(String template, String src, String tgt, String out, Double lenSec) {
        String ckpt = vocoderCkptField.getText().trim();
        String command = template
                .replace("{src}", src)
                .replace("{tgt}", tgt)
                .replace("{out}", out);
        if (!command.contains("--vocoder-ckpt") && !ckpt.isEmpty()) {
            command = command + " --vocoder-ckpt \"" + ckpt + "\"";
        } else {
            command = command.replace("{vocoder_ckpt}", ckpt);
        }
        if (lenSec != null && lenSec > 0) {
            command = command + " --len-sec " + String.format(Locale.US, "%.3f", lenSec);
        }
        return command;
    }

    // ===================================================================
    //  Experiments
    // ===================================================================

    private void runMinLengthExperiment() {
        String src = srcVoiceField.getText().trim();
        String tgt = tgtVoiceField.getText().trim();
        if (src.isEmpty() || tgt.isEmpty()) {
            showError("Выберите source и target WAV на вкладке П.2 VC.");
            return;
        }
        List<Double> lengths = parseLengthList(minLenField.getText());
        if (lengths.isEmpty()) {
            showError("Введите корректные длины, например: 1,2,3,5");
            return;
        }
        refreshExperimentFilesInfo();
        String template = knnVcCommandField.getText().trim();

        SwingWorker<Void, String> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() {
                try {
                    ensureLabDirs();
                    AudioIOUtils.AudioData tgtData = AudioIOUtils.readMonoWav(new File(tgt));
                    for (double lenSec : lengths) {
                        if (lenSec <= 0.0) continue;
                        String tag = "knn_len" + lenSec;
                        Path shortRefPath = tempDir().resolve(tag + "_ref.wav");
                        Path outPath = outputsDir().resolve(tag + "_" + timestampNow() + ".wav");

                        double[] trimmed = trimToSeconds(tgtData.samples(), tgtData.sampleRate(), lenSec);
                        AudioIOUtils.writeMonoWav(shortRefPath.toFile(), trimmed, tgtData.sampleRate());

                        String command = renderKnnCommand(template, src, shortRefPath.toString(),
                                outPath.toString(), lenSec);
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
                                addEvaluationRow("min-length", "kNN-VC", src, shortRefPath.toString(), tgt,
                                        lenSec, secs, outPath.toString()));
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
            showError("Выберите source и target WAV на вкладке П.2 VC.");
            return;
        }
        String srcText = srcTextForCompare.getText().trim();
        String tgtText = tgtTextForCompare.getText().trim();
        if (srcText.isEmpty() || tgtText.isEmpty()) {
            showError("Введите тексты source и target — они нужны CosyVoice zero-shot.");
            return;
        }
        refreshExperimentFilesInfo();

        SwingWorker<Void, String> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() {
                try {
                    ensureLabDirs();
                    Path knnOut = outputsDir().resolve("cmp_knn_" + timestampNow() + ".wav");
                    runOneCompare("kNN-VC",
                            renderKnnCommand(knnVcCommandField.getText().trim(), src, tgt,
                                    knnOut.toString(), null),
                            knnOut,
                            src, tgt, this::publish);

                    String cosyOut = outputsDir().resolve("cmp_cosy_" + timestampNow() + ".wav").toString();
                    String cosyCommand = cosyVcCommandField.getText().trim()
                            .replace("{src}", src)
                            .replace("{tgt}", tgt)
                            .replace("{src_text}", escapeForQuotedArg(srcText))
                            .replace("{tgt_text}", escapeForQuotedArg(tgtText))
                            .replace("{out}", cosyOut);
                    runOneCompare("CosyVoice zero-shot", cosyCommand, Path.of(cosyOut), src, tgt, this::publish);
                } catch (Exception ex) {
                    publish("compare: исключение: " + ex.getMessage());
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

    private void runOneCompare(String model, String command, Path outPath,
                               String src, String tgt,
                               java.util.function.Consumer<String> logSink) {
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
                addEvaluationRow("compare", model, src, tgt, tgt, null, secs, outPath.toString()));
    }

    // ===================================================================
    //  Evaluation / metrics
    // ===================================================================

    private void addEvaluationRow(String experiment, String model, String srcPath, String processingRefPath,
                                  String metricsTargetPath, Double refLenSec, double runtimeSec, String outPath) {
        try {
            AudioIOUtils.AudioData src = AudioIOUtils.readMonoWav(new File(srcPath));
            AudioIOUtils.AudioData tgt = AudioIOUtils.readMonoWav(new File(metricsTargetPath));
            AudioIOUtils.AudioData out = AudioIOUtils.readMonoWav(new File(outPath));

            double simToTarget = similarityProxy(out.samples(), out.sampleRate(), tgt.samples(), tgt.sampleRate());
            double simToSource = similarityProxy(out.samples(), out.sampleRate(), src.samples(), src.sampleRate());

            experimentsModel.addRow(new Object[]{
                    experiment,
                    model,
                    srcPath,
                    processingRefPath,
                    metricsTargetPath,
                    refLenSec == null ? "—" : formatDouble(refLenSec),
                    formatDouble(runtimeSec),
                    formatDouble(simToTarget),
                    formatDouble(simToSource),
                    outPath
            });
            appendLog(String.format(Locale.US,
                    "%s/%s: готово. обработка ref=%s; метрики target=%s; target=%.3f source=%.3f time=%.3fs",
                    experiment, model, processingRefPath, metricsTargetPath, simToTarget, simToSource, runtimeSec));
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

    // ===================================================================
    //  Process execution
    // ===================================================================

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

    private void appendLog(String text) {
        if (!logArea.getText().isEmpty()) logArea.append("\n");
        logArea.append(text);
    }

    private void showError(String message) {
        JOptionPane.showMessageDialog(this, message, "Ошибка", JOptionPane.ERROR_MESSAGE);
    }

    private record ExecResult(int exitCode, String output) {}
}
