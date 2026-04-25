package org.example;

import com.google.gson.Gson;
import com.google.gson.JsonSyntaxException;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public final class AudioMetricsBridge {
    private static final Gson GSON = new Gson();

    private AudioMetricsBridge() {}

    public record Result(
            Double snrDb,
            Double sdrDb,
            Double siSdrDb,
            Double pesqMos,
            String error
    ) {}

    public static Path scriptPath() {
        return Paths.get(System.getProperty("user.dir", "."))
                .resolve("python")
                .resolve("lab3_audio_metrics.py");
    }

    public static Result evaluate(String refWavPath, String degWavPath, Path outJsonPath) {
        try {
            Path script = scriptPath();
            if (!Files.exists(script)) {
                return new Result(null, null, null, null, "Не найден Python-скрипт: " + script.toAbsolutePath());
            }
            Files.createDirectories(outJsonPath.getParent());

            ProcessBuilder pb = new ProcessBuilder(
                    "py", "-3.11",
                    script.toAbsolutePath().toString(),
                    "--ref", refWavPath,
                    "--deg", degWavPath,
                    "--out", outJsonPath.toAbsolutePath().toString()
            );
            pb.redirectErrorStream(true);
            Process p = pb.start();
            String out = readAll(p.getInputStream());
            int code = p.waitFor();
            if (code != 0) {
                return new Result(null, null, null, null, "Метрики: код " + code + "\n" + out);
            }
            if (!Files.exists(outJsonPath)) {
                return new Result(null, null, null, null, "Нет JSON результата: " + outJsonPath + "\n" + out);
            }
            String json = Files.readString(outJsonPath, StandardCharsets.UTF_8);
            Payload payload = GSON.fromJson(json, Payload.class);
            if (payload == null) {
                return new Result(null, null, null, null, "Пустой JSON результата метрик.");
            }
            return new Result(payload.snr_db, payload.sdr_db, payload.si_sdr_db, payload.pesq_mos, payload.error);
        } catch (InterruptedException ie) {
            Thread.currentThread().interrupt();
            return new Result(null, null, null, null, "Прервано: " + ie.getMessage());
        } catch (JsonSyntaxException jse) {
            return new Result(null, null, null, null, "Ошибка JSON: " + jse.getMessage());
        } catch (Exception ex) {
            return new Result(null, null, null, null, "Ошибка метрик: " + ex.getMessage());
        }
    }

    private static String readAll(InputStream is) throws IOException {
        byte[] b = is.readAllBytes();
        return new String(b, StandardCharsets.UTF_8).trim();
    }

    private static final class Payload {
        Double snr_db;
        Double sdr_db;
        Double si_sdr_db;
        Double pesq_mos;
        String error;
    }
}

