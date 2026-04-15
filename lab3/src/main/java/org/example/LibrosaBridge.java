package org.example;

import com.google.gson.Gson;
import com.google.gson.JsonSyntaxException;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public final class LibrosaBridge {
    private static final Gson GSON = new Gson();

    private LibrosaBridge() {}

    public record Result(
            int sampleRate,
            int hopLength,
            double[][] melDb, 
            double centroidHz,
            double rolloffHz,
            double bandwidthHz,
            double zcr,
            double[] centroidSeriesHz,
            double[] rolloffSeriesHz,
            double[] bandwidthSeriesHz,
            double[] zcrSeries,
            String error
    ) {}

    public static Path scriptPath() {
        // Проектный путь: <project>/python/lab3_librosa_analyze.py
        return Paths.get(System.getProperty("user.dir", "."))
                .resolve("python")
                .resolve("lab3_librosa_analyze.py");
    }

    public static Result analyzeWavWithLibrosa(String wavPath, Path cacheJsonPath) {
        try {
            Path script = scriptPath();
            if (!Files.exists(script)) {
                return new Result(0, 0, null, 0, 0, 0, 0, null, null, null, null,
                        "Не найден Python-скрипт: " + script.toAbsolutePath());
            }
            Files.createDirectories(cacheJsonPath.getParent());

            ProcessBuilder pb = new ProcessBuilder(
                    "py", "-3.11",
                    script.toAbsolutePath().toString(),
                    "--in", wavPath,
                    "--out", cacheJsonPath.toAbsolutePath().toString()
            );
            pb.redirectErrorStream(true);
            Process p = pb.start();
            String out = readAll(p.getInputStream());
            int code = p.waitFor();
            if (code != 0) {
                return new Result(0, 0, null, 0, 0, 0, 0, null, null, null, null,
                        "librosa-скрипт завершился с кодом " + code + "\n" + out);
            }
            if (!Files.exists(cacheJsonPath)) {
                return new Result(0, 0, null, 0, 0, 0, 0, null, null, null, null,
                        "Нет JSON-результата: " + cacheJsonPath.toAbsolutePath() + "\n" + out);
            }
            String json = Files.readString(cacheJsonPath, StandardCharsets.UTF_8);
            Payload payload = GSON.fromJson(json, Payload.class);
            if (payload == null || payload.mel_db == null) {
                return new Result(0, 0, null, 0, 0, 0, 0, null, null, null, null,
                        "Пустой/невалидный JSON из librosa.");
            }
            return new Result(
                    payload.sample_rate,
                    payload.hop_length,
                    payload.mel_db,
                    payload.centroid_hz,
                    payload.rolloff_hz,
                    payload.bandwidth_hz,
                    payload.zcr,
                    payload.centroid_series_hz,
                    payload.rolloff_series_hz,
                    payload.bandwidth_series_hz,
                    payload.zcr_series,
                    payload.error
            );
        } catch (InterruptedException ie) {
            Thread.currentThread().interrupt();
            return new Result(0, 0, null, 0, 0, 0, 0, null, null, null, null, "Прервано: " + ie.getMessage());
        } catch (JsonSyntaxException jse) {
            return new Result(0, 0, null, 0, 0, 0, 0, null, null, null, null, "Ошибка JSON: " + jse.getMessage());
        } catch (Exception ex) {
            return new Result(0, 0, null, 0, 0, 0, 0, null, null, null, null, "Ошибка librosa-анализа: " + ex.getMessage());
        }
    }

    private static String readAll(InputStream is) throws IOException {
        byte[] b = is.readAllBytes();
        return new String(b, StandardCharsets.UTF_8).trim();
    }

    // Gson DTO
    private static final class Payload {
        int sample_rate;
        int hop_length;
        double[][] mel_db;
        double centroid_hz;
        double rolloff_hz;
        double bandwidth_hz;
        double zcr;
        double[] centroid_series_hz;
        double[] rolloff_series_hz;
        double[] bandwidth_series_hz;
        double[] zcr_series;
        String error;
    }
}

