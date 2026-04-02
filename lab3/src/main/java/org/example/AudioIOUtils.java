package org.example;

import javax.sound.sampled.*;
import java.io.*;

public class AudioIOUtils {

    public record AudioData(double[] samples, int sampleRate) {}

    public static AudioData readMonoWav(File file) throws UnsupportedAudioFileException, IOException {
        try (AudioInputStream input = AudioSystem.getAudioInputStream(file)) {
            AudioFormat base = input.getFormat();
            AudioFormat target = new AudioFormat(
                    AudioFormat.Encoding.PCM_SIGNED,
                    base.getSampleRate(),
                    16,
                    1,
                    2,
                    base.getSampleRate(),
                    false
            );
            try (AudioInputStream pcmStream = AudioSystem.getAudioInputStream(target, input)) {
                byte[] raw = readAllBytes(pcmStream);
                int samplesCount = raw.length / 2;
                double[] samples = new double[samplesCount];
                for (int i = 0; i < samplesCount; i++) {
                    int lo = raw[2 * i] & 0xff;
                    int hi = raw[2 * i + 1];
                    short s = (short) ((hi << 8) | lo);
                    samples[i] = s / 32768.0;
                }
                return new AudioData(samples, (int) target.getSampleRate());
            }
        }
    }

    public static void writeMonoWav(File file, double[] signal, int sampleRate) throws IOException {
        byte[] pcm = new byte[signal.length * 2];
        for (int i = 0; i < signal.length; i++) {
            double v = Math.max(-1.0, Math.min(1.0, signal[i]));
            short s = (short) Math.round(v * 32767.0);
            pcm[2 * i] = (byte) (s & 0xff);
            pcm[2 * i + 1] = (byte) ((s >> 8) & 0xff);
        }
        AudioFormat format = new AudioFormat(sampleRate, 16, 1, true, false);
        try (ByteArrayInputStream bais = new ByteArrayInputStream(pcm);
             AudioInputStream ais = new AudioInputStream(bais, format, signal.length)) {
            AudioSystem.write(ais, AudioFileFormat.Type.WAVE, file);
        }
    }

    private static byte[] readAllBytes(AudioInputStream stream) throws IOException {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        byte[] buffer = new byte[8192];
        int read;
        while ((read = stream.read(buffer)) != -1) {
            baos.write(buffer, 0, read);
        }
        return baos.toByteArray();
    }
}
