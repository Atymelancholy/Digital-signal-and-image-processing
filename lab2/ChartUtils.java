package org.example;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYBarDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;

public class ChartUtils {

    public static final Color LIGHT_BEIGE = new Color(250, 245, 238);
    public static final Color OFF_WHITE = new Color(252, 250, 245);
    public static final Color CREAM = new Color(255, 253, 248);
    public static final Color LIGHT_GRAY = new Color(240, 238, 235);
    public static final Color MEDIUM_GRAY = new Color(180, 175, 170);
    public static final Color DARK_BROWN = new Color(70, 50, 30);
    public static final Color WARM_GRAY = new Color(100, 90, 80);
    public static final Color FOREST_GREEN = new Color(45, 100, 45);
    public static final Color EARTH_GREEN = new Color(85, 107, 47);
    public static final Color TERRA_COTTA = new Color(204, 85, 0);
    public static final Color SLATE_BLUE = new Color(72, 61, 139);
    public static final Color DARK_SLATE = new Color(47, 79, 79);

    public static void applyChartStyle(JFreeChart chart, Color lineColor) {
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(CREAM);
        plot.setDomainGridlinePaint(new Color(220, 220, 220));
        plot.setRangeGridlinePaint(new Color(220, 220, 220));
        plot.setDomainGridlineStroke(new BasicStroke(0.5f));
        plot.setRangeGridlineStroke(new BasicStroke(0.5f));

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        renderer.setSeriesPaint(0, lineColor);
        renderer.setSeriesStroke(0, new BasicStroke(1.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        plot.setRenderer(renderer);

        plot.setOutlinePaint(MEDIUM_GRAY);
        plot.setOutlineStroke(new BasicStroke(1.0f));

        chart.getTitle().setFont(new Font("Arial", Font.BOLD, 12));
        chart.getTitle().setPaint(DARK_BROWN);

        plot.getDomainAxis().setLabelFont(new Font("Arial", Font.PLAIN, 10));
        plot.getDomainAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 9));
        plot.getDomainAxis().setLabelPaint(DARK_BROWN);
        plot.getDomainAxis().setTickLabelPaint(WARM_GRAY);

        plot.getRangeAxis().setLabelFont(new Font("Arial", Font.PLAIN, 10));
        plot.getRangeAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 9));
        plot.getRangeAxis().setLabelPaint(DARK_BROWN);
        plot.getRangeAxis().setTickLabelPaint(WARM_GRAY);

        if (chart.getLegend() != null) {
            chart.getLegend().setBackgroundPaint(null);
            chart.getLegend().setItemFont(new Font("Arial", Font.PLAIN, 8));
            chart.getLegend().setItemPaint(WARM_GRAY);
        }
    }

    public static void applyScatterStyle(JFreeChart chart, Color pointColor) {
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(CREAM);
        plot.setDomainGridlinePaint(new Color(220, 220, 220));
        plot.setRangeGridlinePaint(new Color(220, 220, 220));
        plot.setDomainGridlineStroke(new BasicStroke(0.5f));
        plot.setRangeGridlineStroke(new BasicStroke(0.5f));
        plot.setOutlinePaint(MEDIUM_GRAY);
        plot.setOutlineStroke(new BasicStroke(1.0f));

        chart.getTitle().setFont(new Font("Arial", Font.BOLD, 12));
        chart.getTitle().setPaint(DARK_BROWN);

        plot.getDomainAxis().setLabelFont(new Font("Arial", Font.PLAIN, 10));
        plot.getDomainAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 9));
        plot.getDomainAxis().setLabelPaint(DARK_BROWN);
        plot.getDomainAxis().setTickLabelPaint(WARM_GRAY);
        plot.getRangeAxis().setLabelFont(new Font("Arial", Font.PLAIN, 10));
        plot.getRangeAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 9));
        plot.getRangeAxis().setLabelPaint(DARK_BROWN);
        plot.getRangeAxis().setTickLabelPaint(WARM_GRAY);

        if (chart.getLegend() != null) {
            chart.getLegend().setBackgroundPaint(null);
            chart.getLegend().setItemFont(new Font("Arial", Font.PLAIN, 8));
            chart.getLegend().setItemPaint(WARM_GRAY);
        }

        XYLineAndShapeRenderer points = new XYLineAndShapeRenderer(false, true);
        points.setSeriesPaint(0, pointColor);
        points.setSeriesShape(0, new Ellipse2D.Double(-2, -2, 4, 4));
        plot.setRenderer(0, points);

        var base = plot.getDataset(0);
        if (base instanceof XYSeriesCollection) {
            XYSeriesCollection coll = (XYSeriesCollection) base;
            double barWidth = 0.5;
            if (coll.getSeriesCount() > 0 && coll.getSeries(0).getItemCount() >= 2) {
                double x0 = coll.getSeries(0).getX(0).doubleValue();
                double x1 = coll.getSeries(0).getX(1).doubleValue();
                barWidth = Math.max(1e-9, (x1 - x0) * 0.1);
            }
            XYBarDataset barsDataset = new XYBarDataset(coll, barWidth);
            XYBarRenderer bars = new XYBarRenderer();
            bars.setShadowVisible(false);
            bars.setBarPainter(new org.jfree.chart.renderer.xy.StandardXYBarPainter());
            bars.setSeriesPaint(0, pointColor);
            bars.setDrawBarOutline(false);
            bars.setMargin(0.0);
            plot.setDataset(1, barsDataset);
            plot.setRenderer(1, bars);
            plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
        }

        if (plot.getRangeAxis() instanceof NumberAxis) {
            ((NumberAxis) plot.getRangeAxis()).setAutoRangeIncludesZero(true);
        }
    }

    public static void applyMiniChartStyle(JFreeChart chart, Color c, boolean isPhasePlot) {
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(CREAM);
        plot.setDomainGridlinePaint(new Color(240, 240, 240));
        plot.setRangeGridlinePaint(new Color(240, 240, 240));
        plot.setDomainGridlineStroke(new BasicStroke(0.2f));
        plot.setRangeGridlineStroke(new BasicStroke(0.2f));

        XYLineAndShapeRenderer renderer = isPhasePlot ? new XYLineAndShapeRenderer(false, true)
                : new XYLineAndShapeRenderer(true, false);
        renderer.setSeriesPaint(0, c);
        if (isPhasePlot) {
            renderer.setSeriesShape(0, new Ellipse2D.Double(-2, -2, 4, 4));
        } else {
            renderer.setSeriesStroke(0, new BasicStroke(0.8f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        }
        plot.setRenderer(renderer);

        plot.setOutlinePaint(MEDIUM_GRAY);
        plot.setOutlineStroke(new BasicStroke(0.5f));

        plot.getDomainAxis().setLabelFont(new Font("Arial", Font.PLAIN, 0));
        plot.getDomainAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 0));
        plot.getDomainAxis().setLabelPaint(new Color(0, 0, 0, 0));
        plot.getDomainAxis().setTickLabelPaint(new Color(0, 0, 0, 0));
        plot.getDomainAxis().setAxisLinePaint(new Color(0, 0, 0, 0));
        plot.getDomainAxis().setTickMarkPaint(new Color(0, 0, 0, 0));

        plot.getRangeAxis().setLabelFont(new Font("Arial", Font.PLAIN, 0));
        plot.getRangeAxis().setTickLabelFont(new Font("Arial", Font.PLAIN, 0));
        plot.getRangeAxis().setLabelPaint(new Color(0, 0, 0, 0));
        plot.getRangeAxis().setTickLabelPaint(new Color(0, 0, 0, 0));
        plot.getRangeAxis().setAxisLinePaint(new Color(0, 0, 0, 0));
        plot.getRangeAxis().setTickMarkPaint(new Color(0, 0, 0, 0));
    }

    public static JPanel createDynamicChartPanel(double[] data, String title, String xLabel, String yLabel,
                                                 Color color, boolean isTimeDomain, boolean isSpectrum) {
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.setBackground(OFF_WHITE);
        mainPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(5, 5, 5, 5)));

        boolean isPhasePlot = (yLabel != null && yLabel.toLowerCase().contains("фаза"))
                || (title != null && title.toLowerCase().contains("фаз"));

        JLabel chartTitle = new JLabel(title);
        chartTitle.setFont(new Font("Arial", Font.BOLD, 11));
        chartTitle.setForeground(DARK_BROWN);
        chartTitle.setBorder(BorderFactory.createEmptyBorder(0, 5, 5, 5));

        XYSeries fullSeries = new XYSeries("Данные");
        double dt = isTimeDomain ? (1.0 / SignalData.SAMPLE_RATE * 1000) : 1.0;
        double scale = isSpectrum ? (SignalData.SAMPLE_RATE / (double) data.length) : 1.0;
        int displayLength = isSpectrum ? data.length / 2 : data.length;
        int step = isPhasePlot ? 1 : Math.max(1, displayLength / 2000);

        for (int i = 0; i < displayLength; i += step) {
            double xValue = i * (isTimeDomain ? dt : scale);
            double yValue = data[i];
            if (!Double.isNaN(yValue) && !Double.isInfinite(yValue)) {
                fullSeries.add(xValue, yValue);
            }
        }

        XYSeriesCollection fullDataset = new XYSeriesCollection();
        fullDataset.addSeries(fullSeries);

        JFreeChart fullChart = ChartFactory.createXYLineChart(title, xLabel, yLabel, fullDataset,
                PlotOrientation.VERTICAL, true, true, false);

        if (isPhasePlot) {
            applyScatterStyle(fullChart, color);
        } else {
            applyChartStyle(fullChart, color);
        }

        ChartPanel fullChartPanel = new ChartPanel(fullChart);
        fullChartPanel.setBackground(CREAM);
        fullChartPanel.setMouseWheelEnabled(true);
        fullChartPanel.setRangeZoomable(true);
        fullChartPanel.setDomainZoomable(true);
        fullChartPanel.setDisplayToolTips(true);

        JPanel controlPanel = createNavigationControlPanel(fullChartPanel);
        // Фиксируем высоту панели управления, чтобы она точно была видна
        controlPanel.setPreferredSize(new Dimension(controlPanel.getPreferredSize().width, 35));
        controlPanel.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));

        mainPanel.add(chartTitle, BorderLayout.NORTH);
        mainPanel.add(fullChartPanel, BorderLayout.CENTER);
        mainPanel.add(controlPanel, BorderLayout.SOUTH);

        return mainPanel;
    }

    public static JPanel createMiniChartPanel(double[] data, String title, String xLabel, String yLabel,
                                              Color color, boolean isTimeDomain, boolean isSpectrum) {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(OFF_WHITE);
        panel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(3, 3, 3, 3)));

        JLabel titleLabel = new JLabel(title);
        titleLabel.setFont(new Font("Arial", Font.BOLD, 8));
        titleLabel.setForeground(DARK_BROWN);
        titleLabel.setHorizontalAlignment(SwingConstants.CENTER);
        titleLabel.setBorder(BorderFactory.createEmptyBorder(0, 0, 2, 0));

        XYSeries series = new XYSeries("Данные");
        double dt = isTimeDomain ? (1.0 / SignalData.SAMPLE_RATE * 1000) : 1.0;
        double scale = isSpectrum ? (SignalData.SAMPLE_RATE / (double) data.length) : 1.0;
        int displayLength = isSpectrum ? data.length / 2 : data.length;
        boolean isPhasePlot = (yLabel != null && yLabel.toLowerCase().contains("фаза"))
                || (title != null && title.toLowerCase().contains("фаз"));
        int step = isPhasePlot ? 1 : Math.max(1, displayLength / 200);

        for (int i = 0; i < displayLength; i += step) {
            double xValue = i * (isTimeDomain ? dt : scale);
            double yValue = data[i];
            if (!Double.isNaN(yValue) && !Double.isInfinite(yValue)) {
                series.add(xValue, yValue);
            }
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createXYLineChart(null, null, null, dataset,
                PlotOrientation.VERTICAL, false, true, false);

        applyMiniChartStyle(chart, color, isPhasePlot);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setBackground(CREAM);
        chartPanel.setMouseWheelEnabled(true);
        chartPanel.setRangeZoomable(true);
        chartPanel.setDomainZoomable(true);
        chartPanel.setDisplayToolTips(true);
        chartPanel.setMinimumDrawWidth(0);
        chartPanel.setMinimumDrawHeight(0);
        chartPanel.setMaximumDrawWidth(Integer.MAX_VALUE);
        chartPanel.setMaximumDrawHeight(Integer.MAX_VALUE);

        panel.add(titleLabel, BorderLayout.NORTH);
        panel.add(chartPanel, BorderLayout.CENTER);

        panel.setPreferredSize(new Dimension(300, 180));
        panel.setMinimumSize(new Dimension(300, 180));

        return panel;
    }

    private static JPanel createNavigationControlPanel(ChartPanel chartPanel) {
        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT, 2, 0));
        panel.setBackground(LIGHT_GRAY);
        panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

        JButton zoomInBtn = createSmallButton("+", "Увеличить масштаб");
        JButton zoomOutBtn = createSmallButton("-", "Уменьшить масштаб");
        JButton resetBtn = createSmallButton("⟲", "Сбросить масштаб");
        JButton leftBtn = createSmallButton("←", "Сдвинуть влево");
        JButton rightBtn = createSmallButton("→", "Сдвинуть вправо");

        zoomInBtn.addActionListener(e -> chartPanel.zoomInBoth(1.5, 1.5));
        zoomOutBtn.addActionListener(e -> chartPanel.zoomOutBoth(1.5, 1.5));
        resetBtn.addActionListener(e -> chartPanel.restoreAutoBounds());
        leftBtn.addActionListener(e -> chartPanel.getChart().getXYPlot().getDomainAxis().pan(-0.1));
        rightBtn.addActionListener(e -> chartPanel.getChart().getXYPlot().getDomainAxis().pan(0.1));

        panel.add(new JLabel("Масштаб:"));
        panel.add(zoomInBtn);
        panel.add(zoomOutBtn);
        panel.add(resetBtn);
        panel.add(Box.createHorizontalStrut(10));
        panel.add(new JLabel("Сдвиг:"));
        panel.add(leftBtn);
        panel.add(rightBtn);

        return panel;
    }

    private static JButton createSmallButton(String text, String tooltip) {
        JButton button = new JButton(text);
        button.setFont(new Font("Arial", Font.BOLD, 10));
        button.setPreferredSize(new Dimension(30, 25));
        button.setBackground(CREAM);
        button.setForeground(DARK_BROWN);
        button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(MEDIUM_GRAY, 1),
                BorderFactory.createEmptyBorder(2, 5, 2, 5)));
        button.setToolTipText(tooltip);
        button.setFocusPainted(false);
        button.setCursor(new Cursor(Cursor.HAND_CURSOR));
        button.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseEntered(java.awt.event.MouseEvent evt) {
                button.setBackground(OFF_WHITE);
            }
            public void mouseExited(java.awt.event.MouseEvent evt) {
                button.setBackground(CREAM);
            }
        });
        return button;
    }
}