// файл: org/example/Complex.java
package org.example;

public class Complex {
    private final double re;
    private final double im;

    public Complex(double real, double imag) {
        this.re = real;
        this.im = imag;
    }

    public double re() { return re; }
    public double im() { return im; }

    public Complex add(Complex b) {
        return new Complex(this.re + b.re, this.im + b.im);
    }

    public Complex subtract(Complex b) {
        return new Complex(this.re - b.re, this.im - b.im);
    }

    public Complex multiply(Complex b) {
        return new Complex(
                this.re * b.re - this.im * b.im,
                this.re * b.im + this.im * b.re);
    }

    public Complex conjugate() {
        return new Complex(re, -im);
    }

    public double abs() {
        return Math.sqrt(re * re + im * im);
    }

    public double phase() {
        return Math.atan2(im, re);
    }
}