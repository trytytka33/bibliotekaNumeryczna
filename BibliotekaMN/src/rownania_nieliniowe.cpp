//
// Created by Szymon Ros on 11/06/2025.
//
#include "../include/rownania_nieliniowe.h"
#include <iostream>
#include <cmath>
using namespace std;

namespace biblioteka_numeryczna {
    const double TOLERANCJA = 0.00000001;
const int MAX_ITER = 500;

double f1(double x) {
    if (x >= 1) return NAN;
    return log(1-x) + 1.0/(x*x + 3);
}

double pochodnaF1(double x) {
    if (x >= 1) return NAN;
    return -1.0/(1-x) - 2*x/((x*x + 3)*(x*x + 3));
}

double f2(double x) {
    if (abs(x) < 0.0000000001) return NAN;
    return x*x*x + 30*sin(x) - 12.0/x - 28;
}

double pochodnaF2(double x) {
    if (abs(x) < 0.0000000001) return NAN;
    return 3*x*x + 30*cos(x) + 12.0/(x*x);
}

double f3(double x) {
    if (abs(x+2) < 0.0000000001 || abs(x+4) < 0.0000000001) return NAN;
    return cos(3*M_PI*x)/(x+2) - 1.0/(x+4);
}

double pochodnaF3(double x) {
    if (abs(x+2) < 0.0000000001 || abs(x+4) < 0.0000000001) return NAN;
    double czlon1 = -3*M_PI*sin(3*M_PI*x)/(x+2);
    double czlon2 = -cos(3*M_PI*x)/((x+2)*(x+2));
    double czlon3 = 1.0/((x+4)*(x+4));
    return czlon1 + czlon2 + czlon3;
}

double metodaBisekcji(double a, double b, double (*funkcja)(double)) {
    double fa = funkcja(a);
    double fb = funkcja(b);

    if (isnan(fa) || isnan(fb) || fa * fb > 0) {
        return NAN;
    }

    double c;
    for (int i = 0; i < MAX_ITER; i++) {
        c = (a + b) / 2.0;
        double fc = funkcja(c);

        if (isnan(fc)) return NAN;

        if (abs(fc) < TOLERANCJA || abs(b - a) < TOLERANCJA) {
            return c;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return c;
}

double metodaNewtona(double x0, double (*funkcja)(double), double (*pochodna)(double)) {
    double x = x0;

    for (int i = 0; i < MAX_ITER; i++) {
        double fx = funkcja(x);
        double dfx = pochodna(x);

        if (isnan(fx) || isnan(dfx) || abs(dfx) < 0.000000000000001) {
            return NAN;
        }

        if (abs(fx) < TOLERANCJA) {
            return x;
        }

        double dx = fx / dfx;
        x = x - dx;

        if (abs(dx) < TOLERANCJA) {
            return x;
        }
    }
    return x;
}

double metodaSiecznych(double x0, double x1, double (*funkcja)(double)) {
    double f0 = funkcja(x0);
    double f1 = funkcja(x1);

    if (isnan(f0) || isnan(f1)) {
        return NAN;
    }

    for (int i = 0; i < MAX_ITER; i++) {
        if (abs(f1) < TOLERANCJA) {
            return x1;
        }

        double mianownik = f1 - f0;
        if (abs(mianownik) < 0.000000000000001) {
            return NAN;
        }

        double x2 = x1 - f1 * (x1 - x0) / mianownik;
        double f2 = funkcja(x2);

        if (isnan(f2)) {
            return NAN;
        }

        if (abs(f2) < TOLERANCJA) {
            return x2;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }
    return x1;
}

vector<double> znajdzWszystkiePierwiastki(double a, double b, double (*funkcja)(double), double (*pochodna)(double), double krok = 0.1) {
    vector<double> pierwiastki;

    for (double x = a; x < b - krok; x += krok) {
        double f_x = funkcja(x);
        double f_x_krok = funkcja(x + krok);

        if (!isnan(f_x) && !isnan(f_x_krok) && f_x * f_x_krok < 0) {
            double pierwiastek = metodaBisekcji(x, x + krok, funkcja);
            if (!isnan(pierwiastek)) {
                bool juzJest = false;
                for (double p : pierwiastki) {
                    if (abs(p - pierwiastek) < 0.001) {
                        juzJest = true;
                        break;
                    }
                }
                if (!juzJest) {
                    pierwiastki.push_back(pierwiastek);
                }
            }
        }
    }


    for (double x = a; x <= b; x += krok * 2) {
        double pierwiastek = metodaNewtona(x, funkcja, pochodna);
        if (!isnan(pierwiastek) && pierwiastek >= a && pierwiastek <= b) {
            bool juzJest = false;
            for (double p : pierwiastki) {
                if (abs(p - pierwiastek) < 0.001) {
                    juzJest = true;
                    break;
                }
            }
            if (!juzJest) {
                pierwiastki.push_back(pierwiastek);
            }
        }
    }

    for (double x1 = a; x1 < b - krok; x1 += krok * 3) {
        for (double x2 = x1 + krok; x2 <= b; x2 += krok * 3) {
            double pierwiastek = metodaSiecznych(x1, x2, funkcja);
            if (!isnan(pierwiastek) && pierwiastek >= a && pierwiastek <= b) {
                bool juzJest = false;
                for (double p : pierwiastki) {
                    if (abs(p - pierwiastek) < 0.001) {
                        juzJest = true;
                        break;
                    }
                }
                if (!juzJest) {
                    pierwiastki.push_back(pierwiastek);
                }
            }
        }
    }

    return pierwiastki;
}

void wypiszWynik(const string& nazwaMetody, double wynik, double (*funkcja)(double)) {
    cout << "  " << nazwaMetody << ": x = " << wynik;
    if (!isnan(wynik)) {
        cout << ", błąd bezwzględny = " << abs(funkcja(wynik));
    }
    cout << endl;
}
}