//
// Created by Szymon Ros on 11/06/2025.
//
#include "../include/calkowanie_numeryczne.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
using namespace std;
using namespace std::chrono;

namespace biblioteka_numeryczna {
    double FunkcjaTrygonometryczna(double x) {
        return x * pow(cos(x), 3);
    }


    double Horner(const vector<double>& wspolczynniki, double x) {
        double wynik = 0.0;
        for (double wsp : wspolczynniki) {
            wynik = wynik * x + wsp;
        }
        return wynik;
    }


    double MetodaProstokatow(const vector<double>& wspolczynniki, double poczatekPrzedzialu, double koniecPrzedzialu, int n) {
        double dx = (koniecPrzedzialu - poczatekPrzedzialu) / n;
        double wynik = 0.0;
        for (int i = 0; i < n; i++) {
            double x = poczatekPrzedzialu + i * dx;
            wynik += Horner(wspolczynniki, x) * dx;
        }
        return wynik;
    }

    double MetodaTrapezow(const vector<double>& wspolczynniki, double poczatekPrzedzialu, double koniecPrzedzialu, int n) {
        double h = (koniecPrzedzialu - poczatekPrzedzialu) / n;
        double wynik = 0.0;
        for (int i = 0; i < n; i++) {
            double x1 = poczatekPrzedzialu + i * h;
            double x2 = poczatekPrzedzialu + (i + 1) * h;
            double f1 = Horner(wspolczynniki, x1);
            double f2 = Horner(wspolczynniki, x2);
            wynik += (f1 + f2) * h / 2.0;
        }
        return wynik;
    }


    double Simpson(const vector<double>& wspolczynniki, double a, double b, int n) {
        if (n % 2 != 0) {
            n += 1;
        }

        double dx = (b - a) / n;
        double wynik = Horner(wspolczynniki, a) + Horner(wspolczynniki, b);

        for (int i = 1; i < n; i++) {
            double x = a + i * dx;
            if (i % 2 == 0) {
                wynik += 2 * Horner(wspolczynniki, x);
            } else {
                wynik += 4 * Horner(wspolczynniki, x);
            }
        }

        wynik *= dx / 3.0;
        return wynik;
    }


    double MetodaProstokatowTryg(double poczatekPrzedzialu, double koniecPrzedzialu, int n) {
        double dx = (koniecPrzedzialu - poczatekPrzedzialu) / n;
        double wynik = 0.0;
        for (int i = 0; i < n; i++) {
            double x = poczatekPrzedzialu + i * dx;
            wynik += FunkcjaTrygonometryczna(x) * dx;
        }
        return wynik;
    }

    // dla funkcji trygonometrycznej
    double MetodaTrapezowTryg(double poczatekPrzedzialu, double koniecPrzedzialu, int n) {
        double h = (koniecPrzedzialu - poczatekPrzedzialu) / n;
        double wynik = 0.0;
        for (int i = 0; i < n; i++) {
            double x1 = poczatekPrzedzialu + i * h;
            double x2 = poczatekPrzedzialu + (i + 1) * h;
            double f1 = FunkcjaTrygonometryczna(x1);
            double f2 = FunkcjaTrygonometryczna(x2);
            wynik += (f1 + f2) * h / 2.0;
        }
        return wynik;
    }

    //  dla funkcji trygonometrycznej
    double SimpsonTryg(double a, double b, int n) {
        if (n % 2 != 0) {
            n += 1;
        }

        double dx = (b - a) / n;
        double wynik = FunkcjaTrygonometryczna(a) + FunkcjaTrygonometryczna(b);

        for (int i = 1; i < n; i++) {
            double x = a + i * dx;
            if (i % 2 == 0) {
                wynik += 2 * FunkcjaTrygonometryczna(x);
            } else {
                wynik += 4 * FunkcjaTrygonometryczna(x);
            }
        }

        wynik *= dx / 3.0;
        return wynik;
    }


    double CalkaAnalitycznaWielomianu(const vector<double>& wspolczynniki, double a, double b) {
        vector<double> calka_wspolczynniki(wspolczynniki.size() + 1, 0.0);


        for (size_t i = 0; i < wspolczynniki.size(); i++) {
            calka_wspolczynniki[i + 1] = wspolczynniki[i] / (i + 1);
        }

        //  F(b) - F(a)
        double wynik_b = 0.0;
        double wynik_a = 0.0;
        double x_pow_b = 1.0;
        double x_pow_a = 1.0;

        for (size_t i = 0; i < calka_wspolczynniki.size(); i++) {
            wynik_b += calka_wspolczynniki[i] * x_pow_b;
            wynik_a += calka_wspolczynniki[i] * x_pow_a;
            x_pow_b *= b;
            x_pow_a *= a;
        }

        return wynik_b - wynik_a;
    }
    void WyswietlZbieznosc(const vector<double>& wspolczynniki, double a, double b) {
        vector<int> liczbyPodzialow = {10, 50, 100, 500, 1000, 5000, 10000};
        double dokladnyWynik = CalkaAnalitycznaWielomianu(wspolczynniki, a, b);

        cout << "\nTest zbieznosci dla wielomianu:" << endl;
        cout << setw(15) << "Liczba podzialow" << setw(20) << "Metoda prostokatow"
             << setw(20) << "Metoda trapezow" << setw(20) << "Metoda Simpsona" << endl;

        for (int n : liczbyPodzialow) {
            double wynikProstokatow = MetodaProstokatow(wspolczynniki, a, b, n);
            double wynikTrapezow = MetodaTrapezow(wspolczynniki, a, b, n);
            double wynikSimpson = Simpson(wspolczynniki, a, b, n);

            cout << setw(15) << n
                 << setw(20) << fixed << setprecision(10) << wynikProstokatow
                 << setw(20) << wynikTrapezow
                 << setw(20) << wynikSimpson << endl;
        }
    }

    void WyswietlZbieznoscTryg(double a, double b) {
        vector<int> liczbyPodzialow = {10, 50, 100, 500, 1000, 5000, 10000};
        double dokladnyWynik = -6599.1667;

        cout << "\nTest zbieznosci dla calki trygonometrycznej:" << endl;
        cout << setw(15) << "Liczba podzialow" << setw(20) << "Metoda prostokatow"
             << setw(20) << "Metoda trapezow" << setw(20) << "Metoda Simpsona" << endl;

        for (int n : liczbyPodzialow) {
            double wynikProstokatow = MetodaProstokatowTryg(a, b, n);
            double wynikTrapezow = MetodaTrapezowTryg(a, b, n);
            double wynikSimpson = SimpsonTryg(a, b, n);

            cout << setw(15) << n
                 << setw(20) << fixed << setprecision(10) << wynikProstokatow
                 << setw(20) << wynikTrapezow
                 << setw(20) << wynikSimpson << endl;
        }
    }

    bool wczytajDane(const string& filename, int& degree, vector<double>& wspolczynniki, double& a, double& b) {
        ifstream file(filename);

        if (!file.is_open()) {
            cerr << "Blad: Nie mozna otworzyc pliku " << filename << endl;
            return false;
        }

        string line;

        getline(file, line);
        getline(file, line);
        size_t pos = line.find(':');
        if (pos != string::npos) {
            degree = stoi(line.substr(pos + 1));
        } else {
            cerr << "Blad formatu przy odczycie stopnia wielomianu" << endl;
            return false;
        }

        getline(file, line);

        getline(file, line);
        istringstream iss(line);
        double wsp;
        while (iss >> wsp) {
            wspolczynniki.push_back(wsp);
        }
        reverse(wspolczynniki.begin(), wspolczynniki.end());

        getline(file, line);
        pos = line.find(':');
        if (pos != string::npos) {
            string interval = line.substr(pos + 1);
            istringstream interval_stream(interval);
            interval_stream >> a >> b;
        } else {
            cerr << "Blad formatu przy odczycie przedziału całkowania" << endl;
            return false;
        }

        if (wspolczynniki.size() != degree + 1) {
            cout << "Uwaga: Liczba wspolczynnikow (" << wspolczynniki.size()
                    << ") nie zgadza się ze stopniem wielomianu (" << degree << ")." << endl;
        }
        file.close();
        return true;
    }
    double wartoscFunkcji1(double argument1) {
        return pow(argument1, 2) * pow(sin(argument1), 3);
    }

    double wartoscFunkcji2(double argument2) {
        return exp(pow(argument2, 2)) * (1 - argument2);
    }

    double wartoscFunkcji3(double argument3) {
        return -14*pow(argument3,6) + 18*pow(argument3,5) + 5*pow(argument3,4) - 2*pow(argument3,3) - 19*pow(argument3,2) + 14*argument3 - 5;
    }

    double wartoscFunkcji4(double argument4) {
        return argument4*pow(cos(argument4), 3);
    }

    // Podstawowa kwadratura GL bez podziału na podprzedziały
    double kwadraturaGL(unordered_map<int, vector<double>>& punkty, unordered_map<int, vector<double>>& wagi,
                        double poczatekPrzedzialu, double koniecPrzedzialu, int liczbaPunktow,
                        double (*funkcja)(double)) {
        double wynik = 0;
        double stala = (koniecPrzedzialu - poczatekPrzedzialu) / 2;
        double stala2 = (poczatekPrzedzialu + koniecPrzedzialu) / 2;

        for (int i = 0; i < liczbaPunktow; i++) {
            double argument = stala * punkty[liczbaPunktow-1][i] + stala2;
            wynik += wagi[liczbaPunktow-1][i] * funkcja(argument);
        }
        return stala * wynik;
    }


    double kwadraturaGLPodzial(unordered_map<int, vector<double>>& punkty, unordered_map<int, vector<double>>& wagi,
                              double poczatekPrzedzialu, double koniecPrzedzialu, int liczbaPunktow,
                              double (*funkcja)(double), int liczbaPodzialow) {
        double wynik = 0;
        double dlugoscPrzedzialu = (koniecPrzedzialu - poczatekPrzedzialu) / liczbaPodzialow;

        for (int i = 0; i < liczbaPodzialow; i++) {
            double poczatekPodprzedzialu = poczatekPrzedzialu + i * dlugoscPrzedzialu;
            double koniecPodprzedzialu = poczatekPodprzedzialu + dlugoscPrzedzialu;

            wynik += kwadraturaGL(punkty, wagi, poczatekPodprzedzialu, koniecPodprzedzialu, liczbaPunktow, funkcja);
        }

        return wynik;
    }


    void analizaZbieznosci(unordered_map<int, vector<double>>& punkty, unordered_map<int, vector<double>>& wagi,
                          double poczatekPrzedzialu, double koniecPrzedzialu, double dokladnaWartosc,
                          double (*funkcja)(double), const string& nazwaFunkcji) {

        vector<int> podzialy = {4, 8, 16, 32, 64, 128};

        cout << "\n===== Analiza zbieznosci dla funkcji " << nazwaFunkcji << " =====" << endl;
        cout << "Dokladna wartosc calki: " << dokladnaWartosc << endl;
        cout << setw(10) << "Podzialy" << setw(15) << "GL(2 wezly)" << setw(15) << "GL(3 wezly)" << setw(15) << "GL(4 wezly)" << endl;

        for (int n : podzialy) {
            double wynik2 = kwadraturaGLPodzial(punkty, wagi, poczatekPrzedzialu, koniecPrzedzialu, 2, funkcja, n);
            double wynik3 = kwadraturaGLPodzial(punkty, wagi, poczatekPrzedzialu, koniecPrzedzialu, 3, funkcja, n);
            double wynik4 = kwadraturaGLPodzial(punkty, wagi, poczatekPrzedzialu, koniecPrzedzialu, 4, funkcja, n);

            double blad2 = fabs(dokladnaWartosc - wynik2);
            double blad3 = fabs(dokladnaWartosc - wynik3);
            double blad4 = fabs(dokladnaWartosc - wynik4);

            cout << setw(10) << n
                 << setw(15) << scientific << setprecision(6) << blad2
                 << setw(15) << scientific << setprecision(6) << blad3
                 << setw(15) << scientific << setprecision(6) << blad4 << endl;
        }
    }
}