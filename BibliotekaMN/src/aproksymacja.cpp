//
// Created by Szymon Ros on 11/06/2025.
//
#include "../include/aproksymacja.h"
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;
namespace biblioteka_numeryczna {
    double funkcjaAproksymowana(double x) {
        return exp(x) * cos(6 * x) - pow(x, 3) + 5 * pow(x, 2) - 10;
    }

    //x^i * x^j = x^(i+j)
    double iloczynSkalarnyElementowZBazy(int pierwszyElement, int drugiElement, const double poczatekPrzedzialu, const double koniecPrzedzialu, int n = 100){
        int sumaPoteg = pierwszyElement + drugiElement;                                                                                                          //sluzy jako wykladnik w wyniku
        if (n % 2 != 0) {
            n += 1;
        }
        double dx = (koniecPrzedzialu - poczatekPrzedzialu) / n;
        double wynik = pow(poczatekPrzedzialu, sumaPoteg) + pow(koniecPrzedzialu, sumaPoteg);                                                                                    //suma pierwszego i ostatniego elementu, bo zaczynamy obliczanie wartosci od 2 elementów, a nie 1 elementu

        for (int i = 1; i < n; i++) {
            double x = poczatekPrzedzialu + i * dx;
            if (i % 2 == 0) {
                wynik += 2 * pow(x, sumaPoteg);
            } else {
                wynik += 4 * pow(x, sumaPoteg);
            }
        }
        wynik *= dx / 3.0;
        return wynik;
    }
    //x^i * f(x)
    double iloczynSkalarnyElementuZBazyIFunkcji(int elementZBazy, const double poczatekPrzedzialu, const double koniecPrzedzialu, int n = 100) {
        if (n % 2 != 0) {
            n+=1;
        }
        double dx = (koniecPrzedzialu - poczatekPrzedzialu) / n;
        double wynik = funkcjaAproksymowana(poczatekPrzedzialu) * pow(poczatekPrzedzialu, elementZBazy) + funkcjaAproksymowana(koniecPrzedzialu) * pow(koniecPrzedzialu, elementZBazy);
        for (int i = 1; i < n; i++) {
            double x = poczatekPrzedzialu + i * dx;
            if (i % 2 == 0) {
                wynik += 2 * funkcjaAproksymowana(x) * pow(x, elementZBazy);
            } else {
                wynik += 4 * funkcjaAproksymowana(x) * pow(x, elementZBazy);
            }
        }
        wynik *= dx / 3.0;
        return wynik;

    }
    void rozkladLU_zPivotingiem(const vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, vector<int>& P);

    vector<double> permutujWektor(const vector<double>& b, const vector<int>& P);

    vector<double> rozwiazLy_b(const vector<vector<double>>& L, const vector<double>& b);

    vector<double> rozwiazUx_y(const vector<vector<double>>& U, const vector<double>& y);
    double sprawdzPoprawnosc(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);
    pair<vector<double>, vector<double>> rozwiazUkladLU(const vector<vector<double>>& A, const vector<double>& b);
    double bladAproksymacji(const vector<double>& wspolczynniki, double poczatekPrzedzialu, double koniecPrzedzialu, int liczbaPunktow = 100) {
        double maxBlad = 0.0;
        double delta = (koniecPrzedzialu - poczatekPrzedzialu) / liczbaPunktow;

        for (int i = 0; i <= liczbaPunktow; i++) {
            double x = poczatekPrzedzialu + i * delta;
            double wartoscFunkcji = funkcjaAproksymowana(x);
            double wartoscWielomianu = 0.0;

            for (int j = 0; j < wspolczynniki.size(); j++) {
                wartoscWielomianu += wspolczynniki[j] * pow(x, j);
            }
            double blad = fabs(wartoscFunkcji - wartoscWielomianu);
            if (blad > maxBlad) {
                maxBlad = blad;
            }
        }
        return maxBlad;
    }


    void pokazAproksymacje(const vector<double>& wspolczynniki, double poczatekPrzedzialu, double koniecPrzedzialu, int liczbaPunktow = 10) {
        double delta = (koniecPrzedzialu - poczatekPrzedzialu) / liczbaPunktow;

        cout << fixed << setprecision(6);
        cout << "\nPorownanie funkcji i jej aproksymacji:\n";
        cout << setw(10) << "x" << setw(20) << "f(x)" << setw(20) << "W(x)" << setw(20) << "Blad" << endl;

        for (int i = 0; i <= liczbaPunktow; i++) {
            double x = poczatekPrzedzialu + i * delta;
            double wartoscFunkcji = funkcjaAproksymowana(x);
            double wartoscWielomianu = 0.0;

            for (int j = 0; j < wspolczynniki.size(); j++) {
                wartoscWielomianu += wspolczynniki[j] * pow(x, j);
            }

            double blad = fabs(wartoscFunkcji - wartoscWielomianu);
            cout << setw(10) << x << setw(20) << wartoscFunkcji << setw(20) << wartoscWielomianu << setw(20) << blad << endl;
        }
    }
}