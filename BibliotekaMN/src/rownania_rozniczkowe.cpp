//
// Created by Szymon Ros on 11/06/2025.
//
#include <iostream>
#include <vector>
#include "../include/rownania_rozniczkowe.h"
using namespace std;

namespace biblioteka_numeryczna {
    const double alfa = 3* pow(10, -12);
    const double Beta = 0;
    const double t0 = 1000.0;
    const double czasKoncowy = 100.0;
    const double T0 = 2973;
    const double czas_chlodzenia = 2973.0;

    double obliczPochodna(double T) {
        return -alfa * (pow(T, 4) - Beta);
    }

    double rozwiazanieDokladne(double t) {          //rozwiazanie policzone metoda rozdzielania zmiennych
        return t0 / pow(1 + 3 * alfa * pow(t0, 3) * t, 1.0/3.0);
    }
    double temperatura_dokladna(double t) {
        double mianownik = 1.0 + 3.0 * alfa * pow(T0, 3) * t;
        return T0 / pow(mianownik, 1.0/3.0);
    }
    vector<pair<double, double>> metodaEulera(double krok) {
        vector<pair<double, double>> wyniki;
        double T = t0;
        double t = 0.0;

        wyniki.push_back({t, T});

        while (t < czasKoncowy) {
            double dT = krok * obliczPochodna(T);
            T += dT;
            t += krok;

            wyniki.push_back({t, T});
        }

        return wyniki;
    }
    double oblicz_sredni_blad_kwadratowy(const vector<double>& wyniki, double krok) {
        double suma = 0.0;
        int liczba_punktow = 0;

        for (size_t i = 0; i < wyniki.size(); i++) {
            double t = i * krok;
            double T_dokladna = temperatura_dokladna(t);
            double roznica = wyniki[i] - T_dokladna;
            suma += roznica * roznica;
            liczba_punktow++;
        }

        return sqrt(suma / liczba_punktow);
    }

    double obliczBladSredniKwadratowy(const vector<pair<double, double>>& rozwiazanieNumeryczne) {
        double sumaKwadratowBledu = 0.0;
        int licznik = 0;

        for (const auto& punkt : rozwiazanieNumeryczne) {
            double t = punkt.first;
            double T_numeryczne = punkt.second;
            double T_dokladne = rozwiazanieDokladne(t);

            double blad = T_numeryczne - T_dokladne;
            sumaKwadratowBledu += blad * blad;
            licznik++;
        }

        return sumaKwadratowBledu / licznik;
    }
    double oblicz_zmiane_temperatury(double T) {
        if (T <= 0.0) return 0.0;
        return -alfa * pow(T, 4);
    }
    vector<double> oblicz_metoda_Eulera(double krok) {
        vector<double> wyniki = {T0};
        double T = T0;

        for (double t = 0; t < czas_chlodzenia; t += krok) {
            double zmiana = oblicz_zmiane_temperatury(T);
            T += zmiana * krok;
            if (T < 0) T = 0;
            wyniki.push_back(T);
        }
        return wyniki;
    }
    vector<double> oblicz_metoda_srodka(double krok) {
        vector<double> wyniki = {T0};
        double T = T0;

        for (double t = 0; t < czas_chlodzenia; t += krok) {
            double zmiana1 = oblicz_zmiane_temperatury(T);
            double T_polowa = T + zmiana1 * krok/2;
            double zmiana2 = oblicz_zmiane_temperatury(T_polowa);
            T += zmiana2 * krok;
            if (T < 0) T = 0;
            wyniki.push_back(T);
        }
        return wyniki;
    }
    vector<double> oblicz_metoda_Heuna(double krok) {
        vector<double> wyniki = {T0};
        double T = T0;

        for (double t = 0; t < czas_chlodzenia; t += krok) {
            double k1 = oblicz_zmiane_temperatury(T);
            double T_predyktor = T + k1 * krok;
            double k2 = oblicz_zmiane_temperatury(T_predyktor);
            T += (k1 + k2) * krok / 2;
            if (T < 0) T = 0;
            wyniki.push_back(T);
        }
        return wyniki;
    }
    vector<double> oblicz_metoda_RK4(double krok) {
        vector<double> wyniki = {T0};
        double T = T0;

        for (double t = 0; t < czas_chlodzenia; t += krok) {
            double k1 = oblicz_zmiane_temperatury(T);
            double k2 = oblicz_zmiane_temperatury(T + k1 * krok/2);
            double k3 = oblicz_zmiane_temperatury(T + k2 * krok/2);
            double k4 = oblicz_zmiane_temperatury(T + k3 * krok);
            T += (k1 + 2*k2 + 2*k3 + k4) * krok / 6;
            if (T < 0) T = 0;
            wyniki.push_back(T);
        }
        return wyniki;
    }
}