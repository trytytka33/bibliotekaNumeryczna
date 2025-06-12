#include "../include/interpolacja.h"
#include <cmath>
#include <stdexcept>
using namespace std;

namespace biblioteka_numeryczna {
    double interpolacjaLagrangea(const vector<double>& wezly_x,
                                const vector<double>& wartosci_y,
                                double punkt_x) {
        if (wezly_x.size() != wartosci_y.size() || wezly_x.empty()) {
            throw invalid_argument("Wektory wezly_x i wartosci_y muszą mieć taki sam rozmiar i nie mogą być puste");
        }

        int n = wezly_x.size();
        double suma = 0;

        for (int i = 0; i < n; i++) {
            double skladnik = wartosci_y[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    if (abs(wezly_x[i] - wezly_x[j]) < 1e-10) {
                        throw invalid_argument("Punkty węzłowe nie mogą być identyczne");
                    }
                    skladnik *= (punkt_x - wezly_x[j]) / (wezly_x[i] - wezly_x[j]);
                }
            }
            suma += skladnik;
        }
        return suma;
    }

    vector<vector<double>> wspolczynnikiNewton(const vector<double>& wezly_x,
                                              const vector<double>& wartosci_y) {
        if (wezly_x.size() != wartosci_y.size() || wezly_x.empty()) {
            throw invalid_argument("Wektory wezly_x i wartosci_y muszą mieć taki sam rozmiar i nie mogą być puste");
        }

        int n = wezly_x.size();
        vector<vector<double>> tablica(n, vector<double>(n, 0));

        // Pierwsza kolumna to wartości funkcji
        for (int i = 0; i < n; i++) {
            tablica[i][0] = wartosci_y[i];
        }

        // Oblicz różnice dzielone
        for (int j = 1; j < n; j++) {
            for (int i = j; i < n; i++) {
                if (abs(wezly_x[i] - wezly_x[i-j]) < 1e-10) {
                    throw invalid_argument("Punkty węzłowe nie mogą być identyczne");
                }
                tablica[i][j] = (tablica[i][j-1] - tablica[i-1][j-1]) / (wezly_x[i] - wezly_x[i-j]);
            }
        }
        return tablica;
    }

    double interpolacjaNewtona(const vector<double>& wezly_x,
                              const vector<vector<double>>& tablica_newton,
                              double punkt_x) {
        if (wezly_x.empty() || tablica_newton.empty()) {
            throw invalid_argument("Wektory nie mogą być puste");
        }

        int n = wezly_x.size();
        if (tablica_newton.size() != n) {
            throw invalid_argument("Rozmiary wektorów są niezgodne");
        }

        double wynik = tablica_newton[0][0];

        for (int i = 1; i < n; i++) {
            double temp = tablica_newton[i][i];
            for (int j = 0; j < i; j++) {
                temp *= (punkt_x - wezly_x[j]);
            }
            wynik += temp;
        }

        return wynik;
    }

    double metodaNaturalna(const vector<double>& wspolczynniki, double x) {
        if (wspolczynniki.empty()) {
            throw invalid_argument("Wektor współczynników nie może być pusty");
        }

        double wynik = 0.0;
        double potega = 1.0;
        int rozmiar = wspolczynniki.size();

        for (int i = 0; i < rozmiar; i++) {
            wynik += wspolczynniki[i] * potega;
            potega *= x;
        }
        return wynik;
    }

    double metodaHornera(const vector<double>& wspolczynniki, double x) {
        if (wspolczynniki.empty()) {
            throw invalid_argument("Wektor współczynników nie może być pusty");
        }

        int rozmiar = wspolczynniki.size();
        double wynik = 0.0;

        for (int i = 0; i < rozmiar; i++) {
            wynik = wynik * x + wspolczynniki[rozmiar - 1 - i];
        }
        return wynik;
    }

    double sredniaBladKwadratowy(const vector<double>& wszystkie_x,
                                const vector<double>& wszystkie_y,
                                const vector<double>& wezly_x,
                                const vector<double>& wezly_y) {
        if (wszystkie_x.size() != wszystkie_y.size()) {
            throw invalid_argument("Wektory wszystkie_x i wszystkie_y muszą mieć taki sam rozmiar");
        }

        double suma_bledow = 0;
        int liczba_punktow = 0;

        for (size_t i = 0; i < wszystkie_x.size(); i++) {
            // Sprawdź czy punkt nie jest węzłem
            bool jest_wezlem = false;
            for (size_t j = 0; j < wezly_x.size(); j++) {
                if (abs(wszystkie_x[i] - wezly_x[j]) < 1e-10) {
                    jest_wezlem = true;
                    break;
                }
            }

            if (!jest_wezlem) {
                double interpolowana_wartosc = interpolacjaLagrangea(wezly_x, wezly_y, wszystkie_x[i]);
                suma_bledow += pow(wszystkie_y[i] - interpolowana_wartosc, 2);
                liczba_punktow++;
            }
        }

        return liczba_punktow > 0 ? suma_bledow / liczba_punktow : 0;
    }

    void wybierzWezly(const vector<double>& dane_x,
                     const vector<double>& dane_y,
                     int krok,
                     vector<double>& wezly_x,
                     vector<double>& wezly_y) {
        wezly_x.clear();
        wezly_y.clear();

        for (size_t i = 0; i < dane_x.size() && i < dane_y.size(); i += krok) {
            wezly_x.push_back(dane_x[i]);
            wezly_y.push_back(dane_y[i]);
        }
    }
}