//
// Created by Szymon Ros on 10/06/2025.
//
#ifndef INTERPOLACJA_H
#define INTERPOLACJA_H
#include <vector>
using namespace std;

namespace biblioteka_numeryczna {
    /**
     * @brief Interpolacja wielomianowa metodą Lagrange'a
     * @param wezly_x Wektor punktów węzłowych (argumenty)
     * @param wartosci_y Wektor wartości funkcji w punktach węzłowych
     * @param punkt_x Punkt, w którym obliczamy wartość wielomianu
     * @return Wartość wielomianu interpolacyjnego w punkcie x
     *
     * Przykład użycia:
     * vector<double> wezly = {0, 1, 2};
     * vector<double> wartosci = {1, 4, 9};
     * double wynik = interpolacjaLagrangea(wezly, wartosci, 1.5);
     */
    double interpolacjaLagrangea(const vector<double>& wezly_x,
                                const vector<double>& wartosci_y,
                                double punkt_x);

    /**
     * @brief Oblicza współczynniki wielomianu interpolacyjnego Newtona
     * @param wezly_x Wektor punktów węzłowych
     * @param wartosci_y Wektor wartości funkcji w punktach węzłowych
     * @return Macierz różnic dzielonych (tablica Newtona)
     *
     * Przykład użycia:
     * vector<double> wezly = {0, 1, 2};
     * vector<double> wartosci = {1, 4, 9};
     * auto tablica = wspolczynnikiNewton(wezly, wartosci);
     */
    vector<vector<double>> wspolczynnikiNewton(const vector<double>& wezly_x,
                                              const vector<double>& wartosci_y);

    /**
     * @brief Interpolacja wielomianowa metodą Newtona
     * @param wezly_x Wektor punktów węzłowych
     * @param tablica_newton Tablica współczynników Newtona (z funkcji wspolczynnikiNewton)
     * @param punkt_x Punkt, w którym obliczamy wartość wielomianu
     * @return Wartość wielomianu interpolacyjnego w punkcie x
     *
     * Przykład użycia:
     * vector<double> wezly = {0, 1, 2};
     * vector<double> wartosci = {1, 4, 9};
     * auto tablica = wspolczynnikiNewton(wezly, wartosci);
     * double wynik = interpolacjaNewtona(wezly, tablica, 1.5);
     */
    double interpolacjaNewtona(const vector<double>& wezly_x,
                              const vector<vector<double>>& tablica_newton,
                              double punkt_x);

    /**
     * @brief Metoda naturalna obliczenia wartości wielomianu (dla porównania wydajności)
     * @param wspolczynniki Współczynniki wielomianu (od a0 do an)
     * @param x Argument
     * @return Wartość wielomianu
     */
    double metodaNaturalna(const vector<double>& wspolczynniki, double x);

    /**
     * @brief Metoda Hornera obliczenia wartości wielomianu
     * @param wspolczynniki Współczynniki wielomianu (od a0 do an)
     * @param x Argument
     * @return Wartość wielomianu
     */
    double metodaHornera(const vector<double>& wspolczynniki, double x);

    /**
     * @brief Oblicza średni błąd kwadratowy interpolacji
     * @param wszystkie_x Wszystkie punkty danych
     * @param wszystkie_y Wszystkie wartości funkcji
     * @param wezly_x Punkty węzłowe użyte do interpolacji
     * @param wezly_y Wartości w punktach węzłowych
     * @return Średni błąd kwadratowy
     */
    double sredniaBladKwadratowy(const vector<double>& wszystkie_x,
                                const vector<double>& wszystkie_y,
                                const vector<double>& wezly_x,
                                const vector<double>& wezly_y);

    /**
     * @brief Wybiera węzły co k-ty element z danych
     * @param dane_x Wszystkie punkty x
     * @param dane_y Wszystkie wartości y
     * @param krok Krok wybierania węzłów (co krok-ty element)
     * @param wezly_x Wynikowe węzły x
     * @param wezly_y Wynikowe węzły y
     */
    void wybierzWezly(const vector<double>& dane_x,
                     const vector<double>& dane_y,
                     int krok,
                     vector<double>& wezly_x,
                     vector<double>& wezly_y);
}

#endif // INTERPOLACJA_H