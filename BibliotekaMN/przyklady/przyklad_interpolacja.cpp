//
// Created by Szymon Ros on 10/06/2025.
//
#include "../include/biblioteka_metody_numeryczne.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
using namespace std;

int main() {
    // Przykład 1: Interpolacja funkcji kwadratowej
    vector<double> wezly_x = {-2, -1, 0, 1, 2};
    vector<double> wartosci_y = {4, 1, 0, 1, 4}; // f(x) = x^2

    cout << "=== Przykład interpolacji funkcji f(x) = x^2 ===" << endl;
    cout << "Punkty węzłowe:" << endl;
    for (size_t i = 0; i < wezly_x.size(); i++) {
        cout << "(" << wezly_x[i] << ", " << wartosci_y[i] << ")" << endl;
    }

    // Interpolacja w różnych punktach
    vector<double> punkty_testowe = {-1.5, -0.5, 0.5, 1.5};
    cout << "\nInterpolacja w punktach testowych:" << endl;
    for (double x : punkty_testowe) {
        double wynik = biblioteka_numeryczna::interpolacjaLagrangea(wezly_x, wartosci_y, x);
        double rzeczywista = x * x; // rzeczywista wartość
        cout << "x = " << x << ", interpolacja = " << wynik
             << ", rzeczywista = " << rzeczywista
             << ", błąd = " << abs(wynik - rzeczywista) << endl;
    }

    // Przykład 2: Wczytywanie danych z pliku (podobnie jak w oryginalnym kodzie)
    cout << "\n=== Przykład z danymi z pliku ===" << endl;

    // Symulujemy dane z pliku
    vector<double> wszystkie_x, wszystkie_y;
    for (int i = 0; i <= 20; i++) {
        wszystkie_x.push_back(i * 0.5);
        wszystkie_y.push_back(sin(i * 0.5)); // f(x) = sin(x)
    }

    cout << "Analiza błędu interpolacji dla funkcji sin(x):" << endl;

    // Testuj różne kroki węzłów (podobnie jak w oryginalnym kodzie)
    for (int krok = 1; krok <= 5; krok++) {
        vector<double> wezly_wybranych_x, wezly_wybranych_y;
        biblioteka_numeryczna::wybierzWezly(wszystkie_x, wszystkie_y, krok,
                                           wezly_wybranych_x, wezly_wybranych_y);

        double sredni_blad = biblioteka_numeryczna::sredniaBladKwadratowy(
            wszystkie_x, wszystkie_y, wezly_wybranych_x, wezly_wybranych_y);

        cout << "Węzły co " << krok << " punkt (łącznie " << wezly_wybranych_x.size()
             << " węzłów), średni błąd kwadratowy = " << sredni_blad << endl;
    }

    // Interpolacja w wybranym punkcie
    cout << "\nPodaj wartość argumentu dla interpolacji: ";
    double argument;
    cin >> argument;

    vector<double> wezly_co_5_x, wezly_co_5_y;
    biblioteka_numeryczna::wybierzWezly(wszystkie_x, wszystkie_y, 5, wezly_co_5_x, wezly_co_5_y);

    double wynik_interpolacji = biblioteka_numeryczna::interpolacjaLagrangea(
        wezly_co_5_x, wezly_co_5_y, argument);

    cout << "Wartość wielomianu interpolacyjnego dla argumentu " << argument
         << " wynosi: " << wynik_interpolacji << endl;
    cout << "Rzeczywista wartość sin(" << argument << ") = " << sin(argument) << endl;
    cout << "Błąd interpolacji = " << abs(wynik_interpolacji - sin(argument)) << endl;

    return 0;
}