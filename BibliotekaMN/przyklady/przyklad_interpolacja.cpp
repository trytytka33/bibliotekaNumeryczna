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
    // interpolacja funkcji kwadratowej
    vector<double> wezly_x = {-2, -1, 0, 1, 2};
    vector<double> wartosci_y = {4, 1, 0, 1, 4}; // f(x) = x^2

    cout << "=== Przykład interpolacji funkcji f(x) = x^2 ===" << endl;
    cout << "Punkty węzłowe:" << endl;
    for (size_t i = 0; i < wezly_x.size(); i++) {
        cout << "(" << wezly_x[i] << ", " << wartosci_y[i] << ")" << endl;
    }

    // interpolacja w roznych punktach
    vector<double> punkty_testowe = {-1.5, -0.5, 0.5, 1.5};
    cout << "\nInterpolacja w punktach testowych:" << endl;
    for (double x : punkty_testowe) {
        double wynik = biblioteka_numeryczna::interpolacjaLagrangea(wezly_x, wartosci_y, x);
        double rzeczywista = x * x; 
        cout << "x = " << x << ", interpolacja = " << wynik
             << ", rzeczywista = " << rzeczywista
             << ", błąd = " << abs(wynik - rzeczywista) << endl;
    }

    //  wczytywanie danych z pliku 
    cout << "\nDane z pliku" << endl;

    vector<double> wszystkie_x, wszystkie_y;
    for (int i = 0; i <= 20; i++) {
        wszystkie_x.push_back(i * 0.5);
        wszystkie_y.push_back(sin(i * 0.5)); // y = sin(x)
    }

    cout << "Analiza bledu interpolacji dla funkcji sin(x):" << endl;

    //różne kroki węzłów 
    for (int krok = 1; krok <= 5; krok++) {
        vector<double> wezly_wybranych_x, wezly_wybranych_y;
        biblioteka_numeryczna::wybierzWezly(wszystkie_x, wszystkie_y, krok,
                                           wezly_wybranych_x, wezly_wybranych_y);

        double sredni_blad = biblioteka_numeryczna::sredniaBladKwadratowy(
            wszystkie_x, wszystkie_y, wezly_wybranych_x, wezly_wybranych_y);

        cout << "Wezly co " << krok << " punkt (lacznie " << wezly_wybranych_x.size()
             << " wezlow), sredni blad kwadratowy = " << sredni_blad << endl;
    }

    // interpolacja w wybranym punkcie
    cout << "\nPodaj wartość argumentu dla interpolacji: ";
    double argument;
    cin >> argument;

    vector<double> wezly_co_5_x, wezly_co_5_y;
    biblioteka_numeryczna::wybierzWezly(wszystkie_x, wszystkie_y, 5, wezly_co_5_x, wezly_co_5_y);

    double wynik_interpolacji = biblioteka_numeryczna::interpolacjaLagrangea(
        wezly_co_5_x, wezly_co_5_y, argument);

    cout << "Wartosc wielomianu interpolacyjnego dla argumentu " << argument
         << " wynosi: " << wynik_interpolacji << endl;
    cout << "Rzeczywista wartosc sin(" << argument << ") = " << sin(argument) << endl;
    cout << "Blad interpolacji = " << abs(wynik_interpolacji - sin(argument)) << endl;

    return 0;
}
