//
// Created by Szymon Ros on 10/06/2025.
//
// BIBLIOTEKA METOD NUMERYCZNYCH
// 
//
// Główny plik nagłówkowy biblioteki metod numerycznych.
// 
//
// Dostępne moduły:
// - uklady_liniowe.h: Rozwiązywanie układów równań liniowych
// - interpolacja.h: Interpolacja wielomianowa (Lagrange, Newton)
// - aproksymacja.h: Aproksymacja funkcji metodą najmniejszych kwadratów
// - calkowanie_numeryczne.h: Całkowanie numeryczne (prostokąty, trapezy, Simpson, Gauss)
// - rownania_rozniczkowe.h: Rozwiązywanie równań różniczkowych zwyczajnych
// - rownania_nieliniowe.h: Znajdowanie pierwiastków równań nieliniowych
//
// Przykład użycia:
// #include "biblioteka_metody_numeryczne.h"
// using namespace biblioteka_numeryczna;
//
// int main() {
//     // Rozwiązywanie układu równań
//     vector<vector<double>> A = {{2, 1}, {1, 3}};
//     vector<double> b = {3, 4};
//     vector<double> x = eliminacjaGaussa(A, b);
//
//     // Interpolacja
//     vector<double> wezly = {0, 1, 2};
//     vector<double> wartosci = {1, 4, 9};
//     double wynik = interpolacjaLagrangea(wezly, wartosci, 1.5);
//
//     return 0;
// }

#ifndef BIBLIOTEKA_NUMERYCZNA_H
#define BIBLIOTEKA_NUMERYCZNA_H

#include "uklady_liniowe.h"
#include "interpolacja.h"
#include "aproksymacja.h"
#include "calkowanie_numeryczne.h"
#include "rownania_rozniczkowe.h"
#include "rownania_nieliniowe.h"

#endif // BIBLIOTEKA_NUMERYCZNA_H
