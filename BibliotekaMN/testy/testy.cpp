//
// Created by Szymon Ros on 11/06/2025.
//
//
// Testy jednostkowe dla biblioteki metod numerycznych
// Autor: Szymon Ros
// Data: 11/06/2025
//

#include "../include/biblioteka_metody_numeryczne.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace biblioteka_numeryczna;

// Pomocnicze funkcje testowe
const double EPSILON = 1e-6;

bool porownajDouble(double a, double b, double eps = EPSILON) {
    return abs(a - b) < eps;
}

void wypiszWynikTestu(const string& nazwaTestu, bool wynik) {
    cout << "[" << (wynik ? "PASS" : "FAIL") << "] " << nazwaTestu << endl;
}

// =============================================================================
// TESTY DLA MODUŁU UKŁADY LINIOWE
// =============================================================================

void testEliminacjaGaussa1() {
    cout << "\n=== TESTY UKŁADÓW LINIOWYCH ===" << endl;

    // Test 1: Prosty układ 2x2
    vector<vector<double>> A = {{2, 1}, {1, 3}};
    vector<double> b = {3, 4};

    vector<double> x = eliminacjaGaussa(A, b);

    // Sprawdzenie: 2*x[0] + 1*x[1] = 3, 1*x[0] + 3*x[1] = 4
    // Rozwiązanie: x[0] = 1, x[1] = 1
    bool test1 = porownajDouble(x[0], 1.0) && porownajDouble(x[1], 1.0);
    wypiszWynikTestu("eliminacjaGaussa - test 1 (układ 2x2)", test1);
}

void testEliminacjaGaussa2() {
    // Test 2: Układ 3x3
    vector<vector<double>> A = {{1, 2, 3}, {2, -1, 1}, {3, 0, -1}};
    vector<double> b = {9, 8, 3};

    vector<double> x = eliminacjaGaussa(A, b);

    // Weryfikacja przez podstawienie
    double sprawdzenie1 = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
    double sprawdzenie2 = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
    double sprawdzenie3 = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];

    bool test2 = porownajDouble(sprawdzenie1, b[0]) &&
                 porownajDouble(sprawdzenie2, b[1]) &&
                 porownajDouble(sprawdzenie3, b[2]);
    wypiszWynikTestu("eliminacjaGaussa - test 2 (układ 3x3)", test2);
}

void testMetodaGaussaSeidela1() {
    // Test 1: Układ z dominującą przekątną
    vector<vector<double>> A = {{4, 1}, {1, 3}};
    vector<double> b = {5, 4};
    vector<double> x0 = {0, 0};

    vector<double> x = metodaGaussaSeidela(A, b, x0);

    // Weryfikacja przez podstawienie
    double sprawdzenie1 = A[0][0]*x[0] + A[0][1]*x[1];
    double sprawdzenie2 = A[1][0]*x[0] + A[1][1]*x[1];

    bool test1 = porownajDouble(sprawdzenie1, b[0]) && porownajDouble(sprawdzenie2, b[1]);
    wypiszWynikTestu("metodaGaussaSeidela - test 1", test1);
}

void testMetodaGaussaSeidela2() {
    // Test 2: Większy układ
    vector<vector<double>> A = {{5, 1, 1}, {1, 4, 1}, {1, 1, 3}};
    vector<double> b = {7, 6, 5};
    vector<double> x0 = {0, 0, 0};

    vector<double> x = metodaGaussaSeidela(A, b, x0);

    // Weryfikacja
    double sprawdzenie1 = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
    double sprawdzenie2 = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
    double sprawdzenie3 = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];

    bool test2 = porownajDouble(sprawdzenie1, b[0], 1e-3) &&
                 porownajDouble(sprawdzenie2, b[1], 1e-3) &&
                 porownajDouble(sprawdzenie3, b[2], 1e-3);
    wypiszWynikTestu("metodaGaussaSeidela - test 2", test2);
}

void testRozkladLU1() {
    // Test 1: Prosta macierz 2x2
    vector<vector<double>> A = {{2, 1}, {1, 1}};
    vector<vector<double>> L, U;
    vector<int> P;

    rozkladLU_zPivotingiem(A, L, U, P);

    // Sprawdzenie wymiarów
    bool test1 = (L.size() == 2 && U.size() == 2 && P.size() == 2);
    wypiszWynikTestu("rozkladLU_zPivotingiem - test 1 (wymiary)", test1);
}

void testRozkladLU2() {
    // Test 2: Macierz 3x3
    vector<vector<double>> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 10}};
    vector<vector<double>> L, U;
    vector<int> P;

    rozkladLU_zPivotingiem(A, L, U, P);

    // Sprawdzenie czy L jest dolna trójkątna i U górna trójkątna
    bool test2 = (L[0][1] == 0 && L[0][2] == 0 && L[1][2] == 0) &&
                 (U[1][0] == 0 && U[2][0] == 0 && U[2][1] == 0);
    wypiszWynikTestu("rozkladLU_zPivotingiem - test 2 (struktura)", test2);
}

// =============================================================================
// TESTY DLA MODUŁU INTERPOLACJA
// =============================================================================

void testInterpolacjaLagrangea1() {
    cout << "\n=== TESTY INTERPOLACJI ===" << endl;

    // Test 1: Interpolacja liniowa
    vector<double> wezly = {0, 1};
    vector<double> wartosci = {1, 3};

    double wynik1 = interpolacjaLagrangea(wezly, wartosci, 0.5);
    bool test1 = porownajDouble(wynik1, 2.0);  // f(0.5) = 2
    wypiszWynikTestu("interpolacjaLagrangea - test 1 (liniowa)", test1);
}

void testInterpolacjaLagrangea2() {
    // Test 2: Interpolacja kwadratowa
    vector<double> wezly = {0, 1, 2};
    vector<double> wartosci = {1, 4, 9};  // f(x) = x^2 + 2x + 1

    double wynik2 = interpolacjaLagrangea(wezly, wartosci, 1.5);
    // f(1.5) = 1.5^2 + 2*1.5 + 1 = 2.25 + 3 + 1 = 6.25
    bool test2 = porownajDouble(wynik2, 6.25);
    wypiszWynikTestu("interpolacjaLagrangea - test 2 (kwadratowa)", test2);
}

void testWspolczynnikiNewton1() {
    // Test 1: Proste dane
    vector<double> wezly = {0, 1, 2};
    vector<double> wartosci = {1, 2, 5};

    auto tablica = wspolczynnikiNewton(wezly, wartosci);

    bool test1 = (tablica.size() == 3 && tablica[0].size() >= 1);
    wypiszWynikTestu("wspolczynnikiNewton - test 1 (wymiary)", test1);
}

void testWspolczynnikiNewton2() {
    // Test 2: Sprawdzenie pierwszego wiersza (powinien być równy wartościom)
    vector<double> wezly = {1, 2, 3};
    vector<double> wartosci = {2, 8, 18};

    auto tablica = wspolczynnikiNewton(wezly, wartosci);

    bool test2 = porownajDouble(tablica[0][0], 2.0) &&
                 porownajDouble(tablica[1][0], 8.0) &&
                 porownajDouble(tablica[2][0], 18.0);
    wypiszWynikTestu("wspolczynnikiNewton - test 2 (pierwszy wiersz)", test2);
}

void testInterpolacjaNewtona1() {
    // Test 1: Porównanie z Lagrangem
    vector<double> wezly = {0, 1, 2};
    vector<double> wartosci = {1, 4, 9};

    auto tablica = wspolczynnikiNewton(wezly, wartosci);
    double newtonWynik = interpolacjaNewtona(wezly, tablica, 1.5);
    double lagrangeWynik = interpolacjaLagrangea(wezly, wartosci, 1.5);

    bool test1 = porownajDouble(newtonWynik, lagrangeWynik);
    wypiszWynikTestu("interpolacjaNewtona - test 1 (zgodność z Lagrangem)", test1);
}

void testInterpolacjaNewtona2() {
    // Test 2: Wartość w węźle
    vector<double> wezly = {1, 2, 3, 4};
    vector<double> wartosci = {1, 8, 27, 64};

    auto tablica = wspolczynnikiNewton(wezly, wartosci);
    double wynik = interpolacjaNewtona(wezly, tablica, 2.0);

    bool test2 = porownajDouble(wynik, 8.0);  // Powinno być dokładnie 8
    wypiszWynikTestu("interpolacjaNewtona - test 2 (wartość w węźle)", test2);
}

void testMetodaHornera1() {
    // Test 1: Wielomian 1 + 2x + 3x^2
    vector<double> wspolczynniki = {1, 2, 3};
    double wynik = metodaHornera(wspolczynniki, 2.0);
    // 1 + 2*2 + 3*4 = 1 + 4 + 12 = 17

    bool test1 = porownajDouble(wynik, 17.0);
    wypiszWynikTestu("metodaHornera - test 1", test1);
}

void testMetodaHornera2() {
    // Test 2: Porównanie z metodą naturalną
    vector<double> wspolczynniki = {5, -3, 2, 1};
    double x = 1.5;

    double horner = metodaHornera(wspolczynniki, x);
    double naturalna = metodaNaturalna(wspolczynniki, x);

    bool test2 = porownajDouble(horner, naturalna);
    wypiszWynikTestu("metodaHornera - test 2 (zgodność z naturalną)", test2);
}

// =============================================================================
// TESTY DLA MODUŁU APROKSYMACJA
// =============================================================================

void testIloczynSkalarnyElementowZBazy1() {
    cout << "\n=== TESTY APROKSYMACJI ===" << endl;

    // Test 1: Iloczyn <1, 1> na przedziale [0,1]
    double wynik = iloczynSkalarnyElementowZBazy(0, 0, 0.0, 1.0, 1000);
    bool test1 = porownajDouble(wynik, 1.0, 1e-3);  // ∫₀¹ 1·1 dx = 1
    wypiszWynikTestu("iloczynSkalarnyElementowZBazy - test 1", test1);
}

void testIloczynSkalarnyElementowZBazy2() {
    // Test 2: Iloczyn <1, x> na przedziale [0,1]
    double wynik = iloczynSkalarnyElementowZBazy(0, 1, 0.0, 1.0, 1000);
    bool test2 = porownajDouble(wynik, 0.5, 1e-3);  // ∫₀¹ 1·x dx = 1/2
    wypiszWynikTestu("iloczynSkalarnyElementowZBazy - test 2", test2);
}

void testFunkcjaAproksymowana1() {
    // Test 1: Sprawdzenie czy funkcja zwraca wartość
    double wynik = funkcjaAproksymowana(1.0);
    bool test1 = !isnan(wynik) && !isinf(wynik);
    wypiszWynikTestu("funkcjaAproksymowana - test 1 (poprawność)", test1);
}

void testFunkcjaAproksymowana2() {
    // Test 2: Sprawdzenie różnych wartości
    double wynik1 = funkcjaAproksymowana(0.0);
    double wynik2 = funkcjaAproksymowana(1.0);
    bool test2 = wynik1 != wynik2;  // Funkcja powinna być nietrywialną
    wypiszWynikTestu("funkcjaAproksymowana - test 2 (nietrywialność)", test2);
}

void testIloczynSkalarnyElementuZBazyIFunkcji1() {
    // Test 1: Podstawowy test
    double wynik = iloczynSkalarnyElementuZBazyIFunkcji(0, 0.0, 1.0, 100);
    bool test1 = !isnan(wynik) && !isinf(wynik);
    wypiszWynikTestu("iloczynSkalarnyElementuZBazyIFunkcji - test 1", test1);
}

void testIloczynSkalarnyElementuZBazyIFunkcji2() {
    // Test 2: Różne elementy bazy
    double wynik1 = iloczynSkalarnyElementuZBazyIFunkcji(0, 0.0, 1.0, 100);
    double wynik2 = iloczynSkalarnyElementuZBazyIFunkcji(1, 0.0, 1.0, 100);
    bool test2 = abs(wynik1 - wynik2) > 1e-6;  // Powinny być różne
    wypiszWynikTestu("iloczynSkalarnyElementuZBazyIFunkcji - test 2", test2);
}

void testBladAproksymacji1() {
    // Test 1: Podstawowy test błędu
    vector<double> wspolczynniki = {1.0, 0.0, 0.0};
    double blad = bladAproksymacji(wspolczynniki, 0.0, 1.0, 100);
    bool test1 = blad >= 0;  // Błąd powinien być nieujemny
    wypiszWynikTestu("bladAproksymacji - test 1 (nieujemność)", test1);
}

void testBladAproksymacji2() {
    // Test 2: Różne współczynniki
    vector<double> wspolczynniki1 = {1.0};
    vector<double> wspolczynniki2 = {1.0, 1.0};
    double blad1 = bladAproksymacji(wspolczynniki1, 0.0, 1.0, 100);
    double blad2 = bladAproksymacji(wspolczynniki2, 0.0, 1.0, 100);
    bool test2 = blad1 != blad2;
    wypiszWynikTestu("bladAproksymacji - test 2 (różne współczynniki)", test2);
}

// =============================================================================
// TESTY DLA MODUŁU CAŁKOWANIE NUMERYCZNE
// =============================================================================

void testFunkcjaTrygonometryczna1() {
    cout << "\n=== TESTY CAŁKOWANIA NUMERYCZNEGO ===" << endl;

    // Test 1: Wartość w zerze
    double wynik = FunkcjaTrygonometryczna(0.0);
    // sin(0) + cos(0) = 0 + 1 = 1
    bool test1 = porownajDouble(wynik, 1.0);
    wypiszWynikTestu("FunkcjaTrygonometryczna - test 1", test1);
}

void testFunkcjaTrygonometryczna2() {
    // Test 2: Wartość w π/2
    double wynik = FunkcjaTrygonometryczna(M_PI/2);
    // sin(π/2) + cos(π/2) = 1 + 0 = 1
    bool test2 = porownajDouble(wynik, 1.0);
    wypiszWynikTestu("FunkcjaTrygonometryczna - test 2", test2);
}

void testHorner1() {
    // Test 1: Wielomian 1 + x
    vector<double> wspolczynniki = {1, 1};
    double wynik = Horner(wspolczynniki, 2.0);
    bool test1 = porownajDouble(wynik, 3.0);  // 1 + 2 = 3
    wypiszWynikTestu("Horner - test 1", test1);
}

void testHorner2() {
    // Test 2: Wielomian x^2
    vector<double> wspolczynniki = {0, 0, 1};
    double wynik = Horner(wspolczynniki, 3.0);
    bool test2 = porownajDouble(wynik, 9.0);  // 3^2 = 9
    wypiszWynikTestu("Horner - test 2", test2);
}

void testMetodaProstokatow1() {
    // Test 1: Całka z wielomianu stałego
    vector<double> wspolczynniki = {2};  // f(x) = 2
    double wynik = MetodaProstokatow(wspolczynniki, 0.0, 1.0, 1000);
    bool test1 = porownajDouble(wynik, 2.0, 1e-3);  // ∫₀¹ 2 dx = 2
    wypiszWynikTestu("MetodaProstokatow - test 1", test1);
}

void testMetodaProstokatow2() {
    // Test 2: Całka z x
    vector<double> wspolczynniki = {0, 1};  // f(x) = x
    double wynik = MetodaProstokatow(wspolczynniki, 0.0, 2.0, 1000);
    bool test2 = porownajDouble(wynik, 2.0, 1e-2);  // ∫₀² x dx = 2
    wypiszWynikTestu("MetodaProstokatow - test 2", test2);
}

void testMetodaTrapezow1() {
    // Test 1: Funkcja liniowa (powinna być dokładna)
    vector<double> wspolczynniki = {1, 2};  // f(x) = 1 + 2x
    double wynik = MetodaTrapezow(wspolczynniki, 0.0, 1.0, 10);
    double dokladna = 1.0 * 1.0 + 2.0 * 0.5;  // x + x² = 1 + 1 = 2
    bool test1 = porownajDouble(wynik, dokladna, 1e-10);
    wypiszWynikTestu("MetodaTrapezow - test 1", test1);
}

void testMetodaTrapezow2() {
    // Test 2: Funkcja kwadratowa
    vector<double> wspolczynniki = {0, 0, 1};  // f(x) = x²
    double wynik = MetodaTrapezow(wspolczynniki, 0.0, 1.0, 1000);
    bool test2 = porownajDouble(wynik, 1.0/3.0, 1e-3);  // ∫₀¹ x² dx = 1/3
    wypiszWynikTestu("MetodaTrapezow - test 2", test2);
}

void testSimpson1() {
    // Test 1: Funkcja kwadratowa (powinna być dokładna)
    vector<double> wspolczynniki = {1, 0, 1};  // f(x) = 1 + x²
    double wynik = Simpson(wspolczynniki, 0.0, 1.0, 100);
    double dokladna = 1.0 + 1.0/3.0;  // ∫₀¹ (1 + x²) dx = 1 + 1/3 = 4/3
    bool test1 = porownajDouble(wynik, dokladna, 1e-10);
    wypiszWynikTestu("Simpson - test 1", test1);
}

void testSimpson2() {
    // Test 2: Funkcja sześcienna
    vector<double> wspolczynniki = {0, 0, 0, 1};  // f(x) = x³
    double wynik = Simpson(wspolczynniki, 0.0, 2.0, 1000);
    bool test2 = porownajDouble(wynik, 4.0, 1e-6);  // ∫₀² x³ dx = 4
    wypiszWynikTestu("Simpson - test 2", test2);
}

void testCalkaAnalitycznaWielomianu1() {
    // Test 1: Wielomian stały
    vector<double> wspolczynniki = {5};
    double wynik = CalkaAnalitycznaWielomianu(wspolczynniki, 0.0, 2.0);
    bool test1 = porownajDouble(wynik, 10.0);  // ∫₀² 5 dx = 10
    wypiszWynikTestu("CalkaAnalitycznaWielomianu - test 1", test1);
}

void testCalkaAnalitycznaWielomianu2() {
    // Test 2: Wielomian x²
    vector<double> wspolczynniki = {0, 0, 1};
    double wynik = CalkaAnalitycznaWielomianu(wspolczynniki, 1.0, 3.0);
    // ∫₁³ x² dx = [x³/3]₁³ = 27/3 - 1/3 = 26/3
    bool test2 = porownajDouble(wynik, 26.0/3.0);
    wypiszWynikTestu("CalkaAnalitycznaWielomianu - test 2", test2);
}

void testWartoscFunkcji1() {
    // Test 1: Pierwsza funkcja testowa
    double wynik = wartoscFunkcji1(1.0);
    bool test1 = !isnan(wynik) && !isinf(wynik);
    wypiszWynikTestu("wartoscFunkcji1 - test 1", test1);
}

void testWartoscFunkcji1_2() {
    // Test 2: Różne wartości
    double wynik1 = wartoscFunkcji1(0.0);
    double wynik2 = wartoscFunkcji1(1.0);
    bool test2 = wynik1 != wynik2;
    wypiszWynikTestu("wartoscFunkcji1 - test 2", test2);
}

void testWartoscFunkcji2_1() {
    // Test 1: Druga funkcja testowa
    double wynik = wartoscFunkcji2(1.0);
    bool test1 = !isnan(wynik) && !isinf(wynik);
    wypiszWynikTestu("wartoscFunkcji2 - test 1", test1);
}

void testWartoscFunkcji2_2() {
    // Test 2: Różne wartości
    double wynik1 = wartoscFunkcji2(0.5);
    double wynik2 = wartoscFunkcji2(1.5);
    bool test2 = wynik1 != wynik2;
    wypiszWynikTestu("wartoscFunkcji2 - test 2", test2);
}

// =============================================================================
// TESTY DLA MODUŁU RÓWNANIA RÓŻNICZKOWE
// =============================================================================

void testObliczPochodna1() {
    cout << "\n=== TESTY RÓWNAŃ RÓŻNICZKOWYCH ===" << endl;

    // Test 1: Podstawowa wartość
    double wynik = obliczPochodna(25.0);
    bool test1 = !isnan(wynik) && !isinf(wynik);
    wypiszWynikTestu("obliczPochodna - test 1", test1);
}

void testObliczPochodna2() {
    // Test 2: Różne temperatury
    double wynik1 = obliczPochodna(20.0);
    double wynik2 = obliczPochodna(30.0);
    bool test2 = wynik1 != wynik2;
    wypiszWynikTestu("obliczPochodna - test 2", test2);
}

void testRozwiazanieDokladne1() {
    // Test 1: Wartość w t=0
    double wynik = rozwiazanieDokladne(0.0);
    bool test1 = !isnan(wynik) && !isinf(wynik);
    wypiszWynikTestu("rozwiazanieDokladne - test 1", test1);
}

void testRozwiazanieDokladne2() {
    // Test 2: Różne czasy
    double wynik1 = rozwiazanieDokladne(0.0);
    double wynik2 = rozwiazanieDokladne(1.0);
    bool test2 = wynik1 != wynik2;
    wypiszWynikTestu("rozwiazanieDokladne - test 2", test2);
}

void testMetodaEulera1() {
    // Test 1: Podstawowy test
    auto wynik = metodaEulera(0.1);
    bool test1 = !wynik.empty() && wynik.size() > 1;
    wypiszWynikTestu("metodaEulera - test 1 (rozmiar)", test1);
}

void testMetodaEulera2() {
    // Test 2: Pierwszy punkt powinien być warunkiem początkowym
    auto wynik = metodaEulera(0.1);
    bool test2 = !wynik.empty() && wynik[0].first == 0.0;
    wypiszWynikTestu("metodaEulera - test 2 (warunek początkowy)", test2);
}

void testObliczBladSredniKwadratowy1() {
    // Test 1: Pusty wektor
    vector<pair<double, double>> pusty;
    double blad = obliczBladSredniKwadratowy(pusty);
    bool test1 = blad == 0.0 || isnan(blad);
    wypiszWynikTestu("obliczBladSredniKwadratowy - test 1", test1);
}

void testObliczBladSredniKwadratowy2() {
    // Test 2: Jeden punkt
    vector<pair<double, double>> jeden = {{0.0, 1.0}};
    double blad = obliczBladSredniKwadratowy(jeden);
    bool test2 = blad >= 0;
    wypiszWynikTestu("obliczBladSredniKwadratowy - test 2", test2);
}

void testObliczZmianeTemperatury1() {
    // Test 1: Podstawowa wartość
    double wynik = oblicz_zmiane_temperatury(80.0);
    bool test1 = !isnan(wynik) && !isinf(wynik);
    wypiszWynikTestu("oblicz_zmiane_temperatury - test 1", test1);
}

void testObliczZmianeTemperatury2() {
    // Test 2: Temperatura otoczenia powinna dać zero
    double temp_otoczenia = 20.0;  // Zakładając, że to temperatura otoczenia
    double wynik = oblicz_zmiane_temperatury(temp_otoczenia);
    bool test2 = abs(wynik) < 1e-6;  // Powinno być blisko zera
       wypiszWynikTestu("oblicz_zmiane_temperatury - test 2", test2);
}

void testObliczMetodaEulera1() {
   // Test 1: Podstawowy test
   auto wynik = oblicz_metoda_Eulera(0.1);
   bool test1 = !wynik.empty() && wynik.size() > 1;
   wypiszWynikTestu("oblicz_metoda_Eulera - test 1 (rozmiar)", test1);
}

void testObliczMetodaEulera2() {
   // Test 2: Pierwszy element powinien być warunkiem początkowym
   auto wynik = oblicz_metoda_Eulera(0.1);
   bool test2 = !wynik.empty() && wynik[0] > 0;  // Zakładając dodatnią temperaturę początkową
   wypiszWynikTestu("oblicz_metoda_Eulera - test 2 (warunek początkowy)", test2);
}

void testObliczMetodaSrodka1() {
   // Test 1: Podstawowy test
   auto wynik = oblicz_metoda_srodka(0.1);
   bool test1 = !wynik.empty() && wynik.size() > 1;
   wypiszWynikTestu("oblicz_metoda_srodka - test 1 (rozmiar)", test1);
}

void testObliczMetodaSrodka2() {
   // Test 2: Porównanie z metodą Eulera (powinna być dokładniejsza)
   auto euler = oblicz_metoda_Eulera(0.1);
   auto srodek = oblicz_metoda_srodka(0.1);
   bool test2 = euler.size() == srodek.size() && !euler.empty() && !srodek.empty();
   wypiszWynikTestu("oblicz_metoda_srodka - test 2 (zgodność rozmiaru)", test2);
}

void testObliczMetodaHeuna1() {
   // Test 1: Podstawowy test
   auto wynik = oblicz_metoda_Heuna(0.1);
   bool test1 = !wynik.empty() && wynik.size() > 1;
   wypiszWynikTestu("oblicz_metoda_Heuna - test 1 (rozmiar)", test1);
}

void testObliczMetodaHeuna2() {
   // Test 2: Sprawdzenie, czy wyniki są rozsądne
   auto wynik = oblicz_metoda_Heuna(0.1);
   bool test2 = true;
   for(double temp : wynik) {
       if(isnan(temp) || isinf(temp)) {
           test2 = false;
           break;
       }
   }
   wypiszWynikTestu("oblicz_metoda_Heuna - test 2 (poprawność wartości)", test2);
}

void testObliczMetodaRK41() {
   // Test 1: Podstawowy test
   auto wynik = oblicz_metoda_RK4(0.1);
   bool test1 = !wynik.empty() && wynik.size() > 1;
   wypiszWynikTestu("oblicz_metoda_RK4 - test 1 (rozmiar)", test1);
}

void testObliczMetodaRK42() {
   // Test 2: RK4 powinna być najdokładniejsza
   auto euler = oblicz_metoda_Eulera(0.1);
   auto rk4 = oblicz_metoda_RK4(0.1);
   bool test2 = euler.size() == rk4.size() && !rk4.empty();
   wypiszWynikTestu("oblicz_metoda_RK4 - test 2 (zgodność rozmiaru)", test2);
}

void testTemperaturaDokladna1() {
   // Test 1: Wartość w t=0
   double wynik = temperatura_dokladna(0.0);
   bool test1 = !isnan(wynik) && !isinf(wynik) && wynik > 0;
   wypiszWynikTestu("temperatura_dokladna - test 1 (t=0)", test1);
}

void testTemperaturaDokladna2() {
   // Test 2: Temperatura powinna maleć z czasem (chłodzenie)
   double temp1 = temperatura_dokladna(0.0);
   double temp2 = temperatura_dokladna(1.0);
   bool test2 = temp1 > temp2;  // Temperatura powinna maleć
   wypiszWynikTestu("temperatura_dokladna - test 2 (chłodzenie)", test2);
}

void testObliczSredniBladKwadratowy1() {
   // Test 1: Pusty wektor
   vector<double> pusty;
   double blad = oblicz_sredni_blad_kwadratowy(pusty, 0.1);
   bool test1 = blad == 0.0 || isnan(blad);
   wypiszWynikTestu("oblicz_sredni_blad_kwadratowy - test 1 (pusty)", test1);
}

void testObliczSredniBladKwadratowy2() {
   // Test 2: Porównanie dokładności metod
   auto euler = oblicz_metoda_Eulera(0.1);
   auto rk4 = oblicz_metoda_RK4(0.1);

   if(!euler.empty() && !rk4.empty()) {
       double bladEuler = oblicz_sredni_blad_kwadratowy(euler, 0.1);
       double bladRK4 = oblicz_sredni_blad_kwadratowy(rk4, 0.1);
       bool test2 = bladRK4 <= bladEuler;  // RK4 powinna być dokładniejsza
       wypiszWynikTestu("oblicz_sredni_blad_kwadratowy - test 2 (porównanie)", test2);
   } else {
       wypiszWynikTestu("oblicz_sredni_blad_kwadratowy - test 2 (porównanie)", false);
   }
}