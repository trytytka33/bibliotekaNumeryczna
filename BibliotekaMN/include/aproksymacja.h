//
// Created by Szymon Ros on 11/06/2025.
//

#ifndef APROKSYMACJA_H
#define APROKSYMACJA_H
#include <vector>
using namespace std;

namespace biblioteka_numeryczna {
    /**
     * @brief Oblicza iloczyn skalarny dwóch elementów bazy w przestrzeni L2
     * @param pierwszyElement Indeks pierwszego elementu bazy (0 dla stałej, 1 dla x, 2 dla x^2, itd.)
     * @param drugiElement Indeks drugiego elementu bazy
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param n Liczba punktów do całkowania numerycznego
     * @return Wartość iloczynu skalarnego <φ_i, φ_j>
     *
     * Przykład użycia:
     * double iloczyn = iloczynSkalarnyElementowZBazy(1, 2, 0.0, 1.0, 1000);
     * // Oblicza iloczyn skalarny x * x^2 na przedziale [0,1]
     */
    double iloczynSkalarnyElementowZBazy(int pierwszyElement, int drugiElement,
                                        const double poczatekPrzedzialu,
                                        const double koniecPrzedzialu, int n);

    /**
     * @brief Funkcja która ma być aproksymowana
     * @param x Argument funkcji
     * @return Wartość funkcji f(x)
     *
     * Przykład użycia:
     * double wartosc = funkcjaAproksymowana(1.5);
     */
    double funkcjaAproksymowana(double x);

    /**
     * @brief Oblicza iloczyn skalarny elementu bazy z funkcją aproksymowaną
     * @param elementZBazy Indeks elementu bazy (0 dla stałej, 1 dla x, itd.)
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param n Liczba punktów do całkowania numerycznego
     * @return Wartość iloczynu skalarnego <φ_i, f>
     *
     * Przykład użycia:
     * double iloczyn = iloczynSkalarnyElementuZBazyIFunkcji(2, -1.0, 1.0, 1000);
     * // Oblicza iloczyn skalarny x^2 z funkcją f na przedziale [-1,1]
     */
    double iloczynSkalarnyElementuZBazyIFunkcji(int elementZBazy,
                                               const double poczatekPrzedzialu,
                                               const double koniecPrzedzialu, int n);

    /**
     * @brief Przeprowadza rozkład LU z pivotingiem dla macierzy Grama
     * @param A Macierz współczynników (macierz Grama)
     * @param L Macierz dolna trójkątna (zostanie wypełniona)
     * @param U Macierz górna trójkątna (zostanie wypełniona)
     * @param P Wektor permutacji (zostanie wypełniony)
     *
     * Przykład użycia:
     * vector<vector<double>> A = {{1, 0.5}, {0.5, 0.33}};
     * vector<vector<double>> L, U;
     * vector<int> P;
     * rozkladLU_zPivotingiem(A, L, U, P);
     */
    void rozkladLU_zPivotingiem(const vector<vector<double>>& A,
                               vector<vector<double>>& L,
                               vector<vector<double>>& U,
                               vector<int>& P);

    /**
     * @brief Permutuje wektor zgodnie z wektorem permutacji
     * @param b Wektor do permutacji
     * @param P Wektor permutacji
     * @return Permutowany wektor
     *
     * Przykład użycia:
     * vector<double> b = {1, 2, 3};
     * vector<int> P = {2, 0, 1};
     * vector<double> permutowany = permutujWektor(b, P);
     */
    vector<double> permutujWektor(const vector<double>& b, const vector<int>& P);

    /**
     * @brief Rozwiązuje układ Ly = b (podstawianie w przód)
     * @param L Macierz dolna trójkątna
     * @param b Wektor wyrazów wolnych
     * @return Wektor rozwiązań y
     *
     * Przykład użycia:
     * vector<vector<double>> L = {{1, 0}, {0.5, 1}};
     * vector<double> b = {1, 2};
     * vector<double> y = rozwiazLy_b(L, b);
     */
    vector<double> rozwiazLy_b(const vector<vector<double>>& L, const vector<double>& b);

    /**
     * @brief Rozwiązuje układ Ux = y (podstawianie wsteczne)
     * @param U Macierz górna trójkątna
     * @param y Wektor wyrazów wolnych
     * @return Wektor rozwiązań x
     *
     * Przykład użycia:
     * vector<vector<double>> U = {{2, 1}, {0, 1.5}};
     * vector<double> y = {3, 2};
     * vector<double> x = rozwiazUx_y(U, y);
     */
    vector<double> rozwiazUx_y(const vector<vector<double>>& U, const vector<double>& y);

    /**
     * @brief Sprawdza poprawność rozwiązania i zwraca maksymalny błąd
     * @param A Macierz współczynników
     * @param x Wektor rozwiązań
     * @param b Wektor wyrazów wolnych
     * @return Maksymalny błąd residuum
     *
     * Przykład użycia:
     * vector<vector<double>> A = {{2, 1}, {1, 3}};
     * vector<double> x = {1, 1};
     * vector<double> b = {3, 4};
     * double blad = sprawdzPoprawnosc(A, x, b);
     */
    double sprawdzPoprawnosc(const vector<vector<double>>& A,
                            const vector<double>& x,
                            const vector<double>& b);

    /**
     * @brief Rozwiązuje układ równań metodą rozkładu LU z pivotingiem
     * @param A Macierz współczynników (macierz Grama)
     * @param b Wektor wyrazów wolnych (iloczyny skalarne z funkcją)
     * @return Para (y, x) gdzie y to rozwiązanie Ly=Pb, x to współczynniki aproksymacji
     *
     * Przykład użycia:
     * vector<vector<double>> gram = {{1, 0.5}, {0.5, 0.33}};
     * vector<double> rhs = {2, 1};
     * auto [y, wspolczynniki] = rozwiazUkladLU(gram, rhs);
     */
    pair<vector<double>, vector<double>> rozwiazUkladLU(const vector<vector<double>>& A,
                                                       const vector<double>& b);

    /**
     * @brief Oblicza błąd aproksymacji metodą najmniejszych kwadratów
     * @param wspolczynniki Współczynniki aproksymacji wielomianowej
     * @param poczatekPrzedzialu Początek przedziału
     * @param koniecPrzedzialu Koniec przedziału
     * @param liczbaPunktow Liczba punktów do obliczenia błędu
     * @return Średni błąd kwadratowy aproksymacji
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1.0, 2.0, -0.5};
     * double blad = bladAproksymacji(wspolczynniki, 0.0, 1.0, 1000);
     */
    double bladAproksymacji(const vector<double>& wspolczynniki,
                           double poczatekPrzedzialu,
                           double koniecPrzedzialu,
                           int liczbaPunktow);

    /**
     * @brief Wyświetla porównanie funkcji oryginalnej z aproksymacją
     * @param wspolczynniki Współczynniki aproksymacji wielomianowej
     * @param poczatekPrzedzialu Początek przedziału wyświetlania
     * @param koniecPrzedzialu Koniec przedziału wyświetlania
     * @param liczbaPunktow Liczba punktów do wyświetlenia
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1.0, 2.0, -0.5};
     * pokazAproksymacje(wspolczynniki, -1.0, 1.0, 20);
     */
    void pokazAproksymacje(const vector<double>& wspolczynniki,
                          double poczatekPrzedzialu,
                          double koniecPrzedzialu,
                          int liczbaPunktow);
}

#endif //APROKSYMACJA_H