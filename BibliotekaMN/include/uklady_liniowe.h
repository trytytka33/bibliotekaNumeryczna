//
// Created by Szymon Ros on 10/06/2025.
//

#ifndef UKLADY_LINIOWE_H
#define UKLADY_LINIOWE_H

#include <vector>
using namespace std;

namespace biblioteka_numeryczna {
    /**
     * @brief Rozwiązuje układ równań liniowych metodą eliminacji Gaussa
     * @param macierz_A Macierz współczynników (zostanie zmodyfikowana)
     * @param wektor_b Wektor wyrazów wolnych (zostanie zmodyfikowany)
     * @return Wektor rozwiązań
     *
     * Przykład użycia:
     * vector<vector<double>> A = {{2, 1}, {1, 3}};
     * vector<double> b = {3, 4};
     * vector<double> rozwiazanie = eliminacjaGaussa(A, b);
     */
    vector<double> eliminacjaGaussa(vector<vector<double>>& macierz_A,
                                   vector<double>& wektor_b);

    /**
     * @brief Rozwiązuje układ równań metodą Gaussa-Seidela
     * @param macierz_A Macierz współczynników
     * @param wektor_b Wektor wyrazów wolnych
     * @param przyblizenie_poczatkowe Przybliżenie początkowe
     * @param tolerancja Tolerancja błędu
     * @param max_iteracji Maksymalna liczba iteracji
     * @return Wektor rozwiązań
     *
     * Przykład użycia:
     * vector<vector<double>> A = {{4, 1}, {1, 3}};
     * vector<double> b = {5, 4};
     * vector<double> x0 = {0, 0};
     * vector<double> rozwiazanie = metodaGaussaSeidela(A, b, x0);
     */
    vector<double> metodaGaussaSeidela(const vector<vector<double>>& macierz_A,
                                      const vector<double>& wektor_b,
                                      const vector<double>& przyblizenie_poczatkowe,
                                      double tolerancja = 1e-6,
                                      int max_iteracji = 1000);

    /**
     * @brief Wczytuje dane układu równań z pliku
     * @param nazwaPliku Nazwa pliku z danymi
     * @param A Macierz współczynników (zostanie wypełniona)
     * @param B Wektor wyrazów wolnych (zostanie wypełniony)
     * @param N Rozmiar układu (zostanie ustawiony)
     *
     * Przykład użycia:
     * vector<vector<double>> A;
     * vector<double> B;
     * int N;
     * wczytajDane("dane.txt", A, B, N);
     */
    void wczytajDane(const string& nazwaPliku,
                    vector<vector<double>>& A,
                    vector<double>& B, int& N);

    /**
     * @brief Tworzy macierz dopełnioną [A|b] z macierzy A i wektora b
     * @param A Macierz współczynników
     * @param B Wektor wyrazów wolnych
     * @param macierzDopelniona Wynikowa macierz dopełniona
     *
     * Przykład użycia:
     * vector<vector<double>> A = {{2, 1}, {1, 3}};
     * vector<double> B = {3, 4};
     * vector<vector<double>> macierzDopelniona;
     * utworzMacierzDopelniona(A, B, macierzDopelniona);
     */
    void utworzMacierzDopelniona(const vector<vector<double>>& A,
                                const vector<double>& B,
                                vector<vector<double>>& macierzDopelniona);

    /**
     * @brief Wyświetla macierz w czytelnej formie
     * @param macierz Macierz do wyświetlenia
     * @param nazwa Nazwa macierzy (do wyświetlenia jako nagłówek)
     *
     * Przykład użycia:
     * vector<vector<double>> A = {{2, 1}, {1, 3}};
     * wyswietlMacierz(A, "Macierz A");
     */
    void wyswietlMacierz(const vector<vector<double>>& macierz,
                        const string& nazwa);

    /**
     * @brief Wyświetla wektor w czytelnej formie
     * @param wektor Wektor do wyświetlenia
     * @param nazwa Nazwa wektora (do wyświetlenia jako nagłówek)
     *
     * Przykład użycia:
     * vector<double> x = {1.5, 2.3};
     * wyswietlWektor(x, "Rozwiązanie");
     */
    void wyswietlWektor(const vector<double>& wektor,
                        const string& nazwa);

    /**
     * @brief Sprawdza poprawność rozwiązania przez podstawienie do równania Ax = b
     * @param A Macierz współczynników
     * @param B Wektor wyrazów wolnych
     * @param X Wektor rozwiązań do sprawdzenia
     *
     * Przykład użycia:
     * vector<vector<double>> A = {{2, 1}, {1, 3}};
     * vector<double> B = {3, 4};
     * vector<double> X = {1, 1};
     * sprawdzWynik(A, B, X);
     */
    void sprawdzWynik(const vector<vector<double>>& A,
                        const vector<double>& B,
                        const vector<double>& X);

    /**
     * @brief Sprowadza macierz dopełnioną do postaci schodkowej
     * @param macierzDopelniona Macierz dopełniona (zostanie zmodyfikowana)
     *
     * Przykład użycia:
     * vector<vector<double>> macierzDopelniona = {{2, 1, 3}, {1, 3, 4}};
     * sprowadzanieDoPostaciSchodkowej(macierzDopelniona);
     */
    void sprowadzanieDoPostaciSchodkowej(vector<vector<double>>& macierzDopelniona);

    /**
     * @brief Rozwiązuje układ równań z macierzy w postaci schodkowej metodą podstawiania wstecznego
     * @param macierzSchodkowa Macierz w postaci schodkowej
     * @return Wektor rozwiązań
     *
     * Przykład użycia:
     * vector<vector<double>> macierzSchodkowa = {{2, 1, 3}, {0, 2.5, 2.5}};
     * vector<double> rozwiazanie = rozwiazUklad(macierzSchodkowa);
     */
    vector<double> rozwiazUklad(const vector<vector<double>>& macierzSchodkowa);
    /**
 * @brief Przeprowadza rozkład LU z pivotingiem
 * @param A Macierz współczynników
 * @param L Macierz dolna trójkątna (zostanie wypełniona)
 * @param U Macierz górna trójkątna (zostanie wypełniona)
 * @param P Wektor permutacji (zostanie wypełniony)
 *
 * Przykład użycia:
 * vector<vector<double>> A = {{2, 1}, {1, 3}};
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
 */
    vector<double> permutujWektor(const vector<double>& b, const vector<int>& P);

/**
 * @brief Rozwiązuje układ Ly = b (podstawianie w przód)
 * @param L Macierz dolna trójkątna
 * @param b Wektor wyrazów wolnych
 * @return Wektor rozwiązań y
 */
    vector<double> rozwiazLy_b(const vector<vector<double>>& L, const vector<double>& b);

/**
 * @brief Rozwiązuje układ Ux = y (podstawianie wsteczne)
 * @param U Macierz górna trójkątna
 * @param y Wektor wyrazów wolnych
 * @return Wektor rozwiązań x
 */
    vector<double> rozwiazUx_y(const vector<vector<double>>& U, const vector<double>& y);

/**
 * @brief Sprawdza poprawność rozwiązania i zwraca maksymalny błąd
 * @param A Macierz współczynników
 * @param x Wektor rozwiązań
 * @param b Wektor wyrazów wolnych
 * @return Maksymalny błąd
 */
    double sprawdzPoprawnosc(const vector<vector<double>>& A,
                        const vector<double>& x,
                        const vector<double>& b);

/**
 * @brief Rozwiązuje układ równań metodą rozkładu LU z pivotingiem
 * @param A Macierz współczynników
 * @param b Wektor wyrazów wolnych
 * @return Para (y, x) gdzie y to rozwiązanie Ly=Pb, x to końcowe rozwiązanie
 */
    pair<vector<double>, vector<double>> rozwiazUkladLU(const vector<vector<double>>& A,
                                                   const vector<double>& b);

/**
 * @brief Testuje rozwiązanie układu równań z pliku
 * @param nazwaPliku Nazwa pliku z danymi
 */
    void testUkladu(const string& nazwaPliku);
}

#endif // UKLADY_LINIOWE_H