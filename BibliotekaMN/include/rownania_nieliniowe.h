//
// Created by Szymon Ros on 11/06/2025.
//

#ifndef ROWNANIA_NIELINIOWE_H
#define ROWNANIA_NIELINIOWE_H
#include <iostream>
#include <vector>
using namespace std;

namespace biblioteka_numeryczna {
    /**
     * @brief Pierwsza funkcja testowa f1(x)
     * @param x Argument funkcji
     * @return Wartość funkcji f1(x)
     *
     * Przykład użycia:
     * double wartosc = f1(2.5);
     */
    double f1(double x);

    /**
     * @brief Pochodna pierwszej funkcji testowej f1'(x)
     * @param x Argument funkcji
     * @return Wartość pochodnej f1'(x)
     *
     * Przykład użycia:
     * double pochodna = pochodnaF1(2.5);
     */
    double pochodnaF1(double x);

    /**
     * @brief Druga funkcja testowa f2(x)
     * @param x Argument funkcji
     * @return Wartość funkcji f2(x)
     *
     * Przykład użycia:
     * double wartosc = f2(1.0);
     */
    double f2(double x);

    /**
     * @brief Pochodna drugiej funkcji testowej f2'(x)
     * @param x Argument funkcji
     * @return Wartość pochodnej f2'(x)
     *
     * Przykład użycia:
     * double pochodna = pochodnaF2(1.0);
     */
    double pochodnaF2(double x);

    /**
     * @brief Trzecia funkcja testowa f3(x)
     * @param x Argument funkcji
     * @return Wartość funkcji f3(x)
     *
     * Przykład użycia:
     * double wartosc = f3(0.5);
     */
    double f3(double x);

    /**
     * @brief Pochodna trzeciej funkcji testowej f3'(x)
     * @param x Argument funkcji
     * @return Wartość pochodnej f3'(x)
     *
     * Przykład użycia:
     * double pochodna = pochodnaF3(0.5);
     */
    double pochodnaF3(double x);

    /**
     * @brief Znajduje pierwiastek równania metodą bisekcji
     * @param a Lewy kraniec przedziału poszukiwań
     * @param b Prawy kraniec przedziału poszukiwań
     * @param funkcja Wskaźnik do funkcji, której pierwiastka szukamy
     * @return Przybliżona wartość pierwiastka
     *
     * Warunki wstępne: f(a) * f(b) < 0 (funkcja zmienia znak na przedziale)
     *
     * Przykład użycia:
     * double pierwiastek = metodaBisekcji(-2.0, 2.0, f1);
     */
    double metodaBisekcji(double a, double b, double (*funkcja)(double));

    /**
     * @brief Znajduje pierwiastek równania metodą Newtona-Raphsona
     * @param x0 Przybliżenie początkowe
     * @param funkcja Wskaźnik do funkcji, której pierwiastka szukamy
     * @param pochodna Wskaźnik do pochodnej funkcji
     * @return Przybliżona wartość pierwiastka
     *
     * Warunki wstępne: pochodna nie może być zero w otoczeniu pierwiastka
     *
     * Przykład użycia:
     * double pierwiastek = metodaNewtona(1.0, f1, pochodnaF1);
     */
    double metodaNewtona(double x0, double (*funkcja)(double), double (*pochodna)(double));

    /**
     * @brief Znajduje pierwiastek równania metodą siecznych
     * @param x0 Pierwsze przybliżenie początkowe
     * @param x1 Drugie przybliżenie początkowe
     * @param funkcja Wskaźnik do funkcji, której pierwiastka szukamy
     * @return Przybliżona wartość pierwiastka
     *
     * Przykład użycia:
     * double pierwiastek = metodaSiecznych(0.0, 1.0, f1);
     */
    double metodaSiecznych(double x0, double x1, double (*funkcja)(double));

    /**
     * @brief Znajduje wszystkie pierwiastki funkcji na danym przedziale
     * @param a Początek przedziału poszukiwań
     * @param b Koniec przedziału poszukiwań
     * @param funkcja Wskaźnik do funkcji, której pierwiastków szukamy
     * @param pochodna Wskaźnik do pochodnej funkcji
     * @param krok Krok skanowania przedziału
     * @return Wektor znalezionych pierwiastków
     *
     * Metoda skanuje przedział [a,b] z krokiem i stosuje metodę Newtona
     * dla każdego punktu startowego gdzie funkcja zmienia znak
     *
     * Przykład użycia:
     * vector<double> pierwiastki = znajdzWszystkiePierwiastki(-5.0, 5.0, f1, pochodnaF1, 0.1);
     */
    vector<double> znajdzWszystkiePierwiastki(double a, double b,
                                             double (*funkcja)(double),
                                             double (*pochodna)(double),
                                             double krok);

    /**
     * @brief Wypisuje wynik obliczeń wraz z weryfikacją
     * @param nazwaMetody Nazwa użytej metody numerycznej
     * @param wynik Znaleziony pierwiastek
     * @param funkcja Wskaźnik do funkcji (do weryfikacji)
     *
     * Funkcja wypisuje nazwę metody, znaleziony pierwiastek
     * oraz wartość funkcji w tym punkcie (powinna być bliska zeru)
     *
     * Przykład użycia:
     * double pierwiastek = metodaNewtona(1.0, f1, pochodnaF1);
     * wypiszWynik("Metoda Newtona", pierwiastek, f1);
     */
    void wypiszWynik(const string& nazwaMetody, double wynik, double (*funkcja)(double));
}

#endif //ROWNANIA_NIELINIOWE_H