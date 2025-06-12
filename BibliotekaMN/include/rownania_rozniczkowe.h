//
// Created by Szymon Ros on 11/06/2025.
//

#ifndef ROWNANIA_ROZNICZKOWE_H
#define ROWNANIA_ROZNICZKOWE_H
#include <iostream>
#include <vector>
using namespace std;

namespace biblioteka_numeryczna {
    /**
     * @brief Oblicza pochodną funkcji temperatury (przykładowe równanie różniczkowe)
     * @param T Aktualna temperatura
     * @return Wartość pochodnej dT/dt
     *
     * Implementuje przykładowe równanie różniczkowe opisujące zmianę temperatury
     *
     * Przykład użycia:
     * double pochodna = obliczPochodna(25.0);
     */
    double obliczPochodna(double T);

    /**
     * @brief Rozwiązanie dokładne równania różniczkowego (jeśli znane)
     * @param t Czas
     * @return Dokładna wartość funkcji w czasie t
     *
     * Przykład użycia:
     * double dokladna_wartosc = rozwiazanieDokladne(1.5);
     */
    double rozwiazanieDokladne(double t);

    /**
     * @brief Rozwiązuje równanie różniczkowe metodą Eulera
     * @param krok Krok całkowania
     * @return Wektor par (czas, wartość) reprezentujących rozwiązanie numeryczne
     *
     * Metoda Eulera: y_{n+1} = y_n + h * f(t_n, y_n)
     *
     * Przykład użycia:
     * vector<pair<double, double>> rozwiazanie = metodaEulera(0.1);
     */
    vector<pair<double, double>> metodaEulera(double krok);

    /**
     * @brief Oblicza średni błąd kwadratowy rozwiązania numerycznego
     * @param rozwiazanieNumeryczne Wektor par (czas, wartość) z rozwiązania numerycznego
     * @return Średni błąd kwadratowy względem rozwiązania dokładnego
     *
     * Przykład użycia:
     * auto rozwiazanie = metodaEulera(0.1);
     * double blad = obliczBladSredniKwadratowy(rozwiazanie);
     */
    double obliczBladSredniKwadratowy(const vector<pair<double, double>>& rozwiazanieNumeryczne);

    /**
     * @brief Oblicza zmianę temperatury zgodnie z prawem chłodzenia Newtona
     * @param T Aktualna temperatura
     * @return Pochodna temperatury dT/dt
     *
     * Implementuje prawo chłodzenia Newtona: dT/dt = -k(T - T_otoczenia)
     *
     * Przykład użycia:
     * double zmiana = oblicz_zmiane_temperatury(80.0);
     */
    double oblicz_zmiane_temperatury(double T);

    /**
     * @brief Rozwiązuje równanie różniczkowe metodą Eulera (uproszczona wersja)
     * @param krok Krok całkowania
     * @return Wektor wartości temperatury w kolejnych krokach czasowych
     *
     * Przykład użycia:
     * vector<double> rozwiazanie = oblicz_metoda_Eulera(0.1);
     */
    vector<double> oblicz_metoda_Eulera(double krok);

    /**
     * @brief Rozwiązuje równanie różniczkowe metodą punktu środkowego
     * @param krok Krok całkowania
     * @return Wektor wartości temperatury w kolejnych krokach czasowych
     *
     * Metoda punktu środkowego: y_{n+1} = y_n + h * f(t_n + h/2, y_n + h*f(t_n, y_n)/2)
     *
     * Przykład użycia:
     * vector<double> rozwiazanie = oblicz_metoda_srodka(0.1);
     */
    vector<double> oblicz_metoda_srodka(double krok);

    /**
     * @brief Rozwiązuje równanie różniczkowe metodą Heuna
     * @param krok Krok całkowania
     * @return Wektor wartości temperatury w kolejnych krokach czasowych
     *
     * Metoda Heuna (ulepszona metoda Eulera):
     * k1 = h * f(t_n, y_n)
     * k2 = h * f(t_n + h, y_n + k1)
     * y_{n+1} = y_n + (k1 + k2)/2
     *
     * Przykład użycia:
     * vector<double> rozwiazanie = oblicz_metoda_Heuna(0.1);
     */
    vector<double> oblicz_metoda_Heuna(double krok);

    /**
     * @brief Rozwiązuje równanie różniczkowe metodą Runge-Kutta 4. rzędu
     * @param krok Krok całkowania
     * @return Wektor wartości temperatury w kolejnych krokach czasowych
     *
     * Klasyczna metoda RK4:
     * k1 = h * f(t_n, y_n)
     * k2 = h * f(t_n + h/2, y_n + k1/2)
     * k3 = h * f(t_n + h/2, y_n + k2/2)
     * k4 = h * f(t_n + h, y_n + k3)
     * y_{n+1} = y_n + (k1 + 2*k2 + 2*k3 + k4)/6
     *
     * Przykład użycia:
     * vector<double> rozwiazanie = oblicz_metoda_RK4(0.1);
     */
    vector<double> oblicz_metoda_RK4(double krok);

    /**
     * @brief Dokładne rozwiązanie równania chłodzenia Newtona
     * @param t Czas
     * @return Dokładna wartość temperatury w czasie t
     *
     * Analityczne rozwiązanie: T(t) = T_otoczenia + (T_0 - T_otoczenia) * exp(-k*t)
     *
     * Przykład użycia:
     * double dokladna_temp = temperatura_dokladna(2.0);
     */
    double temperatura_dokladna(double t);

    /**
     * @brief Oblicza średni błąd kwadratowy metody numerycznej
     * @param wyniki Wektor wyników numerycznych
     * @param krok Krok czasowy użyty w metodzie
     * @return Średni błąd kwadratowy względem rozwiązania dokładnego
     *
     * Przykład użycia:
     * vector<double> wyniki = oblicz_metoda_RK4(0.1);
     * double blad = oblicz_sredni_blad_kwadratowy(wyniki, 0.1);
     */
    double oblicz_sredni_blad_kwadratowy(const vector<double>& wyniki, double krok);
}
#endif //ROWNANIA_ROZNICZKOWE_H