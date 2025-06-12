//
// Created by Szymon Ros on 11/06/2025.
//

#ifndef CALKOWANIE_NUMERYCZNE_H
#define CALKOWANIE_NUMERYCZNE_H
#include <vector>
#include <unordered_map>
using namespace std;

namespace biblioteka_numeryczna {
    /**
     * @brief Przykładowa funkcja trygonometryczna do testowania metod całkowania
     * @param x Argument funkcji
     * @return Wartość funkcji f(x) = sin(x) + cos(x)
     *
     * Przykład użycia:
     * double wartosc = FunkcjaTrygonometryczna(1.5);
     */
    double FunkcjaTrygonometryczna(double x);

    /**
     * @brief Oblicza wartość wielomianu metodą Hornera
     * @param wspolczynniki Współczynniki wielomianu (od a0 do an)
     * @param x Argument wielomianu
     * @return Wartość wielomianu w punkcie x
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1, 2, 3}; // 1 + 2x + 3x^2
     * double wartosc = Horner(wspolczynniki, 2.0);
     */
    double Horner(const vector<double>& wspolczynniki, double x);

    /**
     * @brief Całkowanie numeryczne wielomianu metodą prostokątów
     * @param wspolczynniki Współczynniki wielomianu
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param n Liczba podprzedziałów
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1, 2, 3};
     * double calka = MetodaProstokatow(wspolczynniki, 0.0, 1.0, 1000);
     */
    double MetodaProstokatow(const vector<double>& wspolczynniki, double poczatekPrzedzialu, double koniecPrzedzialu, int n = 100);

    /**
     * @brief Całkowanie numeryczne wielomianu metodą trapezów
     * @param wspolczynniki Współczynniki wielomianu
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param n Liczba podprzedziałów
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1, 2, 3};
     * double calka = MetodaTrapezow(wspolczynniki, 0.0, 1.0, 1000);
     */
    double MetodaTrapezow(const vector<double>& wspolczynniki, double poczatekPrzedzialu, double koniecPrzedzialu, int n = 100);

    /**
     * @brief Całkowanie numeryczne wielomianu metodą Simpsona
     * @param wspolczynniki Współczynniki wielomianu
     * @param a Początek przedziału całkowania
     * @param b Koniec przedziału całkowania
     * @param n Liczba podprzedziałów (musi być parzysta)
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1, 2, 3};
     * double calka = Simpson(wspolczynniki, 0.0, 1.0, 1000);
     */
    double Simpson(const vector<double>& wspolczynniki, double a, double b, int n = 100);

    /**
     * @brief Całkowanie funkcji trygonometrycznej metodą prostokątów
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param n Liczba podprzedziałów
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * double calka = MetodaProstokatowTryg(0.0, 3.14159, 1000);
     */
    double MetodaProstokatowTryg(double poczatekPrzedzialu, double koniecPrzedzialu, int n = 100);

    /**
     * @brief Całkowanie funkcji trygonometrycznej metodą trapezów
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param n Liczba podprzedziałów
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * double calka = MetodaTrapezowTryg(0.0, 3.14159, 1000);
     */
    double MetodaTrapezowTryg(double poczatekPrzedzialu, double koniecPrzedzialu, int n = 100);

    /**
     * @brief Całkowanie funkcji trygonometrycznej metodą Simpsona
     * @param a Początek przedziału całkowania
     * @param b Koniec przedziału całkowania
     * @param n Liczba podprzedziałów (musi być parzysta)
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * double calka = SimpsonTryg(0.0, 3.14159, 1000);
     */
    double SimpsonTryg(double a, double b, int n = 100);

    /**
     * @brief Oblicza dokładną wartość całki wielomianu analitycznie
     * @param wspolczynniki Współczynniki wielomianu
     * @param a Początek przedziału całkowania
     * @param b Koniec przedziału całkowania
     * @return Dokładna wartość całki
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1, 2, 3};
     * double calka = CalkaAnalitycznaWielomianu(wspolczynniki, 0.0, 1.0);
     */
    double CalkaAnalitycznaWielomianu(const vector<double>& wspolczynniki, double a, double b);

    /**
     * @brief Wyświetla analizę zbieżności metod całkowania dla wielomianu
     * @param wspolczynniki Współczynniki wielomianu
     * @param a Początek przedziału całkowania
     * @param b Koniec przedziału całkowania
     *
     * Porównuje różne metody całkowania z dokładnym wynikiem
     *
     * Przykład użycia:
     * vector<double> wspolczynniki = {1, 2, 3};
     * WyswietlZbieznosc(wspolczynniki, 0.0, 1.0);
     */
    void WyswietlZbieznosc(const vector<double>& wspolczynniki, double a, double b);

    /**
     * @brief Wyświetla analizę zbieżności metod całkowania dla funkcji trygonometrycznej
     * @param a Początek przedziału całkowania
     * @param b Koniec przedziału całkowania
     *
     * Porównuje różne metody całkowania funkcji trygonometrycznej
     *
     * Przykład użycia:
     * WyswietlZbieznoscTryg(0.0, 3.14159);
     */
    void WyswietlZbieznoscTryg(double a, double b);

    /**
     * @brief Wczytuje dane wielomianu z pliku
     * @param filename Nazwa pliku z danymi
     * @param degree Stopień wielomianu (zostanie ustawiony)
     * @param wspolczynniki Współczynniki wielomianu (zostaną wczytane)
     * @param a Początek przedziału (zostanie wczytany)
     * @param b Koniec przedziału (zostanie wczytany)
     * @return true jeśli wczytanie się powiodło, false w przeciwnym razie
     *
     * Przykład użycia:
     * int stopien;
     * vector<double> wspolczynniki;
     * double a, b;
     * bool sukces = wczytajDane("dane.txt", stopien, wspolczynniki, a, b);
     */
    bool wczytajDane(const string& filename, int& degree, vector<double>& wspolczynniki, double& a, double& b);

    /**
     * @brief Pierwsza funkcja testowa dla kwadratur Gaussa-Legendre'a
     * @param argument1 Argument funkcji
     * @return Wartość funkcji f1(x)
     *
     * Przykład użycia:
     * double wartosc = wartoscFunkcji1(1.5);
     */
    double wartoscFunkcji1(double argument1);

    /**
     * @brief Druga funkcja testowa dla kwadratur Gaussa-Legendre'a
     * @param argument2 Argument funkcji
     * @return Wartość funkcji f2(x)
     *
     * Przykład użycia:
     * double wartosc = wartoscFunkcji2(2.0);
     */
    double wartoscFunkcji2(double argument2);

    /**
     * @brief Trzecia funkcja testowa dla kwadratur Gaussa-Legendre'a
     * @param argument3 Argument funkcji
     * @return Wartość funkcji f3(x)
     *
     * Przykład użycia:
     * double wartosc = wartoscFunkcji3(0.5);
     */
    double wartoscFunkcji3(double argument3);

    /**
     * @brief Czwarta funkcja testowa dla kwadratur Gaussa-Legendre'a
     * @param argument4 Argument funkcji
     * @return Wartość funkcji f4(x)
     *
     * Przykład użycia:
     * double wartosc = wartoscFunkcji4(1.0);
     */
    double wartoscFunkcji4(double argument4);

    /**
     * @brief Kwadratura Gaussa-Legendre'a
     * @param punkty Mapa punktów kwadratury dla różnych rzędów
     * @param wagi Mapa wag kwadratury dla różnych rzędów
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param liczbaPunktow Liczba punktów kwadratury
     * @param funkcja Wskaźnik do funkcji do całkowania
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * unordered_map<int, vector<double>> punkty, wagi;
     * double calka = kwadraturaGL(punkty, wagi, 0.0, 1.0, 5, wartoscFunkcji1);
     */
    double kwadraturaGL(unordered_map<int, vector<double>>& punkty, unordered_map<int, vector<double>>& wagi,
                    double poczatekPrzedzialu, double koniecPrzedzialu, int liczbaPunktow,
                    double (*funkcja)(double));

    /**
     * @brief Kwadratura Gaussa-Legendre'a z podziałem przedziału
     * @param punkty Mapa punktów kwadratury dla różnych rzędów
     * @param wagi Mapa wag kwadratury dla różnych rzędów
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param liczbaPunktow Liczba punktów kwadratury
     * @param funkcja Wskaźnik do funkcji do całkowania
     * @param liczbaPodzialow Liczba podziałów przedziału
     * @return Przybliżona wartość całki
     *
     * Przykład użycia:
     * unordered_map<int, vector<double>> punkty, wagi;
     * double calka = kwadraturaGLPodzial(punkty, wagi, 0.0, 1.0, 5, wartoscFunkcji1, 10);
     */
    double kwadraturaGLPodzial(unordered_map<int, vector<double>>& punkty, unordered_map<int, vector<double>>& wagi,
                          double poczatekPrzedzialu, double koniecPrzedzialu, int liczbaPunktow,
                          double (*funkcja)(double), int liczbaPodzialow);

    /**
     * @brief Analizuje zbieżność kwadratur Gaussa-Legendre'a
     * @param punkty Mapa punktów kwadratury dla różnych rzędów
     * @param wagi Mapa wag kwadratury dla różnych rzędów
     * @param poczatekPrzedzialu Początek przedziału całkowania
     * @param koniecPrzedzialu Koniec przedziału całkowania
     * @param dokladnaWartosc Dokładna wartość całki
     * @param funkcja Wskaźnik do funkcji do całkowania
     * @param nazwaFunkcji Nazwa funkcji (do wyświetlania)
     *
     * Przykład użycia:
     * unordered_map<int, vector<double>> punkty, wagi;
     * analizaZbieznosci(punkty, wagi, 0.0, 1.0, 2.0, wartoscFunkcji1, "f1(x)");
     */
    void analizaZbieznosci(unordered_map<int, vector<double>>& punkty, unordered_map<int, vector<double>>& wagi,
                      double poczatekPrzedzialu, double koniecPrzedzialu, double dokladnaWartosc,
                      double (*funkcja)(double), const string& nazwaFunkcji);
}
    #endif //CALKOWANIE_NUMERYCZNE_H