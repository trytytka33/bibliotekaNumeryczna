//
// Created by Szymon Ros on 10/06/2025.
//
#include "../include/uklady_liniowe.h"
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace biblioteka_numeryczna {
    vector<double> eliminacjaGaussa(vector<vector<double>>& macierz_A,
                                   vector<double>& wektor_b) {
        int n = macierz_A.size();
        if (n == 0 || wektor_b.size() != n) {
            throw invalid_argument("Nieprawidłowe wymiary macierzy lub wektora");
        }

        // Eliminacja w przód
        for (int i = 0; i < n; i++) {
            // Znajdź element główny
            int max_wiersz = i;
            for (int k = i + 1; k < n; k++) {
                if (abs(macierz_A[k][i]) > abs(macierz_A[max_wiersz][i])) {
                    max_wiersz = k;
                }
            }

            // Zamień wiersze
            swap(macierz_A[i], macierz_A[max_wiersz]);
            swap(wektor_b[i], wektor_b[max_wiersz]);

            // Sprawdź czy element główny nie jest zerem
            if (abs(macierz_A[i][i]) < 1e-10) {
                throw runtime_error("Macierz jest osobliwa");
            }

            // Eliminuj kolumnę
            for (int k = i + 1; k < n; k++) {
                double wspolczynnik = macierz_A[k][i] / macierz_A[i][i];
                for (int j = i; j < n; j++) {
                    macierz_A[k][j] -= wspolczynnik * macierz_A[i][j];
                }
                wektor_b[k] -= wspolczynnik * wektor_b[i];
            }
        }

        // Podstawienie wsteczne
        vector<double> rozwiazanie(n);
        for (int i = n - 1; i >= 0; i--) {
            rozwiazanie[i] = wektor_b[i];
            for (int j = i + 1; j < n; j++) {
                rozwiazanie[i] -= macierz_A[i][j] * rozwiazanie[j];
            }
            rozwiazanie[i] /= macierz_A[i][i];
        }

        return rozwiazanie;
    }

    vector<double> metodaGaussaSeidela(const vector<vector<double>>& macierz_A,
                                      const vector<double>& wektor_b,
                                      const vector<double>& przyblizenie_poczatkowe,
                                      double tolerancja,
                                      int max_iteracji) {
        int n = macierz_A.size();
        if (n == 0 || wektor_b.size() != n || przyblizenie_poczatkowe.size() != n) {
            throw invalid_argument("Nieprawidłowe wymiary");
        }

        vector<double> x = przyblizenie_poczatkowe;
        vector<double> x_nowe(n);

        for (int iteracja = 0; iteracja < max_iteracji; iteracja++) {
            for (int i = 0; i < n; i++) {
                double suma = 0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        suma += macierz_A[i][j] * x[j];
                    }
                }
                x_nowe[i] = (wektor_b[i] - suma) / macierz_A[i][i];
                x[i] = x_nowe[i]; // Aktualizuj od razu (różnica od Jacobiego)
            }

            // Sprawdź zbieżność
            double norma_bledu = 0;
            for (int i = 0; i < n; i++) {
                norma_bledu += abs(x_nowe[i] - x[i]);
            }
            if (norma_bledu < tolerancja) {
                break;
            }
        }

        return x;
    }
    void wczytajDane(const string& nazwaPliku, vector<vector<double>>& A, vector<double>& B, int& N) {
        ifstream plik(nazwaPliku);
        if (!plik) {
            cerr << "Nie można otworzyć pliku: " << nazwaPliku << endl;
            return;
        }
        if (!plik) {
            string sciezkaAlternatywna = "../" + nazwaPliku;
            plik.open(sciezkaAlternatywna);
        }

        string linia;
        getline(plik, linia);

        getline(plik, linia);               //liczba niewiadomych
        stringstream ss(linia);
        string tmp;
        ss >> tmp >> N;

        A.resize(N, vector<double>(N, 0));
        B.resize(N, 0);

        getline(plik, linia);

        getline(plik, linia);                       //wektor wyrazow wolnych
        stringstream ssB(linia);
        for (int i = 0; i < N; i++) {
            ssB >> B[i];
        }

        getline(plik, linia);

        for (int i = 0; i < N; i++) {                       //wektor wspolczynnikow
            getline(plik, linia);
            stringstream ssA(linia);
            for (int j = 0; j < N; j++) {
                ssA >> A[i][j];
            }
        }

        plik.close();
    }
       void rozkladLU_zPivotingiem(const vector<vector<double>>& A,
                               vector<vector<double>>& L,
                               vector<vector<double>>& U,
                               vector<int>& P) {
        int n = A.size();

        U.resize(n, vector<double>(n, 0.0));
        L.resize(n, vector<double>(n, 0.0));
        P.resize(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                U[i][j] = A[i][j];
            }
            P[i] = i;
            L[i][i] = 1.0;
        }

        for (int k = 0; k < n; k++) {
            cout << "Iteracja " << k+1 << endl;

            double max_val = 0.0;
            int max_idx = k;

            for (int i = k; i < n; i++) {
                double abs_val = fabs(U[i][k]);
                if (abs_val > max_val) {
                    max_val = abs_val;
                    max_idx = i;
                }
            }
            if (max_val < 1e-10) {
                cerr << "Macierz jest osobliwa!" << endl;
            }

            if (max_idx != k) {
                cout << "Zamiana wierszy: " << k+1 << " i " << max_idx+1 << endl;

                swap(P[k], P[max_idx]);
                swap(U[k], U[max_idx]);

                for (int j = 0; j < k; j++) {
                    swap(L[k][j], L[max_idx][j]);
                }
            }

            for (int i = k + 1; i < n; i++) {
                L[i][k] = U[i][k] / U[k][k];

                for (int j = k; j < n; j++) {
                    U[i][j] = U[i][j] - L[i][k] * U[k][j];
                }
            }

            cout << "Macierz U po iteracji " << k+1 << ":" << endl;
            wyswietlMacierz(U, "U (aktualna)");

            cout << "Macierz L po iteracji " << k+1 << ":" << endl;
            wyswietlMacierz(L, "L (aktualna)");
        }

        wyswietlMacierz(L, "Macierz L (końcowa)");
        wyswietlMacierz(U, "Macierz U (końcowa)");
    }

    vector<double> permutujWektor(const vector<double>& b, const vector<int>& P) {
        int n = b.size();
        vector<double> pb(n);

        for (int i = 0; i < n; i++) {
            pb[i] = b[P[i]];
        }

        return pb;
    }

    vector<double> rozwiazLy_b(const vector<vector<double>>& L, const vector<double>& b) {
        int n = L.size();
        vector<double> y(n);

        cout << "Ly = b" << endl;
        for (int i = 0; i < n; i++) {
            double suma = 0.0;
            for (int j = 0; j < i; j++) {
                suma += L[i][j] * y[j];
            }
            y[i] = (b[i] - suma) / L[i][i];
        }

        return y;
    }

    vector<double> rozwiazUx_y(const vector<vector<double>>& U, const vector<double>& y) {
        int n = U.size();
        vector<double> x(n);

        cout << "Ux = y" << endl;
        for (int i = n - 1; i >= 0; i--) {
            double suma = 0.0;
            for (int j = i + 1; j < n; j++) {
                suma += U[i][j] * x[j];
            }
            if (fabs(U[i][i]) < 1e-10) {
                cerr << "Wartość zbyt bliska zeru!" << endl;
                exit(1);
            }
            x[i] = (y[i] - suma) / U[i][i];

            cout << n-i << ": " << x[i] << endl;
        }

        return x;
    }

    double sprawdzPoprawnosc(const vector<vector<double>>& A,
                            const vector<double>& x,
                            const vector<double>& b) {
        int n = A.size();
        double maxBled = 0.0;

        cout << "Sprawdzanie poprawności rozwiązania " << endl;
        cout << "A * x = b?" << endl;

        for (int i = 0; i < n; i++) {
            double suma = 0.0;
            for (int j = 0; j < n; j++) {
                suma += A[i][j] * x[j];
            }
            double blad = fabs(suma - b[i]);
            cout << "Wiersz " << i+1 << ": " << suma << " ?= " << b[i]
                 << " (błąd: " << blad << ")" << endl;

            if (blad > maxBled) {
                maxBled = blad;
            }
        }

        cout << "Maksymalny błąd: " << maxBled << endl;

        return maxBled;
    }

    pair<vector<double>, vector<double>> rozwiazUkladLU(const vector<vector<double>>& A,
                                                        const vector<double>& b) {
        int n = A.size();
        vector<vector<double>> L, U;
        vector<int> P;

        rozkladLU_zPivotingiem(A, L, U, P);

        vector<double> pb = permutujWektor(b, P);
        wyswietlWektor(pb, "Permutowany wektor b");

        vector<double> y = rozwiazLy_b(L, pb);

        vector<double> x = rozwiazUx_y(U, y);

        return make_pair(y, x);
    }
    void utworzMacierzDopelniona(const vector<vector<double>>& A, const vector<double>& B, vector<vector<double>>& macierzDopelniona) {
        int N = A.size();
        macierzDopelniona.resize(N, vector<double>(N + 1, 0));
        //macierz dopelniona ma o 1 wiecej kolumn, ostatnia kolumna to wyrazy wolne

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                macierzDopelniona[i][j] = A[i][j];
            }
            macierzDopelniona[i][N] = B[i];             //dodawanie wyrazow wolnych do ostatniej koluny
        }
    }
    void wyswietlMacierz(const vector<vector<double>>& macierz, const string& nazwa) {
        cout << nazwa << ":" << endl;
        for (const auto& wiersz : macierz) {
            for (double element : wiersz) {
                cout << setw(10) << fixed << setprecision(4) << element << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    void wyswietlWektor(const vector<double>& wektor, const string& nazwa) {
        cout << nazwa << ":" << endl;
        for (double element : wektor) {
            cout << setw(10) << fixed << setprecision(4) << element << " ";
        }
        cout << endl << endl;
    }
    void sprawdzWynik(const vector<vector<double>>& A, const vector<double>& B, const vector<double>& X) {
    int N = A.size();
    cout << "sprawdzanie wynikow: :" << endl;

    for (int i = 0; i < N; i++) {
        double suma = 0;
        for (int j = 0; j < N; j++) {
            suma += A[i][j] * X[j];
        }
        cout << "Równanie " << i+1 << ": ";
        cout << suma << " w przyblizeniu rowne: " << B[i] << endl;
    }
    cout << endl;
}
void sprowadzanieDoPostaciSchodkowej(vector<vector<double>>& macierzDopelniona) {
    int N = macierzDopelniona.size();

    for (int i = 0; i < N; i++) {
        if (fabs(macierzDopelniona[i][i]) < 1e-10) {                // element na przekatnej jest rowny 0 (w tym przypadku skorzystalem z precyzji)
            int indeksMax = i;
            double wartoscMax = fabs(macierzDopelniona[i][i]);

            for (int j = i + 1; j < N; j++) {                           // szukanie elementu o najwiekszym module w kolumne gdzie 0 jest na przekatnej macierzy
                if (fabs(macierzDopelniona[j][i]) > wartoscMax) {
                    indeksMax = j;
                    wartoscMax = fabs(macierzDopelniona[j][i]);
                }
            }

            if (indeksMax != i) {
                cout << "Pivoting: zamiana wierszy " << i+1 << " i " << indeksMax+1 << endl;
                macierzDopelniona[i].swap(macierzDopelniona[indeksMax]);
            }
        }

        if (fabs(macierzDopelniona[i][i]) < 1e-10) { //jesli na przekatnej sa 0 to brak rozwiazan
            cout << "brak rozwiazan!" << endl;
            continue;
        }


        for (int j = i + 1; j < N; j++) {
            double mnoznnik = macierzDopelniona[j][i] / macierzDopelniona[i][i];                //zerowanie elementow ponizej przekatnej
            for (int k = i; k <= N; k++) {
                macierzDopelniona[j][k] -= mnoznnik * macierzDopelniona[i][k];
            }
        }

        cout << "Po eliminacji dla wiersza " << i+1 << ":" << endl;                                    //wyswietlanie macierzy co iteracje po wyzerowaniu kolumny ponizej przekatnej
        for (int r = 0; r < N; r++) {
            for (int c = 0; c <= N; c++) {
                cout << setw(10) << fixed << setprecision(4) << macierzDopelniona[r][c] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}
vector<double> rozwiazUklad(const vector<vector<double>>& macierzSchodkowa) {
    int N = macierzSchodkowa.size();
    vector<double> X(N, 0);

    for (int i = N - 1; i >= 0; i--) {
        double suma = 0;
        for (int j = i + 1; j < N; j++) {
            suma += macierzSchodkowa[i][j] * X[j];
        }

        if (fabs(macierzSchodkowa[i][i]) < 1e-10) {
            cout << "dzielenie przez 0 przy x" << i+1 << endl;
            X[i] = 0; // wartosci 0 dla wartosci nieokreslonych
        } else {
            X[i] = (macierzSchodkowa[i][N] - suma) / macierzSchodkowa[i][i];
        }
    }

    return X;
}

    void testUkladu(const string& nazwaPliku) {
        vector<vector<double>> A;
        vector<double> b;
        int N;

        cout << "Testowanie układu z pliku: " << nazwaPliku << endl;

        wczytajDane(nazwaPliku, A, b, N);

        cout << "Liczba niewiadomych N: " << N << endl << endl;
        wyswietlWektor(b, "Wektor wyrazów wolnych b");
        wyswietlMacierz(A, "Macierz współczynników A");

        pair<vector<double>, vector<double>> wynik = rozwiazUkladLU(A, b);
        wyswietlWektor(wynik.first, "Rozwiązanie z (y)");
        wyswietlWektor(wynik.second, "Rozwiązanie x");

        double maxBled = sprawdzPoprawnosc(A, wynik.second, b);
    }

}