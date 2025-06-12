#include "include/biblioteka_metody_numeryczne.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <unordered_map>
using namespace std;
using namespace std::chrono;
using namespace biblioteka_numeryczna;

const double alfa = 3.0e-12;
const double T0 = 2973;
const double czas_chlodzenia = 2973.0;

int main() {
   cout << "========================================" << endl;
   cout << "URUCHAMIANIE TESTÓW JEDNOSTKOWYCH" << endl;
   cout << "Biblioteka Metod Numerycznych" << endl;
   cout << "Autor: Szymon Ros" << endl;
   cout << "Data: 11/06/2025" << endl;
   cout << "========================================" << endl;

   try {
       // Testy układów liniowych
       testEliminacjaGaussa1();
       testEliminacjaGaussa2();
       testMetodaGaussaSeidela1();
       testMetodaGaussaSeidela2();
       testRozkladLU1();
       testRozkladLU2();

       // Testy interpolacji
       testInterpolacjaLagrangea1();
       testInterpolacjaLagrangea2();
       testWspolczynnikiNewton1();
       testWspolczynnikiNewton2();
       testInterpolacjaNewtona1();
       testInterpolacjaNewtona2();
       testMetodaHornera1();
       testMetodaHornera2();

       // Testy aproksymacji
       testIloczynSkalarnyElementowZBazy1();
       testIloczynSkalarnyElementowZBazy2();
       testFunkcjaAproksymowana1();
       testFunkcjaAproksymowana2();
       testIloczynSkalarnyElementuZBazyIFunkcji1();
       testIloczynSkalarnyElementuZBazyIFunkcji2();
       testBladAproksymacji1();
       testBladAproksymacji2();

       // Testy całkowania numerycznego
       testFunkcjaTrygonometryczna1();
       testFunkcjaTrygonometryczna2();
       testHorner1();
       testHorner2();
       testMetodaProstokatow1();
       testMetodaProstokatow2();
       testMetodaTrapezow1();
       testMetodaTrapezow2();
       testSimpson1();
       testSimpson2();
       testCalkaAnalitycznaWielomianu1();
       testCalkaAnalitycznaWielomianu2();
       testWartoscFunkcji1();
       testWartoscFunkcji1_2();
       testWartoscFunkcji2_1();
       testWartoscFunkcji2_2();

       // Testy równań różniczkowych
       testObliczPochodna1();
       testObliczPochodna2();
       testRozwiazanieDokladne1();
       testRozwiazanieDokladne2();
       testMetodaEulera1();
       testMetodaEulera2();
       testObliczBladSredniKwadratowy1();
       testObliczBladSredniKwadratowy2();
       testObliczZmianeTemperatury1();
       testObliczZmianeTemperatury2();
       testObliczMetodaEulera1();
       testObliczMetodaEulera2();
       testObliczMetodaSrodka1();
       testObliczMetodaSrodka2();
       testObliczMetodaHeuna1();
       testObliczMetodaHeuna2();
       testObliczMetodaRK41();
       testObliczMetodaRK42();
       testTemperaturaDokladna1();
       testTemperaturaDokladna2();
       testObliczSredniBladKwadratowy1();
       testObliczSredniBladKwadratowy2();

       cout << "\n========================================" << endl;
       cout << "TESTY ZAKOŃCZONE" << endl;
       cout << "========================================" << endl;

   } catch (const exception& e) {
       cout << "\nBŁĄD PODCZAS TESTÓW: " << e.what() << endl;
       return 1;
   }

   return 0;
}