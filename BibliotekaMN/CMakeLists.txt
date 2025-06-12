cmake_minimum_required(VERSION 3.31)
project(BibliotekaMN)

set(CMAKE_CXX_STANDARD 20)

add_executable(BibliotekaMN main.cpp
        include/biblioteka_metody_numeryczne.h
        include/interpolacja.h
        src/interpolacja.cpp
        include/uklady_liniowe.h
        src/uklady_liniowe.cpp
        przyklady/przyklad_uklady_liniowe.cpp
        include/calkowanie_numeryczne.h
        src/calkowanie_numeryczne.cpp
        include/aproksymacja.h
        src/aproksymacja.cpp
        include/rownania_rozniczkowe.h
        src/rownania_rozniczkowe.cpp
        include/rownania_nieliniowe.h
        src/rownania_nieliniowe.cpp
        testy/testy.cpp
)

