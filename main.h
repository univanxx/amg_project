#ifndef MAIN_H
#define MAIN_H

#include <cmath>

const double double_precision = 1e-9;	//точность сравнения double'ов

enum ElementType {
    Internal = 0,
    FreedomExit = 10,
    HardBorder = 20,
    Symmetry = 30,
    Custom = 40
};

#endif

