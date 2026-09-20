#ifndef EQUATIONS_H
#define EQUATIONS_H

#include "Matrix.h"

#include <vector>

class Equations{
public:

    Matrix A;
    std::vector<double> b;

    void assemblePoissonMatrix(int Nx, int Ny);

};

#endif