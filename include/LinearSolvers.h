#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "Preconditioners.h"
#include "Equations.h"
#include "MathTools.h"
#include "Type.h"

class LinearSolvers{
public:

    static int richardsonIteration
    (
        const Matrix& A,
        const Preconditioners& M,
        const std::vector<double>& b,
        std::vector<double>& x,
        double tolerance,
        Side side
    );
};
#endif
