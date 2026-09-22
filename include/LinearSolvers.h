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

    void directCholesky
    (
        const Matrix& A_sparse, 
        const std::vector<double>& b, 
        std::vector<double>& x
    );

    void directLU
    (
        const Matrix& A_sparse, 
        const std::vector<double>& b, 
        std::vector<double>& x
    );

};
#endif
