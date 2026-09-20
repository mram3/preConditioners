#include "LinearSolvers.h"

#include <iostream>

using namespace std;

int LinearSolvers::richardsonIteration
(
    const Matrix& A,
    const Preconditioners& M,
    const vector<double>& b,
    vector<double>& x,
    double tolrance,
    Side side
)
{
    int max_iter = 10000;
    int iter = 0;
    double residual = 1.0;
    vector<double> z(b.size(), 0.0);//initializing the M^-1 * (b - Ax) vector 

    while(iter < max_iter){

        //A*x
        vector<double> Ax = A.SpMV(x);
        //b-A*x
        vector<double> r = MathTools::vectorSub(b, Ax);

        residual = MathTools::L2Norm(r);

        if(residual < tolrance){
            return iter;
        }
        
        switch (side)
        {
        case Side::Left:
            M.apply(z, r);
            break;

        case Side::Right:

            break;

        case Side::Split:
        
            break;

        case Side::None:
        
            break;
        }

        //update x
        x = MathTools::vectorAdd(z, x);

        iter++;
    }

    cout << "Does not converge\n";
    return iter;
}