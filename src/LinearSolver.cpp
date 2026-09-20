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
    vector<double> z_lower(b.size(), 0.0);

    while(iter < max_iter){

        //A*x
        vector<double> Ax = A.SpMV(x);
        //b-A*x
        vector<double> r = MathTools::vectorSub(b, Ax);

        if(residual < tolrance){
            return iter;
        }
        
        switch (side)
        {
        case Side::Left:
            M.apply(z, r);
            residual = MathTools::L2Norm(r);
            x = MathTools::vectorAdd(z, x);
            break;

        case Side::Right:
            residual = MathTools::L2Norm(r);
            M.apply(z, r);
            x = MathTools::vectorAdd(z, x);
            break;

        case Side::Split:
            
            // 1. Forward sweep only: L * z_lower = r
            M.applyLower(z_lower, r); 
            
            residual = MathTools::L2Norm(z_lower);
            
            // 2. Backward sweep only: U * z = z_lower
            M.applyUpper(z, z_lower); 
            
            x = MathTools::vectorAdd(z, x);
            break;

        case Side::None:
            residual = MathTools::L2Norm(r);
            x = MathTools::vectorAdd(r, x);
            break;
        }

        if(residual < tolrance){
            return iter;
        }
        
        iter++;
    }

    cout << "Does not converge\n";
    return iter;
}