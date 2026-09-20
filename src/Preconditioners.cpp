#include "Preconditioners.h"

#include <vector>
#include <cmath>

using namespace std;

void identityPreconditioner::setup
(
    const Matrix& A
)
{

}

void identityPreconditioner::apply
(
    vector<double>& z,
    const vector<double>& r
) const
{
    int n = r.size();
    for(int i = 0; i<n; ++i){
        z[i] = r[i];
    }
}

void jacobiPreconditioner::setup
(
    const Matrix& A
)
{
    const auto& rowPtr = A.getrowPtr();
    const auto& col = A.getcol();
    const auto& values = A.getvalues();

    int n = rowPtr.size() - 1;

    invDiag.assign(n, 0.0);

    for(int i = 0; i < n; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i == col[j]){
                invDiag[i] = 1/values[j];
            }
        }
    }
}

void jacobiPreconditioner::apply
(
    vector<double>& z,
    const vector<double>& r
) const
{
    int n = r.size();

    for(int i = 0; i < n; ++i){
        z[i] = invDiag[i] * r[i];
    }
}

void sorPreconditioner::setup
(
    const Matrix& A
)
{
    matPtr = &A; //storing the address of Matrix object A inside matPtr
    const auto& rowPtr = A.getrowPtr();
    const auto& col = A.getcol();
    const auto& values = A.getvalues();

    int n = rowPtr.size() - 1;

    invDiag.assign(n, 0.0);

    for(int i = 0; i < n; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i == col[j]){
                invDiag[i] = 1/values[j];
            }
        }
    }
}

void sorPreconditioner::apply
(
    vector<double>& z,
    const vector<double>& r
) const
{
    const auto& rowPtr = matPtr->getrowPtr();//calling the Matrix object's function
    const auto& col = matPtr->getcol();
    const auto& values = matPtr->getvalues();

    int n = r.size();

    //forward substitution
    for(int i = 0; i < n; ++i){
        double sum = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(col[j] < i){
                sum += z[col[j]] * values[j];
            }
        }
        z[i] = w * (r[i] - sum) * invDiag[i];
    }
}

void sgsPreconditioner::setup
(
    const Matrix& A
)
{
    matPtr = &A;

    const auto& rowPtr = A.getrowPtr();
    const auto& col = A.getcol();
    const auto& values = A.getvalues();

    int n = rowPtr.size() - 1;

    invDiag.assign(n, 0.0);

    for(int i = 0; i < n; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i == col[j]){
                invDiag[i] = 1/values[j];
            }
        }
    }
}

void sgsPreconditioner::apply
(
    vector<double>& z,
    const vector<double>& r
) const
{
    const auto& rowPtr = matPtr->getrowPtr();
    const auto& col = matPtr->getcol();
    const auto& values = matPtr->getvalues();

    int n = r.size();

    /* Symmetric Gauss Seidel implementation
    (D+L)*D^-1*(D+U) z = r
    Let D^-1 * (D+U) z = y -> (D+L) * y = r -> forward substitution
    then (D+U) z = D * y -> backward elimination to find z
    */

    //forward substitution 
    for(int i = 0; i < n; ++i){
        double sum = 0.0;
        for(int j = rowPtr[i]; j<rowPtr[i+1]; ++j){
            if(i > col[j]){
                sum += values[j] * z[col[j]];
            }
        }
        z[i] = (r[i] - sum) * invDiag[i];
    }

    //y = D*y
    for(int i = 0; i < n; i++){
        z[i] = z[i]/invDiag[i];
    }

    //backward elimination
    for(int i = n-1; i >=0; --i){
        double sum = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i < col[j]){
                sum += values[j] * z[col[j]];
            }
        }
        z[i] = (z[i]-sum) * invDiag[i];
    }
}

void ilu0Preconditioner::setup
(
    const Matrix& A
)
{
    matPtr = &A;

    const auto& rowPtr = A.getrowPtr();
    const auto& col = A.getcol();
    const auto& values = A.getvalues();

    luvalues = values;

    int n = rowPtr.size() - 1;
    vector<int> diagIdx(n, 0);

    for(int i = 0; i < n; ++i){//to store where the diagonals are
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i == col[j]){
                diagIdx[i] = j;
                break;
            }
        }
    }

    //ILU(0) factorization
    for(int i = 1; i < n; ++i){//we can skip first row because u[0][j] = a[0][j]
        for(int j1 = rowPtr[i]; j1 < rowPtr[i+1]; ++j1){
            int k = col[j1];//current column

            if(i<=k) break;//we're looking for lower triangular elements within the row
                           //if all lower elements are zero no decomposition needs to be done there

            luvalues[j1] /= luvalues[diagIdx[k]];//dividing the lower triangular element
                                                 //with the diagonal element of its own column 
                                                 //Lij = Lij/Uii
            
            for(int j2 = diagIdx[k]+1; j2<rowPtr[k+1]; ++j2){
                //now going to the row of dividing diagonal and searching through it's columns
                int target = col[j2];

                for(int l = j1 +1; l < rowPtr[i+1]; ++l){
                    /*coming back to the original row & checking if there's non zero value in the 
                    column corresponding to non zero value in row of dividing diagonal*/
                    if(col[l] == target){
                        luvalues[l] -= luvalues[j1]*luvalues[j2];
                        //if a match is found, we do Aij = Aij - Lik*Ukj
                    }

                    //if a match is not found we do nothing, eliminating the fill in
                    //and preserving the original CSR (sparsity) structure
                }
            }
        }
    }
}

void ilu0Preconditioner::apply
(
    vector<double>& z,
    const vector<double>& r
) const
{
    const auto& rowPtr = matPtr->getrowPtr();
    const auto& col = matPtr->getcol();

    int n = r.size();
    /*
    ILU(0) Implementation
    to solve: (LU)z = r
    Let Uz = y -> solve Ly = r (forward elimination)
    Then solve U*z = y (Backward substitution)
    */
    for(int i = 0; i < n; ++i){//forward elimination
        double sum = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i > col[j]){
                sum += luvalues[j]*z[col[j]];
            }
        }
        z[i] = r[i] - sum;
    }

    for(int i = n-1; i >=0; --i){//backward substitution 
        double sum = 0.0;
        double diag = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i < col[j]){
                sum += luvalues[j] * z[col[j]];
            }
            else if(i == col[j]){
                diag = luvalues[j];
            }
        }
        z[i] = (z[i]-sum) / diag;
    }
}

void ic0Preconditioner::setup
(
    const Matrix& A
)
{
    matPtr = &A;

    const auto& rowPtr = A.getrowPtr();
    const auto& col = A.getcol();
    const auto& values = A.getvalues();

    lltvalues = values;

    int n = rowPtr.size() - 1;
    vector<int> diagIdx(n, 0);

    for(int i = 0; i < n; ++i){//storing the diagonals' indices
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i == col[j]){
                diagIdx[i] = j;
                break;
            }
        }
    }

    //IC(0) factorization
    for(int i = 0; i < n; ++i){
        for(int j1 = rowPtr[i]; j1 < rowPtr[i+1]; ++j1){
            int k = col[j1];//current column

            if(k>i) break;//looking at only lower triangle

            else if (k<i){
                double sum = 0.0; //to store the dot prod between row i and row k

                //searching within the i row if previous column has any entires
                for(int j2 = rowPtr[i]; j2 < j1; ++j2){
                    int l = col[j2];
                    //going to kth row to check if there's a nonzero match
                    for(int j3 = rowPtr[k]; j3 < rowPtr[k+1]; ++j3){
                        if(l == col[j3]){
                            sum += lltvalues[j2]*lltvalues[j3];
                            //if there's a matching nonzero value in same column in both rows
                            //we do dot product Lil * Lkl
                            break;
                        }
                    }
                }
                //then we apply Lik = (Aik - (Lil*Lkl)) / Aii
                lltvalues[j1] -= sum;
                lltvalues[j1] /= lltvalues[diagIdx[k]];

                //Mirroring L to U(L^T)
                for(int j4 = rowPtr[k]; j4 < rowPtr[k+1]; ++j4){
                    if(col[j4] == i){
                        lltvalues[j4] = lltvalues[j1];
                        break;
                    }
                }
            }

            else if(k==i){
                double sum = 0.0;
                for(int j2 = rowPtr[i]; j2 < j1; ++j2){
                    sum += lltvalues[j2] * lltvalues[j2];
                }

                lltvalues[j1] = sqrt(lltvalues[j1]-sum);
            }
        }
    }
}

void ic0Preconditioner::apply
(
    vector<double>& z,
    const vector<double>& r
) const
{
    const auto& rowPtr = matPtr->getrowPtr();
    const auto& col = matPtr->getcol();

    int n = r.size();
    /*
    IC(0) Implementation
    to solve: (LL^T)z = r
    Let L^Tz = y -> solve Ly = r (forward elimination)
    Then solve L^T*z = y (Backward substitution)
    */
    for(int i = 0; i < n; ++i){//forward elimination
        double sum = 0.0;
        double diag = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i > col[j]){
                sum += lltvalues[j]*z[col[j]];
            }
            else if(i == col[j]){
                diag = lltvalues[j];
            }
        }
        z[i] = (r[i] - sum)/ diag;
    }

    for(int i = n-1; i >=0; --i){//backward substitution 
        double sum = 0.0;
        double diag = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i < col[j]){
                sum += lltvalues[j] * z[col[j]];
            }
            else if(i == col[j]){
                diag = lltvalues[j];
            }
        }
        z[i] = (z[i]-sum) / diag;
    }
}