#include "Equations.h"

using namespace std;

void Equations::assemblePoissonMatrix
(
    int Nx,
    int Ny
) {
    //creating a test matrix to check the preconditioner performance
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
            int cell_P = i * Ny + j; 
            
            double a_P = 0.0; //diagonal co efficient
            double Su = 1.0;  //source term (taken as 1.0)

            // Boundary conditions
            double phi_West  = 0.0;
            double phi_South = 0.0;
            double phi_East  = 1.0;
            double phi_North = 1.0;

            // 1. WEST FACE (i - 1)
            if (i > 0) { // Internal face
                int cell_W = (i - 1) * Ny + j;
                A.addCoeff(cell_P, cell_W, -1.0);
                a_P += 1.0;
            } else {                        // Boundary face
                a_P += 2.0;                 // cell center to boundary is half the cell thickness
                Su += 2.0 * phi_West;       // moving the boundary term to RHS
            }

            // 2. EAST FACE (i + 1)
            if (i < Nx - 1) { 
                int cell_E = (i + 1) * Ny + j;
                A.addCoeff(cell_P, cell_E, -1.0);
                a_P += 1.0;
            } else {     
                a_P += 2.0; 
                Su += 2.0 * phi_East;
            }

            // 3. SOUTH FACE (j - 1)
            if (j > 0) { 
                int cell_S = i * Ny + (j - 1);
                A.addCoeff(cell_P, cell_S, -1.0);
                a_P += 1.0;
            } else {     
                a_P += 2.0;
                Su += 2.0 * phi_South;
            }

            // 4. NORTH FACE (j + 1)
            if (j < Ny - 1) { 
                int cell_N = i * Ny + (j + 1);
                A.addCoeff(cell_P, cell_N, -1.0);
                a_P += 1.0;
            } else {     
                a_P += 2.0;
                Su += 2.0 * phi_North;
            }

            // 5. ADD DIAGONAL AND RHS
            A.addCoeff(cell_P, cell_P, a_P);
            b.push_back(Su);
        }
    }

    //Comprssing COO to CSR
    A.compressToCSR();
}