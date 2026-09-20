/*
Compiling Instruction: g++ main.cpp src/*.cpp -Iinclude -std=c++17 -o main && ./main
*/
#include <iostream>

#include "LinearSolvers.h"

using namespace std;

int main(){
    int Nx = 50; //number of cells along x axis
    int Ny = 50; //number of cells along y axis

    int NCells = Nx*Ny;
    std::vector<double> x(NCells, 0.0); //initialize the solution vector

    Equations eqn; 
    LinearSolvers solve;
 
    eqn.assemblePoissonMatrix(Nx, Ny); //creates a square matrix of NCells*NCells
    
    jacobiPreconditioner jP; //preconditioner
    sorPreconditioner gS(1.0);
    sorPreconditioner sP(1.88);
    sgsPreconditioner sGS;
    ilu0Preconditioner ilu;
    ic0Preconditioner ic;

    jP.setup(eqn.A);//sets up the preconditioner
    gS.setup(eqn.A);
    sP.setup(eqn.A);
    sGS.setup(eqn.A);
    ilu.setup(eqn.A);
    ic.setup(eqn.A);
    
    //apply preconditioner and solve
    int iterJP = solve.richardsonIteration(eqn.A, jP, eqn.b, x, 1e-6, Side::Right);

    std::fill(x.begin(), x.end(), 0.0);
    int iterGS = solve.richardsonIteration(eqn.A, gS, eqn.b, x, 1e-6, Side::Right);

    std::fill(x.begin(), x.end(), 0.0);
    int iterSP = solve.richardsonIteration(eqn.A, sP, eqn.b, x, 1e-6, Side::Right);

    std::fill(x.begin(), x.end(), 0.0);
    int iterSGS = solve.richardsonIteration(eqn.A, sGS, eqn.b, x, 1e-6, Side::Split);

    std::fill(x.begin(), x.end(), 0.0);
    int iterILU = solve.richardsonIteration(eqn.A, ilu, eqn.b, x, 1e-6, Side::Split);

    std::fill(x.begin(), x.end(), 0.0);
    int iterIC = solve.richardsonIteration(eqn.A, ic, eqn.b, x, 1e-6, Side::Split);

    cout << "Jacopi Preconditioner      : "<<iterJP<< std::endl;
    cout << "Gauss Seidel Preconditioner: "<<iterGS<< std::endl;
    cout << "Successive Over Relaxation : "<<iterSP<< std::endl;
    cout << "Symmetic Gauss Seidel      : "<<iterSGS<< std::endl;
    cout << "ILU(0) Preconditioner      : "<<iterILU<< std::endl;
    cout << "IC(0) Preconditioner       : "<<iterIC<< std::endl;
    return 0;
}