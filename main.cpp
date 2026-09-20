/*
Compiling Instruction: g++ main.cpp src/*.cpp -Iinclude -std=c++17 -o main && ./main
*/
#include <iostream>

#include "LinearSolvers.h"

int main(){
    int Nx = 50; //number of cells along x axis
    int Ny = 50; //number of cells along y axis

    int NCells = Nx*Ny;
    std::vector<double> x(NCells, 1.0); //initialize the solution vector

    Equations eqn; 
    LinearSolvers solve;
    ic0Preconditioner iP; //preconditioner

    eqn.assemblePoissonMatrix(Nx, Ny); //creates a square matrix of NCells*NCells
    
    iP.setup(eqn.A);//sets up the preconditioner 
    
    //apply preconditioner and solve
    int iter = solve.richardsonIteration(eqn.A, iP, eqn.b, x, 1e-6, Side::Left);

    std::cout << iter<< std::endl;
    return 0;
}