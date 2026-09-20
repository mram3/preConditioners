#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix{
public:

    void addCoeff(
        int row, 
        int colm,
        double val
    );

    void compressToCSR();

    std::vector<double> SpMV(
        const std::vector<double>& result
    )const;

    const std::vector<int>& getrowPtr() const;
    const std::vector<int>& getcol() const;
    const std::vector<double>& getvalues() const;

private:

    //COO array
    std::vector<int> rowsCOO;

    //CSR arrays
    std::vector<int> rowPtr;
    std::vector<int> col;
    std::vector<double> values;
};

#endif