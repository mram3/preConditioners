#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

#include "Matrix.h"
#include <stdexcept>

class Preconditioners
{
public:

    virtual ~Preconditioners() = default;

    virtual void setup(const Matrix& A) = 0;

    virtual void apply(std::vector<double>& z, const std::vector<double>& r) const = 0;

    virtual void applyLower(std::vector<double>& z, const std::vector<double>& r) const {
        throw std::runtime_error("Error: Split preconditioning not supported for this preconditioner.");
    }
    
    virtual void applyUpper(std::vector<double>& z, const std::vector<double>& r) const {
        throw std::runtime_error("Error: Split preconditioning not supported for this preconditioner.");
    }
};

class identityPreconditioner : public Preconditioners
{
public:

    void setup(const Matrix& A) override;
    void apply(std::vector<double>& z, const std::vector<double>& r)const override;
};

class jacobiPreconditioner : public Preconditioners
{
public:
    void setup(const Matrix& A) override;
    void apply(std::vector<double>& z, const std::vector<double>& r)const override;

private:
    std::vector<double> invDiag;
};

class sorPreconditioner : public Preconditioners
{
public:
    sorPreconditioner(double w): w(w) {};//constructor
    void setup(const Matrix& A) override;
    void apply(std::vector<double>& z, const std::vector<double>& r)const override;

private:
    double w; //relaxation factor
    std::vector<double> invDiag;
    const Matrix* matPtr = nullptr;//pointer to Matrix object
};

class sgsPreconditioner : public Preconditioners
{
public:
    void setup(const Matrix& A) override;
    void apply(std::vector<double>& z, const std::vector<double>& r)const override;

private:
    std::vector<double> invDiag;
    const Matrix* matPtr = nullptr;
};

class ilu0Preconditioner : public Preconditioners
{
public:
    void setup(const Matrix& A) override;
    void apply(std::vector<double>& z, const std::vector<double>& r)const override;

private:
    std::vector<double> luvalues;
    const Matrix* matPtr = nullptr;
};

class ic0Preconditioner : public Preconditioners
{
public:
    void setup(const Matrix& A) override;
    void apply(std::vector<double>& z, const std::vector<double>& r)const override;

private:
    std::vector<double> lltvalues;
    const Matrix* matPtr = nullptr;
};

#endif