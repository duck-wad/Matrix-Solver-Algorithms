#pragma once

#include <iostream>
#include <vector>
#include <cmath>

void GaussianElimination(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
void LUDecomposition(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
void CholeskyDecomposition(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
void GaussSeidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const double errorTol = 0.001, const int maxIter = 500);

//puts the input matrix A into reduced-echelon form
//if matrix L is passed in as argument, it will be filled as the bottom triangle for LU decomposition
void ForwardElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& L);
//dummy function if an L matrix is not passed in
void ForwardElimination(std::vector<std::vector<double>>& A, std::vector<double>& b);

std::vector<double> ForwardSubstitution(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
std::vector<double> BackwardSubstitution(const std::vector<std::vector<double>>& A, const std::vector<double>& b);