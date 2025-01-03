#pragma once

#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> GaussianElimination(const std::vector<std::vector<double>> A, const std::vector<double> b);
std::vector<double> LUDecomposition(const std::vector<std::vector<double>> A, const std::vector<double> b);
std::vector<double> CholeskyDecomposition(const std::vector<std::vector<double>> A, const std::vector<double> b);