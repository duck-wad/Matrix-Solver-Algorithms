#include <iostream>
#include <vector>

#include "Algorithm.h"
#include "MatrixOperations.h"

int main() {
	std::vector<std::vector<double>> A = {
		{25, 15, -5},
		{15, 18, 0},
		{-5, 0, 11}
	};

	std::vector<double> b = { 5.0, -10.0, 14.0};

	std::cout << "A matrix: " << std::endl;
	PrintMatrix(A);
	std::cout << "b vector: " << std::endl;
	PrintVector(b);

	GaussianElimination(A, b);
	LUDecomposition(A, b);
	CholeskyDecomposition(A, b);
	GaussSeidel(A, b, 0.0001, 100);

	return 0;
}