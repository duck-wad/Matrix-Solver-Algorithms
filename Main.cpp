#include <iostream>
#include <vector>

#include "Algorithm.h"
#include "Matrix.h"

void PrintVector(const std::vector<double>& vec) {
	for (size_t i = 0; i < vec.size(); i++) {
		std::cout << vec[i] << std::endl;
	}
}

int main() {
	std::vector<std::vector<double>> A = {
		{1.0, 1.0, 0.0},
		{-2.0, -1.0, 2.0},
		{3.0, 6.0, 7.0}
	};
	std::vector<double> b = { 5.0, -10.0, 14.0 };

	std::vector<double> x = LUDecomposition(A, b);

	PrintVector(x);

	return 0;
}