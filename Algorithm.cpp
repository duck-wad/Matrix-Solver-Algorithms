#include "Algorithm.h"
#include "Matrix.h"

std::vector<double> GaussianElimination(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {

	if (A.size() != b.size() || A[0].size() != b.size()) {
		throw std::invalid_argument("Matrix and vector must have same dimensions");
	}

	std::vector<std::vector<double>> A_REF = A;
	std::vector<double> b_REF = b;
	
	ForwardElimination(A_REF, b_REF);

	std::vector<double> x(b.size());

	//perform backward substitution to solve for x
	for (size_t i = A_REF.size(); i-- > 0; ) {
		if (i == A_REF.size() - 1) {
			x[i] = b_REF[i] / A_REF[i][i];
		}
		else {
			double sum = 0.0;
			for (size_t j = i + 1; j < A_REF.size(); j++) {
				sum += A_REF[i][j] * x[j];
			}
			x[i] = (b_REF[i] - sum) / A_REF[i][i];
		}
	}

	return x;
}

//if L is not passed in, dummy function is called which then calls actual function with an empty L
void ForwardElimination(std::vector<std::vector<double>>& A, std::vector<double>& b) {
	std::vector<std::vector<double>> L = {};
	ForwardElimination(A, b, L);
}

//if L needs to be filled, the function expects L to be sized before being passed in
void ForwardElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& L) {

	if (A.size() != b.size() || A[0].size() != b.size()) {
		throw std::invalid_argument("Matrix and vector must have same dimensions");
	}

	bool fillL = !L.empty();
	if (fillL) {
		if (L.size() != A.size() || L[0].size() != A[0].size()) {
			throw std::invalid_argument("L is not sized correctly");
		}
		for (size_t i = 0; i < L.size(); i++) {
			L[i][i] = 1.0;
		}
	}

	//append b to last column of A matrix to put into augmented form to streamline row operations
	std::vector<std::vector<double>> temp = Transpose(A);
	temp.push_back(b);
	std::vector<std::vector<double>> aug_A = Transpose(temp);
	
	//perform row operations to bring aug_A to REF
	for (size_t piv = 0; piv < aug_A.size(); piv++) {
		for (size_t i = piv + 1; i < aug_A.size(); i++) {
			double f = aug_A[i][piv] / aug_A[piv][piv];
			std::vector<double> pivRow = aug_A[piv] * f;
			aug_A[i] = aug_A[i] - pivRow;
			if (fillL) {
				L[piv][i] = f;
			}
		}
	}
	//set b to last column of augmented matrix
	//un-augment to get the REF for A
	temp = Transpose(aug_A);
	b = temp[A.size()];
	temp.pop_back();
	A = Transpose(temp);
}
