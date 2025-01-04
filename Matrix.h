#pragma once

#include <iostream>
#include <vector>
#include <cassert>

template<typename T>
std::vector<T> operator* (const std::vector<T> row, const T c) {
	std::vector<T> output(row.size());

	for (size_t i = 0; i < row.size(); i++) {
		output[i] = row[i] * c;
	}
	return output;
}

//subtract row2 from row1
template<typename T>
std::vector<T> operator- (const std::vector<T> row1, const std::vector<T> row2) {

	if (row1.size() != row2.size()) {
		throw std::invalid_argument("Rows must be same size to subtract");
	}

	std::vector<T> output(row1.size());

	for (size_t i = 0; i < row1.size(); i++) {
		output[i] = row1[i] - row2[i];
	}
	return output;
}

//tranpose an nxm matrix
template<typename T>
std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>>& matrix) {

	size_t rows = matrix.size();
	assert(rows != 0 && "Matrix must be populated to be transposed");
	size_t cols = matrix[0].size();

	std::vector<std::vector<T>> output(cols, std::vector<T>(rows));

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			output[j][i] = matrix[i][j];
		}
	}
	return output;
}

