#pragma once

#include <vector>

using Matrix = std::vector<std::vector<double>>;

void subtract_rows(
    Matrix& A, int r1, int r2, double factor, double tol
);

void RREF(Matrix& A);

int matrix_rank(Matrix A);