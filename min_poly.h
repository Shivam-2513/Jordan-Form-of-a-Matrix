#pragma once

#include <vector>
#include <utility>

using Matrix = std::vector<std::vector<double>>;

std::vector<std::pair<double, double>>
minimalPolynomialFactors(
    const Matrix& A,
    double eigenMergeTol = 1e-7,
    double rankRelTol = 1e-8,
    double qrTol = 1e-12
);