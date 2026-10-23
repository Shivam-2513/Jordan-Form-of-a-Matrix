
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include "min_poly.h"
using namespace std;
using Matrix = vector<vector<double>>;

Matrix Identity(int n) {
    Matrix I(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        I[i][i] = 1.0;
    return I;
}

double maxAbs(const Matrix& A) {
    double ans = 0.0;
    for (const auto& row : A)
        for (double x : row)
            ans = max(ans, abs(x));
    return ans;
}

// Householder QR decomposition: A = Q * R.
pair<Matrix, Matrix> householderQR(const Matrix& A) {
    int n = (int)A.size();
    Matrix Q = Identity(n), R = A;

    for (int k = 0; k < n - 1; ++k) {
        double normx = 0.0;

        for (int i = k; i < n; ++i)
            normx = hypot(normx, R[i][k]);

        if (normx <= 1e-15)
            continue;

        double alpha = -copysign(normx, R[k][k]);
        vector<double> v(n - k);

        for (int i = k; i < n; ++i)
            v[i - k] = R[i][k];

        v[0] -= alpha;

        double normv = 0.0;
        for (double x : v)
            normv = hypot(normv, x);

        if (normv <= 1e-15)
            continue;

        for (double& x : v)
            x /= normv;

        // R <- P * R
        for (int j = k; j < n; ++j) {
            double dot = 0.0;

            for (int i = k; i < n; ++i)
                dot += v[i - k] * R[i][j];

            for (int i = k; i < n; ++i)
                R[i][j] -= 2.0 * v[i - k] * dot;
        }

        // Q <- Q * P
        for (int i = 0; i < n; ++i) {
            double dot = 0.0;

            for (int j = k; j < n; ++j)
                dot += Q[i][j] * v[j - k];

            for (int j = k; j < n; ++j)
                Q[i][j] -= 2.0 * dot * v[j - k];
        }
    }

    return {Q, R};
}

// Reduce A to upper Hessenberg form using orthogonal similarities.
Matrix hessenberg(Matrix H) {
    int n = (int)H.size();

    for (int k = 0; k < n - 2; ++k) {
        int start = k + 1;
        double normx = 0.0;

        for (int i = start; i < n; ++i)
            normx = hypot(normx, H[i][k]);

        if (normx <= 1e-15)
            continue;

        double alpha = -copysign(normx, H[start][k]);
        vector<double> v(n - start);

        for (int i = start; i < n; ++i)
            v[i - start] = H[i][k];

        v[0] -= alpha;

        double normv = 0.0;
        for (double x : v)
            normv = hypot(normv, x);

        if (normv <= 1e-15)
            continue;

        for (double& x : v)
            x /= normv;

        // H <- P * H
        for (int j = k; j < n; ++j) {
            double dot = 0.0;

            for (int i = start; i < n; ++i)
                dot += v[i - start] * H[i][j];

            for (int i = start; i < n; ++i)
                H[i][j] -= 2.0 * v[i - start] * dot;
        }

        // H <- H * P
        for (int i = 0; i < n; ++i) {
            double dot = 0.0;

            for (int j = start; j < n; ++j)
                dot += H[i][j] * v[j - start];

            for (int j = start; j < n; ++j)
                H[i][j] -= 2.0 * dot * v[j - start];
        }

        for (int i = k + 2; i < n; ++i)
            H[i][k] = 0.0;
    }

    return H;
}

// Approximate all eigenvalues using shifted QR iteration.
// Assumes all eigenvalues are real.
vector<double> realEigenvaluesQR(
    Matrix A, double qrTol = 1e-12
) {
    int n = (int)A.size();
    if (n == 0) return {};

    // Check whether A is upper or lower triangular.
    double scale = maxAbs(A);
    double tol = qrTol * max(1.0, scale);

    bool upper = true, lower = true;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j && abs(A[i][j]) > tol)
                upper = false;

            if (i < j && abs(A[i][j]) > tol)
                lower = false;
        }
    }

    // For a triangular matrix, eigenvalues are its diagonal entries.
    if (upper || lower) {
        vector<double> eig(n);

        for (int i = 0; i < n; ++i)
            eig[i] = A[i][i];

        return eig;
    }
    A = hessenberg(A);

    vector<double> eig(n);
    int m = n, iterations = 0;
    const int maxIterations = 50000;

    while (m > 1) {
        double sub = abs(A[m - 1][m - 2]);

        double localScale =
            1.0 + abs(A[m - 1][m - 1])
                + abs(A[m - 2][m - 2]);

        // Deflate a converged eigenvalue.
        if (sub <= qrTol * localScale) {
            A[m - 1][m - 2] = 0.0;
            eig[m - 1] = A[m - 1][m - 1];
            --m;
            iterations = 0;
            continue;
        }

        if (++iterations > maxIterations)
            throw runtime_error(
                "QR iteration did not converge."
            );

        double a = A[m - 2][m - 2];
        double b = A[m - 2][m - 1];
        double c = A[m - 1][m - 2];
        double d = A[m - 1][m - 1];

        // Wilkinson-style shift when the trailing 2x2
        // block has real eigenvalues.
        double disc = (a - d) * (a - d) + 4.0 * b * c;
        double mu = d;

        if (disc >= 0.0) {
            double root = sqrt(disc);
            double l1 = 0.5 * (a + d + root);
            double l2 = 0.5 * (a + d - root);

            mu = (abs(l1 - d) < abs(l2 - d))
                    ? l1 : l2;
        }

        // Occasional exceptional shift.
        if (iterations % 30 == 0)
            mu = d + 0.75 * abs(c);

        Matrix B(m, vector<double>(m));

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j)
                B[i][j] = A[i][j]
                        - (i == j ? mu : 0.0);

        auto [Q, R] = householderQR(B);

        // A <- Q^T * A * Q on the active block,
        // retaining the couplings to the other rows/columns.
        Matrix temp = A;

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                double s = 0.0;

                for (int k = 0; k < m; ++k)
                    s += Q[k][i] * temp[k][j];

                A[i][j] = s;
            }
        }

        temp = A;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                double s = 0.0;

                for (int k = 0; k < m; ++k)
                    s += temp[i][k] * Q[k][j];

                A[i][j] = s;
            }
        }

        // Remove negligible entries below the subdiagonal.
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i - 1; ++j)
                if (abs(A[i][j]) < 1e-10)
                    A[i][j] = 0.0;
    }

    eig[0] = A[0][0];
    return eig;
}

// Numerical rank using Gaussian elimination and partial pivoting.
int numericalRank(Matrix A, double relTol = 1e-8) {
    int n = (int)A.size();
    if (n == 0)
        return 0;

    double scale = maxAbs(A);
    if (scale == 0.0)
        return 0;

    double tol = relTol * scale;
    int r = 0;

    for (int col = 0; col < n && r < n; ++col) {
        int p = r;

        for (int i = r + 1; i < n; ++i)
            if (abs(A[i][col]) > abs(A[p][col]))
                p = i;

        if (abs(A[p][col]) <= tol)
            continue;

        swap(A[p], A[r]);

        double pivot = A[r][col];

        for (int i = r + 1; i < n; ++i) {
            double factor = A[i][col] / pivot;
            A[i][col] = 0.0;

            for (int j = col + 1; j < n; ++j)
                A[i][j] -= factor * A[r][j];
        }

        ++r;
    }

    return r;
}

// Returns {eigenvalue, exponent in minimal polynomial}.
// Assumes a real spectrum.
//
// Tolerances can be adjusted for the input matrix.
vector<pair<double, double>> minimalPolynomialFactors(
    const Matrix& input,
    double eigenMergeTol,
    double rankRelTol,
    double qrTol
) {
    int n = (int)input.size();

    if (n == 0)
        throw invalid_argument("Matrix must be non-empty.");

    for (const auto& row : input)
        if ((int)row.size() != n)
            throw invalid_argument("Matrix must be square.");

    for (const auto& row : input)
        for (double x : row)
            if (!isfinite(x))
                throw invalid_argument(
                    "Matrix entries must be finite."
                );

    // Normalize to reduce overflow/underflow risks.
    double scale = maxAbs(input);

    if (scale == 0.0)
        return {{0.0, 1.0}};

    Matrix A = input;

    for (auto& row : A)
        for (double& x : row)
            x /= scale;

    vector<double> vals = realEigenvaluesQR(A, qrTol);
    sort(vals.begin(), vals.end());

    // Cluster approximately equal eigenvalues.
    struct Cluster {
        double sum;
        int count;
    };

    vector<Cluster> clusters;

    for (double x : vals) {
        if (clusters.empty()) {
            clusters.push_back({x, 1});
            continue;
        }

        double mean =
            clusters.back().sum / clusters.back().count;

        double tol = eigenMergeTol
                   * max({1.0, abs(mean), abs(x)});

        if (abs(x - mean) <= tol) {
            clusters.back().sum += x;
            clusters.back().count++;
        } else {
            clusters.push_back({x, 1});
        }
    }

    vector<pair<double, double>> ans;

    for (const auto& cl : clusters) {
        double lambda = cl.sum / cl.count;
        int algebraicMultiplicity = cl.count;

        // B = A - lambda*I.
        Matrix B = A;

        for (int i = 0; i < n; ++i)
            B[i][i] -= lambda;

        // The exponent is the smallest k for which
        // nullity((A - lambda*I)^k) equals the
        // algebraic multiplicity of lambda.
        Matrix power = B;
        int exponent = algebraicMultiplicity;

        for (int k = 1; k <= algebraicMultiplicity; ++k) {
            int nullity =
                n - numericalRank(power, rankRelTol);

            if (nullity >= algebraicMultiplicity) {
                exponent = k;
                break;
            }

            // power <- power * B
            Matrix next(n, vector<double>(n, 0.0));

            for (int i = 0; i < n; ++i) {
                for (int t = 0; t < n; ++t) {
                    if (power[i][t] == 0.0)
                        continue;

                    for (int j = 0; j < n; ++j)
                        next[i][j] += power[i][t] * B[t][j];
                }
            }

            power.swap(next);
        }

        ans.push_back({
            lambda * scale,
            (double)exponent
        });
    }

    return ans;
}
