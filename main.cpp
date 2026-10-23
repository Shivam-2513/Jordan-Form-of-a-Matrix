#include <bits/stdc++.h>
#include <iomanip>
#include "min_poly.h"
#include "matrix_rank.h"
#define loop(i, start, end) for(int i = start; i <= end; ++i)
#define rloop(i,end, start) for(int i = end; i >= start ; --i)
using namespace std;
using db = double;
template<typename T>
using V = vector<T>;
using Mat = vector<vector<double>>;
const db EPS = 1e-9;

Mat diag(db x, int n){
    Mat d(n, V<db>(n, 0));
    loop(i, 0, n - 1) d[i][i] = x;
    return d;
}
V<V<db>> transpose(V<V<db>>& A){
    int n = A.size(), m = A[0].size();
    V<V<db>> At(m, V<db>(n, 0.0));
    loop(i, 0, n - 1) loop(j, 0, m - 1) At[j][i] = A[i][j];
    return At;
}
V<V<db>> matmul(V<V<db>>& A, V<V<db>>& B){
    int n = A.size(), m = A[0].size(), r = B[0].size();
    V<V<db>> C(n, V<db>(r, 0.0));
    loop(i, 0, n - 1){
        loop(j, 0, r - 1){
            loop(k, 0, m - 1){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}
Mat matsub(Mat &A, Mat &B){
    int n = A.size(), m = A[0].size();
    Mat C(n, V<db>(m, 0));
    loop(i, 0, n - 1) loop(j, 0, m - 1) C[i][j] = A[i][j] - B[i][j];
    return C;
}
V<V<db>> system_solve(V<V<db>>& A, V<V<db>>& b){
    int n = A.size(), m = A[0].size();
    V<V<db>> x(m, V<db>(1, 0.0));
    auto At = transpose(A);
    auto M = matmul(At, A);
    auto c = matmul(At, b);
    loop(i, 0, m - 1){
        M[i].emplace_back(c[i][0]);
    }
    RREF(M);
    rloop(i, m - 1, 0){
        x[i][0] = M[i][m];
    }
    return x;
}

int main(){
    // Mat A = {{1, 1, -1, 1, -1},
    //         {0, 1, 0, 0, 0},
    //         {0, 0, 1, 2, -1},
    //         {0, 0, 0, 3, 1},
    //         {0, 0, 0, 0, 3}};
    // Mat A = {{2, 0, 0, 0, 0, 0},
    //          {1, 2, 0, 0, 0, 0},
    //          {-1, 0, 2, 0, 0, 0},
    //          {0, 1, 0, 2, 0, 0},
    //          {1, 1, 1, 1, 2, 0},
    //          {0, 0, 0, 0, 1, -1}};
    Mat A = {
            {1, 1, 0,  0, 3, -2},
            {0, 1, 1, -1, 1, -1},
            {0, 0, 1,  0, 0,  0},
            {0, 0, 0,  1, 3, -2},
            {0, 0, 0,  0, 4,  1},
            {0, 0, 0,  0, 0,  4}
        };
    int n = A.size();

    //--these will be used in the future if i want to generate the chain basis also--
    Mat 𝒞(n);
    Mat u(n); loop(i, 0, n - 1) u[i].emplace_back(0.0);
    //--

    V<pair<db, db>> result = minimalPolynomialFactors(A);

    int i = 0;
    Mat Jordan(n, V<db>(n, 0));
    for(auto [𝛌, AM]: result){
        V<int> d = {0}; //di means dim(Ker(A - 𝛌I)^i) = n - rank((A - 𝛌I)^i)
        /*
        AM is the block size,
        we need to find the sizes of the boxes in this block;

        no. of boxes of size >= k is dk - d(k-1)
        therefore, no. of boxes of size exactly k, nb[k] :=  (dk - d(k-1)) - (d(k+1) - dk) 
                                                          = -d(k+1) + 2dk - d(k-1)
        we find di till di == di-1 for the first time. once that happens. di == di-1 for all greater i's
        call this i as H
        for all 1 <= bs <= H-1, 
            loop(_, 1, nb[bs])
                loop(_', 1, bs-1)
                    (i, i)     -> 𝛌
                    (i, i + 1) -> 1
                    i          += 1
                (i, i) -> 𝛌
                i      += 1
        */
        // establishing the d[] vector:
        Mat 𝛌I = diag(𝛌, n);
        Mat N = matsub(A, 𝛌I);
        Mat power = N;
        while(1){
            int rk = matrix_rank(power);
            power = matmul(N, power);
            d.emplace_back(n - rk);
            if(*(d.end() - 1) == *(d.end() - 2)) break;
        }
        // establishing nb[] vector:
        V<int> nb(d.size() - 1, 0);
        loop(i, 1, d.size() - 2){
            nb[i] = -d[i + 1] + 2 * d[i] - d[i - 1];
        }
        //filling the jordan matrix for this 𝛌:
        loop(bs, 1, nb.size() - 1){
            loop(_, 1, nb[bs]){
                loop(__, 1, bs - 1){
                    Jordan[i][i] = 𝛌;
                    Jordan[i][i + 1] = 1;
                    i += 1;
                }
                Jordan[i][i] = 𝛌;
                i += 1;
            }
        }
    }

    //output the Jordan form:
    cout << "Jordan Form: \n";
    loop(i, 0, n - 1){
        loop(j, 0, n - 1) cout << Jordan[i][j] << ' ';
        cout << '\n';
    }
}