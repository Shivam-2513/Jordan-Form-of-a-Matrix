#include <bits/stdc++.h>
#include "matrix_rank.h"
#define loop(i, start, end) for(int i = start; i <= end; ++i)
#define rloop(i,end, start) for(int i = end; i >= start ; --i)
using namespace std;
using db = double;
template<typename T>
using V = vector<T>;
const db EPS = 1e-9;
void subtract_rows(V<V<db>>& A, int i, int j, db p, db q){
    int n = A.size(), m = A[0].size();
    loop(col, 0, m - 1){
        A[i][col] = p * A[i][col] - q * A[j][col];
    }
}
void RREF(V<V<db>>& A){
    int n = A.size(), m = A[0].size();
    int curr_row = 0, pivot_col = 0;
    while(curr_row < n && pivot_col < m){
        int pivot_row = curr_row;
        loop(i, curr_row, n - 1){
            if(abs(A[i][pivot_col]) > abs(A[curr_row][pivot_col])){
                pivot_row = i;
                swap(A[curr_row], A[pivot_row]);
            }
        }
        pivot_row = curr_row;
        if(abs(A[pivot_row][pivot_col]) > EPS /*means non zero*/){
            db x = A[pivot_row][pivot_col];
            loop(j, 0, m - 1) A[pivot_row][j] /= x;
            loop(r, 0, n - 1){
                if(r == pivot_row) continue;
                db m = A[r][pivot_col] / A[pivot_row][pivot_col];
                subtract_rows(A, r, pivot_row, 1, m);
            }
            curr_row++;
        }
        pivot_col++;
    }
}
int matrix_rank(V<V<db>> B){
    int n = B.size(), m = B[0].size();
    V<V<db>> A(n, V<db>(m, 0.0));
    loop(i, 0, n - 1) loop(j, 0, m - 1) A[i][j] = B[i][j];
    RREF(A);
    int rk = 0;
    loop(i, 0, n - 1){
        loop(j, 0, m - 1){
            if(abs(A[i][j]) > EPS){
                rk++;
                break;
            }
        }
    }
    return rk;
}