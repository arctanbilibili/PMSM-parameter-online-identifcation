#include "matrix.h"

void* matrix_alloc(int size) {
    return 0;
}
// 矩阵加法
void matrix_add(const float* A, const float* B, float* C, int Ar, int Ac) {
    int i, j;
    for (i = 0; i < Ar; i++) {
        for (j = 0; j < Ac; j++) {
            int index = Ac * i + j;
            *(C + index) = *(A + index) + *(B + index);
        }
    }
}
// 矩阵减法
void matrix_sub(const float* A, const float* B, float* C, int Ar, int Ac) {
    int i, j;
    for (i = 0; i < Ar; i++) {
        for (j = 0; j < Ac; j++) {
            int index = Ac * i + j;
            *(C + index) = *(A + index) - *(B + index);
        }
    }
}
// 矩阵乘法
void matrix_mul(const float* A,
                const float* B,
                float* C,
                int Ar,
                int AcBr,
                int Bc) {
    int i, j, k;
    for (i = 0; i < Ar; i++) {
        for (j = 0; j < Bc; j++) {
            int Cindex = Bc * i + j;
            *(C + Cindex) = 0.0;
            for (k = 0; k < AcBr; k++) {
                //*(C+Cindex) += A[i][k] * B[k][j];
                *(C + Cindex) += *(A + AcBr * i + k) * *(B + Bc * k + j);
            }
        }
    }
}
// 矩阵乘法
void matrix_const_mul(const float* A, float k, float* C, int Ar, int Ac) {
    int i, j;
    for (i = 0; i < Ar; ++i) {
        for (j = 0; j < Ac; ++j) {
            *(C + Ac * i + j) = *(A + Ac * i + j) * k;
        }
    }
}
// 简易-1^k
static int pow_1(int k) {
    return k % 2 ? 1 : -1;
}
// 矩阵行列式
float matrix_det(const float* A, int n) {
    float result = 0.0;
    if (n == 1) {
        result = A[0];
    } else if (n == 2) {
        result = A[0] * A[n + 1] - A[1] * A[n];
    } else {
        int sign = 1;
        int i, j, k;
        float submatrix[M_N * M_N];
        for (i = 0; i < n; i++) {
            int sub_i = 0, sub_j = 0;
            for (j = 1; j < n; j++) {
                for (k = 0; k < n; k++) {
                    if (k != i) {
                        submatrix[sub_j * (n - 1) + sub_i] = A[j * n + k];
                        sub_i++;
                        if (sub_i == n - 1) {
                            sub_i = 0;
                            sub_j++;
                        }
                    }
                }
            }
            result += sign * A[i] * matrix_det(submatrix, n - 1);
            sign = -sign;
        }
    }
    return result;
}
// 矩阵逆是否有效
static int matrix_inv_valid(const float* B, int Bn) {
    float det = matrix_det(B, Bn);
    if (det > -M_ep && det < M_ep)
        return -1;
    return 0;
}
// A->B
void matrix_copy(const float* A, float* B, int Ar, int Ac) {
    int i, j;
    for (i = 0; i < Ar; ++i) {
        for (j = 0; j < Ac; ++j) {
            *(B + Ac * i + j) = *(A + Ac * i + j);
        }
    }
}
// 矩阵转置
void matrix_tanspose(const float* A, float* B, int Ar, int Ac) {
    int i, j;
    for (i = 0; i < Ar; ++i) {
        for (j = 0; j < Ac; ++j) {
            *(B + Ar * j + i) = *(A + Ac * i + j);
        }
    }
}
// 矩阵求逆
void matrix_inv(const float* A, float* B, int n) {
    // 使用高斯-约旦消元法求解矩阵的逆矩阵
    int i, j, k;
    float temp;
    float Atmp[M_N * M_N];
    matrix_copy(A, Atmp, n, n);
    // 将B初始化为单位矩阵
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            *(B + i * n + j) = (i == j) ? 1.0 : 0.0;
        }
    }
    // 高斯-约旦消元
    for (i = 0; i < n; i++) {
        temp = *(Atmp + i * n + i);
        for (j = 0; j < n; j++) {
            *(Atmp + i * n + j) /= temp;
            *(B + i * n + j) /= temp;
        }
        for (j = 0; j < n; j++) {
            if (i != j) {
                temp = *(Atmp + j * n + i);
                for (k = 0; k < n; k++) {
                    *(Atmp + j * n + k) -= *(Atmp + i * n + k) * temp;
                    *(B + j * n + k) -= *(B + i * n + k) * temp;
                }
            }
        }
    }
}
// A=单位矩阵
void matrix_I(float* A, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j) {
                *(A + i * n + j) = 1.0f;
            } else {
                *(A + i * n + j) = 0.0f;
            }
        }
    }
}
// 0矩阵
void matrix_Zero(float* A, int Ar, int Ac) {
    int i, j;
    for (i = 0; i < Ar; ++i) {
        for (j = 0; j < Ac; ++j) {
            *(A + i * Ac + j) = 0.0f;
        }
    }
}
// 矩阵除法
// 0 正常 -1 无效
int matrix_div(const float* A, const float* B, float* C, int Ar, int AcBn) {
    if (matrix_inv_valid(B, AcBn) != 0) {
        return -1;
    }
    // 将矩阵除法转化为矩阵乘法，即C = A * B^-1
    float inv_B[M_N * M_N];
    matrix_inv(B, inv_B, AcBn);
    matrix_mul(A, inv_B, C, Ar, AcBn, AcBn);
    return 0;
}
