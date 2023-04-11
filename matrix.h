#ifndef _MATRIX_H_
#define _MATRIX_H_

#define M_N 4
#define M_ep 1e-6

void* matrix_alloc(int size);
// 矩阵加法
void matrix_add(const float* A, const float* B, float* C, int Ar, int Ac);
// 矩阵减法
void matrix_sub(const float* A, const float* B, float* C, int Ar, int Ac);
// 矩阵乘法
void matrix_mul(const float* A,
                const float* B,
                float* C,
                int Ar,
                int AcBr,
                int Bc);
// 矩阵乘法
void matrix_const_mul(const float* A, float k, float* C, int Ar, int Ac);
// 矩阵行列式
float matrix_det(const float* A, int n);
// A->B
void matrix_copy(const float* A, float* B, int Ar, int Ac);
// 矩阵转置
void matrix_tanspose(const float* A, float* B, int Ar, int Ac);
// 矩阵求逆
void matrix_inv(const float* A, float* B, int n);
// A=单位矩阵
void matrix_I(float* A, int n);
// 0矩阵
void matrix_Zero(float* A, int Ar, int Ac);
// 矩阵除法
// 0 正常 -1 无效
int matrix_div(const float* A, const float* B, float* C, int Ar, int AcBn);

#endif
