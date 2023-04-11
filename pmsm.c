#include "pmsm.h"
#include <stdarg.h>

static float Ld;
static float Lq;
static float R;
static float Fai;

static float xk6_1[5];  // t-1 时刻值
static float xk6[5];    // t   时刻值
static float output[4000][4];  // 备份池

static FILE* fp = NULL;
int openxk6(const char* fileName) {
    int j;
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Failed to open file.\n");
        return -1;
    }
    // 把第一行读到xk6
    for (j = 0; j < 6; ++j) {
        fscanf(fp, "%f", &xk6_1[j]);
        xk6[j] = xk6_1[j];
    }
    return 0;
}
int readxk6() {
    static int i = 0;
    int j;
    for (j = 0; j < 6; ++j) {
        xk6_1[j] = xk6[j];  // 备份
        fscanf(fp, "%f", &xk6[j]);
    }
    // 定位到第line行
    char line[200];
    fgets(line, 200, fp);
    if (feof(fp)) {
        return -1;  // 读取无效
    }
    return i++;
}
int closexk6() {
    fclose(fp);
    return 0;
}
int save2output(int arg, ...) {
    static int i = 0;
    int j;
    va_list args;
    va_start(args, arg);
    for (j = 0; j < arg; j++) {
        output[i][j] = va_arg(args, double);
    }
    va_end(args);
    ++i;
    return 0;
}
int output2fd(const char* fileName, int r, int c) {
    int i, j;
    FILE* fpout = NULL;
    fpout = fopen(fileName, "w+");
    if (fpout == NULL) {
        printf("save Fault.\n");
        return -1;
    }
    for (i = 0; i < r; ++i) {
        for (j = 0; j < c; ++j) {
            fprintf(fpout, "%.7f ", output[i][j]);
        }
        fprintf(fpout, "\n");
    }
    fclose(fpout);
    return 0;
}

float readId() {
    return xk6[0];
}
float readIq() {
    return xk6[1];
}
float readId_1() {
    return xk6_1[0];
}
float readIq_1() {
    return xk6_1[1];
}
float readUd() {
    return xk6_1[2];
}
float readUq() {
    return xk6_1[3];
}
float readWe() {
    return xk6[4];
}
// LdLq迭代初始值
float readLd() {
    return 0.002f;
}
float readLq() {
    return 0.002f;
}
// 标准值
float readR() {
    return 1.0f;
}
float readFai() {
    return 0.175f;
}
