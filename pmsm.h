#ifndef _PMSM_H_
#define _PMSM_H_
#include <stdio.h>
#include <string.h>
// extern float Ld;
// extern float Lq;
// extern float R;
// extern float Fai;

// fs
int openxk6(const char* fileName);
int readxk6();
int closexk6();
int save2output(int arg, ...);
int output2fd(const char* fileName, int r, int c);

float readId();
float readIq();
float readId_1();
float readIq_1();
float readUd();
float readUq();
float readLd();
float readLq();
float readR();
float readFai();
float readWe();

#endif
