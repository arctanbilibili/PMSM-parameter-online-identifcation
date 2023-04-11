#include "matrix.h"
#include "pmsm.h"

static float x_o[2][1];
static float RDivLq;
static float oneDivLd;
const static float Ki4 = 2454.5;  // Lq
const static float Ki5 = 674.5;   // Ld

void MARS_init() {
    x_o[0][0] = readId();
    x_o[1][0] = readIq();

    RDivLq = readR() / readLq();
    oneDivLd = 1.0f / readLd();
}
float getMARSLd() {
    if (oneDivLd < 10.0f)
        return 0.01;
    return 1.0f / oneDivLd;
}
float getMARSLq() {
    if (RDivLq < 10.0f)
        return 0.01;
    return readR() / RDivLq;
}

// % update
void MARS_update(float To) {
    float id = readId();
    float iq = readIq();
    float ud = readUd();
    float uq = readUq();
    float we = readWe();
    float R = readR();
    float fai = readFai();
    float Ld = getMARSLd();  // convert from oneDivLd
    float Lq = getMARSLq();  // convert from RDivLq

    float A_o[2][2];
    float B_o[2][2];
    float u_o[2];
    float x_o_next[2][1];  //[id_o,iq_o]
    float tmp21_1[2][1];
    float tmp21_2[2][1];
    float tmp21_3[2][1];

    A_o[0][0] = -R / Ld;
    A_o[0][1] = we * Lq / Ld;
    A_o[1][0] = -we * Ld / Lq;
    A_o[1][1] = -R / Lq;

    B_o[0][0] = 1.0f / Ld;
    B_o[0][1] = 0.0f;
    B_o[1][0] = 0.0f;
    B_o[1][1] = 1.0f / Lq;

    u_o[0] = readUd();
    u_o[1] = readUq() - we * fai;

    // RDivLq = RDivLq + Ki4*(iqk - iq_o(n))*iq_o(n)*To;
    // oneDivLd = oneDivLd + Ki5*(idk - id_o(n))*udk_1*To;

    matrix_mul((float*)A_o, (float*)x_o, (float*)tmp21_1, 2, 2, 1);
    matrix_mul((float*)B_o, (float*)u_o, (float*)tmp21_2, 2, 2, 1);
    matrix_add((float*)tmp21_1, (float*)tmp21_2, (float*)tmp21_3, 2, 1);
    matrix_const_mul((float*)tmp21_3, To, (float*)tmp21_1, 2, 1);
    matrix_add((float*)x_o, (float*)tmp21_1, (float*)x_o_next, 2, 1);  // x_oNext = x_o + (A_o*x_o + B_o*u_o)*To;
    RDivLq += Ki4 * (iq - x_o_next[1][0]) * x_o_next[1][0] * To;
    oneDivLd += Ki5 * (id - x_o_next[0][0]) * ud * To;
    matrix_copy((float*)x_o_next, (float*)x_o, 2, 1);
}

// on pc
float testMARS(const char* inFile) {
    int i = 0;
    char outFile[40];
    strcpy(outFile, inFile);
    strcat(outFile, "_MARS.txt");
    openxk6(inFile);
    MARS_init();

    while ((i = readxk6()) < 1400 && i >= 0) {
        printf("%lf ", readId());
        printf("%lf ", readIq());
        printf("%lf ", readUd());
        printf("%lf ", readUq());
        printf("%lf : ", readWe());
        if (i > 1) {
            MARS_update(0.0002f);
        }
        printf("%lf ", getMARSLd());
        printf("%lf\n", getMARSLq());
        save2output(2, getMARSLd(), getMARSLq());
    }
    output2fd(outFile, i, 2);
    closexk6();
    return -1;
}
