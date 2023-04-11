#include "matrix.h"
#include "pmsm.h"

static float xek_1k_1[4][1];
static float Kk[4][2];
static float Pk_1k_1[4][4];
static float Qd[4][4];
static float Rd[2][2];

// void EKF_init(float xek_1k_1[4][1],float Pk_1k_1[4][4],\
//                 float Kk[4][2],float Qd[4][4],float Rd[2][2]){
void EKF_init() {
    xek_1k_1[0][0] = readId();
    xek_1k_1[1][0] = readIq();
    xek_1k_1[2][0] = readLq();
    xek_1k_1[3][0] = readLd();
    matrix_Zero((float*)Kk, 4, 2);
    matrix_Zero((float*)Pk_1k_1, 4, 4);

    matrix_Zero((float*)Qd, 4, 4);
    matrix_Zero((float*)Rd, 2, 2);
    Qd[0][0] = 10.0f;
    Qd[1][1] = 10.0f;
    Qd[2][2] = 0.000002f;
    Qd[3][3] = 0.000002f;
    Rd[0][0] = 50.0f;
    Rd[1][1] = 50.0f;
}
// void EKF_update(float xek_1k_1[4][1],float Pk_1k_1[4][4],float To){
void EKF_update(float To) {
    int status;
    float id = readId();
    float iq = readIq();
    float ud = readUd();
    float uq = readUq();
    float we = readWe();
    float R = readR();
    float fai = readFai();
    float Ld = xek_1k_1[3][0];
    float Lq = xek_1k_1[2][0];

    float Fk_1[4][4];
    float Fk_1T[4][4];
    float Hk[2][4];
    float HkT[4][2];
    float xekk_1[4][1];
    float xekk[4][1];
    float Pkk_1[4][4];

    float tmp44_1[4][4];
    float tmp44_2[4][4];
    float tmp44_3[4][4];
    float tmp21_1[2][1];
    float tmp21_2[2][1];
    float tmp22_1[2][2];
    float tmp22_2[2][2];
    float tmp24_1[2][4];
    float tmp41_1[4][1];
    float tmp42_1[4][2];

    matrix_Zero((float*)Fk_1, 4, 4);
    Fk_1[0][0] = -R / Ld;
    Fk_1[0][1] = Lq * we / Ld;
    Fk_1[0][2] = we * iq / Ld;
    Fk_1[0][3] = (R * id - Lq * we * iq - ud) / Ld / Ld;
    Fk_1[1][0] = -Ld * we / Lq;
    Fk_1[1][1] = -R / Lq;
    Fk_1[1][2] = (R * iq + we * Ld * id + we * fai - uq) / Lq / Lq;
    Fk_1[1][3] = -we * id / Lq;
    matrix_tanspose((float*)Fk_1, (float*)Fk_1T, 4, 4);  // Fk_1T = Fk_1'

    matrix_Zero((float*)Hk, 2, 4);
    Hk[0][0] = 1.0f;
    Hk[1][1] = 1.0f;
    matrix_tanspose((float*)Hk, (float*)HkT, 2, 4);

    xekk_1[0][0] =
        xek_1k_1[0][0] + (-R * id / Ld + Lq * we * iq / Ld + ud / Ld) * To;
    xekk_1[1][0] =
        xek_1k_1[1][0] +
        (-R * iq / Lq - Ld * we * id / Lq - we * fai / Lq + uq / Lq) * To;
    xekk_1[2][0] = xek_1k_1[2][0];
    xekk_1[3][0] = xek_1k_1[3][0];

    // Pkk_1  = Pk_1k_1 + (Fk_1*Pk_1k_1 + Pk_1k_1*Fk_1')*To + Qd;
    matrix_mul((float*)Fk_1, (float*)Pk_1k_1, (float*)tmp44_1, 4, 4, 4);
    matrix_mul((float*)Pk_1k_1, (float*)Fk_1T, (float*)tmp44_2, 4, 4, 4);
    matrix_add((float*)tmp44_1, (float*)tmp44_2, (float*)tmp44_3, 4, 4);
    matrix_const_mul((float*)tmp44_3, To, (float*)tmp44_1, 4,
                     4);  // tmp44_1 = (Fk_1*Pk_1k_1 + Pk_1k_1*Fk_1')*To
    matrix_add((float*)Pk_1k_1, (float*)tmp44_1, (float*)tmp44_2, 4, 4);
    matrix_add((float*)tmp44_2, (float*)Qd, (float*)Pkk_1, 4, 4);

    // xekk   = xekk_1 + Kk*([id;iq] - Hk*xekk_1);
    matrix_mul((float*)Hk, (float*)xekk_1, (float*)tmp21_1, 2, 4, 1);
    tmp21_2[0][0] = id - tmp21_1[0][0];
    tmp21_2[1][0] = iq - tmp21_1[1][0];
    matrix_mul((float*)Kk, (float*)tmp21_2, (float*)tmp41_1, 4, 2, 1);  // tmp41_1 = Kk*([id;iq] - Hk*xekk_1)
    matrix_add((float*)xekk_1, (float*)tmp41_1, (float*)xekk, 4, 1);

    // Pkk    = Pkk_1 - Kk*Hk*Pkk_1;
    // Kk     = Pkk_1*Hk'/(Hk*Pkk_1*Hk' + Rd);
    matrix_mul((float*)Hk, (float*)Pkk_1, (float*)tmp24_1, 2, 4, 4);
    matrix_mul((float*)Kk, (float*)tmp24_1, (float*)tmp44_1, 4, 2, 4);
    matrix_sub((float*)Pkk_1, (float*)tmp44_1, (float*)Pk_1k_1, 4,4);  // Pk_1k_1 = Pkk

    matrix_mul((float*)Hk, (float*)Pkk_1, (float*)tmp24_1, 2, 4, 4);
    matrix_mul((float*)tmp24_1, (float*)HkT, (float*)tmp22_1, 2, 4, 2);
    matrix_add((float*)tmp22_1, (float*)Rd, (float*)tmp22_2, 2, 2);  // tmp22_2 = (Hk*Pkk_1*Hk' + Rd);

    matrix_mul((float*)Pkk_1, (float*)HkT, (float*)tmp42_1, 4, 4, 2);
    status = matrix_div((float*)tmp42_1, (float*)tmp22_2, (float*)Kk, 4, 2);
    if (status < 0) {
        matrix_Zero((float*)Kk, 4, 2);       // clear
        matrix_Zero((float*)Pk_1k_1, 4, 4);  // clear
        printf("Inv fault\n");
        return;
    }  // 数据不对不更新，直接等下一次数据

    // xek_1k_1   = xekk;
    matrix_copy((float*)xekk, (float*)xek_1k_1, 4, 1);
}

float getEKFLd() {
    return xek_1k_1[3][0];
}
float getEKFLq() {
    return xek_1k_1[2][0];
}
float getEKFId() {
    return xek_1k_1[0][0];
}
float getEKFIq() {
    return xek_1k_1[1][0];
}

// on pc
float testEKF(const char* inFile) {
    int i = 0;
    char outFile[40];
    strcpy(outFile, inFile);
    strcat(outFile, "_EKF.txt");
    openxk6(inFile);
    EKF_init();

    while ((i = readxk6()) < 1400 && i >= 0) {
        printf("%lf ", readId());
        printf("%lf ", readIq());
        printf("%lf ", readUd());
        printf("%lf ", readUq());
        printf("%lf : ", readWe());
        if (i > 1) {
            EKF_update(0.0002f);
        }
        printf("%lf ", getEKFLd());
        printf("%lf\n", getEKFLq());
        save2output(4, getEKFLd(), getEKFLq(), getEKFId(), getEKFIq());
    }
    output2fd(outFile, i, 4);
    closexk6();
    return -1;
}