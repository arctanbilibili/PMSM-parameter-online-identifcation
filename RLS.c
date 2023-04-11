#include "matrix.h"
#include "pmsm.h"

const static float lamda = 0.995f;
static float theta[2][1];
static float Pk[2][2];

// 初始化4个参数值
void RLS_init() {
    matrix_Zero((float*)theta, 2, 1);
    matrix_I((float*)Pk, 2);
    theta[0][0] = readLq();
    theta[1][0] = readLd();
}
// 在线辨识更新两个参数
void RLS_update(float To) {
    float id_1 = readId_1();
    float iq_1 = readIq_1();
    float id = readId();
    float iq = readIq();
    float ud = readUd();
    float uq = readUq();
    float we = readWe();
    float R = readR();
    float fai = readFai();

    float yk[2][1];  // 涉及读取
    float eye2[2][2];
    float thetak_1[2][1];
    float Pk_1[2][2];
    float Faik[2][2];  // 涉及读取
    float err[2][1];
    float Kk[2][2];

    int status;
    float tmp22_1[2][2];
    float tmp22_2[2][2];
    float tmp22_3[2][2];
    float tmp21_1[2][1];

    matrix_copy((float*)theta, (float*)thetak_1, 2, 1);
    matrix_copy((float*)Pk, (float*)Pk_1, 2, 2);
    // 赋值Faik，yk
    Faik[0][0] = -we * iq;
    Faik[0][1] = (id - id_1) / To;
    Faik[1][0] = (iq - iq_1) / To;
    Faik[1][1] = we * id;
    yk[0][0] = ud - id * R;
    yk[0][1] = uq - iq * R - fai * we;
    // faik = [-w*iq did;
    //         diq w*id]';
    // yk = [ud-id*R;
    //         uq-iq*R-fai*w];
    matrix_I((float*)eye2, 2);
    matrix_const_mul((float*)eye2, lamda, (float*)tmp22_1, 2, 2);  // tmp22_1 = lamda*eye2
    matrix_tanspose((float*)Faik, (float*)tmp22_2, 2, 2);  // tmp22_2 = faik'
    matrix_mul((float*)tmp22_2, (float*)Pk_1, (float*)tmp22_3, 2, 2, 2);  // tmp22_3 = faik'*Pk_1
    matrix_mul((float*)tmp22_3, (float*)Faik, (float*)tmp22_2, 2, 2, 2);  // tmp22_2 = faik'*Pk_1*faik
    matrix_add((float*)tmp22_1, (float*)tmp22_2, (float*)tmp22_3, 2, 2);  // tmp22_3 = lamda*eye2+faik'*Pk_1*faik
    status = matrix_div((float*)Faik, (float*)tmp22_3, (float*)tmp22_1, 2, 2);  // tmp22_1 = Faik/(lamda*eye2+faik'*Pk_1*fai)
    if (status < 0) {
        printf("Inv fault\n");
        return;
    }  // 数据不对不更新，直接等下一次数据
    matrix_mul((float*)Pk_1, (float*)tmp22_1, (float*)Kk, 2, 2, 2);  // Kk = Pk_1*Faik/(lamda*eye2+faik'*Pk_1*fai)

    matrix_tanspose((float*)Faik, (float*)tmp22_1, 2, 2);  // tmp22_1 = faik'
    matrix_mul((float*)Kk, (float*)tmp22_1, (float*)tmp22_2, 2, 2, 2);  // tmp22_2 = kk*faik'
    matrix_sub((float*)eye2, (float*)tmp22_2, (float*)tmp22_3, 2, 2);  // tmp22_3 = eye2-kk*faik'
    matrix_mul((float*)tmp22_3, (float*)Pk_1, (float*)tmp22_1, 2, 2, 2);  // tmp22_1 = (eye2-kk*faik')*Pk_1
    matrix_const_mul((float*)tmp22_1, 1.0f / lamda, (float*)Pk, 2, 2);  // Pk = (eye2-kk*faik')*Pk_1/a

    matrix_tanspose((float*)Faik, (float*)tmp22_1, 2, 2);  // tmp22_1 = faik'
    matrix_mul((float*)tmp22_1, (float*)thetak_1, (float*)tmp21_1, 2, 2, 1);  // tmp21_1 = faik'*theta_1
    matrix_sub((float*)yk, (float*)tmp21_1, (float*)err, 2, 1);  // err = yk-faik'*theta_1

    matrix_mul((float*)Kk, (float*)err, (float*)tmp21_1, 2, 2, 1);  // tmp21_1 = Kk*err
    matrix_add((float*)thetak_1, (float*)tmp21_1, (float*)theta, 2, 1);  // theta = thetak_1+Kk*err
}

float getRLSLd() {
    return theta[1][0];
}
float getRLSLq() {
    return theta[0][0];
}

// on pc
float testRLS(const char* inFile) {
    int i = 0;
    char outFile[40];
    strcpy(outFile, inFile);
    strcat(outFile, "_RLS.txt");
    openxk6(inFile);
    RLS_init();

    while ((i = readxk6()) < 1400 && i >= 0) {
        printf("%lf ", readId());
        printf("%lf ", readIq());
        printf("%lf ", readUd());
        printf("%lf ", readUq());
        printf("%lf : ", readWe());
        if (i > 1) {
            RLS_update(0.0002f);
        }
        printf("%lf ", getRLSLd());
        printf("%lf\n", getRLSLq());
        save2output(2, getRLSLd(), getRLSLq());
    }
    output2fd(outFile, i, 2);
    closexk6();
    return -1;
}