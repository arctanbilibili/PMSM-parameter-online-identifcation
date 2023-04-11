#include <stdio.h>
#include "matrix.h"

// 1、
// RLS辨识：
// 设置初始值RLS_init();
// 中断内循环调用RLS_update(float To);//To为中断间隔时间

// 2、
// EKF辨识：
// 设置初始值EKF_init();
// 中断内循环调用EKF_update(float To);//To为中断间隔时间

// 3、
// MARS辨识：
// 设置初始值MARS_init();
// 中断内循环调用MARS_update(float To);//To为中断间隔时间

/* 
* 移植只需修改pmsm.c和pmsm.h文件
* （这两个文件模拟电机运行，实际上是读取 xk6.txt文件中的{id ,iq ,ud ,uq ,we}来模拟每个时刻的测量值）
* 以及删除三个方法c文件内的testRLS，testEKF，testMARS函数即可
* 采集数据函数readxk6()到xk6[5]数组以及xk_1[5]数组
* xk6 为当前时刻的 {id ,iq ,ud ,uq ,we} 5个数据
*/

extern float testRLS(const char* inFile);
extern float testEKF(const char* inFile);
extern float testMARS(const char* inFile);

int main() {
    // 测试，xk6.txt 每行前5列代表数据 {id ,iq ,ud ,uq ,we}
    // testRLS("../xk6.txt");
    // testEKF("../xk6.txt");
    // testMARS("../xk6.txt");

    return 0;
}
