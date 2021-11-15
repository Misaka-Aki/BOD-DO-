#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "gauss.h"
#include "display.h"
#include "transpose.h"
#include "mult.h"
using namespace std;
#define MAXSIZE 10   //矩阵最大阶数 为 10


int N  ;// 断面数；
double T;//温度
double Q0 ;// 上游（起始）河段的流量，m3 / s；
double L0 ;// 上游（起始）河段的BOD浓度，mg / L；
double O0;// 上游（起始）河段的DO浓度，mg / L；
double Q[MAXSIZE];// 各断面进入河流的污水流量，m3 / s；
double L[MAXSIZE];// 各断面进入河流的BOD浓度，mg / L；
double O[MAXSIZE];// 各断面进入河流的DO浓度，mg / L；
double Q3[MAXSIZE];// 各断面处的取水量，m3 / s；
double k1[MAXSIZE];// 各断面下游河段的BOD衰减速率常数（耗氧系数），d - 1；
double k2[MAXSIZE];// 各断面下游河段的大气复氧速率常数（复氧系数），d - 1；
double t[MAXSIZE];// 各断面下游河段内的流行时间，d；
double Os;// 饱和溶解氧值，mg / L；
//其他参数
double Q1[MAXSIZE];// 上游进入断面i的流量，m3 / s；
double L1[MAXSIZE];// 上游进入断面i的BOD浓度，mg / L；
double O1[MAXSIZE];// 上游进入断面i的DO浓度，mg / L；
double Q2[MAXSIZE];// 由断面i输出到下游的流量，m3 / s；
double L2[MAXSIZE];// 由断面i输出到下游的BOD浓度，mg / L；
double O2[MAXSIZE];// 由断面i输出到下游的DO浓度，mg / L；
double aa[MAXSIZE];//α[MAXSIZE]；
double bb[MAXSIZE];//β[MAXSIZE]；
double r[MAXSIZE];//γ[MAXSIZE]；
double s[MAXSIZE];//δ[MAXSIZE]；

double kd[MAXSIZE];
double ka[MAXSIZE];
double a[MAXSIZE];
double b[MAXSIZE];
double c[MAXSIZE];
double d[MAXSIZE];
double f[MAXSIZE];
double g[MAXSIZE];
double h[MAXSIZE];
double m[MAXSIZE];
double n[MAXSIZE];

double L20;
double O20;

double A[MAXSIZE][MAXSIZE];
double B[MAXSIZE][MAXSIZE];
double C[MAXSIZE][MAXSIZE];
double D[MAXSIZE][MAXSIZE];
double E[MAXSIZE][MAXSIZE];
double F[MAXSIZE][MAXSIZE];
double V[MAXSIZE][MAXSIZE];
double U[MAXSIZE][MAXSIZE];
double R[MAXSIZE][MAXSIZE];

//转置的矩阵  本程序没有用到 
double _L[MAXSIZE][MAXSIZE];
double _O[MAXSIZE][MAXSIZE];
double _f[MAXSIZE][MAXSIZE];
double _g[MAXSIZE][MAXSIZE];
double _h[MAXSIZE][MAXSIZE];

//逆矩阵
double A1[MAXSIZE][MAXSIZE]; //A矩阵的逆矩阵
double C1[MAXSIZE][MAXSIZE]; //C矩阵的逆矩阵

int main()
{
	int i, j;
	
	printf("请输入断面数x=");
	scanf("%d", &N);
	
	printf("请输入各断面进入河流的流量，m3/s  Q=\n");
	for (i = 1; i <= N; i++)
		scanf("%lf", &Q[i]);

	printf("请输入各断面处的取水量，m3/s  Q3=\n");
	for (i = 1; i <= N; i++)
		scanf("%lf", &Q3[i]);
	
	printf("请输入各断面进入河流的BOD浓度，mg/L  L=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &L[i]);

	Transpose1to2(L, _L);

	printf("请输入各断面进入河流的DO浓度，mg/L  O=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &O[i]);

	Transpose1to2(O, _O);

	printf("各断面下游河段的BOD衰减速率常数（耗氧系数），d-1  k1=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &kd[i]);

	printf("各断面下游河段的大气复氧速率常数（复氧系数），d-1   k2=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &ka[i]);

	printf("各断面下游河段内的流行时间，d  t=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &t[i]);

	printf("上游（起始）河段的流量，m3/s  Q0=\n");
	scanf("%lf", &Q2[0]);

	printf("上游（起始）河段的BOD浓度，mg/L  L0=\n");
	scanf("%lf", &L20);

	printf("上游（起始）河段的DO浓度，mg/L  O0=\n");
	scanf("%lf", &O20);

	printf("温度值，℃  T=\n");
	scanf("%lf", &T);
	
	
	Os = 468 / (31.6 + T);
	printf("饱和溶解氧值，mg/L  Os= %.2lf  \n", Os);


/*1.BOD多河段矩阵模型*/
	for (i = 0; i < N; i++)
	{
		Q2[i + 1] = Q2[i] + Q[i + 1] - Q3[i + 1];
		Q1[i + 1] = Q2[i];
	}


	for ( i = 0; i < N; i++)
	{
		aa[i] = exp(-kd[i] * t[i]);
		a[i] = aa[i] * (Q1[i + 1] - Q3[i + 1]) / Q2[i + 1];
		b[i + 1] = Q[i + 1] / Q2[i + 1];
		
		
		
		/*2.BOD-DO耦合矩阵模型*/
		r[i] = exp(-ka[i] * t[i]);
		bb[i] = kd[i] * (aa[i] - r[i]) / (ka[i] - kd[i]);
		
		c[i] = (Q1[i + 1] - Q3[i + 1]) * r[i] / Q2[i + 1];
		
		d[i] = bb[i]*(Q1[i + 1] - Q3[i + 1]) / Q2[i + 1];
		
		s[i] = Os*(1 - r[i]);
		
		f[i] = s[i]*(Q1[i + 1] - Q3[i + 1]) / Q2[i + 1];
	}

	/*f转置*/ 
	Transpose1to2(f, _f);

	for ( i = 0; i < N; i++)
		g[i] = 0;
		
	/*g转置*/ 
	g[0] = a[0] * L20;
	Transpose1to2(g, _g);

	for ( i = 0; i < N; i++)
		h[i] = 0;
	h[0] = c[0] * O20 - d[0] * L20;

    /*h转置*/
	Transpose1to2(h, _h);


	/*矩阵A*/ 
	for ( i = 0; i < MAXSIZE; i++)
		A[i][i] = 1;
	for (i = 0; i < N; i++)
		A[i + 1][i] = -a[i + 1];
		
		
	/*矩阵B*/
	for (i = 0; i < N; i++)
		B[i][i] = b[i + 1];
		
	/*矩阵C*/
	for (i = 0; i < MAXSIZE; i++)
		C[i][i] = 1;
	for (i = 0; i < N-1; i++)
		C[i + 1][i] = -c[i + 1];

	/*矩阵D*/
	D[N][N] = 0;
	for (i = 0; i < N - 1; i++)
		D[i + 1][i] = d[i + 1];


	/*矩阵A、C求逆*/
	Gauss(A, A1, N);
	Gauss(C, C1, N);
	
	/*矩阵U*/
	Mult(A1, B, U, N, N, N, N);      


	double temp1[MAXSIZE][MAXSIZE] = { 0 };
	double temp2[MAXSIZE][MAXSIZE] = { 0 };
	double temp3[MAXSIZE][MAXSIZE] = { 0 };
	double temp4[MAXSIZE][MAXSIZE] = { 0 };

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			C1[i][j] *= -1;
	}
	

	Mult(C1, D, temp1, N, N, N, N);

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			C1[i][j] *= -1;
	}
	
	Mult(temp1, A1, temp2, N, N, N, N);
	Mult(temp2, B, V, N, N, N, N);   

	/*矩阵m*/
	V_Mult(A1, g, m, N, N, N);          


	double C1B[MAXSIZE][MAXSIZE] = { 0 };
	double C1BO[MAXSIZE] = { 0 };
	double fh[MAXSIZE] = { 0 };
	double C1fh[MAXSIZE] = { 0 };
	double C1D[MAXSIZE][MAXSIZE] = { 0 };
	double C1DA1[MAXSIZE][MAXSIZE] = { 0 };
	double C1DA1g[MAXSIZE] = { 0 };

	/*求各C1B*/
	Mult(C1, B, C1B, N, N, N, N);
	/*求各C1BO*/
	V_Mult(C1B, O, C1BO, N, N, N);
	
	/*求（f+h）*/
	for (i = 0; i < N; i++)
		fh[i] = f[i] + h[i];
	
	/*求各C1fh*/
	V_Mult(C1, fh, C1fh, N, N, N);

	/*分步，求各C1DA1  */
	Mult(C1, D, C1D, N, N, N, N);
	Mult(C1D, A1, C1DA1, N, N, N, N);
	
	/*求C1DA1g*/
	V_Mult(C1DA1, g, C1DA1g, N, N, N);


	/*矩阵n*/
	for (i = 0; i < N; i++)
		n[i] = C1BO[i] + C1fh[i] - C1DA1g[i];     


	double UL[MAXSIZE] = { 0 };

	/*计算UL*/
	V_Mult(U, L, UL, N, N, N);
	
	/*计算L2*/
	for (i = 0; i < N; i++)
		L2[i] = UL[i] + m[i];    

	/*计算VL*/
	double VL[MAXSIZE] = { 0 };
	V_Mult(V, L, VL, N, N, N);
	
	/*计算O2*/
	for (i = 0; i < N; i++)
		O2[i] = VL[i] + n[i];



	printf(" U   河流BOD稳态响应矩阵\n");       
	display(U, N, N);
	
	printf(" V   河流DO稳态响应矩阵\n"); 
	display(V, N, N);
	
	printf(" m\n")  ; 
	display(m, N);
	
	printf(" n\n")  ; 
	display(n, N);
	
	printf(" L2   各断面的BOD浓度\n"); 
	display(L2, N);
	
	printf(" O2   各断面的DO浓度\n"); 
	display(O2, N);


	return 0;
}




