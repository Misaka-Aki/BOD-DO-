#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "gauss.h"
#include "display.h"
#include "transpose.h"
#include "mult.h"
using namespace std;
#define MAXSIZE 10   //���������� Ϊ 10


int N  ;// ��������
double T;//�¶�
double Q0 ;// ���Σ���ʼ���Ӷε�������m3 / s��
double L0 ;// ���Σ���ʼ���Ӷε�BODŨ�ȣ�mg / L��
double O0;// ���Σ���ʼ���Ӷε�DOŨ�ȣ�mg / L��
double Q[MAXSIZE];// ����������������ˮ������m3 / s��
double L[MAXSIZE];// ��������������BODŨ�ȣ�mg / L��
double O[MAXSIZE];// ��������������DOŨ�ȣ�mg / L��
double Q3[MAXSIZE];// �����洦��ȡˮ����m3 / s��
double k1[MAXSIZE];// ���������κӶε�BOD˥�����ʳ���������ϵ������d - 1��
double k2[MAXSIZE];// ���������κӶεĴ����������ʳ���������ϵ������d - 1��
double t[MAXSIZE];// ���������κӶ��ڵ�����ʱ�䣬d��
double Os;// �����ܽ���ֵ��mg / L��
//��������
double Q1[MAXSIZE];// ���ν������i��������m3 / s��
double L1[MAXSIZE];// ���ν������i��BODŨ�ȣ�mg / L��
double O1[MAXSIZE];// ���ν������i��DOŨ�ȣ�mg / L��
double Q2[MAXSIZE];// �ɶ���i��������ε�������m3 / s��
double L2[MAXSIZE];// �ɶ���i��������ε�BODŨ�ȣ�mg / L��
double O2[MAXSIZE];// �ɶ���i��������ε�DOŨ�ȣ�mg / L��
double aa[MAXSIZE];//��[MAXSIZE]��
double bb[MAXSIZE];//��[MAXSIZE]��
double r[MAXSIZE];//��[MAXSIZE]��
double s[MAXSIZE];//��[MAXSIZE]��

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

//ת�õľ���  ������û���õ� 
double _L[MAXSIZE][MAXSIZE];
double _O[MAXSIZE][MAXSIZE];
double _f[MAXSIZE][MAXSIZE];
double _g[MAXSIZE][MAXSIZE];
double _h[MAXSIZE][MAXSIZE];

//�����
double A1[MAXSIZE][MAXSIZE]; //A����������
double C1[MAXSIZE][MAXSIZE]; //C����������

int main()
{
	int i, j;
	
	printf("�����������x=");
	scanf("%d", &N);
	
	printf("�������������������������m3/s  Q=\n");
	for (i = 1; i <= N; i++)
		scanf("%lf", &Q[i]);

	printf("����������洦��ȡˮ����m3/s  Q3=\n");
	for (i = 1; i <= N; i++)
		scanf("%lf", &Q3[i]);
	
	printf("�������������������BODŨ�ȣ�mg/L  L=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &L[i]);

	Transpose1to2(L, _L);

	printf("�������������������DOŨ�ȣ�mg/L  O=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &O[i]);

	Transpose1to2(O, _O);

	printf("���������κӶε�BOD˥�����ʳ���������ϵ������d-1  k1=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &kd[i]);

	printf("���������κӶεĴ����������ʳ���������ϵ������d-1   k2=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &ka[i]);

	printf("���������κӶ��ڵ�����ʱ�䣬d  t=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &t[i]);

	printf("���Σ���ʼ���Ӷε�������m3/s  Q0=\n");
	scanf("%lf", &Q2[0]);

	printf("���Σ���ʼ���Ӷε�BODŨ�ȣ�mg/L  L0=\n");
	scanf("%lf", &L20);

	printf("���Σ���ʼ���Ӷε�DOŨ�ȣ�mg/L  O0=\n");
	scanf("%lf", &O20);

	printf("�¶�ֵ����  T=\n");
	scanf("%lf", &T);
	
	
	Os = 468 / (31.6 + T);
	printf("�����ܽ���ֵ��mg/L  Os= %.2lf  \n", Os);


/*1.BOD��Ӷξ���ģ��*/
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
		
		
		
		/*2.BOD-DO��Ͼ���ģ��*/
		r[i] = exp(-ka[i] * t[i]);
		bb[i] = kd[i] * (aa[i] - r[i]) / (ka[i] - kd[i]);
		
		c[i] = (Q1[i + 1] - Q3[i + 1]) * r[i] / Q2[i + 1];
		
		d[i] = bb[i]*(Q1[i + 1] - Q3[i + 1]) / Q2[i + 1];
		
		s[i] = Os*(1 - r[i]);
		
		f[i] = s[i]*(Q1[i + 1] - Q3[i + 1]) / Q2[i + 1];
	}

	/*fת��*/ 
	Transpose1to2(f, _f);

	for ( i = 0; i < N; i++)
		g[i] = 0;
		
	/*gת��*/ 
	g[0] = a[0] * L20;
	Transpose1to2(g, _g);

	for ( i = 0; i < N; i++)
		h[i] = 0;
	h[0] = c[0] * O20 - d[0] * L20;

    /*hת��*/
	Transpose1to2(h, _h);


	/*����A*/ 
	for ( i = 0; i < MAXSIZE; i++)
		A[i][i] = 1;
	for (i = 0; i < N; i++)
		A[i + 1][i] = -a[i + 1];
		
		
	/*����B*/
	for (i = 0; i < N; i++)
		B[i][i] = b[i + 1];
		
	/*����C*/
	for (i = 0; i < MAXSIZE; i++)
		C[i][i] = 1;
	for (i = 0; i < N-1; i++)
		C[i + 1][i] = -c[i + 1];

	/*����D*/
	D[N][N] = 0;
	for (i = 0; i < N - 1; i++)
		D[i + 1][i] = d[i + 1];


	/*����A��C����*/
	Gauss(A, A1, N);
	Gauss(C, C1, N);
	
	/*����U*/
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

	/*����m*/
	V_Mult(A1, g, m, N, N, N);          


	double C1B[MAXSIZE][MAXSIZE] = { 0 };
	double C1BO[MAXSIZE] = { 0 };
	double fh[MAXSIZE] = { 0 };
	double C1fh[MAXSIZE] = { 0 };
	double C1D[MAXSIZE][MAXSIZE] = { 0 };
	double C1DA1[MAXSIZE][MAXSIZE] = { 0 };
	double C1DA1g[MAXSIZE] = { 0 };

	/*���C1B*/
	Mult(C1, B, C1B, N, N, N, N);
	/*���C1BO*/
	V_Mult(C1B, O, C1BO, N, N, N);
	
	/*��f+h��*/
	for (i = 0; i < N; i++)
		fh[i] = f[i] + h[i];
	
	/*���C1fh*/
	V_Mult(C1, fh, C1fh, N, N, N);

	/*�ֲ������C1DA1  */
	Mult(C1, D, C1D, N, N, N, N);
	Mult(C1D, A1, C1DA1, N, N, N, N);
	
	/*��C1DA1g*/
	V_Mult(C1DA1, g, C1DA1g, N, N, N);


	/*����n*/
	for (i = 0; i < N; i++)
		n[i] = C1BO[i] + C1fh[i] - C1DA1g[i];     


	double UL[MAXSIZE] = { 0 };

	/*����UL*/
	V_Mult(U, L, UL, N, N, N);
	
	/*����L2*/
	for (i = 0; i < N; i++)
		L2[i] = UL[i] + m[i];    

	/*����VL*/
	double VL[MAXSIZE] = { 0 };
	V_Mult(V, L, VL, N, N, N);
	
	/*����O2*/
	for (i = 0; i < N; i++)
		O2[i] = VL[i] + n[i];



	printf(" U   ����BOD��̬��Ӧ����\n");       
	display(U, N, N);
	
	printf(" V   ����DO��̬��Ӧ����\n"); 
	display(V, N, N);
	
	printf(" m\n")  ; 
	display(m, N);
	
	printf(" n\n")  ; 
	display(n, N);
	
	printf(" L2   �������BODŨ��\n"); 
	display(L2, N);
	
	printf(" O2   �������DOŨ��\n"); 
	display(O2, N);


	return 0;
}




