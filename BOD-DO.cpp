#include "stdio.h"
#include "stdlib.h"
#include "iostream"
#include "math.h"

using namespace std;

#define MAXSIZE 10

int N  ;
double T;
double Q0 ;
double L0 ;
double O0;
double Q[MAXSIZE];
double L[MAXSIZE];
double O[MAXSIZE];
double Q3[MAXSIZE];
double k1[MAXSIZE];
double k2[MAXSIZE];
double t[MAXSIZE];
double Os;

double Q1[MAXSIZE];
double L1[MAXSIZE];
double O1[MAXSIZE];
double Q2[MAXSIZE];
double L2[MAXSIZE];
double O2[MAXSIZE];
double aa[MAXSIZE];
double bb[MAXSIZE];
double r[MAXSIZE];
double s[MAXSIZE];

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


double _L[MAXSIZE][MAXSIZE];
double _O[MAXSIZE][MAXSIZE];
double _f[MAXSIZE][MAXSIZE];
double _g[MAXSIZE][MAXSIZE];
double _h[MAXSIZE][MAXSIZE];


double A1[MAXSIZE][MAXSIZE]; 
double C1[MAXSIZE][MAXSIZE]; 


void Transpose2(double a[][MAXSIZE], double b[][MAXSIZE], int row, int col);
void Transpose1to2(double a[MAXSIZE], double b[][MAXSIZE]);


void displayMatrix(double a[][MAXSIZE], int row, int col);

void displayMatrix(double a[MAXSIZE], int col);


bool Gauss(double A[][MAXSIZE], double B[][MAXSIZE], int n);

void Matrix_Mult(double A[][MAXSIZE], double B[][MAXSIZE], double C[][MAXSIZE],
	int row1, int col1, int row2, int col2);

void MatrixV_Mult(double A[][MAXSIZE], double B[MAXSIZE], double C[MAXSIZE],
	int row1, int col1, int row2);


int main()
{
	int i, j;
	printf("请输入断面数x=");
	scanf("%d", &N);
	printf("请输入各断面进入河流的流量，m3/s  Q=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &Q[i]);

	printf("请输入各断面处的取水量，m3/s  Q3=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &Q3[i]);
	
	printf("请输入各断面进入河流的BOD浓度，mg/L  L=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &L[i]);

	Transpose1to2(L, _L);

	printf("请输入各断面进入河流的DO浓度，mg/L  O=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &O[i]);

	Transpose1to2(O, _O);

	printf("请输入各断面下游河段的BOD衰减速率常数（耗氧系数），d-1  k1=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &kd[i]);

	printf("请输入各断面下游河段的大气复氧速率常数（复氧系数），d-1   k2=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &ka[i]);

	printf("请输入各断面下游河段内的流行时间，d  t=\n");
	for (i = 0; i < N; i++)
		scanf("%lf", &t[i]);

	printf("请输入上游（起始）河段的流量，m3/s  Q0=\n");
	scanf("%lf", &Q2[0]);

	printf("请输入上游（起始）河段的BOD浓度，mg/L  L0=\n");
	scanf("%lf", &L20);

	printf("请输入上游（起始）河段的DO浓度，mg/L  O0=\n");
	scanf("%lf", &O20);

	printf("请输入温度值，℃  T=\n");
	scanf("%lf", &T);
	Os = 468 / (31.6 + T);
	printf("饱和溶解氧值，mg/L  Os= %.2lf  \n", Os);

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
		r[i] = exp(-ka[i] * t[i]);
		c[i] = (Q1[i + 1] - Q3[i + 1]) * r[i] / Q2[i + 1];
		bb[i] = kd[i] * (aa[i] - r[i]) / (ka[i] - kd[i]);
		d[i] = bb[i]*(Q1[i + 1] - Q3[i + 1]) / Q2[i + 1];
		s[i] = Os*(1 - r[i]);
		f[i] = s[i]*(Q1[i + 1] - Q3[i + 1]) / Q2[i + 1];
	}


	Transpose1to2(f, _f);

	for ( i = 0; i < N; i++)
		g[i] = 0;
	g[0] = aa[0] * L20;

	Transpose1to2(g, _g);

	for ( i = 0; i < N; i++)
		h[i] = 0;
	h[0] = c[0] * O20 - d[0] * L20;
	h[0] = r[0] * O20 - kd[0] / (ka[0] - kd[0]) * (aa[0] - r[0]);

	Transpose1to2(h, _h);



	for ( i = 0; i < MAXSIZE; i++)
		A[i][i] = 1;
	for (i = 0; i < N; i++)
		A[i + 1][i] = -a[i + 1];

	for (i = 0; i < N; i++)
		B[i][i] = b[i + 1];

	for (i = 0; i < MAXSIZE; i++)
		C[i][i] = 1;
	for (i = 0; i < N-1; i++)
		C[i + 1][i] = -c[i + 1];

	D[N][N] = 0;
	for (i = 0; i < N - 1; i++)
		D[i + 1][i] = d[i + 1];
	for (i = 0; i < N - 1; i++)
		D[i + 1][i] = kd[i+1] / (ka[i + 1] - kd[i + 1]) * (aa[i + 1] - r[i + 1]);


	Gauss(A, A1, N);
	Gauss(C, C1, N);
	Matrix_Mult(A1, B, U, N, N, N, N); 

	double temp1[MAXSIZE][MAXSIZE] = { 0 };
	double temp2[MAXSIZE][MAXSIZE] = { 0 };
	double temp3[MAXSIZE][MAXSIZE] = { 0 };
	double temp4[MAXSIZE][MAXSIZE] = { 0 };

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			C1[i][j] *= -1;
	}
	Matrix_Mult(C1, D, temp1, N, N, N, N);

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			C1[i][j] *= -1;
	}
	Matrix_Mult(temp1, A1, temp2, N, N, N, N);
	Matrix_Mult(temp2, B, V, N, N, N, N); 

	
	MatrixV_Mult(A1, g, m, N, N, N); 

	double C1B[MAXSIZE][MAXSIZE] = { 0 };
	double C1BO[MAXSIZE] = { 0 };
	double fg[MAXSIZE] = { 0 };
	double C1fg[MAXSIZE] = { 0 };
	double C1D[MAXSIZE][MAXSIZE] = { 0 };
	double C1DA1[MAXSIZE][MAXSIZE] = { 0 };
	double C1DA1g[MAXSIZE] = { 0 };

	Matrix_Mult(C1, B, C1B, N, N, N, N);
	MatrixV_Mult(C1B, O, C1BO, N, N, N);

	for (i = 0; i < N; i++)
		fg[i] = f[i] + g[i];
	MatrixV_Mult(C1, fg, C1fg, N, N, N);

	Matrix_Mult(C1, D, C1D, N, N, N, N);
	Matrix_Mult(C1D, A1, C1DA1, N, N, N, N);
	MatrixV_Mult(C1DA1, g, C1DA1g, N, N, N);



	for (i = 0; i < N; i++)
		n[i] = C1BO[i] + C1fg[i] - C1DA1g[i];      
		
	double UL[MAXSIZE] = { 0 };

	MatrixV_Mult(U, L, UL, N, N, N);

	for (i = 0; i < N; i++)
		L2[i] = UL[i] + m[i];            
	double VL[MAXSIZE] = { 0 };

	MatrixV_Mult(V, L, VL, N, N, N);
	for (i = 0; i < N; i++)
		O2[i] = VL[i] + n[i];              
	printf(" U   河流BOD稳态响应矩阵\n");
	displayMatrix(U, N, N);
	
	printf(" V   河流DO稳态响应矩阵\n");
	displayMatrix(V, N, N);
	
	printf(" m\n");
	displayMatrix(m, N);
	
	printf(" n\n");
	displayMatrix(n, N);
	
	printf(" L2   各断面的BOD浓度\n");
	displayMatrix(L2, N);
	
	printf(" O2   各断面的DO浓度\n");
	displayMatrix(O2, N);


	return 0;
}



void Matrix_Mult(double A[][MAXSIZE], double B[][MAXSIZE], double C[][MAXSIZE],
	int row1, int col1, int row2, int col2)
{
	int i = 0, j = 0, k = 0;
	for (i = 0; i < row1; i++)
	{
		for (j = 0; j < col2; ++j)
		{
			C[i][j] = 0.0;
			for (k = 0; k < col1; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void MatrixV_Mult(double A[][MAXSIZE], double B[MAXSIZE], double C[MAXSIZE],
	int row1, int col1, int row2)
{
	int i, j;
	for (i = 0; i < row1; i++)
	{
		C[i] = 0;
		for (j = 0; j < col1; j++)
			C[i] += A[i][j] * B[j];
	}
}



bool Gauss(double A[][MAXSIZE], double B[][MAXSIZE], int n)
{
	int i, j, k;
	double max, temp;
	double t[MAXSIZE][MAXSIZE];                
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			t[i][j] = A[i][j];
		}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			B[i][j] = (i == j) ? (float)1 : 0;
		}
	}
	for (i = 0; i < n; i++)
	{

		max = t[i][i];
		k = i;
		for (j = i + 1; j < n; j++)
		{
			if (fabs(t[j][i]) > fabs(max))
			{
				max = t[j][i];
				k = j;
			}
		}

		if (k != i)
		{
			for (j = 0; j < n; j++)
			{
				temp = t[i][j];
				t[i][j] = t[k][j];
				t[k][j] = temp;

				temp = B[i][j];
				B[i][j] = B[k][j];
				B[k][j] = temp;
			}
		}

		if (t[i][i] == 0)
		{
			printf("There is no inverse matrix!\n");
			return false;
		}
		
		temp = t[i][i];
		for (j = 0; j < n; j++)
		{
			t[i][j] = t[i][j] / temp;        
			B[i][j] = B[i][j] / temp;        
		}
		for (j = 0; j < n; j++)        
		{
			if (j != i)                
			{
				temp = t[j][i];
				for (k = 0; k < n; k++)        
				{
					t[j][k] = t[j][k] - t[i][k] * temp;
					B[j][k] = B[j][k] - B[i][k] * temp;
				}
			}
		}
	}
	return true;
}



void displayMatrix(double a[][MAXSIZE], int row, int col)
{
	
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			printf("\t%.4lf", a[i][j]);
		}
		printf("\n\n");
	}
}


void displayMatrix(double a[MAXSIZE], int col)
{
	
	for (int j = 0; j < col; j++)
	{
			printf("\t%.4lf", a[j]);
	}
	printf("\n\n");
}



void Transpose(double *m1, double *m2, int m, int n)
{
	int i = 0, j;
	for (int i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
			m2[j*m + i] = m1[i*n + j];
	}
}

void Transpose2(double a[][MAXSIZE], double b[][MAXSIZE], int row, int col)
{
	int i = 0, j;
	for (int i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			b[j][i] = a[i][j];
		}
	}
}

void Transpose1to2(double a[MAXSIZE], double b[MAXSIZE][MAXSIZE])
{
	int i = 0;
	for (int i = 0; i < MAXSIZE; i++)
	{
		b[0][i] = a[i];
	}
}


