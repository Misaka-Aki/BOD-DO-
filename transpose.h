#ifndef _TRANSPOSE_H_
#define _TRANSPOSE_H_ 
#define MAXSIZE 10

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

#endif

