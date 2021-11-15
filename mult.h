#ifndef _MULT_H_
#define _MULT_H_ 
#define MAXSIZE 10

//æÿ’Û≥À“‘æÿ’Û
void Mult(double A[][MAXSIZE], double B[][MAXSIZE], double C[][MAXSIZE],
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

//æÿ’Û≥À“‘œÚ¡ø
void V_Mult(double A[][MAXSIZE], double B[MAXSIZE], double C[MAXSIZE],
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

#endif
