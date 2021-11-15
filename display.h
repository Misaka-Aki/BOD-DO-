#ifndef _DISPLAY_H_
#define _DISPLAY_H_ 
#define MAXSIZE 10

//显示矩阵
void display(double a[][MAXSIZE], int row, int col)
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

//显示向量
void display(double a[MAXSIZE], int col)
{
	for (int j = 0; j < col; j++)
	{
			printf("\t%.4lf", a[j]);
	}
	printf("\n\n");
}

#endif
