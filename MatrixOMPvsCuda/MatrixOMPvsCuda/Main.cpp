#include "Main.h"

static double tim1;
static double tim2;

Main *cu = new Main();

double **mathLine(double **Matrix)
{
	double determinant = 1; double t;
	double end, start = omp_get_wtime();
	for (int i = 0; i < SIZE - 1; i++)
	{
		int maxN = i;
		double maxValue = fabs(Matrix[i][i]);
		for (int j = i + 1; j < SIZE; j++) //поиск максимального элемента
		{
			double tmp = fabs(Matrix[j][i]);
			if (tmp > maxValue)
			{
				maxN = j;
				maxValue = tmp;
			}
		}
		if (maxN > i)
		{
			for (int j = 0; j < SIZE; j++)
			{
				double tmp = Matrix[i][j];
				Matrix[i][j] = Matrix[maxN][j];
				Matrix[maxN][j] = tmp;
			}
		}
		else if (maxValue == 0)
		{
			printf("null in matrix %d\n",i);
			//return NULL;
		}

		double val = Matrix[i][i];
		determinant *= val;

		for (int j = i + 1; j < SIZE; j++)
		{
			double k = Matrix[j][i] / val;
			Matrix[j][i] = 0;
			for (int c = i + 1; c < SIZE; c++)
				Matrix[j][c] = Matrix[j][c] - Matrix[i][c] * k;
		}
	}
	end = omp_get_wtime();
	t = (end - start); tim1 = t;
	cu->timelinever = t;
	printf("time: %f\n", t);
	printf("determinant - %f\n", fabs(determinant));
	if (SIZE < 10)
	{
		for (int i = 0; i < SIZE; i++)
		{
			for (int j = 0; j < SIZE; j++)
			{
				printf("%.2f\t", Matrix[i][j]);
			}
			printf("\n");
		}
	}
	return (double**)Matrix;
}

double **mathOmp(double **Matrix)
{
	omp_set_num_threads(20);
	int ndata = SIZE / omp_get_num_threads();
	double determinant = 1, k, t, maxValue, tmp;
	int pol = SIZE / (4), i, j, c, maxN;
	double end, start = omp_get_wtime();
	for (i = 0; i < SIZE - 1; i++)
	{
#pragma omp parallel private (j,tmp) shared(maxN, maxValue)
		{
#pragma omp sections
		{
#pragma omp section
			maxN = i;
#pragma omp section
			maxValue = fabs(Matrix[i][i]);
		}
#pragma omp for schedule(dynamic,pol)
		for (j = i + 1; j < SIZE; j++) //поиск максимального элемента
		{
			tmp = fabs(Matrix[j][i]);
			if (tmp > maxValue)
			{
				maxN = j;
				maxValue = tmp;
			}
		}
		}
#pragma omp parallel for private(j,tmp) if (maxN > i)
		for (j = 0; j < SIZE; j++)
		{
			tmp = Matrix[i][j];
			Matrix[i][j] = Matrix[maxN][j];
			Matrix[maxN][j] = tmp;
		}
		determinant *= Matrix[i][i];
#pragma omp parallel for private(k,c,j)
		for (j = i + 1; j < SIZE; j++)
		{
			k = Matrix[j][i] / Matrix[i][i];
			Matrix[j][i] = 0;
			for (c = i + 1; c < SIZE; c++)
				Matrix[j][c] = Matrix[j][c] - Matrix[i][c] * k;
		}
	}

	end = clock();
	t = (end - start); tim2 = t;
	cu->timeompver = t;
	printf("time: %f\n", t);
	printf("determinant - %3.3f\n", fabs(determinant));
	if (SIZE < 10)
	{
		for (int i = 0; i < SIZE; i++)
		{
			for (int j = 0; j < SIZE; j++)
			{
				printf("%.2f\t", Matrix[i][j]);
			}
			printf("\n");
		}
	}
	return (double**)Matrix;
}

void provMat(double **Matrix, double **Matrix2, double **Matrix3)
{
	printf("-------prov-------\n");
#pragma omp parallel for
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			if (fabs(Matrix[i][j] - Matrix2[i][j]) < 0.01)
			{

			}
			if (Matrix[i][j] != Matrix2[i][j])
			{
				printf("Matrix[%d][%d] - error", i, j);
			}
		}
	}
}

int main()
{
#ifndef _OPENMP
	printf("This compiled code has no OpenMP support:( Check your compiler: if it supports OpenMP you must apply a correspondent compiler key.\n");
	goto exit;
#endif
	double **Matrix = new double*[SIZE];// (double **)calloc(SIZE, sizeof(*Matrix));
	//double **Matrix2 = new double*[SIZE];// (double **)calloc(SIZE, sizeof(*Matrix2));
	double **Matrix3 = new double*[SIZE];// (double **)calloc(SIZE, sizeof(*Matrix3));
	srand(time(NULL));
	for (int i = 0; i < SIZE; i++)
	{
		Matrix[i] = new double[SIZE];// (double*)calloc(SIZE, sizeof(*Matrix));
		//Matrix2[i] = new double[SIZE];// (double*)calloc(SIZE, sizeof(*Matrix2));
		Matrix3[i] = new double[SIZE];// (double*)calloc(SIZE, sizeof(*Matrix3));
		for (int j = 0; j < SIZE; j++)
		{
			double k = (100 + rand() % 10000)*0.01;
			Matrix[i][j] = k;
			//Matrix2[i][j] = k;
			Matrix3[i][j] = k;
		}
	}
	if (SIZE < 10)
	{
		for (int i = 0; i < SIZE; i++)
		{
			for (int j = 0; j < SIZE; j++)
			{
				printf("%.2f\t", Matrix3[i][j]);
			}
			printf("\n");
		}
	}

	printf("------line------\n");
	Matrix = mathLine(Matrix);
	if (Matrix == NULL)
		goto exit;
	/*printf("------omp------\n");
	Matrix2 = mathOmp(Matrix2);*/
	printf("-----cuda-start------\n");
	cu->getMatrixFromCuda(Matrix3,cu);

	printf("-----result-time-----\n");
	printf("line time = %f \n omp time = %f \n cuda time = %f \n", cu->timelinever, cu->timeompver, cu->timecudaver);

	//provMat(Matrix, Matrix2, Matrix3);

exit:
	system("pause");
}