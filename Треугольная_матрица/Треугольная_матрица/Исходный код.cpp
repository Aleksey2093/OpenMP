#include "CudaInfo.cuh";

using namespace::std;

bool printMatrixToConsole(string title, double **matrix, int n, double diagol)
{
	cout << endl << "-----------------------------" << endl;
	cout << title << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << matrix[i][j] << "\t\t";
		cout << endl;
	}
	cout << "-----------------------------" << endl;
	cout << "Определитель " << diagol << endl;
	cout << "-----------------------------" << endl << endl;
	return true;
}

double math_Line(double **matrix, int n)
{
	cout << "Расчет последовательной версии" << endl;
	double t1 = omp_get_wtime(), t2;
	for (int i = 0; i < n - 1; i++)
	{
		int maxN = i;
		double maxValue = fabs(matrix[i][i]);
		for (int j = i + 1; j < n; j++)
		{
			double val = fabs(matrix[j][i]);
			if (val > maxValue)
			{
				maxN = j;
				maxValue = val;
			}
		}
		if (maxN > i)
		{
			double *row = matrix[maxN];
			matrix[maxN] = matrix[i];
			matrix[i] = row;
		}
		else if (maxValue == 0)
			return -1;

		double val = matrix[i][i];

		for (int j = i + 1; j < n; j++)
		{
			double k = matrix[j][i] / val;
			matrix[j][i] = 0;
			for (int c = i + 1; c < n; c++)
			{
				matrix[j][c] = matrix[j][c] - matrix[j][c] * k;
			}
		}
	}
	long double diagol = matrix[0][0];
	for (int i = 1; i < n; i++)
		diagol *= matrix[i][i];
	t2 = omp_get_wtime();
	double start = (t2 - t1);
	printf("Время последовательной версии - %.10f \n", start);
	if (diagol != 0)
		printf("Определитель - %.20f \n", diagol);
	else
		cout << "Определитель - 0 " << endl;	if (n < 10)
		printMatrixToConsole("Последовательная", matrix, n, diagol);
	return start;
}

double math_Openmp(double **matrix, int n)
{
	cout << "Расчет версии openmp" << endl;
	double t1 = omp_get_wtime(), t2;
	long double diagol;
#pragma omp parallel num_threads(20)
	{
#pragma omp for
		for (int i = 0; i < n - 1; i++)
		{
			int maxN = i;
			double maxValue = fabs(matrix[i][i]);
			//#pragma omp parallel for num_threads(20)
			for (int j = i + 1; j < n; j++)
			{
				double val = fabs(matrix[j][i]);
				if (val > maxValue)
				{
					maxN = j;
					maxValue = val;
				}
			}
			if (maxN > i)
			{
				double *row = matrix[maxN];
				matrix[maxN] = matrix[i];
				matrix[i] = row;
			}
			else if (maxValue == 0)
			{
				cout << "Нулевой элемент" << endl;
			}

			double value = matrix[i][i];
			//#pragma omp parallel for num_threads(20)
			for (int j = i + 1; j < n; j++)
			{
				double k = matrix[j][i] / value;
				matrix[j][i] = 0;
				for (int c = i + 1; c < n; c++)
				{
					matrix[j][c] = matrix[j][c] - matrix[j][c] * k;
				}
			}
		}
#pragma omp single
		{
			diagol = matrix[0][0];
		}
#pragma omp for
		for (int i = 1; i < n; i++)
			diagol *= matrix[i][i];

	}
	t2 = omp_get_wtime();
	double start = (t2 - t1);
	printf("Время omp версии - %.10f \n", start);
	if (diagol != 0)
		printf("Определитель - %.20f \n", diagol);
	else
		cout << "Определитель - 0 " << endl;
	if (n < 10)
		printMatrixToConsole("Последовательная", matrix, n, diagol);
	return start;
}

double math_Cuda(double **matrix, int n)
{
	CudaInfo *cuda = new CudaInfo();
	cuda->Info();
	cout << "Расчет версии cuda" << endl;
	cuda->OpredelitUpgrade(matrix, n);
	return 1;
}

double math_Gibrid(double **matrix, int n)
{
	cout << "Расчет гибридной версии" << endl;
	return 1;
}

double** getMatrix(double matrix[SIZE][SIZE])
{
	double **mat = new double*[SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		mat[i] = new double[SIZE];
		for (int j = 0; j < SIZE; j++)
			mat[i][j] = matrix[i][j];
	}
	return mat;
}

int main()
{
	srand(time(NULL));
	setlocale(LC_ALL, "Russian");
	long int n = SIZE;
	//cout << "Введите размер матрицы: ";
	//cin >> n;
	double mat[SIZE][SIZE], **matrix = new double*[n];
	cout << "Исходная матрица: " << endl;
	int num = 1;
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			mat[i][j] = 1 + rand() % 100;
			num++;
		}
	}
	matrix = getMatrix(mat);
	double res0 = math_Line(matrix, n);
	cout << endl;
	if (res0 != -1)
	{
#ifndef _OPENMP
		printf("This compiled code has no OpenMP support:( Check your compiler: if it supports OpenMP you must apply a correspondent compiler key.\n");
		goto exit;
#endif
		matrix = getMatrix(mat);
		double res1 = math_Openmp(matrix, n);
		printf("Производительность omp %f \n", fabs(res0 / res1));// / res1);
		matrix = getMatrix(mat);
		double res2 = math_Cuda(matrix, n);
		matrix = getMatrix(mat);
		double res3 = math_Gibrid(matrix, n);
	}
exit:
	system("pause");
}