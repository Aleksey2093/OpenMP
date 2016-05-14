#include "CudaInfo.cuh";

using namespace::std;

bool printMatrixToConsole(string title, float **matrix, int n, float diagol)
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

float math_Line(float **matrix, int n)
{
	cout << "Расчет последовательной версии" << endl;
	float t1 = omp_get_wtime(), t2;
	for (int i = 0; i < n - 1; i++)
	{
		int maxN = i;
		float maxValue = fabs(matrix[i][i]);
		for (int j = i + 1; j < n; j++)
		{
			float val = fabs(matrix[j][i]);
			if (val > maxValue)
			{
				maxN = j;
				maxValue = val;
			}
		}
		if (maxN > i)
		{
			float *row = matrix[maxN];
			matrix[maxN] = matrix[i];
			matrix[i] = row;
		}
		else if (maxValue == 0)
		{
			cout << "null row - " << i << endl;
		}

		float val = matrix[i][i];

		for (int j = i + 1; j < n; j++)
		{
			float k = matrix[j][i] / val;
			matrix[j][i] = 0;
			for (int c = i+1; c < n; c++)
			{
				matrix[j][c] = matrix[j][c] - matrix[j][c] * k;
			}
		}
	}
	long float diagol = matrix[0][0];
	for (int i = 1; i < n; i++)
		diagol *= matrix[i][i];
	t2 = omp_get_wtime();
	float start = (t2 - t1);
	printf("Время последовательной версии - %.10f \n", start);
	if (diagol != 0)
		printf("Определитель - %.20f \n", diagol);
	else
		cout << "Определитель - 0 " << endl;	if (n < 10)
		printMatrixToConsole("Последовательная", matrix, n, diagol);
	return start;
}

float math_Openmp(float **matrix, int n)
{
	cout << "Расчет версии openmp" << endl;
	float t1 = omp_get_wtime(), t2;
	long float diagol;
#pragma omp parallel num_threads(20)
	{
#pragma omp for
		for (int i = 0; i < n - 1; i++)
		{
			int maxN = i;
			float maxValue = fabs(matrix[i][i]);
			//#pragma omp parallel for num_threads(20)
			for (int j = i + 1; j < n; j++)
			{
				float val = fabs(matrix[j][i]);
				if (val > maxValue)
				{
					maxN = j;
					maxValue = val;
				}
			}
			if (maxN > i)
			{
				float *row = matrix[maxN];
				matrix[maxN] = matrix[i];
				matrix[i] = row;
			}
			else if (maxValue == 0)
			{
				cout << "Нулевой элемент" << endl;
			}

			float value = matrix[i][i];
			//#pragma omp parallel for num_threads(20)
			for (int j = i + 1; j < n; j++)
			{
				float k = matrix[j][i] / value;
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
	float start = (t2 - t1);
	printf("Время omp версии - %.10f \n", start);
	if (diagol != 0)
		printf("Определитель - %.20f \n", diagol);
	else
		cout << "Определитель - 0 " << endl;
	if (n < 10)
		printMatrixToConsole("OMP матрица", matrix, n, diagol);
	return start;
}

float math_Cuda(float matrix[SIZE][SIZE], int n)
{
	CudaInfo *cuda = new CudaInfo();
	cuda->Info();
	cout << "Расчет версии cuda" << endl;
	cuda->OpredelitUpgrade(matrix, n);
	return 1;
}

float math_Gibrid(float matrix[SIZE][SIZE], int n)
{
	cout << "Расчет гибридной версии" << endl;
	return 1;
}

float** getMatrix(float matrix[SIZE][SIZE])
{
	float **mat = new float*[SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		mat[i] = new float[SIZE];
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
	float mat[SIZE][SIZE], **matrix = new float*[n];
	cout << "Исходная матрица: " << endl;
	int num = 1;
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new float[n];
		for (int j = 0; j < n; j++)
		{
			mat[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			mat[i][j] = mat[i][j] + (1 + rand() % 30);
			// (float)(1 + rand() % 1000) + (float)1 / (float)(1 + rand() % 1000);
			num++;
		}
	}
	matrix = getMatrix(mat);
	if (n < 10)
		printMatrixToConsole("Исходная матрица", matrix, n, 0);
	float res0 = math_Line(matrix, n);
	cout << endl;
	if (res0 != -5)
	{
#ifndef _OPENMP
		printf("This compiled code has no OpenMP support:( Check your compiler: if it supports OpenMP you must apply a correspondent compiler key.\n");
		goto exit;
#endif
		matrix = getMatrix(mat);
		float res1 = math_Openmp(matrix, n);
		printf("Производительность omp %f \n", fabs(res0 / res1));// / res1);
		float res2 = math_Cuda(mat, n);
		float res3 = math_Gibrid(mat, n);
	}
	else
	{
		cout << "Обнаружены нули" << endl;
	}
exit:
	system("pause");
}