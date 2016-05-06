#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream> //std
#include <sstream>
#include <fstream> //���� �����
#include <locale.h> //���� � ��� ������, �� ������ ��� ������� ������
#include <time.h>
#include <omp.h>

using namespace::std;

bool printMatrixToConsole(string title, float **matrix, int n, int m)
{
	cout << endl << "-----------------------------" << endl;
	cout << title << endl;
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new float[m];
		for (int j = 0; j < m; j++)
		{
			matrix[i][j] = (i + 1) * (i + 1) + j;
			matrix[i][j] += (((float)(i + 1) * (float)(j + 1)) / (float)rand());// *0.1f;
			cout << matrix[i][j] << "\t\t";
		}
		cout << endl;
	}
	cout << "-----------------------------" << endl << endl;
	return true;
}

int ������_���������������(float **matrix, int n, int m)
{
	int ������������ = 0;
	cout << "������ ���������������� ������";
	for (int i = 0; i < n-1; i++)
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
			������������++;
		}
		else if (maxValue == 0)
			return maxValue;

		float value = matrix[i][i];

		for (int j = i + 1; j < n; j++)
		{
			float k = matrix[j][i];
			for (int c = i+1; c < m; c++)
			{
				matrix[j][c] = matrix[j][c] - matrix[j][c] * k;
			}
		}
	}
	printMatrixToConsole("����������� ����������� �������", matrix, n, m);
	return 1;
}

int ������_openmp(float **matrix, int n, int m)
{
	cout << "������ ������ openmp" << endl;
	return 1;
}

int ������_cuda(float **matrix, int n, int m)
{
	cout << "������ ������ cuda" << endl;
	return 1;
}

int ������_��������(float **matrix, int n, int m)
{
	cout << "������ ��������� ������" << endl;
	return 1;
}

int main()
{
	srand(time(NULL));
	setlocale(LC_ALL, "Russian"); //���� � ��� ������, �� ������ ��� ������� ������
	long int n, m;
	cout << "������� ���������� ����� �������: ";
	cin >> n;
	/*cout << "������� ���������� �������� �������: ";
	cin >> m;*/
	m = n;
	float **matrix = new float*[n];
	cout << "�������� �������: "<< endl;
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new float[m];
		for (int j = 0; j < m; j++)
		{
			matrix[i][j] = (i+1) * (i+1) + j;
			matrix[i][j] += (((float)(i + 1) * (float)(j + 1)) / (float)rand());// *0.1f;
			cout << matrix[i][j] << "\t\t";
		}
		cout << endl;
	}
	int res0 = ������_���������������(matrix, n, m);
	int res1 = ������_��������(matrix, n, m);
	int res2 = ������_openmp(matrix, n, m);
	int res3 = ������_��������(matrix, n, m);

	system("pause");
}