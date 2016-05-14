#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream> //std
#include <sstream>
#include <fstream> //ввод вывод
#include <locale.h> //Если у вас линукс, то можете эту строчку убрать
#include <time.h>
#include <omp.h>

using namespace std;

int main()
{
	double Sum = 0;
#pragma omp parallel
	{
		printf("thad %d \n", omp_get_thread_num());
		//if (omp_get_thread_num() == 0)
		//{
		for (int i = 0; i < 1; i++)
		{
#pragma omp for
			for (int j = 0; j < 16; j++)
			{
				printf("thadddddd %d \n", omp_get_thread_num());
				//cout << "olololo " + omp_get_thread_num() << endl;
				double Sum1 = 0.5 / 0.4 * 2 + 2 + 0.5 * 0.4 * Sum;
				if (Sum1 / 2 == 0)
				{
					Sum += Sum1;
				}
				else if (Sum1 / 2 != 0)
				{
					Sum -= Sum1;
				}
				if (j == 1000000)
				{
					//						cout << omp_get_thread_num << endl;
				}
			}
		}
	}
	system("pause");
}