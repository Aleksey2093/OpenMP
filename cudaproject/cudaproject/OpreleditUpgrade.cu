#include <stdio.h>
#include <stdlib.h>
#include "CudaInfo.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define SIZE 3

__host__
void matrix_view(float(*array)[SIZE], char *q) {                  //вывод матрицы на экран
	int i, j;
	for (i = 0; i<SIZE; i++) {
		for (j = 0; j<SIZE; j++)
			printf(q, array[i][j]);
		puts("");
	}
	puts("");
}

__host__
void matrix_rand(float(*array)[SIZE]) {               //заполнение случайным образом на девайсе
	int i, j;
	for (i = 0; i<SIZE; i++)
		for (j = 0; j<SIZE; j++)
			array[i][j] = 1 + rand() % 16;
}

__device__
void subrow(float(*array)[SIZE], int m, int n, float k) {         //вычитание строк на устройстве
	int g = blockIdx.y*blockDim.y + threadIdx.y;
	if (g<SIZE)
		array[m][g] -= k*array[n][g];
}

__global__
void determinant(float(*mtx)[SIZE]) {             //Ядро: приведение к треугольному виду
	int i = i = blockIdx.x*blockDim.x + threadIdx.x;
	int j;
	float coeff;
	for (j = 0; j<SIZE - 1; j++) {
		if (!mtx[j][j]) subrow(mtx, j, SIZE - 1, 3);
		if (i >= j && i<SIZE - 1) {
			coeff = mtx[i + 1][j] / mtx[j][j];
			subrow(mtx, i + 1, j, coeff);
		}
		__syncthreads();
	}
}

int CudaInfo::OpredelitUpgrade(void)
{
	int i;
	float mtx_h[SIZE][SIZE], (*mtx_d)[SIZE];
	long double det;
	cudaMalloc((void **)&mtx_d, sizeof(float)*SIZE*SIZE);          //выделение памяти на устройстве
	puts("исходная матрица\n");
	matrix_rand(mtx_h);
	cudaMemcpy(mtx_d, mtx_h, sizeof(float)*SIZE*SIZE,       //копирование массива в память видеокарты
		cudaMemcpyHostToDevice);
	matrix_view(mtx_h, "| %.0f ");
	dim3 threadsPerBlock(16, 16);
	dim3 numBlock(SIZE / threadsPerBlock.x, SIZE / threadsPerBlock.y);
	determinant << <numBlock, threadsPerBlock >> >(mtx_d);                      //вызов ядра
	cudaThreadSynchronize();
	cudaMemcpy(mtx_h, mtx_d, sizeof(float)*SIZE*SIZE,       //копирование массива из памяти видеокарты
		cudaMemcpyDeviceToHost);
	puts("трeугольный вид\n");
	matrix_view(mtx_h, "| %.8f ");
	det = 1;
	for (i = 0; i<SIZE; i++) {
		printf("%.8f\n", mtx_h[i][i]);
		det *= mtx_h[i][i];
	}
	printf("det=%.0Lf\n ", det);
	cudaFree(mtx_d);
	getchar();
	return 0;
}