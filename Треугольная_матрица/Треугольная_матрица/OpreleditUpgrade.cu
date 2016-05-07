#include <stdio.h>
#include <stdlib.h>
#include "CudaInfo.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__host__
void matrix_view(float(*array)[SIZE], char *q) {                  //����� ������� �� �����
	int i, j;
	for (i = 0; i<SIZE; i++) {
		for (j = 0; j<SIZE; j++)
			printf(q, array[i][j]);
		puts("");
	}
	puts("");
}

__host__
void matrix_rand(float(*array)[SIZE], double **matrix, int n) {               //���������� ��������� ������� �� �������
	int i, j;
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			array[i][j] = /*matrix[i][j]*/ 1 + rand() % 16;
}

__device__
void subrow(float(*array)[SIZE], int m, int n, float k) {         //��������� ����� �� ����������
	int g = blockIdx.y*blockDim.y + threadIdx.y;
	if (g<SIZE)
		array[m][g] -= k*array[n][g];
}

__global__
void determinant(float(*mtx)[SIZE]) {             //����: ���������� � ������������ ����
	int i = blockIdx.x*blockDim.x + threadIdx.x;
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

int CudaInfo::OpredelitUpgrade(float mtx_h[SIZE][SIZE], int n)//(double **matrix, int n)
{
	int i;
	float /*mtx_h[SIZE][SIZE],*/ (*mtx_d)[SIZE];
	long double det;
	cudaMalloc((void **)&mtx_d, sizeof(float)*SIZE*SIZE);          //��������� ����� �� ����������
	puts("�������� �������\n");
//	matrix_rand(mtx_h,matrix,n);
	time_t t1 = clock();
	cudaMemcpy(mtx_d, mtx_h, sizeof(float)*SIZE*SIZE,       //����������� ������� � ����� ����������
		cudaMemcpyHostToDevice);
	if (SIZE<10)
		matrix_view(mtx_h, "| %.0f ");
	dim3 threadsPerBlock(SIZE, SIZE, 1);
	dim3 numBlock(SIZE / threadsPerBlock.x, SIZE / threadsPerBlock.y, 1);
	determinant << <numBlock, threadsPerBlock >> >(mtx_d);                      //����� ����
	cudaThreadSynchronize();
	cudaMemcpy(mtx_h, mtx_d, sizeof(float)*SIZE*SIZE,       //����������� ������� �� ������ ����������
		cudaMemcpyDeviceToHost);
	puts("����������� ���\n");
	if (SIZE<10)
	matrix_view(mtx_h, "| %.8f ");
	det = 1;
	for (i = 0; i<SIZE; i++) {
		//printf("%.8f\n", mtx_h[i][i]);
		det *= mtx_h[i][i];
	}
	time_t t2 = clock();
	double timerun = ((t2 - t1) / CLOCKS_PER_SEC);
	printf("time work=%f\n",timerun);
	printf("det=%.0Lf\n ", det);
	cudaFree(mtx_d);
	return 0;
}