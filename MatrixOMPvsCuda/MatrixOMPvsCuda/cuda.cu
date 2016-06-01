#include "Main.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>

__device__ void swaprow(double *Matrix, int row1, int row2, int n)
{
	int r = blockIdx.x*blockDim.x;
	int idx = r + threadIdx.x;
	if (idx >= row2*n && idx < (row2+1)*n)
	{
		int id = idx - row2*n + row1*n;
		int id2 = idx;//threadIdx.x + row2*n;
		//printf("id: %d, id2: %d",id,id2);
		double tmp;
		tmp = Matrix[id];
		Matrix[id] = Matrix[id2];
		Matrix[id2] = tmp;
	}
}

__device__ int maxN;
__device__ int *copy;
__device__ double maxValue;
__global__ void seamax(double *Matrix, int i, int n)
{
	if (blockIdx.x*blockDim.x + threadIdx.x == 0)
	{
		maxN = i;
		maxValue = fabs(Matrix[i*n + i]);
		for (int j = i + 1; j < n; j++) //поиск максимального элемента
		{
			double tmp = fabs(Matrix[j*n + i]);
			if (tmp > maxValue)
			{
				maxN = j;
				maxValue = tmp;
			}
		}
		copy = &maxN;
	}
	__syncthreads();
	if (maxN > i)
	{
		swaprow(Matrix, i, maxN, n);
	}
	else if (maxValue == 0)
	{
		if (blockIdx.x*blockDim.x + threadIdx.x == 0)
			printf("nulllll please stop \n");
	}
}

__global__ void proizrow(double *Matrix, int i, int j, int n)
{
	int id2 = blockIdx.x*blockDim.x + threadIdx.x; //последующие
	if (id2 < n*n)
	{
		int id = i*n + threadIdx.x; //первая строка
		int jni = j*n + i;
		double jim = Matrix[jni] / Matrix[i*n + i];
		__syncthreads();
		if (id2 == jni)
		{
			Matrix[id2] = 0;
		}
		else if (id2 > jni && id2 < (-1)*jni - i + n)
		{
			Matrix[id2] = Matrix[id2] - Matrix[id] * (jim);
		}
	}
}

double** Main::getMatrixFromCuda(double** Matrix, Main *cu)
{
	int n = SIZE;
	double end, start;
	double *dev_Matrix, *Matrixline = new double[n*n];
	int g = (n*n + 127) / 128, th = 128;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Matrixline[i*n + j] = Matrix[i][j];
		}
	}
	cudaMalloc((void**)&dev_Matrix, n*n*sizeof(double));
	cudaMemcpy(dev_Matrix, Matrixline, n*n*sizeof(double), cudaMemcpyHostToDevice);
	
	if (n < 10)
	{
		printf("\n");
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				printf("%.2f\t", Matrixline[i*n + j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	dim3 grid(n, (n+127) / 128);
	start = omp_get_wtime();
	for (int i = 0; i < n-1; i++)
	{
		seamax << <n, th >> >(dev_Matrix, i, n);
#pragma omp parallel for
		for (int j = i + 1; j < n; j++)
		{
			proizrow << <n, th >> >(dev_Matrix, i, j, n);
		}
	}
	end = omp_get_wtime();
	cudaMemcpy(Matrixline, dev_Matrix, n*n*sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(dev_Matrix);
	end = end - start;
	cu->timecudaver = end;

	if (n < 10)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				printf("%.2f\t", Matrixline[i*n + j]);
			}
			printf("\n");
		}
	}

	return Matrix;
}