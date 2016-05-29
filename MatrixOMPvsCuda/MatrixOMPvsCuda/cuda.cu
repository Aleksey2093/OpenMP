#include "Main.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cuda_runtime_api.h>

__device__ int maxNcuda;
__device__ float maxValueCuda;

__global__ void determ()
{

}

__global__ void swaprow(float *Matrix, int row1, int row2, int n)
{
	int id = blockIdx.x + row1*n;
	int id2 = blockIdx.x + row2*n;
	__shared__ float tmp[SIZE];
	tmp[blockIdx.x] = Matrix[id];
	Matrix[id] = Matrix[id2];
	Matrix[id2] = tmp[blockIdx.x];
}

__global__ void proizrow(float *Matrix, int i, int j, int n)
{
	int id = i*n + threadIdx.x; //первая строка
	int id2 = blockIdx.x*blockDim.x /*+ j*n*/ + threadIdx.x; //последующие
	if (id2 < n*n)
	{
		float jim = Matrix[j*n + i];
		float val = Matrix[i*n + i];
		__syncthreads();
		if (id2 == j*n + i)
		{
			Matrix[id2] = 0;
		}
		__syncthreads();
		if (id2 > j * n + i)
		{
			float t = Matrix[id2] - Matrix[id] * (jim / val);
			Matrix[id2] = t;
		}
	}
}

__global__ void seamax(float *Matrix, int row1, int n)
{

}

double** Main::getMatrixFromCuda(double** Matrix, Main *cu)
{
	int n = SIZE;
	double end, start;
	float *dev_Matrix, *Matrixline = new float[n*n];
	int g = (n*n + 127) / 128, th = 9;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Matrixline[i*n + j] = Matrix[i][j];
		}
	}
	cudaMalloc((void**)&dev_Matrix, n*n*sizeof(float));
	cudaMemcpy(dev_Matrix, Matrixline, n*n*sizeof(float), cudaMemcpyHostToDevice);

	int maxN; float maxValue;
	
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
	start = omp_get_wtime();
	for (int i = 0; i < n - 1; i++)
	{
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				{
					maxN = i;
				}
#pragma omp section
				{
					maxValue = fabs(Matrixline[i*n + i]);
				}
			}
#pragma omp for
			for (int j = i + 1; j < n; j++) //поиск максимального элемента
			{
				float tmp = fabs(Matrixline[j*n + i]);
				if (tmp > maxValue)
				{
					maxN = j;
					maxValue = tmp;
				}
			}
#pragma omp single
			{
				if (maxN > i)
				{
					swaprow << <n, 1 >> >(dev_Matrix, i, maxN, n);
				}
				else if (maxValue == 0)
				{
					printf("null in matrix\n");
				}
			}
#pragma omp for
			for (int j = i + 1; j < n; j++)
			{
				proizrow << <n, th >> >(dev_Matrix, i, j, n);
			}
#pragma omp single
			{
				cudaMemcpy(Matrixline, dev_Matrix, n*n*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
	}
	end = omp_get_wtime();

	cudaMemcpy(Matrixline, dev_Matrix, n*n*sizeof(float), cudaMemcpyDeviceToHost);

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