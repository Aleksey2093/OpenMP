#include "Main.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cuda_runtime_api.h>

__global__ void determ()
{

}

__global__ void swaprow(float *Matrix, int row1, int row2, int n)
{
	float t;
	int id = blockIdx.x + row1*n;
	int id2 = blockIdx.x + row2*n;
	Matrix[id] = Matrix[id2];
	Matrix[id2] = t;
}


__global__ void proizrow(float *Matrix, int i, int j, int n, float val)
{
	/*__shared__ float koe[SIZE];
	int id = blockIdx.x + i*n;
	float row1[SIZE];
	if (blockIdx.x < SIZE) //первая строчка сетки нужная нам строка идем по ней и собираем нужную информацию
	{
		row1[blockIdx.x] = Matrix[id];
		//printf("%.1f\t", row1[blockIdx.x]);
	}
	if (blockIdx.x >= SIZE)
	{
		for (int p = 0; p < SIZE; p++)
		{
			Matrix[blockIdx.x + j*n] = Matrix[blockIdx.x + j*n] + row1[p];
		}
	}*/
	int id = blockIdx.x + i*n;
	int id2 = blockIdx.x + j*n;
	float t = Matrix[id2] - Matrix[id];// *Matrix[i*SIZE + j] / val;
	if (id2 == i*SIZE + j)
	{
		Matrix[id2] = 0;
	}
	printf("(%.1f)", Matrix[i*SIZE + j]);
	Matrix[id2] = t;
}

/*for (int j = i + 1; j < SIZE; j++)
{
float k = Matrixline[j*SIZE + i] / val;
Matrixline[j*SIZE + i] = 0;
for (int c = i + 1; c < SIZE; c++)
Matrixline[j*SIZE + c] = Matrixline[j*SIZE + c] - Matrixline[i*SIZE + c] * k;
}*/

__global__ void seamax(float *Matrix, int rowstart, int colstart, int n, int g, int th)
{
	/*		
	for (int j = i + 1; j < SIZE; j++) //поиск максимального элемента
		{
			double tmp = fabs(Matrix[j][i]);
			if (tmp > maxValue)
			{
				maxN = j;
				maxValue = tmp;
			}
		}*/
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	__shared__ float maxblock[(SIZE + 127) / 128]; //массив максимумов по блокам
	maxblock[threadIdx.x + blockDim.x * blockIdx.x] = Matrix[rowstart*n + colstart];
	for (int i = 0; i < th; i++)
	{

	}
}

double** Main::getMatrixFromCuda(double** Matrix)
{
	int n = SIZE*SIZE;
	float *dev_Matrix, *Matrixline = new float[n];
	int g = (SIZE + 127) / 128, th = 128;
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			Matrixline[i*SIZE + j] = Matrix[i][j];
		}
	}
	cudaMalloc((void**)&dev_Matrix, n*sizeof(float));
	cudaMemcpy(dev_Matrix, Matrixline, n*sizeof(float), cudaMemcpyHostToDevice);

	int maxN; float maxValue;
	printf("\n");
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			printf("%.2f\t",Matrixline[i*SIZE + j]);
		}
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < SIZE-1; i++)
	{
		//seamax << <g, th >> >(dev_Matrix, i, i + 1, SIZE, g, th);
//#pragma omp parallel
		{
			maxN = i;
			maxValue = fabs(Matrixline[i*SIZE+i]);
			for (int j = i + 1; j < SIZE; j++) //поиск максимального элемента
			{
				float tmp = fabs(Matrixline[j*SIZE+i]);
				if (tmp > maxValue)
				{
					maxN = j;
					maxValue = tmp;
				}
			}
			if (maxN > i)
			{
				cudaMemcpy(dev_Matrix, Matrixline, n*sizeof(float), cudaMemcpyHostToDevice);
				swaprow<<<SIZE,1>>>(dev_Matrix, i, maxN, SIZE);
				cudaMemcpy(Matrixline, dev_Matrix, n*sizeof(float), cudaMemcpyDeviceToHost);
			}
			else if (maxValue == 0)
			{
				printf("null in matrix\n");
			}
			float val = Matrixline[i *SIZE + i];

			//cudaMemcpy(dev_Matrix, Matrixline, n*sizeof(float), cudaMemcpyHostToDevice);
			for (int j = i + 1; j < SIZE; j++)
			{
				//float k = Matrixline[j * SIZE + i] / val;
				cudaMemcpy(dev_Matrix, Matrixline, n*sizeof(float), cudaMemcpyHostToDevice);
				proizrow << <SIZE, 1 >> >(dev_Matrix, i, j, SIZE, val);
				cudaMemcpy(Matrixline, dev_Matrix, n*sizeof(float), cudaMemcpyDeviceToHost);
			}
			//cudaMemcpy(Matrixline, dev_Matrix, n*sizeof(float), cudaMemcpyDeviceToHost);

			for (int i = 0; i < SIZE; i++)
			{
				for (int j = 0; j < SIZE; j++)
				{
					printf("%.2f\t", Matrixline[i*SIZE + j]);
				}
				printf("\n");
			}
			printf("\n");
			//__global__ void proizrow(float *Matrix, int i, int j, int n, float val)
			/*for (int j = i + 1; j < SIZE; j++)
			{
				float k = Matrixline[j*SIZE + i] / val;
				Matrixline[j*SIZE + i] = 0;
				for (int c = i + 1; c < SIZE; c++)
					Matrixline[j*SIZE + c] = Matrixline[j*SIZE + c] - Matrixline[i*SIZE + c] * k;
			}*/
		}
	}

	/*cudaMalloc((void**)&dev_Matrix, n*sizeof(float));
	cudaMemcpy(dev_Matrix, Matrix, n*sizeof(float), cudaMemcpyHostToDevice);
	dim3 grid(n, 1);
	add << <SIZE, 1 >> >(dev_Matrix);
	cudaMemcpy(Matrix, dev_Matrix, n*sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(dev_Matrix);*/

	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			printf("%.2f\t", Matrixline[i*SIZE + j]);
		}
		printf("\n");
	}

	return Matrix;
}