#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream> //std
#include <sstream>
#include <fstream> //in out
#include <locale.h> 
#include <time.h>
#include <omp.h>
#include <cuda_runtime_api.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define SIZE 7

class CudaInfo
{
public:
	CudaInfo();
	~CudaInfo();
	void Info();
	int kernel();
	int Opredelit();
	int OpredelitUpgrade(float matrix[SIZE][SIZE], int n);
private:

};