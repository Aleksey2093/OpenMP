#pragma once
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream> //std
#include <sstream>
#include <fstream> //in out
#include <locale.h> 
#include <time.h>
#include <omp.h>

#define SIZE 9

class Main
{
public:
	double** getMatrixFromCuda(double** Matrix3, Main *cu);
public:
	double timelinever;
	double timeompver;
	double timecudaver;
































































































};

