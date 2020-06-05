#include<iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>
#include<random>
#include<mpi.h>

using namespace std;

const int N = 200; //test value

int ProcNum; // Number of available processes
int ProcRank; // Rank of current process

void init(double** &mx, double* &b, double max, double min, int &m)
{
	if (ProcRank == 0)
	{
		cout << "Input the size of the matrix:\n";
		cin >> m;
	}
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	mx = new double*[m];
	for (int i = 0; i<m; i++)
	    mx[i] = new double[m];
	b = new double[m];
	for (int i = 0; i < m; ++i)
	{
		b[i] = (max - min) * ((double) rand()/(double)RAND_MAX) + min;
	}
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			mx[i][j] = (max - min) * ((double)rand()/(double)RAND_MAX) + min;
		}
	}
	/*MPI_Bcast(mx, m*m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(b, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
}

void printMatrix(double **&mx, int m)
{
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			cout << mx[i][j] << " ";
		}
		cout << '\n';
	}
}

void printVector(double *&v, int m)
{
	for (int i = 0; i < m; ++i)
	{
		cout << v[i] << " ";
	}
	cout << '\n';
}

void matrixToDiagonal(double **&mx, double **&p, int m)
{
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			p[i][j] = mx[i][j];
		}
	}

	for (int i = 0; i < m-1; ++i)
	{
		for (int k = i+1; k < m; ++k)
		{
			double fact = p[k][i]/p[i][i];
			for (int j = i; j < m; ++j)
			{
				p[k][j] -= fact*p[i][j];
			}
		}
	}
}

double getDeterminantDiagonal(double **&mx, int m)
{
	double d = 1;
	for (int i = 0; i < m; ++i)
	{
		d*=mx[i][i];
	}
	return d;
}

double getDeterminant(double **&mx, int m) {
	double d;
	double **p;
	p = new double*[m];
	for (int i = 0; i<m; i++)
	    p[i] = new double[m]; 
	d = 0;
	if (m < 1) cout << "Can't compute detrminant\n";
	if (m == 1) {
	    d = mx[0][0];
	    return d;
	}
	if (m == 2) {
	    d = mx[0][0] * mx[1][1] - (mx[1][0] * mx[0][1]);
	    return d;
	}
	if (m > 2) {
	    matrixToDiagonal(mx,p,m);
	    //printMatrix(p,m);
	    d = getDeterminantDiagonal(p, m);
	    return d;
	}
	return d;
}

void replaceColumn(double **&mx, double **&p, double *&b, int m, int i)
{
	for (int j = 0; j < m; ++j)
	{
		for (int k = 0; k < m; ++k)
		{
			if (k == i) 
			{
				p[j][k] = b[j];
			}
			else p[j][k] = mx[j][k];
		}
	}
}

double getXi(double **&mx, double *&b, int m, int i, double d)
{
	double **p;
	p = new double*[m];
	for (int j = 0; j < m; j++)
	    p[j] = new double[m];
	replaceColumn(mx, p, b, m, i);
	return getDeterminant(p, m)/d;
}

void ParallelResultCalculation(double **&mx, double *&b, int*&mn, int*& ms, double *&x, int m, double d)
{
	for (int i = ms[ProcRank]; i < mn[ProcRank] + ms[ProcRank]; ++i)
	{
		x[i-ms[ProcRank]] = getXi (mx, b, m, i, d);
	}
}

void TestPartialResults(double *&x, int*&mn, int*& ms)
{
	for (int j = 0; j < ProcNum; ++j)
	{
		if (ProcRank == j)
		{
			cout << "In process " << j << endl;
			for (int i = ms[ProcRank]; i < mn[ProcRank] + ms[ProcRank]; ++i)
			{
				cout << x[i-ms[ProcRank]] << " ";
			}
		}
		cout << endl;
	}
}


int main(int argc, char*argv[]) 
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	srand(time(NULL));
	int m = N;
	double** mx;
	double* b;
	double** pProcMatrix; // Stripe of the matrix on current process
	double* pProcResult; // Block of result vector on current process
	
	init(mx, b, 2, 0, m);
	MPI_Barrier(MPI_COMM_WORLD);

	int* matrixNumbers;
	matrixNumbers = new int[ProcNum];
	int* matrixStarts;
	matrixStarts = new int[ProcNum];
	int temp = m%ProcNum;
	int sum = 0;
	for (int i = 0; i < ProcNum; ++i)
	{
		matrixNumbers[i] = m/ProcNum;
		if (temp > 0) 
		{
			matrixNumbers[i]++;
			temp--;
		}
		matrixStarts[i] = sum;
		sum+=matrixNumbers[i];
		if (ProcRank == i)
		{
			pProcResult = new double[matrixNumbers[i]];
		}
	}
	
		
	/*cout << "matrix A:\n";
	printMatrix(mx, m);
	cout << "vector b: \n";
	for (int i = 0; i < m; ++i)
	{
		cout << b[i] << " ";
	}
	cout << '\n';*/
		

	double* x = new double[m];
	MPI_Scatter(x, matrixNumbers[ProcRank], MPI_DOUBLE, pProcResult, matrixNumbers[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	time_t START, FINISH;
	double d;
	if (ProcRank == 0){
		START = clock();
		d = getDeterminant(mx, m);
	}
	MPI_Bcast(&d, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "Determinant is " << d << '\n';
	ParallelResultCalculation(mx, b, matrixNumbers, matrixStarts, pProcResult, m, d);
	MPI_Barrier(MPI_COMM_WORLD);
	//TestPartialResults(pProcResult, matrixNumbers, matrixStarts);
	MPI_Allgather(pProcResult, matrixNumbers[ProcRank], MPI_DOUBLE, x, matrixNumbers[ProcRank], MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		for (int i = 0; i < m; ++i)
		{
			cout << "Value of " << i+1 << "-th unknown: " << x[i] << endl;
		}
		FINISH = clock();
		cout << "Time taken by MPI algorithm: " << double(FINISH - START) / CLOCKS_PER_SEC << "s\n";
	}
	MPI_Finalize();
} 