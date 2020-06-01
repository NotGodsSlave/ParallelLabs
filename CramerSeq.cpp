#include<iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>
#include<random>

using namespace std;

const int N = 200; //test value

void init(double **&mx, double *&b, double max, double min, int m)
{
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


int main() 
{
	srand(time(NULL));
	int m = N;
	cout << "Input the size of the matrix\n";
	cin >> m;
	double** mx;
	mx = new double*[m];
	for (int i = 0; i<m; i++)
	    mx[i] = new double[m];
	double* b = new double[m];
	init(mx, b, 2, 0, m);
	/*cout << "matrix A:\n";
	printMatrix(mx, m);
	cout << "vector b: \n";
	for (int i = 0; i < m; ++i)
	{
		cout << b[i] << " ";
	}
	cout << '\n';*/
	double* x = new double[m];
	time_t START = clock();
	double d = getDeterminant(mx, m);
	//cout << "Determinant is " << d << '\n';
	for (int i = 0; i < m; ++i)
	{
		x[i] = getXi(mx, b, m, i, d);
		cout << "The value of " << i+1 << "-th unknown is " << x[i] << '\n';
	}
	time_t FINISH = clock();
	cout << "Time taken by sequential algorithm: " << double(FINISH - START) / CLOCKS_PER_SEC << "s\n";
} 