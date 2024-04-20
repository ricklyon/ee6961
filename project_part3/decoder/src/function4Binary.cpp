#include "function4Binary.h"
#include <math.h>

void ArrayMultiply4Binary(int *res, const int *a, const int* b, int n, int l)
/**
a: m x n
b: n x l
res: m x l
*/
{

	//static int add[4][4]={{0, 1, 2, 3},{1, 0, 3, 2},{2, 3, 0, 1},{3, 2, 1, 0}};
	//static int mul[4][4]={{0, 0, 0, 0},{0, 1, 2, 3},{0, 2, 3, 1},{0, 3, 1, 2}};
	int j, k;

	//for (i=0; i<m; i++)
	for (j = 0; j<l; j++)
	{
		int s = 0;
		for (k = 0; k<n; k++)
			s = s ^ (a[k] * b[k*l + j]);

		res[j] = s;
	} // for j, i

}

void ArrayMultiply4Binary(int *res, const int *a, const char* b, int n, int l)
/**
a: m x n
b: n x l
res: m x l
*/
{

	//static int add[4][4]={{0, 1, 2, 3},{1, 0, 3, 2},{2, 3, 0, 1},{3, 2, 1, 0}};
	//static int mul[4][4]={{0, 0, 0, 0},{0, 1, 2, 3},{0, 2, 3, 1},{0, 3, 1, 2}};
	int j, k;

	//for (i=0; i<m; i++)
	for (j = 0; j<l; j++)
	{
		int s = 0;
		for (k = 0; k<n; k++)
			s = s ^ (a[k] * (int)b[k*l + j]);

		res[j] = s;
	} // for j, i

}