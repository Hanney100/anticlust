/*
 * copied from:
 * Xiao Yang et al. “A three-phase search approach with dynamic population size for solving
 the maximally diverse grouping problem”. In: European Journal of Operational Research
 302.3 (2022). [SOURCE-CODE: https://raw.githubusercontent.com/toyamaailab/toyamaailab.github.io/main/resource/TPSDP_Code.zip], pp. 925–953. ISSN: 0377-2217. DOI: https://doi.org/
 10.1016/j.ejor.2022.02.003. 
 * */
 
#define WIN32_LEAN_AND_MEAN
#define _CRT_SECURE_NO_WARNINGS
#include <Rcpp.h>
#include <cstdlib>
#include <cstdlib>
#include "stdio.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <algorithm>
 
using namespace std;
using namespace Rcpp;


typedef struct Solution {
	int* p;  		// cluster membership of each element
	int* SizeG;  	// size of each cluster
	double cost;	// global cost funtion value of solution
}Solution;

typedef struct Neighborhood {
	int  type;
	int  v;
	int  g;
	int  x;
	int  y;
}Neighborhood;

char* File_Name;
char* Output_File_Name;
char* Solution_File;
int popSize;
int N, K;  // node number and group number
double f, f_best;
Solution CS, NS, GS, OS;
Solution* Pop;
Solution* Offs;
Neighborhood* Neighbors;
double total_time, starting_time, Time_limit;
int* p; // patition array for each vertex
int* bestp; // patition array for each vertex 
int* SizeG;
double** Delta_Matrix;  // incremental matrix 
double** Delta_Matrix_p1;
double** Delta_Matrix_p2;
double** Delta_Matrix1;
double** Delta;
double** Delta1;
double* gDiv;
double* gDiv_p1;
double* gDiv_p2;
int* SelectEle;
int* SelectEleTemp;
int* SelGroup;
int* p1;
int* p2;

// for crossover
int* vEle;
int* gEle;
int* LBGroup;
int* UBGroup;
int* BigThanLB;
int* ub;


double** D;   // distance matrix between elements
double** DT;
int* LB; // Lower bound of the number of elements in the group i 
int* UB; // Upper bound of the number of elements in the group i
double theta, theta_max, theta_min; 
int beta_min; 
int LMAX;

double** AvgCon;
int* Rd, * UnderLB;

void AssignMemery()
{
	/*  Allocates memory dynamically for various arrays and matrices necessary 
	for the algorithm's execution. This includes structures for population management, 
	distance matrices, diversity measures, and neighborhood exploration.
	*/
	int i, j;

	p = new int[N];
	bestp = new int[N];
	SizeG = new int[K];

	Pop = new Solution[popSize];
	Offs = new Solution[popSize];

	Delta_Matrix = new double* [N];
	for (i = 0; i < N; i++) Delta_Matrix[i] = new double[K];
	Delta_Matrix_p1 = new double* [N];
	for (i = 0; i < N; i++) Delta_Matrix_p1[i] = new double[K];
	Delta_Matrix_p2 = new double* [N];
	for (i = 0; i < N; i++) Delta_Matrix_p2[i] = new double[K];
	Delta_Matrix1 = new double* [N];
	for (i = 0; i < N; i++) Delta_Matrix1[i] = new double[K];
	gDiv = new double[K];
	gDiv_p1 = new double[K];
	gDiv_p2 = new double[K];

	Delta = new double* [N];
	for (i = 0; i < N; i++) Delta[i] = new double[K];

	Delta1 = new double* [N];
	for (i = 0; i < N; i++) Delta1[i] = new double[K];

	for (i = 0; i < popSize; i++) {
		Pop[i].p = new int[N];
		Offs[i].p = new int[N];
	}

	for (i = 0; i < popSize; i++) {
		Pop[i].SizeG = new int[K];
		Offs[i].SizeG = new int[K];
	}

	CS.p = new int[N];
	NS.p = new int[N];
	GS.p = new int[N];
	OS.p = new int[N];

	CS.SizeG = new int[K];
	NS.SizeG = new int[K];
	GS.SizeG = new int[K];
	OS.SizeG = new int[K];

	Neighbors = new Neighborhood[N * (N - 1) / 2 + N * K];

	AvgCon = new double* [K];
	for (i = 0; i < K; i++) AvgCon[i] = new double[K];
	Rd = new int[K];
	for (i = 0; i < K; i++) Rd[i] = 0;
	UnderLB = new int[K];

	ub = new int[K];
	LBGroup = new int[K];
	UBGroup = new int[K];
	BigThanLB = new int[K];
	vEle = new int[N];
	gEle = new int[K];
	SelectEle = new int[N];
	SelGroup = new int[K];
	SelectEleTemp = new int[N];
	p1 = new int[N];
	p2 = new int[N];

}

void ReleaseMemery()
{
	/* responsible for reading the input file, 
	   initializing matrices, and setting constraints on group sizes.
	*/
	int i;

	delete[] p; p = NULL;
	delete[] bestp; bestp = NULL;
	delete[] SizeG; SizeG = NULL;

	delete[] CS.p; CS.p = NULL;
	delete[] CS.SizeG; CS.SizeG = NULL;
	delete[] GS.p; GS.p = NULL;
	delete[] GS.SizeG; GS.SizeG = NULL;
	delete[] NS.p; NS.p = NULL;
	delete[] NS.SizeG; NS.SizeG = NULL;
	delete[] OS.p; OS.p = NULL;
	delete[] OS.SizeG; OS.SizeG = NULL;

	delete[] LB; LB = NULL;
	delete[] UB; UB = NULL;
	delete[] Neighbors; Neighbors = NULL;

	for (i = 0; i < N; i++)
	{
		delete[] Delta_Matrix[i]; Delta_Matrix[i] = NULL;
		delete[] Delta_Matrix1[i]; Delta_Matrix1[i] = NULL;
		delete[] Delta_Matrix_p1[i]; Delta_Matrix_p1[i] = NULL;
		delete[] Delta_Matrix_p2[i]; Delta_Matrix_p2[i] = NULL;
		delete[] Delta[i]; Delta[i] = NULL;
		delete[] Delta1[i]; Delta1[i] = NULL;
		delete[] D[i]; D[i] = NULL;
		delete[] DT[i]; DT[i] = NULL;
	}
	for (i = 0; i < K; i++) {
		delete[] AvgCon[i]; AvgCon[i] = NULL;
	}

	for (i = 0; i < popSize; i++) {
		delete[] Pop[i].p; Pop[i].p = NULL;
		delete[] Offs[i].p; Offs[i].p = NULL;
		delete[] Pop[i].SizeG; Pop[i].SizeG = NULL;
		delete[] Offs[i].SizeG; Offs[i].SizeG = NULL;
	}
	delete[] Rd; Rd = NULL;
	delete[] UnderLB; UnderLB = NULL;
	delete[] SelectEle; SelectEle = NULL;
	delete[] SelGroup; SelGroup = NULL;
	delete[] SelectEleTemp; SelectEleTemp = NULL;
	delete[] p1; p1 = NULL;
	delete[] p2; p2 = NULL;
	delete[] gDiv; gDiv = NULL;
	delete[] gDiv_p1; gDiv_p1 = NULL;
	delete[] gDiv_p2; gDiv_p2 = NULL;

	delete[] ub; ub = NULL;
	delete[] LBGroup; LBGroup = NULL;
	delete[] UBGroup; UBGroup = NULL;
	delete[] BigThanLB; BigThanLB = NULL;
	delete[] vEle; vEle = NULL;
	delete[] gEle; gEle = NULL;
}

int Proof(Solution& S)
{
	/* verifies the feasibility of a solution and calculates its cost
	*/
	int i, j;
	double ff;
	int flag;
	ff = 0.0;
	for (i = 0; i < N; i++)
		for (j = i + 1; j < N; j++)
		{
			if (S.p[i] == S.p[j])
			{
				ff += D[i][j];
			}
		}
	S.cost = ff;
	for (i = 0; i < K; i++) S.SizeG[i] = 0;
	for (i = 0; i < N; i++) S.SizeG[S.p[i]]++;
	flag = 1;
	for (i = 0; i < K; i++)
		if (S.SizeG[i] < LB[i] || S.SizeG[i]> UB[i]) { flag = 0; break; }
	
	return flag;
}

void Outputing(Solution& S, char* filename)
{
	/* writes the solution details to a file if it is feasible */
	int i; int r;
	FILE* fp;
	char buff[80];
	r = rand() % 1000;
	if (Proof(S) == 0) return;
	sprintf(buff, "%s", filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "N = %d  G = %d  f = %lf\n", N, K, S.cost);
	for (i = 0; i < K; i++)
		fprintf(fp, "%5.2d   %5.2d   %5.2d \n", LB[i], UB[i], S.SizeG[i]);
	printf("\n");
	for (i = 0; i < N; i++)
		fprintf(fp, "%5.3d   %5.2d\n", i, S.p[i]);
	fclose(fp);
}

void Out_results(double best, double ave, double worst, char* filename, char instance[])
{
	/* writes the best, average, and worst results to a file */
	int i;
	FILE* fp;
	char buff[80];
	sprintf(buff, "%s", filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s   %lf   %lf   %lf\n", instance, best, ave, worst);
	fclose(fp);
}

void RandomInitiaSol(int p[], int SizeG[])
{
	/* Algorithm 2: initializes a random solution that respects the group size constraints */
	int i, j;
	int p1;
	int count;
	int tot_number;
	int sum = 0;
	int* Flag;
	int* SizeGroup;
	SizeGroup = new int[K];
	for (i = 0; i < K; i++) SizeGroup[i] = 0;
	Flag = new int[N];
	for (i = 0; i < N; i++) Flag[i] = 0;
	for (i = 0; i < K; i++) sum += LB[i];
	tot_number = 0;
	while (tot_number < sum)
	{
		p1 = rand() % N;
		if (Flag[p1] == 0)
		{
			count = 0;
			while (count < K)
			{
				if (SizeGroup[count] < LB[count])
				{
					p[p1] = count;
					Flag[p1] = 1;
					SizeGroup[count]++;
					tot_number++;
					break;
				}
				else count++;
			}
		}

	}

	tot_number = 0;
	while (tot_number < N - sum)
	{
		p1 = rand() % N;
		if (Flag[p1] == 0)
		{

			while (1)
			{
				count = rand() % K;
				if (SizeGroup[count] < UB[count])
				{
					p[p1] = count;
					Flag[p1] = 1;
					SizeGroup[count]++;
					tot_number++;
					break;
				}

			}
		}

	}

	for (i = 0; i < K; i++)  SizeG[i] = SizeGroup[i];
	delete[] SizeGroup; SizeGroup = NULL;
	delete[] Flag; Flag = NULL;
}

void BuildNeighbors()
{
	/* Initializes the neighbor structure for optimization */
	int i, j, g;
	int count;
	int SN = N * (N - 1) / 2 + N * K;
	count = 0;
	
	for (i = 0; i < N; i++)
		for (g = 0; g < K; g++)
		{
			Neighbors[count].type = 1;
			Neighbors[count].v = i;
			Neighbors[count].g = g;
			count++;
		}
	for (i = 0; i < N; i++)
		for (j = i + 1; j < N; j++)
		{
			Neighbors[count].type = 2;
			Neighbors[count].x = i;
			Neighbors[count].y = j;
			count++;
		}
		
}

void Clear_Delta_Matrix()
{
	/* Resets the delta matrix and objective function value */
	int x, g;
	f = 0.0;
	for (x = 0; x < N; x++)
		for (g = 0; g < K; g++)
			Delta_Matrix[x][g] = 0.0;
	return;
}

void Build_Delta_Matrix()
{
	/*  Builds the delta matrix and calculates the objective function value */
	int i, j;
	Clear_Delta_Matrix();

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
		{
			Delta_Matrix[i][p[j]] += D[i][j];
		}

	for (i = 0; i < N; i++)
		for (j = 0; j < K; j++)
			Delta[i][j] = Delta_Matrix[i][j] - Delta_Matrix[i][p[i]];

	f = 0.0;
	for (i = 0; i < N; i++)
		f += Delta_Matrix[i][p[i]];
	f = f / 2;

	return;
}

void Build_GroupDiv_For_Crossover() {
	/*  Builds group diversity values for crossover */
	int i, j;
	for (i = 0; i < K; i++) {
		gDiv[i] = 0.0;
	}
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
		{
			if (p[i] == p[j]) {
				gDiv[p[i]] += D[i][j];
			}

		}
}

void One_Move_Update_Delta_Matrix(int i, int g0, int g1)
{
	/* udates the delta matrix for single-element moves */
	int x, j, k;

	for (j = 0; j < N; j++)
	{
		if (j != i)
		{
			Delta_Matrix[j][g0] -= D[i][j];
			Delta_Matrix[j][g1] += D[i][j];
		}
	}

	for (x = 0; x < N; x++)
	{
		if (x != i)
		{
			Delta[x][g0] = Delta_Matrix[x][g0] - Delta_Matrix[x][p[x]];
			Delta[x][g1] = Delta_Matrix[x][g1] - Delta_Matrix[x][p[x]];

			if (p[x] == g0)
			{
				for (k = 0; k < K; k++)
				{
					Delta[x][k] = Delta_Matrix[x][k] - Delta_Matrix[x][g0];
				}
			}

			if (p[x] == g1)
			{
				for (k = 0; k < K; k++)
				{
					Delta[x][k] = Delta_Matrix[x][k] - Delta_Matrix[x][g1];
				}
			}

		}
	}
	x = i;
	for (k = 0; k < K; k++) Delta[x][k] = Delta_Matrix[x][k] - Delta_Matrix[x][g1];

	return;
}
void One_Move_Update_Delta_Matrix1(int i, int g0, int g1)
{
	int x, j, k;

	for (j = 0; j < N; j++)
	{
		if (j != i)
		{
			Delta_Matrix[j][g0] -= D[i][j];
			Delta_Matrix[j][g1] += D[i][j];
		}
	}
	return;
}

void RandLS(int partition[], int SizeGroup[], double* cost)
{
	/* Algorithm 3: Double-neighborhood local search method.
	 Performs a random local search and updates the cost */
	int i, v, g, x, y;
	int old_g, old_g1, swap;
	double delt;
	int Flag;

	for (i = 0; i < N; i++) p[i] = partition[i];
	Build_Delta_Matrix();
	f_best = f;

	delt = -99999.0;
	do
	{
		Flag = 0;

		for (v = 0; v < N; v++)
			for (g = 0; g < K; g++)
				if ((p[v] != g) && (SizeGroup[p[v]] > LB[p[v]]) && (SizeGroup[g] < UB[g]))
				{
					delt = Delta_Matrix[v][g] - Delta_Matrix[v][p[v]];
					if (delt > 0.0001)
					{
						old_g = p[v];
						One_Move_Update_Delta_Matrix1(v, old_g, g);
						SizeGroup[old_g] = SizeGroup[old_g] - 1;
						SizeGroup[g] = SizeGroup[g] + 1;
						p[v] = g;
						f += delt;
						Flag = 1;
					}

				}

		for (x = 0; x < N; x++)
			for (y = x + 1; y < N; y++)
				if (p[x] != p[y])
				{
					delt = (Delta_Matrix[x][p[y]] - Delta_Matrix[x][p[x]]) + (Delta_Matrix[y][p[x]] - Delta_Matrix[y][p[y]]) - DT[x][y];
					if (delt > 0.0001)
					{
						old_g = p[x];
						old_g1 = p[y];
						One_Move_Update_Delta_Matrix1(x, old_g, old_g1);
						One_Move_Update_Delta_Matrix1(y, old_g1, old_g);

						swap = p[x];
						p[x] = p[y];
						p[y] = swap;

						f += delt;
						Flag = 1;
					}

				}
	} while (Flag == 1);

	for (i = 0; i < N; i++)  partition[i] = p[i];
	*cost = f;
	
}

void StrongPerturbation(int L, int partition[], int SizeGroup[])
{
	/* Algorithm 4: Undirected Perturbation. Applies a strong perturbation to the partition */
	int i, v, g, x, y;
	int NumberNeighbors, old_g, old_g1, swap;
	int iter = 0;
	int cur_index, starting_index;
	double delt;
	double total_time = 0.0, starting_time = 0.0;
	int theta, count = 0;

	theta = L;

	NumberNeighbors = N * (N - 1) / 2 + N * K;
	for (i = 0; i < N; i++) p[i] = partition[i];

	do
	{

		cur_index = rand() % NumberNeighbors;

		if (Neighbors[cur_index].type == 1)
		{
			v = Neighbors[cur_index].v;
			g = Neighbors[cur_index].g;
			if ((p[v] != g) && (SizeGroup[p[v]] > LB[p[v]]) && (SizeGroup[g] < UB[g]))
			{
				delt = Delta_Matrix[v][g] - Delta_Matrix[v][p[v]];

				old_g = p[v];
				
				SizeGroup[old_g] = SizeGroup[old_g] - 1;
				SizeGroup[g] = SizeGroup[g] + 1;
				p[v] = g;
				count++;

			}

		}

		else if (Neighbors[cur_index].type == 2)
		{
			x = Neighbors[cur_index].x;
			y = Neighbors[cur_index].y;
			if (p[x] != p[y])
			{
				delt = (Delta_Matrix[x][p[y]] - Delta_Matrix[x][p[x]]) + (Delta_Matrix[y][p[x]] - Delta_Matrix[y][p[y]]) - 2 * D[x][y];

				old_g = p[x];
				old_g1 = p[y];

				swap = p[x];
				p[x] = p[y];
				p[y] = swap;

				count++;
			}
		}
	} while (count < theta);

	for (i = 0; i < N; i++)  partition[i] = p[i];
}

void DirectPerturbation(int LMAX, int partition[], int SizeGroup[])
{
	/* Algorithm 6: Directed Perturbation. 
	Iteratively refines partitions to balance group sizes and minimize costs */
	int i, j, v, g, x, y, k, L, number, Minsd, MinE, Flag;
	double delt, delt1, delt_max;
	int NumberNeighbors, old_g, old_g1, swap;
	for (i = 0; i < N; i++) p[i] = partition[i];
	for (j = 0; j < K; j++) SizeG[j] = SizeGroup[j];
	Build_Delta_Matrix();
	for (L = 0; L < LMAX; L++)
	{
		number = 0;
		for (i = 0; i < K; i++) {
			UnderLB[i] = 0;
			Rd[i] = -1;
			for (j = 0; j < K; j++)
			{
				AvgCon[i][j] = 0.0;
			}
		}
		
		for (k = 0; k < K; k++) {
			Minsd = 99999999;
			MinE = 0;
			for (i = 0; i < N; i++) {
				if (p[i] == k) {
					if (Delta_Matrix[i][k] < Minsd) {
						Minsd = Delta_Matrix[i][k];
						MinE = i;
					}
				}
			}
			Rd[k] = MinE;
			SizeG[k] -= 1;
			if (SizeG[k] < LB[k]) {
				UnderLB[k] = 1;
				number += 1;
			}
		}
		
		for (i = 0; i < K; i++) {
			for (j = 0; j < K; j++) {
				Delta_Matrix[Rd[i]][p[Rd[j]]] = Delta_Matrix[Rd[i]][p[Rd[j]]] - D[Rd[i]][Rd[j]];
				AvgCon[p[Rd[i]]][p[Rd[j]]] = Delta_Matrix[Rd[i]][p[Rd[j]]] / SizeG[p[Rd[j]]];
			}
		}
		
		int Elep;
		int maxValue;
		int nn = 0;
		int i;
		while (nn < number) {
			maxValue = -9999;
			i = rand() % K;
			do {
				i = (i + 1) % K;
			} while (UnderLB[i] == 0);
			for (j = 0; j < K; j++) {
				if (AvgCon[j][i] > maxValue) {
					maxValue = AvgCon[j][i];
					Elep = j;
				}
			}
			SizeG[i] += 1;
			for (k = 0; k < K; k++) {
				if (Rd[k] != -1) {
					Delta_Matrix[Rd[k]][i] += D[Rd[k]][Rd[Elep]];
					AvgCon[p[Rd[k]]][i] = Delta_Matrix[Rd[k]][i] / SizeG[i];
				}
			}
			for (k = 0; k < K; k++) {
				AvgCon[p[Rd[Elep]]][k] = 0.0;
			}
			p[Rd[Elep]] = i;
			UnderLB[i] = 0;
			Rd[Elep] = -1;
			nn++;
		}
		int GP;
		nn = 0;
		while (nn < K - number) {
			Elep = rand() % K;
			do {
				Elep = (Elep + 1) % K;
			} while (Rd[Elep] == -1);
			maxValue = -9999;
			for (j = 0; j < K; j++) {
				if (AvgCon[Elep][j] > maxValue) {
					maxValue = AvgCon[Elep][j];
					GP = j;
				}
			}
			if (SizeG[GP] < UB[GP]) {
				SizeG[GP] += 1;
				for (k = 0; k < K; k++) {
					if (Rd[k] != -1) {
						Delta_Matrix[Rd[k]][GP] += D[Rd[k]][Rd[Elep]];
						AvgCon[p[Rd[k]]][GP] = Delta_Matrix[Rd[k]][GP] / SizeG[GP];
					}
				}
				for (k = 0; k < K; k++) {
					AvgCon[p[Rd[Elep]]][k] = 0.0;
				}
				p[Rd[Elep]] = GP;
				Rd[Elep] = -1;
				nn += 1;
			}
			else {
				for (k = 0; k < K; k++) {
					AvgCon[k][GP] = 0.0;
				}
			}
		}
		Build_Delta_Matrix();
	}
	for (i = 0; i < N; i++)  partition[i] = p[i];
	for (j = 0; j < K; j++)  SizeGroup[j] = SizeG[j];
}

void Crossover(int partition1[], int partition2[], int sc[], int scSizeGroup[]) {
	/* Algorithm 5: combines partitions in a way that maintains group constraints */
	int i, j, k, gDivMax, g;
	int lengthSE;
	int lengthSG;
	int flag;
	int num;
	int pickG;
	int count;
	int pickV;
	int sum;
	int sumLB;
	int sumLowerThanLB;
	for (i = 0; i < N; i++) {
		p[i] = partition1[i];
		p1[i] = partition1[i];
	}
	Build_Delta_Matrix();
	for (i = 0; i < N; i++) {
		for (j = 0; j < K; j++) {
			Delta_Matrix_p1[i][j] = Delta_Matrix[i][j];
		}
	}
	Build_GroupDiv_For_Crossover();
	for (i = 0; i < K; i++) {
		gDiv_p1[i] = gDiv[i];
	}

	for (i = 0; i < N; i++) {
		p[i] = partition2[i];
		p2[i] = partition2[i];
	}
	Build_Delta_Matrix();
	for (i = 0; i < N; i++) {
		for (j = 0; j < K; j++) {
			Delta_Matrix_p2[i][j] = Delta_Matrix[i][j];
		}
	}
	Build_GroupDiv_For_Crossover();
	for (i = 0; i < K; i++) {
		gDiv_p2[i] = gDiv[i];
	}

	for (i = 0; i < N; i++) {
		vEle[i] = i;
		sc[i] = -1;
	}
	for (i = 0; i < K; i++) {
		LBGroup[i] = 0;
		UBGroup[i] = 0;
		BigThanLB[i] = 0;
		gEle[i] = i;
		ub[i] = UB[i];
		scSizeGroup[i] = 0;
	}

	for (i = 0; i < K; i++) {
		if (rand() / (RAND_MAX + 1.0) < 0.5) {
			gDivMax = -9999;
			for (j = 0; j < K; j++) {
				if (gDiv_p1[j] > gDivMax) {
					gDivMax = gDiv_p1[j];
					g = j;
				}
			}
			lengthSE = 0;
			for (j = 0; j < N; j++) {
				if (p1[j] == g) {
					SelectEle[lengthSE++] = j;
				}
			}
			lengthSG = 0;
			for (j = 0; j < K; j++) {
				if (ub[j] != -1 && ub[j] >= lengthSE) {
					SelGroup[lengthSG++] = j;
				}
			}
			// no group satisfied
			if (lengthSG == 0) {
				num = 999999;
				for (j = 0; j < K; j++) {
					if (ub[j] != -1 && lengthSE - ub[j] < num) {
						num = lengthSE - ub[j];
						pickG = j;
					}
				}

				count = 0;
				while (count < lengthSE - num) {
					pickV = rand() % lengthSE;
					do
					{
						pickV = (pickV + 1) % lengthSE;
					} while (SelectEle[pickV] == -1);
					sc[SelectEle[pickV]] = pickG;
					SelectEleTemp[count++] = SelectEle[pickV];
					vEle[SelectEle[pickV]] = -1;
					SelectEle[pickV] = -1;


				}
				lengthSE = count;
			}
			else {
				pickG = rand() % lengthSG;
				pickG = SelGroup[pickG];
				for (j = 0; j < lengthSE; j++) {
					sc[SelectEle[j]] = pickG;
					vEle[SelectEle[j]] = -1;
					SelectEleTemp[j] = SelectEle[j];
				}
			}
		}
		else {
			gDivMax = -9999;
			for (j = 0; j < K; j++) {
				if (gDiv_p2[j] > gDivMax) {
					gDivMax = gDiv_p2[j];
					g = j;
				}
			}
			lengthSE = 0;
			for (j = 0; j < N; j++) {
				if (p2[j] == g) {
					SelectEle[lengthSE++] = j;
				}
			}
			lengthSG = 0;
			for (j = 0; j < K; j++) {
				if (ub[j] != -1 && ub[j] >= lengthSE) {
					SelGroup[lengthSG++] = j;
				}
			}
			// no group satisfied
			if (lengthSG == 0) {
				num = 999999;
				for (j = 0; j < K; j++) {
					if (ub[j] != -1 && lengthSE - ub[j] < num) {
						num = lengthSE - ub[j];
						pickG = j;
					}
				}

				count = 0;
				while (count < lengthSE - num) {
					pickV = rand() % lengthSE;
					do
					{
						pickV = (pickV + 1) % lengthSE;
					} while (SelectEle[pickV] == -1);
					sc[SelectEle[pickV]] = pickG;
					SelectEleTemp[count++] = SelectEle[pickV];
					vEle[SelectEle[pickV]] = -1;
					SelectEle[pickV] = -1;


				}
				lengthSE = count;
			}
			else {
				pickG = rand() % lengthSG;
				pickG = SelGroup[pickG];
				for (j = 0; j < lengthSE; j++) {
					sc[SelectEle[j]] = pickG;
					vEle[SelectEle[j]] = -1;
					SelectEleTemp[j] = SelectEle[j];
				}
			}
		}
		for (j = 0; j < lengthSE; j++) {
			gDiv_p1[p1[SelectEleTemp[j]]] -= Delta_Matrix_p1[SelectEleTemp[j]][p1[SelectEleTemp[j]]];
			gDiv_p2[p2[SelectEleTemp[j]]] -= Delta_Matrix_p2[SelectEleTemp[j]][p2[SelectEleTemp[j]]];
			p1[SelectEleTemp[j]] = -1;
			p2[SelectEleTemp[j]] = -1;
		}
		ub[pickG] = -1;
		scSizeGroup[pickG] = lengthSE;
	}
	count = 0;
	sumLB = 0;
	sumLowerThanLB = 0;
	for (i = 0; i < K; i++) {
		sumLB += LB[i];
		if (scSizeGroup[i] < LB[i]) {
			count = count + scSizeGroup[i];
			sumLowerThanLB += scSizeGroup[i];
			LBGroup[i] = 1;
		}
		else {
			count = count + LB[i];
		}

		if (scSizeGroup[i] > LB[i]) {
			BigThanLB[i] = 1;
		}
	}

	for (i = 0; i < N; i++) {
		if (vEle[i] != -1) {
			count++;
		}
	}
	while (count < sumLB) {
		pickG = rand() % K;
		do {
			pickG = (pickG + 1) % K;
		} while (BigThanLB[pickG] == 0);
		lengthSE = 0;
		for (j = 0; j < N; j++) {
			if (sc[j] == pickG) {
				SelectEle[lengthSE++] = j;
			}
		}
		pickV = rand() % lengthSE;
		sc[SelectEle[pickV]] = -1;
		vEle[SelectEle[pickV]] = SelectEle[pickV];
		scSizeGroup[pickG]--;
		if (scSizeGroup[pickG] == LB[pickG]) {
			BigThanLB[pickG] = 0;
		}
		count++;
	}
	sum = 0;
	for (i = 0; i < K; i++) {
		if (LBGroup[i] == 1) {
			sum += LB[i];
		}

	}

	while (sumLowerThanLB < sum) {
		pickG = rand() % K;
		do {
			pickG = (pickG + 1) % K;
		} while (LBGroup[pickG] == 0);
		lengthSE = 0;
		for (i = 0; i < N; i++) {
			if (vEle[i] != -1) {
				SelectEle[lengthSE++] = i;
			}
		}
		pickV = rand() % lengthSE;
		sc[SelectEle[pickV]] = pickG;
		vEle[SelectEle[pickV]] = -1;
		scSizeGroup[pickG]++;
		if (scSizeGroup[pickG] == LB[pickG]) {
			LBGroup[pickG] = 0;
		}
		sumLowerThanLB++;
	}
	sum = 0;
	for (i = 0; i < K; i++) {
		sum += scSizeGroup[i];
		if (scSizeGroup[i] < UB[i]) {
			UBGroup[i] = 1;
		}
	}

	while (sum < N) {
		pickG = rand() % K;
		do {
			pickG = (pickG + 1) % K;
		} while (UBGroup[pickG] == 0);
		lengthSE = 0;
		for (i = 0; i < N; i++) {
			if (vEle[i] != -1) {
				SelectEle[lengthSE++] = i;
			}
		}
		pickV = rand() % lengthSE;
		sc[SelectEle[pickV]] = pickG;
		vEle[SelectEle[pickV]] = -1;
		scSizeGroup[pickG]++;
		if (scSizeGroup[pickG] == UB[pickG]) {
			UBGroup[pickG] = 0;
		}
		sum++;
	}


}
double FitRadioAndDis(int partition1[], int partition2[], double cost1, double cost2) {
	/* evaluates the quality and dissimilarity of partitions */
	double radio;
	int i, j;
	int count;
	count = 0;
	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			if (partition1[i] == partition1[j] && partition2[i] != partition2[j] || partition1[i] != partition1[j] && partition2[i] == partition2[j]) {
				count++;
			}
		}
	}
	radio = cost1 / cost2 + 0.05 * (count / (N * N) * K);
	return radio;
}


void InitialSol(Solution& S)
{
	/* Algorithm 2: initializes a solution S */
	int i, j;
	int counter = 0;
	RandomInitiaSol(S.p, S.SizeG);
	RandLS(S.p, S.SizeG, &S.cost);
}

bool Cmpare(const Solution &a, const Solution &b) {
	/* ompares two solutions based on their cost */
	return a.cost > b.cost;
}


void SearchAlgorithm()
{
	/* Algorithm 1: The main procedure of TPSDP. */
	int i, j, k;
	int L;
	
	int radio;
	int pickS;
	starting_time = clock();
	GS.cost = -9999999;
	for (i = 0; i < popSize; i++) {
		InitialSol(CS);
		for (j = 0; j < N; j++) Pop[i].p[j] = CS.p[j];
		for (k = 0; k < K; k++) Pop[i].SizeG[k] = CS.SizeG[k];
		Pop[i].cost = CS.cost;
		if (Pop[i].cost > GS.cost) {
			for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
			for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
			GS.cost = Pop[i].cost;
		}
	}

	while (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC < Time_limit)
	{
		L = int(theta * N / K);
		for (i = 0; i < popSize; i++) {
			for (j = 0; j < N; j++) Offs[i].p[j] = Pop[i].p[j];
			for (k = 0; k < K; k++) Offs[i].SizeG[k] = Pop[i].SizeG[k];
			Offs[i].cost = Pop[i].cost;
		}
		// 1 StrongPerturbation and LocalSearch
		
		for (i = 0; i < popSize; i++) {
			StrongPerturbation(L, Pop[i].p, Pop[i].SizeG);
			RandLS(Pop[i].p, Pop[i].SizeG, &Pop[i].cost);
			if (Pop[i].cost > GS.cost) {
				for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
				for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
				GS.cost = Pop[i].cost;
			}
		}
		
		
		// 2 Crossover and LocalSearch
		
		if(popSize > 1){
			for (i = 0; i < popSize; i++) {
				pickS = rand() % popSize;
				do {
					pickS = (pickS + 1) % popSize;
				} while (pickS == i);
				Crossover(Pop[i].p, Pop[pickS].p, Offs[i].p, Offs[i].SizeG);
				RandLS(Offs[i].p, Offs[i].SizeG, &Offs[i].cost);
			}
			for (i = 0; i < popSize; i++) {
				if (Offs[i].cost >= Pop[i].cost) {
					for (j = 0; j < N; j++) Pop[i].p[j] = Offs[i].p[j];
					for (k = 0; k < K; k++) Pop[i].SizeG[k] = Offs[i].SizeG[k];
					Pop[i].cost = Offs[i].cost;
				}
				else if (FitRadioAndDis(Offs[i].p, Pop[i].p, Offs[i].cost, Pop[i].cost) > 1) {
					for (j = 0; j < N; j++) Pop[i].p[j] = Offs[i].p[j];
					for (k = 0; k < K; k++) Pop[i].SizeG[k] = Offs[i].SizeG[k];
					Pop[i].cost = Offs[i].cost;
				}
				if (Pop[i].cost > GS.cost) {
					for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
					for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
					GS.cost = Pop[i].cost;
				}
			}
		}
		
		
		//3 DirectPerturbation and LocalSearch
		
		for (i = 0; i < popSize; i++) {
			DirectPerturbation(LMAX, Pop[i].p, Pop[i].SizeG);
			RandLS(Pop[i].p, Pop[i].SizeG, &Pop[i].cost);
			if (Pop[i].cost > GS.cost) {
				for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
				for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
				GS.cost = Pop[i].cost;
			}
		}
		
		// 4 liner decrease population
		  sort(Pop, Pop + popSize, Cmpare);
		  popSize = int(beta_min - popSize) / Time_limit * (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC) + popSize;
		  theta = theta_max - (theta_max-theta_min) * (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC) / Time_limit;
		  
		
	}
}

// [[Rcpp::export]]
Rcpp::List three_phase_search_dynamic_population_size(Rcpp::NumericMatrix matrix, 
                                            int N, 
											int K, 
                                            int upper_bound, 
											int lower_bound, 
											int popSize, 
											int time_limit,
											double theta_max,
											double theta_min,
											int beta_min,
											int LMAX
                                           )
{
	/* receives data from r, calls SearchAlgorithm, save results for r
	*/	

	N = N;
    Rcpp::Rcout << "Thanks for viewing my code!" << std::endl;
	K = K;
	popSize = popSize;  //beta_max in algortihm
	theta = theta_max;
	theta_max = theta_max;
	beta_min = beta_min;
	LMAX = LMAX;

	D = new double* [N];
	DT = new double* [N];
	int i;
	for (i = 0; i < N; i++) {
		D[i] = new double[N];
		DT[i] = new double[N];
	} 

	for (i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			D[i][j] = matrix(i, j);
			DT[i][j] = 2* matrix(i,j);
		}
	}
	
	LB = new int[K]; 
	UB = new int[K];
	for (i = 0; i < K; i++) { LB[i] = lower_bound; UB[i] = upper_bound; }
	
	AssignMemery();
	
	BuildNeighbors();

	SearchAlgorithm();

	//save GS -> solution with result
	Rcpp::IntegerVector result(N);
	for (int i = 0; i < N; i++) {
        result[i] = GS.p[i];
    }
	int cost = GS.cost;

	ReleaseMemery();

 	return Rcpp::List::create(Rcpp::Named("cost") = cost, Rcpp::Named("result") = result);
}
