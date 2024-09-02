#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <R.h>


typedef struct {
    int *p;        // cluster membership of each element
    int *SizeG;    // size of each cluster
    double cost;   // global cost function value of solution
} Solution;

typedef struct {
    int type;
    int v;
    int g;
    int x;
    int y;
} Neighborhood;

char *File_Name;
char *Output_File_Name;
char *Solution_File;
int popSize;
int N, K;  // node number and group number
double f, f_best;
Solution CS, NS, GS, OS;
Solution *Pop;
Solution *Offs;
Neighborhood *Neighbors;
double total_time, starting_time, Time_limit;
int *p; // partition array for each vertex
int *bestp; // partition array for each vertex
int *SizeG;
double **Delta_Matrix;  // incremental matrix 
double **Delta_Matrix_p1;
double **Delta_Matrix_p2;
double **Delta_Matrix1;
double **Delta;
double **Delta1;
double *gDiv;
double *gDiv_p1;
double *gDiv_p2;
int *SelectEle;
int *SelectEleTemp;
int *SelGroup;
int *p1;
int *p2;

// for crossover
int *vEle;
int *gEle;
int *LBGroup;
int *UBGroup;
int *BigThanLB;
int *ub;

double** D;   // distance matrix between elements
double** DT;
int* LB; // Lower bound of the number of elements in the group i 
int* UB; // Upper bound of the number of elements in the group i
double theta, theta_max, theta_min; 
int beta_min; 
int LMAX;

double** AvgCon;
int* Rd, * UnderLB;

void Build_Delta_Matrix();
void One_Move_Update_Delta_Matrix1(int v, int old_g, int new_g);
void Build_GroupDiv_For_Crossover();
void Build_Delta_Matrix(void);
void Build_GroupDiv_For_Crossover(void);
void AssignMemory();
void ReleaseMemory();
void InitialSol(Solution *S);
void StrongPerturbation(int L, int partition[], int SizeGroup[]);
void RandLS(int partition[], int SizeGroup[], double* cost);
void Crossover(int partition1[], int partition2[], int sc[], int scSizeGroup[]);
double FitRadioAndDis(int partition1[], int partition2[], double cost1, double cost2);
void RandomInitiaSol(int p[], int SizeG[]);
void DirectPerturbation(int LMAX, int partition[], int SizeGroup[]);
int Cmpare(const void *a, const void *b);

void three_phase_search_dynamic_population_size(double *D, 
                                            double *DT,
                                            int *N,
                                            int *K, 
                                            int *upper_bound, int *lower_bound, 
                                            int *result,
											int	*cost,
											int *popSize, 
											int *time_limit,
											double *theta_max,
											double *theta_min,
											int *beta_min,
											int *LMAX
                                           ) {

}


void SearchAlgorithm() {
    int i, j, k;
    int L;
    int pickS;
    clock_t starting_time = clock();
    GS.cost = -9999999;
    
    // Initial population generation
    for (i = 0; i < popSize; i++) {
        InitialSol(&CS);
        for (j = 0; j < N; j++) Pop[i].p[j] = CS.p[j];
        for (k = 0; k < K; k++) Pop[i].SizeG[k] = CS.SizeG[k];
        Pop[i].cost = CS.cost;
        if (Pop[i].cost > GS.cost) {
            for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
            for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
            GS.cost = Pop[i].cost;
        }
    }

    while (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC < Time_limit) {
        L = (int)(theta * N / K);
        for (i = 0; i < popSize; i++) {
            for (j = 0; j < N; j++) Offs[i].p[j] = Pop[i].p[j];
            for (k = 0; k < K; k++) Offs[i].SizeG[k] = Pop[i].SizeG[k];
            Offs[i].cost = Pop[i].cost;
        }
        // Strong Perturbation and Local Search
        for (i = 0; i < popSize; i++) {
            StrongPerturbation(L, Pop[i].p, Pop[i].SizeG);
            RandLS(Pop[i].p, Pop[i].SizeG, &Pop[i].cost);
            if (Pop[i].cost > GS.cost) {
                for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
                for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
                GS.cost = Pop[i].cost;
            }
        }

        // Crossover and Local Search
        if (popSize > 1) {
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
                } else if (FitRadioAndDis(Offs[i].p, Pop[i].p, Offs[i].cost, Pop[i].cost) > 1) {
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

        // Direct Perturbation and Local Search
        for (i = 0; i < popSize; i++) {
            DirectPerturbation(LMAX, Pop[i].p, Pop[i].SizeG);
            RandLS(Pop[i].p, Pop[i].SizeG, &Pop[i].cost);
            if (Pop[i].cost > GS.cost) {
                for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
                for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
                GS.cost = Pop[i].cost;
            }
        }

        // Linear decrease population size
        // Note: Implement sort function based on the comparison function `Cmpare`
        qsort(Pop, popSize, sizeof(Solution), Cmpare);
        popSize = (int)(beta_min - popSize) / Time_limit * (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC) + popSize;
        theta = theta_max - (theta_max - theta_min) * (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC) / Time_limit;
    }
}

void InitialSol(Solution *S) {
    /* Algorithm 2: initializes a solution S */
    RandomInitiaSol(S->p, S->SizeG);
    RandLS(S->p, S->SizeG, &(S->cost));
}

int Cmpare(const void *a, const void *b) {
    /* Compares two solutions based on their cost */
    Solution *solA = (Solution *)a;
    Solution *solB = (Solution *)b;
    return (solB->cost - solA->cost); // Return positive if b is greater, negative if a is greater
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
	SizeGroup =(int *)malloc(K * sizeof(int));;
	for (i = 0; i < K; i++) SizeGroup[i] = 0;
	Flag = (int *)malloc(N * sizeof(int));;
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
	free(SizeGroup);; SizeGroup = NULL;
	free(Flag);; Flag = NULL;
}

void RandLS(int partition[], int SizeGroup[], double* cost) {
    int i, v, g, x, y;
    int old_g, old_g1, swap;
    double delt;
    int Flag;

    for (i = 0; i < N; i++) p[i] = partition[i];
    Build_Delta_Matrix();
    delt = -99999.0;
    do {
        Flag = 0;
        for (v = 0; v < N; v++) {
            for (g = 0; g < K; g++) {
                if ((p[v] != g) && (SizeGroup[p[v]] > LB[p[v]]) && (SizeGroup[g] < UB[g])) {
                    delt = Delta_Matrix[v][g] - Delta_Matrix[v][p[v]];
                    if (delt > 0.0001) {
                        old_g = p[v];
                        One_Move_Update_Delta_Matrix1(v, old_g, g);
                        SizeGroup[old_g] -= 1;
                        SizeGroup[g] += 1;
                        p[v] = g;
                        *cost += delt;
                        Flag = 1;
                    }
                }
            }
        }
        for (x = 0; x < N; x++) {
            for (y = x + 1; y < N; y++) {
                if (p[x] != p[y]) {
                    delt = (Delta_Matrix[x][p[y]] - Delta_Matrix[x][p[x]]) + (Delta_Matrix[y][p[x]] - Delta_Matrix[y][p[y]]) - DT[x][y];
                    if (delt > 0.0001) {
                        old_g = p[x];
                        old_g1 = p[y];
                        One_Move_Update_Delta_Matrix1(x, old_g, old_g1);
                        One_Move_Update_Delta_Matrix1(y, old_g1, old_g);
                        swap = p[x];
                        p[x] = p[y];
                        p[y] = swap;
                        *cost += delt;
                        Flag = 1;
                    }
                }
            }
        }
    } while (Flag == 1);

    for (i = 0; i < N; i++) partition[i] = p[i];
}

void StrongPerturbation(int L, int partition[], int SizeGroup[]) {
    int i, v, g, x, y;
    int NumberNeighbors, old_g, old_g1, swap;
    int cur_index;
    double delt;
    int theta, count = 0;

    theta = L;
    NumberNeighbors = N * (N - 1) / 2 + N * K;
    for (i = 0; i < N; i++) p[i] = partition[i];

    do {
        cur_index = rand() % NumberNeighbors;
        if (Neighbors[cur_index].type == 1) { // Adjust condition according to the structure of Neighbors
            v = Neighbors[cur_index].v;
            g = Neighbors[cur_index].g;
            if ((p[v] != g) && (SizeGroup[p[v]] > LB[p[v]]) && (SizeGroup[g] < UB[g])) {
                delt = Delta_Matrix[v][g] - Delta_Matrix[v][p[v]];
                old_g = p[v];
                SizeGroup[old_g] -= 1;
                SizeGroup[g] += 1;
                p[v] = g;
                count++;
            }
        } else if (Neighbors[cur_index].type == 2) { // Adjust condition according to the structure of Neighbors
            x = Neighbors[cur_index].x;
            y = Neighbors[cur_index].y;
            if (p[x] != p[y]) {
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

    for (i = 0; i < N; i++) partition[i] = p[i];
}

void DirectPerturbation(int LMAX, int partition[], int SizeGroup[]) {
    int i, j, v, g, x, y, k, L, number, Minsd, MinE, Flag;
    double delt, delt1, delt_max;
    int NumberNeighbors, old_g, old_g1, swap;
    for (i = 0; i < N; i++) p[i] = partition[i];
    for (j = 0; j < K; j++) SizeG[j] = SizeGroup[j];
    Build_Delta_Matrix();

    for (L = 0; L < LMAX; L++) {
        number = 0;
        for (i = 0; i < K; i++) {
            UnderLB[i] = 0;
            Rd[i] = -1;
            for (j = 0; j < K; j++) {
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
                if (AvgCon[j][Elep] > maxValue) {
                    maxValue = AvgCon[j][Elep];
                    GP = j;
                }
            }
            SizeG[Elep] += 1;
            for (k = 0; k < K; k++) {
                if (Rd[k] != -1) {
                    Delta_Matrix[Rd[k]][Elep] += D[Rd[k]][Rd[GP]];
                    AvgCon[p[Rd[k]]][Elep] = Delta_Matrix[Rd[k]][Elep] / SizeG[Elep];
                }
            }
            for (k = 0; k < K; k++) {
                AvgCon[p[Rd[GP]]][k] = 0.0;
            }
            p[Rd[GP]] = Elep;
            Rd[GP] = -1;
            nn++;
        }
    }

    for (i = 0; i < N; i++) partition[i] = p[i];
}

void Crossover(int partition1[], int partition2[], int sc[], int scSizeGroup[]) {
    int i, j, k, gDivMax, g;
    int lengthSE, lengthSG;
    int flag;
    int num;
    int pickG;
    int count;
    int pickV;
    int sum;
    int sumLB;
    int sumLowerThanLB;

    // Initialize p and p1
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

    // Initialize p and p2
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

    // Initialize arrays
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

    // Main crossover process
    for (i = 0; i < K; i++) {
        if ((rand() / (RAND_MAX + 1.0)) < 0.5) {
            // Process for partition1
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
            // No group satisfied
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
                    do {
                        pickV = (pickV + 1) % lengthSE;
                    } while (SelectEle[pickV] == -1);
                    sc[SelectEle[pickV]] = pickG;
                    SelectEleTemp[count++] = SelectEle[pickV];
                    vEle[SelectEle[pickV]] = -1;
                    SelectEle[pickV] = -1;
                }
                lengthSE = count;
            } else {
                pickG = rand() % lengthSG;
                pickG = SelGroup[pickG];
                for (j = 0; j < lengthSE; j++) {
                    sc[SelectEle[j]] = pickG;
                    vEle[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        } else {
            // Process for partition2
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
            // No group satisfied
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
                    do {
                        pickV = (pickV + 1) % lengthSE;
                    } while (SelectEle[pickV] == -1);
                    sc[SelectEle[pickV]] = pickG;
                    SelectEleTemp[count++] = SelectEle[pickV];
                    vEle[SelectEle[pickV]] = -1;
                    SelectEle[pickV] = -1;
                }
                lengthSE = count;
            } else {
                pickG = rand() % lengthSG;
                pickG = SelGroup[pickG];
                for (j = 0; j < lengthSE; j++) {
                    sc[SelectEle[j]] = pickG;
                    vEle[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        }
        // Update group division values
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
            count += scSizeGroup[i];
            sumLowerThanLB += scSizeGroup[i];
            LBGroup[i] = 1;
        } else {
            count += LB[i];
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
    /* Evaluates the quality and dissimilarity of partitions */
    double radio;
    int i, j;
    int count = 0;

    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            if ((partition1[i] == partition1[j] && partition2[i] != partition2[j]) || 
                (partition1[i] != partition1[j] && partition2[i] == partition2[j])) {
                count++;
            }
        }
    }

    radio = cost1 / cost2 + 0.05 * ((double)count / (N * N) * K);
    return radio;
}

void BuildNeighbors()
{
	/* Initializes the neighbor structure for optimization */
	int i, j, g;
	int count;
	// int SN = N * (N - 1) / 2 + N * K; ???
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


void AssignMemory() {
    int i, j;

    p = (int *)malloc(N * sizeof(int));
    bestp = (int *)malloc(N * sizeof(int));
    SizeG = (int *)malloc(K * sizeof(int));

    Pop = (Solution *)malloc(popSize * sizeof(Solution));
    Offs = (Solution *)malloc(popSize * sizeof(Solution));

    Delta_Matrix = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++) {
        Delta_Matrix[i] = (double *)malloc(K * sizeof(double));
    }
    Delta_Matrix_p1 = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++) {
        Delta_Matrix_p1[i] = (double *)malloc(K * sizeof(double));
    }
    Delta_Matrix_p2 = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++) {
        Delta_Matrix_p2[i] = (double *)malloc(K * sizeof(double));
    }
    Delta_Matrix1 = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++) {
        Delta_Matrix1[i] = (double *)malloc(K * sizeof(double));
    }
    Delta = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++) {
        Delta[i] = (double *)malloc(K * sizeof(double));
    }
    Delta1 = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++) {
        Delta1[i] = (double *)malloc(K * sizeof(double));
    }
    gDiv = (double *)malloc(K * sizeof(double));
    gDiv_p1 = (double *)malloc(K * sizeof(double));
    gDiv_p2 = (double *)malloc(K * sizeof(double));

    for (i = 0; i < popSize; i++) {
        Pop[i].p = (int *)malloc(N * sizeof(int));
        Offs[i].p = (int *)malloc(N * sizeof(int));
    }
    for (i = 0; i < popSize; i++) {
        Pop[i].SizeG = (int *)malloc(K * sizeof(int));
        Offs[i].SizeG = (int *)malloc(K * sizeof(int));
    }

    CS.p = (int *)malloc(N * sizeof(int));
    NS.p = (int *)malloc(N * sizeof(int));
    GS.p = (int *)malloc(N * sizeof(int));
    OS.p = (int *)malloc(N * sizeof(int));

    CS.SizeG = (int *)malloc(K * sizeof(int));
    NS.SizeG = (int *)malloc(K * sizeof(int));
    GS.SizeG = (int *)malloc(K * sizeof(int));
    OS.SizeG = (int *)malloc(K * sizeof(int));

    Neighbors = (Neighborhood *)malloc((N * (N - 1) / 2 + N * K) * sizeof(Neighborhood));

    AvgCon = (double **)malloc(K * sizeof(double *));
    for (i = 0; i < K; i++) {
        AvgCon[i] = (double *)malloc(K * sizeof(double));
    }
    Rd = (int *)malloc(K * sizeof(int));
    for (i = 0; i < K; i++) {
        Rd[i] = 0;
    }
    UnderLB = (int *)malloc(K * sizeof(int));

    ub = (int *)malloc(K * sizeof(int));
    LBGroup = (int *)malloc(K * sizeof(int));
    UBGroup = (int *)malloc(K * sizeof(int));
    BigThanLB = (int *)malloc(K * sizeof(int));
    vEle = (int *)malloc(N * sizeof(int));
    gEle = (int *)malloc(K * sizeof(int));
    SelectEle = (int *)malloc(N * sizeof(int));
    SelGroup = (int *)malloc(K * sizeof(int));
    SelectEleTemp = (int *)malloc(N * sizeof(int));
    p1 = (int *)malloc(N * sizeof(int));
    p2 = (int *)malloc(N * sizeof(int));
}

void ReleaseMemory() {
    int i;

    free(p);
    free(bestp);
    free(SizeG);

    free(CS.p);
    free(CS.SizeG);
    free(GS.p);
    free(GS.SizeG);
    free(NS.p);
    free(NS.SizeG);
    free(OS.p);
    free(OS.SizeG);

    free(LB);
    free(UB);
    free(Neighbors);

    for (i = 0; i < N; i++) {
        free(Delta_Matrix[i]);
        free(Delta_Matrix1[i]);
        free(Delta_Matrix_p1[i]);
        free(Delta_Matrix_p2[i]);
        free(Delta[i]);
        free(Delta1[i]);
        free(D[i]);
        free(DT[i]);
    }
    for (i = 0; i < K; i++) {
        free(AvgCon[i]);
    }

    for (i = 0; i < popSize; i++) {
        free(Pop[i].p);
        free(Offs[i].p);
        free(Pop[i].SizeG);
        free(Offs[i].SizeG);
    }
    free(Rd);
    free(UnderLB);
    free(SelectEle);
    free(SelGroup);
    free(SelectEleTemp);
    free(p1);
    free(p2);
    free(gDiv);
    free(gDiv_p1);
    free(gDiv_p2);

    free(ub);
    free(LBGroup);
    free(UBGroup);
    free(BigThanLB);
    free(vEle);
    free(gEle);
}
