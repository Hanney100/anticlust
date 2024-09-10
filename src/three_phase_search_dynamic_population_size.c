#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>


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

int populationSize;
int N, K;  // node number and group number
double objective_function, f_best;
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
double *groupDiversity;
double *groupDiversity_p1;
double *groupDiversity_p2;
int *SelectEle;
int *SelectEleTemp;
int *SelectGroup;
int *p1;
int *p2;

// for crossover
int *vectorElement;
int *groupElement;
int *LBGroup;
int *UBGroup;
int *BigThanLB;
int *ub;

double** Distances;   // distance matrix between elements
double** DistancesT;
int* LB; // Lower bound of the number of elements in the group i 
int* UB; // Upper bound of the number of elements in the group i
double theta, theta_max, theta_min; 
int beta_min; 
int eta_max;

double** AvgCon;
int* Rd, * UnderLB;

int random_int(int max);
double uniform_rnd_number(void);
void ClearDeltaMatrix();
void BuildDeltaMatrix();
void OneMoveUpdateDeltaMatrix1(int i, int oldGroup, int newGroup);
void OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup);
void BuildGroupDiversityForCrossover();
void AssignMemory();
void ReleaseMemory();
void BuildNeighbors();
void SearchAlgorithm();
void InitialSol(Solution *S);
void StrongPerturbation(int L, int partition[], int SizeGroup[]);
void RandLS(int partition[], int SizeGroup[], double* cost);
void Crossover(int partition1[], int partition2[], int score[], int scSizeGroup[]);
double FitRadioAndDis(int partition1[], int partition2[], double cost1, double cost2);
void RandomInitialSol(int p[], int SizeG[]);
void DirectPerturbation(int eta_max, int partition[], int SizeGroup[]);
int Cmpare(const void *a, const void *b);

void three_phase_search_dynamic_population_size(
                      double *distances, 
                      int *N_in,
                      int *K_in, 
                      int *upper_bound, 
                      int *lower_bound, 
											int *Beta_max, 
											int *time_limit,
											double *Theta_max,
											double *Theta_min,
											int *Beta_min,
											int *Lmax,
											int *result,
											double *cost,
											int *mem_error) {
  N = *N_in;
  K = *K_in;
  populationSize = *Beta_max;  //beta_max in algortihm
  theta = *Theta_max;
  theta_max = *Theta_max;
  beta_min = *Beta_min;
  eta_max = *Lmax;
  Time_limit =  *time_limit;

  
  // Allocate memory for Distances and DistancesT arrays
  Distances = (double**)malloc(N * sizeof(double*));
  if (Distances == NULL) { *mem_error = 1; return; }
  DistancesT = (double**)malloc(N * sizeof(double*));
  if (DistancesT == NULL) { *mem_error = 1; return; }
  for (int i = 0; i < N; i++) {
    Distances[i] = (double*)malloc(N * sizeof(double));
    if (Distances[i] == NULL) { *mem_error = 1; return; }
    DistancesT[i] = (double*)malloc(N * sizeof(double));
    if (DistancesT[i] == NULL) { *mem_error = 1; return; }
  }
  
  // Fill Distances and DistancesT with values from input
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Distances[i][j] = distances[i * N + j];
      DistancesT[i][j] = 2 * distances[i * N + j];
    }
  }
  
  // Allocate memory for LB and UB arrays
  LB = (int*)malloc(K * sizeof(int));
  if (LB == NULL) { *mem_error = 1; return; }
  UB = (int*)malloc(K * sizeof(int));
  if (UB == NULL) { *mem_error = 1; return; }
  for (int i = 0; i < K; i++) {
    LB[i] = *lower_bound;  // Assuming lower_bound is a pointer to an int
    UB[i] = *upper_bound;  // Assuming upper_bound is a pointer to an int
  }
  
  AssignMemory();
  if (*mem_error == 1) {
    return;
  }
  
  BuildNeighbors();
  
  SearchAlgorithm();
  
  //save GS -> solution with result
  for (int i = 0; i < N; i++){
    result[i] = GS.p[i];
  }
  *cost = GS.cost;
  
  // Remember to free the allocated memory after use
  for (int i = 0; i < N; i++) {
    free(Distances[i]); Distances[i] = NULL;
    free(DistancesT[i]); DistancesT[i] = NULL;
  }
  free(Distances); Distances = NULL;
  free(DistancesT); DistancesT = NULL;
  free(LB); LB = NULL;
  free(UB); UB = NULL;
  
  ReleaseMemory();
}


void SearchAlgorithm() {
    /* Algorithm 1: The main procedure of TPSDP. */

    int eta;
    int pickedSolution;
    clock_t starting_time = clock();
    GS.cost = -INFINITY;
    
    // Initial population generation
    int i, j, k;
    for (i = 0; i < populationSize; i++) {
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
        
        eta = (int)(theta * N / K);
        for (i = 0; i < populationSize; i++) {
            for (j = 0; j < N; j++) Offs[i].p[j] = Pop[i].p[j];
            for (k = 0; k < K; k++) Offs[i].SizeG[k] = Pop[i].SizeG[k];
            Offs[i].cost = Pop[i].cost;
        }
        // Strong Perturbation and Local Search
        for (i = 0; i < populationSize; i++) {
            StrongPerturbation(eta, Pop[i].p, Pop[i].SizeG);
            RandLS(Pop[i].p, Pop[i].SizeG, &Pop[i].cost);
            if (Pop[i].cost > GS.cost) {
                for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
                for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
                GS.cost = Pop[i].cost;
            }
        }

        // Crossover and Local Search
        if (populationSize > 1) {
            for (i = 0; i < populationSize; i++) {
                pickedSolution = random_int(populationSize);
                do {
                    pickedSolution = (pickedSolution + 1) % populationSize;
                } while (pickedSolution == i);
                Crossover(Pop[i].p, Pop[pickedSolution].p, Offs[i].p, Offs[i].SizeG);
                RandLS(Offs[i].p, Offs[i].SizeG, &Offs[i].cost);
            }
            for (i = 0; i < populationSize; i++) {
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
        for (i = 0; i < populationSize; i++) {
            DirectPerturbation(eta_max, Pop[i].p, Pop[i].SizeG);
            RandLS(Pop[i].p, Pop[i].SizeG, &Pop[i].cost);
            if (Pop[i].cost > GS.cost) {
                for (j = 0; j < N; j++) GS.p[j] = Pop[i].p[j];
                for (k = 0; k < K; k++) GS.SizeG[k] = Pop[i].SizeG[k];
                GS.cost = Pop[i].cost;
            }
        }

        // Linear decrease population size
        // Note: Implement sort function based on the comparison function `Cmpare`
        qsort(Pop, populationSize, sizeof(Solution), Cmpare);
        populationSize = (int)(beta_min - populationSize) / Time_limit * (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC) + populationSize;
        theta = theta_max - (theta_max - theta_min) * (1.0 * (clock() - starting_time) / CLOCKS_PER_SEC) / Time_limit;
    }
}

/* Algorithm 2: initializes a solution S */
void InitialSol(Solution *S) {
    /* Algorithm 2: initializes a solution S */
    RandomInitialSol(S->p, S->SizeG);
    RandLS(S->p, S->SizeG, &(S->cost));
}

int Cmpare(const void *a, const void *b) {
    /* Compares two solutions based on their cost */
    Solution *solA = (Solution *)a;
    Solution *solB = (Solution *)b;
    return (solB->cost - solA->cost); // Return positive if b is greater, negative if a is greater
}

void RandomInitialSol(int p[], int SizeG[]) {
	/* Algorithm 2: initializes a random solution that respects the group size constraints */
	int i, j;

    // Allocates memory
	int* isAssigned = (int *)malloc(N * sizeof(int));  // Tracks if an element is assigned
	int* groupSize =(int *)malloc(K * sizeof(int)); // Stores the size of each group
	
	for (i = 0; i < K; i++) groupSize[i] = 0;
	for (i = 0; i < N; i++) isAssigned[i] = 0;

	
    // Calculate the total number of elements that need to satisfy the lower bounds
    int total_assigned = 0;
    int total_LB = 0;
    for (int i = 0; i < K; i++) {
        total_LB += LB[i];
    }

    // First phase: Assign elements to satisfy lower bound constraints (LB)
    while (total_assigned < total_LB) {
        int selected_element = random_int(N);

        if (!isAssigned[selected_element]) {  // Only assign unassigned elements
            for (int group = 0; group < K; group++) {
                if (groupSize[group] < LB[group]) {
                    p[selected_element] = group;
                    isAssigned[selected_element] = 1;
                    groupSize[group]++;
                    total_assigned++;
                    break;  // Move to the next element once assigned
                }
            }
        }
    }
	
	// Second phase: Assign remaining elements, respecting the upper bound (UB)
    while (total_assigned < N) {
        int selected_element = random_int(N);
        if (!isAssigned[selected_element]) {
            int group;
            do {
                group = random_int(K);  // Randomly select a group
            } while (groupSize[group] >= UB[group]);  // Ensure group doesn't exceed UB

            p[selected_element] = group;
            isAssigned[selected_element] = 1;
            groupSize[group]++;
            total_assigned++;
        }
    }

    // Copy the final group sizes into the output array SizeG
   	for (i = 0; i < K; i++)  SizeG[i] = groupSize[i];

    // Free allocated memory
	free(groupSize); groupSize = NULL;
	free(isAssigned); isAssigned = NULL;
}

void RandLS(int partition[], int SizeGroup[], double* cost) {
    const double DELTA_THRESHOLD = 0.0001;  // Define a constant for comparison threshold
    int i, vertex, group, otherGroup, x, y;
    int oldGroup, oldGroup1, tempSwap;
    double delta;
    int improved;

    // Initialize the partition array
    for (i = 0; i < N; i++) p[i] = partition[i];

    // Build the delta matrix for cost changes
    BuildDeltaMatrix();

    // Initialize the delta value
    delta = -99999.0;

    do {
        improved = 0;  // Reset improvement flag

        // First loop: Move individual elements to improve partition
        for (vertex = 0; vertex < N; vertex++) {
            for (group = 0; group < K; group++) {
                // Check if moving `vertex` to `group` is valid and beneficial
                if ((p[vertex] != group) && (SizeGroup[p[vertex]] > LB[p[vertex]]) && (SizeGroup[group] < UB[group])) {
                    delta = Delta_Matrix[vertex][group] - Delta_Matrix[vertex][p[vertex]];

                    if (delta > DELTA_THRESHOLD) {
                        oldGroup = p[vertex];

                        // Update delta matrix for the move
                        OneMoveUpdateDeltaMatrix1(vertex, oldGroup, group);

                        // Update group sizes
                        SizeGroup[oldGroup] -= 1;
                        SizeGroup[group] += 1;

                        // Assign vertex to new group
                        p[vertex] = group;

                        // Update total cost
                        *cost += delta;

                        // Mark as improved
                        improved = 1;
                    }
                }
            }
        }

        // Second loop: Swap pairs of elements between groups
        for (x = 0; x < N; x++) {
            for (y = x + 1; y < N; y++) {
                // Only swap if nodes are in different groups
                if (p[x] != p[y]) {
                    delta = (Delta_Matrix[x][p[y]] - Delta_Matrix[x][p[x]])
                          + (Delta_Matrix[y][p[x]] - Delta_Matrix[y][p[y]])
                          - DistancesT[x][y];

                    if (delta > DELTA_THRESHOLD) {
                        oldGroup = p[x];
                        oldGroup1 = p[y];

                        // Update delta matrix for the swap
                        OneMoveUpdateDeltaMatrix1(x, oldGroup, oldGroup1);
                        OneMoveUpdateDeltaMatrix1(y, oldGroup1, oldGroup);

                        // Swap the two nodes between groups
                        tempSwap = p[x];
                        p[x] = p[y];
                        p[y] = tempSwap;

                        // Update total cost
                        *cost += delta;

                        // Mark as improved
                        improved = 1;
                    }
                }
            }
        }
    } while (improved == 1);  // Continue until no improvement is made

    // Update the partition array with the final assignments
    for (i = 0; i < N; i++) partition[i] = p[i];
}

void StrongPerturbation(int L, int partition[], int SizeGroup[]) {
    /* Algorithm 4: Undirected Perturbation. Applies a strong perturbation to the partition */

    int NumberNeighbors = N * (N - 1) / 2 + N * K;
    int count = 0;
    int cur_index;
    int v, g, x, y;
    int old_g, old_g1, swap;
    double delt;

    // Copy the current partition
    int* temp_partition = (int*)malloc(N * sizeof(int));
    if (temp_partition == NULL) {
        fprintf(stderr, "Memory allocation failed for temp_partition.\n");
        return;
    }
    for (int i = 0; i < N; i++) {
        temp_partition[i] = partition[i];
    }

    theta = L;
    // Perturbation loop
    while (count < theta) {
        cur_index = random_int(NumberNeighbors);

        if (Neighbors[cur_index].type == 1) { // Type 1 neighbor: (element, group)
            v = Neighbors[cur_index].v;
            g = Neighbors[cur_index].g;

            // Apply perturbation if constraints are met
            if (temp_partition[v] != g &&
                SizeGroup[temp_partition[v]] > LB[temp_partition[v]] &&
                SizeGroup[g] < UB[g]) {

                delt = Delta_Matrix[v][g] - Delta_Matrix[v][temp_partition[v]];
                old_g = temp_partition[v];
                SizeGroup[old_g]--;
                SizeGroup[g]++;
                temp_partition[v] = g;
                count++;
            }
        } else if (Neighbors[cur_index].type == 2) { // Type 2 neighbor: (element x, element y)
            x = Neighbors[cur_index].x;
            y = Neighbors[cur_index].y;

            // Apply perturbation if elements are in different groups
            if (temp_partition[x] != temp_partition[y]) {
                delt = (Delta_Matrix[x][temp_partition[y]] - Delta_Matrix[x][temp_partition[x]]) +
                       (Delta_Matrix[y][temp_partition[x]] - Delta_Matrix[y][temp_partition[y]]) -
                       2 * Distances[x][y];

                old_g = temp_partition[x];
                old_g1 = temp_partition[y];
                swap = temp_partition[x];
                temp_partition[x] = temp_partition[y];
                temp_partition[y] = swap;
                count++;
            }
        }
    }

    // Copy the perturbed partition back to the original partition
    for (int i = 0; i < N; i++) {
        partition[i] = temp_partition[i];
    }

    free(temp_partition); // Free the allocated memory
}

void DirectPerturbation(int eta_max, int partition[], int SizeGroup[]) {
    /* Algorithm 6: Directed Perturbation. 
	Iteratively refines partitions to balance group sizes and minimize costs */

    int i, j, k, eta, number, minDeltaValue, minElement;
    double delt, delt1, delt_max;
    int NumberNeighbors, old_g, old_g1, swap;

    // Initialize the partition and size groups
    for (i = 0; i < N; i++) p[i] = partition[i];
    for (j = 0; j < K; j++) SizeG[j] = SizeGroup[j];
    BuildDeltaMatrix();

    // Main loop for perturbation iterations
    for (eta = 0; eta < eta_max; eta++) {

        // Reset tracking variables
        number = 0;
        for (i = 0; i < K; i++) {
            UnderLB[i] = 0;
            Rd[i] = -1;
            for (j = 0; j < K; j++) {
                AvgCon[i][j] = 0.0;
            }
        }

        // Find the minimum scoring element for each group
        for (k = 0; k < K; k++) {
            minDeltaValue = 99999999;
            minElement = 0;
            for (i = 0; i < N; i++) {
                if (p[i] == k) {
                    if (Delta_Matrix[i][k] < minDeltaValue) {
                        minDeltaValue = Delta_Matrix[i][k];
                        minElement = i;
                    }
                }
            }

            // Record the minimum element for removal
            Rd[k] = minElement;
            SizeG[k] -= 1;

            // If the group size falls below the lower bound, mark it
            if (SizeG[k] < LB[k]) {
                UnderLB[k] = 1;
                number += 1;
            }
        }

        // Rebuild the Delta matrix and average connections after removal
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                Delta_Matrix[Rd[i]][p[Rd[j]]] = Delta_Matrix[Rd[i]][p[Rd[j]]] - Distances[Rd[i]][Rd[j]];
                AvgCon[p[Rd[i]]][p[Rd[j]]] = Delta_Matrix[Rd[i]][p[Rd[j]]] / SizeG[p[Rd[j]]];
            }
        }

        // Handle groups that are under the lower bound (LB)
        int selectedGroup;
        int maxAvgCon;
        int nn = 0;
        while (nn < number) {
            maxAvgCon = -9999;
            i = random_int(K);

            // Find the element with the highest average connection to the group
            do {
                i = (i + 1) % K;
            } while (UnderLB[i] == 0);
            for (j = 0; j < K; j++) {
                if (AvgCon[j][i] > maxAvgCon) {
                    maxAvgCon = AvgCon[j][i];
                    selectedGroup = j;
                }
            }

            // Move the selected element to the group
            SizeG[i] += 1;
            for (k = 0; k < K; k++) {
                if (Rd[k] != -1) {
                    Delta_Matrix[Rd[k]][i] += Distances[Rd[k]][Rd[selectedGroup]];
                    AvgCon[p[Rd[k]]][i] = Delta_Matrix[Rd[k]][i] / SizeG[i];
                }
            }

            // Clear old connections for the moved element and finalize the move
            for (k = 0; k < K; k++) {
                AvgCon[p[Rd[selectedGroup]]][k] = 0.0;
            }
            p[Rd[selectedGroup]] = i;
            UnderLB[i] = 0;
            Rd[selectedGroup] = -1;
            nn++;
        }

        // Handle remaining elements for groups that are above LB
        int groupWithMaxAvgCon;
        nn = 0;
        while (nn < K - number) {
            selectedGroup = random_int(K);
             // Find the group with the highest average connection for the element
            do {
                selectedGroup = (selectedGroup + 1) % K;
            } while (Rd[selectedGroup] == -1);
            maxAvgCon = -9999;
            for (j = 0; j < K; j++) {
                if (AvgCon[j][selectedGroup] > maxAvgCon) {
                    maxAvgCon = AvgCon[j][selectedGroup];
                    groupWithMaxAvgCon = j;
                }
            }
            // Move the selected element to the group
            SizeG[selectedGroup] += 1;
            for (k = 0; k < K; k++) {
                if (Rd[k] != -1) {
                    Delta_Matrix[Rd[k]][selectedGroup] += Distances[Rd[k]][Rd[groupWithMaxAvgCon]];
                    AvgCon[p[Rd[k]]][selectedGroup] = Delta_Matrix[Rd[k]][selectedGroup] / SizeG[selectedGroup];
                }
            }

            // Clear old connections for the moved element and finalize the move
            for (k = 0; k < K; k++) {
                AvgCon[p[Rd[groupWithMaxAvgCon]]][k] = 0.0;
            }
            p[Rd[groupWithMaxAvgCon]] = selectedGroup;
            Rd[groupWithMaxAvgCon] = -1;
            nn++;
        }
    }

    for (i = 0; i < N; i++) partition[i] = p[i];
}

void Crossover(int partition1[], int partition2[], int score[], int scSizeGroup[]) {
    /* Algorithm 5: combines partitions in a way that maintains group constraints */

    int i, j, k, maxGroupDiversity, selectedGroup;
    int elementCount, groupCount;
    int unassignedFlag;
    int surplus;
    int targetGroup = -1;
    int processedCount;
    int selectedElement;
    int totalAssigned, totalLowerBound, totalBelowLowerBound;

    // Initialize p and p1 with partition1
    for (i = 0; i < N; i++) {
        p[i] = partition1[i];
        p1[i] = partition1[i];
    }
    BuildDeltaMatrix();
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            Delta_Matrix_p1[i][j] = Delta_Matrix[i][j];
        }
    }
    BuildGroupDiversityForCrossover();
    for (i = 0; i < K; i++) {
        groupDiversity_p1[i] = groupDiversity[i];
    }

    // Initialize p and p2 with partition2
    for (i = 0; i < N; i++) {
        p[i] = partition2[i];
        p2[i] = partition2[i];
    }
    BuildDeltaMatrix();
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            Delta_Matrix_p2[i][j] = Delta_Matrix[i][j];
        }
    }
    BuildGroupDiversityForCrossover();
    for (i = 0; i < K; i++) {
        groupDiversity_p2[i] = groupDiversity[i];
    }

    // Initialize arrays
    for (i = 0; i < N; i++) {
        vectorElement[i] = i;
        score[i] = -1;
    }
    for (i = 0; i < K; i++) {
        LBGroup[i] = 0;
        UBGroup[i] = 0;
        BigThanLB[i] = 0;
        groupElement[i] = i;
        ub[i] = UB[i];
        scSizeGroup[i] = 0;
    }

    // Main crossover process
    for (i = 0; i < K; i++) {
        if (uniform_rnd_number() < 0.5) {
            // Process partition1
            maxGroupDiversity = -9999;
            for (j = 0; j < K; j++) {
                if (groupDiversity_p1[j] > maxGroupDiversity) {
                    maxGroupDiversity = groupDiversity_p1[j];
                    selectedGroup = j;
                }
            }

            elementCount = 0;
            for (j = 0; j < N; j++) {
                if (p1[j] == selectedGroup) {
                    SelectEle[elementCount++] = j;
                }
            }

            groupCount = 0;
            for (j = 0; j < K; j++) {
                if (ub[j] != -1 && ub[j] >= elementCount) {
                    SelectGroup[groupCount++] = j;
                }
            }

            if (groupCount == 0) { // No valid group found
                int minDiff = 999999;
                for (j = 0; j < K; j++) {
                    if (ub[j] != -1 && elementCount - ub[j] < minDiff) {
                        minDiff = elementCount - ub[j];
                        targetGroup = j;
                    }
                }

                processedCount = 0;
                while (processedCount < elementCount - minDiff) {
                    selectedElement = random_int(elementCount);
                    do {
                        selectedElement = (selectedElement + 1) % elementCount;
                    } while (SelectEle[selectedElement] == -1);

                    score[SelectEle[selectedElement]] = targetGroup;
                    SelectEleTemp[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];
                for (j = 0; j < elementCount; j++) {
                    score[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        } else {
            // Process partition2 (similar to partition1 logic)
            maxGroupDiversity = -9999;
            for (j = 0; j < K; j++) {
                if (groupDiversity_p2[j] > maxGroupDiversity) {
                    maxGroupDiversity = groupDiversity_p2[j];
                    selectedGroup = j;
                }
            }

            elementCount = 0;
            for (j = 0; j < N; j++) {
                if (p2[j] == selectedGroup) {
                    SelectEle[elementCount++] = j;
                }
            }

            groupCount = 0;
            for (j = 0; j < K; j++) {
                if (ub[j] != -1 && ub[j] >= elementCount) {
                    SelectGroup[groupCount++] = j;
                }
            }

            if (groupCount == 0) { // No valid group found
                int minDiff = 999999;
                for (j = 0; j < K; j++) {
                    if (ub[j] != -1 && elementCount - ub[j] < minDiff) {
                        minDiff = elementCount - ub[j];
                        targetGroup = j;
                    }
                }

                processedCount = 0;
                while (processedCount < elementCount - minDiff) {
                    selectedElement = random_int(elementCount);
                    do {
                        selectedElement = (selectedElement + 1) % elementCount;
                    } while (SelectEle[selectedElement] == -1);

                    score[SelectEle[selectedElement]] = targetGroup;
                    SelectEleTemp[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];
                for (j = 0; j < elementCount; j++) {
                    score[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        }

        // Update group diversity
        for (j = 0; j < elementCount; j++) {
            groupDiversity_p1[p1[SelectEleTemp[j]]] -= Delta_Matrix_p1[SelectEleTemp[j]][p1[SelectEleTemp[j]]];
            groupDiversity_p2[p2[SelectEleTemp[j]]] -= Delta_Matrix_p2[SelectEleTemp[j]][p2[SelectEleTemp[j]]];
            p1[SelectEleTemp[j]] = -1;
            p2[SelectEleTemp[j]] = -1;
        }

        ub[targetGroup] = -1;
        scSizeGroup[targetGroup] = elementCount;
    }

    // Adjust assignments to maintain group size constraints
    processedCount = 0;
    totalLowerBound = 0;
    totalBelowLowerBound = 0;
    for (i = 0; i < K; i++) {
        totalLowerBound += LB[i];
        if (scSizeGroup[i] < LB[i]) {
            processedCount += scSizeGroup[i];
            totalBelowLowerBound += scSizeGroup[i];
            LBGroup[i] = 1;
        } else {
            processedCount += LB[i];
        }
        if (scSizeGroup[i] > LB[i]) {
            BigThanLB[i] = 1;
        }
    }

    // Assign unprocessed elements
    for (i = 0; i < N; i++) {
        if (vectorElement[i] != -1) {
            processedCount++;
        }
    }

    while (processedCount < totalLowerBound) {
        targetGroup = random_int(K);
        do {
            targetGroup = (targetGroup + 1) % K;
        } while (BigThanLB[targetGroup] == 0);

        elementCount = 0;
        for (j = 0; j < N; j++) {
            if (score[j] == targetGroup) {
                SelectEle[elementCount++] = j;
            }
        }

        selectedElement = random_int(elementCount);
        score[SelectEle[selectedElement]] = -1;
        vectorElement[SelectEle[selectedElement]] = SelectEle[selectedElement];
        scSizeGroup[targetGroup]--;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            BigThanLB[targetGroup] = 0;
        }
        processedCount++;
    }

    int sumLB = 0;
    for (i = 0; i < K; i++) {
        if (LBGroup[i] == 1) {
            sumLB += LB[i];
        }
    }

    while (totalBelowLowerBound < sumLB) {
        targetGroup = random_int(K);
        do {
            targetGroup = (targetGroup + 1) % K;
        } while (LBGroup[targetGroup] == 0);

        elementCount = 0;
        for (i = 0; i < N; i++) {
            if (vectorElement[i] != -1) {
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = random_int(elementCount);
        score[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            LBGroup[targetGroup] = 0;
        }
        totalBelowLowerBound++;
    }

    int totalSize = 0;
    for (i = 0; i < K; i++) {
        totalSize += scSizeGroup[i];
        if (scSizeGroup[i] < UB[i]) {
            UBGroup[i] = 1;
        }
    }

    while (totalSize < N) {
        targetGroup = random_int(K);
        do {
            targetGroup = (targetGroup + 1) % K;
        } while (UBGroup[targetGroup] == 0);

        elementCount = 0;
        for (i = 0; i < N; i++) {
            if (vectorElement[i] != -1) {
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = random_int(elementCount);
        score[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == UB[targetGroup]) {
            UBGroup[targetGroup] = 0;
        }
        totalSize++;
    }
}

double FitRadioAndDis(int partition1[], int partition2[], double cost1, double cost2) {
    /* 
     * Evaluates the quality and dissimilarity of partitions.
     * It calculates a 'radio' value that combines the ratio of costs and a
     * dissimilarity factor between partition1 and partition2.
     */
    
    // Handle potential division by zero
    if (cost2 == 0.0) {
        fprintf(stderr, "Error: Division by zero (cost2 is zero).\n");
        return -1;  
    }

    int i, j;
    int count = 0;
    int totalPairs = (N * (N - 1)) / 2;  // Number of unique pairs (i, j) where i < j

    // Loop over all pairs of elements to count dissimilar pairs
    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            int sameInPartition1 = (partition1[i] == partition1[j]);
            int sameInPartition2 = (partition2[i] == partition2[j]);

            // Count cases where elements are grouped differently in the two partitions
            if (sameInPartition1 != sameInPartition2) {
                count++;
            }
        }
    }

    // Calculate the radio value: ratio of costs + dissimilarity factor
    double dissimilarityFactor = 0.05 * ((double)count / totalPairs) * K;
    double radio = cost1 / cost2 + dissimilarityFactor;

    return radio;
}

void BuildNeighbors() {
	/* Initializes the neighbor structure for optimization */
	int i, j;
	int count = 0;
	
    // Type 1 neighbors: (i, j) where each element i can be in group j
	for (i = 0; i < N; i++)
		for (j = 0; j < K; j++)
		{
			Neighbors[count].type = 1;
			Neighbors[count].v = i;
			Neighbors[count].g = j;
			count++;
		}

    // Type 2 neighbors: (i, j) where each pair of elements (i, j) are neighbors    
	for (i = 0; i < N; i++)
		for (j = i + 1; j < N; j++)
		{
			Neighbors[count].type = 2;
			Neighbors[count].x = i;
			Neighbors[count].y = j;
			count++;
		}
		
}

void ClearDeltaMatrix() {
	/* Resets the delta matrix and objective function value */
    objective_function = 0.0;
    for (int x = 0; x < N; ++x) {
        for (int g = 0; g < K; ++g) {
            Delta_Matrix[x][g] = 0.0;
        }
    }
}

void BuildDeltaMatrix() {
	/*  Builds the delta matrix and calculates the objective function value */
	int i, j;
	ClearDeltaMatrix();

	// Update Delta_Matrix based on distances
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Delta_Matrix[i][p[j]] += Distances[i][j];
        }
    }

    // Compute Delta values based on Delta_Matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < K; ++j) {
            Delta[i][j] = Delta_Matrix[i][j] - Delta_Matrix[i][p[i]];
        }
    }

    // Calculate the objective function value
    objective_function = 0.0;
    for (int i = 0; i < N; ++i) {
        objective_function += Delta_Matrix[i][p[i]];
    }
    objective_function /= 2.0;
}

void BuildGroupDiversityForCrossover() {
	/*  Builds group diversity values for crossover */
	
    // Initialize group diversity values to zero
	for (int i = 0; i < K; i++) groupDiversity[i] = 0.0;
	
	// Compute group diversity based on distances
    for (int i = 0; i < N; ++i) {
        int group_i = p[i];
        for (int j = 0; j < N; ++j) {
            if (group_i == p[j]) {
                groupDiversity[group_i] += Distances[i][j];
            }
        }
    }
}

void OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup) {
	/* udates the delta matrix for single-element moves */
	int x, j, k;

    // Update Delta_Matrix for all elements affected by the move
	for (j = 0; j < N; j++) {
		if (j != i) {
			Delta_Matrix[j][oldGroup] -= Distances[i][j];
			Delta_Matrix[j][newGroup] += Distances[i][j];
		}
	}

    // Recompute Delta values for all elements except the moved element
	for (x = 0; x < N; x++) {
		if (x != i) {
			Delta[x][oldGroup] = Delta_Matrix[x][oldGroup] - Delta_Matrix[x][p[x]];
			Delta[x][newGroup] = Delta_Matrix[x][newGroup] - Delta_Matrix[x][p[x]];

			if (p[x] == oldGroup) {
				for (k = 0; k < K; k++) {
					Delta[x][k] = Delta_Matrix[x][k] - Delta_Matrix[x][oldGroup];
				}
			}

			if (p[x] == newGroup) {
				for (k = 0; k < K; k++) {
					Delta[x][k] = Delta_Matrix[x][k] - Delta_Matrix[x][newGroup];
				}
			}

		}
	}

    // Update Delta values for the moved element
	x = i; 
	for (k = 0; k < K; k++) Delta[x][k] = Delta_Matrix[x][k] - Delta_Matrix[x][newGroup];
}

void OneMoveUpdateDeltaMatrix1(int i, int oldGroup, int newGroup) {
	for (int j = 0; j < N; j++) {
		if (j != i) {
			Delta_Matrix[j][oldGroup] -= Distances[i][j];
			Delta_Matrix[j][newGroup] += Distances[i][j];
		}
	}
}

void AssignMemory() {
    /*  Allocates memory dynamically for various arrays and matrices necessary 
	for the algorithm's execution. This includes structures for population management, 
	distance matrices, diversity measures, and neighborhood exploration.
	*/

    int i, j;
    
    p = (int*)malloc(N * sizeof(int));
    bestp = (int*)malloc(N * sizeof(int));
    SizeG = (int*)malloc(K * sizeof(int));
    
    Pop = (Solution*)malloc(populationSize * sizeof(Solution));
    Offs = (Solution*)malloc(populationSize * sizeof(Solution));
    
    Delta_Matrix = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix[i] = (double*)malloc(K * sizeof(double));
    Delta_Matrix_p1 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p1[i] = (double*)malloc(K * sizeof(double));
    Delta_Matrix_p2 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p2[i] = (double*)malloc(K * sizeof(double));
    Delta_Matrix1 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix1[i] = (double*)malloc(K * sizeof(double));
    groupDiversity = (double*)malloc(K * sizeof(double));
    groupDiversity_p1 = (double*)malloc(K * sizeof(double));
    groupDiversity_p2 = (double*)malloc(K * sizeof(double));
    
    Delta = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta[i] = (double*)malloc(K * sizeof(double));
    
    Delta1 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta1[i] = (double*)malloc(K * sizeof(double));
    
    for (i = 0; i < populationSize; i++) {
        Pop[i].p = (int*)malloc(N * sizeof(int));
        Offs[i].p = (int*)malloc(N * sizeof(int));
    }
    
    for (i = 0; i < populationSize; i++) {
        Pop[i].SizeG = (int*)malloc(K * sizeof(int));
        Offs[i].SizeG = (int*)malloc(K * sizeof(int));
    }
    
    CS.p = (int*)malloc(N * sizeof(int));
    NS.p = (int*)malloc(N * sizeof(int));
    GS.p = (int*)malloc(N * sizeof(int));
    OS.p = (int*)malloc(N * sizeof(int));
    
    CS.SizeG = (int*)malloc(K * sizeof(int));
    NS.SizeG = (int*)malloc(K * sizeof(int));
    GS.SizeG = (int*)malloc(K * sizeof(int));
    OS.SizeG = (int*)malloc(K * sizeof(int));
    
    Neighbors = (Neighborhood*)malloc((N * (N - 1) / 2 + N * K) * sizeof(Neighborhood));
    
    AvgCon = (double**)malloc(K * sizeof(double*));
    for (i = 0; i < K; i++) AvgCon[i] = (double*)malloc(K * sizeof(double));
    Rd = (int*)malloc(K * sizeof(int));
    for (i = 0; i < K; i++) Rd[i] = 0;
    UnderLB = (int*)malloc(K * sizeof(int));
    
    ub = (int*)malloc(K * sizeof(int));
    LBGroup = (int*)malloc(K * sizeof(int));
    UBGroup = (int*)malloc(K * sizeof(int));
    BigThanLB = (int*)malloc(K * sizeof(int));
    vectorElement = (int*)malloc(N * sizeof(int));
    groupElement = (int*)malloc(K * sizeof(int));
    SelectEle = (int*)malloc(N * sizeof(int));
    SelectGroup = (int*)malloc(K * sizeof(int));
    SelectEleTemp = (int*)malloc(N * sizeof(int));
    p1 = (int*)malloc(N * sizeof(int));
    p2 = (int*)malloc(N * sizeof(int));
}

void ReleaseMemory() {
    /* responsible for reading the input file, 
    initializing matrices, and setting constraints on group sizes. */ 
    int i;
    
    free(p); p = NULL;
    free(bestp); bestp = NULL;
    free(SizeG); SizeG = NULL;
    
    free(CS.p); CS.p = NULL;
    free(CS.SizeG); CS.SizeG = NULL;
    free(GS.p); GS.p = NULL;
    free(GS.SizeG); GS.SizeG = NULL;
    free(NS.p); NS.p = NULL;
    free(NS.SizeG); NS.SizeG = NULL;
    free(OS.p); OS.p = NULL;
    free(OS.SizeG); OS.SizeG = NULL;
    
    free(LB); LB = NULL;
    free(UB); UB = NULL;
    free(Neighbors); Neighbors = NULL;
    
    for (i = 0; i < N; i++) {
        free(Delta_Matrix[i]); Delta_Matrix[i] = NULL;
        free(Delta_Matrix1[i]); Delta_Matrix1[i] = NULL;
        free(Delta_Matrix_p1[i]); Delta_Matrix_p1[i] = NULL;
        free(Delta_Matrix_p2[i]); Delta_Matrix_p2[i] = NULL;
        free(Delta[i]); Delta[i] = NULL;
        free(Delta1[i]); Delta1[i] = NULL;
    }
    free(Delta_Matrix); Delta_Matrix = NULL;
    free(Delta_Matrix_p1); Delta_Matrix_p1 = NULL;
    free(Delta_Matrix_p2); Delta_Matrix_p2 = NULL;
    free(Delta_Matrix1); Delta_Matrix1 = NULL;
    free(Delta); Delta = NULL;
    free(Delta1); Delta1 = NULL;
    free(groupDiversity); groupDiversity = NULL;
    free(groupDiversity_p1); groupDiversity_p1 = NULL;
    free(groupDiversity_p2); groupDiversity_p2 = NULL;
    free(AvgCon); AvgCon = NULL;
    free(Rd); Rd = NULL;
    free(UnderLB); UnderLB = NULL;
    free(ub); ub = NULL;
    free(LBGroup); LBGroup = NULL;
    free(UBGroup); UBGroup = NULL;
    free(BigThanLB); BigThanLB = NULL;
    free(vectorElement); vectorElement = NULL;
    free(groupElement); groupElement = NULL;
    free(SelectEle); SelectEle = NULL;
    free(SelectGroup); SelectGroup = NULL;
    free(SelectEleTemp); SelectEleTemp = NULL;
    free(p1); p1 = NULL;
    free(p2); p2 = NULL;
}

// Generate a random integer from zero to max-1
int random_int(int max) {
GetRNGstate();
  double my_number = unif_rand();
  PutRNGstate();
  return (int) floor(my_number * max);
}
