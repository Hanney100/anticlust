#ifndef THREE_PHASE_HEADER_H 
#define THREE_PHASE_HEADER_H

typedef struct {
    int *s;        // cluster membership of each element
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

int random_int(int max);
double uniform_rnd_number(void);
void swap_elements(int* a, int* b);
void fisher_yates_shuffle(int arr[], int n);
void ClearDeltaMatrix();
void BuildDeltaMatrix();
void OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup);
void BuildGroupDiversityForCrossover();
void AssignMemory();
void ReleaseMemory();
void BuildNeighbors();
void SearchAlgorithm();
void InitialSol(Solution *S);
void UndirectedPerturbation(int L, int partition[], int SizeGroup[]);
void DoubleNeighborhoodLocalSearch(int partition[], int SizeGroup[], double* cost);
void Crossover(int partition1[], int partition2[], int score[], int scSizeGroup[]);
double LocalSearchCriterionCalcutlation(int partition1[], int partition2[], double cost1, double cost2);
void RandomInitialSol(int s[], int SizeG[]);
void DirectPerturbation(int eta_max, int partition[], int SizeGroup[]);
int Cmpare(const void *a, const void *b);

#endif // THREE_PHASE_HEADER_H