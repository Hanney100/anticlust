#ifndef THREE_PHASE_HEADER_H 
#define THREE_PHASE_HEADER_H

typedef struct {
    int *s;        // cluster membership of each element
    int *SizeG;    // size of each cluster
    double cost;   // global cost function value of solution
} Solution;

int random_int(int max);
double uniform_rnd_number(void);
void swap_elements(int* a, int* b);
void fisher_yates_shuffle(int arr[], int n);
double LocalSearchCriterionCalcutlation(int partition1[], int partition2[], double cost1, double cost2);
void RandomInitialSol(int s[], int SizeG[]);
int Cmpare(const void *a, const void *b);

#endif // THREE_PHASE_HEADER_H