#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);
double cheapInsertion(double **distanceMatrix, int *besttour, int startVertex, int N);
double farInsertion(double **distanceMatrix, int *besttour, int startVertex, int N);
float nearAddition(double **distanceMatrix, int *besttour, int startVertex,int N);

struct Tour {
  double cost;
  int tour[1000];
};

double **calculateDistance(double **coordinates, int verticesCount) {
    double (*distanceMatrix)[verticesCount] = malloc(verticesCount * verticesCount * sizeof(double));
    if (distanceMatrix == NULL) {
        perror("Memory Allocation Failed for distanceMatrix");
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < verticesCount; i++) {
        for (int j = 0; j < verticesCount; j++) {
            distanceMatrix[i][j] = sqrt(((coordinates[i][0] - coordinates[j][0]) * (coordinates[i][0] - coordinates[j][0])) +
                                        ((coordinates[i][1] - coordinates[j][1]) * (coordinates[i][1] - coordinates[j][1])));
        }
    }

    double **distanceMatrixPtrs = malloc(verticesCount * sizeof(double*));
    if (distanceMatrixPtrs == NULL) {
        perror("Memory Allocation Failed for distanceMatrixPtrs");
        free(distanceMatrix);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < verticesCount; i++) {
        distanceMatrixPtrs[i] = distanceMatrix[i];
    }

    return distanceMatrixPtrs;
}

int main(char argc, char *argv[]){

    #ifdef USE_MPI
    MPI_Init(NULL, NULL);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #else
    int size = 1;
    int rank = 0;
    #endif


  
    char *fileName = "9_coords.coord";
#ifdef USE_MPI


    double TimeStart , TimeEnd;

    int N = readNumOfCoords(fileName);

    double  **coords_2od_ary = readCoords(fileName, N);
    double **distanceMatrix= calculate_distance(coords_2od_ary, N);


    double localBCCI = DBL_MAX, localBCFI = DBL_MAX, localBCNA = DBL_MAX;
    int localBTCI[N + 1], localBTFI[N + 1], localBTNA[N + 1];
    int tour[N + 1], tourFI[N+1], tourNA[N+1];

    // Allocate memory to gather all best tours and costs on rank 0
    double *allBCCI = NULL, *allBCFI = NULL, *allBCNA = NULL;
    // Array of arrays (each element is an int[N + 1])
    int (*aBTCI)[N + 1] = NULL;
    int (*aBTFI)[N + 1] = NULL;
    int (*aBTNA)[N + 1] = NULL;


    //Initialize all the global best costs and tours
    if (rank == 0) {
        allBCCI = malloc(size * sizeof(double));
        aBTCI = malloc(size * sizeof(*aBTCI));
        if (allBCCI == NULL || aBTCI == NULL) {
            perror("Memory Allocation Failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        allBCFI = malloc(size * sizeof(double));
        aBTFI = malloc(size * sizeof(*aBTFI));
        if (allBCFI == NULL || aBTFI == NULL) {
            perror("Memory Allocation Failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        allBCNA = malloc(size * sizeof(double));
        aBTNA = malloc(size * sizeof(*aBTNA));
        if (allBCNA == NULL || aBTNA == NULL) {
            perror("Memory Allocation Failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

    }

    // Start the timer
    TimeStart = MPI_Wtime();


    // Loop over each vertex as the starting point, distributed across ranks
    for (int i = rank; i < N; i += size) {

        // Cheapest Insertion TSP
        double costCI = cheapestInsertionTSP(distanceMatrix, tour, i, N);
        if (costCI < localBCCI) {
            localBCCI = costCI;
            memcpy(localBTCI, tour, (N + 1) * sizeof(int));
        }

        // Farthest Insertion TSP
       double costFI = farthestInsertionTSP(distanceMatrix, tourFI, i, N);
        if (costFI < localBCFI){
            localBCFI = costFI;
            memcpy(localBTFI, tourFI, (N + 1) * sizeof(int));
}
        // Nearest Addition TSP
        double costNA = nearestAdditionTSP(distanceMatrix, tour, i, N);
        if (costNA < localBCNA ) {
            localBCNA = costNA;
            memcpy(localBTNA, tour, (N + 1) * sizeof(int));
        }

}


    // End the timer
    TimeEnd = MPI_Wtime();

    // Gather all local best costs and tours at rank 0
    MPI_Gather(&localBCCI, 1, MPI_DOUBLE, allBCCI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(localBTCI, N + 1, MPI_INT, aBTCI, N + 1, MPI_INT, 0, MPI_COMM_WORLD);


    MPI_Gather(&localBCFI, 1, MPI_DOUBLE, allBCFI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(localBTFI, N + 1, MPI_INT, aBTFI, N + 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gather(&localBCNA, 1, MPI_DOUBLE, allBCNA, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(localBTNA, N + 1, MPI_INT, aBTNA, N + 1, MPI_INT, 0, MPI_COMM_WORLD);



 if (rank == 0) {

        // Initialize the global best tour and cost
        double globalBestCostCI = DBL_MAX;
        int globalBestTourCI[N + 1];

        for (int i = 0; i < size; i++) {
            if (allBCCI[i] < globalBestCostCI ) {
                globalBestCostCI = allBCCI[i];
                memcpy(globalBestTourCI, aBTCI[i], (N + 1) * sizeof(int));
            }
        }

        double globalBestCostFI = DBL_MAX;
        int globalBestTourFI[N + 1];

        for (int i = 0; i < size; i++) {
            if (allBCFI[i] < globalBestCostFI) {
                globalBestCostFI = allBCFI[i];
                memcpy(globalBestTourFI, aBTFI[i], (N + 1) * sizeof(int));
            }
        }



        double globalBestCostNA = DBL_MAX;
        int globalBestTourNA[N + 1];

        for (int i = 0; i < size; i++) {
            if (allBCNA[i] < globalBestCostNA){
                globalBestCostNA = allBCNA[i];
                memcpy(globalBestTourNA, aBTNA[i], (N + 1) * sizeof(int));
            }
        }


       // Write the global best tours to files
       writeTourToFile(globalBestTourCI, N + 1, "best_ci.dat");
       writeTourToFile(globalBestTourFI, N + 1, "best_fi.dat");
       writeTourToFile(globalBestTourNA, N + 1, "best_na.dat");

        // Free the memory allocated for all best tours and costs
        free(allBCCI);
        free(aBTCI);

        free(allBCFI);
        free(aBTFI);

        free(allBCNA);
        free(aBTNA);


     }

    free(distanceMatrix);
    free(coords_2od_ary);
  MPI_Finalize();
    #endif

}


