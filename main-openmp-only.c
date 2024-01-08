#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <string.h>

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);
double cheapInsertion(double **distanceMatrix, int *besttour, int startVertex, int N);
double farInsertion(double **distanceMatrix, int *besttour, int startVertex, int N);
float nearAddition(double **distanceMatrix, int *besttour, int startVertex, int N);

struct Tour {
  double cost;
  int tour[1000];
};

double **calculateDistance(double **coordinates, int verticesCount) {
    double (*distanceMatrix)[verticesCount] = malloc(verticesCount * verticesCount * sizeof(double));
    if (!distanceMatrix) {
        perror("Memory Allocation Failed");
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
    if (!distanceMatrixPtrs) {
        perror("Memory Allocation Failed for distanceMatrixPtrs");
        free(distanceMatrix);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < verticesCount; i++) {
        distanceMatrixPtrs[i] = distanceMatrix[i];
    }

    return distanceMatrixPtrs;
}

void getBestTourStartVertexZero(double cost, struct Tour *localBestTour, int *tourTSP, int N) {
    if (cost < localBestTour->cost) {
        localBestTour->cost = cost;
        memcpy(localBestTour->tour, tourTSP, (N + 1) * sizeof(int));
    }
}

void getBestTourRest(double cost, struct Tour *localBestTour, int *tourTSP, int N) {
    long int_cost = (long)(cost * 10000);
    long int_localBestCost = (long)(localBestTour->cost * 10000);

    if (int_cost < int_localBestCost) {
        localBestTour->cost = cost;
        memcpy(localBestTour->tour, tourTSP, (N + 1) * sizeof(int));
    }
}


int main(int argc, char *argv[]) {
    char *fileName = argv[1];
    char *outFileNameC = argv[2];
    char *outFileNameF = argv[3];
    char *outFileNameN = argv[4];

    int N = readNumOfCoords(fileName);
    double **coordsTwodArray = readCoords(fileName, N);
    double **distanceMatrix = calculateDistance(coordsTwodArray, N);

    struct Tour localBestTourC = {DBL_MAX, {-1}};
    struct Tour localBestTourF = {DBL_MAX, {-1}};
    struct Tour localBestTourN = {DBL_MAX, {-1}};
    int *tour = (int *)malloc((N + 1) * sizeof(int));

    for (int i = 0; i < N; i++) {
        int tourC[N + 1], tourF[N + 1], tourN[N + 1];
        double costC = cheapInsertion(distanceMatrix, tourC, i, N);
        double costF = farInsertion(distanceMatrix, tourF, i, N);
        double costN = nearAddition(distanceMatrix, tourN, i, N);

        // Compare and store best tours
        if (i == 0) {
            getBestTourStartVertexZero(costC, &localBestTourC, tourC, N);
            getBestTourStartVertexZero(costF, &localBestTourF, tourF, N);
            getBestTourStartVertexZero(costN, &localBestTourN, tourN, N);
        } else {
            getBestTourRest(costC, &localBestTourC, tourC, N);
            getBestTourRest(costF, &localBestTourF, tourF, N);
            getBestTourRest(costN, &localBestTourN, tourN, N);
        }
    }

    writeTourToFile(localBestTourC.tour, N + 1, outFileNameC);
    writeTourToFile(localBestTourF.tour, N + 1, outFileNameF);
    writeTourToFile(localBestTourN.tour, N + 1, outFileNameN);

    free(distanceMatrix);
    free(tour);
    free(coordsTwodArray);

    return 0;
}