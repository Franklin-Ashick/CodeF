#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <omp.h>

// Cheapest Insertion Algorithm for TSP
double cheapInsertion(double **distanceMatrix, int *bestTour, int startVertex, int verticesCount) {
    bool *visited = (bool *)calloc(verticesCount, sizeof(bool));

    // Initialization
    int tourSize = 1;
    bestTour[0] = startVertex;
    visited[startVertex] = true;

    // Main loop to build the tour
    while (tourSize < verticesCount) {
        double minCost = DBL_MAX;
        int minPosition = -1;
        int minVertex = -1;

        // Parallel loop for finding the cheapest unvisited vertex
        #pragma omp parallel for
        for (int i = 0; i < verticesCount; i++) {
            if (!visited[i]) {
                double localMinCost = DBL_MAX;
                int localMinPosition = -1;
                int localMinVertex = i;

                for (int j = 0; j < tourSize; j++) {
                    int nextVertex = (j + 1) % tourSize;
                    double insertionCost = distanceMatrix[bestTour[j]][i] + distanceMatrix[i][bestTour[nextVertex]] - distanceMatrix[bestTour[j]][bestTour[nextVertex]];

                    if (insertionCost < localMinCost) {
                        localMinCost = insertionCost;
                        localMinPosition = j + 1; // Insert after j
                    }
                }

                // Update global minimum in a critical section
                #pragma omp critical
                {
                    if (localMinCost < minCost) {
                        minCost = localMinCost;
                        minVertex = localMinVertex;
                        minPosition = localMinPosition;
                    }
                }
            }
        }

        // Insert the vertex into the tour
        for (int i = tourSize; i > minPosition; i--) {
            bestTour[i] = bestTour[i - 1];
        }
        bestTour[minPosition] = minVertex;

        // Update tour size and visited status
        tourSize++;
        visited[minVertex] = true;
    }

    // Closing the tour loop
    bestTour[verticesCount] = bestTour[0];

    // Calculate the total cost of the tour
    double totalTourCost = 0.0;
    for (int i = 0; i < verticesCount; i++) {
        totalTourCost += distanceMatrix[bestTour[i]][bestTour[(i + 1) % verticesCount]];
    }

    free(visited);
    return totalTourCost; // Return the total cost of the tour
}
