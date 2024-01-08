#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <omp.h>

// Farthest Insertion Algorithm for TSP
double farInsertion(double **distanceMatrix, int *bestTour, int startVertex, int verticesCount) {
    bool *visited = (bool *)calloc(verticesCount, sizeof(bool));

    // Initialization
    int tourSize = 1;
    bestTour[0] = startVertex;
    visited[startVertex] = true;

    // Main loop to build the tour
    while (tourSize < verticesCount) {
        double maxDistance = -1.0;
        int farthestVertex = -1;

        // Parallel loop for finding the farthest unvisited vertex
        #pragma omp parallel
        {
            double localMaxDistance = -1.0;
            int localFarthestVertex = -1;

            #pragma omp for nowait
            for (int i = 0; i < verticesCount; i++) {
                if (!visited[i]) {
                    for (int j = 0; j < tourSize; j++) {
                        if (distanceMatrix[bestTour[j]][i] > localMaxDistance) {
                            localMaxDistance = distanceMatrix[bestTour[j]][i];
                            localFarthestVertex = i;
                        }
                    }
                }
            }

            // Critical section to update the global maximum
            #pragma omp critical
            {
                if (localMaxDistance > maxDistance) {
                    maxDistance = localMaxDistance;
                    farthestVertex = localFarthestVertex;
                }
            }
        }

        // Determine the best position to insert the farthest vertex
        double minInsertionCost = DBL_MAX;
        int bestPosition = -1;

        for (int i = 0; i < tourSize; i++) {
            int next = (i + 1) % tourSize;
            double insertionCost = distanceMatrix[bestTour[i]][farthestVertex] + distanceMatrix[farthestVertex][bestTour[next]] - distanceMatrix[bestTour[i]][bestTour[next]];

            if (insertionCost < minInsertionCost) {
                minInsertionCost = insertionCost;
                bestPosition = i + 1;
            }
        }

        // Insert the farthest vertex at the best position
        for (int i = tourSize; i > bestPosition; i--) {
            bestTour[i] = bestTour[i - 1];
        }
        bestTour[bestPosition] = farthestVertex;

        // Update tour size and visited status
        tourSize++;
        visited[farthestVertex] = true;
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
