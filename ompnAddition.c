#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define EPSILON 1e-9 // For floating-point comparisons

float nearAddition(double **distanceMatrix, int *bestTour, int startVertex, int verticesCount) {
    bool *visited = (bool *)calloc(verticesCount, sizeof(bool));
    if (!visited) {
        perror("Memory Allocation Failed for Visited Array");
        exit(EXIT_FAILURE);
    }

    // Initialize the tour with the start vertex
    bestTour[0] = startVertex;
    visited[startVertex] = true;

    // Build the tour by adding the nearest unvisited vertex
    for (int tourSize = 1; tourSize < verticesCount; tourSize++) {
        double minDistance = DBL_MAX;
        int nearestVertex = -1;

        // Parallel region to find the nearest unvisited vertex
        #pragma omp parallel for reduction(min:minDistance)
        for (int i = 0; i < verticesCount; i++) {
            if (!visited[i]) {
                double distance = distanceMatrix[bestTour[tourSize - 1]][i];
                if (distance < minDistance) {
                    #pragma omp critical
                    {
                        if (distance < minDistance) {
                            minDistance = distance;
                            nearestVertex = i;
                        }
                    }
                }
            }
        }

        // Add the nearest unvisited vertex to the tour
        visited[nearestVertex] = true;
        bestTour[tourSize] = nearestVertex;
    }

    // Close the tour by returning to the start vertex
    bestTour[verticesCount] = bestTour[0];

    // Calculate the total cost of the tour
    double totalTourCost = 0.0;
    for (int i = 0; i < verticesCount; i++) {
        totalTourCost += distanceMatrix[bestTour[i]][bestTour[(i + 1) % verticesCount]];
    }

    free(visited);
    return totalTourCost;
}
