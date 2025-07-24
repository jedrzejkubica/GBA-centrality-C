#include <stdio.h>
#include <stddef.h>

#include "compactAdjacency.h"
#include "adjacency.h"
#include "mem.h"


compactAdjacencyMatrix *adjacency2compact(adjacencyMatrix *A) {
    int checkA = checkAdjacency(A);
    if (checkA == 1) {
        fprintf(stderr, "E: weights are not in [0, 1], please fix the network\n");
        exit(1);
    } else if (checkA == 2) {
        fprintf(stderr, "W: your network has self-loops, they have been removed\n");
    }

    compactAdjacencyMatrix *compact = mallocOrDie(sizeof(compactAdjacencyMatrix), "E: OOM for compact\n");

    compact->nbNodes = A->nbCols;
    compact->offsets = mallocOrDie(sizeof(size_t) * (A->nbCols + 1), "E: OOM for offsets\n");
    
    compact->offsets[0] = 0;
    size_t sumInDegrees = 0;
    for (size_t j = 0; j < A->nbCols; j++) {
        // in-degree of node j
        size_t inDegree = 0;
        for (size_t i = 0; i < A->nbCols; i++) {
            if (A->weights[i * A->nbCols + j] > 0) {
                inDegree++;
            }
        }
        sumInDegrees += inDegree;
        compact->offsets[j+1] = sumInDegrees;
    }

    compact->predecessors = mallocOrDie(sizeof(unsigned int) * sumInDegrees, "E: OOM for predecessors\n");
    compact->weights = mallocOrDie(sizeof(float) * sumInDegrees, "E: OOM for weights in compact\n");

    size_t idxPred = 0;
    for (size_t j = 0; j < A->nbCols; j++) {
        for (size_t i = 0; i < A->nbCols; i++) {
            if (A->weights[i * A->nbCols + j] > 0) {
                compact->predecessors[idxPred] = i;
                compact->weights[idxPred] = A->weights[i * A->nbCols + j];
                idxPred++;
            }
        }
    }

    compact->offsetsReverseEdge = mallocOrDie(sizeof(size_t) * sumInDegrees,
											  "E: OOM for offsetsReverseEdge in compact\n");
    for (unsigned int j = 0; j < A->nbCols; j++) {
		for (size_t offset = compact->offsets[j]; offset < compact->offsets[j + 1]; offset++) {
			unsigned int p = compact->predecessors[offset];
			// find offset of j->p edge if it exists
			size_t offsetReverseEdge = sumInDegrees;
			if (A->weights[j * A->nbCols + p] > 0) {
				// look for j among the predecessors of p
				for (size_t offsetRev = compact->offsets[p]; offsetRev < compact->offsets[p + 1]; offsetRev++) {
					if (compact->predecessors[offsetRev] == j) {
						offsetReverseEdge = offsetRev;
						break;
					}
				}
			}
			compact->offsetsReverseEdge[offset] = offsetReverseEdge;
		}
	}
    return compact;
}

void freeCompactAdjacency(compactAdjacencyMatrix *compact) {
    free(compact->offsets);
    free(compact->predecessors);
    free(compact->weights);
    free(compact->offsetsReverseEdge);
    free(compact);
}
