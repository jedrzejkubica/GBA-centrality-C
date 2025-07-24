#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "gbaCentrality.h"
#include "compactAdjacency.h"
#include "pathCounts.h"
#include "mem.h"


void gbaCentrality(adjacencyMatrix *A, geneScores *causal, float alpha, geneScores *scores) {
    // sanity check:
    if (A->nbCols != causal->nbGenes) {
        fprintf(stderr, "E: gbaCentrality() called with adjacency matrix and causal genes of different sizes");
        exit(1);
    }
	unsigned int nbGenes = causal->nbGenes;

	// start by copying causal scores, ie scores = alpha**0 * I * casual
	memcpy(scores->scores, causal->scores, nbGenes * sizeof(float));

    unsigned int maxDistance = 5;

    compactAdjacencyMatrix *interactomeComp = adjacency2compact(A);
    fprintf(stderr, "I: calculating B_1\n");
    pathCountsWithPredMatrix *pathCountsCurrent = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *pathCountsNext = NULL;
	float alphaPowK = 1;
    for (size_t k = 1; k <= maxDistance; k++) {
		// scores += alpha**k * B_k * causal
		alphaPowK *= alpha;
        pathCountsMatrix *interactomePathCounts = countPaths(pathCountsCurrent, interactomeComp);
        for (size_t i = 0; i < nbGenes; i++) {
            for (size_t j = 0; j < nbGenes; j++) {
                scores->scores[j] += alphaPowK * interactomePathCounts->data[i * nbGenes + j] * causal->scores[i];
            }
        }

		if (k < maxDistance) {
			// build B_(k+1) for next iteration
			fprintf(stderr, "I: calculating B_%ld\n", k+1);
			pathCountsNext = buildNextPathCounts(pathCountsCurrent, interactomePathCounts, interactomeComp);
			freePathCountsWithPred(pathCountsCurrent);
			pathCountsCurrent = pathCountsNext;
        }
        freePathCounts(interactomePathCounts);
    }

    freePathCountsWithPred(pathCountsCurrent);
    freeCompactAdjacency(interactomeComp);
}
