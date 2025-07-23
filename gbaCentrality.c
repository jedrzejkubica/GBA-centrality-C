#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "gbaCentrality.h"
#include "compactAdjacency.h"
#include "pathCountsWithPredecessors.h"
#include "pathCounts.h"
#include "mem.h"


geneScores *gbaCentrality(adjacencyMatrix *A, geneScores *causal, float alpha) {
    // sanity check:
    if (A->nbCols != causal->nbGenes) {
        fprintf(stderr, "E: gbaCentrality() called with adjacency matrix and causal genes of different sizes");
        exit(1);
    }
	unsigned int nbGenes = causal->nbGenes;

    geneScores *scores = mallocOrDie(sizeof(geneScores), "E: OOM for scores");
    scores->nbGenes = nbGenes;
    scores->scores = mallocOrDie(nbGenes * sizeof(float), "E: OOM for scores->scores");
	// start by copying causal scores, ie scores = alpha**0 * I * casual
	memcpy(scores->scores, causal->scores, nbGenes * sizeof(float));

    unsigned int maxDistance = 5;

    compactAdjacencyMatrix *interactomeComp = adjacency2compact(A);
    pathCountsWithPredMatrix *pathCountsCurrent = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *pathCountsNext = NULL;
	float alphaPowK = 1;
    for (size_t k = 1; k <= maxDistance; k++) {
		// scores += alpha**k * B_k * causal
		alphaPowK *= alpha;
        pathCountsMatrix *interactomePathCounts = countPaths(pathCountsCurrent, interactomeComp);
        for (size_t i = 0; i < nbGenes; i++) {
            for (size_t j = 0; j < nbGenes; j++) {
                scores->scores[i] += alphaPowK * interactomePathCounts->data[i * nbGenes + j] * causal->scores[i];
            }
        }
        freePathCounts(interactomePathCounts);

		if (k < maxDistance) {
			// build B_(k+1) for next iteration
			fprintf(stderr, "I: calculating B_%ld\n", k+1);
			pathCountsNext = buildNextPathCounts(pathCountsCurrent, interactomeComp);
			freePathCountsWithPred(pathCountsCurrent);
			pathCountsCurrent = pathCountsNext;
        }
    }

    freePathCountsWithPred(pathCountsCurrent);
    freeCompactAdjacency(interactomeComp);
    return(scores);
}

void printScores(geneScores *scores) {
    float *currentScoreP = scores->scores;
    for (unsigned int i = 0; i < scores->nbGenes; i++) {
        printf("%f\n", *currentScoreP);
		currentScoreP++;
    }
}

void freeScores(geneScores *scores) {
    free(scores->scores);
    free(scores);
}

