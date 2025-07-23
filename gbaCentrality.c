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

    geneScores *scores = mallocOrDie(sizeof(geneScores), "E: OOM for scores");
    scores->nbGenes = causal->nbGenes;
    scores->scores = mallocOrDie(causal->nbGenes * sizeof(float), "E: OOM for scores->scores");
	// start by copying causal scores, ie scores = alpha**0 * I * casual
	memcpy(scores->scores, causal->scores, nbGenes * sizeof(float));

    unsigned int maxDistance = 5;

    compactAdjacencyMatrix *interactomeComp = adjacency2compact(A);
    pathCountsWithPredMatrix *pathCountsCurrent = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *pathCountsNext = NULL;
	float alphaPowK = alpha;
    for (size_t k = 1; k < maxDistance; k++) {
        fprintf(stderr, "I: calculating B_%ld\n", k+1);
        pathCountsNext = buildNextPathCounts(pathCountsCurrent, interactomeComp);
        
        pathCountsMatrix *interactomePathCounts = countPaths(pathCountsNext, interactomeComp);

        for (size_t i = 0; i < interactomePathCounts->nbCols; i++) {
            for (size_t j = 0; j < interactomePathCounts->nbCols; j++) {
                scores->scores[i] += alphaPowK * interactomePathCounts->data[i * interactomePathCounts->nbCols + j] * causal->scores[i];
            }
        }
		alphaPowK *= alpha;
        freePathCounts(interactomePathCounts);
        freePathCountsWithPred(pathCountsCurrent);
        pathCountsCurrent = pathCountsNext;
    }

    freePathCountsWithPred(pathCountsNext);
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

