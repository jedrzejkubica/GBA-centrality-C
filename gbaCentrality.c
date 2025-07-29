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

    unsigned int maxDistance = 2;

    compactAdjacencyMatrix *interactomeComp = adjacency2compact(A);
    fprintf(stderr, "INFO gbaCentrality.so: calculating B_1\n");
    pathCountsWithPredMatrix *pathCountsCurrent = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *pathCountsNext = NULL;
	float alphaPowK = 1;
    for (size_t k = 1; k <= maxDistance; k++) {
		// scores += alpha**k * B_k * causal
		alphaPowK *= alpha;
        pathCountsMatrix *interactomePathCounts = countPaths(pathCountsCurrent, interactomeComp);
        rowSums *sums = sumRowElements(interactomePathCounts);  // for row-wise normalization of interactomePathCounts
        for (size_t i = 0; i < nbGenes; i++) {
            for (size_t j = 0; j < nbGenes; j++) {
                if (sums->data[i] == 0) {
                    scores->scores[j] += alphaPowK * interactomePathCounts->data[i * nbGenes + j] * causal->scores[i];
                } else {
                    scores->scores[j] += alphaPowK * interactomePathCounts->data[i * nbGenes + j] / sums->data[i] * causal->scores[i];
                }
            }
        }
        freeRowSums(sums);

		if (k < maxDistance) {
			// build B_(k+1) for next iteration
			fprintf(stderr, "INFO gbaCentrality.so: calculating B_%ld\n", k+1);
			pathCountsNext = buildNextPathCounts(pathCountsCurrent, interactomePathCounts, interactomeComp);
			freePathCountsWithPred(pathCountsCurrent);
			pathCountsCurrent = pathCountsNext;
        }
        freePathCounts(interactomePathCounts);
    }

    freePathCountsWithPred(pathCountsCurrent);
    freeCompactAdjacency(interactomeComp);
}

rowSums *sumRowElements(pathCountsMatrix *pathCounts) {
    rowSums *sums = mallocOrDie(sizeof(rowSums), "E: OOM for row sums\n");\
    unsigned int nbNodes = pathCounts->nbCols;
    sums->nbNodes = nbNodes;
    sums->data = mallocOrDie(sizeof(unsigned int) * nbNodes, "E: OOM for row sums data\n");
    for (size_t i = 0; i < sums->nbNodes; i++) {
        sums->data[i] = 0;
        for (size_t j = 0; j < pathCounts->nbCols; j++) {
            sums->data[i] += pathCounts->data[i * pathCounts->nbCols + j];
        }
    }
    return sums;
}

float calculateScoresDiff(geneScores *scores, geneScores *scoresPrev) {
    float scoresDiff = 0;
    for (size_t i = 0; i < scores->nbGenes; i++) {
        float diff = scores->scores[i] - scoresPrev->scores[i];
        scoresDiff += diff * diff;
    }
    scoresDiff = sqrt(scoresDiff);  // L2 norm

    return scoresDiff;
}

void freeRowSums(rowSums *sums) {
    if (sums) {
        free(sums->data);
        free(sums);
    }
}