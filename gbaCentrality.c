#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gbaCentrality.h"
#include "compactAdjacency.h"
#include "pathCounts.h"
#include "scores.h"
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

    compactAdjacencyMatrix *interactomeComp = adjacency2compact(A);
    fprintf(stderr, "INFO gbaCentrality.so: calculating B_1\n");

    pathCountsWithPredMatrix *pathCountsCurrent = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *pathCountsNext = NULL;

    geneScores *scoresPrev = mallocOrDie(sizeof(geneScores), "E: OOM for scoresPrev\n");
    scoresPrev->nbGenes = nbGenes;
    scoresPrev->scores = mallocOrDie(sizeof(float) * nbGenes, "E: OOM for scoresPrev scores");

	float alphaPowK = 1;
    size_t k = 1;
    float threshold = 10E-6;
    float scoresDiff = 1;

    while (scoresDiff > threshold) {
        // save scores from B_k-1
        memcpy(scoresPrev->scores, scores->scores, nbGenes * sizeof(float));

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

        // calculate the difference between scores for B_k-1 and B_k
        scoresDiff = calculateScoresDiff(scores, scoresPrev);
        fprintf(stderr, "INFO gbaCentrality.so: scoresDiff = %f\n", scoresDiff);

		if (scoresDiff > threshold) {
			// build B_(k+1) for next iteration
			fprintf(stderr, "INFO gbaCentrality.so: calculating B_%ld\n", k+1);
			pathCountsNext = buildNextPathCounts(pathCountsCurrent, interactomePathCounts, interactomeComp);
			freePathCountsWithPred(pathCountsCurrent);
			pathCountsCurrent = pathCountsNext;
            k++;
        }
        freePathCounts(interactomePathCounts);
    }
    freeScores(scoresPrev);
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