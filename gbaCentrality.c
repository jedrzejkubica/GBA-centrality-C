#include <stdlib.h>
#include <stdio.h>
#include <math.h> 

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
    for (size_t i = 0; i < scores->nbGenes; i++) {
        scores->scores[i] = 0.0;
    }

    unsigned int maxDistance = 2;

    compactAdjacencyMatrix *interactomeComp = adjacency2compact(A);
    pathCountsWithPredMatrix *interactomePathCountsWithPred = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *interactomeNext = NULL;
    for (size_t k = 1; k < maxDistance; k++) {
        printf("calculating A**%ld\n", k+1);
        interactomeNext = buildNextPathCounts(interactomePathCountsWithPred, interactomeComp);
        
        pathCountsMatrix *interactomePathCounts = countPaths(interactomeNext, interactomeComp);

        for (size_t i = 0; i < interactomePathCounts->nbCols; i++) {
            for (size_t j = 0; j < interactomePathCounts->nbCols; j++) {
                scores->scores[i] += pow(alpha, k) * interactomePathCounts->data[i * interactomePathCounts->nbCols + j] * causal->scores[i];
            }
        }
        freePathCounts(interactomePathCounts);
        freePathCountsWithPred(interactomePathCountsWithPred);
        interactomePathCountsWithPred = interactomeNext;
    }

    freePathCountsWithPred(interactomeNext);
    freeCompactAdjacency(interactomeComp);
    freeAdjacency(A);

    return(scores);
}

void printScores(geneScores *scores) {
    float *currentScoreP = scores->scores;
    for (unsigned int i = 0; i < scores->nbGenes; i++) {
        printf("%f\n", *currentScoreP);currentScoreP++;
    }
}

void freeScores(geneScores *scores) {
    free(scores->scores);
    free(scores);
}

