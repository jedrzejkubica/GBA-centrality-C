#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gbaCentrality.h"
#include "compactAdjacency.h"
#include "pathCounts.h"
#include "scores.h"
#include "mem.h"


/*
  Return sqrt[sum of square of (diff between same-index elements)] 
*/
float calculateScoresDiff(geneScores *scores, geneScores *scoresPrev);


void gbaCentrality(adjacencyMatrix *A, geneScores *causal, float alpha, geneScores *scores) {
    // sanity check:
    if (A->nbCols != causal->nbGenes) {
        fprintf(stderr, "E: gbaCentrality() called with adjacency matrix and causal genes of different sizes");
        exit(1);
    }
    unsigned int nbGenes = causal->nbGenes;

    // start by copying causal scores, ie scores = alpha**0 * I * casual
    memcpy(scores->scores, causal->scores, nbGenes * sizeof(SCORETYPE));

    compactAdjacencyMatrix *interactomeComp = adjacency2compact(A);
    fprintf(stderr, "INFO gbaCentrality(): calculating B_1\n");

    pathCountsWithPredMatrix *pathCountsCurrent = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *pathCountsNext = NULL;

    geneScores *scoresPrev = mallocOrDie(sizeof(geneScores), "E: OOM for scoresPrev\n");
    scoresPrev->nbGenes = nbGenes;
    scoresPrev->scores = mallocOrDie(sizeof(SCORETYPE) * nbGenes, "E: OOM for scoresPrev scores");

    float alphaPowK = 1;
    size_t k = 1;
    
    // for convergence test
    float threshold = 10E-4;
    float scoresDiff = 1;

    while (scoresDiff > threshold) {
        // save scores from B_k-1
        memcpy(scoresPrev->scores, scores->scores, nbGenes * sizeof(SCORETYPE));

        // scores += alpha**k * B_k * causal
        alphaPowK *= alpha;
        pathCountsMatrix *interactomePathCounts = countPaths(pathCountsCurrent, interactomeComp);

        for (size_t i = 0; i < nbGenes; i++) {
            // calculate sum of row i
            double rowSum = 0;
            for (size_t j = 0; j < nbGenes; j++) {
                rowSum += interactomePathCounts->data[i * nbGenes + j];
            }
            // updates scores
            if (rowSum != 0) {
                double scoreSum = 0;
                for (size_t j = 0; j < nbGenes; j++) {
                    scoreSum += interactomePathCounts->data[i * nbGenes + j] * causal->scores[j];
                }
                scores->scores[i] += alphaPowK * scoreSum / rowSum;
            }
        }
        
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

float calculateScoresDiff(geneScores *scores, geneScores *scoresPrev) {
    double scoresDiff = 0;
    for (size_t i = 0; i < scores->nbGenes; i++) {
        float diff = scores->scores[i] - scoresPrev->scores[i];
        scoresDiff += diff * diff;
    }
    scoresDiff = sqrt(scoresDiff);  // L2 norm

    return((float)scoresDiff);
}
