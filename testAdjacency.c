#include <stdio.h>

#include "adjacency.h"
#include "compactAdjacency.h"
#include "pathCountsWithPredecessors.h"
#include "pathCounts.h"
#include "gbaCentrality.h"
#include "mem.h"


adjacencyMatrix *diamond4(void) {
    adjacencyMatrix *diamond4 = mallocOrDie(sizeof(adjacencyMatrix), "E: OOM for diamond4");
    diamond4->nbCols = 4;
    diamond4->weights = mallocOrDie(16 * sizeof(float), "E: OOM for diamond4 weights");
    for (size_t i = 0; i < 16; i++)
        diamond4->weights[i] = 0;

    diamond4->weights[1] = 0.2;
    diamond4->weights[2] = 0.5;
    diamond4->weights[4] = 0.4;
    diamond4->weights[7] = 1;
    diamond4->weights[8] = 1;
    diamond4->weights[11] = 1;
    diamond4->weights[13] = 1;
    diamond4->weights[14] = 1;
    // add a self-interaction 4->4 for testing
    diamond4->weights[15] = 1;
    return(diamond4);
}

geneScores *causalGenes(void) {
    geneScores *causal = mallocOrDie(sizeof(geneScores), "E: OOM for causal genes");
    causal->nbGenes = 4;
    causal->scores = mallocOrDie(4 * sizeof(float), "E: OOM for known causal genes");

    for (size_t i = 0; i < 4; i++)
        causal->scores[i] = 0;

    causal->scores[0] = 1;
    causal->scores[3] = 1;

    return(causal);
}

int main(void) {
    adjacencyMatrix *diam4 = diamond4();
    printf("diam4 with its self-loop\n");
    printAdjacency(diam4);

    geneScores *causal = causalGenes();
    printf("causal genes\n");
    printScores(causal);

    geneScores *scores = gbaCentrality(diam4, causal, 0.5);
    printf("scores\n");
    printScores(scores);

    return(0);
}

