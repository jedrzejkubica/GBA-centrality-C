#include <stdio.h>

#include "adjacency.h"
#include "scores.h"
#include "gbaCentrality.h"
#include "mem.h"


adjacencyMatrix *diamond4(void) {
    // diamond matrix for testing
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

adjacencyMatrix *asymmetric1and4(void) {
    adjacencyMatrix *asymmetric1and4 = mallocOrDie(sizeof(adjacencyMatrix), "E: OOM for asymmetric1and4");
    asymmetric1and4->nbCols = 7;
    asymmetric1and4->weights = mallocOrDie(7 * 7 * sizeof(float), "E: OOM for asymmetric1and4 weights");
    for (size_t i = 0; i < 49; i++)
        asymmetric1and4->weights[i] = 0;

    asymmetric1and4->weights[1] = 1.0;
    asymmetric1and4->weights[7] = 1.0;
    asymmetric1and4->weights[9] = 1.0;
    asymmetric1and4->weights[15] = 1.0;
    asymmetric1and4->weights[17] = 1.0;
    asymmetric1and4->weights[18] = 1.0;
    asymmetric1and4->weights[19] = 1.0;
    asymmetric1and4->weights[20] = 1.0;
    asymmetric1and4->weights[23] = 1.0;
    asymmetric1and4->weights[30] = 1.0;
    asymmetric1and4->weights[37] = 1.0;
    asymmetric1and4->weights[44] = 1.0;
    return(asymmetric1and4);
}

geneScores *causalGenes(void) {
    geneScores *causal = mallocOrDie(sizeof(geneScores), "E: OOM for causal genes");
    causal->nbGenes = 7;
    causal->scores = mallocOrDie(7 * sizeof(SCORETYPE), "E: OOM for known causal genes");

    for (size_t i = 0; i < 7; i++)
        causal->scores[i] = 0;

    //causal->scores[1] = 1;
    causal->scores[2] = 1;

    return(causal);
}

int main(void) {
    /*
    adjacencyMatrix *diam4 = diamond4();
    printf("diam4 with its self-loop\n");
    printAdjacency(diam4);
    */
    adjacencyMatrix *asymm1to4 = asymmetric1and4();
    printAdjacency(asymm1to4);

    geneScores *causal = causalGenes();
    printf("causal genes\n");
    printScores(causal);

    geneScores *result = mallocOrDie(sizeof(geneScores), "E: OOM for result scores");
    result->nbGenes = 7;
    result->scores = mallocOrDie(7 * sizeof(SCORETYPE), "E: OOM for result scores");

    gbaCentrality(asymm1to4, causal, 0.5, result);
    printf("scores\n");
    printScores(result);
    freeAdjacency(asymm1to4);
    freeScores(result);
    freeScores(causal);

    return(0);
}

