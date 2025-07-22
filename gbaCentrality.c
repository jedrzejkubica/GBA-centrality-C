#include <stdlib.h>
#include <stdio.h>

#include "gbaCentrality.h"
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

    // TODO: code GBA-centrality algorithm
    
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

