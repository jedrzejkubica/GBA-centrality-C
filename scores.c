#include <stdlib.h>
#include <stdio.h>

#include "scores.h"


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

