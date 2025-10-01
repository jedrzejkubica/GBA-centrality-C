#ifndef _SCORES_H_
#define _SCORES_H_


// type used for each geneScores->scores value
#define SCORETYPE double

/*
    Store one score per gene, each score is a float >=0.
    scores MUST be large enough to store nbGenes floats.
 */
typedef struct {
    unsigned int nbGenes;
    SCORETYPE *scores;
} geneScores;


void printScores(geneScores *scores);

void freeScores(geneScores *scores);

#endif
