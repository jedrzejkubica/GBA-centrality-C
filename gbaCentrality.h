#ifndef _GBACENTRALITY_H_
#define _GBACENTRALITY_H_

#include "adjacency.h"


/*
    Store one score per gene, each score is a float >=0.
    scores MUST be large enough to store nbGenes floats.
 */
typedef struct {
    unsigned int nbGenes;
    float *scores;
} geneScores;


/*
    Given a network represented by A and some seed nodes (eg known causal genes),
    apply the GBA-centrality algorithm to calculate a score for each node in
    the network.
    In "causal", each score must be in [0,1] and most should be 0.
    "alpha" is the GBA-centrality parameter, it must be in ]0,1[, typically 0.5 is good.
    Returns a pointer to a newly allocated structure.
*/
geneScores *gbaCentrality(adjacencyMatrix *A, geneScores *causal, float alpha);

void printScores(geneScores *scores);

void freeScores(geneScores *scores);

#endif
