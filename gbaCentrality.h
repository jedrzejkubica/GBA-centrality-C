#ifndef _GBACENTRALITY_H_
#define _GBACENTRALITY_H_

#include "pathCounts.h"
#include "adjacency.h"
#include "scores.h"

/*
  gbaCentrality() is the only symbol we want to export in the shared library
  -> we will compile with -fvisibility=hidden and changed to "default" here
*/
#pragma GCC visibility push(default)

/*
  rowSums[i] is the sum of elements in the i-th row of a pathCountsMatrix
*/
typedef struct {
    unsigned int nbNodes;
    float *data;
} rowSums;

/*
    Given a network represented by A and some seed nodes (eg known causal genes),
    apply the GBA-centrality algorithm to calculate a score for each node in
    the network.
    In "causal", each score must be in [0,1] and most should be 0.
    "alpha" is the GBA-centrality parameter, it must be in ]0,1[, typically 0.5 is good.
    "scores" must be allocated and will be filled in-place.
*/
void gbaCentrality(adjacencyMatrix *A, geneScores *causal, float alpha, geneScores *scores);

/*
    this is used to normalize the rows of interactomePathCounts later
*/
rowSums *sumRowElements(pathCountsMatrix *pathCounts);

float calculateScoresDiff(geneScores *scores, geneScores *scoresPrev);

void freeRowSums(rowSums *sums);

#endif
