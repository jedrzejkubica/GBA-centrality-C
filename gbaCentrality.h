#ifndef _GBACENTRALITY_H_
#define _GBACENTRALITY_H_

#include "adjacency.h"
#include "scores.h"

/*

/*
    Given a network represented by A and some seed nodes (eg known causal genes),
    apply the GBA-centrality algorithm to calculate a score for each node in
    the network.
    In "causal", each score must be in [0,1] and most should be 0.
    "alpha" is the GBA-centrality parameter, it must be in ]0,1[, typically 0.5 is good.
    Returns a pointer to a newly allocated structure.
*/
geneScores *gbaCentrality(adjacencyMatrix *A, geneScores *causal, float alpha);

#endif
