#ifndef _GBACENTRALITY_H_
#define _GBACENTRALITY_H_

#include "network.h"
#include "pathCounts.h"
#include "scores.h"

/*
  gbaCentrality() is the only symbol we want to export in the shared library
  -> we will compile with -fvisibility=hidden and changed to "default" here
*/
#pragma GCC visibility push(default)

/*
    Given a network represented by A and some seed nodes (eg known causal genes),
    apply the GBA-centrality algorithm to calculate a score for each node in
    the network.
    In "causal", each score must be in [0,1] and most should be 0.
    "alpha" is the GBA-centrality parameter, it must be in ]0,1[, typically 0.5 is good.
    "scores" must be allocated and will be filled in-place.
*/
void gbaCentrality(network *N, geneScores *causal, float alpha, geneScores *scores);

#endif
