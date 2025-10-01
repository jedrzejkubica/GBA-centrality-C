#ifndef _PATHCOUNTS_H_
#define _PATHCOUNTS_H_

#include "compactAdjacency.h"


/*
  type to store the sum of path weights of a given length, floating point number
  of some precision (eg maybe double)
*/
#define PATHCOUNTSTYPE float


/*
    this structure is dependent on a compactAdjacencyMatrix;
    data is of size offsets[nbNodes]*nbNodes;
    for any k in [0, deg(j)-1]:
    data[i*offsets[nbNodes] + offsets[j] + k]
    is the number of paths (or the sum of path weights for a weighted network)
    between nodes i and j whose penultimate node is predecessors[offsets[j] + k]
*/
typedef struct {
    PATHCOUNTSTYPE *data;
} pathCountsWithPredMatrix;


/*
    data[i*nbCols + j] contains the number of paths
    (or the sum of path weights for a weighted network)
    of a given length from node i to j
*/
typedef struct {
    unsigned int nbCols;
    PATHCOUNTSTYPE *data;
} pathCountsMatrix;


/*
    build pathCountsWithPredMatrix for paths of length 1,
    return a pointer to a freshly allocated structure;
*/
pathCountsWithPredMatrix *buildFirstPathCounts(compactAdjacencyMatrix *compact);

/*
    build a pathCountsWithPredMatrix for paths of length k+1
    given pathCountsWithPredMatrix for paths of length k
    (excluding paths looping back to the starting node),
    return a pointer to a freshly allocated structure;
*/
pathCountsWithPredMatrix *buildNextPathCounts(pathCountsWithPredMatrix *pathCountsWithPred, pathCountsMatrix *pathCounts,
                                              compactAdjacencyMatrix *compact);

void freePathCountsWithPred(pathCountsWithPredMatrix *pathCounts);


/*  
    produce pathCountsMatrix corresponding to pathCountsWithPredMatrix
*/
pathCountsMatrix *countPaths(pathCountsWithPredMatrix *pathCountsWithPred, compactAdjacencyMatrix *compact);

void printPathCounts(pathCountsMatrix *pathCounts);

void freePathCounts(pathCountsMatrix *pathCounts);

#endif
