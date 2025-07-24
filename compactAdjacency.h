#ifndef _COMPACTADJACENCY_H_
#define _COMPACTADJACENCY_H_

#include <stddef.h>
#include "adjacency.h"

/*
    this module defines the compact representation of adjacency matrix
*/

/*
    compactAdjacencyMatrix represents an adjacency matrix of size nbNodes x nbNodes;

    offsets has nbNodes+1 elements,
    offsets[j] is the first index corresponding to node j in predecessors and weights,
    offsets[nbNodes] is the sum of degrees, used to avoid overflowing;

    nodes with an edge going into node j are nodes predecessors[offsets[j]]
    up to predecessors[offsets[j+1]-1];

    the corresponding edge weights (coming into node j) are weights[offsets[j]]
    up to weights[offsets[j+1]-1]
*/
typedef struct {
    unsigned int nbNodes;
    size_t *offsets;
    unsigned int *predecessors;
    float *weights;
} compactAdjacencyMatrix;

compactAdjacencyMatrix *adjacency2compact(adjacencyMatrix *A);

void freeCompactAdjacency(compactAdjacencyMatrix *compactA);

#endif
