#ifndef _COMPACTADJACENCY_H_
#define _COMPACTADJACENCY_H_

#include <stddef.h>
#include "network.h"

/*
    this module defines the compact representation of an adjacency matrix
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

    offsetsReverseEdge has sumOfDegrees elements (same as predecessors and weights),
    for any i < sumOfDegrees offsetsReverseEdge[i] is the offset corresponding to 
    edge j->p (assuming i is the offset of edge p->j) if j->p exists, sumOfDegrees
    otherwise
*/
typedef struct {
    unsigned int nbNodes;
    size_t *offsets;
    unsigned int *predecessors;
    float *weights;
    size_t *offsetsReverseEdge;
} compactAdjacencyMatrix;

/*
  Allocate and populate a compactAdjacencyMatrix from a network.
  Pre-conditions: the network's edges must be sorted by increasing
  dest then increasing source, and must not contain self-interactions;
  these are true if the network went through checkNetwork()
 */
compactAdjacencyMatrix *network2compact(network *N);

void freeCompactAdjacency(compactAdjacencyMatrix *compactA);

#endif
