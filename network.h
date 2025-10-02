#ifndef _NETWORK_H_
#define _NETWORK_H_

/*
    this module defines data structure for representing a network,
    ie a list of edges
*/

/*
  An edge of given weight goes from source to dest
  NOTE: weights must be 0/1 for unweighted, or in [0, 1] for weighted
*/
typedef struct {
    unsigned int source;
    unsigned int dest;
    float weight;
} edge;

/*
  nbNodes = total number of nodes in the network
  nbEdges = total number of edges in the network
  edges = an array of nbEdges pointers to edges

  NOTE: if the network is undirected each edge must be present twice (A->B and B->A),
  with the same weight
 */
typedef struct {
    unsigned int nbNodes;
    unsigned int nbEdges;
    edge *edges;
} network;


/*
    check if weights are in ]0, 1],
    set weight=0 for edges that connect a node to itself,
    sort edges by increasing dest then increasing source,
    modifying N in-place;
    return -1 if any weight is not in [0, 1], otherwise the number
    of self-interactions that were set to zero-weight (0 if there were none)
*/
int checkNetwork(network *N);

void freeNetwork(network *N);

void printNetwork(network *N);

#endif
