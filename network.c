#include <stdio.h>
#include <stdlib.h>

#include "network.h"
#include "mem.h"


/*
  compare two edges so they can be sorted by increasing dest then
  increasing source
 */
int compEdges(const void *elem1, const void *elem2) {
    int retVal = 0;
    edge *edge1 = (edge*) elem1;
    edge *edge2 = (edge*) elem2;

    if (edge1->dest < edge2->dest)
        retVal = -1;
    else if (edge1->dest > edge2->dest)
        retVal = 1;
    else if (edge1->source < edge2->source)
        retVal = -1;
    else if (edge1->source > edge2->source)
        retVal = 1;
    return(retVal);
}

int checkNetwork(network *N) {
    int nbSelfLoops = 0;

    edge *edgeP = N->edges;
    for (size_t i = 0; i < N->nbEdges; i++) {
        // check if weights are in ]0, 1]
        if ((edgeP->weight <= 0) || (edgeP->weight > 1)) {
            // immediately return, we can't do anything with this network
            return(-1);
        }
        // remove self-loops
        if (edgeP->source == edgeP->dest) {
            edgeP->weight = 0;
            nbSelfLoops++;
        }
        edgeP++;
    }
    // sort the edges in-place
    qsort(N->edges, N->nbEdges, sizeof(edge), compEdges);
    return(nbSelfLoops);
}

void freeNetwork(network *N) {
    free(N->edges);
    free(N);
}

void printNetwork(network *N) {
    printf("Network has %u nodes and %u edges:\n", N->nbNodes, N->nbEdges);
    edge *edgeP = N->edges;
    for (size_t i = 0; i < N->nbEdges; i++) {
        printf("%u -> %u w=%0.2f\n", edgeP->source, edgeP->dest, edgeP->weight);
        edgeP++;
    }
}
