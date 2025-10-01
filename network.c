#include <stdio.h>
#include <string.h>

#include "network.h"
#include "mem.h"

int checkNetwork(network *N) {
    int retVal = 0;

    edge *edgeP = N->edges;
    for (size_t i = 0; i < N->nbEdges; i++) {
        // remove self-loops
        if (edgeP->source == edgeP->dest) {
            edgeP->weight = 0;
            retVal = 2;
        }
        // check if weights are in ]0, 1]
        if ((edgeP->weight <= 0) || (edgeP->weight > 1)) {
            // immediately return, we can't do anything with this network
            return(1);
        }
        edgeP++;
    }
    return(retVal);
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
    }
}
