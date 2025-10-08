/*
  Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2025

  This file was written by Jędrzej Kubica and Nicolas Thierry-Mieg
  (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with this program.
  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "compactAdjacency.h"
#include "network.h"
#include "mem.h"


compactAdjacencyMatrix *network2compact(network *N) {
    int nbSelfLoops = checkNetwork(N);
    if (nbSelfLoops == -1) {
        fprintf(stderr, "E: weights are not in ]0, 1], please fix the network\n");
        exit(1);
    } else if (nbSelfLoops > 0) {
        fprintf(stderr, "WARNING: your network has %i self-loops, they have been removed\n", nbSelfLoops);
    }

    fprintf(stderr, "network2compact: %u edges between %u nodes\n", N->nbEdges, N->nbNodes);
    
    compactAdjacencyMatrix *compact = mallocOrDie(sizeof(compactAdjacencyMatrix), "E: OOM for compact\n");

    compact->nbNodes = N->nbNodes;
    compact->offsets = mallocOrDie(sizeof(size_t) * (N->nbNodes + 1), "E: OOM for offsets\n");
    compact->predecessors = mallocOrDie(sizeof(unsigned int) * (N->nbEdges - nbSelfLoops),
                                        "E: OOM for predecessors\n");
    compact->weights = mallocOrDie(sizeof(float) * (N->nbEdges - nbSelfLoops), "E: OOM for weights\n");
    compact->offsetsReverseEdge = mallocOrDie(sizeof(size_t) * (N->nbEdges - nbSelfLoops),
                                              "E: OOM for offsetsReverseEdge\n");

    size_t sumInDegrees = 0;
    unsigned int currentDest = 0;
    edge *currentEdgeP = N->edges;
    compact->offsets[0] = 0;

    for (size_t ei = 0; ei < N->nbEdges; ei++) {
        // ignore zero-weight self-loops
        if (currentEdgeP->weight == 0) {
            currentEdgeP++;
            continue;
        }
        while (currentDest < currentEdgeP->dest)
            compact->offsets[++currentDest] = sumInDegrees;
        compact->predecessors[sumInDegrees] = currentEdgeP->source;
        compact->weights[sumInDegrees] = currentEdgeP->weight;

        sumInDegrees++;
        currentEdgeP++;
    }
    while (currentDest < compact->nbNodes)
        // set offsets[nbNodes] & also if some of the last nodes have no incoming edges
        compact->offsets[++currentDest] = sumInDegrees;

    // sanity check
    assert(sumInDegrees + nbSelfLoops == N->nbEdges);

    // fill offsetsReverseEdge
    for (unsigned int j = 0; j < N->nbNodes; j++) {
        for (size_t offset = compact->offsets[j]; offset < compact->offsets[j + 1]; offset++) {
            unsigned int p = compact->predecessors[offset];
            // find offset of j->p edge if it exists, ie look for j among the predecessors of p
            size_t offsetReverseEdge = sumInDegrees;
            for (size_t offsetRev = compact->offsets[p]; offsetRev < compact->offsets[p + 1]; offsetRev++) {
                if (compact->predecessors[offsetRev] == j) {
                    offsetReverseEdge = offsetRev;
                    break;
                }
            }
            compact->offsetsReverseEdge[offset] = offsetReverseEdge;
        }
    }
    return compact;
}

void freeCompactAdjacency(compactAdjacencyMatrix *compact) {
    free(compact->offsets);
    free(compact->predecessors);
    free(compact->weights);
    free(compact->offsetsReverseEdge);
    free(compact);
}
