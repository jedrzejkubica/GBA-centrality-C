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

#include "pathCounts.h"
#include "compactAdjacency.h"
#include "mem.h"


pathCountsWithPredMatrix *buildFirstPathCounts(compactAdjacencyMatrix *compact) {
    pathCountsWithPredMatrix *pathCounts = mallocOrDie(sizeof(pathCountsWithPredMatrix), "E: OOM for path counts\n");
    
    unsigned int nbNodes = compact->nbNodes;
    unsigned int sumDegrees = compact->offsets[nbNodes];
    
    pathCounts->data = mallocOrDie(sizeof(PATHCOUNTSTYPE) * sumDegrees * nbNodes, "E: OOM for path counts data\n");
    // set to 0.0 (all-zeroes is not guaranteed to be 0.0)
    for (size_t i = 0; i < sumDegrees * nbNodes; i++)
        pathCounts->data[i] = 0.0;

    // if i->j is an edge of weight w, then there is a path from i to j
    // with penultimate node i and weight w
    for (size_t offset = 0; offset < sumDegrees; offset++)
        pathCounts->data[compact->predecessors[offset] * sumDegrees + offset] = compact->weights[offset];

    return pathCounts;
}

pathCountsWithPredMatrix *buildNextPathCounts(pathCountsWithPredMatrix *pathCountsWithPred, pathCountsMatrix *pathCounts,
                                              compactAdjacencyMatrix *compact) {
    unsigned int nbNodes = compact->nbNodes;
    unsigned int sumDegrees = compact->offsets[nbNodes];
    
    pathCountsWithPredMatrix *nextPathCounts = mallocOrDie(sizeof(pathCountsWithPredMatrix), "E: OOM for next path counts\n");
    nextPathCounts->data = mallocOrDie(sizeof(PATHCOUNTSTYPE) * sumDegrees * nbNodes, "E: OOM for next path counts data\n");
    
    /*
      Ideally we want to count paths, ie walks that don't contain loops. But this is hard
      to code efficiently...
      We therefore have an approximate solution, that counts all walks except those that:
      - contain a loop (of any length) back to the starting node (with i!=j below);
      - contain a loop of length 2, ie with 2 steps back-and-forth, anywhere in
        the walk (using offsetReverseEdge below).
      Therefore we are still counting walks that contain an inner loop of length >= 3  (eg
      i->j->k->m->j will be counted as a length 4 path between i and j)...
      But our approximate solution is efficient, and should largely reduce the
      "hub" issue, ie the explosion of number of walks when going through a hub.
    */
    for (size_t i = 0; i < nbNodes; i++) {
        for (size_t j = 0; j < nbNodes; j++) {
            for (size_t offset = compact->offsets[j]; offset < compact->offsets[j + 1]; offset++) {
                double sum = 0;
                // ignore loops back to starting node => path count stays zero if i==j
                if (i != j) {
                    unsigned int p = compact->predecessors[offset];
                    sum = pathCounts->data[i * nbNodes + p];
                    // ignore walks whose last step is backtracking on the previous step
                    if (compact->offsetsReverseEdge[offset] < sumDegrees) {
                        sum -= pathCountsWithPred->data[i * sumDegrees + compact->offsetsReverseEdge[offset]];
                    }
                    sum *= compact->weights[offset];
                }
                nextPathCounts->data[i * sumDegrees + offset] = (PATHCOUNTSTYPE)sum;
            }
        }
    }
    
    return nextPathCounts;
}

void freePathCountsWithPred(pathCountsWithPredMatrix *pathCounts) {
    free(pathCounts->data);
    free(pathCounts);
}

pathCountsMatrix *countPaths(pathCountsWithPredMatrix *pathCountsWithPred, compactAdjacencyMatrix *compact) {
    pathCountsMatrix *pathCounts = mallocOrDie(sizeof(pathCountsMatrix), "E: OOM for path counts matrix\n");
    
    unsigned int nbNodes = compact->nbNodes;
    unsigned int sumDegrees = compact->offsets[nbNodes];
    
    pathCounts->nbCols = nbNodes;
    pathCounts->data = mallocOrDie(sizeof(PATHCOUNTSTYPE) * nbNodes * nbNodes, "E: OOM for path counts data\n");
    
    for (size_t i = 0; i < nbNodes; i++) {
        for (size_t j = 0; j < nbNodes; j++) {
            double sum = 0;
            for (size_t k = compact->offsets[j]; k < compact->offsets[j + 1]; k++)
                sum += pathCountsWithPred->data[i * sumDegrees + k];

            pathCounts->data[i * nbNodes + j] = (PATHCOUNTSTYPE)sum;
        }
    }
    
    return pathCounts;
}

void printPathCounts(pathCountsMatrix *pathCounts) {
    for (size_t i = 0; i < pathCounts->nbCols; i++) {
        for (size_t j = 0; j < pathCounts->nbCols; j++)
            printf("%0.2f ", pathCounts->data[i * pathCounts->nbCols + j]);

        printf("\n");
    }
}

void freePathCounts(pathCountsMatrix *pathCounts) {
    free(pathCounts->data);
    free(pathCounts);
}
