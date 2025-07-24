#include "pathCountsWithPredecessors.h"
#include "pathCounts.h"
#include "compactAdjacency.h"
#include "adjacency.h"
#include "mem.h"


pathCountsWithPredMatrix *buildFirstPathCounts(compactAdjacencyMatrix *compact) {
    pathCountsWithPredMatrix *pathCounts = mallocOrDie(sizeof(pathCountsWithPredMatrix), "E: OOM for path counts\n");
    
    unsigned int nbNodes = compact->nbNodes;
    unsigned int sumDegrees = compact->offsets[nbNodes];
    
    pathCounts->data = mallocOrDie(sizeof(float) * sumDegrees * nbNodes, "E: OOM for path counts data\n");
    // set to 0.0 float (all-zeroes is not guaranteed to be 0.0)
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
    nextPathCounts->data = mallocOrDie(sizeof(float) * sumDegrees * nbNodes, "E: OOM for next path counts data\n");
    
    for (size_t i = 0; i < nbNodes; i++) {
        for (size_t j = 0; j < nbNodes; j++) {
            for (size_t offset = compact->offsets[j]; offset < compact->offsets[j + 1]; offset++) {
                float sum = 0;
                // we want to ignore loops so path count is zero if i==j
                if (i != j) {
                    unsigned int p = compact->predecessors[offset];
					sum = pathCounts->data[i * nbNodes + p];
					if (compact->offsetsReverseEdge[offset] < sumDegrees) {
						sum -= compact->offsetsReverseEdge[offset];
					}
                    sum *= compact->weights[offset];
                }
                nextPathCounts->data[i * sumDegrees + offset] = sum;
            }
        }
    }
    
    return nextPathCounts;
}

void freePathCountsWithPred(pathCountsWithPredMatrix *pathCounts) {
    free(pathCounts->data);
    free(pathCounts);
}
