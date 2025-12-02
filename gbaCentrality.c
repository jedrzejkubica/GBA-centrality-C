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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gbaCentrality.h"
#include "network.h"
#include "compactAdjacency.h"
#include "pathCounts.h"
#include "scores.h"
#include "mem.h"


/*
  Return sqrt[sum of square of (diff between same-index elements)] 
*/
SCORETYPE calculateScoresDiff(geneScores *scores, geneScores *scoresPrev);


void gbaCentrality(network *N, geneScores *causal, float alpha, geneScores *scores) {
    // sanity check:
    if (N->nbNodes != causal->nbGenes) {
        fprintf(stderr, "ERROR: gbaCentrality() called with network and causal genes of different sizes");
        exit(1);
    }
    size_t nbGenes = causal->nbGenes;

    // start by copying causal scores, ie scores = alpha**0 * I * casual
    memcpy(scores->scores, causal->scores, nbGenes * sizeof(SCORETYPE));

    compactAdjacencyMatrix *interactomeComp = network2compact(N);
   
    fprintf(stderr, "INFO gbaCentrality(): calculating B_1\n");

    pathCountsWithPredMatrix *pathCountsCurrent = buildFirstPathCounts(interactomeComp);
    pathCountsWithPredMatrix *pathCountsNext = NULL;

    geneScores *scoresPrev = mallocOrDie(sizeof(geneScores), "ERROR: OOM for scoresPrev\n");
    scoresPrev->nbGenes = nbGenes;
    scoresPrev->scores = mallocOrDie(sizeof(SCORETYPE) * nbGenes, "ERROR: OOM for scoresPrev scores");

    float alphaPowK = 1;
    size_t k = 1;
    
    // for convergence test
    SCORETYPE threshold = 1E-4;
    SCORETYPE scoresDiff = 1;

    while (scoresDiff > threshold) {
        // save scores from B_k-1
        memcpy(scoresPrev->scores, scores->scores, nbGenes * sizeof(SCORETYPE));

        // scores += alpha**k * causal * B_k
        alphaPowK *= alpha;
        pathCountsMatrix *interactomePathCounts = countPaths(pathCountsCurrent, interactomeComp);

        for (size_t j = 0; j < nbGenes; j++) {
            // calculate sum of column j
            double colSum = 0;
            for (size_t i = 0; i < nbGenes; i++) {
                colSum += interactomePathCounts->data[i * nbGenes + j];
            }
            // updates scores
            if (colSum != 0) {
                double scoreSum = 0;
                for (size_t i = 0; i < nbGenes; i++) {
                    scoreSum += interactomePathCounts->data[i * nbGenes + j] * causal->scores[i];
                }
                scores->scores[j] += alphaPowK * scoreSum / colSum;
            }
        }

        // calculate the difference between scores for B_k-1 and B_k
        scoresDiff = calculateScoresDiff(scores, scoresPrev);
        fprintf(stderr, "INFO gbaCentrality(): scoresDiff = %f\n", scoresDiff);

        if (scoresDiff > threshold) {
            // build B_(k+1) for next iteration
            fprintf(stderr, "INFO gbaCentrality(): calculating B_%ld\n", k+1);
            pathCountsNext = buildNextPathCounts(pathCountsCurrent, interactomePathCounts, interactomeComp);
            freePathCountsWithPred(pathCountsCurrent);
            pathCountsCurrent = pathCountsNext;
            k++;
        }
        freePathCounts(interactomePathCounts);
    }
    freeScores(scoresPrev);
    freePathCountsWithPred(pathCountsCurrent);
    freeCompactAdjacency(interactomeComp);
}

/*
  Return L2-norm of (scores - scoresPrev)
*/
SCORETYPE calculateScoresDiff(geneScores *scores, geneScores *scoresPrev) {
    double scoresDiff = 0; // double == high precision for the running sum
    for (size_t i = 0; i < scores->nbGenes; i++) {
        double diff = scores->scores[i] - scoresPrev->scores[i];
        scoresDiff += diff * diff;
    }
    scoresDiff = sqrt(scoresDiff);
    return(scoresDiff);
}
