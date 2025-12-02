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
#include "signal.h"
#include "scores.h"
#include "mem.h"


/*
  Private function: return sqrt[sum of square of (diff between same-index elements)] 
*/
static SIGNALTYPE calculateScoresDiff(signalMatrix *sumOfSignal, size_t nbNodes);


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
    size_t sumDegrees = interactomeComp->offsets[nbGenes];

    // calculate normalization factors (used in each iteration)
    fprintf(stderr, "INFO gbaCentrality(): calculating normalization factors\n");
    normFactorVector *normFactVec = buildNormFactorVector(interactomeComp, alpha);
    for (size_t i = 0; i < sumDegrees; i++)
        fprintf(stderr, "norm factor = %f\n", normFactVec->data[i]);
    
    fprintf(stderr, "INFO gbaCentrality(): calculating B_1\n");
    signalWithPredMatrix *signalCurrent = buildFirstSignal(interactomeComp, normFactVec);
    signalWithPredMatrix *signalNext = NULL;

    size_t k = 1;
    // for convergence test
    SCORETYPE threshold = 1E-4;
    SCORETYPE scoresDiff = 1;

    while (scoresDiff > threshold) {
        // scores += causal * B_k
        signalMatrix *sumOfSignal = signalSum(signalCurrent, interactomeComp);

        for (size_t j = 0; j < nbGenes; j++) {
            for (size_t i = 0; i < nbGenes; i++) {
                // updates scores: vector * B_k
                scores->scores[j] += sumOfSignal->data[i* sumDegrees + j] * causal->scores[i];
            }
        }

        // calculate the difference between scores for B_k-1 and B_k,
        scoresDiff = calculateScoresDiff(sumOfSignal, nbGenes);
        fprintf(stderr, "INFO gbaCentrality(): scoresDiff = %f\n", scoresDiff);

        if (scoresDiff > threshold) {
            // build B_(k+1) for next iteration
            fprintf(stderr, "INFO gbaCentrality(): calculating B_%ld\n", k+1);
            signalNext = buildNextSignal(signalCurrent, sumOfSignal, interactomeComp, normFactVec);
            freeSignalWithPred(signalCurrent);
            signalCurrent = signalNext;
            k++;
        }
        freeSignal(sumOfSignal);
    }
    freeNormFactorVector(normFactVec);
    freeSignalWithPred(signalCurrent);
    freeCompactAdjacency(interactomeComp);
}

/*
  Return L2-norm of sumOfSignal
*/
static SIGNALTYPE calculateScoresDiff(signalMatrix *sumOfSignal, size_t nbNodes) {
    SIGNALTYPE scoresDiff = 0; // double == high precision for the running sum
    for (size_t i = 0; i < nbNodes * nbNodes; i++) {
        scoresDiff += sumOfSignal->data[i] * sumOfSignal->data[i];
    }
    scoresDiff = sqrt(scoresDiff);
    return(scoresDiff);
}
