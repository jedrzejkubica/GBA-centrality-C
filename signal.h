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

#ifndef _SIGNAL_H_
#define _SIGNAL_H_

#include <stddef.h>
#include "compactAdjacency.h"


/*
  type to store the sum of path weights of a given length, floating point number
  of some precision (eg maybe double)
*/
#define SIGNALTYPE float


/*
    this structure is dependent on a compactAdjacencyMatrix;
    data is of size offsets[nbNodes]*nbNodes;
    for any k in [0, inDeg(j)-1]:
    data[i*offsets[nbNodes] + offsets[j] + k]
    is the "signal" that is propagated from node i to node j via the
    k'th predecessor of node j, ie. predecessors[offsets[j] + k].

    Ideally this "signal" would be propagated along all "paths" from i to j,
    normalized by the in-degree of j, and multiuplying by an attenuation factor
    alpha at each step. However counting paths is hard, so instead we propagate
    along walks but avoiding backtracks and loops back to the starting node i.
*/
typedef struct {
    SIGNALTYPE *data;
} signalWithPredMatrix;


/*
    data[i*nbCols + j] contains the number of paths
    (or the sum of path weights for a weighted network)
    of a given length from node i to j
*/
typedef struct {
    size_t nbCols;
    SIGNALTYPE *data;
} signalMatrix;


/*
    build signalWithPredMatrix for paths of length 1,
    return a pointer to a freshly allocated structure;
*/
signalWithPredMatrix *buildFirstSignal(compactAdjacencyMatrix *compact);

/*
    build a signalWithPredMatrix for paths of length k+1
    given signalWithPredMatrix for paths of length k
    (excluding paths looping back to the starting node),
    return a pointer to a freshly allocated structure;
*/
signalWithPredMatrix *buildNextSignal(signalWithPredMatrix *signalWithPred, signalMatrix *signal,
                                              compactAdjacencyMatrix *compact);

void freeSignalWithPred(signalWithPredMatrix *signal);


/*  
    produce signalMatrix corresponding to signalWithPredMatrix
*/
signalMatrix *countPaths(signalWithPredMatrix *signalWithPred, compactAdjacencyMatrix *compact);

void printSignal(signalMatrix *signal);

void freeSignal(signalMatrix *signal);

#endif
