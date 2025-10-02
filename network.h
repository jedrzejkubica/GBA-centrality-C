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
