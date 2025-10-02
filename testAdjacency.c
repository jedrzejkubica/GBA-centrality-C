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

#include "network.h"
#include "scores.h"
#include "gbaCentrality.h"
#include "mem.h"


network *diamond4(void) {
    // diamond network for testing
    network *diamond4 = mallocOrDie(sizeof(network), "E: OOM for diamond4");
    diamond4->nbNodes = 4;
    diamond4->nbEdges = 9;
    diamond4->edges = mallocOrDie(diamond4->nbEdges * sizeof(edge), "E: OOM for diamond4 edges");
    edge *edgeP = diamond4->edges;

    edgeP->source = 0;
    edgeP->dest = 1;
    edgeP->weight = 0.2;
    edgeP++;
    edgeP->source = 0;
    edgeP->dest = 2;
    edgeP->weight = 0.5;
    edgeP++;

    edgeP->source = 1;
    edgeP->dest = 0;
    edgeP->weight = 0.4;
    edgeP++;
    edgeP->source = 1;
    edgeP->dest = 3;
    edgeP->weight = 1;
    edgeP++;

    edgeP->source = 2;
    edgeP->dest = 0;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 2;
    edgeP->dest = 3;
    edgeP->weight = 1;
    edgeP++;

    edgeP->source = 3;
    edgeP->dest = 1;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 3;
    edgeP->dest = 2;
    edgeP->weight = 1;
    edgeP++;
    
    // add a self-interaction 3->3 for testing
    edgeP->source = 3;
    edgeP->dest = 3;
    edgeP->weight = 1;
    return(diamond4);
}

geneScores *causalGenesDiam4(void) {
    geneScores *causal = mallocOrDie(sizeof(geneScores), "E: OOM for causal genes");
    causal->nbGenes = 4;
    causal->scores = mallocOrDie(causal->nbGenes * sizeof(SCORETYPE), "E: OOM for causal genes");

    for (size_t i = 0; i < causal->nbGenes; i++)
        causal->scores[i] = 0;

    //causal->scores[1] = 1;
    causal->scores[2] = 1;

    return(causal);
}


/*
  7 nodes: node 1 is connected to nodes 0 and 2, node 2 is a hub connected to 1 and 3-6,
  and we create this as an undirected unweighted network
*/
network *asymmetric(void) {
    network *asym = mallocOrDie(sizeof(network), "E: OOM for asym");
    asym->nbNodes = 7;
    asym->nbEdges = 12;
    asym->edges = mallocOrDie(asym->nbEdges * sizeof(edge), "E: OOM for asym edges");
    edge *edgeP = asym->edges;

    edgeP->source = 0;
    edgeP->dest = 1;
    edgeP->weight = 1;
    edgeP++;
    
    edgeP->source = 1;
    edgeP->dest = 0;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 1;
    edgeP->dest = 2;
    edgeP->weight = 1;
    edgeP++;

    edgeP->source = 2;
    edgeP->dest = 1;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 2;
    edgeP->dest = 3;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 2;
    edgeP->dest = 4;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 2;
    edgeP->dest = 5;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 2;
    edgeP->dest = 6;
    edgeP->weight = 1;
    edgeP++;

    edgeP->source = 3;
    edgeP->dest = 2;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 4;
    edgeP->dest = 2;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 5;
    edgeP->dest = 2;
    edgeP->weight = 1;
    edgeP++;
    edgeP->source = 6;
    edgeP->dest = 2;
    edgeP->weight = 1;
    
    return(asym);
}

geneScores *causalGenesAsym(void) {
    geneScores *causal = mallocOrDie(sizeof(geneScores), "E: OOM for causal genes");
    causal->nbGenes = 7;
    causal->scores = mallocOrDie(causal->nbGenes * sizeof(SCORETYPE), "E: OOM for causal genes");

    for (size_t i = 0; i < causal->nbGenes; i++)
        causal->scores[i] = 0;

    //causal->scores[1] = 1;
    causal->scores[2] = 1;

    return(causal);
}

int main(void) {
    {
        // dimaond4 network
        network *diam4 = diamond4();
        printf("diam4 with its self-loop\n");
        printNetwork(diam4);
    
        geneScores *causal = causalGenesDiam4();
        printf("Diam4 causal genes:\n");
        printScores(causal);

        geneScores *result = mallocOrDie(sizeof(geneScores), "E: OOM for result scores");
        result->nbGenes = causal->nbGenes;
        result->scores = mallocOrDie(causal->nbGenes * sizeof(SCORETYPE), "E: OOM for result scores");

        gbaCentrality(diam4, causal, 0.5, result);
        printf("Diam4 scores\n");
        printScores(result);
        freeNetwork(diam4);
        freeScores(causal);
        freeScores(result);
    }

    {
        // asymmetric network
        network *asym = asymmetric();
        printf("asymmetric with one hub\n");
        printNetwork(asym);
    
        geneScores *causal = causalGenesAsym();
        printf("Asym causal genes:\n");
        printScores(causal);

        geneScores *result = mallocOrDie(sizeof(geneScores), "E: OOM for result scores");
        result->nbGenes = causal->nbGenes;
        result->scores = mallocOrDie(causal->nbGenes * sizeof(SCORETYPE), "E: OOM for result scores");

        gbaCentrality(asym, causal, 0.5, result);
        printf("Asym scores\n");
        printScores(result);
        freeNetwork(asym);
        freeScores(causal);
        freeScores(result);
    }

    return(0);
}

