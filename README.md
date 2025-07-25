# GBA-centrality-C

This repository contains the code for GBA-centrality.

run `make` to generate a library gbaCentrality.so

It contains the gbaCentrality() function, defined in gbaCentrality.h

It takes as input:
- pointer to memory allocated for a adjacencyMatrix structure as defined in adjacency.h
- pointer to memory allocated for a geneScores structure (seed scores) as defined in scores.h
- alpha (attenuation parameter)
- pointer to memory allocated for a geneScores structure (for final scores)

Required dependencies:
- bear
