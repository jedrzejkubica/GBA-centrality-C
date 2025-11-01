# GBA-centrality-C

This repository contains the code for the GBA centrality scoring algorithm.

### How to use this repository

```
git clone https://github.com/jedrzejkubica/GBA-centrality-C.git
cd GBA-centrality-C
make
```

This generates an executable testAdjacency and a library gbaCentrality.so. The library contains the gbaCentrality() function, which takes as input:
- pointer to memory allocated for an adjacencyMatrix structure (as defined in adjacency.h)
- pointer to memory allocated for a geneScores structure (as defined in scores.h)
- alpha (attenuation parameter set by user)
- pointer to memory allocated for a geneScores structure (to allocate final scores)

### Dependencies:

- bear
