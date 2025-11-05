# GBA-centrality-C

This repository contains the C code for the GBA centrality scoring algorithm.

It does the heavy lifting behind our [python GBA centrality](https://github.com/jedrzejkubica/GBA-centrality) software.

It can also be used stand-alone by calling the gbaCentrality() function in the
shared library gbaCentrality.so.


### How to use this repository

```
git clone https://github.com/jedrzejkubica/GBA-centrality-C.git
cd GBA-centrality-C
make all
```

This generates an executable testAdjacency and a library gbaCentrality.so. The library
contains the gbaCentrality() function (see gbaCentrality.h), which takes as input:
- pointer to a network structure defining your network, ie graph (see network.h)
- pointer to a geneScores structure defining your causal genes, ie seeds (see scores.h)
- alpha attenuation parameter, must be in ]0,1[, typically 0.5 is good
- pointer to allocated memory for another geneScores structure, the content will be
  filled by gbaCentrality() with the GBA-centrality scores


### Dependencies:
- None for production.
- For development you should install "bear", which allows to build a compilation database
  for LSP (usable by emacs and other IDEs). Then call "make" rather than "make all".
