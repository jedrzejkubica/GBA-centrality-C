# GBA-centrality-C

This repository contains the C code for the GBA centrality scoring algorithm. It does the heavy lifting behind our [python GBA centrality](https://github.com/jedrzejkubica/GBA-centrality) software, or can be used stand-alone.


## Install GBA-centrality-C

```
git clone https://github.com/jedrzejkubica/GBA-centrality-C.git
cd GBA-centrality-C
make
```

This generates an executable `testAdjacency` and a shared library `gbaCentrality.so`.

To make sure the installation is correct, run `testAdjacency`, which computes and prints GBA centrality scores for four example networks.

## Use GBA-centrality-C

GBA-centrality-C can be used stand-alone by calling the `gbaCentrality()` function in the shared library `gbaCentrality.so`. The function arguments are as follows (for details see `gbaCentrality.h`):
- pointer to a network structure defining a network (see network.h)
- pointer to a geneScores structure defining seeds (see scores.h)
- alpha attenuation parameter (0 < alpha < 1, typically 0.5 is good)
- pointer to allocated memory for another geneScores structure, the content will be filled with GBA centrality scores


## Dependencies

- None for production
- For development you can install "bear", which allows to build a compilation database for LSP (usable by emacs and other IDEs). Then call `make all` instead of `make`
