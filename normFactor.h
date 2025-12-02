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

#ifndef _NORMFACTOR_H_
#define _NORMFACTOR_H_

#include "compactAdjacency.h"

/*
  For each edge i->j defined by offset (in a compactAdjacency), store
  alpha * w(i->j) / sumOfWeights(edges incoming to j)
  in normFactorVector->data[offset]
*/
typedef struct {
    float *data;
} normFactorVector;


/*
  Build and return the normFactorVector corresponding to the provided compact
*/
normFactorVector *buildNormFactorVector(compactAdjacencyMatrix *compact, float alpha);


void freeNormFactorVector(normFactorVector *normFactVec);

#endif