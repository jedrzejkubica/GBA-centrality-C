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

#ifndef _SCORES_H_
#define _SCORES_H_


// type used for each geneScores->scores value
#define SCORETYPE double

/*
    Store one score per gene, each score is a float >=0.
    scores MUST be large enough to store nbGenes floats.
 */
typedef struct {
    unsigned int nbGenes;
    SCORETYPE *scores;
} geneScores;


void printScores(geneScores *scores);

void freeScores(geneScores *scores);

#endif
