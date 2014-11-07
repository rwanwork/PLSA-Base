/*
**  Probabilistic latent semantic analysis (PLSA, baseline version)
**  Copyright (C) 2009-2010  by Raymond Wan (r.wan@aist.go.jp)
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>  /*  log10 function  */
#include <stdbool.h>
#include <time.h>

#include "wmalloc.h"
#include "plsa-defn.h"
#include "em-estep.h"

void initEM (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  register unsigned int i;  /*  Index into w1  */
  register unsigned int j;  /*  Index into w2  */
  register unsigned int k;  /*  Index into clusters  */
  register PROBNODE sum;
  time_t start;
  time_t end;

  time (&start);
  PROGRESS_MSG ("Begin initialization...");

  /*  Assign probabilities to probz  */
  sum = 0.0;
  for (k = 0; k < num_clusters; k++) {
    GET_PROBZ (k) = RANDOM_FLOAT;
    sum += GET_PROBZ (k);
  }
  for (k = 0; k < num_clusters; k++) {
    GET_PROBZ (k) = DOLOG (GET_PROBZ (k) / sum);
  }

  /*  Assign probabilities to probw1_z  */
  for (i = 0; i < (info -> num_clusters * info -> m); i++) {
    info -> probw1_z[i] = RANDOM_FLOAT;
  }
  for (k = 0; k < num_clusters; k++) {
    sum = 0.0;
    for (i = 0; i < info -> m; i++) {
      sum += GET_PROBW1_Z (k, i);
    }
    for (i = 0; i < info -> m; i++) {
      GET_PROBW1_Z (k, i) = DOLOG (GET_PROBW1_Z (k, i) / sum);
    }
  }

  /*  Assign probabilities to probw2_z  */
  for (i = 0; i < (info -> num_clusters * info -> n); i++) {
    info -> probw2_z[i] = RANDOM_FLOAT;
  }
  for (k = 0; k < num_clusters; k++) {
    sum= 0.0;
    for (j = 0; j < info -> n; j++) {
      sum += GET_PROBW2_Z (k, j);
    }
    for (j = 0; j < info -> n; j++) {
      GET_PROBW2_Z (k, j) = DOLOG (GET_PROBW2_Z (k, j) / sum);
    }
  }

  PROGRESS_MSG ("Initialization complete...");
  time (&end);
  info -> initEM_time += difftime (end, start);

  return;
}


void applyEStep (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  unsigned int i = 0;  /*  Index into w1  */
  unsigned int j = 0;  /*  Index into w2  */
  unsigned int k = 0;  /*  Index into clusters  */
  PROBNODE sum = 0.0;
  time_t start;
  time_t end;

  time (&start);
  for (i = 0; i < info -> m; i++) {
    for (j = 0; j < info -> n; j++) {
      GET_PROBZ_W1W2 (0, i, j) = GET_PROBW1_Z (0, i) + GET_PROBW2_Z (0, j) + GET_PROBZ (0);
      sum = GET_PROBZ_W1W2 (0, i, j);
      for (k = 1; k < num_clusters; k++) {
        GET_PROBZ_W1W2 (k, i, j) = GET_PROBW1_Z (k, i) + GET_PROBW2_Z (k, j) + GET_PROBZ (k);
        logSumsInline (sum, GET_PROBZ_W1W2 (k, i, j));
      }

      /*  Divide through by the denominator  */
      for (k = 0; k < num_clusters; k++) {
        GET_PROBZ_W1W2 (k, i, j) = GET_PROBZ_W1W2 (k, i, j) - sum;
      }
    }
  }

  time (&end);
  info -> applyEStep_time += difftime (end, start);

  return;
}


PROBNODE calculateML (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  register unsigned int i;  /*  Index into w1  */
  register unsigned int j;  /*  Index into w2  */
  register unsigned int k;  /*  Index into clusters  */
  register unsigned int pos_j;  /*  Actual position in the cooccurrence array  */
  register unsigned int cos_count;  /*  Number of cooccurrences in each row  */
  register PROBNODE total = 0.0;
  PROBNODE temp;
  unsigned int count = 0;
  time_t start;
  time_t end;

  time (&start);

  for (i = 0; i < info -> m; i++) {
    cos_count = GET_COS_POSITION (i, 0);
    for (pos_j = 1; pos_j <= cos_count; pos_j++) {
      j = GET_COS_POSITION (i, pos_j);

      /*  Initialize with cluster 0  */
      temp = (GET_PROBZ (0) + GET_PROBW1_Z (0, i) + GET_PROBW2_Z (0, j));
      /*  Log-likelihood for the co-occurrence of two words  */
      for (k = 1; k < num_clusters; k++) {
        /*  temp stores log values  */
        logSumsInline (temp, GET_PROBZ (k) + GET_PROBW1_Z (k, i) + GET_PROBW2_Z (k, j));
      }

      /*  Log-likelihood across all examples  */
      total += (temp * DOEXP (GET_COS (i, pos_j)));
      count++;
    }
  }

  time (&end);
  info -> calculateML_time += difftime (end, start);

  return (total);
}


