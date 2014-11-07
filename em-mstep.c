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
#include "em-mstep.h"

void applyMStep (INFO *info) {
  register unsigned int i;  /*  Index into w1  */
  register unsigned int j;  /*  Index into w2  */
  register unsigned int k;  /*  Index into clusters  */
  register unsigned int pos_j;  /*  Actual position in the cooccurrence array  */
  register unsigned int cos_count;  /*  Number of cooccurrences in each row  */
  register PROBNODE cos;
  bool *flag_z = NULL;
  bool **flag_w1_z = NULL;
  bool **flag_w2_z = NULL;

  time_t start;
  time_t end;

  time (&start);

  /*******************************************************/
  /*  Initialize flags that indicate whether the cell is so far untouched   */

  flag_z = wmalloc (info -> block_size * sizeof (bool));
  for (k = 0; k < info -> block_size; k++) {
    flag_z[k] = false;
  }

  flag_w1_z = wmalloc (info -> block_size * sizeof (bool*));
  for (k = 0; k < info -> block_size; k++) {
    flag_w1_z[k] = wmalloc (info -> m * sizeof (bool));
    for (i = 0; i < info -> m; i++) {
      flag_w1_z[k][i] = false;
    }
  }

  flag_w2_z = wmalloc (info -> block_size * sizeof (bool*));
  for (k = 0; k < info -> block_size; k++) {
    flag_w2_z[k] = wmalloc (info -> n * sizeof (bool));
    for (j = 0; j < info -> n; j++) {
      flag_w2_z[k][j] = false;
    }
  }

  /*******************************************************/
  /*  Update probabilities */

  for (k = 0; k < info -> block_size; k++) {
    for (i = 0; i < info -> m; i++) {
      cos_count = GET_COS_POSITION (i, 0);
      for (pos_j = 1; pos_j <= cos_count; pos_j++) {
        j = GET_COS_POSITION (i, pos_j);
        cos = GET_COS (i, pos_j);

        /*  probz  */
        if (flag_z[k]) {
          logSumsInline (GET_PROBZ (k), cos + GET_PROBZ_W1W2 (k, i, j));
        }
        else {
          GET_PROBZ (k) = cos + GET_PROBZ_W1W2 (k, i, j);
          flag_z[k] = true;
        }

        /*  probw1_z  */
        if (flag_w1_z[k][i]) {
          logSumsInline (GET_PROBW1_Z (k, i), cos + GET_PROBZ_W1W2 (k, i, j));
        }
        else {
          GET_PROBW1_Z (k, i) = cos + GET_PROBZ_W1W2 (k, i, j);
          flag_w1_z[k][i] = true;
        }

        /*  probw2_z  */
        if (flag_w2_z[k][j]) {
          logSumsInline (GET_PROBW2_Z (k, j), cos + GET_PROBZ_W1W2 (k, i, j));
        }
        else {
          GET_PROBW2_Z (k, j) = cos + GET_PROBZ_W1W2 (k, i, j);
          flag_w2_z[k][j] = true;
        }
      }
    }
  }

  /*******************************************************/
  /*  Uninitialize flags  */

  wfree (flag_z);
  for (k = 0; k < info -> block_size; k++) {
    wfree (flag_w1_z[k]);
    wfree (flag_w2_z[k]);
  }
  wfree (flag_w1_z);
  wfree (flag_w2_z);

  time (&end);
  info -> applyMStep_time += difftime (end, start);

  return;
}


/*!  Normalize probabilities  */
void normalizeProbs (INFO *info) {
  unsigned int i;  /*  Index into w1  */
  unsigned int j;  /*  Index into w2  */
  unsigned int k;  /*  Index into clusters  */
  PROBNODE sum;
  PROBNODE norm;
  time_t start;
  time_t end;

  time (&start);

  for (k = 0; k < info -> num_clusters; k++) {
    norm = GET_PROBZ (k);

    /*  probw1_z  */
    for (i = 0; i < info -> m; i++) {
      GET_PROBW1_Z (k, i) = GET_PROBW1_Z (k, i) - norm;
    }

    /*  probw2_z  */
    for (j = 0; j < info -> n; j++) {
      GET_PROBW2_Z (k, j) = GET_PROBW2_Z (k, j) - norm;
    }
  }

  /*  probz  */
  sum = GET_PROBZ (0);
  for (k = 1; k < info -> num_clusters; k++) {
//     fprintf (stderr, "\t%f\n", GET_PROBZ (k));
    logSumsInline (sum, GET_PROBZ (k));
  }
  for (k = 0; k < info -> num_clusters; k++) {
    GET_PROBZ (k) = GET_PROBZ (k) - sum;
  }

  time (&end);
  info -> normalizeProbs_time += difftime (end, start);

  return;
}


