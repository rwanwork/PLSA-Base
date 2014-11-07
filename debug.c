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
#include <limits.h>                                  /*  UINT_MAX  */
#include <stdbool.h>
#include <math.h>
#include <signal.h>
#include <time.h>

#include "wmalloc.h"
#include "plsa-defn.h"
#include "debug.h"

void handler_sigfpe () {
  signal (SIGFPE, handler_sigfpe); /* reset signal */
  fprintf (stderr, "-->  I have received a SIGFPE!\n");

  return;
}


/*  Used by input.c  */
void debugCheckCo (INFO *info) {
  unsigned int j = 0;
  unsigned int cos_count = 0;
  unsigned int curr_j;

  /*  Check flags in co-occurrence table  */
  for (unsigned int i = 0; i < info -> m; i++) {
    cos_count = GET_COS_POSITION (i, 0);
    curr_j = 0;
    for (unsigned int pos_j = 1; pos_j <= cos_count; pos_j++) {
      j = GET_COS_POSITION (i, pos_j);
      while (curr_j < j) {
        fprintf (stderr, "X");
        curr_j++;
      }
      fprintf (stderr, "O");
    }

    while (curr_j < info -> n) {
      fprintf (stderr, "X");
      curr_j++;
    }

    fprintf (stderr, "\n");
  }

  fprintf (stderr, "-----\n");

  return;
}


/*  Used by output.c  */
void checkCoProb (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  unsigned int i = 0;  /*  Index into w1  */
  unsigned int j = 0;  /*  Index into w2  */
  unsigned int k = 0;  /*  Index into clusters  */
  PROBNODE temp;
  PROBNODE tempsum = 0.0;
  unsigned int nonprob = 0;

  for (i = 0; i < info -> m; i++) {
    for (j = 0; j < info -> n; j++) {
      temp = info -> probz[0] + GET_PROBW1_Z (0, i) + GET_PROBW2_Z (0, j);
      for (k = 1; k < num_clusters; k++) {
        /*  temp stores logarithms  */
        logSumsInline (temp, (GET_PROBZ (k) + GET_PROBW1_Z (k, i) + GET_PROBW2_Z (k, j)));
      }
      if (temp > 0) {
        nonprob++;
      }
      tempsum += DOEXP (temp);
    }
  }

  if (info -> verbose) {
    fprintf (stderr, "**\t%u : %f\n", nonprob, tempsum);
  }

  return;
}

