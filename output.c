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
#include <limits.h>                                  /*  UINT_MAX  */
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "wmalloc.h"
#include "plsa-defn.h"
#include "output.h"


void printCoProb (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  unsigned int i = 0;  /*  Index into w1  */
  unsigned int j = 0;  /*  Index into w2  */
  unsigned int k = 0;  /*  Index into clusters  */
  PROBNODE temp;
  PROBNODE tempsum = 0.0;
  unsigned int nonprob = 0;
  FILE *fp = NULL;
  char *fn = wmalloc (sizeof (char) * (strlen (info -> base_fn) + 10));
  static unsigned int snapshot_count = 0;

  time_t start;
  time_t end;

  time (&start);
  snapshot_count++;

  sprintf (fn, "%s.plsa", info -> base_fn);
  if (info -> textio) {
    FOPEN (fn, fp, "w");
    fprintf (fp, "%u\t", info -> m);
    fprintf (fp, "%u\t", info -> n);
    for (i = 0; i < info -> m; i++) {
      fprintf (fp, "%u\t", info -> row_ids[i]);
    }
    for (j = 0; j < info -> n; j++) {
      fprintf (fp, "%u\t", info -> column_ids[j]);
    }
  }
  else {
    FOPEN (fn, fp, "wb");
    fwrite (&info -> m, sizeof (unsigned int), 1, fp);
    fwrite (&info -> n, sizeof (unsigned int), 1, fp);
    fwrite (info -> row_ids, sizeof (unsigned int), info -> m, fp);
    fwrite (info -> column_ids, sizeof (unsigned int), info -> n, fp);
  }

  for (i = 0; i < info -> m; i++) {
    for (j = 0; j < info -> n; j++) {
      temp = info -> probz[0] + GET_PROBW1_Z (0, i) + GET_PROBW2_Z (0, j);
      for (k = 1; k < num_clusters; k++) {
        /*  temp stores logarithms  */
        logSumsInline (temp, (GET_PROBZ (k) + GET_PROBW1_Z (k, i) + GET_PROBW2_Z (k, j)));
      }

      /*  temp stores logarithms; round it to ROUND_DIGITS  */
      if (info -> rounding) {
        temp = (round (temp * ROUND_DIGITS)) / ROUND_DIGITS;
      }

      /*  temp stores logarithms  */
      if (info -> textio) {
        fprintf (fp, "%lf\t", temp);
      }
      else {
        fwrite (&temp, sizeof (PROBNODE), 1, fp);
      }
      if (temp > 0) {
        nonprob++;
      }
      tempsum += DOEXP (temp);
    }
  }

  FCLOSE (fp);
  wfree (fn);

  if ((info -> verbose) && (info -> iter == UINT_MAX)) {
    fprintf (stderr, "==\tNon-probabilities:                              %u\n", nonprob);
    fprintf (stderr, "==\tSum of p(x,y):                                  %f\n", tempsum);
    fprintf (stderr, "==\tTotal output files printed                      %u\n", snapshot_count);
  }

  time (&end);
  info -> printCoProbs_time += difftime (end, start);

  return;
}


