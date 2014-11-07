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
#include <time.h>  /*  time  */

#include "wmalloc.h"
#include "plsa-defn.h"
#include "debug.h"
#include "input.h"


/*!  Initialization that depends on the input file or parameters  */
void initializePostInput (INFO *info) {
  unsigned int i = 0;
  unsigned int size = 0;
  unsigned int temp = 0;

  info -> cos = wmalloc (info -> m * sizeof (COOCCUR*));

  /*  All processes read the co-occurrence data, so all must initialize  */
  /*  Cannot allocate more space since we don't know the number of values in each row  */
  for (i = 0; i < info -> m; i++) {
    info -> cos[i] = NULL;
  }

  size = info -> num_clusters;

  /*  Allocate space  */
  info -> probw1_z = wmalloc (size * info -> m * sizeof (PROBNODE));
  info -> probw2_z = wmalloc (size * info -> n * sizeof (PROBNODE));
  info -> probz = wmalloc (size * sizeof (PROBNODE));
  info -> probz_w1w2 = wmalloc (size * sizeof (PROBNODE*));
  for (i = 0; i < size; i++) {
    info -> probz_w1w2[i] = wmalloc (info -> m * info -> n * sizeof (PROBNODE));
  }

  /*  Set seed if given as an argument, otherwise use the time  */
  if (info -> seed == UINT_MAX) {
    temp = time (NULL);
    info -> seed = temp;
    if (info -> verbose) {
      fprintf (stderr, "==\tApplying seed from time:                        %u\n", temp);
    }
  }
  srand (info -> seed);

  return;
}


/*!
**  Read the co-occurrence data from file.  The format of the file is:
**
**  [rows][columns][row id+][column id+][w1 cos_count (w21 c21) ... (w2n c2n)]+**
**
**  row and column ids are integer values that map to the original
**  vocabulary.  The number of values should be (info -> m) and
**  (info -> n), respectively.
**
**  Every value is an unsigned integer in binary format, unless
**  textmode is TRUE -- if so, values are in text, separated
**  by white space (tab).
**
**  Note:  i indexes for rows (w1); j indexes for columns (w2)
*/
bool readCO (INFO *info) {
  FILE *fp = NULL;
  unsigned int w1 = 0;
  unsigned int w2 = 0;
  unsigned int freq = 0;

  unsigned int rows = 0;
  unsigned int cols = 0;
  unsigned int cos_count = 0;
  unsigned int found_pairs = 0;
  unsigned int found_w1 = 0;

  unsigned int sum_freq = 0;
  unsigned int nonzero_count = 0;
  time_t start;
  time_t end;

  time (&start);

  PROGRESS_MSG ("Reading from co-occurrence file...");

  /*  Open the file; read the number of rows and columns and check them  */
  if (info -> textio) {
    FOPEN (info -> co_fn, fp, "r");
    fscanf (fp, "%u", &rows);
    fscanf (fp, "%u", &cols);
  }
  else {
    FOPEN (info -> co_fn, fp, "rb");
    fread (&rows, sizeof (unsigned int), 1, fp);
    fread (&cols, sizeof (unsigned int), 1, fp);
  }

  info -> m = rows;
  info -> n = cols;

  initializePostInput (info);

  info -> row_ids = wmalloc (info -> m * sizeof (unsigned int));
  info -> column_ids = wmalloc (info -> n * sizeof (unsigned int));

  if (info -> textio) {
    for (unsigned int i = 0; i < info -> m; i++) {
      fscanf (fp, "%u", &(info -> row_ids[i]));
    }
    for (unsigned int j = 0; j < info -> n; j++) {
      fscanf (fp, "%u", &(info -> column_ids[j]));
    }
  }
  else {
    fread (info -> row_ids, sizeof (unsigned int), info -> m, fp);
    fread (info -> column_ids, sizeof (unsigned int), info -> n, fp);
  }

  found_pairs = 0;
  found_w1 = 0;
  for (unsigned int i = 0; i < info -> m; i++) {
    if (info -> textio) {
      fscanf (fp, "%u", &w1);
    }
    else {
      fread (&w1, sizeof (unsigned int), 1, fp);
    }

    if (feof (fp)) {
      break;
    }
    found_w1++;
    if (info -> textio) {
      fscanf (fp, "%u", &cos_count);
    }
    else {
      fread (&cos_count, sizeof (unsigned int), 1, fp);
    }

    /*  Allocate space for the row  */
    info -> cos[i] = wmalloc (sizeof (COOCCUR) * (cos_count + 1));

    /*  Position 0 of each row is cos_count    */
    info ->  cos[i][0].x = 0.0;
    info ->  cos[i][0].column = cos_count;

    /*  Term found is a query term  */
    for (unsigned int j = 1; j <= cos_count; j++) {
      if (info -> textio) {
        fscanf (fp, "%u", &w2);
        fscanf (fp, "%u", &freq);
      }
      else {
        fread (&w2, sizeof (unsigned int), 1, fp);
        fread (&freq, sizeof (unsigned int), 1, fp);
      }

      if (freq != 0) {
        nonzero_count++;
      }

      if (w2 > info -> n) {
        fprintf (stderr, "Word 2 (%u) is out of range (%u).\n", w2, info -> n);
        exit (EXIT_FAILURE);
      }

      if (info -> debug) {
        fprintf (stderr, "==\t\tRead (%u, %u) --> %u\n", i, w2, freq);
      }

      SET_COS (i, j, w2, DOLOG (freq));

      sum_freq += freq;

      found_pairs++;
    }
  }
  FCLOSE (fp);

  /*  Check if the header of the file matches reality  */
  if (found_w1 != info -> m) {
    fprintf (stderr, "Not all query terms found!  (%u, %u)\n", found_w1, info -> m);
    exit (EXIT_FAILURE);
  }

  if (info -> verbose) {
    unsigned int zero_count = (info -> m * info -> n) - nonzero_count;
    fprintf (stderr, "==\tMaximum number of pairs:                        %u\n", info -> m * info -> n);
    fprintf (stderr, "==\tActual number of pairs in data file:            %u\n", found_pairs);
    fprintf (stderr, "==\tPercentage of zeroes:                           %.2f %% (%u)\n", (double) zero_count / (double) ((info -> m * info -> n)) * 100, zero_count);
    fprintf (stderr, "==\tSum of co-occurrence counts:                    %u\n", sum_freq);
  }

#if DEBUG
    debugCheckCo (info);
#endif

  time (&end);
  info -> readCO_time += difftime (end, start);

  return (true);
}


