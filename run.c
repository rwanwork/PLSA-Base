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
#include <math.h>  /*  fabs  */
#include <float.h>  /*  DBL_EPSILON  */
#include <time.h>
#include <signal.h>

#include "wmalloc.h"
#include "plsa-defn.h"
#include "em-estep.h"
#include "em-mstep.h"
#include "input.h"
#include "output.h"
#include "parameters.h"
#include "debug.h"
#include "run.h"


INFO *initialize () {
  INFO *info = wmalloc (sizeof (INFO));

  time (&(info -> program_start));
  info -> run_time = 0;
  info -> readCO_time = 0;
  info -> initEM_time = 0;
  info -> calculateML_time = 0;
  info -> applyEStep_time = 0;
  info -> applyMStep_time = 0;
  info -> normalizeProbs_time = 0;
  info -> printCoProbs_time = 0;

  /*  MPI unavailable in this version  */
  info -> world_id = MAINPROC;
  info -> world_size = 1;

  /*  Set a handler for floating point exceptions  */
  info -> sigfpe_count = 0;
  signal (SIGFPE, handler_sigfpe);

  return info;
}


void uninitialize (INFO *info) {
  double total_time = 0;
  unsigned int k = 0;

  wfree (info -> cos);
  wfree (info -> probw1_z);
  wfree (info -> probw2_z);
  wfree (info -> probz);

  for (k = 0; k < info -> num_clusters; k++) {
    wfree (info -> probz_w1w2[k]);
  }
  wfree (info -> probz_w1w2);
  wfree (info -> base_fn);
  wfree (info -> co_fn);
  wfree (info -> row_ids);
  wfree (info -> column_ids);

  time (&(info -> program_end));

  if (info -> verbose) {
    total_time = difftime (info -> program_end, info -> program_start);
    if (total_time > 60) {
      fprintf (stderr, "==\tProgram execution:                              %.3f mins\n", total_time / 60);
    }
    else {
      fprintf (stderr, "==\tProgram execution:                              %.3f secs\n", total_time);
    }
    if (total_time > 1) {
      fprintf (stderr, "==\t  run() time:                                   %6.2f %%\n", info -> run_time / total_time * 100);

      fprintf (stderr, "==\t    Read data in:                               %6.2f %%\n", info -> readCO_time / total_time * 100);
      fprintf (stderr, "==\t    EM initialization:                          %6.2f %%\n", info -> initEM_time / total_time * 100);
      fprintf (stderr, "==\t    Calculate ML:                               %6.2f %%\n", info -> calculateML_time / total_time * 100);
      fprintf (stderr, "==\t    Apply E step:                               %6.2f %%\n", info -> applyEStep_time / total_time * 100);
      fprintf (stderr, "==\t    Apply M step:                               %6.2f %%\n", info -> applyMStep_time / total_time * 100);
      fprintf (stderr, "==\t    Normalize probabilities:                    %6.2f %%\n", info -> normalizeProbs_time / total_time * 100);
      fprintf (stderr, "==\t    Print probabilities:                        %6.2f %%\n", info -> printCoProbs_time / total_time * 100);
    }
  }

  wfree (info);

  return;
}


bool run (INFO *info) {
  PROBNODE curr_ML = 0;
  PROBNODE prev_ML = 0;
  PROBNODE diff = 0.0;

  info -> iter = 0;

  time_t start;
  time_t end;
  time_t loop_start;
  time_t loop_end;
  double timediff = 0.0;

  time (&start);

  /*  All processes read in co-occurrence data  */
  if (!readCO (info)) {
    /*  If there is an error, all processes are terminated  */
    fprintf (stderr, "Error reading co-occurrence data.\n");
    return false;
  }

  /*  Only the main process initializes to ensure the random seed affects it only  */
  initEM (info);
  if (info -> verbose) {
    fprintf (stderr, "==\tm = %u; n = %u\n", info -> m, info -> n);
  }


  time (&loop_start);
  while (true) {
    curr_ML = calculateML (info);

    if (info -> iter == 0) {
      if (info -> verbose) {
        fprintf (stderr, "[---]  Initial = %f\n", curr_ML);
      }
    }
    else {
      diff = (curr_ML - prev_ML) / prev_ML * 100 * -1;
      if (info -> verbose) {
        fprintf (stderr, "[%3u]  %f --> %f\t[%f, %2.4f %%]\n", info -> iter, prev_ML, curr_ML, (curr_ML - prev_ML), diff);
      }
      if ((curr_ML < prev_ML) || (DBL_LESS (fabs (diff), ML_DELTA))) {
        info -> iter = UINT_MAX;  /*  Set an indicator to leave loop  */
      }
    }

    prev_ML = curr_ML;

#if DEBUG
    checkCoProb (info);
#endif

    if (info ->  iter != UINT_MAX) {
      info -> iter++;
    }

    if (info -> iter > (info -> maxiter)) {
      info -> iter = UINT_MAX;  /*  Set an indicator to leave loop  */
    }

    /*  Check if we are suppose to exit this loop  */
    if (info -> iter == UINT_MAX) {
      break;
    }

    applyEStep (info);

    /*  Calculate M-step  */
    applyMStep (info);

    normalizeProbs (info);
  }
  time (&loop_end);
  timediff += difftime (loop_end, loop_start);

  if (info -> maxiter == 1) {
    fprintf (stderr, "==\t  Main loop [one iteration only!]:             %6.2f %% (%f)\n", 0.0, timediff);
  }

  if (!info -> no_output) {
    printCoProb (info);
  }

  time (&end);
  info -> run_time += difftime (end, start);

  return (true);
}

