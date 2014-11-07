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

#ifndef PLSA_DEFN_H
#define PLSA_DEFN_H

/*
**  Define the type of floating point to use; highly recommend to use only double
*/
#if 1
/*!  Data type to use for probabilities  */
typedef double PROBNODE;
#else
/*!  Data type to use for probabilities  */
typedef float PROBNODE;
#endif


/*
**  e^(-87.49823353) = 1.0E-38
**  e^(-73.682723)   = 1.0E-32
**  e^(-55.26204223) = 1.0E-24
**  e^(-25.32843602) = 0.00000000001
**  e^(-23.02585093) = 0.0000000001
**  e^(-20.72326584) = 0.000000001
**  e^(-18.42068074) = 0.00000001
**  e^(-16.11809565) = 0.0000001
**  e^(-13.81551056) = 0.000001
**  e^(-11.51292547) = 0.00001
**  e^(-9.210340372) = 0.0001
*/
/*!  Accuracy of floating point values as a log (base e) value, multiplied by -1  */
#define LN_LIMIT 23.02585093

/*!  Minimum probability  */
#define MIN_PROB (1.0E-24)

/*!  Macro to perform a log  */
#define DOLOG(X) (logf (X))

/*!  Macro to perform the exp function  */
#define DOEXP(X) (expf (X))

/*!  Macro to perform log (1 + x)  */
#define DOLOGONE(X) (log1pf (X))

/*!  Macro to perform log (1 + expt(x))  */
#define DOLOG1PEXP(x) DOLOGONE(DOEXP(x))

/*!  Generate a random number between [0, 1); cast to floating point first to prevent overflow  */
#define RANDOM_FLOAT ((PROBNODE)rand () / ((PROBNODE)RAND_MAX + (PROBNODE)1.0))

/*!  Test if two double values are close to each other  */
#define DBL_LESS(A,B) ((B - A) > DBL_EPSILON)

/*!  Minimum difference between two maximum likelihoods  */
#define ML_DELTA 0.001

/*!  ID of the main processor is always 0  */
#define MAINPROC 0

/*!  Number of digits to round; used when outputting to binary only  */
#define ROUND_DIGITS 100000000

/********************************************************************/
/*  Definitions and functions related to MPI  */

/*
**  Definitions from Quinn (2003),  pg.120 [BLOCK_SIZE corrected]  */
/*  id = process rank;
**  p = total number of process;
**  n = number of items
**  index = position of the item to see who is responsible for it
*/
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW ((id) + 1, p, n)-BLOCK_LOW(id, p, n))

/********************************************************************/
/*   Inline functions  */

/*!  Define'd function to indicate program progress  */
#define PROGRESS_MSG(A) \
if (info -> verbose) { \
  fprintf (stderr, "==\t%s\n", A); \
}

#define FOPEN(FILENAME,FP,MODE) \
  FP = fopen ((char*) FILENAME, MODE); \
  if (FP == NULL) { \
    fprintf (stderr, "Error %s %s.\n", (strcmp (MODE, "w") == 0) ? "creating" : "opening", FILENAME); \
    exit (EXIT_FAILURE); \
  }

#define FCLOSE(FP) \
  (void) fclose (FP);

/********************************************************************/
/*  Functions for accessing cooccurrence structure  */

/*!  Function to retrieve from the cooccurrence array  */
#define SET_COS(W,X,Y,Z) \
{ \
  info -> cos[W][X].column = Y; \
  info -> cos[W][X].x = Z; \
}

/*!  Function to retrieve the position from the cooccurrence array  */
#define GET_COS(W,X) (info -> cos[W][X].x)

/*!  Function to retrieve the cooccurrence count from the cooccurrence array  */
#define GET_COS_POSITION(W,X) (info -> cos[W][X].column)

/********************************************************************/
/*  Functions for accessing probabilities  */
/*!  Function to retrieve from P(w1|z); translate 2D to 1D co-ordinates  */
#define GET_PROBW1_Z(X,Y) (info -> probw1_z[X * info -> m + Y])

/*!  Function to retrieve from P(w2|z); translate 2D to 1D co-ordinates  */
#define GET_PROBW2_Z(X,Y) (info -> probw2_z[X * info -> n + Y])

/*!  Function to retrieve from P(z)  */
#define GET_PROBZ(X) (info -> probz[X])

/*!  Function to retrieve from P(z|w1w2); translate 3D to 1D co-ordinates  */
#define GET_PROBZ_W1W2(W,X,Y) (info -> probz_w1w2[W][X * info -> n + Y])

#define logSumsInline(A,B) \
{                          \
  register PROBNODE x, y;  \
  if (A > B) {             \
    x = A;  y = B;         \
  }                        \
  else {                   \
    x = B;  y = A;         \
  }                        \
                           \
  /*  a > b  */            \
                           \
  A = (fabs (y - x) > LN_LIMIT) ? x : x + DOLOG1PEXP (y - x);   \
}

/********************************************************************/
typedef struct cooccur {
  /*!  The co-occurrence count, as a log value  */
  PROBNODE x;
  /*!  Column position of this value  */
  unsigned int column;
} COOCCUR;


typedef struct info {
  /*!  Verbose output?  */
  bool verbose;
  /*!  Debugging output?  */
  bool debug;
  /*!  Text I/O  */
  bool textio;
  /*!  Should the output values be rounded?  */
  bool rounding;
  /*!  Suppress output  */
  bool no_output;

  /*!  Random seed  */
  unsigned int seed;
  /*!  Number of clusters  */
  unsigned int num_clusters;
  /*!  Base filename for the output file  */
  char *base_fn;
  /*!  Maximum number of iterations  */
  unsigned int maxiter;
  /*!  Number of unique query terms  */
  unsigned int m;
  /*!  Number of terms in the document collection  */
  unsigned int n;

  /*!  Co-occurrence filename  */
  char *co_fn;
  /*!  Co-occurrence counts in a COOCCUR data structure  */
  COOCCUR **cos;
  /*!  List of row identifiers (m of them)  */
  unsigned int *row_ids;
  /*!  List of column identifiers (m of them)  */
  unsigned int *column_ids;

  /*!  Iteration; only calculated by the main process and broadcasted to others  */
  unsigned int iter;

  /*!  P(w1|z) of size (k * m)  */
  PROBNODE *probw1_z;
  /*!  P(w2|z) of size (k * n)  */
  PROBNODE *probw2_z;
  /*!  P(z) of size (k); one-dimensional array does not need a pointer  */
  PROBNODE *probz;
  /*!  P(z|w1w2) of size (k * m * n)  */
  PROBNODE **probz_w1w2;

  /*  Variables specific to MPI  */
  /*!  ID of this process  */
  signed int world_id;
  /*!  Number of processes total  */
  signed int world_size;
  /*!  Starting block (cluster) for this process to handle  */
  unsigned int block_start;
  /*!  Ending block (cluster) for this process to handle  */
  unsigned int block_end;
  /*!  Size of the block for this process to handle  */
  unsigned int block_size;

  /*!  Number of floating point exception errors  */
  unsigned int sigfpe_count;

  /*  Various times  */
  time_t program_start;
  double run_time;
  double readCO_time;
  double initEM_time;
  double calculateML_time;
  double applyEStep_time;
  double applyMStep_time;
  double normalizeProbs_time;
  double printCoProbs_time;
  time_t program_end;
} INFO;

#endif
