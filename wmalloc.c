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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wmalloc.h"

static unsigned int inuse_malloc = 0;
static unsigned int max_malloc = 0;
static WMSTRUCT **wm_array;
static char *tempstr;

void *wmalloc (size_t y_arg) {
  void *x_arg = malloc (y_arg);
  if (x_arg == NULL) {
    fprintf (stderr, "Error in malloc while allocating %u bytes in [%s, %u].\n", (unsigned int) y_arg, __FILE__, __LINE__);
    exit (EXIT_FAILURE);
  }
#ifdef COUNT_MALLOC
  countMalloc (x_arg, y_arg, __FILE__, __LINE__);
#endif

  return (x_arg);
}


void *wrealloc (void *x_arg, size_t y_arg) {
#ifdef COUNT_MALLOC
  countFree ((void*) x_arg);
#endif
  x_arg = realloc (x_arg, y_arg);

  if (x_arg == NULL) {
    fprintf (stderr, "Error in realloc while allocating %u bytes in [%s, %u].\n", (unsigned int) y_arg, __FILE__, __LINE__);
    exit (EXIT_FAILURE);
  }
#ifdef COUNT_MALLOC
  countMalloc ((void*) x_arg, y_arg, __FILE__, __LINE__);
#endif

  return (x_arg);
}

void wfree (void *x_arg) {
#ifdef COUNT_MALLOC
  countFree (x_arg);
#endif
  free (x_arg);
}

/*  Function adapted from Algorithms in C (Third edition) by Robert Sedgewick
**  (page 578) */
static unsigned int hash (char *v, signed int M) {
  signed int h = 0;
  signed int a = 127;

  /*  Skip the 0x at the beginning  */
  v += 2;
  for (; *v != '\0'; v++) {
    h = (a*h + (signed int) *v) % M;
  }

  return ((unsigned int) h);
}

void initWMalloc () {
  WMSTRUCT *node = NULL;
  unsigned int i = 0;

  inuse_malloc = 0;
  max_malloc = 0;

  tempstr = malloc (sizeof (char) * TEMPSTRLEN);

  wm_array = malloc (sizeof (WMSTRUCT*) * WM_SIZE);
  for (i = 0; i < WM_SIZE; i++) {
    /*  Add sentinel  */
    node = malloc (sizeof (WMSTRUCT));
    node -> ptr = NULL;
    node -> size = 0;
    node -> file = NULL;
    node -> line = 0;
    node -> next = NULL;
    wm_array[i] = node;
  }

  return;
}


void printWMalloc () {
  WMSTRUCT *prev = NULL;
  WMSTRUCT *curr = NULL;
  unsigned int i = 0;

  fprintf (stderr, "\tMemory used at exit:  %u\n", inuse_malloc);
  fprintf (stderr, "\tMaximum memory used at once:  %u\n", max_malloc);
  fprintf (stderr, "\tMaximum memory (MB):  %.1f\n", (double) max_malloc / (double) (1024 * 1024));

  for (i = 0; i < WM_SIZE; i++) {
    prev = NULL;
    curr = wm_array[i];

    while (curr != NULL) {
      prev = curr;
      curr = curr -> next;
      free (prev -> file);
      free (prev);
    }
  }

  free (tempstr);
  free (wm_array);

  return;
}


void printInUseWMalloc (void) {
  unsigned int i = 0;
  WMSTRUCT *curr = NULL;

  for (i = 0; i < WM_SIZE; i++) {
    if (wm_array[i] -> ptr != NULL) {
      curr = wm_array[i];
      while (curr -> ptr != NULL) {
        fprintf (stderr, "%p\t(%u)\t[%s]\t[%u]\n", curr -> ptr, (unsigned int) curr -> size, curr -> file, curr -> line);
        curr = curr -> next;
      }
    }
  }

  return;
}


void countMalloc (void *ptr, size_t amount, const char *file, const unsigned int line) {
  WMSTRUCT *node = NULL;
  unsigned int pos = 0;

  node = malloc (sizeof (WMSTRUCT));
  node -> ptr = ptr;
  node -> size = amount;
  node -> file = malloc (sizeof (char) * (strlen (file) + 1));
  node -> file = strcpy (node -> file, file);
  node -> line = line;

  (void) snprintf (tempstr, TEMPSTRLEN, "%p", ptr);
  pos = hash (tempstr, WM_SIZE);
  node -> next = wm_array[pos];
  wm_array[pos] = node;

  inuse_malloc += amount;

  if (inuse_malloc > max_malloc) {
    max_malloc = inuse_malloc;
  }

  return;
}


void countFree (void *ptr) {
  WMSTRUCT *prev = NULL;
  WMSTRUCT *curr = NULL;
  unsigned int pos = 0;

  (void) snprintf (tempstr, TEMPSTRLEN, "%p", ptr);
  pos = hash (tempstr, WM_SIZE);

  /*  Locate the node  */
  curr = wm_array[pos];
  while (curr -> ptr != NULL) {
    if (curr -> ptr == ptr) {
      break;
    }
    prev = curr;
    curr = curr -> next;
  }
  if (curr -> ptr == NULL) {
    fprintf (stderr, "(F) Fatal error.  Node %p could not be found.\n", ptr);
    exit (EXIT_FAILURE);
  }

  inuse_malloc -= curr -> size;

  if (inuse_malloc > max_malloc) {
    max_malloc = inuse_malloc;
  }

  /*  Remove the node  */
  if (prev == NULL) {
    wm_array[pos] = curr -> next;
  }
  else {
    prev -> next = curr -> next;
  }
  free (curr -> file);
  free (curr);

  return;
}


