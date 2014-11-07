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

#ifndef WMALLOC_H
#define WMALLOC_H

#define WM_SIZE 65536
#define TEMPSTRLEN 80

typedef struct wmstruct {
  void *ptr;
  size_t size;
  char *file;
  unsigned int line;
  struct wmstruct *next;
} WMSTRUCT;

void *wmalloc (size_t y_arg);
void *wrealloc (void *x_arg, size_t y_arg);
void wfree (void *x_arg);

void initWMalloc (void);
void printWMalloc (void);
void printInUseWMalloc (void);
void countMalloc (void *ptr, size_t amount, const char *file, unsigned int line);
void countFree (void *ptr);

#endif

