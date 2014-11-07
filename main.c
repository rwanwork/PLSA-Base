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
#include <time.h>

#include "plsa-defn.h"
#include "wmalloc.h"
#include "parameters.h"
#include "run.h"


/*!  Main function  */
int main (int argc, char *argv[]) {
  INFO *info;
  bool result = false;

  info = initialize ();

  /*  Process the command line parameters and then check them;
  **  if either fail then print usage information  */
  if ((!processOptions (argc, argv, info)) || (!checkSettings (info))) {
    usage (argv[0]);
  }
  else {
    result = run (info);
  }

  uninitialize (info);

  return (EXIT_SUCCESS);
}

