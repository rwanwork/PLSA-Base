###########################################################################
##  PLSA-Base
##  Probabilistic latent semantic analysis (PLSA, baseline version)
##  
##  Copyright (C) 2009-2020 by Raymond Wan, All rights reserved.
##  Contact:  rwan.work@gmail.com
##  
##  This file is part of PLSA-Base.
##  
##  PLSA-Base is free software; you can redistribute it and/or 
##  modify it under the terms of the GNU General Public License 
##  as published by the Free Software Foundation; either version 
##  3 of the License, or (at your option) any later version.
##  
##  PLSA-Base is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU Lesser General Public License for more details.
##  
##  You should have received a copy of the GNU General Public 
##  License along with PLSA-Base; if not, see 
##  <http://www.gnu.org/licenses/>.
###########################################################################


project ("PLSA-Base" LANGUAGES C)

cmake_minimum_required (VERSION 3.2 FATAL_ERROR)


########################################
##  Define some variables

##  Testing compilation flags, some of which are suggested by the Valgrind 3.3 >
# set (MY_CXX_FLAGS "-pedantic -Wno-long-long -g -fno-inline -O0 -Wall")
##  Release compilation flags, suggested by the Valgrind 3.3 book
set (MY_CXX_FLAGS "-O3 -Wall")

set (TARGET_NAME_EXEC "plsa")
set (CURR_PROJECT_NAME "PLSA-Base")

########################################
##  Define the source files
##  Source files for both the test executable and library
set (SRC_FILES
  debug.c
  em-estep.c
  em-mstep.c
  input.c
  main.c
  output.c
  parameters.c
  run.c
  wmalloc.c
)

add_executable (${TARGET_NAME_EXEC} ${SRC_FILES})
target_link_libraries (${TARGET_NAME_EXEC} m)

########################################
##  Set initial compilation flags

##  Set compiler flags based on global variable
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")

