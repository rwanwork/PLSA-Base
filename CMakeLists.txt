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


##  Set the minimum required CMake version
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)


########################################
##  Define some variables

##  Testing compilation flags, some of which are suggested by the Valgrind 3.3 >
# set (MY_C_FLAGS "-pedantic -Wno-long-long -g -fno-inline -O0 -Wall")
##  Release compilation flags, suggested by the Valgrind 3.3 book
set (MY_C_FLAGS "-O3 -Wall -Wno-unused-result -Wno-unused-but-set-variable")

set (TARGET_NAME_EXEC "plsa")
set (CURR_PROJECT_NAME "PLSA-Base")

##  Define the project
project (${CURR_PROJECT_NAME} LANGUAGES C)


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


########################################
##  Define the executable

add_executable (${TARGET_NAME_EXEC} ${SRC_FILES})

##  Link the executable to the math library
target_link_libraries (${TARGET_NAME_EXEC} m)


########################################
##  Set initial compilation flags

##  Set compiler flags based on global variable
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")


