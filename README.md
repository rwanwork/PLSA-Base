Probabilistic latent semantic analysis (PLSA, baseline version)
===============================================================


Introduction
------------

Probabilistic latent semantic analysis is a method for analyzing co-occurrence data.  This source code is a baseline version of PLSA which makes no use of any parallelization.  This work is described by the poster paper:

    R. Wan, V. N. Anh, and H. Mamitsuka. Efficient Probabilistic Latent Semantic Analysis Through Parallelization. In Proc. 5th Asia Information Retrieval Symposium (poster session), volume 5839 of LNCS, pages 432-443, October 2009.

In that paper, there are five systems.  This source code generates one of them:  Baseline.  However, in the description below, we refer to this system as "PLSA".



Requirements
------------

| Software        | Minimum version | Tested version | Required? | Web site                              |
| --------------- | :-------------: | :------------: | :-------: | ------------------------------------- |
| gcc             | 4.3.2           | 10.2.0         | Yes       | http://gcc.gnu.org/                   |
| CMake           | 3.2             | 3.16.3         | Yes       | http://www.cmake.org/                 |


The PLSA software is written in C and has been compiled using the GNU gcc compiler v4.3.2 under Linux.


Compiling
---------

The PLSA-Base software is written in C and has been compiled using v4.3.2 of gcc. The system has been tested on both 32-bit and 64-bit systems, but it does not make use of any features from 64-bit architectures.

CMake is used to compile the software and it is recommended that an "out-of-source" build is performed so as not to clutter the original source directories. We give some brief instructions below on how to do this:

  1. You should have cloned this software from `GitHub`.  Within the main directory, create a directory called `build`.  Then enter it.  (Actually, build/ can be anywhere since it will be deleted later; this is just an example.)
  2. Then run

           cmake ..
           
  where ".." represents the location of the top-level `CMakeLists.txt`.
  3. Type `make` to compile the C source code of PLSA-Base. If this succeeds, then the executable `plsa` will exist in your current directory.


Running PLSA
------------

Run `plsa` without any arguments to get a list of possible options:

    Probabilistic Latent Semantic Analysis (baseline)
    =================================================

    Usage:  ./plsa [options]

    Options:
    --base <file>      :  Base filename for output file.
    --cooccur <file>   :  Co-occurrence filename.
    --clusters <int>   :  Number of clusters.
    --seed <int>       :  Random seed.
                       :    (Default:  current time).
    --maxiter <int>    :  Maximum iterations.
    --text             :  Text mode (I/O is in text, not binary).
    --verbose          :  Verbose mode.
    --debug            :  Debugging output.
    --rounding         :  Round using 100000000 as the multiplication factor.
    --nooutput         :  Suppress outputting p(x,y) to file.

    PLSA version:  Mar  7 2010 (15:10:57)


From top-to-bottom, the options mean:

* --base:      The filename, before the extension, of the output file.  The extension is fixed as ".plsa".
* --cooccur:   The input co-occurrence file, whose format is described below.
* --clusters:  The number of latent states.
* --seed:      The random seed to use.  If none is provided, the current system time is used.
* --maxiter:   The maximum number of iterations of the EM algorithm to perform.  One of two stopping criteria.
* --text:      Indicate that the input file is in text and not binary; useful for debugging.
* --verbose:   Verbose output.
* --debug:     Debugging output.  Output is generated as each value is read from the input file.  (Note that a lot of output will be generated.)
* --rounding:  Round the output values in p(x,y) using the specified rounding factor.  That is, if the factor is "1000", then three decimal places are used.  Useful for comparing methods due to the problem with floating point arithmetic (details below).
* --nooutput:  Do not produce the final output file.  Eliminates the creation of a fairly large file.

Many of these parameters have no defaults (such as `--maxiter` and  `--clusters`), so they will have to be explicitly given.


Data format
-----------

The input file is provided to the system using the --cooccur switch.  It can be in one of two formats:  text or binary.  Binary mode is the default and suggested for large files to save space; text mode is good for testing.

In either mode, the format of the file is:

    [rows][columns][row id+][column id+][w1 cos_count (w21 c21) ... (w2n c2n)]+

which means:

- rows --  The number of rows in the matrix (unsigned int)
- columns --  The number of columns in the matrix (unsigned int)
- row id --  A sequence of row id's.  The number of them should be exactly "rows".
- column id --  A sequence of column id's.  The number of them should be exactly "columns".
- [w1 cos_count (w21 c21) ... (w2n c2n)] -- A row in the matrix which
    describes only the non-zero positions.
    - w1 --  The ID of the row.
    - cos_count --  The number of non-zero values in the row.
    - w21 -- Column ID of the value.
    - c21 -- The value itself.

All values are unsigned integers.  The program does not make use of either "row id" or "column id" (it was included for a future feature of the program that has not yet been implemented).  Any value is fine (i.e., 0 or sequential integers).

In text format, the values are separated by a single whitespace (usually the tab character).  In binary mode,  they are unsigned integers (usually 4 bytes in size each).

Remember to specify whether it is a text or binary file by using (or not using) the `--text` switch.  The program does NOT check for binary or text mode and problems will occur if the mode does not match the file.

Please see the source in input.c for further details on the file format.


Sample run
----------

Suppose we have a 3 by 4 matrix like this:

    0     0     3     4
    0     0     0     6
    7     2     1     0

then the test file will appear like this:

    3     4
    0     1     2
    0     1     2     3
    0     2     2     3     3     4
    1     1     3     6
    2     3     0     7     1     2     2     1


The first 3 rows are self-explanatory.  The fourth row is interpretted as follows.  In row 0, there are 2 non-zero co-occurrences.  They are at position 2 (with a value of 3) and position 3 (with a value of 4).  Therefore, we have:  (0, 2) = 3 and (0, 3) = 4.

We execute `plsa` as follows:

    ./plsa --cooccur test.cooccur --maxiter 30 --clusters 4 --text --verbose --debug --base out

Without `--verbose` and `--debug`, there would be no output generated.  Output will vary from system to system, but some sample output would look like the following:


    Settings
    --------
    ==  Base filename:                                  out
    ==  Co-occurrence filename:                         test.cooccur
    ==  Probability data type:                          double
    ==  Clusters:                                       4
    ==  Random seed:                                    [from time]
    ==  Exponent difference [utils.h::addLogsFloat]:    23.02585093
    ==  Termination conditions
    ==    Maximum EM iterations:                        30
    ==    Percentage difference:                        0.001000
    ==  Text mode:                                      yes
    ==  Rounding:                                       no
    ==  Suppress output to file:                        no


    ==  Reading from co-occurrence file...
    ==  Applying seed from time:                        1267944876
    ==    Read (0, 2) --> 3
    ==    Read (0, 3) --> 4
    ==    Read (1, 3) --> 6
    ==    Read (2, 0) --> 7
    ==    Read (2, 1) --> 2
    ==    Read (2, 2) --> 1
    ==  Maximum number of pairs:                        12
    ==  Actual number of pairs in data file:            6
    ==  Percentage of zeroes:                           50.00 % (6)
    ==  Sum of co-occurrence counts:                    23
    ==  Begin initialization...
    ==  Initialization complete...
    ==  m = 3; n = 4
    [---]  Initial = -61.076698
    [  1]  -61.076698 --> -52.328544  [8.748153, 14.3232 %]
    [  2]  -52.328544 --> -50.542405  [1.786140, 3.4133 %]
    [  3]  -50.542405 --> -47.318958  [3.223446, 6.3777 %]
    [  4]  -47.318958 --> -42.657292  [4.661667, 9.8516 %]
    [  5]  -42.657292 --> -38.997613  [3.659679, 8.5793 %]
    [  6]  -38.997613 --> -37.772321  [1.225292, 3.1420 %]
    [  7]  -37.772321 --> -37.573313  [0.199008, 0.5269 %]
    [  8]  -37.573313 --> -37.535241  [0.038072, 0.1013 %]
    [  9]  -37.535241 --> -37.523889  [0.011352, 0.0302 %]
    [ 10]  -37.523889 --> -37.519834  [0.004055, 0.0108 %]
    [ 11]  -37.519834 --> -37.518251  [0.001583, 0.0042 %]
    [ 12]  -37.518251 --> -37.517603  [0.000648, 0.0017 %]
    [ 13]  -37.517603 --> -37.517333  [0.000270, 0.0007 %]
    ==  Non-probabilities:                              0
    ==  Sum of p(x,y):                                  1.000000
    ==  Total output files printed                      1
    ==  Program execution:                              0.000 secs


The resulting matrix of p(x,y) is stored as out.plsa and look like the following:

    -----
    3  4  0  1  2  0  1  2  3  
    -20.596084  -21.298047  -2.036939  -1.749157  
    -102.984165  -103.935308  -11.724635  -1.343766  
    -1.189584  -2.442347  -3.135510  -14.177547
    -----

Whose output format is the same as the input format, except that the integral co-occurrence counts are replaced with probabilities in log-space as floating point values.


Other issues
------------

This section lists some issues relevant to running or extending the program, in no particular order.

Applicable to both versions:

1.  All values are stored in log-space to prevent underflow of floating point values.

2.  Probabilities are multiplied by adding values in log-space.  To prevent underflow, a macro called `logSumsInline` is used (defined in `plsa-defn.h`).  It is provided a cut-off called `LN_LIMIT` which stipulates what cut-off (as a ln value) to use when adding two numbers.  If the smallest number is "too small", then no addition is performed and the larger value is taken.  As a result, the order in which numbers are added matters and if   MPI is used a small, negligible difference will appear.  In other words, if

    x = a + b + c + d

then when MPI with two processors is used, then

    x = y + w
    y = a + b
    w = c + d

If c is very small compared to (a + b) but not so when compared to d, then in the first case, it might be dropped.  Generally, the percentage difference in maximum likelihood will not change enough to matter.

3.  Also related to the previous point, p(x,y) values written to the 
output file may be different if the binary files are compared directly due 
to the number of digits of precision.  If sent to a text file with only 6 
digits of precision, there may not be any problems.  To solve this, the 
`--rounding` option was added.  It and the cut-off defined in plsa-defn.h 
will need to be tweaked if two runs are to be compared.

4.  Generating the output file could take half of the total execution time, depending on the matrix size.  To suppress printing this file, use the `--nooutput` switch.

5.  To get the time required for a single iteration of the loop (as reported in the paper cited in Section 1), use the `--maxiter 1` option.

6.  The four random seeds used for the paper cited in Section 1 were:  "20444 3612 31325 17062".  Of course, this information alone is not enough to generate the exact same results since every system's random number generator might be slightly different.

7.  The number of latent states must be larger than the number of processors under MPI.  If this is not the case, then the number of latent states is increased automatically.


Applicable to this version only:

1.  There is no maximum number of latent states.


About PLSA
----------

This software was implemented while I was at Kyoto University (around 2009).  My contact details:

     E-mail:  rwan.work@gmail.com 

My homepage is [here](http://www.rwanwork.info/).

The latest version of PLSA can be downloaded from [GitHub](https://github.com/rwanwork/PLSA-Base).

If you have any information about bugs, suggestions for the documentation or just have some general comments, feel free to contact me via e-mail or GitHub.


Copyright and license
---------------------

    Probablistic latent semantic analysis (PLSA, baseline version)
    Copyright (C) 2009-2021 by Raymond Wan

This software is distributed under the terms of the GNU General Public License (GPL, version 3 or later) -- see the file LICENSE for details.

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts and no Back-Cover Texts. A copy of the license is included with the archive as LICENSE.

