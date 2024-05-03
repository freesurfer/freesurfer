/**
 * @brief Misc utility class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#ifndef MyUtils_h
#define MyUtils_h

#include <math.h>
#include <vector>

class MyUtils
{
public:
  MyUtils()
  {}

  static void FloodFill( char** data, int x, int y, int min_x, int min_y, int max_x, int max_y, int fill_value, int border_value );
  static void FloodFill( char* data, int x, int y, int dim_x, int dim_y, int fill_value, int border_value );
};

#endif
