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

#include "MyUtils.h"
#include <math.h>
#include <stddef.h>
#include <iostream>


using namespace std;

typedef struct
{
  int xl, xr, y, dy;
}
LINESEGMENT;

#define MAXDEPTH 10000

#define PUSH(XL, XR, Y, DY) \
  if( sp < stack+MAXDEPTH && Y+(DY) >= min_y && Y+(DY) <= max_y ) \
{ sp->xl = XL; sp->xr = XR; sp->y = Y; sp->dy = DY; sp++; }

#define POP(XL, XR, Y, DY) \
{ sp--; XL = sp->xl; XR = sp->xr; Y = sp->y + (DY = sp->dy); }

void MyUtils::FloodFill(char** data, int x, int y,
                        int min_x, int min_y,
                        int max_x, int max_y,
                        int fill_value, int border_value)
{
  int left, x1, x2, dy;
  LINESEGMENT stack[MAXDEPTH], *sp = stack;

  if (data[y][x] == border_value || data[y][x] == fill_value)
  {
    return;
  }

  if (x < min_x || x > max_x || y < min_y || y > max_y)
  {
    return;
  }

  PUSH(x, x, y, 1);        /* needed in some cases */
  PUSH(x, x, y+1, -1);    /* seed segment (popped 1st) */

  while (sp > stack )
  {
    POP(x1, x2, y, dy);

    for (x = x1; x >= min_x &&
         data[y][x] != border_value && data[y][x] != fill_value; x--)
    {
      data[y][x] = fill_value;
    }

    if ( x >= x1 )
    {
      goto SKIP;
    }

    left = x+1;
    if ( left < x1 )
    {
      PUSH(left, x1-1, y, -dy);  /* leak on left? */
    }

    x = x1+1;

    do
    {
      for (; x<=max_x &&
           data[y][x] != border_value && data[y][x] != fill_value; x++)
      {
        data[y][x] = fill_value;
      }

      PUSH(left, x-1, y, dy);

      if (x > x2+1)
      {
        PUSH(x2+1, x-1, y, -dy);  /* leak on right? */
      }

SKIP:
      for (x++; x <= x2 &&
           (data[y][x] == border_value || data[y][x] == fill_value); x++)
      {
        ;
      }

      left = x;
    }
    while (x <= x2);
  }
}

void MyUtils::FloodFill(char* data, int x, int y,
                        int dim_x, int dim_y,
                        int fill_value, int border_value)
{
  int left, x1, x2, dy;
  LINESEGMENT stack[MAXDEPTH], *sp = stack;

  if (data[y*dim_x+x] == border_value || data[y*dim_x+x] == fill_value)
  {
    return;
  }

  int min_x = 0, max_x = dim_x-1, min_y = 0, max_y = dim_y-1;
  if (x < min_x || x > max_x || y < min_y || y > max_y)
  {
    return;
  }

  PUSH(x, x, y, 1);        /* needed in some cases */
  PUSH(x, x, y+1, -1);    /* seed segment (popped 1st) */

  while (sp > stack )
  {
    POP(x1, x2, y, dy);

    for (x = x1; x >= min_x &&
         data[y*dim_x+x] != border_value && data[y*dim_x+x] != fill_value; x--)
    {
      data[y*dim_x+x] = fill_value;
    }

    if ( x >= x1 )
    {
      goto SKIP;
    }

    left = x+1;
    if ( left < x1 )
    {
      PUSH(left, x1-1, y, -dy);  /* leak on left? */
    }

    x = x1+1;

    do
    {
      for (; x<=max_x &&
           data[y*dim_x+x] != border_value && data[y*dim_x+x] != fill_value; x++)
      {
        data[y*dim_x+x] = fill_value;
      }

      PUSH(left, x-1, y, dy);

      if (x > x2+1)
      {
        PUSH(x2+1, x-1, y, -dy);  /* leak on right? */
      }

SKIP:
      for (x++; x <= x2 &&
           (data[y*dim_x+x] == border_value || data[y*dim_x+x] == fill_value); x++)
      {
        ;
      }

      left = x;
    }
    while (x <= x2);
  }
}

