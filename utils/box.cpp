/*
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
 */

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "box.h"
#include "error.h"
#include "macros.h"
#include "proto.h"

int BoxPrint(BOX *box, FILE *fp)
{
  fprintf(fp, "BOX: x: %d --> %d, y: %d --> %d, z: %d --> %d\n", box->x0, box->x1, box->y0, box->y1, box->z0, box->z1);
  return (NO_ERROR);
}
int BoxExpand(BOX *box_src, BOX *box_dst, int dx, int dy, int dz)
{
  box_dst->x0 = MAX(box_src->x0 - dx, 0);
  box_dst->y0 = MAX(box_src->y0 - dy, 0);
  box_dst->z0 = MAX(box_src->z0 - dz, 0);
  box_dst->x1 = MIN(box_src->x1 + dx, box_src->width - 1);
  box_dst->y1 = MIN(box_src->y1 + dy, box_src->height - 1);
  box_dst->z1 = MIN(box_src->z1 + dz, box_src->depth - 1);
  box_dst->width = box_src->width;
  box_dst->height = box_src->height;
  box_dst->depth = box_src->depth;
  return (NO_ERROR);
}
