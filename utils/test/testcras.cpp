/**
 * @brief testing sampled volume cras calculation
 *
 */
/*
 * Original Author: Y. Tosa
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

#include <iostream>
#include <iomanip>
#include <cstdlib>

extern "C"
{

#include "error.h"
#include "mri.h"

const char *Progname = "testcras";
}

using namespace std;

int PrettyMatrixPrint(MATRIX *mat)
{
  int row;

  if (mat == NULL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat = NULL!")) ;

  if (mat->type != MATRIX_REAL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not Real type")) ;

  if (mat->rows != 4 || mat->cols != 4)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not of 4 x 4")) ;

  for (row=1; row < 5; ++row)
    printf("              %8.4f %8.4f %8.4f %10.4f\n",
           mat->rptr[row][1], mat->rptr[row][2], mat->rptr[row][3], mat->rptr[row][4]);
  return (NO_ERROR);
}

void printInfo(MRI *mri)
{
  if (mri->nframes > 1)
    printf("    dimensions: %d x %d x %d x %d\n", mri->width, mri->height, mri->depth, mri->nframes) ;
  else
    printf("    dimensions: %d x %d x %d\n", mri->width, mri->height, mri->depth) ;
  printf("   voxel sizes: %6.4f, %6.4f, %6.4f\n", mri->xsize, mri->ysize, mri->zsize) ;
  printf("          type: %s (%d)\n",
         mri->type == MRI_UCHAR   ? "UCHAR" :
         mri->type == MRI_SHORT   ? "SHORT" :
         mri->type == MRI_INT     ? "INT" :
         mri->type == MRI_LONG    ? "LONG" :
         mri->type == MRI_BITMAP  ? "BITMAP" :
         mri->type == MRI_TENSOR  ? "TENSOR" :
         mri->type == MRI_FLOAT   ? "FLOAT" : "UNKNOWN", mri->type) ;
  printf("           fov: %2.3f\n", mri->fov) ;
  printf("        xstart: %2.1f, xend: %2.1f\n", mri->xstart*mri->xsize, mri->xend*mri->xsize) ;
  printf("        ystart: %2.1f, yend: %2.1f\n", mri->ystart*mri->ysize, mri->yend*mri->ysize) ;
  printf("        zstart: %2.1f, zend: %2.1f\n", mri->zstart*mri->zsize, mri->zend*mri->zsize) ;
  printf("            TR: %2.2f msec, TE: %2.2f msec, TI: %2.2f msec, flip angle: %2.2f degrees\n",
         mri->tr, mri->te, mri->ti, DEGREES(mri->flip_angle)) ;
  printf("       nframes: %d\n", mri->nframes) ;
  printf("ras xform %spresent\n", mri->ras_good_flag ? "" : "not ") ;
  printf("    xform info: x_r = %8.4f, y_r = %8.4f, z_r = %8.4f, c_r = %10.4f\n",
         mri->x_r, mri->y_r, mri->z_r, mri->c_r);
  printf("              : x_a = %8.4f, y_a = %8.4f, z_a = %8.4f, c_a = %10.4f\n",
         mri->x_a, mri->y_a, mri->z_a, mri->c_a);
  printf("              : x_s = %8.4f, y_s = %8.4f, z_s = %8.4f, c_s = %10.4f\n",
         mri->x_s, mri->y_s, mri->z_s, mri->c_s);

  printf("\nvoxel to ras transform:\n");
  PrettyMatrixPrint(mri->i_to_r__);
  printf("\nras to voxel transform:\n");
  PrettyMatrixPrint(mri->r_to_i__);
}

int main(int argc, char *argv[])
{
  if (argc <= 1)
  {
    cout << "Usage: testcras srcvolname" << endl;
    return -1;
  }
  MRI *src = MRIread(argv[1]);
  if (!src)
  {
    cerr << "could not load the volume" << endl;
    return -1;
  }

  MRI *dst = MRIreduce(src, 0);

  cout << "Src volume info: " << endl;
  printInfo(src);
  cout << "\nDst volume info: " << endl;
  printInfo(dst);

  if (argc==3)
  {
    cout << "writing dst volume as " << argv[2] << endl;
    MRIwrite(dst, argv[2]);
  }
  MRIfree(&src);
  MRIfree(&dst);
}
