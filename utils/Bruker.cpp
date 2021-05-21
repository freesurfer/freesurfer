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

/////////////////////////////////////////////////////////////////////////
/* Bruker.c                               */
/* created by : y.tosa                    */
/* date       :8/27/2003                  */
// Warning: Do not edit the following four lines.  CVS maintains them.

// there are many files present in Bruker directory
//
// fid             ... raw data
// method          ... store similar info as acqp  <--- used
// acqp            ... used for reconstruction     <--- used
// pulseprogram    ...
// spnam()         ...
// gdprog.ax, ay, az .
// pdata ... reconstructed image strorage
//   |
//   |-> 1
//       |-> 2dseq  ... bshort image                                   <-- used
//           reco   ... created after reconstruction                   <-- used
//           d3proc ... reconstruted image info (width, height, depth) <-- used
//           procs  ... ???
//           roi    ... ???
//
/**********************************************/
/* email from Maritin Horrman <mah@bruker.de> */
/* on how to reconstruct voxelToRAS transform */
/**********************************************/
/* Dear Dr. Tosa, */

/* here I deliver you a suggestion, to get the transformation matrix: */

/* I think, this should work, if you have a single 3D slab or */
/* a 2D data set with 1 slice package and 1 image per slice. */

/* If you have data sets with several 3D slabs or 2D image files with */
/* several slice packages or several images per slice, it is not that easy */
/* to make the correct conversion. */

/* Start with the reconstructed image. */
/* Transpose as per RECO_transposition: eg is RECO_transposition is 1 start with
 */

/* 0 1 0 0 */
/* 1 0 0 0 */
/* 0 0 1 0 */
/* 0 0 0 1 */

/* or is RECO_transposition = 0 start with a 4x4 identity matrix. */

/* Then convert to coordinates from the centre of the ft'ed image, in mm; with
 */
/* first voxel or coordinates 0,0,0, this gives a vox -> mm matrix: */

/* vox(dim)  = RECO_fov(dim) * 10 / RECO_size(dim) */
/* off(dim) = -vox(dim) * (RECO_ft_size(dim)-1)/2; */
/* (for 2d data sets, vox(3) = ACQ_slice_sepn, */
/*                          off(3)  = -vox(3) * (NSLICES-1)/2) */

/* vox(1)   0 0 off(1) */
/* 0 vox(2) 0 0 off(2) */
/* 0 0   vox(3) off(3) */
/* 0 0 0 1 */

/* swapping matrix to express reversal of x and y is given by */

/* -1  0  0 0 */
/*  0 -1  0 0 */
/*  0  0  1 0   (may be -1 in case of 3d) */
/*  0  0  0 1 */

/* translation matrix is given by */

/* 1 0 0 ACQ_read_offset; */
/* 0 1 0 ACQ_phase1_offset; */
/* 0 0 1 ACQ_slice_offset; */
/* 0 0 0 1 */
/* ( for 2D image with 1 slice package use ACQ_slice_offset of the middle slice
 */
/*   = IMND_slicepack_position[0] or PVM_SPackArrSliceOffset[0].) */

/* rotation matrix is given by */

/* gm = ACQ_grad_matrix for slice / echo, assuming this is the same for all */
/* slices */

/* gm(1) gm(4)  gm(7) 0 */
/* gm(2) gm(5)  gm(8) 0 */
/* gm(3) gm(6)  gm(9) 0 */
/* 0     0      0     1 */

/* If RAS means the following coordinate system: */
/* R->L, A->P, H->F */
/* you must also apply the ras matrix. */

/* The ras matrix is given by */

/* -1 0 0 0 */
/* 0 -1 0 0 */
/* 0  0 1 0 */
/* 0  0 0 1 */

/* The ParaVision patient coordinate system is: */
/* L->R, */

/* The whole 4x4 matrix to give mm coordinates in terms of isocentre, from voxel
 */
/* coordinates, is */

/* [ras matrix *] rotation matrix * translation matrix * swap matrix * vox->mm
 */
/* matrix * transposition matrix */

/* I hope, this helps you a little bit, with best regards */

/* Martin Hoerrmann. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "matrix.h"
#include "mri.h"
#include "utils.h"

#include "Bruker.h"

#define V4_LOAD(v, x, y, z, r) (VECTOR_ELT(v, 1) = x, VECTOR_ELT(v, 2) = y, VECTOR_ELT(v, 3) = z, VECTOR_ELT(v, 4) = r);

MRI *brukerRead(char *fname, int read_volume)
{
  char methodFile[1024];
  char acqpFile[1024];
  char dataFile[1024];
  char d3procFile[1024];
  char recoFile[1024];
  BrukerTransform bTran;
  double TR, TE, TI, flip_angle;

  int succeed = 0;
  MRI *mri = 0;

  int width, height, depth;
  int nframes;
  int type;

  // use stat to do this.
  if ((succeed = checkBrukerFiles(fname, methodFile, acqpFile, dataFile, d3procFile, recoFile, 1)) == 0)
    return (MRI *)0;

  ///////////////////////////////////////////////////////////////////////////
  // Bruker pdata/1/reco
  //        pdata/1/d3proc
  // contain the 2dseq bshort or bfloat data information
  ///////////////////////////////////////////////////////////////////////////
  // read reco to get reconstruction info
  if ((succeed = readBrukerReco(recoFile, &bTran)) == 0) return (MRI *)0;
  // read reco width, height depth, type, and nframes
  if ((succeed = readBrukerD3proc(d3procFile, &width, &height, &depth, &type, &nframes)) == 0) return (MRI *)0;

  // read acqp to get actual acquisition information
  if ((succeed = readBrukerAcqp(acqpFile, &TR, &TE, &TI, &flip_angle, &bTran)) == 0) return (MRI *)0;

  // now allocate MRI
  type = bTran.type;
  mri = MRIalloc(width, height, depth, type);
  mri->tr = TR;
  mri->te = TE;
  mri->ti = TI;
  mri->flip_angle = flip_angle;

  // create voxToRas transform
  buildVoxToRASTransform(mri, &bTran);

  // now ready to read volume
  if (read_volume) {
    if ((succeed = readBrukerVolume(mri, dataFile)) == 0) {
      if (mri) MRIfree(&mri);
      return (MRI *)0;
    }
  }
  return mri;
}

int checkBrukerFiles(
    char *fname, char *methodFile, char *acqpFile, char *dataFile, char *d3procFile, char *recoFile, int flag)
{
  struct stat stat_buf;

  if (stat(fname, &stat_buf) < 0) {
    if (flag) fprintf(stderr, "ERROR: could not stat %s.\n", fname);
    return 0;
  }
  // first check fname is a directory
  if (S_ISDIR(stat_buf.st_mode) == 0) {
    if (flag) fprintf(stderr, "ERROR: %s is not a directory.\n", fname);
    return 0;
  }
  if (flag) printf("INFO: directory %s \n", fname);
  // build up methodFile name
  strcpy(methodFile, fname);
  strcat(methodFile, "/method");
  if (stat(methodFile, &stat_buf) < 0) {
    if (flag) fprintf(stderr, "ERROR: could not stat %s.\n", methodFile);
    return 0;
  }
  if (flag) printf("INFO:    method %s\n", methodFile);
  // build up acqpFile name
  strcpy(acqpFile, fname);
  strcat(acqpFile, "/acqp");
  if (stat(acqpFile, &stat_buf) < 0) {
    if (flag) fprintf(stderr, "ERROR: could not stat %s.\n", acqpFile);
    return 0;
  }
  if (flag) printf("INFO:      acqp %s\n", acqpFile);
  // build up dataFile name

  strcpy(dataFile, fname);
  strcat(dataFile, "/pdata/1/2dseq");
  if (stat(dataFile, &stat_buf) < 0) {
    if (flag) fprintf(stderr, "ERROR: could not stat %s.\n", dataFile);
    return 0;
  }
  if (flag) printf("INFO:     2dseq %s \n", dataFile);
  // build up d3proc file
  strcpy(d3procFile, fname);
  strcat(d3procFile, "/pdata/1/d3proc");
  if (stat(d3procFile, &stat_buf) < 0) {
    if (flag) fprintf(stderr, "ERROR: could not stat %s.\n", d3procFile);
    return 0;
  }
  if (flag) printf("INFO:    d3proc %s \n", d3procFile);
  // build up reco file
  strcpy(recoFile, fname);
  strcat(recoFile, "/pdata/1/reco");
  if (stat(recoFile, &stat_buf) < 0) {
    if (flag) fprintf(stderr, "ERROR: could not stat %s.\n", recoFile);
    return 0;
  }
  if (flag) printf("INFO:    d3proc %s \n", recoFile);

  return 1;
}

int splitParameterValue(char *sWholeLine, char *sParameter, char *sValue)
{
  char *p0, *p1;

  // check to make sure that ## is at the beginning
  p0 = strstr(sWholeLine, "##");
  if (!p0) return (0); /* ignore line */
  //
  p0 += 2; /* advance past ## */
  // locate '='
  p1 = strchr(p0, '=');
  if (!p1) return (0);
  *p1 = '\0';  // mark end of Parameter
  // copy parameter
  strcpy(sParameter, p0);
  // now points to Value
  p1++;
  // copy Value
  strcpy(sValue, p1); /* initialize sValue */
  /* Eliminate parentheses and CR in the value string. */
  for (p0 = sValue; *p0; p0++) {
    if (*p0 == '(' || *p0 == ')')
      *p0 = ' ';
    else if (*p0 == '\n')
      *p0 = '\0';
  }
  return 1;
}

int readBrukerD3proc(char *d3procFile, int *px, int *py, int *pz, int *ptype, int *pnframes)
{
  FILE *fp = 0;
  char line[512];
  char Value[128];
  char Parameter[256];
  // int lRead;

  fp = fopen(d3procFile, "r");
  if (fp == 0) {
    fprintf(stderr, "ERROR: could not open d3proc %s", d3procFile);
    return 0;
  }
  while (fgets(line, sizeof(line), fp)) {
    if (!splitParameterValue(line, Parameter, Value)) continue;

    // now gets the values
    if (!strcmp(Parameter, "END")) break;
    // get volume size
    else if (!strcmp(Parameter, "$IM_SIX"))
      // lRead = 
      sscanf(Value, "%d", px);
    else if (!strcmp(Parameter, "$IM_SIY"))
      // lRead = 
      sscanf(Value, "%d", py);
    else if (!strcmp(Parameter, "$IM_SIZ"))
      // lRead = 
      sscanf(Value, "%d", pz);
    else if (!strcmp(Parameter, "$IM_SIT")) {
      sscanf(Value, "%d", pz);
    }
    else if (!strcmp(Parameter, "$IM_SIT")) {
      sscanf(Value, "%d", pnframes);
      if (*pnframes > 1) {
        fprintf(stderr, "ERROR: nframes %d but one is supported.\n", *pnframes);
        return 0;
      }
    }
    else if (!strcmp(Parameter, "$DATTYPE")) {
      if (strcmp(Value, "ip_short") == 0)
        *ptype = MRI_SHORT;
      else if (strcmp(Value, "ip_int") == 0)
        *ptype = MRI_INT;
      else {
        fprintf(stderr,
                "ERROR: currently only short and int types are "
                "supported (Value=%s).\n",
                Value);
        return 0;
      }
    }
  }
  fclose(fp);
  return 1;
}

int buildVoxToRASTransform(MRI *mri, BrukerTransform *pTran)
{
  MATRIX *transposMatrix;
  MATRIX *voxmmMatrix;
  MATRIX *swapMatrix;
  MATRIX *translationMatrix;
  MATRIX *rotationMatrix;
  MATRIX *tmp;
  MATRIX *rottranMatrix;
  MATRIX *voxInvMatrix;
  VECTOR *rcs, *ctr;

  transposMatrix = MatrixAlloc(4, 4, MATRIX_REAL);
  voxmmMatrix = MatrixAlloc(4, 4, MATRIX_REAL);
  swapMatrix = MatrixAlloc(4, 4, MATRIX_REAL);
  translationMatrix = MatrixAlloc(4, 4, MATRIX_REAL);
  rotationMatrix = MatrixAlloc(4, 4, MATRIX_REAL);
  voxInvMatrix = MatrixAlloc(4, 4, MATRIX_REAL);
  //
  if (pTran->transposition)
    stuff_four_by_four(transposMatrix, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  else
    stuff_four_by_four(transposMatrix, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  //
  stuff_four_by_four(voxmmMatrix,
                     pTran->vox_size[0],
                     0,
                     0,
                     pTran->offset[0],
                     0,
                     pTran->vox_size[1],
                     0,
                     pTran->offset[1],
                     0,
                     0,
                     pTran->vox_size[2],
                     pTran->offset[2],
                     0,
                     0,
                     0,
                     1);
  //
  stuff_four_by_four(swapMatrix, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  //
  stuff_four_by_four(translationMatrix,
                     1,
                     0,
                     0,
                     pTran->read_offset,
                     0,
                     1,
                     0,
                     pTran->phase1_offset,
                     0,
                     0,
                     1,
                     pTran->slice_offset,
                     0,
                     0,
                     0,
                     1);
  //
  stuff_four_by_four(rotationMatrix,
                     pTran->grad_matrix[0],
                     pTran->grad_matrix[3],
                     pTran->grad_matrix[6],
                     0,
                     pTran->grad_matrix[1],
                     pTran->grad_matrix[4],
                     pTran->grad_matrix[7],
                     0,
                     pTran->grad_matrix[2],
                     pTran->grad_matrix[5],
                     pTran->grad_matrix[8],
                     0,
                     0,
                     0,
                     0,
                     1);

  tmp = MatrixMultiply(voxmmMatrix, transposMatrix, NULL);
  MatrixMultiply(swapMatrix, tmp, tmp);
  MatrixMultiply(translationMatrix, tmp, tmp);
  MatrixMultiply(rotationMatrix, tmp, tmp);
  //
  ctr = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(ctr, (double)mri->width / 2.0, (double)mri->height / 2.0, (double)mri->depth / 2.0, 1.);
  rcs = MatrixMultiply(tmp, ctr, NULL);

  mri->c_r = V3_X(rcs);
  mri->c_a = V3_Y(rcs);
  mri->c_s = V3_Z(rcs);

  //
  if (pTran->transposition)
    stuff_four_by_four(voxInvMatrix,
                       1. / pTran->vox_size[1],
                       0,
                       0,
                       0,
                       0,
                       1. / pTran->vox_size[0],
                       0,
                       0,
                       0,
                       0,
                       1. / pTran->vox_size[2],
                       0,
                       0,
                       0,
                       0,
                       1);
  else
    stuff_four_by_four(voxInvMatrix,
                       1. / pTran->vox_size[0],
                       0,
                       0,
                       0,
                       0,
                       1. / pTran->vox_size[1],
                       0,
                       0,
                       0,
                       0,
                       1. / pTran->vox_size[2],
                       0,
                       0,
                       0,
                       0,
                       1);

  // now obtain only the rotation and the translation part in MGH way
  rottranMatrix = MatrixMultiply(voxInvMatrix, tmp, NULL);
  if (pTran->dim == 3) {
    mri->x_r = *MATRIX_RELT(rottranMatrix, 1, 1);
    mri->x_a = *MATRIX_RELT(rottranMatrix, 2, 1);
    mri->x_s = *MATRIX_RELT(rottranMatrix, 3, 1);

    mri->y_r = *MATRIX_RELT(rottranMatrix, 1, 2);
    mri->y_a = *MATRIX_RELT(rottranMatrix, 2, 2);
    mri->y_s = *MATRIX_RELT(rottranMatrix, 3, 2);

    mri->z_r = *MATRIX_RELT(rottranMatrix, 1, 3);
    mri->z_a = *MATRIX_RELT(rottranMatrix, 2, 3);
    mri->z_s = *MATRIX_RELT(rottranMatrix, 3, 3);
  }
  else  // 2d case just fake
  {
    mri->x_r = -1;
    mri->x_a = 0;
    mri->x_s = 0;
    mri->y_r = 0;
    mri->y_a = 0;
    mri->y_s = -1;
    mri->z_r = 0;
    mri->z_a = 1;
    mri->z_s = 0;
  }

  if (pTran->transposition == 1) {
    mri->xsize = pTran->vox_size[1];
    mri->ysize = pTran->vox_size[0];
    mri->zsize = pTran->vox_size[2];
  }
  else {
    mri->xsize = pTran->vox_size[0];
    mri->ysize = pTran->vox_size[1];
    mri->zsize = pTran->vox_size[2];
  }

  mri->ras_good_flag = 1;

  MatrixFree(&transposMatrix);
  MatrixFree(&voxmmMatrix);
  MatrixFree(&swapMatrix);
  MatrixFree(&translationMatrix);
  MatrixFree(&rotationMatrix);
  MatrixFree(&tmp);
  MatrixFree(&rottranMatrix);
  MatrixFree(&voxInvMatrix);
  VectorFree(&rcs);
  VectorFree(&ctr);

  return 1;
}

int readBrukerAcqp(char *acqpFile, double *pTR, double *pTE, double *pTI, double *pflip_angle, BrukerTransform *bTran)
{
  FILE *fp = 0;
  char line[512];
  char Parameter[256];
  char Value[128];
  // int lRead = 0;
  int dim = 0;

  fp = fopen(acqpFile, "r");
  if (fp == 0) {
    fprintf(stderr, "ERROR: could not open acqpFile %s", acqpFile);
    return 0;
  }
  while (fgets(line, sizeof(line), fp)) {
    if (!splitParameterValue(line, Parameter, Value)) continue;

    // now gets the values
    if (!strcmp(Parameter, "END"))
      break;
    else if (!strcmp(Parameter, "OWNER")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "warning, no lines read by fgets\n");
      }
      printf("INFO: %s", line);
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "warning, no lines read by fgets\n");
      }
      printf("INFO: %s", line);
    }
    // another check for dimension
    else if (!strcmp(Parameter, "$ACQ_dim")) {
      // lRead =
      sscanf(Value, "%d", &dim);
      if (dim != 3) {
        printf("INFO: acqp tells the dimension is not 3 but %d\n", dim);
        printf("INFO: usually the volume has multiple 2d slices\n");
        printf("INFO: direction cosine info is meaningless.\n");
      }
    }
    // flip angle
    else if (!strcmp(Parameter, "$ACQ_flip_angle")) {
      //  lRead =
      sscanf(Value, "%lf", pflip_angle);
      // convert into radians
      *pflip_angle = *pflip_angle * 3.141592653589793 / 180.;
    }
    // TR
    else if (!strcmp(Parameter, "$ACQ_repetition_time")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: float value must follow ACQ_repetition_time");
        fclose(fp);
        return 0;
      }
      sscanf(line, "%lf", pTR);
    }
    // TE
    else if (!strcmp(Parameter, "$ACQ_echo_time")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: float value must follow ACQ_echo_time");
        fclose(fp);
        return 0;
      }
      sscanf(line, "%lf", pTE);
    }
    // TI
    else if (!strcmp(Parameter, "$ACQ_inversion_time")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: float value must follow ACQ_inversion_time");
        fclose(fp);
        return 0;
      }
      sscanf(line, "%lf", pTI);
    }
    // ACQ_read_offset
    else if (!strcmp(Parameter, "$ACQ_read_offset")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: float value must follow ACQ_read_offset");
        fclose(fp);
        return 0;
      }
      sscanf(line, "%lf", &bTran->read_offset);
    }
    else if (!strcmp(Parameter, "$ACQ_phase1_offset")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: float value must follow ACQ_phase1_offset");
        fclose(fp);
        return 0;
      }
      sscanf(line, "%lf", &bTran->phase1_offset);
    }
    else if (!strcmp(Parameter, "$ACQ_slice_offset")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: float value must follow ACQ_slice_offset");
        fclose(fp);
        return 0;
      }
      sscanf(line, "%lf", &bTran->slice_offset);
    }
    else if (!strcmp(Parameter, "$ACQ_grad_matrix")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: float value must follow ACQ_slice_offset");
        fclose(fp);
        return 0;
      }
      // for 2d image it may contain many 9 element arrays
      // for 3d image it is just one.
      sscanf(line,
             "%lf %lf %lf %lf %lf %lf %lf %lf %lf ",
             &bTran->grad_matrix[0],
             &bTran->grad_matrix[1],
             &bTran->grad_matrix[2],
             &bTran->grad_matrix[3],
             &bTran->grad_matrix[4],
             &bTran->grad_matrix[5],
             &bTran->grad_matrix[6],
             &bTran->grad_matrix[7],
             &bTran->grad_matrix[8]);
    }
  }
  fclose(fp);

  return 1;
}

int readBrukerReco(char *recoFile, BrukerTransform *pTran)
{
  FILE *fp = 0;
  char line[512];
  char Parameter[256];
  char Value[128];
  // int lRead = 0;
  int dim = 0;
  int i;

  fp = fopen(recoFile, "r");
  if (fp == 0) {
    fprintf(stderr, "ERROR: could not open recoFile %s", recoFile);
    return 0;
  }
  while (fgets(line, sizeof(line), fp)) {
    if (!splitParameterValue(line, Parameter, Value)) continue;

    // now gets the values
    if (!strcmp(Parameter, "END")) break;

    // get volume size
    else if (!strcmp(Parameter, "$RECO_transposition")) {
      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: value must follow RECO_transposition");
        fclose(fp);
        return 0;
      }

      sscanf(line, "%d", &pTran->transposition);
    }
    else if (!strcmp(Parameter, "$RECO_fov")) {
      // lRead =
      sscanf(Value, "%d", &dim);
      pTran->dim = dim;

      if (dim != 3) {
        fprintf(stderr, "INFO: fov dimension is %d. The data is not a 3D volume.\n", dim);
      }

      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: value must follow RECO_fov");
        fclose(fp);
        return 0;
      }

      if (dim == 3) {
        sscanf(line, "%lf %lf %lf", &pTran->fov[0], &pTran->fov[1], &pTran->fov[2]);
      }
      else if (dim == 2) {
        sscanf(line, "%lf %lf", &pTran->fov[0], &pTran->fov[1]);

        pTran->fov[2] = pTran->fov[0];
      }
    }
    else if (!strcmp(Parameter, "$RECO_size")) {
      // lRead =
      sscanf(Value, "%d", &dim);
      if (dim != 3) fprintf(stderr, "INFO: size dimension is %d. The data is not a 3D volume.\n", dim);

      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: value must follow RECO_size");
        fclose(fp);
        return 0;
      }
      if (dim == 3)
        sscanf(line, "%d %d %d", &pTran->size[0], &pTran->size[1], &pTran->size[2]);
      else if (dim == 2) {
        sscanf(line, "%d %d", &pTran->size[0], &pTran->size[1]);
        pTran->size[2] = pTran->size[0];  // just fake
      }
    }
    else if (!strcmp(Parameter, "$RECO_ft_size")) {
      // lRead =
      sscanf(Value, "%d", &dim);
      if (dim != 3) fprintf(stderr, "INFO: ft_size dimension is %d. The data is not a 3D volume.\n", dim);

      if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: value must follow RECO_size");
        fclose(fp);
        return 0;
      }
      if (dim == 3)
        sscanf(line, "%d %d %d", &pTran->ft_size[0], &pTran->ft_size[1], &pTran->ft_size[2]);
      else if (dim == 2) {
        sscanf(line, "%d %d", &pTran->ft_size[0], &pTran->ft_size[1]);
        pTran->ft_size[2] = pTran->ft_size[0];  // just fake
      }
    }
    else if (!strcmp(Parameter, "$RECO_wordtype")) {
      if ((strncmp(Value, "_16BIT_SGN_INT", 13) == 0))
        pTran->type = MRI_SHORT;
      else if ((strncmp(Value, "_32BIT_SGN_INT", 13) == 0))
        pTran->type = MRI_INT;
      else {
        fprintf(stderr, "INFO: unsupported data type %s\n", Value);
        return 0;
      }
    }
    else if (!strcmp(Parameter, "$RECO_mode")) {
      printf("INFO: reconstruction mode was %s\n", Value);
    }
    else if (!strcmp(Parameter, "ORIGIN")) {
      printf("INFO: software by %s\n", Value);
    }
  }
  // calculate vox_size and offset
  for (i = 0; i < 3; ++i) {
    pTran->vox_size[i] = pTran->fov[i] * 10. / pTran->size[i];
    pTran->offset[i] = -pTran->vox_size[i] * (pTran->ft_size[i] - 1) / 2.;
  }
  return 1;
}

int readBrukerVolume(MRI *mri, char *dataFile)
{
  FILE *fp = 0;
  int k, j;
  int nread;
  int swap_bytes_flag = 0;
  //  int size;

  if (!mri) {
    fprintf(stderr, "ERROR: readBrukerMethod() must be called before readBrukerVolume");
    return 0;
  }
  // save the data filename
  strcpy(mri->fname, dataFile);
  fp = fopen(dataFile, "r");
  if (fp == 0) {
    fprintf(stderr, "ERROR: could not open dataFile %s", dataFile);
    return 0;
  }
  switch (mri->type) {
    case MRI_UCHAR:
      //  size = sizeof(unsigned char);
      break;
    case MRI_SHORT:
      //   size = sizeof(short);
      break;
    case MRI_FLOAT:
      //   size = sizeof(float);
      break;
    case MRI_INT:
      //   size = sizeof(int);
      break;
    case MRI_LONG:
      //   size = sizeof(long);
      break;
    default:
      fprintf(stderr, "INFO: unknown size.  bail out\n");
      return (0);
  }

  for (k = 0; k < mri->depth; ++k) {
    for (j = 0; j < mri->height; ++j) {
      nread = fread(mri->slices[k][j], sizeof(short), mri->width, fp);
      if (nread != mri->width) {
        fclose(fp);
        MRIfree(&mri);
        mri = 0;
        return 0;
      }
      // this was not needed
      if (swap_bytes_flag) {
	std::vector<short> tmp(mri->width);
        swab(mri->slices[k][j], tmp.data(), mri->width * sizeof(short));
	memcpy(mri->slices[k][j], tmp.data(), mri->width * sizeof(short));
      }
    }
    exec_progress_callback(k, mri->depth, 0, 1);
  }

  fclose(fp);
  return 1;
}

int is_bruker(char *fname)
{
  struct stat stat_buf;
  char methodFile[512];
  char acqpFile[512];
  char dataFile[512];
  char d3procFile[512];
  char recoFile[512];

  if (stat(fname, &stat_buf) < 0) return (0);

  /* if it's a directory, it's a COR dir. */
  if (!S_ISDIR(stat_buf.st_mode)) return 0;
  // must check all these files exist or not
  return checkBrukerFiles(fname, methodFile, acqpFile, dataFile, d3procFile, recoFile, 0);
}
