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

#include "registerio.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fio.h"
#include "macros.h"  // FEQUAL
#include "proto.h"   // nint
#include "resample.h"
#include "timer.h"

extern const char *Progname;

/* ----------------------------------------------------------
  Name: regio_read_register()
  Reads a registration file.

  subject -- name of subject as found in the data base
  inplaneres -- in-plane resolution
  betplaneres -- between-plane resolution
  intensity -- for the register program
  R - matrix to convert from xyz in COR space to xyz in Volume space,
      ie, xyzVol = R*xyzCOR
  float2int - if the regfile has a line after the matrix, the string
      is passed to float2int_code(), the result of which is passed
      back as float2int. If there is no extra line, FLT2INT_TKREG
      is returned (indicating that the regfile was created by
      tkregister).
  -------------------------------------------------------------*/
int regio_read_register(const char *regfile,
                        char **subject,
                        float *inplaneres,
                        float *betplaneres,
                        float *intensity,
                        MATRIX **R,
                        int *float2int)
{
  FILE *fp;
  char tmp[1000];
  int r, c, n;
  float val;

  if (!stricmp(FileNameExtension(regfile, tmp), "LTA")) {
    LTA *lta;
    printf("regio_read_register: loading lta\n");
    lta = LTAread(regfile);
    if (lta == NULL) return (1);
    if (lta->subject[0] == 0) strcpy(lta->subject, "subject-unknown");
    *subject = (char *)calloc(strlen(lta->subject) + 2, sizeof(char));
    strcpy(*subject, lta->subject);

    *intensity = lta->fscale;
    *float2int = FLT2INT_ROUND;
    *inplaneres = lta->xforms[0].src.xsize;
    *betplaneres = lta->xforms[0].src.zsize;
    *R = TransformLTA2RegDat(lta);
    LTAfree(&lta);
    return (0);
  }

  fp = fopen(regfile, "r");
  if (fp == NULL) {
    perror("regio_read_register()");
    fprintf(stderr, "Could not open %s\n", regfile);
    return (1);
  }

  /* subject name */
  n = fscanf(fp, "%s", tmp);
  if (n != 1) {
    perror("regio_read_register()");
    fprintf(stderr, "Error reading subject from %s\n", regfile);
    fclose(fp);
    return (1);
  }
  *subject = (char *)calloc(strlen(tmp) + 2, sizeof(char));
  sprintf(*subject, "%s", tmp);

  /* in-plane resolution */
  n = fscanf(fp, "%f", inplaneres);
  if (n != 1) {
    perror("regio_read_register()");
    fprintf(stderr, "Error reading inplaneres from %s\n", regfile);
    fclose(fp);
    return (1);
  }

  /* between-plane resolution */
  n = fscanf(fp, "%f", betplaneres);
  if (n != 1) {
    perror("regio_read_register()");
    fprintf(stderr, "Error reading betplaneres from %s\n", regfile);
    fclose(fp);
    return (1);
  }

  /* intensity*/
  n = fscanf(fp, "%f", intensity);
  if (n != 1) {
    perror("regio_read_register()");
    fprintf(stderr, "Error reading intensity from %s\n", regfile);
    fclose(fp);
    return (1);
  }

  *R = MatrixAlloc(4, 4, MATRIX_REAL);
  if (*R == NULL) {
    fprintf(stderr, "regio_read_register(): could not alloc R\n");
    fclose(fp);
    return (1);
  }

  /* registration matrix */
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      n = fscanf(fp, "%f", &val);
      if (n != 1) {
        perror("regio_read_register()");
        fprintf(stderr, "Error reading R[%d][%d] from %s\n", r, c, regfile);
        fclose(fp);
        return (1);
      }
      (*R)->rptr[r + 1][c + 1] = val;
    }
  }

  /* Get the float2int method string */
  n = fscanf(fp, "%s", &tmp[0]);
  fclose(fp);

  if (n == EOF)
    *float2int = FLT2INT_TKREG;
  else {
    *float2int = float2int_code(tmp);
    if (*float2int == -1) {
      printf(
          "ERROR: regio_read_register(): float2int method %s from file %s,"
          " match not found\n",
          tmp,
          regfile);
      return (1);
    }
  }

  return (0);
}
/* -------------------------------------------------------------- */
int regio_print_register(
    FILE *fp, const char *subject, float inplaneres, float betplaneres, float intensity, const MATRIX *R, int float2int)
{
  int r, c;
  const char *f2imethod;

  if (subject == NULL)
    fprintf(fp, "subject-unknown\n");
  else if (strlen(subject) == 0)
    fprintf(fp, "subject-unknown\n");
  else
    fprintf(fp, "%s\n", subject);

  fprintf(fp, "%f\n", inplaneres);
  fprintf(fp, "%f\n", betplaneres);
  fprintf(fp, "%f\n", intensity);

  for (r = 0; r < 3; r++) {
    for (c = 0; c < 4; c++) {
      fprintf(fp, "%18.15e ", R->rptr[r + 1][c + 1]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "0 0 0 1\n");

  switch (float2int) {
    case FLT2INT_TKREG:
      f2imethod = "tkregister";
      break;
    case FLT2INT_FLOOR:
      f2imethod = "floor";
      break;
    case FLT2INT_ROUND:
      f2imethod = "round";
      break;
    default:
      printf("ERROR: regio_print_register: f2i code = %d, unrecoginized\n", float2int);
      return (1);
  }

  fprintf(fp, "%s\n", f2imethod);

  return (0);
}

/* -------------------------------------------------------------- */
int regio_write_register(const char *regfile,
                         const char *subject,
                         float inplaneres,
                         float betplaneres,
                         float intensity,
                         const MATRIX *R,
                         int float2int)
{
  FILE *fp;

  fp = fopen(regfile, "w");
  if (fp == NULL) {
    perror("regio_write_register");
    return (1);
  }

  regio_print_register(fp, subject, inplaneres, betplaneres, intensity, R, float2int);
  fclose(fp);

  return (0);
}
/* --------------------------------------------------------------
   regio_read_mincxfm() - reads a 3x4 transform as the last three
   lines of the xfmfile. Blank lines at the end will defeat it.
   If fileinfo != NULL, reads in the "fileinfo". This is the 3rd
   line in the minc xfm file. It will contain information about
   the center of the transform (-center). This is used by
   ltaMNIreadEx() to account for the non-zero center of the MNI
   talairach template. If the center infomation is not there,
   ltaMNIreadEx() assumes that the center is 0 and will produce
   the wrong transform. Thus, if one is going to write out the
   xfm, then one needs to keep track of this information when
   reading it in. regio_write_mincxfm() takes fileinfo as an
   argument. If one is not going to write out the xfm, then
   simply set fileinfo to NULL.
   -------------------------------------------------------------- */
int regio_read_mincxfm(const char *xfmfile, MATRIX **R, char **fileinfo)
{
  FILE *fp;
  char tmpstr[1000];
  int r, c, n, nlines;
  float val;

  memset(tmpstr, '\0', 1000);

  fp = fopen(xfmfile, "r");
  if (fp == NULL) {
    perror("regio_read_mincxfm");
    printf("ERROR: could not read %s\n", xfmfile);
    return (1);
  }

  fgetl(tmpstr, 900, fp); /* MNI Transform File */
  if (strncmp("MNI Transform File", tmpstr, 18)) {
    printf("ERROR: %s does not start as 'MNI Transform File'", xfmfile);
    return (1);
  }

  fgetl(tmpstr, 900, fp); /* fileinfo */
  if (fileinfo != NULL) {
    *fileinfo = strcpyalloc(tmpstr);
    printf("\n%s\n\n", *fileinfo);
  }
  else
    printf("Not reading in xfm fileinfo\n");

  // Close it and open it up again to rewind it
  fclose(fp);
  fp = fopen(xfmfile, "r");

  /* Count the number of lines */
  nlines = 0;
  while (fgets(tmpstr, 1000, fp) != NULL) nlines++;
  rewind(fp);

  /* skip all but the last 3 lines */
  for (n = 0; n < nlines - 3; n++) {
    if (!fgets(tmpstr, 1000, fp)) {
      fprintf(stderr, "regio_read_mincxfm(): could not read file\n");
    }
  }

  *R = MatrixAlloc(4, 4, MATRIX_REAL);
  if (*R == NULL) {
    fprintf(stderr, "regio_read_mincxfm(): could not alloc R\n");
    fclose(fp);
    return (1);
  }
  MatrixClear(*R);

  /* registration matrix */
  for (r = 0; r < 3; r++) { /* note: upper limit = 3 for xfm */
    for (c = 0; c < 4; c++) {
      n = fscanf(fp, "%f", &val);
      if (n != 1) {
        perror("regio_read_mincxfm()");
        fprintf(stderr, "Error reading R[%d][%d] from %s\n", r, c, xfmfile);
        fclose(fp);
        return (1);
      }
      (*R)->rptr[r + 1][c + 1] = val;
      /*printf("%7.4f ",val);*/
    }
    /*printf("\n");*/
  }
  (*R)->rptr[3 + 1][3 + 1] = 1.0;

  fclose(fp);
  return (0);
}
/* --------------------------------------------------------------
   regio_write_mincxfm() - writes a 3x4 transform in something
   like a minc xfm file. See regio_read_mincxfm() for docs on
   fileinfo.
   -------------------------------------------------------------- */
int regio_write_mincxfm(const char *xfmfile, const MATRIX *R, const char *fileinfo)
{
  FILE *fp;
  int r, c;

  fp = fopen(xfmfile, "w");
  if (fp == NULL) {
    perror("regio_write_mincxfm");
    printf("Could open %s for writing\n", xfmfile);
    return (1);
  }
  fprintf(fp, "MNI Transform File\n");
  if (fileinfo) fprintf(fp, "%% %s\n", fileinfo);
  fprintf(fp, "%% This file was created by %s\n", Progname);
  fprintf(fp, "%% %s\n", currentDateTime().c_str());
  fprintf(fp, "\n");
  fprintf(fp, "Transform_Type = Linear;\n");
  fprintf(fp, "Linear_Transform =\n");

  for (r = 0; r < 3; r++) { /* note: upper limit = 3 for xfm */
    for (c = 0; c < 4; c++) fprintf(fp, "%e ", R->rptr[r + 1][c + 1]);
    if (r != 2)
      fprintf(fp, "\n");
    else
      fprintf(fp, ";\n");
  }

  fclose(fp);
  return (0);
}
/* --------------------------------------------------------------
   regio_read_xfm4() - reads a 4x4 transform as the last four
   lines of the xfmfile. Blank lines at the end will defeat it.
   -------------------------------------------------------------- */
int regio_read_xfm4(const char *xfmfile, MATRIX **R)
{
  FILE *fp;
  char tmpstr[1000];
  int r, c, n, nlines;
  float val;

  memset(tmpstr, '\0', 1000);

  fp = fopen(xfmfile, "r");
  if (fp == NULL) {
    perror("regio_read_xfm4");
    fprintf(stderr, "Could read %s\n", xfmfile);
    return (1);
  }

  /* Count the number of lines */
  nlines = 0;
  while (fgets(tmpstr, 1000, fp) != NULL) nlines++;
  rewind(fp);

  /* skip all but the last 3 lines */
  for (n = 0; n < nlines - 4; n++) {
    if (!fgets(tmpstr, 1000, fp)) {
      fprintf(stderr, "regio_read_mincxfm(): could not read file\n");
    }
  }

  *R = MatrixAlloc(4, 4, MATRIX_REAL);
  if (*R == NULL) {
    fprintf(stderr, "regio_read_xfm4(): could not alloc R\n");
    fclose(fp);
    return (1);
  }
  MatrixClear(*R);

  /* registration matrix */
  for (r = 0; r < 3; r++) { /* note: upper limit = 3 for xfm */
    for (c = 0; c < 4; c++) {
      n = fscanf(fp, "%f", &val);
      if (n != 1) {
        perror("regio_read_xfm4()");
        fprintf(stderr, "Error reading R[%d][%d] from %s\n", r, c, xfmfile);
        fclose(fp);
        return (1);
      }
      (*R)->rptr[r + 1][c + 1] = val;
      /*printf("%7.4f ",val);*/
    }
    /*printf("\n");*/
  }
  (*R)->rptr[3 + 1][3 + 1] = 1.0;

  fclose(fp);
  return (0);
}
/* --------------------------------------------------------------
   regio_read_xfm() - reads a transform. If the file ends in .xfm,
   then the last three lines are read, and the last line of the
   transform is set to 0 0 0 1. Otherwise, the last four lines
   are read. Blank lines at the end will defeat it. This should
   be able to read tlas properly.
   -------------------------------------------------------------- */
int regio_read_xfm(const char *xfmfile, MATRIX **R)
{
  char *ext, *fileinfo;
  int err = 0;

  ext = fio_extension(xfmfile);

  if (strcmp(ext, "xfm") == 0) {
    err = regio_read_mincxfm(xfmfile, R, &fileinfo);
    free(fileinfo);
  }
  else
    err = regio_read_xfm4(xfmfile, R);

  free(ext);

  return (err);
}
/*
  write a transform that takes standard surfaces into hires
  volumes (in surface coords) into a register.dat file. See
  mris_register_to_volume.
*/
#include "error.h"
#include "mri_circulars.h"
int regio_write_surfacexform_to_register_dat(
    const MATRIX *B, const char *fname, const MRI_SURFACE *mris, const MRI *mri, const char *subject, int float2int)
{
  MATRIX *Ta, *Sa, *invTa, *A, *R, *S, *invS, *T, *m1, *m2;
  MRI *mri_surf = MRIallocHeader(mris->vg.width, mris->vg.height, mris->vg.depth, MRI_UCHAR, 1);

  MRIcopyVolGeomToMRI(mri_surf, &mris->vg);

  T = MRIxfmCRS2XYZtkreg(mri);
  S = MRIgetVoxelToRasXform(mri);
  invS = MatrixInverse(S, NULL);
  Ta = MRIxfmCRS2XYZtkreg(mri_surf);
  Sa = MRIgetVoxelToRasXform(mri_surf);
  invTa = MatrixInverse(Ta, NULL);
  A = MatrixMultiply(Sa, invTa, NULL);

  m1 = MatrixMultiply(A, B, NULL);
  m2 = MatrixMultiply(invS, m1, NULL);
  R = MatrixMultiply(T, m2, NULL);
  regio_write_register(fname, subject, mri->xsize, mri->zsize, 1, R, float2int);
  MatrixFree(&A);
  MatrixFree(&Ta);
  MatrixFree(&Sa);
  MatrixFree(&invTa);
  MatrixFree(&R);
  MatrixFree(&m1);
  MatrixFree(&m2);
  MatrixFree(&S);
  MatrixFree(&invS);
  MatrixFree(&T);
  MRIfree(&mri_surf);
  return (NO_ERROR);
}

MATRIX *regio_read_surfacexform_from_register_dat(const char *fname,
                                                  const MRI_SURFACE *mris,
                                                  const MRI *mri,
                                                  char **subject)
{
  MATRIX *Ta, *Sa, *invT, *A, *R, *S, *invSa, *T, *m1, *m2, *B;
  float pres, bres, intensity;
  int float2int;
  MRI *mri_surf = MRIallocHeader(mris->vg.width, mris->vg.height, mris->vg.depth, MRI_UCHAR, 1);

  if (regio_read_register(fname, subject, &pres, &bres, &intensity, &B, &float2int) != 0)
    ErrorReturn(NULL, (ERROR_NOFILE, "regio_read_surfacexform_from_register_dat(%s) failed", fname));

  MRIcopyVolGeomToMRI(mri_surf, &mris->vg);

  T = MRIxfmCRS2XYZtkreg(mri);
  S = MRIgetVoxelToRasXform(mri);
  Ta = MRIxfmCRS2XYZtkreg(mri_surf);
  Sa = MRIgetVoxelToRasXform(mri_surf);
  invSa = MatrixInverse(Sa, NULL);
  invT = MatrixInverse(T, NULL);
  A = MatrixMultiply(S, invT, NULL);

  m1 = MatrixMultiply(A, B, NULL);
  m2 = MatrixMultiply(invSa, m1, NULL);
  R = MatrixMultiply(Ta, m2, NULL);
  MatrixFree(&A);
  MatrixFree(&Ta);
  MatrixFree(&Sa);
  MatrixFree(&invT);
  MatrixFree(&B);
  MatrixFree(&m1);
  MatrixFree(&m2);
  MatrixFree(&S);
  MatrixFree(&invSa);
  MatrixFree(&T);
  MRIfree(&mri_surf);
  return (R);
}

/*
  MATRIX *regio_read_registermat(char *regfile)
  Just reads in the matrix and leaves the rest of the
  crap in the file.
*/
MATRIX *regio_read_registermat(const char *regfile)
{
  char *subject;
  float inplaneres, betplaneres, intensity;
  int float2int, err;
  MATRIX *R;
  err = regio_read_register(regfile, &subject, &inplaneres, &betplaneres, &intensity, &R, &float2int);
  if (err) return (NULL);
  free(subject);
  return (R);
}

char *regio_read_subject(const char *regfile)
{
  float inplaneres, betplaneres, intensity;
  int float2int;
  MATRIX *R;
  char *subject;
  regio_read_register(regfile, &subject, &inplaneres, &betplaneres, &intensity, &R, &float2int);
  MatrixFree(&R);
  return (subject);
}
