/**
 * @brief utilities handling control points
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

#include "ctrpoints.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pwd.h>
#include <time.h>
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "label.h"
#include "mri.h"
#include "proto.h"
#include "transform.h"
#include "timer.h"
#include "utils.h"  //  fgetl


MPoint *MRIreadControlPoints(const char *fname, int *count, int *useRealRAS)
{
  FILE *fp;
  char *cp, line[STRLEN];
  int i = 0;
  float xw, yw, zw;
  char text[256];
  int val;
  // int numpoints;
  int num_control_points, nargs;
  MPoint *pointArray = 0;
  char extension[STRLEN];

  FileNameExtension(fname, extension);
  if (!stricmp(extension, "label")) {
    LABEL *area;

    area = LabelRead(NULL, fname);
    if (area == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "MRIreadControlPoints: could not load label file %s", fname));

    *count = num_control_points = area->n_points;

    // allocate memory
    pointArray = (MPoint *)malloc(num_control_points * sizeof(MPoint));
    if (!pointArray)
      ErrorExit(ERROR_NOMEMORY, "MRIreadControlPoints could not allocate %d-sized array", num_control_points);
    *useRealRAS = area->coords == LABEL_COORDS_SCANNER_RAS;
    for (i = 0; i < num_control_points; i++) {
      pointArray[i].x = (double)area->lv[i].x;
      pointArray[i].y = (double)area->lv[i].y;
      pointArray[i].z = (double)area->lv[i].z;
    }
    LabelFree(&area);
    return (pointArray);
  }

  *useRealRAS = 0;
  if (Gdiag & DIAG_SHOW) fprintf(stderr, "reading control points from %s...\n", fname);

  //
  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_BADPARM, "cannot open control point file", fname));

  // get number of points
  num_control_points = 0;
  do {
    cp = fgetl(line, 199, fp);

    if (cp) {
      i = sscanf(cp, "%f %f %f", &xw, &yw, &zw);
      if (i == 3) {
        num_control_points++;
      }
      else  // new format
      {
        i = sscanf(cp, "%s %d", text, &val);
        // if (strcmp("numpoints", text) == 0 && i == 2) {
          // numpoints = val;
        // }
        // else 
        if (strcmp("useRealRAS", text) == 0 && i == 2) {
          *useRealRAS = val;
        }
      }
    }
  } while (cp);

  fprintf(stderr, "Reading %d control points...\n", num_control_points);
  *count = num_control_points;

  // allocate memory
  pointArray = (MPoint *)malloc(num_control_points * sizeof(MPoint));
  if (!pointArray)
    ErrorExit(ERROR_NOMEMORY, "MRIreadControlPoints could not allocate %d-sized array", num_control_points);
  rewind(fp);
  for (i = 0; i < num_control_points; i++) {
    cp = fgetl(line, 199, fp);
    nargs = sscanf(cp, "%f %f %f", &xw, &yw, &zw);
    if (nargs != 3) {
      i--;  // not a control point
      continue;
    }
    pointArray[i].x = (double)xw;
    pointArray[i].y = (double)yw;
    pointArray[i].z = (double)zw;
  }
  fclose(fp);

  return pointArray;
}

int MRIwriteControlPoints(const MPoint *pointArray, int count, int useRealRAS, const char *fname)
{
  FILE *fp;
  int i;
  int res;

  if (Gdiag & DIAG_SHOW) fprintf(stderr, "Writing control points to %s...\n", fname);

  if (useRealRAS > 1 || useRealRAS < 0)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "MRIwriteControlPoints useRealRAS must"
                 " be 0 (surfaceRAS) or 1 (scannerRAS) but %d\n",
                 useRealRAS))

        fp = fopen(fname, "w");
  if (!fp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "MRIwriteControlPoints(%s): could not"
                 " open file",
                 fname));
  for (i = 0; i < count; ++i) {
    if ((res = fprintf(fp, "%f %f %f\n", pointArray[i].x, pointArray[i].y, pointArray[i].z)) < 0)
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "MRIwriteControlPoints(%s): could not"
                   " write file",
                   fname));
  }
  // if res < 0, then error
  res = fprintf(fp, "info\n");
  res = fprintf(fp, "numpoints %d\n", count);
  res = fprintf(fp, "useRealRAS %d\n", useRealRAS);
  // get user id
  char user[1024];
  struct passwd *pw = getpwuid(geteuid());
  if ((pw != NULL) && (pw->pw_name != NULL)) {
    strcpy(user, pw->pw_name);
  } else {
    strcpy(user, "unknown user");
  }
  res = fprintf(fp, "written by %s on %s\n", user, currentDateTime().c_str());
  res = fclose(fp);

  return (NO_ERROR);
}

// (mr) uses LTA to map point array:
MPoint *MRImapControlPoints(const MPoint *pointArray, int count, int useRealRAS, MPoint *trgArray, LTA *lta)
{
  if (trgArray == NULL) trgArray = (MPoint *)malloc(count * sizeof(MPoint));

  if (!lta->xforms[0].src.valid) ErrorExit(ERROR_BADPARM, "MRImapControlPoints LTA src geometry not valid!\n");
  if (!lta->xforms[0].dst.valid) ErrorExit(ERROR_BADPARM, "MRImapControlPoints LTA dst geometry not valid!\n");

  // create face src and target mri from lta:
  MRI *mri_src = MRIallocHeader(1, 1, 1, MRI_UCHAR, 1);
  useVolGeomToMRI(&lta->xforms[0].src, mri_src);
  MRI *mri_trg = MRIallocHeader(1, 1, 1, MRI_UCHAR, 1);
  useVolGeomToMRI(&lta->xforms[0].dst, mri_trg);

  // set vox ras transforms depending on flag:
  MATRIX *src_ras2vox, *src_vox2ras, *trg_vox2ras;
  switch (useRealRAS) {
    case 0:
      src_vox2ras = MRIxfmCRS2XYZtkreg(mri_src);
      src_ras2vox = MatrixInverse(src_vox2ras, NULL);
      MatrixFree(&src_vox2ras);
      trg_vox2ras = MRIxfmCRS2XYZtkreg(mri_trg);
      break;
    case 1:
      src_ras2vox = extract_r_to_i(mri_src);
      trg_vox2ras = extract_i_to_r(mri_trg);
      break;
    default:
      ErrorExit(ERROR_BADPARM, "MRImapControlPoints has bad useRealRAS flag %d\n", useRealRAS);
  }

  // make vox2vox:
  lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);

  // concatenate transforms:
  MATRIX *M = NULL;
  M = MatrixMultiply(lta->xforms[0].m_L, src_ras2vox, M);
  M = MatrixMultiply(trg_vox2ras, M, M);

  if (Gdiag_no > 0) {
    printf("MRImapControlPoints(): M -------------------\n");
    MatrixPrint(stdout, M);
    printf("------------------\n");
  }

  // clenup some stuff:
  MRIfree(&mri_src);
  MRIfree(&mri_trg);
  MatrixFree(&src_ras2vox);
  MatrixFree(&trg_vox2ras);

  // map point array
  VECTOR *p = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(p, 4) = 1.0;
  VECTOR *p2 = VectorAlloc(4, MATRIX_REAL);
  int i;
  for (i = 0; i < count; i++) {
    VECTOR_ELT(p, 1) = pointArray[i].x;
    VECTOR_ELT(p, 2) = pointArray[i].y;
    VECTOR_ELT(p, 3) = pointArray[i].z;
    MatrixMultiply(M, p, p2);
    trgArray[i].x = VECTOR_ELT(p2, 1);
    trgArray[i].y = VECTOR_ELT(p2, 2);
    trgArray[i].z = VECTOR_ELT(p2, 3);
  }

  // cleanup rest
  MatrixFree(&M);
  MatrixFree(&p);
  MatrixFree(&p2);

  return trgArray;
}

/*!
  \fn MPoint *ControlPoints2Vox(MPoint *ras, int npoints, int UseRealRAS, MRI
  *vol)
  \brief Converts control points from RAS to Voxel (col, row, slice)
  in the given volume. The col, row slice remain floating point.
  Note: UseRealRAS has not really been tested.
 */
MPoint *ControlPoints2Vox(MPoint *ras, int npoints, int UseRealRAS, MRI *vol)
{
  MPoint *crs;
  int n;
  MATRIX *vox2ras, *ras2vox, *vras, *vcrs = NULL;

  crs = (MPoint *)calloc(sizeof(MPoint), npoints);

  if (UseRealRAS)
    vox2ras = MRIxfmCRS2XYZ(vol, 0);
  else
    vox2ras = MRIxfmCRS2XYZtkreg(vol);
  ras2vox = MatrixInverse(vox2ras, NULL);
  MatrixFree(&vox2ras);

  vras = MatrixAlloc(4, 1, MATRIX_REAL);
  vras->rptr[4][1] = 1;

  for (n = 0; n < npoints; n++) {
    vras->rptr[1][1] = ras[n].x;
    vras->rptr[1][2] = ras[n].y;
    vras->rptr[1][3] = ras[n].z;
    vcrs = MatrixMultiply(ras2vox, vras, vcrs);
    crs[n].x = vcrs->rptr[1][1];
    crs[n].y = vcrs->rptr[2][1];
    crs[n].z = vcrs->rptr[3][1];
    // printf("%4.1f %4.1f %4.1f   %4.1f %4.1f %4.1f\n",
    //	   ras[n].x,ras[n].y,ras[n].z,crs[n].x,crs[n].y,crs[n].z);
  }

  MatrixFree(&ras2vox);
  MatrixFree(&vras);
  MatrixFree(&vcrs);

  return (crs);
}

/*!
  \fn MPoint *ControlPointsApplyMatrix(MPoint *srcctr, int nctrpoints, MATRIX
  *M, MPoint *outctr)
  \brief Multiplies the xyz of the control points by the given matrix
 */
MPoint *ControlPointsApplyMatrix(MPoint *srcctr, int nctrpoints, MATRIX *M, MPoint *outctr)
{
  int n;
  MATRIX *srcras = NULL, *outras = NULL;

  if (outctr == NULL) outctr = (MPoint *)calloc(sizeof(MPoint), nctrpoints);

  srcras = MatrixAlloc(4, 1, MATRIX_REAL);
  srcras->rptr[4][1] = 1;

  for (n = 0; n < nctrpoints; n++) {
    srcras->rptr[1][1] = srcctr[n].x;
    srcras->rptr[1][2] = srcctr[n].y;
    srcras->rptr[1][3] = srcctr[n].z;
    outras = MatrixMultiply(M, srcras, outras);
    outctr[n].x = outras->rptr[1][1];
    outctr[n].y = outras->rptr[2][1];
    outctr[n].z = outras->rptr[3][1];
    // printf("%4.1f %4.1f %4.1f   %4.1f %4.1f %4.1f\n",
    //	   srcctr[n].x,srcctr[n].y,srcctr[n].z,
    //	   outctr[n].x,outctr[n].y,outctr[n].z);
  }

  return (outctr);
}

/*!
  \fn MATRIX *ControlPoints2TalMatrix(char *subject)
  \brief Computes the matrix that when multiplied by the control
  points in control.dat bring the contol points into the "talairach"
  (ie, fsaverage) space.
 */
MATRIX *ControlPoints2TalMatrix(char *subject)
{
  MATRIX *Ta, *Na, *Tv, *Nv, *invNa, *invTv, *XFM, *M;
  char *SUBJECTS_DIR, tmpstr[2000];
  MRI *fsa, *vol;
  LTA *lta;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr, "%s/fsaverage/mri/orig.mgz", SUBJECTS_DIR);
  printf("%s\n", tmpstr);
  fsa = MRIreadHeader(tmpstr, MRI_VOLUME_TYPE_UNKNOWN);
  if (fsa == NULL) return (NULL);
  Na = MRIxfmCRS2XYZ(fsa, 0);
  Ta = MRIxfmCRS2XYZtkreg(fsa);
  invNa = MatrixInverse(Na, NULL);

  sprintf(tmpstr, "%s/%s/mri/orig.mgz", SUBJECTS_DIR, subject);
  printf("%s\n", tmpstr);
  vol = MRIreadHeader(tmpstr, MRI_VOLUME_TYPE_UNKNOWN);
  if (vol == NULL) return (NULL);
  Nv = MRIxfmCRS2XYZ(vol, 0);
  Tv = MRIxfmCRS2XYZtkreg(vol);
  invTv = MatrixInverse(Tv, NULL);

  sprintf(tmpstr, "%s/%s/mri/transforms/talairach.xfm", SUBJECTS_DIR, subject);
  printf("%s\n", tmpstr);
  lta = LTAreadEx(tmpstr);
  if (lta == NULL) {
    printf("ERROR: reading %s\n", tmpstr);
    return (NULL);
  }
  XFM = MatrixCopy(lta->xforms[0].m_L, NULL);
  LTAfree(&lta);

  M = MatrixMultiply(Ta, invNa, NULL);
  M = MatrixMultiply(M, XFM, M);
  M = MatrixMultiply(M, Nv, M);
  M = MatrixMultiply(M, invTv, M);

  // printf("\n");  MatrixPrint(stdout,M);  printf("\n");

  MatrixFree(&Ta);
  MatrixFree(&Na);
  MatrixFree(&Tv);
  MatrixFree(&Nv);
  MatrixFree(&invTv);
  MatrixFree(&invNa);
  MatrixFree(&XFM);
  MRIfree(&fsa);
  MRIfree(&vol);
  return (M);
}

/*!
  \fn MPoint *GetTalControlPoints(char **subjectlist, int nsubjects, int
  *pnctrtot)
  \brief Loads in the control points for each subject in the list and
  converts them to talairach (fsaverage) space. If a subject does not
  have a control.dat then it is skipped.
 */
MPoint *GetTalControlPoints(char **subjectlist, int nsubjects, int *pnctrtot)
{
  int nthsubject, nc, nctot, CPUseRealRAS, k;
  MPoint *subjctr, *fsactr, *ctr = NULL;
  char *SUBJECTS_DIR, tmpstr[2000];
  MATRIX *M;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  // First, count the total number of control points so structure can be
  // allocated
  nctot = 0;
  for (nthsubject = 0; nthsubject < nsubjects; nthsubject++) {
    sprintf(tmpstr, "%s/%s/tmp/control.dat", SUBJECTS_DIR, subjectlist[nthsubject]);
    if (!fio_FileExistsReadable(tmpstr)) continue;
    subjctr = MRIreadControlPoints(tmpstr, &nc, &CPUseRealRAS);
    nctot += nc;
    free(subjctr);
  }
  printf("  GetTalControlPoints(): nsubjects = %d, nctot = %d\n", nsubjects, nctot);
  ctr = (MPoint *)calloc(sizeof(MPoint), nctot);

  // Now load the control points for each subject and convert them "talairach",
  // ie, fsaverage space.
  nctot = 0;
  for (nthsubject = 0; nthsubject < nsubjects; nthsubject++) {
    sprintf(tmpstr, "%s/%s/tmp/control.dat", SUBJECTS_DIR, subjectlist[nthsubject]);
    if (!fio_FileExistsReadable(tmpstr)) continue;
    subjctr = MRIreadControlPoints(tmpstr, &nc, &CPUseRealRAS);
    M = ControlPoints2TalMatrix(subjectlist[nthsubject]);
    fsactr = ControlPointsApplyMatrix(subjctr, nc, M, NULL);
    for (k = 0; k < nc; k++) {
      ctr[nctot + k].x = fsactr[k].x;
      ctr[nctot + k].y = fsactr[k].y;
      ctr[nctot + k].z = fsactr[k].z;
    }
    nctot += nc;
    free(subjctr);
    free(fsactr);
    MatrixFree(&M);
  }

  // for(nc = 0; nc < nctot; nc++)
  // printf("%g %g %g\n",ctr[nc].x,ctr[nc].y,ctr[nc].z);
  // MRIwriteControlPoints(ctr, nctot, 0, "mycontrol.dat");
  *pnctrtot = nctot;
  return (ctr);
}

/*!
  \fn MPoint *GetTalControlPointsSFile(char *subjectlistfile, int *pnctrtot)
  \brief Loads in the control points for each subject in the file and
  converts them to talairach (fsaverage) space. If a subject does not
  have a control.dat then it is skipped.
 */
MPoint *GetTalControlPointsSFile(const char *subjectlistfile, int *pnctrtot)
{
  MPoint *ctr;
  FILE *fp;
  int nsubjects, r;
  char tmpstr[2000], **subjectlist;

  printf("  GetTalControlPointsSFile(): opening %s\n", subjectlistfile);
  fp = fopen(subjectlistfile, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", subjectlistfile);
    return (NULL);
  }
  nsubjects = 0;
  while (1) {
    r = fscanf(fp, "%s", tmpstr);
    if (r == EOF) break;
    nsubjects++;
  }
  fclose(fp);

  printf("   GetTalControlPointsSFile(): Found %d subjects\n", nsubjects);
  if (nsubjects == 0) {
    printf(" ERROR: no subjects found in %s\n", subjectlistfile);
    return (NULL);
  }
  subjectlist = (char **)calloc(sizeof(char *), nsubjects);

  fp = fopen(subjectlistfile, "r");
  nsubjects = 0;
  while (1) {
    r = fscanf(fp, "%s", tmpstr);
    if (r == EOF) break;
    subjectlist[nsubjects] = strcpyalloc(tmpstr);
    nsubjects++;
  }
  fclose(fp);

  ctr = GetTalControlPoints(subjectlist, nsubjects, pnctrtot);
  free(subjectlist);
  return (ctr);
}
