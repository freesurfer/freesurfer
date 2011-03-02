/**
 * @file  ctrpoints.c
 * @brief utilities handling control points
 *
 */
/*
 * Original Author: Y. Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:42 $
 *    $Revision: 1.9 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "diag.h"
#include "error.h"
#include "mri.h"
#include "utils.h" //  fgetl
#include "ctrpoints.h"
#include "transform.h"

extern char *cuserid(char *);

MPoint *MRIreadControlPoints(const char *fname, int *count, int *useRealRAS)
{
  FILE *fp ;
  char *cp, line[STRLEN] ;
  int  i = 0 ;
  float xw, yw, zw ;
  char text[256];
  int val;
  int numpoints;
  int num_control_points, nargs;
  MPoint *pointArray = 0;

  *useRealRAS = 0;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading control points from %s...\n", fname) ;

  //
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "cannot open control point file", fname)) ;

  // get number of points
  num_control_points = 0;
  do
  {
    cp = fgetl(line, 199, fp) ;

    if (cp)
    {
      i = sscanf(cp, "%f %f %f", &xw, &yw, &zw);
      if (i == 3)
      {
        num_control_points++ ;
      }
      else // new format
      {
        i = sscanf(cp, "%s %d", text, &val);
        if (strcmp("numpoints", text) == 0 && i==2)
        {
          numpoints = val;
        }
        else if (strcmp("useRealRAS", text) == 0 && i==2)
        {
          *useRealRAS = val;
        }
      }
    }
  }
  while (cp) ;

  fprintf(stderr, "Reading %d control points...\n", num_control_points) ;
  *count = num_control_points;

  // allocate memory
  pointArray=(MPoint*) malloc(num_control_points* sizeof(MPoint));
  if (!pointArray)
    ErrorExit(ERROR_NOMEMORY,
              "MRIreadControlPoints could not allocate %d-sized array",
              num_control_points) ;
  rewind(fp) ;
  for (i = 0 ; i < num_control_points ; i++)
  {
    cp = fgetl(line, 199, fp) ;
    nargs = sscanf(cp, "%f %f %f", &xw, &yw, &zw) ;
    if (nargs != 3)
    {
      i-- ;    // not a control point
      continue ;
    }
    pointArray[i].x = (double) xw;
    pointArray[i].y = (double) yw;
    pointArray[i].z = (double) zw;
  }
  fclose(fp) ;

  return pointArray;
}

int MRIwriteControlPoints(const MPoint *pointArray,
                          int count,
                          int useRealRAS,
                          const char *fname)
{
  FILE *fp;
  int i;
  int res;
  time_t time;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "Writing control points to %s...\n", fname) ;

  if (useRealRAS > 1 || useRealRAS < 0)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "MRIwriteControlPoints useRealRAS must"
                 " be 0 (surfaceRAS) or 1 (scannerRAS) but %d\n",
                 useRealRAS))

    fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "MRIwriteControlPoints(%s): could not"
                 " open file", fname)) ;
  for (i=0 ; i < count; ++i)
  {
    if ((res = fprintf(fp, "%f %f %f\n",
                       pointArray[i].x,
                       pointArray[i].y,
                       pointArray[i].z))< 0)
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM, "MRIwriteControlPoints(%s): could not"
                   " write file", fname)) ;
  }
  // if res < 0, then error
  res=fprintf(fp, "info\n");
  res=fprintf(fp, "numpoints %d\n", count);
  res=fprintf(fp, "useRealRAS %d\n", useRealRAS);
  res=fprintf(fp, "written by %s on %s\n",
              cuserid(0), asctime(localtime(&time)));
  res=fclose(fp);

  return (NO_ERROR);
}

// (mr) uses LTA to map point array:
MPoint *MRImapControlPoints(const MPoint *pointArray, int count, int useRealRAS,
                            MPoint *trgArray, LTA* lta)
{
	if (trgArray == NULL)
     trgArray=(MPoint*) malloc(count* sizeof(MPoint));

  if (! lta->xforms[0].src.valid )
    ErrorExit(ERROR_BADPARM,"MRImapControlPoints LTA src geometry not valid!\n");
  if (! lta->xforms[0].dst.valid )
    ErrorExit(ERROR_BADPARM,"MRImapControlPoints LTA dst geometry not valid!\n");

  // create face src and target mri from lta:
	MRI* mri_src = MRIallocHeader(1,1,1,MRI_UCHAR);
	useVolGeomToMRI(&lta->xforms[0].src,mri_src);
	MRI* mri_trg = MRIallocHeader(1,1,1,MRI_UCHAR);
	useVolGeomToMRI(&lta->xforms[0].dst,mri_trg);

  // set vox ras transforms depending on flag:
  MATRIX * src_ras2vox, *trg_vox2ras;
  switch (useRealRAS)
  {
    case 0:
		{
      MATRIX * src_vox2ras = MRIxfmCRS2XYZtkreg(mri_src);
      src_ras2vox = MatrixInverse(src_vox2ras,NULL);	
	  	MatrixFree(&src_vox2ras);
      trg_vox2ras = MRIxfmCRS2XYZtkreg(mri_trg);
      break;
		}
    case 1:
		{
      src_ras2vox = extract_r_to_i(mri_src);
      trg_vox2ras = extract_i_to_r(mri_trg);
      break;
		}
    default:
      ErrorExit(ERROR_BADPARM,
                "MRImapControlPoints has bad useRealRAS flag %d\n",
                useRealRAS) ;
  }
	
	// make vox2vox:
	lta = LTAchangeType(lta,LINEAR_VOX_TO_VOX);

  // concatenate transforms:
  MATRIX *M = NULL;
	M = MatrixMultiply(lta->xforms[0].m_L, src_ras2vox, M);
	M = MatrixMultiply(trg_vox2ras, M, M);
	
	// clenup some stuff:
	MRIfree(&mri_src);
	MRIfree(&mri_trg);
	MatrixFree(&src_ras2vox);
	MatrixFree(&trg_vox2ras);
	
	// map point array
  VECTOR * p  = VectorAlloc(4,MATRIX_REAL);
	VECTOR_ELT(p,4)=1.0;
  VECTOR * p2 = VectorAlloc(4,MATRIX_REAL);
	int i;
	for (i=0;i<count;i++)
	{
    VECTOR_ELT(p,1)=pointArray[i].x;
		VECTOR_ELT(p,2)=pointArray[i].y;
		VECTOR_ELT(p,3)=pointArray[i].z;
		MatrixMultiply(M, p, p2) ;
    trgArray[i].x=VECTOR_ELT(p2,1);
    trgArray[i].y=VECTOR_ELT(p2,2);
    trgArray[i].z=VECTOR_ELT(p2,3);
	}
	
	// cleanup rest
	MatrixFree(&M);
	MatrixFree(&p);
	MatrixFree(&p2);
	
  return trgArray;
}
