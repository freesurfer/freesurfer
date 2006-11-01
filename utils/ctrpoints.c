///////////////////////////////////////////////////////////////////
// ctrpoints.c
//
// originally written by y.tosa
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2006/11/01 20:17:52 $
// Revision       : $Revision: 1.4 $
//
////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "diag.h"
#include "error.h"
#include "mri.h"
#include "utils.h" //  fgetl
#include "ctrpoints.h"
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
  int num_control_points;
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
  } while (cp) ;

  fprintf(stderr, "Reading %d control points...\n", num_control_points) ;
  *count = num_control_points;

  // allocate memory
  pointArray=(MPoint*) calloc(num_control_points, sizeof(MPoint));
  if (!pointArray)
    ErrorExit(ERROR_NOMEMORY, 
              "MRIreadControlPoints could not allocate %d-sized array",
              num_control_points) ;
  rewind(fp) ;
  for (i = 0 ; i < num_control_points ; i++)
  {
    cp = fgetl(line, 199, fp) ;
    sscanf(cp, "%f %f %f", &xw, &yw, &zw) ;
    pointArray[i].x = (Real) xw; pointArray[i].y = (Real) yw; pointArray[i].z = (Real) zw;
  }
  fclose(fp) ;

  return pointArray;
}

int MRIwriteControlPoints(MPoint *pointArray, int count, int useRealRAS, char *fname)
{
  FILE *fp;
  int i;
  int res;
  time_t time;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "Writing control points to %s...\n", fname) ;

  if (useRealRAS > 1 || useRealRAS < 0)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "MRIwriteControlPoints useRealRAS must be 0 (surfaceRAS) or 1 (scannerRAS) but %d\n",
		useRealRAS))

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "MRIwriteControlPoints(%s): could not"
                 " open file", fname)) ;
  for (i=0 ; i < count; ++i)
  {
    if ((res = fprintf(fp, "%f %f %f\n", pointArray[i].x, pointArray[i].y, pointArray[i].z))< 0)
	ErrorReturn(ERROR_BADPARM, 
		    (ERROR_BADPARM, "MRIwriteControlPoints(%s): could not"
		     " writing file", fname)) ;
  }
  // if res < 0, then error
  res=fprintf(fp, "info\n");
  res=fprintf(fp, "numpoints %d\n", count);
  res=fprintf(fp, "useRealRAS %d\n", useRealRAS);
  res=fprintf(fp, "written by %s on %s\n",cuserid(0), asctime(localtime(&time)));
  res=fclose(fp);

  return (NO_ERROR);
}

