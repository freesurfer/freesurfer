/**
 * @brief program for computing a distance transform on the surface
 *
 */
/*
 * Original Author: Bruce Fischl
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "label.h"
#include "annotation.h"
#include "MARS_DT_Boundary.h"


int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static float anterior_dist = -1 ;
static float posterior_dist = -1 ;
static int divide = 1 ;
static char *divide_surf_name ;
static int output_label = 0 ;
static float normalize = -1 ;
static int vol = 0 ;   // do distance from surface into volume instead of label on surface

int
main(int argc, char *argv[])
{
  char          **av, *output_fname ;
  int           ac, nargs, msec, mode=-1 ;
  LABEL         *area = NULL ;
  MRI_SURFACE   *mris ;
  Timer then ;
  MRI           *mri_dist ;

  nargs = handleVersionOption(argc, argv, "mris_distance_transform");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;


  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  then.reset() ;
  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s",
              Progname, argv[1]) ;

  if (vol)
  {
/*
    mri_template = MRIread(argv[2]) ;
    if (!mri_template)
      ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume from %s\n", Progname, argv[2]) ;
*/
  }
  else
  {
    area = LabelRead(NULL, argv[2]) ;
    if (area == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s",
		Progname, argv[2]) ;
    
    if (anterior_dist > 0)
      LabelCropAnterior(area, anterior_dist) ;
    if (posterior_dist > 0)
      LabelCropPosterior(area, posterior_dist) ;
  }
  
  if (stricmp(argv[3], "signed") == 0)
    mode = DTRANS_MODE_SIGNED ;
  else if (stricmp(argv[3], "unsigned") == 0)
    mode = DTRANS_MODE_UNSIGNED ;
  else if (stricmp(argv[3], "outside") == 0)
    mode = DTRANS_MODE_OUTSIDE ;
  else
  {
    print_usage() ;
    ErrorExit(ERROR_BADPARM, "unrecognized mode choice %s\n", argv[3]) ;
  }
  output_fname = argv[4] ;

  MRIScomputeMetricProperties(mris) ;
  if (vol)
  {
    mri_dist = MRIScomputeDistanceToSurface(mris, NULL, 0.25) ;
    MRIwrite(mri_dist, argv[4]) ;
  }
  else
  {
    MRIScomputeSecondFundamentalForm(mris) ;
    if (normalize > 0)
    {
      normalize = sqrt(mris->total_area) ;
      printf("normalizing surface distances by sqrt(%2.1f) = %2.1f\n", mris->total_area,normalize) ;
    }
    if (divide > 1)
    {
      int  i ;
      char fname[STRLEN], ext[STRLEN], base_name[STRLEN] ;
      LABEL *area_division ;
      
      FileNameExtension(output_fname, ext) ;
      FileNameRemoveExtension(output_fname, base_name) ;
      LabelMark(area, mris) ;
      MRIScopyMarksToAnnotation(mris) ;
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      if (MRISreadVertexPositions(mris, divide_surf_name) != NO_ERROR)
	ErrorExit(ERROR_BADPARM, "%s: could not read vertex coords from %s", Progname, divide_surf_name) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISdivideAnnotationUnit(mris, 1, divide) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      
      
      // MRISdivideAnnotationUnit sets the marked to be in [0,divide-1], make it [1,divide]
      // make sure they are oriented along original a/p direction
#define MAX_UNITS 100    
      {
	double cx[MAX_UNITS], cy[MAX_UNITS], cz[MAX_UNITS], min_a ;
	int    index, num[MAX_UNITS], new_index[MAX_UNITS], j, min_i ;
	VERTEX *v ;
	
	memset(num, 0, sizeof(num[0])*divide) ;
	memset(cx, 0, sizeof(cx[0])*divide) ;
	memset(cy, 0, sizeof(cy[0])*divide) ;
	memset(cz, 0, sizeof(cz[0])*divide) ;
	for (i = 0 ; i < area->n_points ; i++)
	{
	  if (area->lv[i].vno < 0 || area->lv[i].deleted > 0)
	    continue ;
	  v = &mris->vertices[area->lv[i].vno] ;
	  v->marked++ ;
	  index = v->marked ;
	  cx[index] += v->x ;
	  cy[index] += v->y ;
	  cz[index] += v->z ;
	  num[index]++ ;
	}
	memset(new_index, 0, sizeof(new_index[0])*divide) ;
	for (i = 1 ; i <= divide ; i++)
	  cy[i] /= num[i] ;
	
	// order them from posterior to anterior
	for (j = 1 ; j <= divide ; j++)
	{
	  min_a = 1e10 ; min_i = 0 ;
	  for (i = 1 ; i <= divide ; i++)
	  {
	    if (cy[i] < min_a)
	    {
	      min_a = cy[i] ;
	      min_i = i ;
	    }
	  }
	  cy[min_i] = 1e10 ;  // make it biggest so it won't be considered again
	  new_index[j] = min_i ;
	}
	for (i = 0 ; i < area->n_points ; i++)
	{
	  if (area->lv[i].vno < 0 || area->lv[i].deleted > 0)
	    continue ;
	  v = &mris->vertices[area->lv[i].vno] ;
	  v->marked = new_index[v->marked] ;
	}
      }
      for (i = 1 ; i <= divide ; i++)
      {
	area_division = LabelFromMarkValue(mris, i) ;
	
	printf("performing distance transform on division %d with %d vertices\n", 
	       i, area_division->n_points) ;
	if (output_label)
	{
	  sprintf(fname, "%s%d.label", base_name, i) ;
	  printf("writing %dth subdivision to %s\n", i, fname) ;
	  LabelWrite(area_division, fname);
	}
	MRISdistanceTransform(mris, area_division, mode) ;
	sprintf(fname, "%s%d.%s", base_name, i, ext) ;
	if (normalize > 0)
	  MRISmulVal(mris, 1.0/normalize) ;
	MRISwriteValues(mris, fname) ;
      }
    }
    else
    {
      MRISdistanceTransform(mris, area, mode) ;
      if (normalize > 0)
	MRISmulVal(mris, 1.0/normalize) ;
      MRISwriteValues(mris, output_fname) ;
    }
  }

  msec = then.milliseconds() ;
  fprintf(stderr,"distance transform took %2.1f minutes\n", (float)msec/(60*1000.0f));

  exit(0) ;
  return(0) ;  /* for ansi */
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "anterior"))
  {
    anterior_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("using anterior-most %2.3f mm of label only\n", anterior_dist) ;
  }
  else if (!stricmp(option, "vol"))
  {
    vol = 1 ;
    printf("computing distance to surface in volume instead of label distance\n") ;
  }
  else if (!stricmp(option, "normalize"))
  {
    normalize = 1 ;
    printf("normalizing distances with surface area\n") ;
  }
  else if (!stricmp(option, "olabel"))
  {
    output_label = 1 ;
    printf("writing out label subdivisions\n") ;
  }
  else if (!stricmp(option, "divide"))
  {
    divide = atoi(argv[2]) ;
    divide_surf_name = argv[3] ;
    nargs = 2 ;
    printf("dividing label into %d parts along primary eigendirection using coords in %s\n", 
           divide, divide_surf_name) ;
  }
  else if (!stricmp(option, "posterior"))
  {
    posterior_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("using posterior-most %2.3f mm of label only\n", posterior_dist) ;
  }
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_help() ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  printf("usage: %s [options] <surface> <label> <mode> <output file>\n",
         Progname) ;
  printf("where mode is one of 'signed', 'unsigned', 'outside'\n") ;
  
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "This program computes the distance transform of a label on "
          "the surface\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "\t-anterior <dist>   - only use anteriormost <dist> portion of label\n") ;
  fprintf(stderr,
          "\t-posterior <dist>  - only use posteriormost <dist> portion of label\n") ;
  fprintf(stderr,
          "\t-divide <n>        - divide label into <n> units along primary eigendirection\n") ;
  fprintf(stderr,
          "\t-olabel            - output label subdivisions\n") ;
  exit(1) ;
}


static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

#if 0
static int
crop_anterior_label(LABEL *area, float anterior_dist)
{
  int    n ;
  float  amax ;

  amax = -100000;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y > amax)
      amax = area->lv[n].y ;
  }
  
  amax -= anterior_dist ;
  printf("cropping all vertices with Y > %2.0f\n", amax) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y < amax)
      area->lv[n].deleted = 1 ;
  }
  return(NO_ERROR) ;
}
static int
crop_posterior_label(LABEL *area, float anterior_dist)
{
  int    n ;
  float  amax ;

  amax = -100000;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y > amax)
      amax = area->lv[n].y ;
  }
  
  amax += anterior_dist ;
  printf("cropping all vertices with Y > %2.0f\n", amax) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y > amax)
      area->lv[n].deleted = 1 ;
  }
  return(NO_ERROR) ;
}

#endif
