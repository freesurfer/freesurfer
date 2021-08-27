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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static int reductions = 1 ;
static int verbose = 0 ;
static int pad = 0 ;
static float thresh = 0.0 ;
static MRI *mri_like = NULL ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs, x0, y0, z0, dx, dy, dz ;
  MRI    *mri_src, *mri_dst = NULL ;
  char   *in_dir, *out_dir = NULL;
  MRI_REGION box ;

  nargs = handleVersionOption(argc, argv, "mri_extract");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    if (isdigit(*(argv[1]+1)))
      break ; // not an option - a negative number
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  /*
     command line:

     mri_extract <src_dir> x0 y0 z0 dx dy dz <dst_dir>
  */

  in_dir = argv[1] ;
  if (mri_like != NULL)
  {
    if (argc < 2)
      ErrorExit(ERROR_BADPARM,
		"usage: %s -like <template vol> <src volume>  <dst volume>", Progname) ;

    
  }
  else
  {
    if (argc < 8)
    {
      printf("usage: %s -like <template vol> <src volume>  <dst volume>\n", Progname) ;
      printf("\nor\n\n");
      printf("usage: %s <src volume> x0 y0 z0 dx dy dz <dst volume>", Progname) ;
      ErrorExit(ERROR_BADPARM, "", Progname) ;
    }

    if (sscanf(argv[2], "%d", &x0) != 1)
      ErrorExit(ERROR_BADPARM,
		"%s: could not scan x0 from '%s'", Progname, argv[2]) ;
    if (sscanf(argv[3], "%d", &y0) != 1)
      ErrorExit(ERROR_BADPARM,
		"%s: could not scan y0 from '%s'", Progname, argv[3]) ;
    if (sscanf(argv[4], "%d", &z0) != 1)
      ErrorExit(ERROR_BADPARM,
		"%s: could not scan z0 from '%s'", Progname, argv[4]) ;
    if (sscanf(argv[5], "%d", &dx) != 1)
      ErrorExit(ERROR_BADPARM,
		"%s: could not scan dx from '%s'", Progname, argv[5]) ;
    if (sscanf(argv[6], "%d", &dy) != 1)
      ErrorExit(ERROR_BADPARM,
		"%s: could not scan dy from '%s'", Progname, argv[6]) ;
    if (sscanf(argv[7], "%d", &dz) != 1)
      ErrorExit(ERROR_BADPARM,
		"%s: could not scan dz from '%s'", Progname, argv[7]) ;
    out_dir = argv[8] ;
  }


  if (verbose)
    fprintf(stderr, "reading from %s...", in_dir) ;
  mri_src = MRIread(in_dir) ;
  if (verbose)
    fprintf(stderr, "done\n") ;

  if (!mri_src)
    exit(1) ;
  
  if (mri_like)
  {
    double xv0, yv0, zv0, xv1, yv1, zv1 ;
    out_dir = argv[2] ;
    MRIvoxelToVoxel(mri_like, mri_src, 0, 0, 0, &xv0, &yv0, &zv0);
    MRIvoxelToVoxel(mri_like, mri_src, mri_like->width-1, mri_like->height-1, mri_like->depth-1, &xv1, &yv1, &zv1);
    printf("corners of template MRI point to (%2.1f, %2.1f, %2.1f) and (%2.1f, %2.1f, %2.1f)\n", xv0, yv0, zv0, xv1, yv1, zv1);

    x0 = MIN(nint(xv0), nint(xv1)) ;
    y0 = MIN(nint(yv0), nint(yv1)) ;
    z0 = MIN(nint(zv0), nint(zv1)) ;
    dx = MAX(nint(xv0), nint(xv1)) - x0 + 1;
    dy = MAX(nint(yv0), nint(yv1)) - y0 + 1;
    dz = MAX(nint(zv0), nint(zv1)) - z0 + 1;
    printf("extracting region (%d, %d, %d) size (%d, %d, %d)\n",x0,y0,z0,dx,dy,dz);
    printf("template volume size (%d, %d, %d)\n", mri_like->width,mri_like->height,mri_like->depth) ;
  }

  MRIboundingBox(mri_src, thresh, &box) ;
  if (x0 < 0)
    x0 = box.x-pad ;
  if (y0 < 0)
    y0 = box.y-pad ;
  if (z0 < 0)
    z0 = box.z-pad ;
  if (dx < 0)
    dx = box.dx+2*pad ;
  if (dy < 0)
    dy = box.dy+2*pad ;
  if (dz < 0)
    dz = box.dz+2*pad ;
  printf("using bounding box (%d, %d, %d) -> (%d, %d, %d)\n", x0, y0,z0, x0+dx-1, y0+dy-1,z0+dz-1) ;

  mri_dst = MRIextract(mri_src, NULL, x0, y0, z0, dx, dy, dz) ;

  if (!mri_dst)
    exit(1) ;

  if (verbose)
    fprintf(stderr, "\nwriting to %s", out_dir) ;
  MRIwrite(mri_dst, out_dir) ;
  if (verbose)
    fprintf(stderr, "\n") ;
  exit(0) ;
  return(0) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (stricmp(option, "like") == 0)
  {
    mri_like = MRIread(argv[2]) ;
    printf("using RAS coordinates and size of %s to compute bounding box\n",argv[2]) ;
    nargs = 1 ;
  }
  else switch (toupper(*option)) {
  case 'V':
    verbose = 1 ;
    break ;
  case 'P':
    pad = atoi(argv[2]) ;
    printf("using pad = %d\n", pad) ;
    nargs = 1 ;
    break ;
  case 'T':
    thresh = atof(argv[2]) ;
    printf("using threshold = %f\n", thresh) ;
    nargs = 1 ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &reductions) ;
    fprintf(stderr, "reducing %d times\n", reductions) ;
    nargs = 1 ;
    break ;
  case 'H':
  case '?':
  case 'U':
    printf("usage: %s -like <template vol> <src volume>  <dst volume>\n", Progname) ;
    printf("\nor\n\n");
    printf("usage: %s <src volume> x0 y0 z0 dx dy dz <dst volume>", Progname) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
