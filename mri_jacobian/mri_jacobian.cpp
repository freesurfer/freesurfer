/**
 * @brief Computes the (determinant of the) jacobian of the input non-linear morph.
 */
/*
 * Original Author: Fischl, B.
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
#include "mrimorph.h"
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
                

const char *Progname ;

static int use_log = 0 ;
static float sigma = 0 ;
static int write_areas = 0 ;
static int init = 0 ;
static int atlas = 0 ;
static int zero_mean = 0 ;
static LTA *lta = NULL ;

static void usage_exit(int code) ;
static int find_debug_node(GCA_MORPH *gcam, int origx, int origy, int origz) ;
static int init_gcam_areas(GCA_MORPH *gcam) ;
static int 	mask_invalid(GCA_MORPH *gcam,  MRI *mri) ;

int tm3dfile = 0;

int
main(int argc, char *argv[])
{
  char      **av, *out_fname ;
  int       ac, nargs ;
  GCA_MORPH *gcam ;
  int       msec, minutes, seconds ;
  Timer start ;
  MRI       *mri, *mri_jacobian, *mri_area, *mri_orig_area ;

  nargs = handleVersionOption(argc, argv, "mri_jacobian");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
    {
      nargs = get_option(argc, argv) ;
      argc -= nargs ;
      argv += nargs ;
    }

  if (argc < 4)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;
  gcam = GCAMread(argv[1]) ;
  
  if (gcam == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read input morph %s\n", Progname,argv[1]);
  if (Gx >= 0 && atlas == 0)
    find_debug_node(gcam, Gx, Gy, Gz) ;

  mri = MRIread(argv[2]) ;
  if (gcam == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read template volume %s\n", Progname,argv[2]);

  GCAMrasToVox(gcam, mri) ;
  if (init || tm3dfile)
    init_gcam_areas(gcam) ;
  if (atlas)
    {
      mri_area = GCAMwriteMRI(gcam, NULL, GCAM_AREA);
      mri_orig_area = GCAMwriteMRI(gcam, NULL, GCAM_ORIG_AREA);
    }
  else
    {
      mri_area = GCAMmorphFieldFromAtlas(gcam, mri, GCAM_AREA, 0, 0);
      mri_orig_area = GCAMmorphFieldFromAtlas(gcam, mri, GCAM_ORIG_AREA, 0, 0);
    }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_orig_area, "o.mgz") ;
      MRIwrite(mri_area, "a.mgz") ;
    }
  if (Gx > 0)
    printf("area = %2.3f, orig = %2.3f\n", 
	   MRIgetVoxVal(mri_area, Gx, Gy, Gz,0),MRIgetVoxVal(mri_orig_area,Gx,Gy,Gz,0)) ;
  if (lta)
    {
      double det ;
      if  (lta->type == LINEAR_RAS_TO_RAS)
	LTArasToVoxelXform(lta, mri, mri) ;
      det = MatrixDeterminant(lta->xforms[0].m_L) ;
      printf("correcting transform with det=%2.3f\n", det) ;
      MRIscalarMul(mri_orig_area, mri_orig_area, 1/det) ;
    }
  
  if (! FZERO(sigma))
    {
      MRI *mri_kernel, *mri_smooth ;
      mri_kernel = MRIgaussian1d(sigma, 100) ;
      mri_smooth = MRIconvolveGaussian(mri_area, NULL, mri_kernel) ;
      MRIfree(&mri_area) ; mri_area = mri_smooth ;
      mri_smooth = MRIconvolveGaussian(mri_orig_area, NULL, mri_kernel) ;
      MRIfree(&mri_orig_area) ; mri_orig_area = mri_smooth ;

      MRIfree(&mri_kernel) ; 
    }
  if (Gx > 0)
    printf("after smoothing area = %2.3f, orig = %2.3f\n", 
	   MRIgetVoxVal(mri_area, Gx, Gy, Gz,0),MRIgetVoxVal(mri_orig_area,Gx,Gy,Gz,0)) ;
  mri_jacobian = MRIdivide(mri_area, mri_orig_area, NULL) ;
  if (Gx > 0)
    printf("jacobian = %2.3f\n", MRIgetVoxVal(mri_jacobian, Gx, Gy, Gz,0)) ;
  if (atlas)
    mask_invalid(gcam, mri_jacobian) ;
  if (use_log)
    {
      MRIlog10(mri_jacobian, NULL, mri_jacobian, 0) ;
      if (zero_mean)
	MRIzeroMean(mri_jacobian, mri_jacobian) ;
      if (Gx > 0)
	printf("log jacobian = %2.3f\n", MRIgetVoxVal(mri_jacobian, Gx, Gy, Gz,0)) ;
    }
  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_jacobian, out_fname) ;
  MRIfree(&mri_jacobian) ;
  if (write_areas)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s_area.mgz", out_fname) ;
      printf("writing area to %s\n", fname) ;
      MRIwrite(mri_area, fname) ;
      sprintf(fname, "%s_orig_area.mgz", out_fname) ;
      printf("writing orig area to %s\n", fname) ;
      MRIwrite(mri_orig_area, fname) ;
    }
  if (atlas && Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      char fname[STRLEN] ;
      FileNameRemoveExtension(out_fname, out_fname) ;
      mri_area = GCAMwriteMRI(gcam, mri_area, GCAM_MEANS);
      sprintf(fname, "%s_means.mgz", out_fname) ;
      printf("writing means to %s\n", fname) ;
      MRIwrite(mri_area, fname) ;
      sprintf(fname, "%s_labels.mgz", out_fname) ;
      mri_area = GCAMwriteMRI(gcam, mri_area, GCAM_LABEL);
      printf("writing labels to %s\n", fname) ;
      MRIwrite(mri_area, fname) ;
    }
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ; minutes = seconds / 60 ;seconds = seconds % 60 ;
  fprintf(stderr, "jacobian calculation took %d minutes and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;
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
  if (!stricmp(option, "dt"))
    {
    }
  else if (!stricmp(option, "debug_voxel"))
    {
      Gx = atoi(argv[2]) ; Gy = atoi(argv[3]) ; Gz = atoi(argv[4]) ;
      nargs = 3 ;
      printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    }
  else if (!stricmp(option, "remove"))
    {
      lta = LTAread(argv[2]) ;
      if (lta == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s\n", Progname, argv[2]) ;
      printf("removing determinant of transform %s\n", argv[2]) ;
      nargs = 1 ;
    }
  else if (!stricmp(option, "tm3d"))
    {
      tm3dfile = 1;
      printf("The input morph originated from a tm3d (mri_cvs_register file).\n") ;
    }
  else switch (toupper(*option))
    {
    case 'A':
      atlas = 1 ;
      printf("outputing in atlas coords\n") ;
      break ;
    case 'W':
      write_areas = 1 ;
      printf("writing area volumes\n") ;
      break ;
    case 'L':
      use_log = 1 ;
      printf("taking log of jacobian values before saving\n") ;
      break ;
    case 'S':
      sigma = atof(argv[2]) ;
      printf("smoothing jacobian volume with sigma=%2.2f\n", sigma) ;
      nargs = 1 ;
      break ;
    case 'Z':
      zero_mean = 1 ;
      use_log = 1 ;
      printf("making log jacobian zero mean\n") ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }
  
  return(nargs) ;
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code)
{
  printf("usage: %s [options] <3d morph> <template volume> <output volume>\n", Progname) ;
  
  printf("\n") ;
  printf("Options:\n\n") ;
  printf("  -dt\n");
  printf("  -debug_voxel Gx Gy Gz\n");
  printf("  -remove\n");
  printf("  -a output is written in atlas coordinate system\n");
  printf("  -w writing area volumes \n");
  printf("  -l taking log of jacoian values before saving\n");
  printf("  -s <sigma>: smoothing jacobian volume with sigma\n");
  printf("  -z making log jacobian zero mean\n");
  printf("  -tm3d: the input morph (m3z) originated from tm3d (mri_cvs_register)\n");
  printf("  -? / -u: writing out help \n");

  exit(code) ;
}

static int
find_debug_node(GCA_MORPH *gcam, int origx, int origy, int origz)
{
  int            x, y, z, xmin, ymin, zmin ;
  double         d, dmin ;
  GCA_MORPH_NODE *gcamn ;


  dmin = 1e10 ; xmin = ymin = zmin = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
	{	
	  for (z = 0 ; z < gcam->depth ; z++)
	    {
	      gcamn = &gcam->nodes[x][y][z] ;
	      d = sqrt(SQR(gcamn->origx-origx) + SQR(gcamn->origy-origy) + SQR(gcamn->origz-origz)) ;
	      if (d < dmin)
		{
		  dmin = d ; xmin = x ; ymin = y ; zmin = z ;
		}
	    }
	}
    }
  gcamn = &gcam->nodes[xmin][ymin][zmin] ;
  printf("Talairach voxel (%d, %d, %d) maps to node (%d, %d, %d) %s --> (%2.1f, %2.1f, %2.1f)\n",
	 origx, origy, origz, xmin, ymin, zmin, cma_label_to_name(gcamn->label), gcamn->x, gcamn->y, gcamn->z) ;
		
  return(NO_ERROR) ;
}

static int
init_gcam_areas(GCA_MORPH *gcam)
{
  int            x, y, z ;
  double         orig_area ;
  GCA_MORPH_NODE *gcamn ;

  orig_area = gcam->spacing * gcam->spacing * gcam->spacing ;
  for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
	{	
	  for (z = 0 ; z < gcam->depth ; z++)
	    {
	      gcamn = &gcam->nodes[x][y][z] ;
	      gcamn->orig_area = gcamn->orig_area1 = gcamn->orig_area2 = orig_area ;
	    }
	}
    }
  return(NO_ERROR) ;
}

static int
mask_invalid(GCA_MORPH *gcam,  MRI *mri)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
	{	
	  for (z = 0 ; z < gcam->depth ; z++)
	    {
	      gcamn = &gcam->nodes[x][y][z] ;
	      if (gcamn->invalid || gcamn->area <= 0)
		MRIsetVoxVal(mri, x, y, z, 0, 1) ;
	    }
	}
    }
  return(NO_ERROR) ;
}
