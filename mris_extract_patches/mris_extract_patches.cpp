/**
 * @brief extract volume patches from a surface and a label file
 *
 * Extract volumetric patches around each labeled vertex and the corresponding vertex in the other hemi
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "fsinit.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static const char *surf_name = "white" ;
static const char *sphere_name = "sphere.d1.left_right";
static const char *hemi_name = "lh" ;
static const char *ohemi_name = "rh" ;
static int hemi = LEFT_HEMISPHERE ;
static int wsize = 32 ; 
static int nbrs = 3 ;
static const char *label_name = "FCD";
static const char *vol_name = "norm.mgz" ;
static const char *ovol_name = "norm.mgz" ;
static char sdir[STRLEN] = "" ;
static float random_patch_pct = 0.0 ;
static int augment = 0 ;

MRI *MRISextractVolumeWindow(MRI_SURFACE *mris, MRI *mri, int wsize, int vno, double theta) ;

int
main(int argc, char *argv[])
{
  int          nargs, a ;
  char         *subject, fname[STRLEN], *out_dir, fname_only[STRLEN] ;
  int          msec, minutes, seconds, n ;
  Timer start ;
  MRI_SURFACE  *mris, *mris_ohemi ;
  MRI          *mri_norm, *mri_patches, *mri_labels, *mri_onorm ;
  LABEL        *area_tmp, *area ;
  int          random_patches, npoints, *non_fcd_vertices = NULL, augment_patches ;

  nargs = handleVersionOption(argc, argv, "mris_extract_patches");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (strlen(sdir) == 0)
  {
    char *env = getenv("SUBJECTS_DIR") ;
    if (env == NULL)
      ErrorExit(ERROR_UNSUPPORTED, "%s: SUBJECTS_DIR must be specified on command line with -sdir or in env", Progname) ;
    strcpy(sdir, env) ;
  }
  Progname = argv[0] ;
  FSinit() ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;
  setRandomSeed(-1L) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc != 3)
    usage_exit(1) ;

  subject = argv[1] ;
  out_dir = argv[2] ;

  printf("processing subject %s hemi %s, label %s and writing results to %s\n", subject, hemi_name, label_name, out_dir) ;
  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, subject, hemi_name, surf_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: MRISread(%s) failed", Progname, fname);

  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, subject, ohemi_name, surf_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris_ohemi = MRISread(fname) ;
  if (!mris_ohemi)
    ErrorExit(ERROR_NOFILE, "%s: MRISread(%s) failed", Progname, fname);
  
  MRIScomputeMetricProperties(mris) ;
  if (MRISreadCanonicalCoordinates(mris, sphere_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: MRISreadCanonicalCoordinates(%s) failed", Progname, sphere_name);

  MRIScomputeMetricProperties(mris_ohemi) ;
  if (MRISreadCanonicalCoordinates(mris_ohemi, sphere_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: MRISreadCanonicalCoordinates(%s) failed", Progname, sphere_name);
  MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
  MRISsetNeighborhoodSizeAndDist(mris_ohemi, nbrs) ;
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeMetricProperties(mris_ohemi) ;
  MRIScomputeSecondFundamentalForm(mris);
  MRIScomputeSecondFundamentalForm(mris_ohemi);

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, subject, vol_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri_norm = MRIread(fname) ;
  if (mri_norm == NULL)
    ErrorExit(ERROR_NOFILE, "%s: MRIread(%s) failed", Progname, fname);

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, subject, ovol_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri_onorm = MRIread(fname) ;
  if (mri_onorm == NULL)
    ErrorExit(ERROR_NOFILE, "%s: MRIread(%s) failed", Progname, fname);

  req = snprintf(fname, STRLEN, "%s/%s/label/%s.%s.label",
		 sdir, subject, hemi_name, label_name);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  area_tmp = LabelRead(subject, fname) ;
  if (area_tmp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: LabelRead(%s) failed", Progname, fname) ;
  
  LabelUnassign(area_tmp) ;
  area = LabelSampleToSurface(mris, area_tmp, mri_norm, CURRENT_VERTICES) ;
  LabelFree(&area_tmp) ;
  random_patches = round(area->n_points*random_patch_pct) ;
  npoints = area->n_points + random_patches ;
  augment_patches = nint(npoints*augment) ;
  
  mri_patches = MRIallocSequence(wsize, wsize, wsize,MRI_FLOAT, npoints+augment_patches);
  mri_labels = MRIallocSequence(npoints+augment_patches, 2, 1, MRI_INT, 1) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    MRI *mri_tmp ;
    mri_tmp = MRISextractVolumeWindow(mris, mri_norm,  wsize, area->lv[n].vno, 0) ;
    MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
    MRIfree(&mri_tmp) ;
    MRIsetVoxVal(mri_labels, n, 0, 0, 0, area->lv[n].vno) ;
    MRIsetVoxVal(mri_labels, n, 1, 0, 0, 1) ; // mark it as an FCD
  }
#define MAX_ITER 1000000
  if (random_patches > 0)
  {
    int    i, vno, iter ;
    VERTEX *v ;

    if (MRISreadCurvatureFile(mris, "sulc") != NO_ERROR)
      ErrorExit(Gerror, "");
    LabelMarkSurface(area, mris) ;
    non_fcd_vertices = (int *)calloc(random_patches, sizeof(int)) ;
    for (i = iter = 0 ; iter < MAX_ITER ; iter++)
    {
      vno = nint(randomNumber(0, mris->nvertices-1)) ;
      v = &mris->vertices[vno] ;
      if (v->curv < 0 || v->marked)
	continue ;
      v->marked = 1 ;
      non_fcd_vertices[i] = vno ;
      i++ ;
       if (i >= random_patches)
	break ;
    }
    if (iter >= MAX_ITER)
      ErrorExit(ERROR_NOFILE, "%s: max iterations %d reached without finding enough random patches", Progname, iter);
    for (i = 0 ; n < npoints ; n++, i++)
    {
      MRI *mri_tmp ;
      mri_tmp = MRISextractVolumeWindow(mris, mri_norm,  wsize,  non_fcd_vertices[i], 0) ;
      MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
      MRIfree(&mri_tmp) ;
      MRIsetVoxVal(mri_labels, n, 0, 0, 0, non_fcd_vertices[i]) ;
      MRIsetVoxVal(mri_labels, n, 1, 0, 0, 0) ; // mark it as not an FCD
    }
  }
  
  for (a = 0 ; a < augment ; a++)
  {
    double theta ;
    int n1, i ;

    theta = randomNumber(0.0, M_PI) ;
    for (n1 = 0 ; n1 < area->n_points ; n1++, n++)
    {
      MRI *mri_tmp ;

      mri_tmp = MRISextractVolumeWindow(mris, mri_norm,  wsize, area->lv[n1].vno, theta) ;
      MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
      MRIfree(&mri_tmp) ;
      MRIsetVoxVal(mri_labels, n, 0, 0, 0, area->lv[n1].vno) ;
      MRIsetVoxVal(mri_labels, n, 1, 0, 0, 1) ; // mark it as an FCD
    }
    if (random_patches > 0)
    {
      for (i = 0 ; n1 < npoints ; n1++, n++, i++)
      {
        MRI *mri_tmp ;
	mri_tmp = MRISextractVolumeWindow(mris, mri_norm,  wsize,  non_fcd_vertices[i], theta) ;
	MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
	MRIfree(&mri_tmp) ;
	MRIsetVoxVal(mri_labels, n, 0, 0, 0, non_fcd_vertices[i]) ;
	MRIsetVoxVal(mri_labels, n, 1, 0, 0, 0) ; // mark it as not an FCD
      }
    }
  }

  FileNameRemoveExtension(FileNameOnly(vol_name, fname_only), fname_only) ;
  req = snprintf(fname, STRLEN, "%s/%s.patches.%s.mgz", out_dir, hemi_name, fname_only) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("writing output file %s\n", fname) ;
  MRIwrite(mri_patches, fname) ;
  req = snprintf(fname, STRLEN, "%s/%s.labels.%s.mgz", out_dir, hemi_name, fname_only) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("writing output file %s\n", fname) ;
  MRIwrite(mri_labels, fname) ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    MRI    *mri_tmp ;
    int    ovno ;
    VERTEX *v ;

    v = &mris->vertices[area->lv[n].vno] ;
    ovno = MRISfindClosestCanonicalVertex(mris_ohemi, v->cx, v->cy, v->cz) ;
    if (ovno < 0)
      ErrorExit(ERROR_BADPARM, "%s: could not find closest vertex to %d\n", area->lv[n].vno) ;
    mri_tmp = MRISextractVolumeWindow(mris_ohemi, mri_onorm,  wsize,  ovno, 0) ;
    MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
    MRIfree(&mri_tmp) ;
    MRIsetVoxVal(mri_labels, n, 0, 0, 0, ovno) ;
    MRIsetVoxVal(mri_labels, n, 1, 0, 0, 0) ; // mark it as an FCD
  }

  if (random_patches > 0)
  {
    int    i ;

    for (i = 0 ; n < npoints ; n++, i++)
    {
      MRI    *mri_tmp ;
      VERTEX *v ;
      int    ovno ;

      v = &mris->vertices[non_fcd_vertices[i]] ;
      ovno = MRISfindClosestCanonicalVertex(mris_ohemi, v->cx, v->cy, v->cz) ;
      if (ovno < 0)
	ErrorExit(ERROR_BADPARM, "%s: could not find closest vertex to %d\n", area->lv[n].vno) ;

      mri_tmp = MRISextractVolumeWindow(mris_ohemi, mri_onorm,  wsize,  ovno, 0) ;
      MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
      MRIfree(&mri_tmp) ;
      MRIsetVoxVal(mri_labels, n, 0, 0, 0, ovno) ;
      MRIsetVoxVal(mri_labels, n, 1, 0, 0, 0) ; // mark it as not an FCD
    }
  }

  for (a = 0 ; a < augment ; a++)
  {
    double theta ;
    int n1 ;

    theta = randomNumber(0.0, M_PI) ;
    for (n1 = 0 ; n1 < area->n_points ; n1++, n++)
    {
      MRI    *mri_tmp ;
      VERTEX *v ;
      int    ovno ;

      v = &mris->vertices[area->lv[n1].vno] ;
      ovno = MRISfindClosestCanonicalVertex(mris_ohemi, v->cx, v->cy, v->cz) ;
      if (ovno < 0)
	ErrorExit(ERROR_BADPARM, "%s: could not find closest vertex to %d\n", area->lv[n1].vno) ;
      mri_tmp = MRISextractVolumeWindow(mris_ohemi, mri_onorm,  wsize,  ovno, theta) ;
      MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
      MRIfree(&mri_tmp) ;
      MRIsetVoxVal(mri_labels, n, 0, 0, 0, ovno) ;
      MRIsetVoxVal(mri_labels, n, 1, 0, 0, 0) ; // mark it as an FCD
    }
    if (random_patches > 0)
    {
      int i ;

      for (i = 0 ; n1 < npoints ; n1++, n++, i++)
      {
        MRI *   mri_tmp ;
	VERTEX  *v ;
	int     ovno ;

	v = &mris->vertices[non_fcd_vertices[i]] ;
	ovno = MRISfindClosestCanonicalVertex(mris_ohemi, v->cx, v->cy, v->cz) ;
	if (ovno < 0)
	  ErrorExit(ERROR_BADPARM, "%s: could not find closest vertex to %d\n", non_fcd_vertices[i]) ;
	
	mri_tmp = MRISextractVolumeWindow(mris_ohemi, mri_onorm,  wsize, ovno, theta) ;
	MRIcopyFrame(mri_tmp, mri_patches, 0, n) ;
	MRIfree(&mri_tmp) ;
	MRIsetVoxVal(mri_labels, n, 0, 0, 0, ovno) ;
	MRIsetVoxVal(mri_labels, n, 1, 0, 0, 0) ; // mark it as not an FCD
      }
    }
  }
  
  FileNameRemoveExtension(FileNameOnly(ovol_name, fname_only), fname_only) ;
  req = snprintf(fname, STRLEN, "%s/%s.patches.%s.mgz", out_dir, ohemi_name, fname_only) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("writing output file with %d patches to %s\n", mri_patches->nframes,fname) ;
  MRIwrite(mri_patches, fname) ;
  req = snprintf(fname, STRLEN, "%s/%s.labels.%s.mgz", out_dir, ohemi_name, fname_only) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("writing output file %s\n", fname) ;
  MRIwrite(mri_labels, fname) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stdout, "patch extraction took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;

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
  if (!stricmp(option, "sd") || !stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]);
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "sphere_name"))
  {
    sphere_name = argv[2] ;
    printf("using sphere file %s\n", sphere_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "vol_name") || !stricmp(option, "vol"))
  {
    vol_name = argv[2] ;
    printf("using volume file %s\n", vol_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ovol_name") || !stricmp(option, "ovol"))
  {
    ovol_name = argv[2] ;
    printf("using contralateral volume file %s\n", ovol_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "hemi"))
  {
    hemi_name = argv[2] ;
    if (stricmp(hemi_name, "lh")  && stricmp(hemi_name, "rh"))
      ErrorExit(ERROR_UNSUPPORTED, "%s: hemi (%s) must be either 'lh' or 'rh'\n", Progname, hemi) ;
    if (!stricmp(hemi_name, "lh"))
    {
      hemi = LEFT_HEMISPHERE ;
      ohemi_name = "rh" ;
    }
    else
    {
      hemi = RIGHT_HEMISPHERE ;
      ohemi_name = "lh" ;
    }

    nargs = 1 ;
  }
  else switch (toupper(*option))
    {
    case 'A':
      augment = atoi(argv[2]) ;
      nargs = 1 ;
      printf("augmenting the data %2d times using planar rotations\n", augment) ;
      break ;
    case 'S':
      surf_name = argv[2] ;
      nargs = 1 ;
      printf("reading surface from %s\n", surf_name) ;
      break ;
    case 'L':
      label_name = argv[2] ;
      nargs = 1 ;
      printf("reading label from %s\n", label_name) ;
      break ;
    case 'R':
      random_patch_pct = atof(argv[2]);
      nargs = 1 ;
      printf("creating %2.2f%% randomly chosen patches to augment training set\n", 100*random_patch_pct);
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    case 'W':
      wsize = atoi(argv[2]) ;
      nargs = 1 ;
      printf("setting window size to %d\n", wsize) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("debugging vertex %d\n", Gdiag_no) ;
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
  printf("Usage: %s [options] <subject> <output dir>\n",
         Progname) ;
  printf("Example: mris_extract_patches bruce $SUBJECTS_DIR/bruce/patches\n");
  exit(code) ;
}

MRI *
MRISextractVolumeWindow(MRI_SURFACE *mris, MRI *mri, int wsize, int vno, double theta)
{
  MRI    *mri_vol ;
  VERTEX *v ;
  double x0, y0, z0, whalf, xs, ys, zs, xv, yv, zv, val, e1c, e2c, nc, ctheta, stheta, x1, y1 ;
  int   xi, yi, zi ;

  ctheta = cos(theta) ;
  stheta = sin(theta) ;

  v = &mris->vertices[vno] ;
  mri_vol = MRIalloc(wsize, wsize, wsize, MRI_FLOAT) ;


  // form a window that has the vertex centered 2/3 of the way down so that more
  // of the window extends 'outwards' (in the surface normal direction), than 'inwards
  whalf = (wsize-1)/2.0;
  x0 = v->x ;
  y0 = v->y ;
  z0 = v->z ;
  for (xi = 0 ; xi < wsize ; xi++)
    for (yi = 0 ; yi < wsize ; yi++)
      for (zi = 0 ; zi < wsize ; zi++)
      {
	e1c = (xi-whalf) ;
	e2c = (yi-whalf) ;
	x1 = e1c*ctheta - e2c*stheta ;
	y1 = e1c*stheta + e2c*ctheta ;
	e1c = x1 ; e2c = y1 ;
	nc = (zi-(whalf/2.0)) ;  // have more of the window extend outwards than inwards
	xs = x0 + e1c*v->e1x + e2c*v->e2x + nc*v->nx ;
	ys = y0 + e1c*v->e1y + e2c*v->e2y + nc*v->ny;
	zs = z0 + e1c*v->e1z + e2c*v->e2z + nc*v->nz ;

	MRISsurfaceRASToVoxel(mris, mri, xs, ys, zs, &xv, &yv, &zv) ;
	MRIsampleVolume(mri, xv, yv, zv, &val) ;
	MRIsetVoxVal(mri_vol, xi, yi, zi, 0, val) ;
      }

  
  return(mri_vol) ;
}
