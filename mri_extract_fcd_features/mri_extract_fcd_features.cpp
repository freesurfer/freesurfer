/**
 * @brief program to extract features from a recon dir to use as input to an FCD classifier.
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
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "histo.h"
#include "cma.h"
#include "label.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;


static char sdir[STRLEN] = "" ;
static const char *white_name = "white" ;
static const char *pial_name = "pial" ;
static const char *vol_name = "norm.mgz" ;
static const char *ribbon_name = "ribbon.mgz" ;
static const char *aparc_name = "aparc+aseg.mgz" ;
//static char *annot_name = "aparc" ;
static const char *aseg_name = "aseg.mgz" ;
static const char *sphere_name = "sphere.d1.left_right";
static const char *cortex_label = "cortex" ;
static const char *flair_name = "FLAIR.masked.mgz" ;



static  MRI *MRIcomputeSurfaceDistanceProbabilities(MRI_SURFACE *mris,  MRI *mri_ribbon, MRI *mri, MRI *mri_aseg) ;
static  MRI *MRIcomputeSurfaceDistanceIntensities(MRI_SURFACE *mris, MRI *mri_ribbon, MRI *mri_aparc, MRI *mri, MRI *mri_aseg, int whalf) ;
static  MRI *MRIcomputeFlairRatio(MRI_SURFACE *mris, MRI *mri_ribbon, MRI *mri_aparc, MRI *mri, MRI *mri_flair, MRI *mri_aseg, int whalf) ;

static int whalf = 5 ;
static int navgs = 0 ;

int
main(int argc, char *argv[]) 
{
  char   fname[STRLEN], *cp ;
  int    nargs ;
  char   *subject, *out_fname, *hemi;
  const char *ohemi ;
  int    msec, minutes, seconds ;
  Timer start ;
  MRI          *mri, *mri_features, *mri_ribbon, *mri_aseg, *mri_aparc, *mri_flair ;
  MRI_SURFACE  *mris, *mris_contra ;
  LABEL        *cortex ;

  nargs = handleVersionOption(argc, argv, "mri_extract_fcd_features");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  subject = argv[1] ;
  hemi = argv[2] ;
  if (strcmp(hemi, "lh") == 0)
    ohemi = "rh" ;
  else
    ohemi = "lh" ;

  out_fname = argv[3] ;
  printf("reading data for subject %s and writing output to %s\n", subject, out_fname) ;

  if (!strlen(sdir))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, white_name) ;
  printf("reading %s\n", fname) ;
  mris  = MRISread(fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s\n", Progname, fname) ;
  MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
  

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, ohemi, white_name) ;
  printf("reading %s\n", fname) ;
  mris_contra  = MRISread(fname) ;
  if (mris_contra == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s\n", Progname, fname) ;
  MRISsaveVertexPositions(mris_contra, WHITE_VERTICES) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, subject, ribbon_name) ;
  printf("reading %s\n", fname) ;
  mri_ribbon  = MRIread(fname) ;
  if (mri_ribbon == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read ribbon from %s\n", Progname, fname) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, subject, aparc_name) ;
  printf("reading %s\n", fname) ;
  mri_aparc  = MRIread(fname) ;
  if (mri_aparc == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read ribbon from %s\n", Progname, fname) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, subject, aseg_name) ;
  printf("reading %s\n", fname) ;
  mri_aseg  = MRIread(fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read aseg from %s\n", Progname, fname) ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, pial_name) ;
  if (MRISreadPialCoordinates(mris, fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read pial coordinates from %s\n", Progname, fname) ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, sphere_name) ;
  if (MRISreadCanonicalCoordinates(mris, fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read left/right spherical coordinates from %s\n", Progname, fname) ;
  

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, ohemi, pial_name) ;
  if (MRISreadPialCoordinates(mris_contra, fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read pial coordinates from %s\n", Progname, fname) ;

  sprintf(fname, "%s/%s/label/%s.%s", sdir, subject, hemi, cortex_label) ;
  cortex = LabelRead(NULL, fname) ;
  if (cortex == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read cortical label from %s\n", Progname, fname) ;
  LabelRipRestOfSurface(cortex, mris) ;
  LabelFree(&cortex) ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, ohemi, sphere_name) ;
  if (MRISreadCanonicalCoordinates(mris_contra, fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read left/right spherical coordinates from %s\n", Progname, fname) ;
  

  sprintf(fname, "%s/%s/mri/%s", sdir, subject, vol_name) ; 
  printf("reading %s\n", fname) ;
  mri  = MRIread(fname) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s\n", Progname, fname) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, subject, flair_name) ; 
  printf("reading %s\n", fname) ;
  mri_flair  = MRIread(fname) ;
  if (mri_flair == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s\n", Progname, fname) ;

  if (0)
    mri_features = MRIcomputeSurfaceDistanceProbabilities(mris, mri_ribbon, mri, mri_aseg) ;
  else
  {
    MRI *mri_ohemi_features, *mri_ohemi_mapped_to_hemi_features, *mri_flair_features, *mri_ohemi_flair_features, *mri_ohemi_mapped_to_hemi_flair_features ;
    char fname[STRLEN], ext[STRLEN], fname_no_ext[STRLEN] ;

    mri_features = MRIcomputeSurfaceDistanceIntensities(mris, mri_ribbon, mri_aparc, mri, mri_aseg, whalf) ;
    mri_ohemi_features = MRIcomputeSurfaceDistanceIntensities(mris_contra, mri_ribbon, mri_aparc, mri, mri_aseg, whalf) ;
    mri_ohemi_mapped_to_hemi_features = MRISmapToSurface(mris_contra, mris, mri_ohemi_features, NULL) ; // map contra feature to this surface
    MRIsubtract(mri_features, mri_ohemi_mapped_to_hemi_features, mri_features) ;
//    MRIwrite(mri_ohemi_mapped_to_hemi_features, "test.mgz") ;

    mri_flair_features = MRIcomputeFlairRatio(mris, mri_ribbon, mri_aparc, mri, mri_flair, mri_aseg, whalf) ;
    mri_ohemi_flair_features = MRIcomputeFlairRatio(mris_contra, mri_ribbon, mri_aparc, mri, mri_flair, mri_aseg, whalf) ;
    mri_ohemi_mapped_to_hemi_flair_features = MRISmapToSurface(mris_contra, mris, mri_ohemi_flair_features, NULL) ; // map contra feature to this s
    MRIsubtract(mri_flair_features, mri_ohemi_mapped_to_hemi_flair_features, mri_flair_features) ;
    strcpy(fname, out_fname) ;
    FileNameExtension(fname, ext) ;
    FileNameRemoveExtension(fname, fname_no_ext) ;
    sprintf(fname, "%s.flair.%s", fname_no_ext, ext) ;
    if (Gdiag_no >= 0)
      printf("feature(%d) = %f\n", Gdiag_no, MRIgetVoxVal(mri_flair_features, Gdiag_no, 0, 0, 0)) ;
    printf("DISABLED: writing output to %s\n", fname) ;
//    MRIwrite(mri_flair_features, fname) ;
  }
 
  if (navgs > 0)
  {
    MRI *mri_tmp ;
    mri_tmp = MRISsmoothMRI(mris, mri_features, navgs, NULL, NULL);
    MRIfree(&mri_features) ;
    mri_features = mri_tmp ;
  }
  printf("writing output to %s\n", out_fname) ;
  if (Gdiag_no >= 0)
    printf("feature(%d) = %f\n", Gdiag_no, MRIgetVoxVal(mri_features, Gdiag_no, 0, 0, 0)) ;


  MRIwrite(mri_features, out_fname) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "feature extraction took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "SDIR"))
  {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  }
  else  if (!stricmp(option, "VNO"))
  {
    Gdiag_no = atoi(argv[2]) ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    nargs = 1 ;
  }
  else  if (!stricmp(option, "PIAL"))
  {
    pial_name = argv[2] ;
    printf("using pial surface '%s'\n", pial_name) ;
    nargs = 1 ;
  }
  else  if (!stricmp(option, "WHITE"))
  {
    white_name = argv[2] ;
    printf("using white surface '%s'\n", white_name) ;
    nargs = 1 ;
  }
  else  if (!stricmp(option, "VOL"))
  {
    vol_name = argv[2] ;
    printf("using intensity volume '%s'\n", vol_name) ;
    nargs = 1 ;
  }
  else switch (toupper(*option)) {
    case 'A':
      navgs = atof(argv[2]) ;
      nargs = 1 ;
      printf("using %d smoothing iterations after feature extraction\n", navgs) ;
      break ;
    case 'W':
      whalf = atoi(argv[2]) ;
      nargs = 1 ;
      printf("setting whalf = %d\n", whalf) ;
      break ;
    case 'V':
      vol_name = argv[2] ;
      nargs = 1 ;
      printf("using %s as vol name\n", vol_name) ;
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
usage_exit(int code) {
  printf("usage: %s [options] <subject> <hemi> <output file>\n", Progname) ;
  printf("\tsdir <subjects dir> - specify SUBJECTS_DIR on the command line instead of in env\n") ;
  exit(code) ;
}



static double max_white_dist = 5 ;
static double min_white_dist = -5 ;

static double max_pial_dist = 3 ;
static double min_pial_dist = -2 ;

static double white_bin_size = .5 ;
static double pial_bin_size = .5 ;

static int NBINS = 256 ;

#define SKIP_LABEL(l)   (IS_CEREBELLAR_GM(l) || IS_VENTRICLE(l) || l == Brain_Stem || IS_CC(l) || IS_WMSA(l) || IS_CAUDATE(l) || IS_AMYGDALA(l) || IS_HIPPO(l))

static int close_order = 3 ;

static MRI *
MRIcomputeSurfaceDistanceProbabilities(MRI_SURFACE *mris,  MRI *mri_ribbon, MRI *mri, MRI *mri_aseg) 
{
  MRI          *mri_features, *mri_binary, *mri_white_dist, *mri_pial_dist, *mri_mask ;
  int          nwhite_bins, npial_bins, x, y, z, label, i ;
  HISTOGRAM2D *hw, *hp ; 
  float        pdist, wdist, val ;
  double       wval, pval ;

  mri_features = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 2) ;
  MRIcopyHeader(mri, mri_features) ;

  mri_binary = MRIcopy(mri_ribbon, NULL) ;
  mri_binary = MRIbinarize(mri_ribbon, NULL, 1, 0, 1) ;
  mri_pial_dist = MRIdistanceTransform(mri_binary, NULL, 1, max_pial_dist+1, DTRANS_MODE_SIGNED,NULL);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_pial_dist, "pd.mgz") ;

  MRIclear(mri_binary) ;
  MRIcopyLabel(mri_ribbon, mri_binary, Left_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_ribbon, mri_binary, Right_Cerebral_White_Matter) ;
  MRIbinarize(mri_binary, mri_binary, 1, 0, 1) ;
  mri_white_dist = MRIdistanceTransform(mri_binary, NULL, 1, max_white_dist+1, DTRANS_MODE_SIGNED,NULL);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_white_dist, "wd.mgz") ;

  nwhite_bins = ceil((max_white_dist - min_white_dist) / white_bin_size)+1 ;
  npial_bins = ceil((max_pial_dist - min_pial_dist) / pial_bin_size)+1 ;
  hw = HISTO2Dalloc(NBINS, nwhite_bins) ;
  hp = HISTO2Dalloc(NBINS, npial_bins) ;

  HISTO2Dinit(hw, NBINS, nwhite_bins, 0, NBINS-1, min_white_dist, max_white_dist) ;
  HISTO2Dinit(hp, NBINS, npial_bins, 0, NBINS-1, min_pial_dist, max_pial_dist) ;

  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
	label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	if (IS_CEREBELLAR_GM(label) || IS_VENTRICLE(label) || label == Brain_Stem || IS_CC(label))
	  continue ;
	pdist = MRIgetVoxVal(mri_pial_dist, x, y, z, 0) ;
	wdist = MRIgetVoxVal(mri_pial_dist, x, y, z, 0) ;
	val  = MRIgetVoxVal(mri, x, y, z, 0) ;
	if (pdist >= min_pial_dist && pdist <= max_pial_dist)
	  HISTO2DaddSample(hp, val, pdist, 0, 0, 0, 0) ;
	if (wdist >= min_white_dist && wdist <= max_white_dist)
	  HISTO2DaddSample(hw, val, wdist, 0, 0, 0, 0) ;
      }
  HISTO2DmakePDF(hp, hp) ;
  HISTO2DmakePDF(hw, hw) ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
	label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	if (IS_CEREBELLAR_GM(label) || IS_VENTRICLE(label) || label == Brain_Stem || IS_CC(label))
	  continue ;
	pdist = MRIgetVoxVal(mri_pial_dist, x, y, z, 0) ;
	wdist = MRIgetVoxVal(mri_pial_dist, x, y, z, 0) ;
	val  = MRIgetVoxVal(mri, x, y, z, 0) ;
	wval = HISTO2DgetCount(hw, val, wdist);
	if (DZERO(wval) == 0)
	  MRIsetVoxVal(mri_features, x, y, z, 1, -log10(wval)) ;
	pval = HISTO2DgetCount(hp, val, pdist);
	if (DZERO(pval) == 0)
	  MRIsetVoxVal(mri_features, x, y, z, 0, -log10(pval)) ;
      }

  MRIclear(mri_binary) ;
  MRIbinarize(mri_ribbon, mri_binary, 1, 0, 1) ;
  mri_mask = MRIcopy(mri_binary, NULL) ;
  for (i = 0 ; i < close_order ; i++)
  {
    MRIdilate(mri_binary, mri_mask) ;
    MRIcopy(mri_mask, mri_binary) ;
  }
  
  for (i = 0 ; i < close_order ; i++)
  {
    MRIerode(mri_binary, mri_mask) ;
    MRIcopy(mri_mask, mri_binary) ;
  }

  MRIwrite(mri_mask, "m.mgz") ;
  MRImask(mri_features, mri_mask, mri_features, 0, 0) ;
  HISTO2Dfree(&hw) ;
  HISTO2Dfree(&hp) ;
  MRIfree(&mri_white_dist) ; MRIfree(&mri_pial_dist) ; MRIfree(&mri_binary) ;
  return(mri_features) ;
}


#define MAX_SEARCH_DIST  50

#define angle_threshold RADIANS(20.0)

static MRI *
MRIcomputeSurfaceDistanceIntensities(MRI_SURFACE *mris,  MRI *mri_ribbon, MRI *mri_aparc, MRI *mri, MRI *mri_aseg, int whalf) 
{
  MRI          *mri_features, *mri_binary, *mri_white_dist, *mri_pial_dist ;
  int          vno, ngm, outside_of_ribbon, label0, label, ohemi_label, xi, yi, zi, xk, yk, zk, x0, y0, z0, hemi_label, assignable ;
  double       xv, yv, zv, step_size, dist, thickness, wdist, pdist, snx, sny, snz, nx, ny, nz, xl, yl, zl, x, y, z, dot, angle ;
  VERTEX       *v ;

  mri_features = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 1) ;  // one samples inwards, one in ribbon, and one outside
  MRIcopyHeader(mri, mri_features) ;

  mri_binary = MRIcopy(mri_ribbon, NULL) ;
  mri_binary = MRIbinarize(mri_ribbon, NULL, 1, 0, 1) ;
  mri_pial_dist = MRIdistanceTransform(mri_binary, NULL, 1, max_pial_dist+1, DTRANS_MODE_SIGNED,NULL);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_pial_dist, "pd.mgz") ;

  MRIclear(mri_binary) ;
  MRIcopyLabel(mri_ribbon, mri_binary, Left_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_ribbon, mri_binary, Right_Cerebral_White_Matter) ;
  MRIbinarize(mri_binary, mri_binary, 1, 0, 1) ;
  mri_white_dist = MRIdistanceTransform(mri_binary, NULL, 1, max_white_dist+1, DTRANS_MODE_SIGNED,NULL);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_white_dist, "wd.mgz") ;

  if (mris->hemisphere == LEFT_HEMISPHERE)
  {
    ohemi_label = Right_Cerebral_Cortex ;
    hemi_label = Left_Cerebral_Cortex ;
  }
  else
  {
    hemi_label = Right_Cerebral_Cortex ;
    ohemi_label = Left_Cerebral_Cortex ;
  }

  step_size = mri->xsize/2 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;  // not cortex
    nx = v->pialx - v->whitex ; ny = v->pialy - v->whitey ; nz = v->pialz - v->whitez ;
    thickness = sqrt(nx*nx + ny*ny + nz*nz) ;
    if (FZERO(thickness))
      continue ;   // no  cortex here


    x = (v->pialx + v->whitex)/2 ; y = (v->pialy + v->whitey)/2 ; z = (v->pialz + v->whitez)/2 ;  // halfway between white and pial is x0
    MRISsurfaceRASToVoxelCached(mris, mri_aseg, x, y, z, &xl, &yl, &zl) ;
    x0 = nint(xl); y0 = nint(yl) ; z0 = nint(zl) ;
    label0 = MRIgetVoxVal(mri_aparc, x0, y0, z0,0) ;

    // compute surface normal in voxel coords
    MRISsurfaceRASToVoxelCached(mris, mri_aseg, x+v->nx, y+v->ny, z+v->nz, &snx, &sny, &snz) ;
    snx -= xl ; sny -= yl ; snz -= zl ;

    for (ngm = 0, xk = -whalf ; xk <= whalf ; xk++)
    {
      xi = mri_aseg->xi[x0+xk] ;
      for (yk = -whalf ; yk <= whalf ; yk++)
      {
	yi = mri_aseg->yi[y0+yk] ;
	for (zk = -whalf ; zk <= whalf ; zk++)
	{
	  zi = mri_aseg->zi[z0+zk] ;
	  label = MRIgetVoxVal(mri_aseg, xi, yi, zi,0) ;
	  if (xi == Gx && yi == Gy && zi == Gz)
	    DiagBreak() ;
	  if (label != hemi_label)
	    continue ;
	  label = MRIgetVoxVal(mri_aparc, xi, yi, zi,0) ;
	  if (label && label != label0)  // if  outside the ribbon it won't be assigned to a parcel
	    continue ;  // constrain it to be in the same cortical parcel

	  // search along vector connecting x0 to this point to make sure it is we don't perforate wm or leave and re-enter cortex
	  nx = xi-x0 ; ny = yi-y0 ; nz = zi-z0 ;
	  thickness = sqrt(nx*nx + ny*ny + nz*nz) ;
	  assignable = 1 ;  // assume this point should be counted
	  if (thickness > 0)
	  {
	    nx /= thickness ; ny /= thickness ; nz /= thickness ;
	    dot = nx*snx + ny*sny + nz*snz ; angle = acos(dot) ;
	    if (fabs(angle) > angle_threshold)
	      assignable = 0 ;
	    outside_of_ribbon = 0 ;
	    for (dist = 0 ; assignable && dist <= thickness ; dist += step_size) 
	    {
	      xv = x0+nx*dist ;  yv = y0+ny*dist ;  zv = z0+nz*dist ; 
	      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
		DiagBreak() ;
	      MRIsampleVolume(mri_pial_dist, xv, yv, zv, &pdist) ;
	      MRIsampleVolume(mri_white_dist, xv, yv, zv, &wdist) ;
	      label = MRIgetVoxVal(mri_aseg, xi, yi, zi,0) ;
	      if (SKIP_LABEL(label) || label == ohemi_label)
		assignable = 0 ;
	      if (wdist < 0)  // entered wm - not assignable
		assignable = 0 ;
	      else
	      {
		if (pdist > 0)  // outside pial surface
		  outside_of_ribbon = 1 ;
		else
		{
		  if (outside_of_ribbon) // left ribbon and reentered
		    assignable = 0 ;
		}
	      }
	    }
	  }  // close of thickness > 0
	  if (assignable)
	    ngm++ ;
	  else
	    DiagBreak() ;
	}
      }
    }
    
    MRIsetVoxVal(mri_features, vno, 0, 0, 0, ngm) ;
  }

  MRIfree(&mri_white_dist) ; MRIfree(&mri_pial_dist) ; MRIfree(&mri_binary) ;
  return(mri_features) ;
}

#define SAMPLE_DIST 3.0
static  MRI *
MRIcomputeFlairRatio(MRI_SURFACE *mris, MRI *mri_ribbon, MRI *mri_aparc, MRI *mri, MRI *mri_flair, MRI *mri_aseg, int whalf) 
{
  MRI    *mri_features, *mri_dist_white, *mri_tmp, *mri_filled ;
  int     xv, yv, zv, vno, label ;
  double  xs, ys, zs, xf, yf, zf, flair, T1, dist, xd, yd, zd, dx, dy, dz, nx, ny, nz, x0, y0, z0, len, dot, theta ;
  VERTEX  *v ;

  MRISsaveVertexPositions(mris, TMP2_VERTICES) ;
  MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
  mri_tmp = MRISfillInterior(mris, mri->xsize/2, NULL) ;
  mri_filled = MRIextractRegionAndPad(mri_tmp, NULL, NULL, nint(30/mri->xsize)) ; MRIfree(&mri_tmp) ;
  mri_dist_white = MRIcloneDifferentType(mri_filled, MRI_FLOAT) ;
  MRIdistanceTransform(mri_filled, mri_dist_white, 1, 100, DTRANS_MODE_SIGNED, NULL) ;
  MRISrestoreVertexPositions(mris, TMP2_VERTICES) ;
  MRIfree(&mri_filled) ;

  mri_features = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 1) ;  // one samples inwards, one in ribbon, and one outside
  MRIcopyHeader(mri, mri_features) ;

  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    xs = v->whitex - v->nx * SAMPLE_DIST ; ys = v->whitey - v->ny * SAMPLE_DIST ; zs = v->whitez - v->nz * SAMPLE_DIST ;
    MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xf, &yf, &zf) ;
    xv = nint(xf) ; yv = nint(yf) ; zv = nint(zf) ;
    label = MRIgetVoxVal(mri_aseg, xv, yv, zv, 0) ;
    if (IS_WHITE_MATTER(label) == 0)
      continue ;

    MRISsurfaceRASToVoxelCached(mris, mri_dist_white, v->whitex, v->whitey, v->whitez, &x0, &y0, &z0) ;
    MRISsurfaceRASToVoxelCached(mris, mri_dist_white, xs, ys, zs, &xd, &yd, &zd) ;
    MRIsampleVolumeGradient(mri_dist_white, xd, yd, zd, &dx, &dy, &dz) ; dx *= -1 ; dy *= -1 ; dz *= -1 ;
    MRIsampleVolume(mri_dist_white, xd, yd, zd, &dist) ;

    nx = xd-x0 ;  ny = yd-y0 ; nz = zd-z0 ; 
    len = sqrt(nx*nx + ny*ny + nz*nz) ;  nx /= len ; ny /= len ; nz /= len ;
    dot = nx*dx + ny*dy + nz*dz ;
    theta = acos(dot) * 360 / (2*M_PI) ;
    if (theta < 0 || theta > 90)
    {
      if (vno == Gdiag_no)
	printf("v %d: angle %2.3f too large\n", vno, theta) ;
//      MRIsetVoxVal(mri_features, vno, 0, 0, 0, 0) ;
      continue ;
    }

    if (dist > -SAMPLE_DIST/2)   // white matter is too thin here
    {
      if (vno == Gdiag_no)
	printf("vno %d: white dist %2.2f not close enough to sample dist %2.2f\n", vno,dist, SAMPLE_DIST) ;
//      MRIsetVoxVal(mri_features, vno, 0, 0, 0, 0) ;
      continue ;
    }

    MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xf, &yf, &zf) ;
    MRIsampleVolume(mri_flair, xf, yv, zf, &flair) ;
    MRIsampleVolume(mri, xf, yv, zf, &T1) ; T1 = 1 ;
    if (FZERO(T1))
      continue ;

    v->marked = 1 ;
    MRIsetVoxVal(mri_features, vno, 0, 0, 0, flair/T1) ;
  }

  MRISreadFromVolume(mri_features, mris, VERTEX_VALS) ;
  MRISsoapBubbleVals(mris, 100) ;
  MRISwriteIntoVolume(mris, mri_features, VERTEX_VALS) ;
  
  MRIfree(&mri_dist_white) ;
  return(mri_features) ;
}

