/**
 * @file  mri_auto_fill.c
 * @brief usage: ./mri_auto_fill filledvol T1vol m3dxform templvol outvol
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:13 $
 *    $Revision: 1.14 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "const.h"
#include "transform.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *MRIcombineHemispheres(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                                  MRI *mri_lh_prob, MRI *mri_rh_prob) ;
static int MRIfillVolume(MRI *mri) ;
static MRI *MRIcheckHemisphereOverlap(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                                      MRI *mri_lh_prob, MRI *mri_rh_prob) ;
static MRI *MRIthresholdFilled(MRI *mri_src, MRI *mri_T1, MRI *mri_mask,
                               MRI *mri_inv_T1, MRI *mri_inv_T1_std,
                               MRI *mri_dst,
                               float threshold, float nsigma, int out_label) ;
static MRI *MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_T1, float thresh,
                             int out_label, MRI *mri_dst);
#if 0
static MRI *MRIfillVentricles(MRI *mri_src, MRI *mri_T1, MRI *mri_mask,
                              MRI *mri_inv_T1, MRI *mri_inv_T1_std,
                              MRI *mri_dst,
                              int out_label) ;
#endif

char *Progname ;

static void usage_exit(int code) ;

#define LH_FILLED_VOLUME 4
#define RH_FILLED_VOLUME 6

static float pct = 95.0f ;
static float nsigma = 0.1f ;
static int fix = 0 ;
static int dilate = 0 ;
int old_flag = 0; // use M3Dmorph instead of GCAmorph

int
main(int argc, char *argv[])
{
  char   **av, path[STRLEN], ventricle_fname[500] ;
  int    ac, nargs ;
  MRI    *mri_filled, *mri_lh_template, *mri_lh_inverse_template, *mri_lh,
  *mri_rh, *mri_lv, *mri_rv, *mri_inv_lv, *mri_inv_rv,
  *mri_rh_template, *mri_rh_inverse_template, *mri_T1, *mri_inv_T1,
  *mri_inv_T1_std, *mri_tmp ;
  char   *filled_fname, *template_fname, *out_fname, *xform_fname, fname[100],
  *T1_fname ;
  M3D    *m3d ;
  GCAM   *gcam;
  struct  timeb start ;
  int     msec ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mri_auto_fill.c,v 1.14 2011/03/02 00:04:13 nicks Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  TimerStart(&start) ;

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

  if (argc < 5)
    usage_exit(1) ;

  filled_fname = argv[1] ;
  T1_fname = argv[2] ;
  xform_fname = argv[3] ;
  template_fname = argv[4] ;
  out_fname = argv[5] ;

  /////////////////////////////////////////////////////////////////////
  fprintf(stderr, "reading transform %s...", xform_fname) ;
  if (old_flag)
  {
    m3d = MRI3DreadSmall(xform_fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, xform_fname) ;
  }
  else
  {
    gcam = GCAMread(xform_fname);
    if (!gcam)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, xform_fname) ;
  }
  fprintf(stderr, "done.\n") ;
  /////////////////////////////////////////////////////////////////////

  /* left hemisphere filled volume */
  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, LH_FILLED_VOLUME);
  else
    strcpy(fname, template_fname) ;
  fprintf(stderr, "reading lh filled template %s...\n", fname) ;
  mri_lh_template = MRIread(fname) ;
  if (!mri_lh_template)
    ErrorExit(ERROR_NOFILE, "%s: could not read lh template volume %s.\n",
              Progname, template_fname) ;
  fprintf(stderr, "applying inverse transform...\n") ;
  mri_lh_inverse_template = MRIapplyInverse3DMorph(mri_lh_template,m3d,NULL);
  MRIfree(&mri_lh_template) ;

  /* right hemisphere filled volume */
  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, RH_FILLED_VOLUME);
  else
    strcpy(fname, template_fname) ;
  fprintf(stderr, "reading rh filled template %s...\n", fname) ;
  mri_rh_template = MRIread(fname) ;
  if (!mri_rh_template)
    ErrorExit(ERROR_NOFILE, "%s: could not read rh template volume %s.\n",
              Progname, template_fname) ;
  /////////////////////////////////////////////////////////////////////////
  fprintf(stderr, "applying inverse transform...\n") ;
  if (old_flag)
    mri_rh_inverse_template = MRIapplyInverse3DMorph(mri_rh_template,m3d,NULL);
  else
  {
    mri_rh_inverse_template = MRIclone(mri_rh_template, NULL);
    mri_rh_inverse_template =
      GCAMmorphFromAtlas(mri_rh_template, gcam, mri_rh_inverse_template,
                         SAMPLE_TRILINEAR);
  }
  MRIfree(&mri_rh_template) ;
  /////////////////////////////////////////////////////////////////////////
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    fprintf(stderr, "writing inverse image to lh and rh_inverse.mgh\n") ;
    MRIwrite(mri_lh_inverse_template, "lh_inverse.mgh@mgh") ;
    MRIwrite(mri_rh_inverse_template, "rh_inverse.mgh@mgh") ;
  }


  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, 0);
  mri_tmp = MRIread(fname) ;
  if (!mri_tmp)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 template volume %s.\n",
              Progname, fname) ;
  ///////////////////////////////////////////////////////////////////////////
  if (old_flag)
    mri_inv_T1 = MRIapplyInverse3DMorph(mri_tmp, m3d, NULL);
  else
  {
    mri_inv_T1 = MRIclone(mri_tmp, NULL);
    mri_inv_T1 = GCAMmorphFromAtlas(mri_tmp, gcam, mri_inv_T1,SAMPLE_TRILINEAR);
  }
  MRIfree(&mri_tmp) ;
  //////////////////////////////////////////////////////////////////////////
  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, 1);
  mri_tmp = MRIread(fname) ;
  if (!mri_tmp)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 template volume %s.\n",
              Progname, fname) ;
  //////////////////////////////////////////////////////////////////////////
  if (old_flag)
    mri_inv_T1_std = MRIapplyInverse3DMorph(mri_tmp, m3d, NULL);
  else
  {
    mri_inv_T1_std = MRIclone(mri_tmp, NULL);
    mri_inv_T1_std = GCAMmorphFromAtlas(mri_tmp, gcam, mri_inv_T1_std,SAMPLE_TRILINEAR);
  }
  MRIfree(&mri_tmp) ;
  //////////////////////////////////////////////////////////////////////////
  mri_T1 = MRIread(T1_fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume %s.\n",
              Progname, T1_fname) ;
  //////////////////////////
  if (old_flag)
    MRI3DmorphFree(&m3d) ;
  else
    GCAMfree(&gcam);
  //////////////////////////

  mri_filled = MRIread(filled_fname) ;
  if (!mri_filled)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s.\n",
              Progname, filled_fname) ;


  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "using inverse image to fill lh volume...\n") ;
  mri_lh = MRIthresholdFilled(mri_filled, mri_T1, mri_lh_inverse_template,
                              mri_inv_T1, mri_inv_T1_std, NULL, pct,nsigma,
                              MRI_LEFT_HEMISPHERE);
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "using inverse image to fill rh volume...\n") ;
  mri_rh = MRIthresholdFilled(mri_filled, mri_T1, mri_rh_inverse_template,
                              mri_inv_T1, mri_inv_T1_std, NULL, pct,nsigma,
                              MRI_LEFT_HEMISPHERE);
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,
            "applying inverse morph to left ventricle representation\n");
  FileNamePath(template_fname, path) ;
  sprintf(ventricle_fname, "%s/left_ventricle.mgh#0", path) ;
  mri_lv = MRIread(ventricle_fname) ;
  if (!mri_lv)
    ErrorExit(ERROR_NOFILE,"%s: could not read %s",Progname,ventricle_fname);


  fprintf(stderr, "reading transform %s...", xform_fname) ;
  //////////////////////////////////////////////////////////////////////////
  if (old_flag)
  {
    m3d = MRI3DreadSmall(xform_fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, xform_fname) ;
    mri_inv_lv = MRIapplyInverse3DMorph(mri_lv,m3d,NULL);
    MRI3DmorphFree(&m3d) ;
  }
  else
  {
    gcam = GCAMread(xform_fname);
    if (!gcam)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, xform_fname) ;
    mri_inv_lv = MRIclone(mri_lv, NULL);
    mri_inv_lv = GCAMmorphFromAtlas(mri_lv, gcam, mri_inv_lv,SAMPLE_TRILINEAR);
    GCAMfree(&gcam);
  }
  fprintf(stderr, "done.\n") ;
  ///////////////////////////////////////////////////////////////////////////

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MRIwrite(mri_inv_lv, "inv_lv.mgh@mgh") ;
  MRIfree(&mri_lv) ;

  mri_lv = MRIfillVentricle(mri_inv_lv,mri_T1,100,MRI_LEFT_HEMISPHERE,NULL);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MRIwrite(mri_lv, "lv.mgh@mgh") ;
  MRIfree(&mri_inv_lv) ;
  MRIunion(mri_lv, mri_lh, mri_lh) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MRIwrite(mri_lh, "lv_union.mgh@mgh") ;
  MRIfree(&mri_lv) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,
            "applying inverse morph to right ventricle representation\n");

  sprintf(ventricle_fname, "%s/right_ventricle.mgh#0", path) ;
  mri_rv = MRIread(ventricle_fname) ;
  if (!mri_rv)
    ErrorExit(ERROR_NOFILE,"%s: could not read %s",Progname,ventricle_fname);

  /////////////////////////////////////////////////////////
  fprintf(stderr, "reading transform %s...", xform_fname) ;
  if (old_flag)
  {
    m3d = MRI3DreadSmall(xform_fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, xform_fname) ;
    fprintf(stderr, "done.\n") ;
    mri_inv_rv = MRIapplyInverse3DMorph(mri_rv,m3d,NULL);
    MRI3DmorphFree(&m3d) ;
  }
  else
  {
    m3d = MRI3DreadSmall(xform_fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, xform_fname) ;
    mri_inv_rv = MRIclone(mri_rv, NULL);
    mri_inv_rv = GCAMmorphFromAtlas(mri_rv, gcam, mri_inv_rv,SAMPLE_TRILINEAR);
    GCAMfree(&gcam);
  }
  MRIfree(&mri_rv) ;
  /////////////////////////////////////////////////////////
  mri_rv = MRIfillVentricle(mri_inv_rv,mri_T1,100,MRI_RIGHT_HEMISPHERE,NULL);
  MRIfree(&mri_inv_rv) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MRIwrite(mri_rv, "rv.mgh@mgh") ;
  MRIunion(mri_rv, mri_rh, mri_rh) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MRIwrite(mri_rh, "rv_union.mgh@mgh") ;
  MRIfree(&mri_rv) ;

  MRIfree(&mri_filled) ;

  /*
     or the labeled left (mri_lh) and right (mri_rh) volumes
     together. In cases where both hemispheres label is on,
     choose the one with the higher probability as defined
     by mri_lh_inverse_template and mri_rh_inverse_template
  */
  mri_filled = MRIcombineHemispheres(mri_lh, mri_rh, NULL,
                                     mri_lh_inverse_template,
                                     mri_rh_inverse_template) ;
  fprintf(stderr, "done.\n") ;

  /*
    change the label of points that are on, but are significantly more
    likely to belong to the other hemisphere (fixes failure of cutting
    planes).
  */
  if (fix)
    MRIcheckHemisphereOverlap(mri_lh, mri_rh, mri_filled,
                              mri_lh_inverse_template,mri_rh_inverse_template);

  /* make sure there are no interior holes or exterior islands */
  MRIfillVolume(mri_filled) ;
  MRIfree(&mri_lh) ;
  MRIfree(&mri_rh) ;
  MRIfree(&mri_lh_inverse_template) ;
  MRIfree(&mri_rh_inverse_template) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing filled volume to %s...\n", out_fname) ;
  if (MRIwrite(mri_filled, out_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not write output volume to %s.\n",
              Progname, out_fname) ;

  msec = TimerStop(&start) ;
  fprintf(stderr, "auto filling took %2.2f minutes\n",
          (float)msec/(1000.0f*60.0f));
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
  if (!stricmp(option, "dilate"))
  {
    dilate = atoi(argv[2]) ;
    fprintf(stderr, "dilating ventricles %d times before filling\n",dilate);
    nargs = 1 ;
  }
  if (!strcmp(option, "old"))
  {
    fprintf(stderr, "using old m3d morph routine\n");
    old_flag = 1;
  }
  else switch (toupper(*option))
    {
    case 'N':
      nsigma = atof(argv[2]) ;
      fprintf(stderr, "using nsigma = %2.2f\n", nsigma) ;
      nargs = 1 ;
      break ;
    case 'F':
      fix = !fix ;
      fprintf(stderr, "%scorrecting hemisphere overlap...\n",
              fix ? "not " : "") ;
      break ;
    case 'T':
    case 'P':
      pct = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using threshold = %2.1f%%\n", pct) ;
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
  printf("usage: %s <filled volume> <T1 volume> <m3d transform> "
         "<template volume> <output volume>\n",
         Progname) ;
  exit(code) ;
}
static MRI *
MRIcombineHemispheres(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                      MRI *mri_lh_prob, MRI *mri_rh_prob)
{
  int       x, y, z, width, height, depth ;
  BUFTYPE   *plh_prob, *prh_prob, *plh, *prh, *pdst,
  lh_prob, rh_prob, lh_label, rh_label, out_label ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_lh, NULL) ;

  width = mri_dst->width ;
  height = mri_dst->height ;
  depth = mri_dst->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      plh = &MRIvox(mri_lh, 0, y, z) ;
      prh = &MRIvox(mri_rh, 0, y, z) ;
      plh_prob = &MRIvox(mri_lh_prob, 0, y, z) ;
      prh_prob = &MRIvox(mri_rh_prob, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        lh_label = *plh++ ;
        rh_label = *prh++ ;
        lh_prob =  *plh_prob++ ;
        rh_prob =  *prh_prob++ ;
        if (lh_label && !rh_label)
          out_label = lh_label ;
        else if (rh_label && !lh_label)
          out_label = rh_label ;
        else if (rh_label && lh_label)
        {
          if (rh_prob > lh_prob)
            out_label = rh_label ;
          else
            out_label = lh_label ;
        }
        else
          out_label = 0 ;
        *pdst++ = out_label ;
      }
    }
  }

  return(mri_dst) ;
}
/* Talairach seed points for white matter in left and right hemispheres */
#define SEED_SEARCH_SIZE            9
#define MIN_NEIGHBORS               3
#define CORPUS_CALLOSUM_TAL_X       0.0
#define CORPUS_CALLOSUM_TAL_Y       0.0
#define CORPUS_CALLOSUM_TAL_Z       27.0

static int neighbors_on(MRI *mri, int x0, int y0, int z0, int label) ;

static int
MRIfillVolume(MRI *mri)
{
  Real     wm_rh_tal_x = 29.0 ;
  Real     wm_rh_tal_y = -12.0 ;
  Real     wm_rh_tal_z = 28.0 ;
  Real     wm_lh_tal_x = -29.0 ;
  Real     wm_lh_tal_y = -12.0 ;
  Real     wm_lh_tal_z = 28.0 ;
  int      wm_lh_x, wm_lh_y, wm_lh_z, wm_rh_x, wm_rh_y, wm_rh_z,xi,yi,zi,
  x, y, z, xnew, ynew,znew;
  Real     xr, yr, zr, dist, min_dist, xd, yd, zd ;
  MRI      *mri_fg, *mri_bg ;
#if 0
  BUFTYPE  val ;
  int      xlim0, xlim1, ylim0, ylim1, zlim0, zlim1 ;

  /* find bounding box for valid data */
  zlim0 = mri->depth ;
  ylim0 = mri->height ;
  xlim0 = mri->width ;
  zlim1 = ylim1 = xlim1 = 0;
  for (z = 0 ; z < mri->depth ; z++)
  {
    for (y = 0 ; y < mri->height; y++)
    {
      for (x = 0 ; x < mri->width ; x++)
      {
        val = MRIvox(mri, x, y, z) ;
        if (val == MRI_RIGHT_HEMISPHERE || val == MRI_LEFT_HEMISPHERE) ;
        {
          if (z<zlim0) zlim0=z;
          if (z>zlim1) zlim1=z;
          if (y<ylim0) ylim0=y;
          if (y>ylim1) ylim1=y;
          if (x<xlim0) xlim0=x;
          if (x>xlim1) xlim1=x;
        }
      }
    }
  }
#endif

  /* find white matter seed point for the left hemisphere */
  MRItalairachToVoxel(mri,CORPUS_CALLOSUM_TAL_X+(SEED_SEARCH_SIZE+1),
                      CORPUS_CALLOSUM_TAL_Y, CORPUS_CALLOSUM_TAL_Z,
                      &xr, &yr, &zr);
  MRItalairachToVoxel(mri, wm_lh_tal_x,wm_lh_tal_y,wm_lh_tal_z,&xr,&yr,&zr);
  wm_lh_x = nint(xr) ;
  wm_lh_y = nint(yr) ;
  wm_lh_z = nint(zr) ;
  if ((MRIvox(mri, wm_lh_x, wm_lh_y, wm_lh_z) != MRI_LEFT_HEMISPHERE) ||
      (neighbors_on(mri, wm_lh_x, wm_lh_y, wm_lh_z, MRI_LEFT_HEMISPHERE)
       < MIN_NEIGHBORS))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for lh wm seed...\n") ;
    for (z = wm_lh_z-SEED_SEARCH_SIZE ; z <= wm_lh_z+SEED_SEARCH_SIZE ; z++)
    {
      zi = mri->zi[z] ;
      for (y = wm_lh_y-SEED_SEARCH_SIZE ;
           y <= wm_lh_y+SEED_SEARCH_SIZE ;
           y++)
      {
        yi = mri->yi[y] ;
        for (x = wm_lh_x-SEED_SEARCH_SIZE;
             x <= wm_lh_x+SEED_SEARCH_SIZE ;
             x++)
        {
          xi = mri->xi[x] ;
          if ((MRIvox(mri, xi, yi, zi) == MRI_LEFT_HEMISPHERE) &&
              neighbors_on(mri, xi, yi, zi, MRI_LEFT_HEMISPHERE) >=
              MIN_NEIGHBORS)
          {
            xd = (xi - wm_lh_x) ;
            yd = (yi - wm_lh_y) ;
            zd = (zi - wm_lh_z) ;
            dist = xd*xd + yd*yd + zd*zd ;
            if (dist < min_dist)
            {
              xnew = xi ;
              ynew = yi ;
              znew = zi ;
              min_dist = dist ;
            }
          }
        }
      }
    }
    wm_lh_x = xnew ;
    wm_lh_y = ynew ;
    wm_lh_z = znew ;
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "lh seed point at (%d, %d, %d): %d neighbors on.\n",
            wm_lh_x, wm_lh_y, wm_lh_z,
            neighbors_on(mri, wm_lh_x, wm_lh_y, wm_lh_z,MRI_LEFT_HEMISPHERE));

  /* find white matter seed point for the right hemisphere */
  MRItalairachToVoxel(mri,CORPUS_CALLOSUM_TAL_X+(SEED_SEARCH_SIZE+1),
                      CORPUS_CALLOSUM_TAL_Y, CORPUS_CALLOSUM_TAL_Z,
                      &xr, &yr, &zr);
  MRItalairachToVoxel(mri, wm_rh_tal_x,wm_rh_tal_y,wm_rh_tal_z,&xr,&yr,&zr);
  wm_rh_x = nint(xr) ;
  wm_rh_y = nint(yr) ;
  wm_rh_z = nint(zr) ;
  if ((MRIvox(mri, wm_rh_x, wm_rh_y, wm_rh_z) == MRI_RIGHT_HEMISPHERE) ||
      (neighbors_on(mri, wm_rh_x, wm_rh_y, wm_rh_z, MRI_RIGHT_HEMISPHERE)
       < MIN_NEIGHBORS))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for rh wm seed...\n") ;
    for (z = wm_rh_z-SEED_SEARCH_SIZE ;
         z <= wm_rh_z+SEED_SEARCH_SIZE ;
         z++)
    {
      zi = mri->zi[z] ;
      for (y = wm_rh_y-SEED_SEARCH_SIZE ;
           y <= wm_rh_y+SEED_SEARCH_SIZE ;
           y++)
      {
        yi = mri->yi[y] ;
        for (x = wm_rh_x-SEED_SEARCH_SIZE ;
             x <= wm_rh_x+SEED_SEARCH_SIZE;
             x++)
        {
          xi = mri->xi[x] ;
          if ((MRIvox(mri, xi, yi, zi) == MRI_RIGHT_HEMISPHERE) &&
              (neighbors_on(mri, xi, yi, zi, MRI_RIGHT_HEMISPHERE)
               >= MIN_NEIGHBORS))
          {
            xd = (xi - wm_rh_x) ;
            yd = (yi - wm_rh_y) ;
            zd = (zi - wm_rh_z) ;
            dist = xd*xd + yd*yd + zd*zd ;
            if (dist < min_dist)
            {
              xnew = xi ;
              ynew = yi ;
              znew = zi ;
              min_dist = dist ;
            }
          }
        }
      }
    }
    wm_rh_x = xnew ;
    wm_rh_y = ynew ;
    wm_rh_z = znew ;
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "rh seed point at (%d, %d, %d): %d neighbors on.\n",
            wm_rh_x, wm_rh_y, wm_rh_z,
            neighbors_on(mri, wm_rh_x, wm_rh_y,wm_rh_z,MRI_RIGHT_HEMISPHERE));

  /* initialize the fill with the detected seed points */
  mri_fg = MRIclone(mri, NULL) ;
  MRIvox(mri_fg, wm_rh_x, wm_rh_y, wm_rh_z) = MRI_RIGHT_HEMISPHERE ;
  MRIvox(mri_fg, wm_lh_x, wm_lh_y, wm_lh_z) = MRI_LEFT_HEMISPHERE ;
  fprintf(stderr, "filling left hemisphere...\n") ;
  MRIgrowLabel(mri, mri_fg, MRI_LEFT_HEMISPHERE, MRI_LEFT_HEMISPHERE) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_fg, "lh.mgh@mgh") ;
  fprintf(stderr, "filling right hemisphere...\n") ;
  MRIgrowLabel(mri, mri_fg, MRI_RIGHT_HEMISPHERE, MRI_RIGHT_HEMISPHERE) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_fg, "rh.mgh@mgh") ;

  mri_bg = MRIclone(mri, NULL) ;
  MRIvox(mri_bg, 0, 0, 0) = 1 ;
  fprintf(stderr, "filling background...\n") ;
  MRIgrowLabel(mri, mri_bg, 0, 1) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_bg, "bg.mgh@mgh") ;
  MRIturnOnFG(mri, mri_fg, mri_bg) ;
  MRIturnOffBG(mri, mri_bg) ;

  MRIfree(&mri_fg) ;
  MRIfree(&mri_bg) ;
  return(NO_ERROR) ;
}

static int
neighbors_on(MRI *mri, int x0, int y0, int z0, int label)
{
  int   nbrs = 0 ;

  if (MRIvox(mri,x0-1,y0,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0+1,y0,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0+1,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0-1,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0+1) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0-1) == label)
    nbrs++ ;
  return(nbrs) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static MRI *
MRIcheckHemisphereOverlap(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                          MRI *mri_lh_prob, MRI *mri_rh_prob)
{
  int       x, y, z, width, height, depth ;
  BUFTYPE   *plh_prob, *prh_prob, *plh, *prh, *pdst, *psrc,
  lh_prob, rh_prob, lh_label, rh_label, label ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_lh, NULL) ;

  width = mri_dst->width ;
  height = mri_dst->height ;
  depth = mri_dst->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      plh = &MRIvox(mri_lh, 0, y, z) ;
      prh = &MRIvox(mri_rh, 0, y, z) ;
      plh_prob = &MRIvox(mri_lh_prob, 0, y, z) ;
      prh_prob = &MRIvox(mri_rh_prob, 0, y, z) ;
      psrc = &MRIvox(mri_dst, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        if (x == 131 && y == 135 && z == 76)
          DiagBreak() ;  /* marked as 255, should be 127 */
        if (x == 125 && y == 148 && z == 100)
          DiagBreak() ;  /* marked as 255, should be 127 */

        label = *psrc++ ;
        lh_label = *plh++ ;
        rh_label = *prh++ ;
        lh_prob =  *plh_prob++ ;
        rh_prob =  *prh_prob++ ;
        if (label)  /* white matter of one hemi or the other */
        {
#if 0
          if (lh_prob < rh_prob)
            label = MRI_RIGHT_HEMISPHERE ;
          else
            label = MRI_LEFT_HEMISPHERE ;
#else
          if (lh_prob < 10 && rh_prob > 2*lh_prob)
            label = MRI_RIGHT_HEMISPHERE ;
          else if ((label == rh_label) && rh_prob < 10 && lh_prob > 2*rh_prob)
            label = MRI_LEFT_HEMISPHERE ;
#endif
        }
        *pdst++ = label ;
      }
    }
  }
  return(mri_dst) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static MRI *
MRIthresholdFilled(MRI *mri_src, MRI *mri_T1, MRI *mri_mask, MRI *mri_inv_T1,
                   MRI *mri_inv_T1_std, MRI *mri_dst, float threshold,
                   float nsigma, int out_label)
{
  BUFTYPE   *pmask, *pdst, *psrc, out_val, mask, in_val,
  T1_val, T1_inv_val, sigma, *pT1, *pinvT1, *psigma ;
  int       width, height, depth, x, y, z, nchanged, noff, non ;
  float     nvox, sdist ;

  if (mri_mask->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "MRI3Dthreshold: mask must be MRI_FLOAT")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  /* now apply the inverse morph to build an average wm representation
     of the input volume
     */

  noff = non = nchanged = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pmask = &MRIvox(mri_mask, 0, y, z) ;
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      pT1 = &MRIvox(mri_T1, 0, y, z) ;
      pinvT1 = &MRIvox(mri_inv_T1, 0, y, z) ;
      psigma = &MRIvox(mri_inv_T1_std, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        T1_val = *pT1++ ;
        T1_inv_val = *pinvT1++ ;
        sigma = *psigma++ ;
        if (!FZERO(sigma))
          sdist = (float)abs(T1_val - T1_inv_val)/sigma ;
        else
          sdist = 0 ;
        if (x == 125 && y == 148 && z == 100)
          DiagBreak() ;  /* marked as 255, should be 127 */
        out_val = 0 ;
        mask = *pmask++ ;   /* value from inverse morphed volume */
        in_val = *psrc++ ;
        if (sdist > nsigma)  /* intensities are too different - don't change */
          out_val = in_val ;
        else
        {
          if (mask < 100-threshold)        /* probably off */
            out_val = 0 ;
          else  if (mask > threshold)      /* probably on */
            out_val = out_label ;
          else              /* not sure, use original fill val */
            out_val = in_val ;
        }
        if (out_val != in_val)
        {
          if (out_val)
            non++ ;
          else
            noff++ ;
          nchanged++ ;
        }
        *pdst++ = out_val ;
      }
    }
  }

  nvox = (float)(width * height * depth) ;
  fprintf(stderr, "%8d of %8d voxels changed - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)nchanged/nvox) ;
  fprintf(stderr, "%8d of %8d voxels off     - %2.1f%%.\n",
          noff, (int)nvox,   100.0f*(float)noff/nvox) ;
  fprintf(stderr, "%8d of %8d voxels on      - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)non/nvox) ;
  return(mri_dst) ;
}
#if 1
MRI *
MRIfillVentricle(MRI *mri_inv_ventricle, MRI *mri_T1, float thresh,
                 int out_label, MRI *mri_dst)
{
  BUFTYPE   *pdst, *pinv_ventricle, out_val, T1_val, inv_ventricle_val, *pT1 ;
  int       width, height, depth, x, y, z,
  ventricle_voxels, dno ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_T1, NULL) ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;
  /* now apply the inverse morph to build an average wm representation
     of the input volume
     */

  for (dno = dilate ; dno > 0 ; dno--)
    MRIdilate(mri_inv_ventricle, mri_inv_ventricle) ;

  ventricle_voxels = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      pT1 = &MRIvox(mri_T1, 0, y, z) ;
      pinv_ventricle = &MRIvox(mri_inv_ventricle, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        T1_val = *pT1++ ;
        inv_ventricle_val = *pinv_ventricle++ ;
        out_val = 0 ;
        if (inv_ventricle_val >= thresh)
        {
          ventricle_voxels++ ;
          out_val = out_label ;
        }
        *pdst++ = out_val ;
      }
    }
  }

#if 0
  MRIfillRegion(mri_T1, mri_dst, 30, out_label, 2*ventricle_voxels) ;
  MRIdilate(mri_dst, mri_dst) ;
  MRIdilate(mri_dst, mri_dst) ;
#endif
  return(mri_dst) ;
}
#else
#if 0
static MRI *
MRIfillVentricles(MRI *mri_src, MRI *mri_T1, MRI *mri_mask,
                  MRI *mri_inv_T1, MRI *mri_inv_T1_std, MRI *mri_dst,
                  int out_label)
{
  BUFTYPE   *pmask, *pdst, *psrc, out_val, mask, in_val,
  T1_val, T1_inv_val, sigma, *pT1, *pinvT1, *psigma ;
  int       width, height, depth, x, y, z, nchanged, noff, non,
  ventricle_voxels;
  float     nvox, sdist ;

  if (mri_mask->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "MRIfillVentricles: mask must be MRI_FLOAT")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  /* now apply the inverse morph to build an average wm representation
  of the input volume
  */

  ventricle_voxels = noff = non = nchanged = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pmask = &MRIvox(mri_mask, 0, y, z) ;
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      pT1 = &MRIvox(mri_T1, 0, y, z) ;
      pinvT1 = &MRIvox(mri_inv_T1, 0, y, z) ;
      psigma = &MRIvox(mri_inv_T1_std, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        T1_val = *pT1++ ;
        T1_inv_val = *pinvT1++ ;
        sigma = *psigma++ ;
        if (!FZERO(sigma))
          sdist = (float)abs(T1_val - T1_inv_val)/sigma ;
        else
          sdist = 0 ;
        if (x == 125 && y == 148 && z == 100)
          DiagBreak() ;  /* marked as 255, should be 127 */
        out_val = 0 ;
        mask = *pmask++ ;   /* value from inverse morphed volume */
        in_val = *psrc++ ;
        if (in_val != out_label && mask > 70 && T1_val < 45 &&
            T1_inv_val < 60)
        {
          ventricle_voxels++ ;
          out_val = out_label ;
        }
        else
          out_val = in_val ;
        if (out_val != in_val)
        {
          if (out_val)
            non++ ;
          else
            noff++ ;
          nchanged++ ;
        }
        *pdst++ = out_val ;
      }
    }
  }

  MRIfillRegion(mri_T1, mri_dst, 30, out_label, 2*ventricle_voxels) ;
  MRIdilate(mri_dst, mri_dst) ;
  MRIdilate(mri_dst, mri_dst) ;
  nvox = (float)(width * height * depth) ;
  fprintf(stderr, "%8d of %8d voxels changed - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)nchanged/nvox) ;
  fprintf(stderr, "%8d of %8d voxels off     - %2.1f%%.\n",
          noff, (int)nvox,   100.0f*(float)noff/nvox) ;
  fprintf(stderr, "%8d of %8d voxels on      - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)non/nvox) ;
  return(mri_dst) ;
}
#endif
#endif
