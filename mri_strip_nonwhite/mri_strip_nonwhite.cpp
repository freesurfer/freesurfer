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
#include "mrimorph.h"
#include "transform.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *MRImaskThresholdNeighborhood(MRI *mri_src, MRI *mri_mask,
    MRI *mri_dst,
    float threshold, int nsize) ;

const char *Progname ;

static void usage_exit(int code) ;

#define T1_VOLUME     0
#define WM_VOLUME     1
#define FILLED_VOLUME 2
#define EDIT_VOLUME   3
#define MAX_VOLUMES   4

static float pct = 0.0f ;
static int nsize = 0 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  MRI    *mri, *mri_template, *mri_inverse_template ;
  char   *in_fname, *template_fname, *out_fname, *xform_fname, fname[100] ;
  M3D    *m3d ;
  int     msec ;

  nargs = handleVersionOption(argc, argv, "mri_strip_nonwhite");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Timer start;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  in_fname = argv[1] ;
  xform_fname = argv[2] ;
  template_fname = argv[3] ;
  out_fname = argv[4] ;

  mri = MRIread(in_fname) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s.\n",
              Progname, in_fname) ;

  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, WM_VOLUME*2); /* means and stds */
  else
    strcpy(fname, template_fname) ;
  mri_template = MRIread(fname) ;
  if (!mri_template)
    ErrorExit(ERROR_NOFILE, "%s: could not read template volume %s.\n",
              Progname, template_fname) ;

  fprintf(stderr, "reading transform %s...", xform_fname) ;
  m3d = MRI3DreadSmall(xform_fname) ;
  if (!m3d)
    ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
              Progname, xform_fname) ;
  fprintf(stderr, "done.\n") ;
  fprintf(stderr,
          "template has %d degrees of freedom - setting thresh = %2.1f\n",
          mri_template->dof, pct) ;
  fprintf(stderr, "applying inverse transform...") ;
  mri_inverse_template = MRIapplyInverse3DMorph(mri_template, m3d, NULL) ;
  MRIfree(&mri_template) ;
  MRIwrite(mri_inverse_template, "inverse.mgh") ;
  if (!nsize)  /* don't erase anything close to white */
    nsize = nint(m3d->node_spacing)+1 ;
  MRI3DmorphFree(&m3d) ;
  fprintf(stderr, "done.\nthresholding inverse image (spacing=%d)...",nsize) ;
  MRImaskThresholdNeighborhood(mri, mri_inverse_template, mri, pct, nsize) ;
  fprintf(stderr, "done.\n") ;

  if (MRIwrite(mri, out_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not write output volume to %s.\n",
              Progname, out_fname) ;

  msec = start.milliseconds() ;
  fprintf(stderr, "skull stripping took %2.2f minutes\n",
          (float)msec/(1000.0f*60.0f));
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
  switch (toupper(*option)) {
  case 'N':
    nsize = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "white matter buffer zone set to %d voxels\n",nsize) ;
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
usage_exit(int code) {
  printf("usage: %s <input volume> <transform> <template volume> <output volume>\n",
         Progname) ;
  exit(code) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static MRI *
MRImaskThresholdNeighborhood(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                             float threshold, int nsize) {
  BUFTYPE   *pmask, out_val ;
  int       width, height, depth, x, y, z, x1, y1, z1, xi, yi, zi,
  xmin, xmax, ymin, ymax, zmin, zmax ;

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

  xmin = width ;
  ymin = height ;
  zmin = depth ;
  xmax = ymax = zmax = 0 ;
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pmask = &MRIvox(mri_mask, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        if (*pmask++ > threshold) {
          if (x < xmin)
            xmin = x ;
          if (x > xmax)
            xmax = x ;
          if (y < ymin)
            ymin = y ;
          if (y > ymax)
            ymax = y ;
          if (z < zmin)
            zmin = z ;
          if (z > zmax)
            zmax = z ;
        }
      }
    }
  }
  xmin = MAX(xmin-nsize, 0) ;
  ymin = MAX(ymin-nsize, 0) ;
  zmin = MAX(zmin-nsize, 0) ;
  xmax = MIN(xmax+nsize, width-1) ;
  ymax = MIN(ymax+nsize, height-1) ;
  zmax = MIN(zmax+nsize, depth-1) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "bounding box = (%d:%d, %d:%d, %d:%d).\n",
            xmin, xmax, ymin, ymax, zmin, zmax) ;

  /* remove stuff outside bounding box */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if ((x < xmin || x > xmax) ||
            (y < ymin || y > ymax) ||
            (z < zmin || z > zmax))
          MRIvox(mri_dst, x, y, z) = 0 ;
      }
    }
  }

  for (z = zmin ; z <= zmax ; z++) {
    for (y = ymin ; y <= ymax ; y++) {
      for (x = xmin ; x <= xmax ; x++) {
        out_val = 0 ;
        for (z1 = -nsize ; z1 <= nsize ; z1++) {
          zi = mri_src->zi[z+z1] ;
          for (y1 = -nsize ; y1 <= nsize ; y1++) {
            yi = mri_src->yi[y+y1] ;
            for (x1 = -nsize ; x1 <= nsize ; x1++) {
              xi = mri_src->xi[x+x1] ;
              if (MRIvox(mri_mask, xi, yi, zi) > threshold) {
                out_val = MRIvox(mri_src, x, y, z) ;
                break ;
              }
            }
            if (out_val > 0)
              break ;
          }
          if (out_val > 0)
            break ;
        }
        MRIvox(mri_dst, x, y, z) = out_val ;
      }
    }
  }

  return(mri_dst) ;
}

