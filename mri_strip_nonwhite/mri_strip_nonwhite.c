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

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *MRImaskThreshold(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, 
                             float threshold) ;

char *Progname ;

static void usage_exit(int code) ;

#define T1_VOLUME     0
#define WM_VOLUME     1
#define FILLED_VOLUME 2
#define EDIT_VOLUME   3
#define MAX_VOLUMES   4

static float threshold = 1.0f ;


int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, i ;
  MRI    *mri, *mri_template, *mri_inverse_template ;
  char   *in_fname, *template_fname, *out_fname, *xform_fname, fname[100] ;
  M3D    *m3d ;

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
    usage_exit(1) ;

  in_fname = argv[1] ; xform_fname = argv[2] ;
  template_fname = argv[3] ; out_fname = argv[4] ;

  mri = MRIread(in_fname) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s.\n",
              Progname, in_fname) ;

  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, WM_VOLUME*2) ; /* means and stds */
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
  fprintf(stderr, "done.\napplying inverse transform...") ;
  mri_inverse_template = MRIapplyInverse3DMorph(mri_template, m3d, NULL) ;
  MRIfree(&mri_template) ;
  MRIwrite(mri_inverse_template, "inverse.mgh") ;
  MRI3DmorphFree(&m3d) ;
  fprintf(stderr, "done.\nthresholding inverse image...") ;
  MRImaskThreshold(mri, mri_inverse_template, mri, threshold) ;
  fprintf(stderr, "done.\n") ;

  if (MRIwrite(mri, out_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not write output volume to %s.\n",
              Progname, out_fname) ;

  fprintf(stderr, "done.\n") ;
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
  switch (toupper(*option))
  {
  case 'T':
    threshold = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using threshold %2.1f\n", threshold) ;
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
  printf("usage: %s <input volume> <transform> <template volume> <output volume>\n", 
         Progname) ;
  exit(code) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static MRI *
MRImaskThreshold(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, float threshold)
{
  BUFTYPE   *psrc, *pdst, *pmask ;
  int       width, height, depth, x, y, z ;

  if (mri_mask->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, 
                       "MRI3Dthreshold: mask must be MRI_FLOAT")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 


  /* now apply the inverse morph to build an average wm representation
     of the input volume 
     */

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pmask = &MRIvox(mri_mask, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        if (*pmask++ > threshold)
          *pdst++ = *psrc++ ;
        else
        {
          *pdst++ = 0 ;
          psrc++ ;
        }
      }
    }
  }

  return(mri_dst) ;
}

