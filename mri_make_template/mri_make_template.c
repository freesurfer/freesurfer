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

char *Progname ;

static void usage_exit(int code) ;
int MRIaccumulateMeansAndVariances(MRI *mri, MRI *mri_mean, MRI *mri_std) ;
int MRIcomputeMeansAndStds(MRI *mri_mean, MRI *mri_std, int ndof) ;
MRI *MRIfloatToChar(MRI *mri_src, MRI *mri_dst) ;

static int use_wm = 0 ;
static int use_filled = 0 ;
static int make_edit_volume = 0 ;
static char *transform_fname = NULL ;

#define T1_VOLUME     0
#define WM_VOLUME     1
#define FILLED_VOLUME 2
#define EDIT_VOLUME   3
#define MAX_VOLUMES   4

static int first_transform = 0 ;

int
main(int argc, char *argv[])
{
  char   **av, *volume_name, subjects_dir[100], *cp ;
  int    ac, nargs, i, which, dof, no_transform ;
  MRI    *mri, *mri_mean = NULL, *mri_std ;
  char   *subject_name, *out_fname, fname[100] ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment.\n",
              Progname) ;
  strcpy(subjects_dir, cp) ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;

  for (which = 0 ; which < MAX_VOLUMES ; which++) /* 1st T1, then wm volumes */
  {
    if ((which > 0 && !use_wm) || (which >= FILLED_VOLUME && !use_filled) ||
        (which >= EDIT_VOLUME  && !make_edit_volume))
      break ;
    volume_name = which == T1_VOLUME ? "T1" : 
      which == FILLED_VOLUME ? "filled" : "wm" ;

    /* for each subject specified on cmd line */
    no_transform = first_transform ;
    for (dof = 0, i = 1 ; i < argc-1 ; i++) 
    {
      if (*argv[i] == '-')   /* don't do transform for next subject */
      { no_transform = 1 ; continue ; }
      dof++ ;
      subject_name = argv[i] ;
      sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, volume_name) ;
      fprintf(stderr, "%d of %d: reading %s...", i, argc-2, fname) ;
      mri = MRIread(fname) ;
      if (!mri)
        ErrorExit(ERROR_NOFILE, "%s: could not open volume %s",Progname,fname);
      if (which == WM_VOLUME)  /* labeled white matter volume - binarize it */
        MRIbinarize(mri, mri, 5, 0, 10) ;
      else if (which == FILLED_VOLUME)
      {
        MRIreplaceValues(mri, mri, 127, 1) ;
        MRIreplaceValues(mri, mri, 255, 1) ;
        MRIreplaceValues(mri, mri, 128, 1) ;
      }
      else if (which == EDIT_VOLUME)
        MRIbinarize(mri, mri, 255, 0, 1) ; /* only edits will survive */
      fprintf(stderr, "done.\n") ;
      if (transform_fname && no_transform-- <= 0)
      {
        int       type ;
        MORPH_3D  *m3d ;
        LTA       *lta ;
        MRI       *mri_tmp ;
        
        sprintf(fname, "%s/%s/mri/transforms/%s", 
                subjects_dir, subject_name, transform_fname) ;

        fprintf(stderr, "reading transform %s...", fname) ;
        type = TransformFileNameType(fname) ;
        switch (type)
        {
        default:
        case TRANSFORM_ARRAY_TYPE:
          lta = LTAread(fname) ;
          if (!lta)
            ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                      Progname, fname) ;
          mri_tmp = LTAtransform(mri, NULL, lta) ;
          LTAfree(&lta) ;
          break ;
        case MORPH_3D_TYPE:
          m3d = MRI3DreadSmall(fname) ;
          if (!m3d)
            ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                      Progname, transform_fname) ;
          fprintf(stderr, "done.\napplying transform...") ;
          mri_tmp = MRIapply3DMorph(mri, m3d, NULL) ;
          MRI3DmorphFree(&m3d) ;
          break ;
        }
        MRIfree(&mri) ; mri = mri_tmp ;
        fprintf(stderr, "done.\n") ;
      }
      if (!mri_mean)
      {
        mri_mean = MRIalloc(mri->width, mri->height, mri->depth, MRI_FLOAT) ;
        mri_std = MRIalloc(mri->width, mri->height, mri->depth, MRI_FLOAT) ;
        if (!mri_mean || !mri_std)
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate templates.\n",
                    Progname) ;
      }

      MRIaccumulateMeansAndVariances(mri, mri_mean, mri_std) ;
      MRIfree(&mri) ;
    }
    MRIcomputeMeansAndStds(mri_mean, mri_std, dof) ;
#if 1
    fprintf(stderr, "\nwriting %s means to %s...", volume_name, out_fname) ;
    mri_mean->dof = dof ;
    mri = MRIfloatToChar(mri_mean, NULL) ;
    if (!which)
      MRIwrite(mri, out_fname) ;
    else
      MRIappend(mri, out_fname) ;
    MRIfree(&mri_mean) ;
    fprintf(stderr, "\nwriting %s variances to %s...", volume_name, out_fname);
    mri = MRIfloatToChar(mri_std, NULL) ;
    if (dof <= 1) /* can't calulate variances - set them to reasonable val */
      MRIreplaceValues(mri, mri, 0, 1) ;
    MRIappend(mri, out_fname) ;
    fprintf(stderr, "done.\n") ;
#endif
    MRIfree(&mri_std) ;
  }

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
  case 'E':
    use_wm = make_edit_volume = 1 ;
    fprintf(stderr, 
            "appending wm volume to template and building auto-edit volue.\n");
    break ;
  case 'F':
    use_filled = 1 ;
    use_wm = 1 ;
    fprintf(stderr, "appending wm and filled volumes to template.\n") ;
    break ;
  case 'W':
    use_wm = 1 ;
    fprintf(stderr, "appending wm volume to template.\n") ;
    break ;
  case 'N':
    first_transform = 1 ;  /* don't use transform on first volume */
    break ;
  case 'T':
    transform_fname = argv[2] ;
    fprintf(stderr, "applying transformation %s to each volume\n", 
            transform_fname) ;
    nargs = 1 ;
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
  printf("usage: %s <subject> <subject> ... <output volume>\n", Progname) ;
  exit(code) ;
}

int
MRIaccumulateMeansAndVariances(MRI *mri, MRI *mri_mean, MRI *mri_std)
{
  int    x, y, z, width, height, depth ;
  float  val ;
  
  width = mri->width ; height = mri->height ; depth = mri->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == 128 && y == 128 && z == 128)
          DiagBreak() ;
        val = MRIvox(mri,x,y,z) ;
        MRIFvox(mri_mean,x,y,z) += val ;
        MRIFvox(mri_std,x,y,z) += val*val ;
      }
    }
  }
  return(NO_ERROR) ;
}

int
MRIcomputeMeansAndStds(MRI *mri_mean, MRI *mri_std, int ndof)
{
  int    x, y, z, width, height, depth ;
  float  sum, sum_sq, mean, var ;
  
  width = mri_std->width ; height = mri_std->height ; depth = mri_std->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == 128 && y == 128 && z == 128)
          DiagBreak() ;
        sum = MRIFvox(mri_mean,x,y,z) ;
        sum_sq = MRIFvox(mri_std,x,y,z) / ndof ;
        mean = MRIFvox(mri_mean,x,y,z) = sum / ndof ; 
        var = sum_sq - mean*mean; 
        MRIFvox(mri_std,x,y,z) = sqrt(var) ;
      }
    }
  }
  return(NO_ERROR) ;
}

MRI *
MRIfloatToChar(MRI *mri_src, MRI *mri_dst)
{
  int   width, height, depth/*, x, y, z, out_val*/ ;
  /*  float fmax, fmin ;*/
  
  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 
  mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
#if 1
  MRIcopy(mri_src, mri_dst) ;
#else
  MRIvalRange(mri_src, &fmin, &fmax) ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
      }
    }
  }
#endif
  return(mri_dst) ;
}

#if 0
MRI *
MRIbinarizeEditting(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, val ;
  BUFTYPE *
  
  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 
  mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIvalRange(mri_src, &fmin, &fmax) ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
      }
    }
  }

  return(mri_dst) ;
}
#endif
