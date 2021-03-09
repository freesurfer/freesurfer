/**
 * @brief create an average of a set of volumes
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
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *align_with_average(MRI *mri_src, MRI *mri_avg) ;
static MATRIX *align_pca(MRI *mri_src, MRI *mri_avg) ;
static MATRIX *pca_matrix(MATRIX *m_in_evectors, double in_means[3],
                          MATRIX *m_ref_evectors, double ref_means[3]) ;

const char *Progname ;
static int align = 0 ;
static int window_flag = 0 ;
static MORPH_PARMS  parms ;

static void usage_exit(int code) ;
static int pct = 0 ;
static int fromFile = 0 ;
static int use_abs = 0 ;


static int sqr_images = 0 ;
static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;
static int thresh_low = 0 ;
static int nreductions = 2 ;
static int conform = 1 ;
static int sinc_flag = 1;
static int sinchalfwindow = 3;
static float scale_factor = 0.0 ;

static float binarize_thresh = 0 ;

int  MRIsqrtAndNormalize(MRI *mri, float num) ;
MRI  *MRIsumSquare(MRI *mri1, MRI *mri2, MRI *mri_dst) ;

int
MRIsqrtAndNormalize(MRI *mri, float num)
{
  int    x, y, z ;
  double val ;

  for (x = 0 ; x < mri->width ; x++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri, x, y, z,0) ;
        val = sqrt(val/num) ;
        MRIsetVoxVal(mri, x, y, z, 0, val);
      }
    }
  }

  return(NO_ERROR) ;
}

MRI  *
MRIsumSquare(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int    x, y, z, f ;
  double val1, val2, val_dst ;

  if (!mri_dst)
    mri_dst = MRIclone(mri1, NULL) ;

  for (f=0 ; f < mri1->nframes ; f++)
  {
    for (x = 0 ; x < mri1->width ; x++)
    {
      for (y = 0 ; y < mri1->height ; y++)
      {
        for (z = 0 ; z < mri1->depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          val1 = MRIgetVoxVal(mri1, x, y, z,f) ;
          if (mri2)
            val2 = MRIgetVoxVal(mri2, x, y, z, f) ;
          else
            val2 = 0 ;
          val_dst = val1*val1 + val2 ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val_dst);
        }
      }
    }
  }

  return(mri_dst) ;
}

MRI  *
MRIsumFrames(MRI *mri, MRI *mri_dst, int squares)
{
  int   x, y, z, f ;
  double  val,val_dst ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(mri->width, mri->height, mri->depth, mri->type);
    MRIcopyHeader(mri,mri_dst);
  }

  for (x = 0 ; x < mri->width ; x++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val_dst = 0.0;
        for (f=0 ; f < mri->nframes ; f++)
        {
          val = MRIgetVoxVal(mri, x, y, z,f) ;
          if (squares)
            val_dst += val*val;
          else
            val_dst += val;
        }
        val_dst /= (double)(mri->nframes);
        if (squares) val_dst = sqrt(val_dst);          
        MRIsetVoxVal(mri_dst, x, y, z, 0, val_dst);
      }
    }
  }
  return(mri_dst) ;
}

int
main(int argc, char *argv[])
{
  char   **av, fname[STRLEN] ;
  int    ac, nargs, i, num = 0, filecount;
  MRI    *mri_src, *mri_avg = NULL, *mri_tmp ;
  char   *in_fname, *out_fname, *list_fname ;
  int          msec, minutes, seconds, skipped = 0 ;
  Timer start ;

  nargs = handleVersionOption(argc, argv, "mri_average");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;
  parms.dt = 1e-6 ;
  parms.tol = 1e-5 ;
  parms.momentum = 0.0 ;
  parms.niterations = 20 ;

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
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;

  if (fromFile)
  {

    list_fname = argv[1] ;
    if (! FileExists(list_fname))
      ErrorExit(Gerror, 
                "%s: fopen(%s) for input volume filelist failed",
                Progname, list_fname) ;
    else
    {
      FILE           *fp ;
      fprintf(stderr,
              "reading input volume filenames from %s...\n",
              list_fname) ;
      fp = fopen(list_fname, "r") ;
      if (!fp)
        ErrorExit(Gerror, "Volumelist file cannot be opened\n") ;

      filecount = 0;
      while (fscanf(fp,"%s",fname) != EOF)
      {
        MRI *mri_src_old = NULL;

        fprintf(stderr, "%d of list: reading %s...\n", filecount+1, fname) ;

        /*Pretty ugly but just pasted asthe one below ... (LZ)*/
        mri_src = MRIread(fname) ;
        if (!mri_src)
          ErrorExit(Gerror, "%s: MRIread(%s) failed", Progname, fname) ;

        float src_min, src_max;
        MRIlimits(mri_src, &src_min, &src_max);
        if (src_min >= src_max)
          continue;

        if (binarize_thresh > 0)
          MRIbinarize(mri_src, mri_src, binarize_thresh, 0, 100) ;
        if (pct)
          MRIbinarize(mri_src, mri_src, 1, 0, 100) ;

        if (scale_factor > 0)
          MRIscalarMul(mri_src, mri_src, scale_factor) ;

        if (conform)
        {
          MRI *mri_tmp ;

          fprintf(stderr, "embedding and interpolating volume\n") ;
          mri_tmp = MRIconform(mri_src) ;
          /*      MRIfree(&mri_src) ;*/
          mri_src_old = mri_src; // free it later...
          mri_src = mri_tmp ;
        }

        if (filecount == 1)
        {
          if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz))
          {
            MRI *mri_tmp ;

            fprintf
            (stderr,
             "translating second volume by (%2.1f, %2.1f, %2.1f)\n",
             tx, ty, tz) ;
            mri_tmp = MRItranslate(mri_src, NULL, tx, ty, tz) ;
            MRIfree(&mri_src) ;
            mri_src = mri_tmp ;
          }
          if (!FZERO(rzrot))
          {
            MRI *mri_tmp ;

            fprintf
            (stderr,
             "rotating second volume by %2.1f degrees around Z axis\n",
             (float)DEGREES(rzrot)) ;
            mri_tmp = MRIrotateZ_I(mri_src, NULL, rzrot) ;
            MRIfree(&mri_src) ;
            mri_src = mri_tmp ;
          }
          if (!FZERO(rxrot))
          {
            MRI *mri_tmp ;

            fprintf
            (stderr,
             "rotating second volume by %2.1f degrees around X axis\n",
             (float)DEGREES(rxrot)) ;
            mri_tmp = MRIrotateX_I(mri_src, NULL, rxrot) ;
            MRIfree(&mri_src) ;
            mri_src = mri_tmp ;
          }
          if (!FZERO(ryrot))
          {
            MRI *mri_tmp ;

            fprintf
            (stderr,
             "rotating second volume by %2.1f degrees around Y axis\n",
             (float)DEGREES(ryrot)) ;
            mri_tmp = MRIrotateY_I(mri_src, NULL, ryrot) ;
            MRIfree(&mri_src) ;
            mri_src = mri_tmp ;
          }
        }
        if (align && mri_avg)  /* don't align the first time */
        {
          mri_tmp = align_with_average(mri_src, mri_avg) ;
          MRIfree(&mri_src) ;
          mri_src = mri_tmp ;
        }

        if (mri_avg &&
            ((mri_src->width != mri_avg->width) ||
             (mri_src->height != mri_avg->height) ||
             (mri_src->depth != mri_avg->depth))
           )
        {
          printf("src image (%d, %d, %d) not compatible with "
                 "avg image (%d, %d, %d)\n",
                 mri_src->width, mri_src->height, mri_src->depth,
                 mri_avg->width, mri_avg->height, mri_avg->depth) ;
          skipped++ ;
          continue ;
        }
        if (mri_avg == NULL)
        {
          mri_avg = MRIallocSequence(mri_src->width,
                                     mri_src->height,
                                     mri_src->depth,
                                     MRI_FLOAT,
                                     mri_src->nframes) ;
          MRIcopyHeader(mri_src, mri_avg) ;
        }
        if (sqr_images)
          MRIsumSquare(mri_src, mri_avg, mri_avg) ;
        else
          MRIaverage(mri_src, filecount-skipped, mri_avg) ;
        MRIfree(&mri_src) ;
        if (mri_src_old) MRIfree(&mri_src_old) ;

        filecount ++;
        /* End of reading all the files in */
      }
      fclose(fp);
    }

  }
  else
  {
    for (num = 0, i = 1 ; i < argc-1 ; i++)
    {
      MRI *mri_src_old = NULL;

      in_fname = argv[i] ;
      fprintf(stderr, "%d of %d: reading %s...\n", num+1, argc-2, in_fname) ;

      mri_src = MRIread(in_fname) ;
      if (!mri_src)
        ErrorExit(Gerror, "%s: MRIread(%s) failed", Progname, in_fname) ;
      if (binarize_thresh > 0)
        MRIbinarize(mri_src, mri_src, binarize_thresh, 0, 100) ;
      if (pct)
        MRIbinarize(mri_src, mri_src, 1, 0, 100) ;

      if (scale_factor > 0)
        MRIscalarMul(mri_src, mri_src, scale_factor) ;

      if (conform)
      {
        MRI *mri_tmp ;

        fprintf(stderr, "embedding and interpolating volume\n") ;
        mri_tmp = MRIconform(mri_src) ;
        /*      MRIfree(&mri_src) ;*/
        mri_src_old = mri_src; // free it later...
        mri_src = mri_tmp ;
      }

      if (i == 2)
      {
        if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz))
        {
          MRI *mri_tmp ;

          fprintf
          (stderr,
           "translating second volume by (%2.1f, %2.1f, %2.1f)\n",
           tx, ty, tz) ;
          mri_tmp = MRItranslate(mri_src, NULL, tx, ty, tz) ;
          MRIfree(&mri_src) ;
          mri_src = mri_tmp ;
        }
#if 1
        if (!FZERO(rzrot))
        {
          MRI *mri_tmp ;

          fprintf
          (stderr,
           "rotating second volume by %2.1f degrees around Z axis\n",
           (float)DEGREES(rzrot)) ;
          mri_tmp = MRIrotateZ_I(mri_src, NULL, rzrot) ;
          MRIfree(&mri_src) ;
          mri_src = mri_tmp ;
        }
        if (!FZERO(rxrot))
        {
          MRI *mri_tmp ;

          fprintf
          (stderr,
           "rotating second volume by %2.1f degrees around X axis\n",
           (float)DEGREES(rxrot)) ;
          mri_tmp = MRIrotateX_I(mri_src, NULL, rxrot) ;
          MRIfree(&mri_src) ;
          mri_src = mri_tmp ;
        }
        if (!FZERO(ryrot))
        {
          MRI *mri_tmp ;

          fprintf
          (stderr,
           "rotating second volume by %2.1f degrees around Y axis\n",
           (float)DEGREES(ryrot)) ;
          mri_tmp = MRIrotateY_I(mri_src, NULL, ryrot) ;
          MRIfree(&mri_src) ;
          mri_src = mri_tmp ;
        }
#else
        if (!FZERO(ryrot) || !FZERO(rxrot) || !FZERO(rzrot))
        {
          MRI *mri_tmp ;
          MATRIX *mX, *mY, *mZ, *mRot, *mTmp ;

          mX = MatrixAllocRotation(3, x_angle, X_ROTATION) ;
          mY = MatrixAllocRotation(3, y_angle, Y_ROTATION) ;
          mZ = MatrixAllocRotation(3, z_angle, Z_ROTATION) ;
          mTmp = MatrixMultiply(mX, mZ, NULL) ;
          mRot = MatrixMultiply(mY, mTmp, NULL);
          fprintf
          (stderr,
           "rotating second volume by (%2.1f, %2.1f, %2.1f) degrees\n",
           (float)DEGREES(rxrot), (float)DEGREES(ryrot)
           (float)DEGREES(rzrot)) ;

          mri_tmp = MRIrotate_I(mri_src, NULL, mRot, NULL) ;
          MRIfree(&mri_src) ;
          mri_src = mri_tmp ;

          MatrixFree(&mX) ;
          MatrixFree(&mY) ;
          MatrixFree(&mZ) ;
          MatrixFree(&mTmp) ;
          MatrixFree(&mRot) ;
        }
#endif

#if 0
        if (!FZERO(rxrot) || !FZERO(ryrot) || !FZERO(rzrot))
          MRIwrite(mri_src, "/disk2/mri/tamily/mri/tmp") ;
#endif
      }
#if 0
      mri_src->xsize =
        mri_src->ysize =
          mri_src->zsize =
            mri_src->thick = 1.0f ;
      mri_src->imnr0 = 1 ;
      mri_src->imnr1 = mri_src->depth ;
#endif
      if (align && mri_avg)  /* don't align the first time */
      {
        mri_tmp = align_with_average(mri_src, mri_avg) ;
        MRIfree(&mri_src) ;
        mri_src = mri_tmp ;
      }

      num++ ;
      if (mri_avg &&
          ((mri_src->width != mri_avg->width) ||
           (mri_src->height != mri_avg->height) ||
           (mri_src->depth != mri_avg->depth))
         )
      {
        fprintf(stderr,"src image (%d, %d, %d) not compatible with "
               "avg image (%d, %d, %d)\n",
               mri_src->width, mri_src->height, mri_src->depth,
               mri_avg->width, mri_avg->height, mri_avg->depth) ;
        skipped++ ;
        continue ;
      }
      
      /* if single input, average across frames */
      if (argc ==3)
      {
        fprintf(stderr, "Single input, working on %d input frames instead.\n", mri_src->nframes) ;
        mri_avg = MRIsumFrames(mri_src,NULL,sqr_images);        
      }
      else
      {
      
        if (mri_avg == NULL)
        {
          mri_avg = MRIallocSequence(mri_src->width,
                                     mri_src->height,
                                     mri_src->depth,
                                     MRI_FLOAT,
                                     mri_src->nframes) ;
          MRIcopyHeader(mri_src, mri_avg) ;
        }
        if (sqr_images)
          MRIsumSquare(mri_src, mri_avg, mri_avg) ;
        else
          MRIaverage(mri_src, (i-1)-skipped, mri_avg) ;
      }
      MRIfree(&mri_src) ;
      if (mri_src_old) MRIfree(&mri_src_old) ;
    }
  }

  if (sqr_images && argc > 3)
    MRIsqrtAndNormalize(mri_avg, num) ;

  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_avg, out_fname) ;
  MRIfree(&mri_avg) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "alignment and averaging took %d minutes and %d seconds.\n",
          minutes, seconds) ;
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
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using dt = %2.3e\n", parms.dt) ;
  }
  else if (!stricmp(option, "tol"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.3e\n", parms.tol) ;
  }
  else if (!stricmp(option, "sqr") || !stricmp(option, "rms"))
  {
    sqr_images = 1 ;
    fprintf
    (stderr, "computing sqrt of sum of squares instead of average (RMS)...\n") ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    fprintf(stderr, "interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "reduce"))
  {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "reducing input images %d times before aligning...\n",
            nreductions) ;
  }
  else if (!stricmp(option, "sinc"))
  {
    sinchalfwindow = atoi(argv[2]);
    sinc_flag = 1;
    nargs = 1;
    fprintf(stderr,"using sinc interpolation with windowwidth of %d\n",
            2*sinchalfwindow);
  }
  else if (!stricmp(option, "trilinear"))
  {
    sinc_flag = 0;
    fprintf(stderr,"using trilinear interpolation\n");
  }
  else if (!stricmp(option, "abs"))
  {
    use_abs = 1 ;
    fprintf(stderr,"taking abs value of input volumes\n");
  }
  else if (!stricmp(option, "window"))
  {
    window_flag = 1 ;
    fprintf(stderr, "applying hanning window to volumes...\n") ;
  }
  else if (!stricmp(option, "help"))
  {
    usage_exit(0) ;
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    fprintf(stderr, "inhibiting isotropic volume interpolation\n") ;
  }
  else switch (toupper(*option))
    {
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      Gdiag |= DIAG_WRITE ;
      nargs = 1 ;
      fprintf(stderr, "writing snapshots every %d iterations\n",
              parms.write_iterations) ;
      break ;
    case 'S':
      scale_factor = atof(argv[2]) ;
      nargs = 1 ;
      printf("scaling all volumes by %f\n", scale_factor) ;
      break ;
    case 'T':
      tx = atof(argv[2]) ;
      ty = atof(argv[3]) ;
      tz = atof(argv[4]) ;
      nargs = 3 ;
      break ;
    case 'P':
      pct = 1 ;
      printf("binarizing images to compute pct at each voxel\n") ;
      break ;
    case 'F':
      fromFile = 1 ;
      printf("read input volumes from file\n") ;
      break ;
    case 'B':
      binarize_thresh = atof(argv[2]) ;
      printf("binarizing images with thresh %2.1f to compute pct at each voxel\n", binarize_thresh) ;
      nargs = 1 ;
      break ;
    case 'R':
      rxrot = RADIANS(atof(argv[2])) ;
      ryrot = RADIANS(atof(argv[3])) ;
      rzrot = RADIANS(atof(argv[4])) ;
      nargs = 3 ;
      break ;
    case 'M':
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using momentum = %2.3f\n", parms.momentum) ;
      break ;
    case 'A':
      align = 1 ;
      fprintf(stderr, "aligning volumes before averaging...\n") ;
      break ;
    case 'H':
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
  printf("usage: %s [options] <volume> ... <output volume>\n", Progname) ;
  printf("\t-a              rigid alignment of input "
         "volumes before averaging\n") ;
  printf("\t-F              read volumes from an input file (first argument is the input filename)\n") ;
  printf("\t-dt <float n>   set dt to n (default=1e-6)\n") ;
  printf("\t-tol <float n>  set tol to n (default=1e-5)\n") ;
  printf("\t-conform        interpolate volume to be isotropic 1mm^3 (this option is on by default)\n") ;
  printf("\t-noconform      inhibit isotropic volume interpolation\n");
  printf("\t-reduce <int n> reduce input images n (default=2) times\n") ;
  printf("\t-sinc <int n>   using sinc "
         "interpolation with windowwidth of 2*n (default=3)\n") ;
  printf("\t-trilinear      use trilinear interpolation\n");
  printf("\t-window         apply hanning window to volumes\n");
  printf("\t-w <int n>      write snapshots every n iterations\n");
  printf("\t-t <x> <y> <z>  translation of second volume\n");
  printf("\t-r <x> <y> <z>  rotation of "
         "second volume around each axis in degrees\n");
  printf("\t-m <float n>    use momentum n (default=0)\n");
  printf("\t-sqr            compute sqrt of average of sum of squares (RMS, same as -rms)\n") ;
  printf("\t-rms            compute sqrt of average of sum of squares (RMS, same as -sqr)\n") ;
  printf("\t-u              print usage\n");
  printf("\t-p              compute %% \n");
  printf("\t-b <float th>   binarize the input volumes using threshold th \n");
  printf("\t-abs            take abs value of volume \n");
  exit(code) ;
}

static MRI *
align_with_average(MRI *mri_src, MRI *mri_avg)
{
  MRI     *mri_aligned, *mri_in_red, *mri_ref_red ;
  MRI     *mri_in_windowed, *mri_ref_windowed, *mri_in_tmp, *mri_ref_tmp ;
  int     i ;
  MATRIX  *m_L ;

  fprintf(stderr, "initializing alignment using PCA...\n") ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwriteImageViews(mri_avg, "ref", 400) ;
    MRIwriteImageViews(mri_src, "before_pca", 400) ;
  }

  m_L = align_pca(mri_src, mri_avg) ;
  if (Gdiag & DIAG_SHOW)
  {
    printf("initial transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }
  if (Gdiag & DIAG_WRITE)
  {
    if (sinc_flag)
      mri_aligned = MRIsincTransform(mri_src, NULL, m_L,sinchalfwindow) ;
    else
      mri_aligned = MRIlinearTransform(mri_src, NULL, m_L) ;
    MRIwriteImageViews(mri_aligned, "after_pca", 400) ;
    MRIfree(&mri_aligned) ;
  }

  fprintf(stderr, "aligning volume with average...\n") ;

  if (window_flag)
  {
    mri_in_windowed =
      MRIwindow(mri_src, NULL, WINDOW_HANNING,127,127,127,100.0f);
    mri_ref_windowed =
      MRIwindow(mri_avg,NULL,WINDOW_HANNING,127,127,127,100.0f);
    mri_src = mri_in_windowed ;
    mri_avg = mri_ref_windowed ;
  }

  MRIscaleMeanIntensities(mri_src, mri_avg, mri_src);

  mri_in_red = mri_in_tmp = MRIcopy(mri_src, NULL) ;
  mri_ref_red = mri_ref_tmp = MRIcopy(mri_avg, NULL) ;
  for (i = 0 ; i < nreductions ; i++)
  {
    mri_in_red = MRIreduceByte(mri_in_tmp, NULL) ;
    mri_ref_red = MRIreduceByte(mri_ref_tmp,NULL);
    MRIfree(&mri_in_tmp);
    MRIfree(&mri_ref_tmp) ;
    mri_in_tmp = mri_in_red ;
    mri_ref_tmp = mri_ref_red ;
  }
  parms.mri_ref = mri_avg ;
  parms.mri_in = mri_src ;  /* for diagnostics */
  MRIrigidAlign(mri_in_red, mri_ref_red, &parms, m_L) ;

  fprintf(stderr, "transforming input volume...\n") ;
  MatrixPrint(stderr, parms.lta->xforms[0].m_L) ;
  fprintf(stderr, "\n") ;

  if (sinc_flag)
    mri_aligned =
      MRIsincTransform
      (mri_src, NULL, parms.lta->xforms[0].m_L,sinchalfwindow) ;
  else
    mri_aligned = MRIlinearTransform(mri_src, NULL, parms.lta->xforms[0].m_L) ;
  if (Gdiag & DIAG_WRITE)
    MRIwriteImageViews(mri_aligned, "after_alignment", 400) ;
  MRIfree(&mri_in_red) ;
  MRIfree(&mri_ref_red) ;

  return(mri_aligned) ;
}


static MATRIX *
align_pca(MRI *mri_in, MRI *mri_ref)
{
  int    row, col, i ;
  float  dot ;
  MATRIX *m_ref_evectors = NULL, *m_in_evectors = NULL ;
  float  in_evalues[3], ref_evalues[3] ;
  double  ref_means[3], in_means[3] ;
#if 0
  MRI     *mri_in_windowed, *mri_ref_windowed ;

  mri_in_windowed = MRIwindow(mri_in, NULL, WINDOW_HANNING,127,127,127,100.0f);
  mri_ref_windowed = MRIwindow(mri_ref,NULL,WINDOW_HANNING,127,127,127,100.0f);
  if (Gdiag & DIAG_WRITE)
  {
    MRIwriteImageViews(mri_in_windowed, "in_windowed", 400) ;
    MRIwriteImageViews(mri_ref_windowed, "ref_windowed", 400) ;
  }
#endif

  if (!m_ref_evectors)
    m_ref_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;
  if (!m_in_evectors)
    m_in_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;

  MRIprincipleComponents(mri_ref, m_ref_evectors, ref_evalues,
                         ref_means, thresh_low);
  MRIprincipleComponents(mri_in,m_in_evectors,in_evalues,in_means,thresh_low);

  /* check to make sure eigenvectors aren't reversed */
  for (col = 1 ; col <= 3 ; col++)
  {
#if 0
    float theta ;
#endif

    for (dot = 0.0f, row = 1 ; row <= 3 ; row++)
      dot += m_in_evectors->rptr[row][col] * m_ref_evectors->rptr[row][col] ;

    if (dot < 0.0f)
    {
      fprintf(stderr,
              "WARNING: mirror image detected in eigenvector #%d\n",
              col) ;
      dot *= -1.0f ;
      for (row = 1 ; row <= 3 ; row++)
        m_in_evectors->rptr[row][col] *= -1.0f ;
    }
#if 0
    theta = acos(dot) ;
    fprintf(stderr, "angle[%d] = %2.1f\n", col, DEGREES(theta)) ;
#endif
  }
  fprintf(stderr, "ref_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_ref_evectors->rptr[i][1],
            m_ref_evectors->rptr[i][2],
            m_ref_evectors->rptr[i][3]) ;

  fprintf(stderr, "\nin_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_in_evectors->rptr[i][1],
            m_in_evectors->rptr[i][2],
            m_in_evectors->rptr[i][3]) ;

  return(pca_matrix(m_in_evectors, in_means,m_ref_evectors, ref_means)) ;
}

static MATRIX *
pca_matrix(MATRIX *m_in_evectors, double in_means[3],
           MATRIX *m_ref_evectors, double ref_means[3])
{
  float   dx, dy, dz ;
  MATRIX  *mRot, *m_in_T, *mOrigin, *m_L, *m_R, *m_T, *m_tmp ;
  double  x_angle, y_angle, z_angle, r11, r21, r31, r32, r33, cosy ;
  int     row, col ;

  m_in_T = MatrixTranspose(m_in_evectors, NULL) ;
  mRot = MatrixMultiply(m_ref_evectors, m_in_T, NULL) ;

  r11 = mRot->rptr[1][1] ;
  r21 = mRot->rptr[2][1] ;
  r31 = mRot->rptr[3][1] ;
  r32 = mRot->rptr[3][2] ;
  r33 = mRot->rptr[3][3] ;
  y_angle = atan2(-r31, sqrt(r11*r11+r21*r21)) ;
  cosy = cos(y_angle) ;
  z_angle = atan2(r21 / cosy, r11 / cosy) ;
  x_angle = atan2(r32 / cosy, r33 / cosy) ;

#define MAX_ANGLE  (RADIANS(30))
  if (fabs(x_angle) > MAX_ANGLE || fabs(y_angle) > MAX_ANGLE ||
      fabs(z_angle) > MAX_ANGLE)
  {
    MatrixFree(&m_in_T) ;
    MatrixFree(&mRot) ;
    fprintf(stderr, "eigenvector swap detected: ignoring PCA...\n") ;
    return(MatrixIdentity(4, NULL)) ;
  }

  mOrigin = VectorAlloc(3, MATRIX_REAL) ;
  mOrigin->rptr[1][1] = ref_means[0] ;
  mOrigin->rptr[2][1] = ref_means[1] ;
  mOrigin->rptr[3][1] = ref_means[2] ;

  fprintf(stderr, "reference volume center of mass at (%2.1f,%2.1f,%2.1f)\n",
          ref_means[0], ref_means[1], ref_means[2]) ;
  fprintf(stderr, "input volume center of mass at     (%2.1f,%2.1f,%2.1f)\n",
          in_means[0], in_means[1], in_means[2]) ;
  dx = ref_means[0] - in_means[0] ;
  dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;

  fprintf(stderr, "translating volume by %2.1f, %2.1f, %2.1f\n",
          dx, dy, dz) ;
  fprintf(stderr, "rotating volume by (%2.2f, %2.2f, %2.2f)\n",
          DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;

  /* build full rigid transform */
  m_R = MatrixAlloc(4,4,MATRIX_REAL) ;
  m_T = MatrixAlloc(4,4,MATRIX_REAL) ;
  for (row = 1 ; row <= 3 ; row++)
  {
    for (col = 1 ; col <= 3 ; col++)
    {
      *MATRIX_RELT(m_R,row,col) = *MATRIX_RELT(mRot, row, col) ;
    }
    *MATRIX_RELT(m_T,row,row) = 1.0 ;
  }
  *MATRIX_RELT(m_R, 4, 4) = 1.0 ;

  /* translation so that origin is at ref eigenvector origin */
  dx = -ref_means[0] ;
  dy = -ref_means[1] ;
  dz = -ref_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ;
  *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ;
  *MATRIX_RELT(m_T, 4, 4) = 1 ;
  m_tmp = MatrixMultiply(m_R, m_T, NULL) ;
  *MATRIX_RELT(m_T, 1, 4) = -dx ;
  *MATRIX_RELT(m_T, 2, 4) = -dy ;
  *MATRIX_RELT(m_T, 3, 4) = -dz ;
  MatrixMultiply(m_T, m_tmp, m_R) ;

  /* now apply translation to take in centroid to ref centroid */
  dx = ref_means[0] - in_means[0] ;
  dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ;
  *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ;
  *MATRIX_RELT(m_T, 4, 4) = 1 ;

  m_L = MatrixMultiply(m_R, m_T, NULL) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("m_T:\n") ;
    MatrixPrint(stdout, m_T) ;
    printf("m_R:\n") ;
    MatrixPrint(stdout, m_R) ;
    printf("m_L:\n") ;
    MatrixPrint(stdout, m_L) ;
  }
  MatrixFree(&m_R) ;
  MatrixFree(&m_T) ;

  MatrixFree(&mRot) ;
  VectorFree(&mOrigin) ;
  return(m_L) ;
}
