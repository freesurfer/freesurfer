#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "transform.h"
#include "mrimorph.h"
#include "fio.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

static int verbose = 0 ;
static int xdim = XDIM, ydim = YDIM, zdim = ZDIM ;

static int brueker = 0 ;
static float scale = 1.0 ;
static float blur_sigma = 0.0f ;
static char *inverse_transform_fname = NULL ;
static char *transform_fname = NULL ;


MRI   *MRIreadBrueker(char *fname) ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  MRI    *mri, *mri_kernel ;
  char   *in_fname, *out_fname ;

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

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  if (argc < 3)
    ErrorExit(ERROR_BADPARM, "%s: no output name specified", Progname) ;
  if (argc > 3)
    ErrorExit(ERROR_BADPARM, "%s: too many command line parameters (%d)", Progname,argc) ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  if (verbose)
    fprintf(stderr, "reading from %s...", in_fname) ;
  if (!brueker)
    mri = MRIread(in_fname) ;
  else
    mri = MRIreadBrueker(in_fname) ;

  if (!mri)
    ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s", 
              Progname, in_fname) ;

  if (!mri->imnr0)
  {
    mri->imnr0++ ;
    mri->imnr1++ ;
  }
  if (xdim != XDIM || ydim != YDIM || zdim != ZDIM)
  {
    MRI  *mri_tmp ;

    fprintf(stderr,"reordering dimensions to (%d, %d, %d)\n",xdim,ydim,zdim);
    mri_tmp = MRIreorder(mri, NULL, xdim, ydim, zdim) ;
    MRIfree(&mri) ;
    mri = mri_tmp ;
  }

  if (blur_sigma > 0.0f)
  {
    MRI *mri_smooth ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    fprintf(stderr, "smoothing volume with sigma = %2.3f\n", blur_sigma) ;
    mri_smooth = MRIclone(mri, NULL) ;
    MRIconvolveGaussian(mri, mri_smooth, mri_kernel) ;
    MRIfree(&mri_kernel) ; MRIfree(&mri) ; 
    mri = mri_smooth ;
  }
  if (verbose)
    fprintf(stderr, "done.\nwriting to %s...", out_fname) ;

  if (transform_fname || inverse_transform_fname)
  {
    MRI *mri_tmp ;
    int type, inverse ;
    char *fname ;
    LTA  *lta ;
    M3D  *m3d ;

    if (inverse_transform_fname)
    {
      inverse = 1 ;
      fname = inverse_transform_fname ;
    }
    else
    {
      inverse = 0 ;
      fname = transform_fname ;
    }
    type = TransformFileNameType(fname) ;
    switch (type)
    {
    default:
    case MNI_TRANSFORM_TYPE:
    case TRANSFORM_ARRAY_TYPE:
      lta = LTAread(fname) ;
      if (!lta)
        ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n", Progname, fname) ;
      if (inverse)
      {
        MATRIX *m_L_inv ;
        m_L_inv = MatrixInverse(lta->xforms[0].m_L, NULL) ;
        MatrixFree(&lta->xforms[0].m_L) ;
        lta->xforms[0].m_L = m_L_inv ;
      }
      mri_tmp = LTAtransform(mri, NULL, lta) ;
      LTAfree(&lta) ;
      break ;
    case MORPH_3D_TYPE:
      m3d = MRI3DreadSmall(fname) ;
      if (!m3d)
        ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                  Progname, fname) ;
      if (inverse)
        mri_tmp = MRIapplyInverse3DMorph(mri, m3d, NULL) ;
      else
        mri_tmp = MRIapply3DMorph(mri, m3d, NULL) ;
      MRI3DmorphFree(&m3d) ;
      break ;
    }
    MRIfree(&mri) ;
    mri = mri_tmp ;
  }

  fprintf(stderr, "writing output to '%s'...\n", out_fname) ;
  MRIwrite(mri, out_fname) ;
  if (verbose)
    fprintf(stderr, "done.\n") ;
  MRIfree(&mri) ;
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
  int  nargs = 0, sign = 1 ;
  char *option, *cp ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "blur"))
  {
    blur_sigma = atof(argv[2]) ;
    fprintf(stderr, "applying %2.2f standard deviation Gaussian kernel\n",
            blur_sigma) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
  case 'T':
    transform_fname = argv[2] ;
    fprintf(stderr, "applying transformation from %s before writing image.\n",
            transform_fname) ;
    nargs = 1 ;
    break ;
  case 'I':
    inverse_transform_fname = argv[2] ;
    fprintf(stderr, "inverting and applying transformation from %s before writing image.\n",
            inverse_transform_fname) ;
    nargs = 1 ;
    break ;
  case 'B':
    brueker = 1 ;
    break ;
  case 'V':
    verbose = !verbose ;
    break ;
  case 'N':
#if 0
    sscanf(argv[2], "%d", &reductions) 
    fprintf(stderr, "reducing %d times\n", reductions) ;
#endif
    nargs = 1 ;
    break ;
  case 'R':
    xdim = atoi(argv[2]) ;
    ydim = atoi(argv[3]) ;
    zdim = atoi(argv[4]) ;
    nargs = 3 ;
    break ;
  case 'X':
    cp = argv[2] ;
    if (*cp == '-')
    {
      cp++ ;  /* skip minus sign */
      sign = -1 ;
    }
    switch (toupper(*cp))
    {
    default:
    case 'X': xdim = XDIM ; break ;
    case 'Y': xdim = YDIM ; break ;
    case 'Z': xdim = ZDIM ; break ;
    }
    xdim *= sign ;
    nargs = 1 ;
    break ;
  case 'Y':
    cp = argv[2] ;
    if (*cp == '-')
    {
      cp++ ;  /* skip minus sign */
      sign = -1 ;
    }
    switch (toupper(*cp))
    {
    default:
    case 'X': ydim = XDIM ; break ;
    case 'Y': ydim = YDIM ; break ;
    case 'Z': ydim = ZDIM ; break ;
    }
    ydim *= sign ;
    nargs = 1 ;
    break ;
  case 'Z':
    cp = argv[2] ;
    if (*cp == '-')
    {
      cp++ ;  /* skip minus sign */
      sign = -1 ;
    }
    switch (toupper(*cp))
    {
    default:
    case 'X': zdim = XDIM ; break ;
    case 'Y': zdim = YDIM ; break ;
    case 'Z': zdim = ZDIM ; break ;
    }
    zdim *= sign ;
    nargs = 1 ;
    break ;
  case 'S':
    scale = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "scaling intensities by %2.2f\n", scale) ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s [input directory] [output directory]\n", argv[0]) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
MRI *
MRIreadBrueker(char *fname)
{
  MRI   *mri_tmp, *mri_short, *mri ;
  FILE  *fp ;
  int   nslices, imgsize = 256 * 256 * 2 ;  /* 16 bits per voxel */
  int   filesize, x, y, z, vox, pad ;
  float fmin, fmax, in_val, out_val ;

  if (!(fp = fopen(fname, "rb"))) 
    ErrorReturn(NULL, 
                (ERROR_NOFILE, "%s: could not read Brueker image.\n", fname)) ;

  fseek(fp, 0, SEEK_END) ;
  filesize = ftell(fp) ;
  fseek(fp, 0, SEEK_SET) ;

  nslices = filesize / imgsize ;
  fprintf(stderr, "reading %d slices...\n", nslices) ;
  mri_short = MRIalloc(256, 256, nslices, MRI_SHORT) ;

  for (z = 0 ; z < nslices ; z++)
  {
    for (y = 0 ; y < 256 ; y++)
    {
      for (x = 0 ; x < 256 ; x++)
      {
        fread2(&vox, fp) ;
        MRISvox(mri_short, x, y, z) = (short)vox ;
      }
    }
  }

  MRIvalRange(mri_short, &fmin, &fmax) ;
  mri_tmp = MRIalloc(256, 256, nslices, MRI_UCHAR) ;
  for (z = 0 ; z < nslices ; z++)
  {
    for (y = 0 ; y < 256 ; y++)
    {
      for (x = 0 ; x < 256 ; x++)
      {
        in_val = (float)MRISvox(mri_short, x, y, z) ;
        out_val = scale * (in_val - fmin) * 255.0 / (fmax-fmin) ;
        if (out_val > 255.0f)
          out_val = 255.0f ;
        MRIvox(mri_tmp, x, y, z) = (BUFTYPE)nint(out_val) ;
      }
    }
  }
  MRIfree(&mri_short) ;

  mri = MRIalloc(256, 256, 256, MRI_UCHAR) ;
  if (nslices != 256)
  {
    pad = (256-nslices)/2 ;
    MRIextractInto(mri_tmp, mri, 0, 0, 0, 256, 256, nslices, 0, 0, pad) ;
  }
  MRIfree(&mri_tmp) ;    
  fclose(fp);
  return(mri) ;
}
