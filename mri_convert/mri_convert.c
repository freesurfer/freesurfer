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
#include "mri_conform.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void usage(int exit_val);

char *Progname ;

static int verbose = 1 ;
static int conform = 0;
static int xdim = XDIM, ydim = YDIM, zdim = ZDIM ;
static int raw_flag;
static int raw_width, raw_height, raw_depth, raw_type;

static float blur_sigma = 0.0f ;
static float scale = 1.0 ;
static char *inverse_transform_fname = NULL ;
static char *transform_fname = NULL ;


int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  MRI    *mri, *mri2, *mri_kernel ;
  char   *in_fname, *out_fname ;
  FILE   *fin;

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
    fprintf(stderr, "reading from %s...\n", in_fname) ;
  if(raw_flag)
  {
    if((fin = fopen(in_fname, "r")) == NULL)
      ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s", 
              Progname, in_fname) ;
    mri = MRIreadRaw(fin, raw_width, raw_height, raw_depth, raw_type);
    fclose(fin);
  }
  else
    mri = MRIread(in_fname) ;

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

    if(verbose)
      printf("smoothing volume with sigma = %2.3f\n", blur_sigma) ;
    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_smooth = MRIclone(mri, NULL) ;
    MRIconvolveGaussian(mri, mri_smooth, mri_kernel) ;
    MRIfree(&mri_kernel) ; MRIfree(&mri) ; 
    mri = mri_smooth ;
  }

  if (conform)
  {
    if(verbose)
      printf("conforming volume...\n");
    mri2 = MRIconform(mri);
    mri = MRIcopy(mri2, NULL);
    MRIfree(&mri2);
  }

  if (verbose)
    fprintf(stderr, "writing to %s...\n", out_fname) ;

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

  if(verbose)
    printf("writing output to '%s'...\n", out_fname) ;
  MRIwrite(mri, out_fname) ;
  MRIfree(&mri) ;
  exit(0) ;
  return(0) ;
}

static void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s [options] input_file output_file\n", Progname) ;
  fprintf(fout, "  options are:\n");
  fprintf(fout, "  -u, -?          display usage and exit\n");
  fprintf(fout, "  -v              verbose\n");
  fprintf(fout, "  -conform        conform the volume to 256x256x256, uchar\n");
  fprintf(fout, "  -r x y z        reorder dimensions\n");
  fprintf(fout, "  -blur sigma     blur the volume\n");
  fprintf(fout, "  -raw x y z type read a raw data file; type is one of:\n");
  fprintf(fout, "                  uchar, int, long, float, short\n");
  fprintf(fout, "  -x xdim -y ydim -z zdim\n");
  fprintf(fout, "                  reorder the axes\n");
  fprintf(fout, "  -T transform    apply a transform\n");

  exit(exit_val) ;

} /* end usage() */

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0, sign = 1 ;
  char *option, *cp ;
  int i;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "blur"))
  {
    blur_sigma = atof(argv[2]) ;
    fprintf(stderr, "applying %2.2f standard deviation Gaussian kernel\n",
            blur_sigma) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "conform"))
    conform = 1 ;
  else if (!stricmp(option, "raw"))
  {
    for(i = 2;i <= 5;i++)
    {
      if(!argv[i])
        usage(1);
    }

    raw_width = atoi(argv[2]);
    raw_height = atoi(argv[3]);
    raw_depth = atoi(argv[4]);

    if(!stricmp(argv[5], "uchar"))
      raw_type = MRI_UCHAR;
    if(!stricmp(argv[5], "int"))
      raw_type = MRI_INT;
    if(!stricmp(argv[5], "long"))
      raw_type = MRI_LONG;
    if(!stricmp(argv[5], "float"))
      raw_type = MRI_FLOAT;
    if(!stricmp(argv[5], "short"))
      raw_type = MRI_SHORT;
    else
      usage(1);

    raw_flag = 1;
    nargs = 4;
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

/* eof */
