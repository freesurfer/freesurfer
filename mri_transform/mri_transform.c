#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "transform.h"

static char vcid[] = "$Id: mri_transform.c,v 1.2 2002/04/22 21:59:00 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static char *out_like_fname = NULL ;
static int invert_flag = 0 ;

int
main(int argc, char *argv[])
{
  char        **av, *in_vol, *out_vol, *xform_fname ;
  int         ac, nargs, i, ras_flag = 0 ;
  MRI         *mri_in, *mri_out ;
  LTA         *lta ;
  MATRIX      *m, *m_total ;
  TRANSFORM   *transform ;

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
    usage_exit() ;

  in_vol = argv[1] ;
  out_vol = argv[argc-1] ;

  printf("reading volume from %s...\n", in_vol) ;
  mri_in = MRIread(in_vol) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, 
              in_vol) ;
  if (out_like_fname)
  {
    MRI *mri_tmp = MRIread(out_like_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read template volume from %s",out_like_fname) ;
    mri_out = MRIalloc(mri_tmp->width, mri_tmp->height, mri_tmp->depth, mri_tmp->type) ;
    /*    MRIcopyHeader(mri_tmp, mri_out) ;*/
    MRIfree(&mri_tmp) ;
  }
  else
    mri_out = MRIalloc(256, 256, 256, mri_in->type) ;

  MRIcopyPulseParameters(mri_in, mri_out) ;
  m_total = MatrixIdentity(4, NULL) ;
  for (i = 2 ; i < argc-1 ; i++)
  {
    xform_fname = argv[i] ;
    if (strcmp(xform_fname, "-I") == 0)
    {
      invert_flag = 1 ;
      continue ;
    }
    printf("reading transform %s...\n", xform_fname) ;
    transform = TransformRead(xform_fname) ;
    if (!transform)
      ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", 
                Progname, xform_fname) ;

    if (transform->type != MORPH_3D_TYPE)
    {
      lta = (LTA *)(transform->xform) ;
      m = MatrixCopy(lta->xforms[0].m_L, NULL) ;
      
      if (lta->type == LINEAR_RAS_TO_RAS)  /* convert it to a voxel transform */
      {
        ras_flag = 1 ;
      }
      else if (ras_flag)
        ErrorExit(ERROR_UNSUPPORTED, "%s: transforms must be all RAS or all voxel",Progname) ;
      
      if (invert_flag)
      {
        MATRIX *m_tmp ;
        printf("inverting transform...\n") ;
        m_tmp = MatrixInverse(m, NULL) ;
        if (!m_tmp)
          ErrorExit(ERROR_BADPARM, "%s: transform is singular!") ;
        MatrixFree(&m) ; m = m_tmp ;
        invert_flag = 0 ;
      }
      MatrixMultiply(m, m_total, m_total) ;
      LTAfree(&lta) ;
      if (ras_flag)  /* convert it to a voxel transform */
      {
        MATRIX *m_tmp ;
        printf("converting RAS xform to voxel xform...\n") ;
        m_tmp = MRIrasXformToVoxelXform(mri_in, mri_out, m_total, NULL) ;
        MatrixFree(&m_total) ; m_total = m_tmp ;
      }
      MRIlinearTransform(mri_in, mri_out, m_total) ;
    }
    else
    {
      if (invert_flag)
        mri_out = TransformApplyInverse(transform, mri_in, NULL) ;
      else
        mri_out = TransformApply(transform, mri_in, NULL) ;
    }
    invert_flag = 0 ;
  }
  
  printf("writing output to %s.\n", out_vol) ;
  MRIwrite(mri_out, out_vol) ;

  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "out_like") || !stricmp(option, "ol"))
  {
    out_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
  }
  else switch (toupper(*option))
  {
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'I':
    invert_flag = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input volume> <input surface> <registration file> <output .float file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
     "\nThis program will paint a average Talairach stats onto a surface\n");
  fprintf(stderr, "-imageoffset <image offset> - set offset to use\n") ;
  fprintf(stderr, "-S                          - paint using surface "
          "coordinates\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

