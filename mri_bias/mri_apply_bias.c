#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "mrinorm.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *apply_bias(MRI *mri_orig, MRI *mri_norm, MRI *mri_bias) ;

char *Progname ;

static void usage_exit(int code) ;

int
main(int argc, char *argv[])
{
  char         **av ;
  int          ac, nargs ;
  MRI          *mri_orig, *mri_norm, *mri_bias ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mri_apply_bias.c,v 1.1 2006/07/13 19:24:15 fischl Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

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

	mri_orig = MRIread(argv[1]) ;
	mri_bias = MRIread(argv[2]) ;

	mri_norm = apply_bias(mri_orig, NULL, mri_bias);

  fprintf(stderr, "writing to %s...\n", argv[3]) ;
  MRIwrite(mri_norm, argv[3]) ;
  MRIfree(&mri_norm) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "bias correction took %d minutes and %d seconds.\n",
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
  if (!stricmp(option, "sdir"))
	{
	}
  else switch (toupper(*option))
	{
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
  printf("usage: %s [options] <input vol> <bias vol> <output volume>\n", Progname) ;
  exit(code) ;
}
static MRI *
apply_bias(MRI *mri_orig, MRI *mri_norm, MRI *mri_bias)
{
	MATRIX   *m_vox2vox;
	VECTOR   *v1, *v2;
	int      x, y, z ;
	double   xd, yd, zd, bias, val_orig, val_norm ;

	if (mri_norm == NULL)
		mri_norm = MRIclone(mri_orig, NULL) ;

	m_vox2vox = MRIgetVoxelToVoxelXform(mri_orig, mri_bias) ;
	v1 = VectorAlloc(4, MATRIX_REAL);
	v2 = VectorAlloc(4, MATRIX_REAL);
	VECTOR_ELT(v1, 4) = 1.0 ;
	VECTOR_ELT(v2, 4) = 1.0 ;

	for (x = 0 ; x < mri_orig->width ; x++)
	{
		V3_X(v1) = x ;
		for (y = 0 ; y < mri_orig->height ; y++)
		{
			V3_Y(v1) = y ;
			for (z = 0 ; z < mri_orig->depth ; z++)
			{
				V3_Z(v1) = z ;
				val_orig = MRIgetVoxVal(mri_orig, x, y, z, 0) ;
				MatrixMultiply(m_vox2vox, v1, v2) ;
				xd = V3_X(v2) ; yd = V3_Y(v2) ;zd = V3_Z(v2);
				MRIsampleVolume(mri_bias, xd, yd, zd, &bias) ;
				val_norm = val_orig * bias ;
				if (mri_norm->type == MRI_UCHAR)
				{
					if (val_norm > 255)
						val_norm = 255 ;
					else if (val_norm < 0)
						val_norm = 0 ;
				}
				MRIsetVoxVal(mri_norm, x, y, z, 0, val_norm) ;
			}
		}
	}

	MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ;
	return(mri_norm) ;
}
