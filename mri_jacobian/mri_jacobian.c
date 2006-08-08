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
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
                

char *Progname ;

static int use_log = 0 ;
static float sigma = 0 ;
static int write_areas = 0 ;
static int init = 1 ;
static int atlas = 0 ;

static void usage_exit(int code) ;
static int find_debug_node(GCA_MORPH *gcam, int origx, int origy, int origz) ;
static int init_gcam_areas(GCA_MORPH *gcam) ;
int
main(int argc, char *argv[])
{
  char      **av, *out_fname ;
  int       ac, nargs ;
	GCA_MORPH *gcam ;
  int       msec, minutes, seconds ;
  struct timeb start ;
	MRI       *mri, *mri_jacobian, *mri_area, *mri_orig_area ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_jacobian.c,v 1.2 2006/08/08 13:12:45 fischl Exp $", "$Name:  $");
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

  if (argc < 4)
    usage_exit(1) ;

	out_fname = argv[argc-1] ;
	gcam = GCAMread(argv[1]) ;
	if (gcam == NULL)
		ErrorExit(ERROR_BADPARM, "%s: could not read input morph %s\n", Progname,argv[1]);
	if (Gx >= 0)
		find_debug_node(gcam, Gx, Gy, Gz) ;

	if (init)
		init_gcam_areas(gcam) ;
	mri = MRIread(argv[2]) ;
	if (gcam == NULL)
		ErrorExit(ERROR_BADPARM, "%s: could not read template volume %s\n", Progname,argv[2]);

	if (atlas)
	{
		mri_area = GCAMwriteMRI(gcam, NULL, GCAM_AREA);
		mri_orig_area = GCAMwriteMRI(gcam, NULL, GCAM_ORIG_AREA);
	}
	else
	{
		mri_area = GCAMmorphFieldFromAtlas(gcam, mri, GCAM_AREA, 0, 0);
		mri_orig_area = GCAMmorphFieldFromAtlas(gcam, mri, GCAM_ORIG_AREA, 0, 0);
	}
	if (FZERO(sigma) == 0)
	{
		MRI *mri_kernel, *mri_smooth ;
		mri_kernel = MRIgaussian1d(sigma, 100) ;
		mri_smooth = MRIconvolveGaussian(mri_area, NULL, mri_kernel) ;
		MRIfree(&mri_area) ; mri_area = mri_smooth ;
		mri_smooth = MRIconvolveGaussian(mri_orig_area, NULL, mri_kernel) ;
		MRIfree(&mri_orig_area) ; mri_orig_area = mri_smooth ;

		MRIfree(&mri_kernel) ; 
	}
	mri_jacobian = MRIdivide(mri_area, mri_orig_area, mri_jacobian) ;
	if (use_log)
		MRIlog10(mri_jacobian, mri_jacobian, 0) ;
  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_jacobian, out_fname) ;
  MRIfree(&mri_jacobian) ;
	if (write_areas)
	{
		MRIwrite(mri_area, "area.mgz") ;
		MRIwrite(mri_orig_area, "orig_area.mgz") ;
	}
	if (atlas && DIAG_WRITE && DIAG_VERBOSE_ON)
	{
		char fname[STRLEN] ;
		FileNameRemoveExtension(out_fname, out_fname) ;
		mri_area = GCAMwriteMRI(gcam, mri_area, GCAM_MEANS);
		sprintf(fname, "%s_intensity.mgz", out_fname) ;
		printf("writing means to %s\n", fname) ;
		MRIwrite(mri_area, fname) ;
		sprintf(fname, "%s_labels.mgz", out_fname) ;
		mri_area = GCAMwriteMRI(gcam, mri_area, GCAM_LABEL);
		printf("writing labels to %s\n", fname) ;
		MRIwrite(mri_area, fname) ;
	}
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ; minutes = seconds / 60 ;seconds = seconds % 60 ;
  fprintf(stderr, "jacobian calculation took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  }
	else if (!stricmp(option, "debug_voxel"))
	{
		Gx = atoi(argv[2]) ; Gy = atoi(argv[3]) ; Gz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
	}
  else switch (toupper(*option))
  {
	case 'A':
		atlas = 1 ;
		printf("outputing in atlas coords\n") ;
		break ;
	case 'W':
		write_areas = 1 ;
		printf("writing area volumes\n") ;
		break ;
	case 'L':
		use_log = 1 ;
		printf("taking log of jacobian values before saving\n") ;
		break ;
	case 'S':
		sigma = atof(argv[2]) ;
		printf("smoothing jacobian volume with sigma=%2.2f\n", sigma) ;
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
  printf("usage: %s [options] <3d morph> <template volume> <output volume>\n", Progname) ;
  exit(code) ;
}

static int
find_debug_node(GCA_MORPH *gcam, int origx, int origy, int origz)
{
	int            x, y, z, xmin, ymin, zmin ;
	double         d, dmin ;
	GCA_MORPH_NODE *gcamn ;


	dmin = 1e10 ; xmin = ymin = zmin = 0 ;
	for (x = 0 ; x < gcam->width ; x++)
	{
		for (y = 0 ; y < gcam->height ; y++)
		{	
			for (z = 0 ; z < gcam->depth ; z++)
			{
				gcamn = &gcam->nodes[x][y][z] ;
				d = sqrt(SQR(gcamn->origx-origx) + SQR(gcamn->origy-origy) + SQR(gcamn->origz-origz)) ;
				if (d < dmin)
				{
					dmin = d ; xmin = x ; ymin = y ; zmin = z ;
				}
			}
		}
	}
	gcamn = &gcam->nodes[xmin][ymin][zmin] ;
	printf("Talairach voxel (%d, %d, %d) maps to node (%d, %d, %d) %s --> (%2.1f, %2.1f, %2.1f)\n",
				 origx, origy, origz, xmin, ymin, zmin, cma_label_to_name(gcamn->label), gcamn->x, gcamn->y, gcamn->z) ;
		
	return(NO_ERROR) ;
}

static int
init_gcam_areas(GCA_MORPH *gcam)
{
	int            x, y, z ;
	double         orig_area ;
	GCA_MORPH_NODE *gcamn ;

	orig_area = gcam->spacing * gcam->spacing * gcam->spacing ;
	for (x = 0 ; x < gcam->width ; x++)
	{
		for (y = 0 ; y < gcam->height ; y++)
		{	
			for (z = 0 ; z < gcam->depth ; z++)
			{
				gcamn = &gcam->nodes[x][y][z] ;
				gcamn->orig_area = gcamn->orig_area1 = gcamn->orig_area2 = orig_area ;
			}
		}
	}
	return(NO_ERROR) ;
}

