
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "cma.h"
#include "colortab.h"
#include "gca.h"

static char vcid[] = "$Id: mri_edit_segmentation_with_surfaces.c,v 1.2 2003/04/09 20:28:22 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int relabel_hypointensities(MRI *mri, MRI_SURFACE *mris, int right) ;
static int relabel_hypointensities_neighboring_gray(MRI *mri) ;
static int edit_hippocampal_complex(MRI *mri, MRI_SURFACE *mris, int right, char *annot_name) ;
static int edit_hippocampus(MRI *mri) ;
int MRIneighbors(MRI *mri, int x0, int y0, int z0, int val) ;

static char *annot_name = "aparc.annot" ;


char *Progname ;

static char *label_name = NULL ;
static char *annotation_name = NULL ;

static char *surf_name = "white" ;

int
main(int argc, char *argv[])
{
  char          **av, *hemi, fname[STRLEN],
		            *in_aseg_name, *out_aseg_name, *surf_dir ;
  int           ac, nargs, h ;
  MRI_SURFACE   *mris ;
  MRI           *mri_aseg ;

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

  if (argc < 3)
    usage_exit() ;

	in_aseg_name = argv[1] ;
  surf_dir = argv[2] ;
	out_aseg_name = argv[3] ;

  mri_aseg = MRIread(in_aseg_name) ;
  if (!mri_aseg)
    ErrorExit(ERROR_NOFILE, "%s: could not read input segmentation %s", Progname, in_aseg_name) ;


	for (h = 0 ; h <= 1 ; h++)
	{
		if (h == 0)
			hemi = "lh" ;
		else
			hemi = "rh" ;
		sprintf(fname, "%s/%s.%s", surf_dir, hemi, surf_name)  ;
		printf("reading input surface %s...\n", fname) ;
		mris = MRISread(fname) ;
		if (!mris)
			ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
								Progname, fname) ;
		MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
		MRIScomputeMetricProperties(mris) ;
		MRISsmoothSurfaceNormals(mris, 10) ;   /* remove kinks in surface */

		printf("%s: editing hippocampal complex...\n", hemi) ;
		edit_hippocampal_complex(mri_aseg, mris, h, annot_name) ;
		printf("%s: relabeling hypointensities...\n", hemi) ;
		relabel_hypointensities(mri_aseg, mris, h) ;
		MRISfree(&mris) ;
	}
	relabel_hypointensities_neighboring_gray(mri_aseg) ;

	edit_hippocampus(mri_aseg) ;
	printf("writing modified segmentation to %s...\n", out_aseg_name) ;
	MRIwrite(mri_aseg, out_aseg_name) ;
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
  else if (!stricmp(option, "-annot"))
	{
		annot_name = argv[2] ;
		nargs = 2 ;
		printf("using annotation file %s...\n", annot_name) ;
	}
  else if (!stricmp(option, "debug_voxel"))
	{
		Gx = atoi(argv[2]) ;
		Gy = atoi(argv[3]) ;
		Gz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
	}
  else switch (toupper(*option))
  {
  case 'L':
    label_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "limiting computations to label %s.\n", label_name) ;
    break ;
  case 'A':
    annotation_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "computing statistics for each annotation in %s.\n", 
            annotation_name) ;
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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <subject name> <hemi> [<surface name>]\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program measures a variety of anatomical properties\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "-l <label file>              - limit calculations to specified "
          "label\n") ;
  fprintf(stderr,
          "-a <annotation file>         - compute properties for each label\n"
          "                               in the annotation file separately"
          "\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
relabel_hypointensities(MRI *mri, MRI_SURFACE *mris, int right)
{
	int   x, y, z, label, changed ;
	MRIS_HASH_TABLE *mht ;
	VERTEX           *v ;
	float            dx, dy, dz, dot, dist ;
	Real             xw, yw, zw ;

	mht = MHTfillVertexTableRes(mris,NULL, CURRENT_VERTICES, 8.0f) ;
	for (changed = x = 0 ; x < mri->width ; x++)
	{
		for (y = 0 ; y < mri->height ; y++)
		{
			for (z = 0 ; z < mri->depth ; z++)
			{
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				label = MRIvox(mri, x, y, z) ;
				if (label == Left_WM_hypointensities)
					MRIvox(mri, x, y, z) = WM_hypointensities ;
				else if (label == Right_WM_hypointensities)
					MRIvox(mri, x, y, z) = WM_hypointensities ;
				if ((!right && (label != Left_Cerebral_Cortex)) ||
						(right && (label != Right_Cerebral_Cortex)))
					continue ;
				if (MRIneighbors(mri, x, y, z, Unknown) >= 3) /* avoid stuff outside of brain in a sea of unknown */
					continue ;

				MRIvoxelToWorld(mri, x, y, z, &xw, &yw, &zw) ;
				v = MHTfindClosestVertexInTable(mht, mris, xw, yw, zw) ;
				if (v == NULL)  /* no vertices within range - assume it is hypointensity */
				{
					dot = -1 ;
					dist = 1000 ;
				}
				else
				{
					dx = xw - v->x ; dy = yw - v->y ; dz = zw - v->z ; 
					dot = v->nx*dx + v->ny*dy + v->nz*dz ;
					dist = sqrt(dx*dx+dy*dy+dz*dz) ;
				}
				if (dot < 0 && dist > 1)
				{
					changed++ ;
					MRIvox(mri, x, y, z) = WM_hypointensities ;
				}
			}
		}
	}

	printf("%d voxels changed to hypointensity...\n", changed) ;
	MHTfree(&mht) ;
	return(NO_ERROR) ;
}
int
relabel_hypointensities_neighboring_gray(MRI *mri)
{
	int    x, y, z, label, changed, i ;
	MRI    *mri_tmp = NULL ;

	for (changed = i = 0 ; i < 2 ; i++)
	{
		mri_tmp = MRIcopy(mri, mri_tmp) ;
		for (x = 0 ; x < mri->width ; x++)
		{
			for (y = 0 ; y < mri->height ; y++)
			{
				for (z = 0 ; z < mri->depth ; z++)
				{
					label = MRIvox(mri_tmp, x, y, z) ;
					if (label != WM_hypointensities)
						continue ;
					if (MRIneighbors(mri_tmp, x, y, z, Left_Cerebral_Cortex) > 0)
					{
						MRIvox(mri, x, y, z) = Left_Cerebral_Cortex ;
						changed++ ;
					}
					else  if (MRIneighbors(mri_tmp, x, y, z, Right_Cerebral_Cortex) > 0)
					{
						MRIvox(mri, x, y, z) = Right_Cerebral_Cortex ;
						changed++ ;
					}
				}
			}
		}
	}
	printf("%d hypointense voxels neighboring cortex changed\n", changed) ;
	return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIneighbors(MRI *mri, int x0, int y0, int z0, int val)
{
  int   nbrs = 0 ;

  if (MRIvox(mri,mri->xi[x0-1],y0,z0) == val)
    nbrs++ ;
  if (MRIvox(mri,mri->xi[x0+1],y0,z0) == val)
    nbrs++ ;
  if (MRIvox(mri,x0,mri->yi[y0+1],z0) == val)
    nbrs++ ;
  if (MRIvox(mri,x0,mri->yi[y0-1],z0) == val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,mri->zi[z0+1]) == val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,mri->zi[z0-1]) == val)
    nbrs++ ;
  return(nbrs) ;
}
static int
edit_hippocampal_complex(MRI *mri, MRI_SURFACE *mris, int right, char *annot_name)
{
	MRIS_HASH_TABLE *mht ;
	int              x, y, z, label, changed, index ;
	Real             xw, yw, zw ;
	char             fname[STRLEN] ;
	VERTEX           *v ;
	float            dx, dy, dz, dot, dist ;
	MRI              *mri_changed ;

	mri_changed = MRIclone(mri, NULL) ;

	if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
		ErrorExit(ERROR_NOFILE, "%s: could not read annotation file %s for hemi %s\n",
							annot_name, right ? "rh" : "lh") ;

	if (mris->ct == NULL)  /* color table not in annotation file */
	{
		char *cp ;
		cp = getenv("CSURF_DIR") ;
		sprintf(fname, "%s/Simple_surface_labels2002.txt", cp) ;
		printf("reading colortable from %s...\n", fname) ;
		mris->ct = CTABread(fname) ;
		if (!mris->ct)
			ErrorExit(ERROR_NOFILE, "%s: could not read color table from %s",Progname, fname) ;
	}
	mht = MHTfillVertexTableRes(mris,NULL, CURRENT_VERTICES, 8.0f) ;
	for (changed = x = 0 ; x < mri->width ; x++)
	{
		for (y = 0 ; y < mri->height ; y++)
		{
			for (z = 0 ; z < mri->depth ; z++)
			{
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				label = MRIvox(mri, x, y, z) ;
				/* only process amygdala and hippocampal voxels for the correct hemisphere */
				if ((right && ((label != Right_Hippocampus) && (label != Right_Amygdala) && (label != Right_Inf_Lat_Vent))) ||
						(!right && ((label != Left_Hippocampus) && (label != Left_Amygdala) && (label != Left_Inf_Lat_Vent))))
					continue ;
				MRIvoxelToWorld(mri, x, y, z, &xw, &yw, &zw) ;
				v = MHTfindClosestVertexInTable(mht, mris, xw, yw, zw) ;
				if (v == NULL)
					continue ;
				index = CTABannotationToIndex(mris->ct, v->annotation) ;
#define MEDIAL_WALL 41  /* only for Simple_surface_labels2002.txt - will have to modify for CMA */
#define LINGUAL_SULCUS 18
#define LINGUAL_SULCUS2 66
#define PARRAHIPPOCAMPAL_GYRUS 19
				if (index == MEDIAL_WALL || index  < 0)
					continue ;
				if (index == PARRAHIPPOCAMPAL_GYRUS)
				{
#define MIN_DIST 0.5
#if 1
					/* don't change voxels where the wm surface normal is pointing superiorly */
					if (v->nz >= 0)
					{
						dist = sqrt(SQR(v->x-xw)+SQR(v->y-yw)+SQR(v->z-zw)) ;
						if (dist > 2)  /* don't process amygdala that is far from surface */
							continue ;
						if (((fabs(v->nz) > fabs(v->nx)) || (fabs(v->nz) > fabs(v->ny))))
							continue ;   /* this voxel is superior to wm surface */
						/* at this point, the wm vertex is running nearly inferior-superior, so use
							 laterality to determine which side of parahippo it is on */
						if (fabs(xw) - fabs(v->x) < MIN_DIST)  /* if it isn't clearly lateral to parahippo, don't change it */
							continue ;
						
					}
					else if (((fabs(v->nz) < fabs(v->nx)) || (fabs(v->nz) < fabs(v->ny)))) /* not clearly inf or sup */
					{
						/* at this point, the wm vertex is running nearly inferior-superior, so use
							 laterality to determine which side of parahippo it is on */
						if (fabs(xw) - fabs(v->x) < MIN_DIST)  /* if it isn't clearly lateral to parahippo, don't change it */
							continue ;
						
					}
#else
					if (zw > (v->z) && (fabs(zw-v->z) > MIN_DIST))
						continue ;  /* don't change wm that is superior to parahippocampal gyrus */
					if (fabs(zw-v->z) < MIN_DIST)  /* too close - make sure it is medial of wm */
					{
						if (fabs(xw) - fabs(v->x) < MIN_DIST)  /* if it isn't clearly lateral to parahippo, don't change it */
							continue ;
					}
#endif
				}
				if (x == Gx && y == Gy && z == Gz)
					printf("voxel (%d, %d, %d): label %d (%s), parc %d, changing to cortex...\n",
								 Gx, Gy, Gz, label, cma_label_to_name(label), index) ;
				if (right)
					MRIvox(mri, x, y, z) = Right_Cerebral_Cortex ;
				else
					MRIvox(mri, x, y, z) = Left_Cerebral_Cortex ;
				changed++ ;
				MRIvox(mri_changed, x, y, z) = 1 ;
			}
		}
	}

 	for (x = 0 ; x < mri->width ; x++)
	{
		for (y = 0 ; y < mri->height ; y++)
		{
			for (z = 0 ; z < mri->depth ; z++)
			{
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				label = MRIvox(mri, x, y, z) ;
				/* only process amygdala and hippocampal voxels for the correct hemisphere */
				if ((right && ((label != Right_Hippocampus) && (label != Right_Amygdala) && (label != Right_Inf_Lat_Vent))) ||
						(!right && ((label != Left_Hippocampus) && (label != Left_Amygdala) && (label != Left_Inf_Lat_Vent))))
					continue ;
				MRIvoxelToWorld(mri, x, y, z, &xw, &yw, &zw) ;
				v = MHTfindClosestVertexInTable(mht, mris, xw, yw, zw) ;
				if (v == NULL)
					continue ;
				index = CTABannotationToIndex(mris->ct, v->annotation) ;
#define MEDIAL_WALL 41  /* only for Simple_surface_labels2002.txt - will have to modify for CMA */
#define LINGUAL_SULCUS 18
#define LINGUAL_SULCUS2 66
#define PARRAHIPPOCAMPAL_GYRUS 19
				if (index == MEDIAL_WALL || index  < 0)
					continue ;
				if ((index == PARRAHIPPOCAMPAL_GYRUS) && ((v->nz > 0) || (MRIneighborsOn(mri_changed, x, y, z,1) == 0)))
					continue ;  /* don't change wm that is superior to parahippocampal gyrus */
				if (x == Gx && y == Gy && z == Gz)
					printf("voxel (%d, %d, %d): label %d (%s), parc %d, changing to cortex...\n",
								 Gx, Gy, Gz, label, cma_label_to_name(label), index) ;
				if (right)
					MRIvox(mri, x, y, z) = Right_Cerebral_Cortex ;
				else
					MRIvox(mri, x, y, z) = Left_Cerebral_Cortex ;
				changed++ ;
				/*				MRIvox(mri_changed, x, y, z) = 1 ;*/
			}
		}
	}

	/* look for points inside wm labeled as cortex  and change them to wm */
	for (x = 0 ; x < mri->width ; x++)
	{
		for (y = 0 ; y < mri->height ; y++)
		{
			for (z = 0 ; z < mri->depth ; z++)
			{
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				label = MRIvox(mri, x, y, z) ;

				/* only process cortex voxels for the correct hemisphere */
				if ((right && ((label != Right_Cerebral_Cortex) && (label != Right_Hippocampus))) ||
						(!right && ((label != Left_Cerebral_Cortex) && (label != Left_Hippocampus))))
					continue ;
				MRIvoxelToWorld(mri, x, y, z, &xw, &yw, &zw) ;
				v = MHTfindClosestVertexInTable(mht, mris, xw, yw, zw) ;
				if (v == NULL)
					continue ;
				index = CTABannotationToIndex(mris->ct, v->annotation) ;
				/* only thin temporal wm */
				if ((index != MEDIAL_WALL) && (index != LINGUAL_SULCUS) && (index != LINGUAL_SULCUS2) &&
						(index != PARRAHIPPOCAMPAL_GYRUS))
					continue ;
				dx = xw - v->x ; dy = yw - v->y ; dz = zw - v->z ; 
				dot = v->nx*dx + v->ny*dy + v->nz*dz ;
				dist = sqrt(dx*dx+dy*dy+dz*dz) ;
				if (dot < 0 && dist < 2)
				{
					if (x == Gx && y == Gy && z == Gz)
						printf("voxel (%d, %d, %d): label %d (%s), parc %d, changing to wm...\n",
									 Gx, Gy, Gz, label, cma_label_to_name(label), index) ;
					changed++ ;
					MRIvox(mri, x, y, z) = right ? Right_Cerebral_White_Matter : Left_Cerebral_White_Matter ;
				}
			}
		}
	}
	printf("%d voxels changed in medial temporal lobe...\n", changed) ;
	MHTfree(&mht) ; MRIfree(&mri_changed) ;
	return(NO_ERROR) ;
}


static int
edit_hippocampus(MRI *mri)
{
	int   x, y, z, right, changed, total = 0, yi, yk, found, label ;

	do
	{
		changed =  0 ;
		for (x = 0 ; x < mri->width ; x++)
		{
			for (y = 0 ; y < mri->height ; y++)
			{
				for (z = 0 ; z < mri->depth ; z++)
				{
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak() ;
					label = MRIvox(mri, x, y, z) ;
					if (IS_CORTEX(label) == 0)
						continue ;
					right = (label == Right_Cerebral_Cortex) ;
					yi = mri->yi[y-1] ;  /* voxel immediately superior */
					if (IS_HIPPO(MRIvox(mri, x, yi, z)) == 0)
						continue ;

					/* check for wm within 3 mm inferior */
					for (found = 0, yk = 1 ; yk <= 3 ; yk++)
					{
						yi = mri->yi[y+yk] ;  /* inferior voxel */
						if (IS_WM(MRIvox(mri, x, yi, z)) != 0)
						{
							found = 1 ;
							break ;
						}
					}
					if (!found)
						continue ;
					
					if (x == Gx && y == Gy && z == Gz)
						printf("changing voxel (%d, %d, %d) from %s to ",Gx, Gy, Gz, cma_label_to_name(label)) ;
					changed++ ;
					MRIvox(mri, x, y, z) = right ? Right_Hippocampus : Left_Hippocampus ;
					if (x == Gx && y == Gy && z == Gz)
						printf("%s\n", cma_label_to_name(MRIvox(mri, x, y, z))) ;
				}
			}
		}
		total += changed ;
	} while (changed > 0) ;

	printf("%d cortex voxels changed to hippocampus...\n", total) ;
	return(NO_ERROR)  ;
}

