//
// mri_hires_register.c
//
// written by Bruce Fischl
// Nov. 9th ,2000
// 
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2005/01/31 21:16:25 $
// Revision       : $Revision: 1.5 $
//
////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gcamorph.h"

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "version.h"
#include "transform.h"
#include "nr.h"

typedef struct
{
	int *xi ;
	int *yi ;
	int *zi ;
	MRI *mri ;
	int nvox ;
} VOXEL_LIST ;

#define MAX_ANGLE       RADIANS(25)

#define PAD             2

static MRI *mri_from_voxel_list(VOXEL_LIST *vl, MRI *mri) ;
static int  free_voxel_list(VOXEL_LIST **pvoxel_list) ;
static VOXEL_LIST *create_voxel_list(MRI *mri, float low_val, float hi_val , 
																		 VOXEL_LIST *vl, int skip) ;
static int write_snapshot(MRI *mri_lowres, MRI *mri_hires, 
													MATRIX *m_vox_xform, INTEGRATION_PARMS *parms, 
													int fno, int conform, char *fname) ;

static double MAX_TRANS = 30 ;

static float compute_powell_sse(float *p) ;
static int    powell_minimize(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, MATRIX *mat);

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;

static TRANSFORM *compute_optimal_transform(VOXEL_LIST *vl_lowres, 
																						VOXEL_LIST *vl_hires, 
																						INTEGRATION_PARMS *parms,
																						TRANSFORM *transform) ;
static double compute_overlap(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, MATRIX *m_L) ;
static double find_optimal_translation(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, 
																			 MATRIX *m_L, float min_trans, 
																			 float max_trans, 
																			 float trans_steps, int nreductions);

static double find_optimal_linear_xform(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires,
																				MATRIX *m_L, MATRIX *m_origin,
																				float min_angle, float max_angle,
																				float min_scale, float max_scale,
																				float min_trans, float max_trans,
																				float angle_steps, float scale_steps, 
																				float trans_steps,
																				int nreductions);
char *Progname ;
static int target_label = Left_Hippocampus ;

static int skip = 2 ;

#define nothing                    0
#define alveus                     1
#define perforant_pathway	         2
#define parasubiculum		           3
#define presubiculum			         4
#define subiculum			             5
#define CA1				                 6
#define CA2				                 7
#define CA3				                 8
#define CA4				                 9
#define GC_DG			                 10
#define HATA				               11
#define fimbria                    12
#define lateral_ventricle		       13
#define molecular_layer_HP         14
#define molecular_layer_subiculum	 15

static MRI *mri_hires_intensity = NULL ;
static char *hires_intensity_fname = NULL ;

#if 0
static int gray_labels[] =
	{
		parasubiculum, presubiculum, subiculum, CA1, CA2, CA3, CA4, HATA,
		molecular_layer_subiculum, molecular_layer_HP, GC_DG
	} ;

static int white_labels[] =
	{
		alveus, perforant_pathway, fimbria
	} ;
#endif

static INTEGRATION_PARMS parms ;
static TRANSFORM  *transform = NULL ;
static GCA_MORPH_PARMS mp ;

int
main(int argc, char *argv[])
{
	char       **av, *hires_fname, *aseg_fname, *intensity_fname, *out_fname, fname[STRLEN] ;
  int        ac, nargs, i ;
	MRI        *mri_intensity, *mri_lowres, *mri_hires, *mri_tmp ;
	VOXEL_LIST *vl_lowres, *vl_hires ;
	MRI_REGION  box ;
  struct timeb start ;
  int          msec, minutes, seconds ;

	parms.write_iterations = 0 ;
	parms.start_t = 0 ;

  TimerStart(&start) ;
  setRandomSeed(-1L) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

	Progname = argv[0] ;
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

	hires_fname = argv[1] ;
	intensity_fname = argv[2] ;
	aseg_fname = argv[3] ;
	out_fname = argv[4] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
	mri_hires = MRIread(hires_fname) ;
	if (!mri_hires)
		ErrorExit(ERROR_NOFILE, "%s: could not read hires label volume %s",
							Progname, hires_fname) ;
	if (mri_hires->type != MRI_UCHAR)
	{
		mri_tmp = MRIchangeType(mri_hires, MRI_UCHAR, 0, 255, 1) ;
		MRIfree(&mri_hires) ;
		mri_hires = mri_tmp ;
	}

	MRIboundingBox(mri_hires, 0, &box) ;
	box.x -= PAD ; box.y -= PAD ; box.z -= PAD ; 
	box.dx += 2*PAD ; box.dy += 2*PAD ; box.dz += 2*PAD ; 
	mri_tmp = MRIextractRegion(mri_hires, NULL, &box) ;
	MRIfree(&mri_hires) ;
	mri_hires = mri_tmp ;

	mri_intensity = MRIread(intensity_fname) ;
	if (!mri_intensity)
		ErrorExit(ERROR_NOFILE, "%s: could not read intensity label volume %s",
							Progname, intensity_fname) ;
	mri_lowres = MRIread(aseg_fname) ;
	if (!mri_lowres)
		ErrorExit(ERROR_NOFILE, "%s: could not read aseg label volume %s",
							Progname, aseg_fname) ;
	mri_tmp = MRIclone(mri_lowres, NULL) ; MRIcopyLabel(mri_lowres, mri_tmp, target_label) ;
	MRIfree(&mri_lowres) ; mri_lowres = mri_tmp ;

  if (Gdiag & DIAG_WRITE && parms.write_iterations > 0)
  {
		sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_lowres, fname, IMAGE_SIZE) ;	
		sprintf(fname, "intensity_%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_intensity, fname, IMAGE_SIZE) ;
		MRIwrite(mri_intensity, "intensity_target.mgz") ;
		MRIwrite(mri_lowres, "aseg_target.mgz") ;
	}

	if (hires_intensity_fname)
	{
		MRI *mri_tmp2 = MRIread(hires_intensity_fname) ;
		if (!mri_tmp2)
			ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s...\n", hires_intensity_fname) ;
		if (mri_tmp2->type != MRI_UCHAR && mri_tmp2->type != MRI_FLOAT)
		{
			mri_hires_intensity = MRIchangeType(mri_tmp2, MRI_FLOAT, 0, 2000, 1) ;
			MRIfree(&mri_tmp2) ;
			mri_tmp2 = mri_hires_intensity ;
		}
	}

	if (transform == NULL)   /* compute optimal linear transform */
	{
		vl_lowres = create_voxel_list(mri_lowres,target_label,target_label,NULL,0);
		for (i = 0 ; i < 3 ; i++)
		{
			printf("------------- outer loop iteration %d ---------------\n",i) ;
			vl_hires = create_voxel_list(mri_hires, 1, 255, NULL, skip) ;
			
			transform = compute_optimal_transform(vl_lowres, vl_hires, &parms,
																						transform) ;
			free_voxel_list(&vl_hires) ;
			vl_hires = create_voxel_list(mri_hires, 1, 255, NULL, skip/4) ;
			powell_minimize(vl_lowres, vl_hires, ((LTA *)(transform->xform))->xforms[0].m_L) ;
			free_voxel_list(&vl_hires) ;
		}
		free_voxel_list(&vl_lowres) ;

		printf("final vox2vox matrix:\n") ;
		MatrixPrint(stdout, ((LTA *)(transform->xform))->xforms[0].m_L) ;
		{
			MRI *mri_aligned, *mri_filtered ;
			char fname[STRLEN] ;
			int  i ;

			mri_aligned = MRIsrcTransformedCentered(mri_hires, mri_lowres, ((LTA *)(transform->xform))->xforms[0].m_L, SAMPLE_NEAREST) ;
			sprintf(fname, "%sfinal.mgz", parms.base_name) ;
			MRIwrite(mri_aligned, fname) ;

			for (i = 1 ; i <= 10 ; i++)
			{
				sprintf(fname, "%sfiltered%d.mgz", parms.base_name, i) ;
				mri_filtered = MRImodeFilter(mri_aligned, NULL, i) ;
				printf("writing filtered image to %s\n", fname) ;
				MRIwrite(mri_filtered, fname) ;
				MRIfree(&mri_filtered) ;
			}
			MRIfree(&mri_aligned) ;
		}

		LTAvoxelToRasXform((LTA *)(transform->xform), mri_hires, mri_lowres) ;
		TransformWrite(transform, out_fname) ;
	}
	else   /* linear transform already computed - compute 3d morph */
	{
		GCA_MORPH *gcam ;
		MATRIX    *m_L ;
		LTA       *lta ;

		lta = ((LTA *)(transform->xform)) ;
    m_L = MRIrasXformToVoxelXform(mri_hires, mri_lowres, lta->xforms[0].m_L, NULL) ;
    MatrixFree(&lta->xforms[0].m_L) ;

		if (MatrixDeterminant(m_L) < 0)
		{
			MATRIX *m_I ;

			printf("transform det < 0 --> replacing with identity and transforming hires label volume\n");
			mri_tmp = MRIsrcTransformedCentered(mri_hires, mri_lowres, m_L, SAMPLE_NEAREST) ;
			MRIfree(&mri_hires) ;
			mri_hires = mri_tmp ;
			m_I = MatrixIdentity(4, NULL) ;
			MRIrasXformToVoxelXform(mri_hires, mri_lowres, m_I, m_L);
			MatrixFree(&m_I) ;
			
			/* make sure none of the labels are on the border */
			mri_tmp = MRIalloc(mri_hires->width+2, mri_hires->height+2, mri_hires->depth+2, MRI_UCHAR) ;
			MRIextractInto(mri_hires, mri_tmp, 0, 0, 0, mri_hires->width, mri_hires->width, mri_hires->depth, 1, 1, 1);
			MRIfree(&mri_hires) ; mri_hires = mri_tmp ;
			MRIwrite(mri_hires, "hires_xformed.mgz") ;
		}

    lta->xforms[0].m_L = m_L ;
		printf("initializing GCAM with vox->vox matrix:\n") ;
		MatrixPrint(stdout, m_L) ;
		memset(&mp, 0, sizeof(mp)) ;
		strcpy(mp.base_name, parms.base_name) ;
		mp.write_iterations = parms.write_iterations ;
		mp.l_jacobian = 1 ;
		mp.l_binary = 1 ;
		mp.dt = 0.01 ;
		mp.noneg = True ;
		mp.exp_k = 20 ;
		mp.max_grad = 0.3 ;
		mp.momentum = 0.0 ;
		mp.l_smoothness = 1 ;
		mp.sigma = 8 ;
		mp.relabel_avgs = -1 ;
		mp.navgs = 64 ;
		mp.levels = 4 ;
		mp.integration_type = GCAM_INTEGRATE_BOTH ;
		mp.nsmall = 1 ;
		mp.reset_avgs = -1 ;
		mp.tol = 0.0001 ;
		mp.niterations = 1000 ;
		gcam = GCAMalloc(mri_hires->width, mri_hires->height, mri_hires->depth);
    GCAMinit(gcam, mri_lowres, NULL, transform, 0) ;
		GCAMinitLabels(gcam, mri_hires) ;
    GCAMsetStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ; /* disable everything */
		GCAMregister(gcam, mri_lowres, &mp) ;
		GCAMwrite(gcam, out_fname) ;
	}

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("registration took %d minutes and %d seconds.\n", 
				 minutes, seconds) ;
	exit(0) ;
	return(0) ;
}

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!stricmp(option, "debug_voxel"))
	{
		Gx = atoi(argv[2]) ;
		Gy = atoi(argv[3]) ;
		Gz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
	}
  else if (!stricmp(option, "view"))
	{
		Gsx = atoi(argv[2]) ;
		Gsy = atoi(argv[3]) ;
		Gsz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("viewing voxel (%d, %d, %d)\n", Gsx, Gsy, Gsz) ;
	}
	else if (!stricmp(option, "dt"))
	{
		mp.dt = atof(argv[2]) ;
		nargs = 1 ;
		printf("using dt = %2.3f\n", mp.dt) ;
	}
	else switch (*option)
	{
  case 'M':
    mp.momentum = atof(argv[2]) ;
    nargs = 1 ;
    printf("momentum = %2.2f\n", mp.momentum) ;
    break ;
	case 'S':
		skip = atoi(argv[2]);
		printf("skipping %d voxels in hires data...\n", skip) ;
		nargs = 1 ;
		break ;
	case 'T':
		printf("reading transform from %s...\n", argv[2]) ;
		transform = TransformRead(argv[2]) ;
		if (transform == NULL)
			ErrorExit(ERROR_NOFILE,"%s: could not read transform from %s\n",argv[2]);
		nargs = 1 ;
		break ;
	case 'I':
		hires_intensity_fname = argv[2] ;
		nargs = 1 ;
		printf("reading intensity image from %s for debugging...\n", hires_intensity_fname) ;
		break ;
	case 'W':
		mp.write_iterations = parms.write_iterations = atoi(argv[2]) ;
		Gdiag |= DIAG_WRITE ;
		nargs = 1 ;
		printf("setting write iterations = %d\n", parms.write_iterations) ;
		break ;
  case '?':
  case 'U':
    usage_exit(1);
    break ;
	default:
    printf("unknown option %s\n", argv[1]) ;
		usage_exit(1) ;
    break ;
	}
	return(nargs) ;
}

static void 
usage_exit(int ecode)
{
	printf(
				 "usage: %s <hires labeling> <input intensity> <input aseg>"
				 " <output xform>\n",
				 Progname) ;
	exit(ecode) ;
}


static TRANSFORM *
compute_optimal_transform(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, 
													INTEGRATION_PARMS *parms, TRANSFORM *transform)
{
	MATRIX    *m_vox_xform, *m_origin, *m_inv_origin, *m_hires_vox2ras, 
		        *m_lowres_ras2vox, *m_trans, *m_tmp ;
	VECTOR    *v_cl, *v_ch ;
	double    hires_cent[3], lowres_cent[3], dx, dy, dz, scale,min_search_scale ;
	int       niter, nscales, good_step, done, trans ;
	double    old_max_overlap, max_overlap ;
	MRI       *mri_lowres, *mri_hires ;

	mri_lowres = vl_lowres->mri ; mri_hires = vl_hires->mri ;

#define MIN_SEARCH_SCALE 0.01
  min_search_scale = MIN_SEARCH_SCALE ;
	m_origin = MatrixIdentity(4, NULL) ;
	MRIcenterOfMass(mri_hires, hires_cent, 0) ;
	MRIcenterOfMass(mri_lowres, lowres_cent, 0) ;
	*MATRIX_RELT(m_origin, 1, 4) = lowres_cent[0] ; 
	*MATRIX_RELT(m_origin, 2, 4) = lowres_cent[1] ;
	*MATRIX_RELT(m_origin, 3, 4) = lowres_cent[2] ; 
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;
	m_inv_origin = MatrixInverse(m_origin, NULL) ;
	if (transform == NULL)
	{
		transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
		m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;
		
		m_lowres_ras2vox = MRIgetRasToVoxelXform(mri_lowres) ;
		m_hires_vox2ras = MRIgetVoxelToRasXform(mri_hires) ;
		MatrixMultiply(m_lowres_ras2vox, m_hires_vox2ras, m_vox_xform) ;
		printf("initial transform from direction cosines:\n") ;
		MatrixPrint(stdout, m_vox_xform) ;
		MatrixFree(&m_lowres_ras2vox) ; MatrixFree(&m_hires_vox2ras) ;
		
		dx = lowres_cent[0] - hires_cent[0] ; dy = lowres_cent[1] - hires_cent[1]  ;
		dz = lowres_cent[2] - hires_cent[2] ;
		
		v_cl = VectorAlloc(4, MATRIX_REAL) ; v_ch = VectorAlloc(4, MATRIX_REAL) ; 
		*MATRIX_RELT(v_cl,4,1) = 1.0 ; *MATRIX_RELT(v_ch,4,1) = 1.0 ;
		V3_X(v_ch) = hires_cent[0] ; V3_Y(v_ch) = hires_cent[1] ; V3_Z(v_ch) = hires_cent[2] ;
		MatrixMultiply(m_vox_xform, v_ch, v_cl) ;
		dx = V3_X(v_cl) - lowres_cent[0] ;
		dy = V3_Y(v_cl) - lowres_cent[1] ;
		dz = V3_Z(v_cl) - lowres_cent[2] ;
		m_trans = MatrixIdentity(4, NULL) ;
		*MATRIX_RELT(m_trans, 1, 4) = -dx ; *MATRIX_RELT(m_trans, 2, 4) = -dy ; *MATRIX_RELT(m_trans, 3, 4) = -dz ;
		
		m_tmp = MatrixCopy(m_vox_xform, NULL) ;
		MatrixMultiply(m_trans, m_tmp, m_vox_xform) ;
		printf("after aligning centroids:\n") ;
		MatrixPrint(stdout, m_vox_xform) ;
		max_overlap = compute_overlap(vl_lowres, vl_hires, m_vox_xform) ;
		printf("initial overlap = %2.2f...\n", max_overlap) ;
		
		if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
		{
			write_snapshot(mri_lowres, mri_hires, m_vox_xform, parms, parms->start_t, 1,NULL);
		}
		parms->start_t++ ;
		MatrixFree(&m_trans) ; MatrixFree(&m_tmp) ; VectorFree(&v_cl) ; VectorFree(&v_ch) ;
		if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
		{
			write_snapshot(mri_lowres, mri_hires, m_vox_xform, parms, parms->start_t,1,NULL);
		}
		parms->start_t++ ;
	}
	else
		m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;


	trans = MAX(MAX(mri_hires->width,mri_hires->height),mri_hires->depth)/8 ;
	max_overlap = find_optimal_translation(vl_lowres, vl_hires, m_vox_xform, 
																				 -trans, trans, 5, 4) ;
		
#define MIN_SCALES 3
  /////////////////////////// loop here ////////////////////////////////////////////
  niter = 0 ; nscales = 1 ; scale = 1.0 ; good_step = 0 ; done = 0 ;
  do
  {
    old_max_overlap = max_overlap ;
    printf("****************************************\n");
    printf("Nine parameter search.  iteration %d nscales = %d ...\n", 
					 niter+1, nscales);
    printf("****************************************\n");
    max_overlap = find_optimal_linear_xform(vl_lowres, vl_hires, 
																						m_vox_xform, m_origin,
																						-MAX_ANGLE*scale,
																						MAX_ANGLE*scale,
																						1-0.5*scale, 
																						1+0.5*scale, 
																						-scale*MAX_TRANS, 
																						scale*MAX_TRANS,
																						3, 3, 3, 2);
    
    if (parms->write_iterations != 0)
    {
			write_snapshot(mri_lowres, mri_hires, ((LTA *)(transform->xform))->xforms[0].m_L, 
										 parms, parms->start_t+niter, 1, NULL) ;

    }
    printf("Result so far: scale %2.3f: max overlap=%2.2f, old max overlap=%2.2f\n",
					 scale,max_overlap, old_max_overlap) ;
    MatrixPrint(stderr, m_vox_xform);
    /* search a finer nbhd (if do-while continues) */
    if ((max_overlap <= old_max_overlap)) /* couldn't take a step */
    {
      if (good_step)
      {
				scale *= 0.25 ;
				if (scale < min_search_scale)
					break ;
				good_step = 0 ;
				printf("reducing scale to %2.4f\n", scale) ;
				nscales++ ;
      }
      else
				done = 1 ;
    }
    else
      good_step = 1 ; /* took at least one good step at this scale */
    
    niter++ ;
  } while (nscales < MIN_SCALES || (done == FALSE)) ;



	/*	m_vox_xform = compute_pca(mri_hires, mri_lowres) ;*/
  parms->start_t += niter ;
	MatrixFree(&m_origin) ;
	MatrixFree(&m_inv_origin) ;
	return(transform) ;
}

#if 1

/* compute intersection of transformed hires with lowres divided by union */
static double
compute_overlap(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, MATRIX *m_L)
{
	int     intersection, Union, x, y, z, width, height, depth, xd, yd, zd, label,
		hwidth, hheight, hdepth, i ;
	VECTOR  *v1, *v2 ;
	MRI     *mri_lowres, *mri_hires ;
	static MRI *mri_intersection = NULL;

	mri_lowres = vl_lowres->mri ; mri_hires = vl_hires->mri ; 

	v1 = VectorAlloc(4, MATRIX_REAL) ;
	v2 = VectorAlloc(4, MATRIX_REAL) ;
	*MATRIX_RELT(v1, 4, 1) = 1.0 ; *MATRIX_RELT(v2, 4, 1) = 1.0 ;

	width = mri_lowres->width ; 
	height = mri_lowres->height; 
	depth = mri_lowres->depth;
	hwidth = mri_hires->width ; 
	hheight = mri_hires->height ;
	hdepth = mri_hires->depth;

	mri_intersection = mri_from_voxel_list(vl_lowres, mri_intersection) ;

	/* first go through lowres volume and for every voxel that is on in it,
		 map it to the hires, and if the hires is on, add one to the overlap
	*/

	/* go  through lowres volume and for every voxel that is on in it,
		 map it to the hires, and if the hires hasn't been counted yet, count it.
	*/
#define IN_INTERSECTION 255
#define IN_UNION        254

	for (intersection = Union = i = 0 ; i < vl_hires->nvox ; i++)
	{
		x = vl_hires->xi[i] ; y = vl_hires->yi[i] ; z = vl_hires->zi[i] ; 

		V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
		MatrixMultiply(m_L, v1, v2) ;
		xd = nint(V3_X(v2)) ; yd = nint(V3_Y(v2)) ; zd = nint(V3_Z(v2)) ; 
		if (xd >= 0 && xd < width &&
				yd >= 0 && yd < height &&
				zd >= 0 && zd < depth)
		{
			label = MRIvox(mri_intersection, xd, yd, zd) ;
			switch (label)
			{
			case 0:  /* hires mapping outside of lowres label*/
				Union++ ;
				MRIvox(mri_intersection, xd, yd, zd) = IN_UNION ;
				break ;
			case IN_UNION:             /* already processed one way or the other */
			case IN_INTERSECTION:
				break ;
			default:                   /* hires mapping into lowres label */
				intersection++ ;
				Union++ ;
				MRIvox(mri_intersection, xd, yd, zd) = IN_INTERSECTION ;
				break ;
			}
		}
		else
			Union++ ; /* penalize for mapping out of FOV */
	}

	/* now count lowres voxels that weren't mapped to in union */
	for (i = 0 ; i < vl_lowres->nvox ; i++)
	{
		x = vl_lowres->xi[i] ; y = vl_lowres->yi[i] ; z = vl_lowres->zi[i] ; 
		label = MRIvox(mri_intersection, x, y, z) ;
		switch (label)
		{
		case 0:  /* hires mapping outside of lowres label*/
		case IN_UNION:             /* already processed one way or the other */
		case IN_INTERSECTION:
			break ;
		default:                   /* hires mapping into lowres label */
			Union++ ;
			break ;
		}
	}

	/* reset intersection volume */
	for (i = 0 ; i < vl_hires->nvox ; i++)
	{
		x = vl_hires->xi[i] ; y = vl_hires->yi[i] ; z = vl_hires->zi[i] ; 

		V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
		MatrixMultiply(m_L, v1, v2) ;
		xd = nint(V3_X(v2)) ; yd = nint(V3_Y(v2)) ; zd = nint(V3_Z(v2)) ; 
		if (xd >= 0 && xd < width &&
				yd >= 0 && yd < height &&
				zd >= 0 && zd < depth)
		{
			MRIvox(mri_intersection, xd, yd, zd) = 0 ;
		}
	}

	VectorFree(&v1) ; VectorFree(&v2) ;
	return((double)intersection/(double)Union) ;
}
#else
static int
compute_overlap(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, MATRIX *m_L)
{
	int     overlap, x, y, z, width, height, depth, xd, yd, zd, label,
      		hwidth, hheight, hdepth, i ;
	MATRIX  *m_inv ;
	VECTOR  *v1, *v2 ;
	MRI     *mri_lowres_orig, *mri_hires_orig ;
	static MRI *mri_lowres = NULL, *mri_hires = NULL ;

	mri_lowres_orig = vl_lowres->mri ; 
	mri_hires_orig = vl_hires->mri ; 
	m_inv = MatrixInverse(m_L, NULL) ;
	if (m_inv == NULL)
		ErrorExit(ERROR_BADPARM, "compute_overlap: singular matrix") ;

	v1 = VectorAlloc(4, MATRIX_REAL) ;
	v2 = VectorAlloc(4, MATRIX_REAL) ;
	*MATRIX_RELT(v1, 4, 1) = 1.0 ; *MATRIX_RELT(v2, 4, 1) = 1.0 ;

	width = mri_lowres_orig->width ; 
	height = mri_lowres_orig->height; 
	depth = mri_lowres_orig->depth;
	hwidth = mri_hires_orig->width ; 
	hheight = mri_hires_orig->height ;
	hdepth = mri_hires_orig->depth;
	overlap = 0 ;
	mri_hires = mri_from_voxel_list(vl_hires, mri_hires) ;
	mri_lowres = mri_from_voxel_list(vl_lowres, mri_lowres) ;

	/* first go through lowres volume and for every voxel that is on in it,
		 map it to the hires, and if the hires is on, add one to the overlap
	*/
	for (i = 0 ; i < vl_lowres->nvox ; i++)
	{
		x = vl_lowres->xi[i] ; y = vl_lowres->yi[i] ; z = vl_lowres->zi[i] ; 

		V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
		MatrixMultiply(m_inv, v1, v2) ;
		xd = nint(V3_X(v2)) ; yd = nint(V3_Y(v2)) ; zd = nint(V3_Z(v2)) ; 
		if (xd >= 0 && xd < hwidth &&
				yd >= 0 && yd < hheight &&
				zd >= 0 && zd < hdepth)
		{
			label = MRIvox(mri_hires, xd, yd, zd) ;
			if (label > 0 && label != 255)
			{
				overlap++ ;
				MRIvox(mri_hires, xd, yd, zd) = 255 ; /* only count it once */
			}
			else if (label == 0)
				overlap-- ;
		}
		else
			overlap-- ;
	}

	/* now go through lowres volume and for every voxel that is on in it,
		 map it to the hires, and if the hires hasn't been counted yet, count it.
	*/
	for (i = 0 ; i < vl_hires->nvox ; i++)
	{
		x = vl_hires->xi[i] ; y = vl_hires->yi[i] ; z = vl_hires->zi[i] ; 

		V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
		MatrixMultiply(m_L, v1, v2) ;
		xd = nint(V3_X(v2)) ; yd = nint(V3_Y(v2)) ; zd = nint(V3_Z(v2)) ; 
		if (xd >= 0 && xd < width &&
				yd >= 0 && yd < height &&
				zd >= 0 && zd < depth)
		{
			label = MRIvox(mri_lowres, xd, yd, zd) ;
			if (label > 0 && label != 255)
			{
				MRIvox(mri_lowres, xd, yd, zd) = 255 ;
				overlap++ ;
			}
			else if (label == 0)
				overlap-- ;
		}
		else
			overlap-- ;
	}

	if (overlap == 0)
		DiagBreak() ;
	MatrixFree(&m_inv) ; VectorFree(&v1) ; VectorFree(&v2) ;
	return(overlap) ;
}
#endif
static double
find_optimal_translation(VOXEL_LIST *vl_lowres,  VOXEL_LIST *vl_hires,
                         MATRIX *m_L, float min_trans, float max_trans, 
                         float trans_steps, int nreductions)
{
	MRI      *mri_lowres, *mri_hires ;
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta, 
           overlap, max_overlap, mean_trans ;
  int      i ;

	mri_lowres = vl_lowres->mri ;
	mri_hires = vl_hires->mri ;
  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_overlap = compute_overlap(vl_lowres, vl_hires, m_L) ;

  for (i = 0 ; i <= nreductions ; i++)
  {
    delta = (max_trans-min_trans) / trans_steps ;
    if (FZERO(delta))
      return(max_overlap) ;
    if (Gdiag & DIAG_SHOW)
    {
      printf(
						 "scanning translations %2.2f->%2.2f (step %2.1f) ",
						 min_trans,max_trans, delta) ;
      fflush(stdout) ;
    }
    for (x_trans = min_trans ; x_trans <= max_trans ; x_trans += delta)
    {
      *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
      for (y_trans = min_trans ; y_trans <= max_trans ; y_trans += delta)
      {
        *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
        for (z_trans= min_trans ; z_trans <= max_trans ; z_trans += delta)
        {
          *MATRIX_RELT(m_trans, 3, 4) = z_trans ;
          if (nint((x_trans)) == -9 && nint((y_trans)) == -5 &&
              nint((z_trans)) == -7)
            DiagBreak() ;

					// get the transform
          m_L_tmp = MatrixMultiply(m_trans, m_L, m_L_tmp) ;
					// calculate the overlap
          overlap = compute_overlap(vl_lowres, vl_hires, m_L_tmp) ;
          if (overlap > max_overlap)
          {
            max_overlap = overlap ;
            x_max = x_trans ; y_max = y_trans ; z_max = z_trans ;
#if 1
            printf("new max overlap %2.2f found at (%2.1f, %2.1f, %2.1f)\n",
									 max_overlap, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
      printf(
						 "max overlap = %2.2f @ (%2.1f, %2.1f, %2.1f)\n", 
						 max_overlap, x_max, y_max, z_max) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max ;
    // create a new transform by multiplying the previous one.
    MatrixMultiply(m_trans, m_L, m_L_tmp) ;
    MatrixCopy(m_L_tmp, m_L) ;
    x_max = y_max = z_max = 0.0 ;  /* we've translated transform by old maxs */

    mean_trans = (max_trans + min_trans) / 2 ;
    delta = (max_trans-min_trans)/4 ;
    min_trans = mean_trans - delta ;
    max_trans = mean_trans + delta ;
  }

  printf("\n") ;

  MatrixFree(&m_trans) ;
  return(max_overlap) ;
}
static double
find_optimal_linear_xform(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires,
													MATRIX *m_L, MATRIX *m_origin,
													float min_angle, float max_angle,
													float min_scale, float max_scale,
													float min_trans, float max_trans,
													float angle_steps, float scale_steps, 
													float trans_steps,
													int nreductions)
{
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
           *m_tmp2, *m_scale, *m_trans, *m_tmp3 = NULL ;
  double   x_angle, y_angle, z_angle, x_max_rot, y_max_rot, z_max_rot, 
           delta_rot, x_max_scale, y_max_scale, z_max_scale, delta_scale, 
           x_trans, delta_trans, y_trans, z_trans,
           mean_angle, x_scale, y_scale, z_scale, mean_scale, x_max_trans,
           y_max_trans, overlap, max_overlap, z_max_trans, mean_trans ;
  int      i ;
	MRI      *mri_lowres, *mri_hires ;

	mri_lowres = vl_lowres->mri ;
	mri_hires = vl_hires->mri ;
  m_trans = MatrixIdentity(4, NULL) ;
  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max_trans = y_max_trans = z_max_trans = x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_overlap = compute_overlap(vl_lowres, vl_hires, m_L) ;
  for (i = 0 ; i < nreductions ; i++)
  {
    delta_trans = (max_trans-min_trans) / (trans_steps-1) ;
    delta_scale = (max_scale-min_scale) / (scale_steps-1) ;
    delta_rot = (max_angle-min_angle) / (angle_steps-1) ;
    if (Gdiag & DIAG_SHOW)
    {
      printf("  scanning %2.2f degree nbhd (%2.1f)\n"
             "  scale %2.3f->%2.3f (step %2.3f), trans %2.2f->%2.2f (step %2.2f)\n",
						 (float)DEGREES(max_angle), (float)DEGREES(delta_rot),
             min_scale,max_scale, delta_scale, min_trans, max_trans, delta_trans);
      fflush(stdout) ;
    }

    // scale /////////////////////////////////////////////////////////////
    for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta_scale)
    {
      /*      printf("x_scale = %2.3f\n", x_scale) ;*/
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ; y_scale <= max_scale ; y_scale += delta_scale)
      {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ; z_scale <= max_scale; z_scale += delta_scale)
        {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) = 
            *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

					// angle /////////////////////////////////////////////////////////////
          for (x_angle = min_angle ; x_angle <= max_angle ; x_angle += delta_rot)
          {
            m_x_rot = MatrixReallocRotation(4, x_angle, X_ROTATION, m_x_rot) ;
            for (y_angle = min_angle ; y_angle <= max_angle ; y_angle += delta_rot)
            {
              m_y_rot = MatrixReallocRotation(4, y_angle, Y_ROTATION, m_y_rot);
              m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
              for (z_angle= min_angle; z_angle <= max_angle; z_angle += delta_rot)
              {
                m_z_rot = MatrixReallocRotation(4, z_angle,Z_ROTATION,m_z_rot);
                m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
                m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
                MatrixMultiply(m_origin, m_tmp2, m_rot) ;

                
                m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
								m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3) ;

								// translation //////////////////////////////////////////////////////
								for (x_trans = min_trans ; x_trans <= max_trans ; x_trans += delta_trans)
								{
									*MATRIX_RELT(m_trans, 1, 4) = x_trans ;
									for (y_trans = min_trans ; y_trans <= max_trans ; y_trans += delta_trans)
									{
										*MATRIX_RELT(m_trans, 2, 4) = y_trans ;
										for (z_trans= min_trans ; z_trans <= max_trans ; z_trans += delta_trans)
										{
											*MATRIX_RELT(m_trans, 3, 4) = z_trans ;

											m_L_tmp = MatrixMultiply(m_trans, m_tmp3, m_L_tmp) ;
											overlap = compute_overlap(vl_lowres, vl_hires, m_L_tmp) ;
		      
											if (overlap > max_overlap)
											{
												max_overlap = overlap ;
												x_max_scale = x_scale ; y_max_scale = y_scale ; 
												z_max_scale = z_scale ;
												x_max_rot = x_angle ; y_max_rot = y_angle ; 
												z_max_rot = z_angle ;
												x_max_trans = x_trans ; y_max_trans = y_trans ; z_max_trans = z_trans ;
											}
										}
									}
								}
              }
            }
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
    {
      printf("  max overlap = %2.2f @ R=(%2.3f,%2.3f,%2.3f),"
						 "S=(%2.3f,%2.3f,%2.3f), T=(%2.1f,%2.1f,%2.1f)\n", 
						 max_overlap, DEGREES(x_max_rot), DEGREES(y_max_rot),
						 DEGREES(z_max_rot),x_max_scale, y_max_scale, z_max_scale,
						 x_max_trans, y_max_trans,z_max_trans) ;
    }

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_scale, 1, 4) = 
      *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
    *MATRIX_RELT(m_scale,1,1) = x_max_scale ;
    *MATRIX_RELT(m_scale,2,2) = y_max_scale ;
    *MATRIX_RELT(m_scale,3,3) = z_max_scale ;
    m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
    MatrixMultiply(m_origin, m_tmp, m_scale) ;


    x_max_scale = y_max_scale = z_max_scale = 1.0 ;

    mean_scale = (max_scale + min_scale) / 2 ;
    delta_scale = (max_scale-min_scale)/4 ;
    min_scale = mean_scale - delta_scale ;
    max_scale = mean_scale + delta_scale ;

    /* update L to reflect new maximum and search around it */
    MatrixReallocRotation(4, x_max_rot, X_ROTATION, m_x_rot) ;
    MatrixReallocRotation(4, y_max_rot, Y_ROTATION, m_y_rot) ;
    MatrixReallocRotation(4, z_max_rot, Z_ROTATION, m_z_rot) ;
    MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
    MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
    m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
    MatrixMultiply(m_origin, m_tmp2, m_rot) ;

    m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
    m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max_trans ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max_trans ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max_trans ;
    MatrixMultiply(m_trans, m_tmp3, m_L_tmp) ;

    MatrixCopy(m_L_tmp, m_L) ;


    x_max_trans = y_max_trans = z_max_trans = 0.0 ;  /* we've translated transform by old maxs */
    mean_trans = (max_trans + min_trans) / 2 ;
    delta_trans = (max_trans-min_trans)/4 ;
    min_trans = mean_trans - delta_trans ;
    max_trans = mean_trans + delta_trans ;

    /* we've rotated transform to old max */
    x_max_rot = y_max_rot = z_max_rot = 0.0 ;

    mean_angle = (max_angle + min_angle) / 2 ;
    delta_rot = (max_angle-min_angle)/4 ;
    min_angle = mean_angle - delta_rot ;
    max_angle = mean_angle + delta_rot ;
  }
  MatrixFree(&m_x_rot) ; MatrixFree(&m_y_rot) ; MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;   MatrixFree(&m_tmp) ; MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ; MatrixFree(&m_trans) ; MatrixFree(&m_tmp3) ;
  return(max_overlap) ;
}
#define NPARMS (4*4)
#ifdef TOL
#undef TOL
#endif
#define TOL 1e-12

#include "nrutil.h"
static VOXEL_LIST *Gvl_lowres, *Gvl_hires ;

static int
powell_minimize(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, MATRIX *mat)
{
	float *p, **xi, fret, fstart;
	int   i, r, c, iter ;

	p = vector(1, NPARMS) ;
	xi = matrix(1, NPARMS, 1, NPARMS) ;
	for (i = r = 1 ; r <= 4 ; r++)
	{
		for (c = 1 ; c <= 4 ; c++)
		{
			p[i++] = *MATRIX_RELT(mat, r, c) ;
		}
	}

	Gvl_lowres = vl_lowres ; Gvl_hires = vl_hires ;
	for (r = 1 ; r <= NPARMS ; r++)
	{
		for (c = 1 ; c <= NPARMS ; c++)
		{
			xi[r][c] = r == c ? 1 : 0 ;
		}
	}

	powell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
	do
	{
		for (r = 1 ; r <= NPARMS ; r++)
		{
			for (c = 1 ; c <= NPARMS ; c++)
			{
				xi[r][c] = r == c ? 1 : 0 ;
			}
		}

		fstart = fret ;
		powell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
		for (i = r = 1 ; r <= 4 ; r++)
		{
			for (c = 1 ; c <= 4 ; c++)
			{
				*MATRIX_RELT(mat, r, c) = p[i++] ;
			}
		}
		*MATRIX_RELT(mat, 4, 1) = 0.0 ; *MATRIX_RELT(mat, 4, 2) = 0.0 ; 
		*MATRIX_RELT(mat, 4, 3) = 0.0 ; *MATRIX_RELT(mat, 4, 4) = 1.0 ; 
		printf("%3.3d: best alignment at after powell: %2.3f (%d steps)\n", parms.start_t,fret, iter) ;
		write_snapshot(vl_lowres->mri, vl_hires->mri, mat, &parms, parms.start_t++,1,NULL);
	} while (fret < fstart) ;

	free_matrix(xi, 1, NPARMS, 1, NPARMS) ;
	free_vector(p, 1, NPARMS) ;
	return(NO_ERROR) ;
}

static float
compute_powell_sse(float *p)
{
	static MATRIX *mat = NULL ;
	float  error ;
	int    i, r, c ;

	if (mat == NULL)
		mat = MatrixAlloc(4, 4, MATRIX_REAL) ;
	for (i = r = 1 ; r <= 4 ; r++)
	{
		for (c = 1 ; c <= 4 ; c++)
		{
			*MATRIX_RELT(mat, r, c) = p[i++] ;
		}
	}
	*MATRIX_RELT(mat, 4, 1) = 0.0 ; *MATRIX_RELT(mat, 4, 2) = 0.0 ; 
	*MATRIX_RELT(mat, 4, 3) = 0.0 ; *MATRIX_RELT(mat, 4, 4) = 1.0 ; 
	error = -compute_overlap(Gvl_lowres, Gvl_hires, mat) ;
	return(error) ;
}

static int
write_snapshot(MRI *mri_lowres, MRI *mri_hires, MATRIX *m_vox_xform, 
							 INTEGRATION_PARMS *parms, int fno, int conform, char *in_fname)
{
	MRI *mri_aligned ;
	char fname[STRLEN] ;

	if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
	{
		printf("hires->lowres vox->vox transform:\n") ;
		MatrixPrint(stdout, m_vox_xform) ;
	}
	if (conform)
	{
		mri_aligned = MRIclone(mri_lowres, NULL) ;
		MRIlinearTransformInterp(mri_hires, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
	}
	else
	{
		MATRIX *m_ras ;

    m_ras = MRIvoxelXformToRasXform(mri_hires, mri_lowres, m_vox_xform, NULL) ;
		
		mri_aligned = MRIapplyRASlinearTransformInterp(mri_hires, NULL, m_ras, SAMPLE_NEAREST);
		MatrixFree(&m_ras) ;
	}
	sprintf(fname, "%s%03d", parms->base_name, fno) ;
	MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
	if (in_fname)
		strcpy(fname, in_fname) ;
	else
		sprintf(fname, "%s%03d.mgz", parms->base_name, fno) ;
	printf("writing snapshot to %s...\n", fname) ;
	MRIwrite(mri_aligned, fname) ;
	MRIfree(&mri_aligned) ;

	{
		mri_aligned = MRIsrcTransformedCentered(mri_hires, mri_lowres, m_vox_xform,
																						SAMPLE_NEAREST) ;
		sprintf(fname, "orig_%s%03d.mgz", parms->base_name, fno) ;
		printf("writing snapshot to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}

	if (mri_hires_intensity)
	{
		mri_aligned = MRIsrcTransformedCentered(mri_hires_intensity, mri_lowres, m_vox_xform, SAMPLE_TRILINEAR) ;
		sprintf(fname, "intensity_%s%03d", parms->base_name, fno) ;
		MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
		sprintf(fname, "intensity_%s%03d.mgz", parms->base_name, fno) ;
		printf("writing snapshot to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}

	return(NO_ERROR) ;
}

static VOXEL_LIST *
create_voxel_list(MRI *mri, float low_val, float hi_val, VOXEL_LIST *vl, 
									int skip)
{
	int   x, y, z, nvox, i ;
	Real  val ;

	skip++ ;  /* next voxel + amount to skip */
	for (nvox = x = 0 ; x < mri->width ; x+=skip)
	{
		for (y = 0 ; y < mri->height ; y+=skip)
		{
			for (z = 0 ; z < mri->depth ; z+=skip)
			{
				val = MRIgetVoxVal(mri, x, y, z, 0) ;
				if (val > 0)
					DiagBreak() ;
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				if (val >= low_val && val <= hi_val)
					nvox++ ;
			}
		}
	}

	printf("allocating %d voxel indices...\n", nvox) ;
	if (vl == NULL)
	{
		vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
		vl->nvox = nvox ;
		vl->xi = (int *)calloc(nvox, sizeof(int)) ;
		vl->yi = (int *)calloc(nvox, sizeof(int)) ;
		vl->zi = (int *)calloc(nvox, sizeof(int)) ;
		if (!vl || !vl->xi || !vl->yi || !vl->zi)
			ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
								Progname, nvox) ;
	}
	for (nvox = x = 0 ; x < mri->width ; x+=skip)
	{
		for (y = 0 ; y < mri->height ; y+=skip)
		{
			for (z = 0 ; z < mri->depth ; z+=skip)
			{
				val = MRIgetVoxVal(mri, x, y, z, 0) ;
				if (val >= low_val && val <= hi_val)
				{
					i = nvox++ ;
					vl->xi[i] = x ; vl->yi[i] = y ; vl->zi[i] = z ;
				}
			}
		}
	}
	vl->mri = mri ;
	return(vl) ;
}

static int
free_voxel_list(VOXEL_LIST **pvl)
{
	VOXEL_LIST *vl = *pvl ;
	*pvl = NULL ;

	if (!vl)
		return(ERROR_BADPARM) ;
	free(vl->xi) ;
	free(vl->yi) ;
	free(vl->zi) ;
	free(vl) ;
	return(NO_ERROR) ;
}
static MRI *
mri_from_voxel_list(VOXEL_LIST *vl, MRI *mri)
{
	int   i ;
	Real  val ;

	if (mri == NULL)
		mri = MRIclone(vl->mri, NULL) ;

	for (i = 0 ; i < vl->nvox ; i++)
	{
		val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0) ;
		MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val) ;
	}
	return(mri) ;
}

