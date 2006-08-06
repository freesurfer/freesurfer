//
// mri_hires_register.c
//
// written by Bruce Fischl
// Nov. 9th ,2000
// 
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2006/08/06 20:55:50 $
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
#include "nr_wrapper.h"
#include "fastmarching.h"
#include "voxlist.h"

#define DEFAULT_MAX_ANGLE       RADIANS(25)
static double MAX_ANGLE = DEFAULT_MAX_ANGLE ;
static double MAX_SCALE = 0.5 ;

#define PAD       10


static int find_label = -1 ;
static float x_vox = 0.0 ;
static float y_vox = 0.0 ;
static float z_vox = 0.0 ;

static int apply_transform = 1 ;

static int write_snapshot(MRI *mri_target, MRI *mri_source, 
													MATRIX *m_vox_xform, MORPH_PARMS *parms, 
													int fno, int conform, char *fname) ;

static double MAX_TRANS = 30 ;

static float compute_powell_sse(float *p) ;
static int    powell_minimize(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *mat);

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;

static TRANSFORM *compute_optimal_transform(VOXEL_LIST *vl_target, 
																						VOXEL_LIST *vl_source, 
																						MORPH_PARMS *parms,
																						TRANSFORM *transform) ;

static double compute_likelihood(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *m_L, float intensity_scale) ;

static double (*pf_likelihood)(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *m_L, float intensity_scale) = compute_likelihood ;

static double find_optimal_translation(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, 
																			 MATRIX *m_L, float min_trans, 
																			 float max_trans, 
																			 float trans_steps, int nreductions);

static double find_optimal_linear_xform(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source,
																				MATRIX *m_L, MATRIX *m_origin,
																				float min_angle, float max_angle,
																				float min_scale, float max_scale,
																				float min_trans, float max_trans,
																				float angle_steps, float scale_steps, 
																				float trans_steps,
																				int nreductions,
																				int rigid);
char *Progname ;

static int skip = 2 ;
static double distance = 1.0 ;



static MORPH_PARMS parms ;
static TRANSFORM  *transform = NULL ;

int
main(int argc, char *argv[])
{
	char       **av, *source_fname, *target_fname, *out_fname, fname[STRLEN] ;
  int        ac, nargs, i ;
	MRI        *mri_target, *mri_source, *mri_tmp ;
	VOXEL_LIST *vl_target, *vl_source ;
	MRI_REGION  box ;
  struct timeb start ;
  int          msec, minutes, seconds ;
	float        fmin, fmax ;

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

	source_fname = argv[1] ;
	target_fname = argv[2] ;
	out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
	mri_source = MRIread(source_fname) ;
	if (!mri_source)
		ErrorExit(ERROR_NOFILE, "%s: could not read source label volume %s",
							Progname, source_fname) ;

	
	pf_likelihood = compute_likelihood ;

	MRIboundingBox(mri_source, 0, &box) ;
	box.x -= PAD ; box.y -= PAD ; box.z -= PAD ; 
	box.dx += 2*PAD ; box.dy += 2*PAD ; box.dz += 2*PAD ; 
	MRIcropBoundingBox(mri_source, &box) ;
	mri_tmp = MRIextractRegion(mri_source, NULL, &box) ;
	MRIfree(&mri_source) ;
	mri_source = mri_tmp ;

	mri_target = MRIread(target_fname) ;
	if (!mri_target)
		ErrorExit(ERROR_NOFILE, "%s: could not read target label volume %s",
							Progname, target_fname) ;

  if (Gdiag & DIAG_WRITE && parms.write_iterations > 0)
  {
		sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_target, fname, IMAGE_SIZE) ;	
		sprintf(fname, "%s_target.mgz", parms.base_name) ;
		MRIwrite(mri_target, fname) ;
	}

	MRImatchMeanIntensity(mri_source, mri_target, mri_source) ;
	
	/* compute optimal linear transform */
	MRIvalRange(mri_target, &fmin, &fmax) ;
	vl_target = VLSTcreate(mri_target,1,fmax+1,NULL,skip,0);
	vl_target->mri2 = mri_target ;
	for (i = 0 ; i < 1 ; i++)
	{
		printf("------------- outer loop iteration %d ---------------\n",i) ;
		MRIvalRange(mri_target, &fmin, &fmax) ;
		vl_source = VLSTcreate(mri_source, 1, fmax+1, NULL, skip, 0) ;
		vl_source->mri2 = mri_source ;
		
		transform = compute_optimal_transform(vl_target, vl_source, &parms,
																					transform) ;
		VLSTfree(&vl_source) ;
		vl_source = VLSTcreate(mri_source, 1, fmax+1, NULL, skip/4, 0) ;
		vl_source->mri2 = mri_source ;
		if (parms.rigid == 0)
		{
			powell_minimize(vl_target, vl_source, ((LTA *)(transform->xform))->xforms[0].m_L) ;
			if (parms.rigid)
				MatrixOrthonormalizeTransform(((LTA *)(transform->xform))->xforms[0].m_L) ;
		}
		VLSTfree(&vl_source) ;
		if (apply_transform)
		{
			MRI *mri_aligned ;
			MATRIX *m_vox_xform ;
			char   fname[STRLEN] ;
			
			FileNameRemoveExtension(out_fname, fname) ;
			strcat(fname, ".mgz") ;
			m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;
			mri_aligned = MRIclone(mri_target, NULL) ;
			MRIlinearTransformInterp(mri_source, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
			printf("writing transformed output volume to %s...\n", fname) ;
			MRIwrite(mri_aligned, fname) ;
			MRIfree(&mri_aligned) ;
		}
		LTAvoxelToRasXform((LTA *)(transform->xform), mri_source, mri_target) ;
		TransformWrite(transform, out_fname) ;
		LTArasToVoxelXform((LTA *)(transform->xform), mri_source, mri_target) ;
		skip /= 2 ;
	}
	VLSTfree(&vl_target) ;
	
	printf("final vox2vox matrix:\n") ;
	MatrixPrint(stdout, ((LTA *)(transform->xform))->xforms[0].m_L) ;
	if (apply_transform)
	{
		MRI *mri_aligned ;
		MATRIX *m_vox_xform ;
		char   fname[STRLEN] ;
		
		FileNameRemoveExtension(out_fname, fname) ;
		strcat(fname, ".mgz") ;
		m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;
		mri_aligned = MRIclone(mri_target, NULL) ;
		MRIlinearTransformInterp(mri_source, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
		printf("writing transformed output volume to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}
	LTAvoxelToRasXform((LTA *)(transform->xform), mri_source, mri_target) ;
	TransformWrite(transform, out_fname) ;

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
  else if (!stricmp(option, "find_label"))
  {
		find_label = atoi(argv[2]) ;
		x_vox = atof(argv[3]) ;
		y_vox = atof(argv[4]) ;
		z_vox = atof(argv[5]) ;
		nargs = 4 ;
    printf("finding label %s (%d) at (%2.1f, %2.1f, %2.1f)\n",
					 cma_label_to_name(find_label), find_label, x_vox, y_vox,z_vox) ;
  }
  else if (!stricmp(option, "distance"))
  {
    distance = atof(argv[2]) ;
		nargs = 1 ;
    printf("expanding border by %2.1f mm every outer cycle\n", distance);
  }
  else if (!stricmp(option, "view"))
	{
		Gsx = atoi(argv[2]) ;
		Gsy = atoi(argv[3]) ;
		Gsz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("viewing voxel (%d, %d, %d)\n", Gsx, Gsy, Gsz) ;
	}
  else if (!stricmp(option, "trans"))
  {
		MAX_TRANS = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting MAX_TRANS = %2.1f\n", MAX_TRANS) ;
  }
  else if (!stricmp(option, "max_angle"))
  {
		MAX_ANGLE = RADIANS(atof(argv[2])) ;
    nargs = 1 ;
    printf("setting max angle for search to %2.1f degrees\n", DEGREES(MAX_ANGLE)) ;
  }
  else if (!stricmp(option, "max_scale"))
  {
		MAX_SCALE = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting max scale for search to %2.1f --> %2.1f\n", 
					 1-MAX_SCALE, 1+MAX_SCALE) ;
  }
	else if (!stricmp(option, "skip"))
	{
		skip = atoi(argv[2]);
		printf("skipping %d voxels in source data...\n", skip) ;
		nargs = 1 ;
	}
	else switch (*option)
	{
	case 'W':
		parms.write_iterations = atoi(argv[2]) ;
		Gdiag |= DIAG_WRITE ;
		nargs = 1 ;
		printf("setting write iterations = %d\n", parms.write_iterations) ;
		break ;
	case 'R':
		parms.rigid = 1 ;
		printf("constraining transform to be rigid...\n") ;
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
	printf("usage: %s <source> <target> <output xform>\n",Progname) ;
	exit(ecode) ;
}


static TRANSFORM *
compute_optimal_transform(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, 
													MORPH_PARMS *parms, TRANSFORM *transform)
{
	MATRIX    *m_vox_xform, *m_origin, *m_inv_origin, *m_source_vox2ras, 
		        *m_target_ras2vox, *m_trans, *m_tmp ;
	VECTOR    *v_cl, *v_ch ;
	double    source_cent[3], target_cent[3], dx, dy, dz, scale,min_search_scale ;
	int       niter, nscales, good_step, done, trans ;
	double    old_max_likelihood, max_likelihood ;
	MRI       *mri_target, *mri_source ;

	mri_target = vl_target->mri ; mri_source = vl_source->mri ;

#define MIN_SEARCH_SCALE 0.01
  min_search_scale = MIN_SEARCH_SCALE ;
	m_origin = MatrixIdentity(4, NULL) ;
	MRIcenterOfMass(mri_source, source_cent, 0) ;
	MRIcenterOfMass(mri_target, target_cent, 0) ;
	*MATRIX_RELT(m_origin, 1, 4) = target_cent[0] ; 
	*MATRIX_RELT(m_origin, 2, 4) = target_cent[1] ;
	*MATRIX_RELT(m_origin, 3, 4) = target_cent[2] ; 
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;
	m_inv_origin = MatrixInverse(m_origin, NULL) ;
	if (transform == NULL)
	{
		transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
		m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;
		
		m_target_ras2vox = MRIgetRasToVoxelXform(mri_target) ;
		m_source_vox2ras = MRIgetVoxelToRasXform(mri_source) ;
		MatrixMultiply(m_target_ras2vox, m_source_vox2ras, m_vox_xform) ;
		printf("initial transform from direction cosines:\n") ;
		MatrixPrint(stdout, m_vox_xform) ;
		MatrixFree(&m_target_ras2vox) ; MatrixFree(&m_source_vox2ras) ;
		
		dx = target_cent[0] - source_cent[0] ; dy = target_cent[1] - source_cent[1]  ;
		dz = target_cent[2] - source_cent[2] ;
		
		v_cl = VectorAlloc(4, MATRIX_REAL) ; v_ch = VectorAlloc(4, MATRIX_REAL) ; 
		*MATRIX_RELT(v_cl,4,1) = 1.0 ; *MATRIX_RELT(v_ch,4,1) = 1.0 ;
		V3_X(v_ch) = source_cent[0] ; V3_Y(v_ch) = source_cent[1] ; V3_Z(v_ch) = source_cent[2] ;
		MatrixMultiply(m_vox_xform, v_ch, v_cl) ;
		dx = V3_X(v_cl) - target_cent[0] ;
		dy = V3_Y(v_cl) - target_cent[1] ;
		dz = V3_Z(v_cl) - target_cent[2] ;
		m_trans = MatrixIdentity(4, NULL) ;
		*MATRIX_RELT(m_trans, 1, 4) = -dx ; *MATRIX_RELT(m_trans, 2, 4) = -dy ; *MATRIX_RELT(m_trans, 3, 4) = -dz ;
		
		m_tmp = MatrixCopy(m_vox_xform, NULL) ;
		MatrixMultiply(m_trans, m_tmp, m_vox_xform) ;
		printf("after aligning centroids:\n") ;
		MatrixPrint(stdout, m_vox_xform) ;
		max_likelihood = (*pf_likelihood)(vl_target, vl_source, m_vox_xform, 1.0) ;
		printf("initial likelihood = %2.4f...\n", max_likelihood) ;
		
		MatrixFree(&m_trans) ; MatrixFree(&m_tmp) ; VectorFree(&v_cl) ; VectorFree(&v_ch) ;
		if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
		{
			write_snapshot(mri_target, mri_source, m_vox_xform, parms, parms->start_t,1,NULL);
		}
		parms->start_t++ ;
	}
	else
		m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;


	trans = MAX(MAX_TRANS, MAX(MAX(mri_source->width,mri_source->height),mri_source->depth)/8) ;
	max_likelihood = find_optimal_translation(vl_target, vl_source, m_vox_xform, 
																				 -trans, trans, 5, 4) ;
		
	if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
	{
		write_snapshot(mri_target, mri_source, m_vox_xform, parms, parms->start_t,1,NULL);
	}
	parms->start_t++ ;
#define MIN_SCALES 3
  /////////////////////////// loop here ////////////////////////////////////////////
  niter = 0 ; nscales = 1 ; scale = parms->rigid ? .25: 1.0 ; 
	good_step = 0 ; done = 0 ;
  do
  {
    old_max_likelihood = max_likelihood ;
    printf("****************************************\n");
    printf("Nine parameter search.  iteration %d nscales = %d ...\n", 
					 niter+1, nscales);
    printf("****************************************\n");
    max_likelihood = find_optimal_linear_xform(vl_target, vl_source, 
																							 m_vox_xform, m_origin,
																							 -MAX_ANGLE*scale,
																							 MAX_ANGLE*scale,
																							 1-MAX_SCALE*scale, 
																							 1+MAX_SCALE*scale, 
																							 -scale*MAX_TRANS, 
																							 scale*MAX_TRANS,
																							 3, 3, 3, 2, parms->rigid);
    
    if (parms->write_iterations != 0)
    {
			write_snapshot(mri_target, mri_source, ((LTA *)(transform->xform))->xforms[0].m_L, 
										 parms, parms->start_t+niter, 1, NULL) ;

    }
    printf("Result so far: scale %2.3f: max likelihood = %2.4f, old max likelihood=%2.4f\n",
					 scale,max_likelihood, old_max_likelihood) ;
    MatrixPrint(stderr, m_vox_xform);
    /* search a finer nbhd (if do-while continues) */
    if ((max_likelihood <= old_max_likelihood)) /* couldn't take a step */
    {
			scale *= 0.25 ;
			if (scale < min_search_scale)
				break ;
			good_step = 0 ;
			printf("reducing scale to %2.4f\n", scale) ;
			nscales++ ;
			//			done = (good_step == 0) ;
    }
    else
      good_step = 1 ; /* took at least one good step at this scale */
    
    niter++ ;
  } while (nscales < MIN_SCALES || (done == FALSE)) ;



	/*	m_vox_xform = compute_pca(mri_source, mri_target) ;*/
  parms->start_t += niter ;
	MatrixFree(&m_origin) ;
	MatrixFree(&m_inv_origin) ;
	return(transform) ;
}


/* compute intersection of transformed source with target divided by union */
static double
compute_likelihood(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *m_L, float intensity_scale)
{
	int     x, y, z, width, height, depth,
		      hwidth, hheight, hdepth, i ;
	VECTOR  *v1, *v2 ;
	MRI     *mri_target, *mri_source ;
	double  sse, error ;
	Real    d1, d2, xd, yd, zd ;
	MATRIX  *m_L_inv ;

	m_L_inv = MatrixInverse(m_L, NULL) ;
	if (m_L_inv == NULL)
		ErrorExit(ERROR_BADPARM, "compute_distance_transform_sse: singular matrix.") ;

	mri_target = vl_target->mri2 ; mri_source = vl_source->mri2 ; 

	v1 = VectorAlloc(4, MATRIX_REAL) ;
	v2 = VectorAlloc(4, MATRIX_REAL) ;
	*MATRIX_RELT(v1, 4, 1) = 1.0 ; *MATRIX_RELT(v2, 4, 1) = 1.0 ;

	width = mri_target->width ; height = mri_target->height; depth = mri_target->depth;
	hwidth = mri_source->width ; hheight = mri_source->height ;hdepth = mri_source->depth;


	/* go through both voxel lists and compute the sse
		 map it to the source, and if the source hasn't been counted yet, count it.
	*/

	sse = 0.0 ;
	for (i = 0 ; i < vl_source->nvox ; i++)
	{
		x = vl_source->xi[i] ; y = vl_source->yi[i] ; z = vl_source->zi[i] ; 

		V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
		MatrixMultiply(m_L, v1, v2) ;
		d1 = MRIgetVoxVal(vl_source->mri2, x, y, z, 0) ;
		xd = V3_X(v2) ; yd = V3_Y(v2) ; zd = V3_Z(v2) ; 
		if (xd < 0)
			xd = 0 ;
		else if (xd >= width-1)
			xd = width-1 ;
		if (yd < 0)
			yd = 0 ;
		else if (yd >= height-1)
			yd = height-1 ;
		if (zd < 0)
			zd = 0 ;
		else if (zd >= depth-1)
			zd = depth-1 ;
		MRIsampleVolume(vl_target->mri2, xd, yd, zd, &d2) ;
		error = d1*intensity_scale-d2 ;
		sse += error*error ;
	}

	/* now count target voxels that weren't mapped to in union */
	for (i = 0 ; i < vl_target->nvox ; i++)
	{
		x = vl_target->xi[i] ; y = vl_target->yi[i] ; z = vl_target->zi[i] ; 
		V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
		MatrixMultiply(m_L_inv, v1, v2) ;
		d1 = MRIgetVoxVal(vl_target->mri2, x, y, z, 0) ;

		xd = V3_X(v2) ; yd = V3_Y(v2) ; zd = V3_Z(v2) ; 
		if (xd < 0)
			xd = 0 ;
		else if (xd >= hwidth-1)
			xd = hwidth-1 ;
		if (yd < 0)
			yd = 0 ;
		else if (yd >= hheight-1)
			yd = hheight-1 ;
		if (zd < 0)
			zd = 0 ;
		else if (zd >= hdepth-1)
			zd = hdepth-1 ;
		MRIsampleVolume(vl_source->mri2, xd, yd, zd, &d2) ;
		error = d1*intensity_scale-d2 ;
		sse += error*error ;
	}

	VectorFree(&v1) ; VectorFree(&v2) ;
	MatrixFree(&m_L_inv) ;
	return(-sqrt(sse / (double)(vl_target->nvox + vl_source->nvox))) ;
}

static double
find_optimal_translation(VOXEL_LIST *vl_target,  VOXEL_LIST *vl_source,
                         MATRIX *m_L, float min_trans, float max_trans, 
                         float trans_steps, int nreductions)
{
	MRI      *mri_target, *mri_source ;
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta, 
           likelihood, max_likelihood, mean_trans ;
  int      i ;

	mri_target = vl_target->mri ;
	mri_source = vl_source->mri ;
  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_likelihood = (*pf_likelihood)(vl_target, vl_source, m_L, 1.0) ;

  for (i = 0 ; i <= nreductions ; i++)
  {
    delta = (max_trans-min_trans) / trans_steps ;
    if (FZERO(delta))
      return(max_likelihood) ;
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
					// calculate the likelihood
          likelihood = (*pf_likelihood)(vl_target, vl_source, m_L_tmp, 1.0) ;
          if (likelihood > max_likelihood)
          {
            max_likelihood = likelihood ;
            x_max = x_trans ; y_max = y_trans ; z_max = z_trans ;
#if 1
            printf("new max likelihood %2.4f found at (%2.1f, %2.1f, %2.1f)\n",
									 max_likelihood, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
      printf(
						 "max likelihood = %2.4f @ (%2.1f, %2.1f, %2.1f)\n", 
						 max_likelihood, x_max, y_max, z_max) ;

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
  return(max_likelihood) ;
}
static double
find_optimal_linear_xform(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source,
													MATRIX *m_L, MATRIX *m_origin,
													float min_angle, float max_angle,
													float min_scale, float max_scale,
													float min_trans, float max_trans,
													float angle_steps, float scale_steps, 
													float trans_steps,
													int nreductions,
													int rigid)
{
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
           *m_tmp2, *m_scale, *m_trans, *m_tmp3 = NULL ;
  double   x_angle, y_angle, z_angle, x_max_rot, y_max_rot, z_max_rot, 
           delta_rot, x_max_scale, y_max_scale, z_max_scale, delta_scale, 
           x_trans, delta_trans, y_trans, z_trans,
           mean_angle, x_scale, y_scale, z_scale, mean_scale, x_max_trans,
           y_max_trans, likelihood, max_likelihood, z_max_trans, mean_trans ;
  int      i ;
	MRI      *mri_target, *mri_source ;

	mri_target = vl_target->mri ;
	mri_source = vl_source->mri ;
  m_trans = MatrixIdentity(4, NULL) ;
  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max_trans = y_max_trans = z_max_trans = x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_likelihood = (*pf_likelihood)(vl_target, vl_source, m_L, 1.0) ;
  for (i = 0 ; i < nreductions ; i++)
  {
    delta_trans = (max_trans-min_trans) / (trans_steps-1) ;
		if (rigid)
		{
			max_scale = min_scale = 1 ;
			delta_scale = max_scale ;
		}
		else
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
											likelihood = (*pf_likelihood)(vl_target, vl_source, m_L_tmp, 1.0) ;
		      
											if (likelihood > max_likelihood)
											{
												max_likelihood = likelihood ;
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
      printf("  max likelihood = %2.4f @ R=(%2.3f,%2.3f,%2.3f),"
						 "S=(%2.3f,%2.3f,%2.3f), T=(%2.1f,%2.1f,%2.1f)\n", 
						 max_likelihood, DEGREES(x_max_rot), DEGREES(y_max_rot),
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
  return(max_likelihood) ;
}
#define NPARMS ((3*4)+1)
#ifdef TOL
#undef TOL
#endif
#define TOL 1e-12

static VOXEL_LIST *Gvl_target, *Gvl_source ;

static int
powell_minimize(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *mat)
{
	float *p, **xi, fret, fstart, intensity_scale;
	int   i, r, c, iter, diag ;

	p = vector(1, NPARMS) ;
	xi = matrix(1, NPARMS, 1, NPARMS) ;
	for (i = r = 1 ; r <= 3 ; r++)
	{
		for (c = 1 ; c <= 4 ; c++)
		{
			p[i++] = *MATRIX_RELT(mat, r, c) ;
		}
	}
	p[NPARMS] = 1.0 ;  // intensity scaling

	Gvl_target = vl_target ; Gvl_source = vl_source ;
	for (r = 1 ; r <= NPARMS ; r++)
	{
		for (c = 1 ; c <= NPARMS ; c++)
		{
			xi[r][c] = r == c ? 1 : 0 ;
		}
	}

	diag = Gdiag ; Gdiag |= DIAG_VERBOSE ;
	powell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
	Gdiag = diag ;
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
		for (i = r = 1 ; r <= 3 ; r++)
		{
			for (c = 1 ; c <= 4 ; c++)
			{
				*MATRIX_RELT(mat, r, c) = p[i++] ;
			}
		}
		intensity_scale = p[NPARMS] ;
		*MATRIX_RELT(mat, 4, 1) = 0.0 ; *MATRIX_RELT(mat, 4, 2) = 0.0 ; 
		*MATRIX_RELT(mat, 4, 3) = 0.0 ; *MATRIX_RELT(mat, 4, 4) = 1.0 ; 
		printf("%3.3d: best alignment at after powell: %2.3f (%d steps, iscale=%2.3f)\n", parms.start_t,fret, iter, intensity_scale) ;
		write_snapshot(vl_target->mri, vl_source->mri, mat, &parms, parms.start_t++,1,NULL);
	} while (fret < fstart) ;

	free_matrix(xi, 1, NPARMS, 1, NPARMS) ;
	free_vector(p, 1, NPARMS) ;
	return(NO_ERROR) ;
}

static float
compute_powell_sse(float *p)
{
	static MATRIX *mat = NULL ;
	float  error, intensity_scale ;
	int    i, r, c ;

	if (mat == NULL)
		mat = MatrixAlloc(4, 4, MATRIX_REAL) ;
	for (i = r = 1 ; r <= 3 ; r++)
	{
		for (c = 1 ; c <= 4 ; c++)
		{
			*MATRIX_RELT(mat, r, c) = p[i++] ;
		}
	}
	intensity_scale = p[NPARMS] ;
	*MATRIX_RELT(mat, 4, 1) = 0.0 ; *MATRIX_RELT(mat, 4, 2) = 0.0 ; 
	*MATRIX_RELT(mat, 4, 3) = 0.0 ; *MATRIX_RELT(mat, 4, 4) = 1.0 ; 
	error = -(*pf_likelihood)(Gvl_target, Gvl_source, mat, intensity_scale) ;
	return(error) ;
}

static int
write_snapshot(MRI *mri_target, MRI *mri_source, MATRIX *m_vox_xform, 
							 MORPH_PARMS *parms, int fno, int conform, char *in_fname)
{
	MRI *mri_aligned ;
	char fname[STRLEN] ;
	LTA  *lta ;

	if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
	{
		printf("source->target vox->vox transform:\n") ;
		MatrixPrint(stdout, m_vox_xform) ;
	}
	if (conform)
	{
		mri_aligned = MRIclone(mri_target, NULL) ;
		MRIlinearTransformInterp(mri_source, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
	}
	else
	{
		lta = LTAalloc(1, NULL) ;
		MatrixCopy(m_vox_xform, lta->xforms[0].m_L) ;
		//		mri_aligned = MRITransformedCenteredMatrix(mri_source_intensity, mri_target, m_vox_xform) ;
		mri_aligned = MRITransformedCenteredMatrix(mri_source, mri_target, m_vox_xform) ;
	}
	if (in_fname)
		sprintf(fname, "%s_%s", parms->base_name, in_fname) ;
	else
		sprintf(fname, "%s_%03d", parms->base_name, fno) ;
	MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
	if (in_fname)
		sprintf(fname, "%s_%s.mgz", parms->base_name, in_fname) ;
	else
		sprintf(fname, "%s_%03d.mgz", parms->base_name, fno) ;
	printf("writing snapshot to %s...\n", fname) ;
	MRIwrite(mri_aligned, fname) ;
	MRIfree(&mri_aligned) ;

	{
#if 0
		mri_aligned = MRIsrcTransformedCentered(mri_source, mri_target, m_vox_xform,
																						SAMPLE_NEAREST) ;
#else
		mri_aligned = MRITransformedCenteredMatrix(mri_source, mri_target, m_vox_xform) ;
#endif
		if (in_fname)
			sprintf(fname, "orig_%s_%s.mgz", parms->base_name, in_fname) ;
		else
			sprintf(fname, "orig_%s_%03d.mgz", parms->base_name, fno) ;
		printf("writing snapshot to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}

	return(NO_ERROR) ;
}
