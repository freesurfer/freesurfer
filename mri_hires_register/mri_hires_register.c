//
// mri_hires_register.c
//
// written by Bruce Fischl
// Nov. 9th ,2000
// 
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2005/01/25 17:07:17 $
// Revision       : $Revision: 1.1 $
//
////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

#define DEFAULT_MAX_STEPS 5
#define MAX_ANGLE       RADIANS(15)
#define MIN_ANGLE       RADIANS(2)

#define MAX_SCALE       2.0
#define MIN_SCALE       0.5

static int max_angles = DEFAULT_MAX_STEPS ;


static double MAX_TRANS = 30 ;

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;
#if 0
static MATRIX *pca_matrix(MATRIX *m_in_evectors, double in_means[3],
                         MATRIX *m_ref_evectors, double ref_means[3]) ;
static MATRIX *compute_pca(MRI *mri_in, MRI *mri_ref) ;
static int order_eigenvectors(MATRIX *m_src_evectors, MATRIX *m_dst_evectors);
#endif

static TRANSFORM *compute_optimal_transform(MRI *mri_aseg, MRI *mri_hires, 
																 INTEGRATION_PARMS *parms) ;
static int compute_overlap(MRI *mri_aseg, MRI *mri_hires, MATRIX *m_L) ;
static int find_optimal_translation(MRI *mri_aseg, MRI *mri_hires, 
																		MATRIX *m_L, float min_trans, 
																		float max_trans, 
																		float trans_steps, int nreductions);

static int find_optimal_linear_xform(MRI *mri_lowres, MRI *mri_hires,
																		 MATRIX *m_L, MATRIX *m_origin,
																		 float min_angle, float max_angle,
																		 float min_scale, float max_scale,
																		 float min_trans, float max_trans,
																		 float angle_steps, float scale_steps, 
																		 float trans_steps,
																		 int nreductions);
char *Progname ;
static int target_label = Left_Hippocampus ;

#define PAD 2


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
#define molecular_layer_subiculum	 12
#define lateral_ventricle		       13
#define molecular_layer_HP         14
#define fimbria                    15

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

int
main(int argc, char *argv[])
{
	char       **av, *hires_fname, *aseg_fname, *intensity_fname, *out_fname, fname[STRLEN] ;
  int        ac, nargs ;
	MRI        *mri_intensity, *mri_aseg, *mri_hires, *mri_aseg_orig, 
		         *mri_hires_orig, *mri_tmp ;
	MRI_REGION hires_box, lowres_box ;
	TRANSFORM  *transform ;

	parms.write_iterations = 1 ;
	parms.start_t = 0 ;
	

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
	mri_hires_orig = MRIread(hires_fname) ;
	if (!mri_hires_orig)
		ErrorExit(ERROR_NOFILE, "%s: could not read hires label volume %s",
							Progname, hires_fname) ;
	if (mri_hires_orig->type != MRI_UCHAR)
	{
		mri_tmp = MRIchangeType(mri_hires_orig, MRI_UCHAR, 0, 255, 1) ;
		MRIfree(&mri_hires_orig) ;
		mri_hires_orig = mri_tmp ;
	}
	mri_intensity = MRIread(intensity_fname) ;
	if (!mri_intensity)
		ErrorExit(ERROR_NOFILE, "%s: could not read intensity label volume %s",
							Progname, intensity_fname) ;
	mri_aseg_orig = MRIread(aseg_fname) ;
	if (!mri_aseg_orig)
		ErrorExit(ERROR_NOFILE, "%s: could not read aseg label volume %s",
							Progname, aseg_fname) ;

	/* crop lowres volume */
	mri_aseg = MRIclone(mri_aseg_orig, NULL) ;
	MRIcopyLabel(mri_aseg_orig, mri_aseg, target_label) ;
	MRIbinarize(mri_aseg, mri_aseg, 1, 0, 10) ;
	MRIwrite(mri_aseg, "b.mgz") ;
	MRIboundingBox(mri_aseg, 0, &lowres_box) ;
	printf("aseg bounding box (%d, %d, %d) + (%d, %d, %d)\n",
				 lowres_box.x, lowres_box.y, lowres_box.z, 
				 lowres_box.dx, lowres_box.dy, lowres_box.dz) ;
	lowres_box.x -= PAD ; lowres_box.y -= PAD ; lowres_box.z -= PAD ;
	lowres_box.dx += 2*PAD ; lowres_box.dy += 2*PAD ; lowres_box.dz += 2*PAD ;
	mri_tmp = MRIextractRegion(mri_aseg, NULL, &lowres_box) ;
	MRIfree(&mri_aseg) ; mri_aseg = mri_tmp ;
	mri_tmp = MRIextractRegion(mri_intensity, NULL, &lowres_box) ;
	MRIfree(&mri_intensity) ; mri_intensity = mri_tmp ;
	MRIwrite(mri_aseg, "a.mgz") ; MRIwrite(mri_intensity, "i.mgz") ;
  if (Gdiag & DIAG_WRITE && parms.write_iterations > 0)
  {
    MRIwriteImageViews(mri_aseg, "target", IMAGE_SIZE) ;
	}

	/* crop hires volume */
	mri_hires = MRIclone(mri_hires_orig, NULL) ;
	setDirectionCosine(mri_hires_orig, MRI_CORONAL) ;
	MRIbinarize(mri_hires_orig, mri_hires, 1, 0, 10) ;
	MRIwrite(mri_hires, "hb.mgz") ;
	MRIboundingBox(mri_hires, 0, &hires_box) ;
	printf("hires bounding box (%d, %d, %d) + (%d, %d, %d)\n",
				 hires_box.x, hires_box.y, hires_box.z, hires_box.dx, 
				 hires_box.dy, hires_box.dz) ;
	mri_tmp = MRIextractRegion(mri_hires, NULL, &hires_box) ;
	MRIfree(&mri_hires) ; mri_hires = mri_tmp ;
	MRIwrite(mri_hires, "h.mgz") ;

	transform = compute_optimal_transform(mri_aseg, mri_hires, &parms) ;

	printf("initial matrix:\n") ;
	MatrixPrint(stdout, ((LTA *)(transform->xform))->xforms[0].m_L) ;
#if 0
  if (Gdiag & DIAG_WRITE && parms.write_iterations > 0)
  {
    MRI *mri_aligned ;

    mri_aligned = MRIclone(mri_aseg, NULL) ;
    MRIlinearTransformInterp(mri_hires, mri_aligned, 
														 ((LTA *)(transform->xform))->xforms[0].m_L, 
														 SAMPLE_NEAREST);
    MRIwriteImageViews(mri_aligned, "after_pca", IMAGE_SIZE) ;
		MRIwrite(mri_aligned, "pca_aligned.mgz") ;
    MRIfree(&mri_aligned) ;
  }
#endif

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
		Gx = Gsx = atoi(argv[2]) ;
		Gy = Gsy = atoi(argv[3]) ;
		Gz = Gsz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
	}
	else switch (*option)
	{
	case 'W':
		parms.write_iterations = atoi(argv[2]) ;
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

#if 0
static MATRIX *
compute_pca(MRI *mri_in, MRI *mri_ref)
{
  int    row, col, i ; 
  float  dot ;
  MATRIX *m_ref_evectors = NULL, *m_in_evectors = NULL ;
  float  in_evalues[3], ref_evalues[3] ;
  double  ref_means[3], in_means[3] ;

  if (!m_ref_evectors)
    m_ref_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;
  if (!m_in_evectors)
    m_in_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;

	MRIbinaryPrincipleComponents(mri_ref, m_ref_evectors, ref_evalues, 
															 ref_means, 1);
	MRIbinaryPrincipleComponents(mri_in,m_in_evectors,in_evalues,in_means,
															 1);

  order_eigenvectors(m_in_evectors, m_in_evectors) ;
  order_eigenvectors(m_ref_evectors, m_ref_evectors) ;

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
      fprintf(stderr, "WARNING: mirror image detected in eigenvector #%d\n",
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

static int
order_eigenvectors(MATRIX *m_src_evectors, MATRIX *m_dst_evectors)
{
  int    row, col, xcol, ycol, zcol ;
  double mx ;

  if (m_src_evectors == m_dst_evectors)
    m_src_evectors = MatrixCopy(m_src_evectors, NULL) ;

  /* find columx with smallest dot product with unit x vector */
  mx = fabs(*MATRIX_RELT(m_src_evectors, 1, 1)) ; xcol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (fabs(*MATRIX_RELT(m_src_evectors, 1, col)) > mx)
    {
      xcol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 1, col)) ;
    }

  mx = fabs(*MATRIX_RELT(m_src_evectors, 2, 1)) ; ycol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (*MATRIX_RELT(m_src_evectors, 2, col) > mx)
    {
      ycol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 2, col)) ;
    }

  mx = fabs(*MATRIX_RELT(m_src_evectors, 3, 1)) ; zcol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (fabs(*MATRIX_RELT(m_src_evectors, 3, col)) > mx)
    {
      zcol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 3, col)) ;
    }

  for (row = 1 ; row <= 3 ; row++)
  {
    *MATRIX_RELT(m_dst_evectors,row,1) = *MATRIX_RELT(m_src_evectors,row,xcol);
    *MATRIX_RELT(m_dst_evectors,row,2) = *MATRIX_RELT(m_src_evectors,row,ycol);
    *MATRIX_RELT(m_dst_evectors,row,3) = *MATRIX_RELT(m_src_evectors,row,zcol);
  }
  return(NO_ERROR) ;
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

#if 0
#define MAX_X_ANGLE  (RADIANS(35))
#define MAX_Y_ANGLE  (RADIANS(15))
#define MAX_Z_ANGLE  (RADIANS(15))
  if (fabs(x_angle) > MAX_X_ANGLE || fabs(y_angle) > MAX_Y_ANGLE || 
      fabs(z_angle) > MAX_Z_ANGLE)
#endif
  {
    MATRIX *m_I ;

    /*    MatrixFree(&m_in_T) ; MatrixFree(&mRot) ;*/
    fprintf(stderr, 
       "eigenvector swap detected (%2.0f, %2.0f, %2.0f): ignoring rotational PCA...\n",
            DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;

    m_I = MatrixIdentity(3, NULL) ;
    MatrixCopy(m_I, mRot) ;
    MatrixFree(&m_I) ;
    x_angle = y_angle = z_angle = 0.0 ;
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
  dx = -ref_means[0] ; dy = -ref_means[1] ; dz = -ref_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ; *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ; *MATRIX_RELT(m_T, 4, 4) = 1 ;
  m_tmp = MatrixMultiply(m_R, m_T, NULL) ;
  *MATRIX_RELT(m_T, 1, 4) = -dx ; *MATRIX_RELT(m_T, 2, 4) = -dy ;
  *MATRIX_RELT(m_T, 3, 4) = -dz ; 
  MatrixMultiply(m_T, m_tmp, m_R) ;

  /* now apply translation to take in centroid to ref centroid */
  dx = ref_means[0] - in_means[0] ; dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ; *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ; *MATRIX_RELT(m_T, 4, 4) = 1 ;

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
  MatrixFree(&m_R) ; MatrixFree(&m_T) ;

  MatrixFree(&mRot) ;
  VectorFree(&mOrigin) ;
  return(m_L) ;
}
#endif

static TRANSFORM *
compute_optimal_transform(MRI *mri_lowres, MRI *mri_hires, 
													INTEGRATION_PARMS *parms)
{
	TRANSFORM *transform ;
	MATRIX    *m_scale, *m_L, *m_tmp, *m_origin, *m_inv_origin ;
	double    hires_cent[3], lowres_cent[3], dx, dy, dz, scale,
            min_scale, max_scale, max_angle, scale_steps, mean,
            angle_steps, delta, min_search_scale ;
	int       initial_overlap, niter, nscales, good_step, done,
            old_max_overlap, max_overlap ;

#define MIN_SEARCH_SCALE 0.1
  min_search_scale = MIN_SEARCH_SCALE ;
	transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
	m_scale = MatrixIdentity(4, NULL) ;
	m_L = ((LTA *)(transform->xform))->xforms[0].m_L ;
	*MATRIX_RELT(m_scale, 1,1) = mri_hires->xsize / mri_lowres->xsize ;
	*MATRIX_RELT(m_scale, 2,2) = mri_hires->xsize / mri_lowres->ysize ;
	*MATRIX_RELT(m_scale, 3,3) = mri_hires->zsize / mri_lowres->zsize ;

  m_origin = MatrixIdentity(4, NULL) ;

  MRIcenterOfMass(mri_hires, hires_cent, 0) ;
  *MATRIX_RELT(m_origin, 1, 4) = hires_cent[0] ; 
  *MATRIX_RELT(m_origin, 2, 4) = hires_cent[1] ;
  *MATRIX_RELT(m_origin, 3, 4) = hires_cent[2] ; 
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;
	m_inv_origin = MatrixInverse(m_origin, NULL) ;

  MRIcenterOfMass(mri_lowres, lowres_cent, 0) ;
	dx = hires_cent[0] - lowres_cent[0] ; dy = hires_cent[1] - lowres_cent[1] ;
	dz = hires_cent[2] - lowres_cent[2] ;
	*MATRIX_RELT(m_L, 1, 4) = dx ; *MATRIX_RELT(m_L, 2, 4) = dy ; *MATRIX_RELT(m_L, 3, 4) = dz ;
	
	m_tmp = MatrixMultiply(m_L, m_inv_origin, NULL) ;
	MatrixMultiply(m_scale, m_tmp, m_L) ;
	MatrixMultiply(m_origin, m_L, m_tmp) ;
	MatrixCopy(m_tmp, m_L) ;

	initial_overlap = compute_overlap(mri_lowres, mri_hires, m_L) ;
	printf("initial overlap = %d...\n", initial_overlap) ;
	max_overlap = find_optimal_translation(mri_lowres, mri_hires, m_L, 
																					-100, 100, 3, 4) ;


  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
  {
    MRI *mri_aligned ;
		char fname[STRLEN] ;

    mri_aligned = MRIclone(mri_lowres, NULL) ;
    MRIlinearTransformInterp(mri_hires, mri_aligned, 
														 ((LTA *)(transform->xform))->xforms[0].m_L, 
														 SAMPLE_NEAREST);
		sprintf(fname, "%s%03d", parms->base_name, parms->start_t+1) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
		sprintf(fname, "%s%03d.mgz", parms->base_name, parms->start_t+1) ;
		printf("writing snapshot to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
    MRIfree(&mri_aligned) ;
  }
  max_angle = MAX_ANGLE ; angle_steps = max_angles ;
  max_scale = MAX_SCALE ; min_scale = MIN_SCALE ; scale_steps = max_angles ;
#define MIN_SCALES 2
  /////////////////////////// loop here ////////////////////////////////////////////
  niter = 1 ; nscales = 0 ; scale = 1.0 ; good_step = 0 ; done = 0 ;
  do
  {
    old_max_overlap = max_overlap ;
    printf("****************************************\n");
    printf("Nine parameter search.  iteration %d nscales = %d ...\n", niter, nscales);
    printf("****************************************\n");
    max_overlap = find_optimal_linear_xform(mri_lowres, mri_hires, 
																						m_L, m_origin,
																						-RADIANS(2*scale),
																						RADIANS(2*scale),
																						1-.25/16.0*scale, 
																						1+.25/16.0*scale, 
																						-scale/16.0*MAX_TRANS, 
																						scale/16.0*MAX_TRANS,
																						3, 3, 3, 2);
    
    if (parms->write_iterations != 0)
    {
      char fname[STRLEN] ;
      MRI *mri_aligned ;
      
      mri_aligned = MRIlinearTransform(mri_hires, NULL, m_L) ;
      sprintf(fname, "%s%03d", parms->base_name, parms->start_t+niter+1) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s%03d.mgz", parms->base_name, parms->start_t+niter+1) ;
			printf("writing snapshot to %s...\n", fname) ;
			MRIwrite(mri_aligned, fname) ;

      MRIfree(&mri_aligned) ;
    }
    printf("Result so far: scale %2.3f: max_overlap=%d, old_max_overlap =%d\n",
						scale,max_overlap, old_max_overlap) ;
    MatrixPrint(stderr, m_L);
    /* search a finer nbhd (if do-while continues) */
    if ((max_overlap <= old_max_overlap)) /* couldn't take a step */
    {
      if (good_step)
      {
				scale *= 0.25 ;
				if (scale < min_search_scale)
					break ;
				mean = (max_scale + min_scale)/2 ;
				delta = (max_scale - min_scale)/2 ;
				max_scale = 1.0 + delta*scale ;
				min_scale = 1.0 - delta*scale ;
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
  } while (nscales++ < MIN_SCALES || (done == FALSE)) ;


	/*	m_L = compute_pca(mri_hires, mri_lowres) ;*/
	MatrixFree(&m_scale) ; MatrixFree(&m_tmp) ; MatrixFree(&m_origin) ;
	MatrixFree(&m_inv_origin) ;
	return(transform) ;
}

static int
compute_overlap(MRI *mri_lowres_orig, MRI *mri_hires_orig, MATRIX *m_L)
{
	int     overlap, x, y, z, width, height, depth, xd, yd, zd, label,
          hwidth, hheight, hdepth ;
	MATRIX  *m_inv ;
	VECTOR  *v1, *v2 ;
	MRI     *mri_lowres, *mri_hires ;

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
	mri_hires = MRIcopy(mri_hires_orig, NULL) ;
	mri_lowres = MRIcopy(mri_lowres_orig, NULL) ;

	/* first go through lowres volume and for every voxel that is on in it,
		 map it to the hires, and if the hires is on, add one to the overlap
	*/
	for (x = 0 ; x < width ; x++)
	{
		for (y = 0 ; y < height ; y++)
		{
			for (z = 0 ; z < depth ; z++)
			{
				if (MRIvox(mri_lowres, x, y, z) == 0)
					continue ;

				V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
				MatrixMultiply(m_inv, v1, v2) ;
				xd = nint(V3_X(v2)) ; yd = nint(V3_Y(v2)) ; zd = nint(V3_Z(v2)) ; 
				if (xd >= 0 && xd < hwidth &&
						yd >= 0 && yd < hheight &&
						zd >= 0 && zd < hdepth)
				{
					label = MRIvox(mri_hires, xd, yd, zd) ;
					if (label > 0)
					{
						overlap++ ;
						MRIvox(mri_hires, xd, yd, zd) = 0 ; /* only count it once */
					}
				}
			}
		}
	}

	/* first go through hires volume and for every voxel that is on in it,
		 map it to the lowres, and if the lowres hasn't been counted yet, count it.
	*/
	for (x = 0 ; x < hwidth ; x++)
	{
		for (y = 0 ; y < hheight ; y++)
		{
			for (z = 0 ; z < hdepth ; z++)
			{
				if (MRIvox(mri_hires_orig, x, y, z) == 0)
					continue ;

				V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
				MatrixMultiply(m_L, v1, v2) ;
				xd = nint(V3_X(v2)) ; yd = nint(V3_Y(v2)) ; zd = nint(V3_Z(v2)) ; 
				if (xd >= 0 && xd < width &&
						yd >= 0 && yd < height &&
						zd >= 0 && zd < depth)
				{
					label = MRIvox(mri_lowres, xd, yd, zd) ;
					if (label)
					{
						MRIvox(mri_lowres, xd, yd, zd) = 0 ;
						overlap++ ;
					}
				}
			}
		}
	}

	MRIfree(&mri_lowres) ; MRIfree(&mri_hires) ;
	MatrixFree(&m_inv) ; VectorFree(&v1) ; VectorFree(&v2) ;
	return(overlap) ;
}
static int
find_optimal_translation(MRI *mri_lowres, MRI *mri_hires, 
                         MATRIX *m_L, float min_trans, float max_trans, 
                         float trans_steps, int nreductions)
{
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta, 
           mean_trans ;
  int      overlap, max_overlap, i ;

  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_overlap = compute_overlap(mri_lowres, mri_hires, m_L) ;

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
          overlap = compute_overlap(mri_lowres, mri_hires, m_L_tmp) ;
          if (overlap > max_overlap)
          {
            max_overlap = overlap ;
            x_max = x_trans ; y_max = y_trans ; z_max = z_trans ;
#if 1
            printf("new max overlap %d found at (%2.1f, %2.1f, %2.1f)\n",
									 max_overlap, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
      printf(
						 "max overlap = %d @ (%2.1f, %2.1f, %2.1f)\n", 
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
static int
find_optimal_linear_xform(MRI *mri_lowres, MRI *mri_hires,
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
           y_max_trans, z_max_trans, mean_trans ;
  int      i, overlap, max_overlap ;

  m_trans = MatrixIdentity(4, NULL) ;
  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max_trans = y_max_trans = z_max_trans = x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_overlap = compute_overlap(mri_lowres, mri_hires, m_L) ;
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
											overlap = compute_overlap(mri_lowres, mri_hires, m_L_tmp) ;
		      
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
      printf("  max overlap = %d @ R=(%2.3f,%2.3f,%2.3f),"
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
