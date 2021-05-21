/**
 * @brief program for computing a nonlinear alignment using a discrete cosine transform
 *
 * compute the discrete cosine transform to align two binary images, typically
 * initializing with an affine xform.
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


//
// mri_nl_align.c
//
// written by Bruce Fischl
// Mar, 2007
// 
// Warning: Do not edit the following four lines.  CVS maintains them.
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
#include "numerics.h"
#include "version.h"
#include "transform.h"
#include "fastmarching.h"
#include "gcamorph.h"
#include "voxlist.h"
#include "dct.h"
#include "hippo.h"

static MRI *Gmri_source, *Gmri_target, *Gmri_orig_source, *Gmri_aseg,
           *Gmri_orig_target ;

static void dfp_step_func(int itno, float sse, void *vparms, float *p) ;
static void compute_gradient(float *p, float *g) ;
static float compute_sse(float *p)  ;
static int    powell_minimize(DCT *dct, MRI *mri_source, MRI *mri_target) ;
static int    quasi_newton_minimize(DCT *dct,MRI *mri_source, MRI *mri_target, GCA_MORPH_PARMS *gmp) ;

static float smooth_intensities = -1.0 ;
static int morph_to = 0 ;
static int find_label = -1 ;
static float x_ras = 0.0 ;
static float y_ras = 0.0 ;
static float z_ras = 0.0 ;
static int mode_filters = 0 ;

static double sigma = 4 ;
static int upsample = 0 ;
static int apply_transform = 1 ;
static int use_powell = 0 ;

#define FLAGS_SOURCESPACE 1

static int flags =  0 ;
static int write_snapshot(DCT *dct, MRI *mri_source, MRI *mri_target, char *fname, int flags) ;
static DCT *find_optimal_dct(DCT *dct, MRI *mri_source, MRI *mri_target, 
                             int ncoef,int skip) ;
int MRImapRegionToTargetMRI(MRI *mri_src, MRI *mri_dst, MRI_REGION *box) ;

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;

static MRI_REGION Gtarget_box, Gsource_box ;

static char *source_intensity_fname = NULL ;
const char *Progname ;

static int skip = 4 ;

static int Gncoef = 5 ;

static TRANSFORM  *transform = NULL ;
static GCA_MORPH_PARMS mp ;


#define NONE  0
#define ANGIO 1
#define HIPPO 2
#define WM    3
#define LABEL 4

static int which = NONE ;
int
main(int argc, char *argv[])
{
	char         **av, *source_fname, *target_fname, *out_fname, fname[STRLEN] ;
  int          ac, nargs, i, new_transform = 0, orig_skip ;
	MRI          *mri_target, *mri_source, *mri_kernel, *mri_smooth_source,
               *mri_smooth_target ;
  Timer start ;
  int          msec, hours, minutes, seconds ;
	MATRIX       *m_L/*, *m_I*/ ;
	LTA          *lta ;
  DCT          *dct = NULL ;

  mp.npasses = skip+1 ;

  start.reset() ;
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

  if (argc < 3)
		usage_exit(1) ;

  mp.diag_morph_from_atlas = !morph_to ;

	source_fname = argv[1] ;
	target_fname = argv[2] ;
	out_fname = argv[3] ;
	printf("source = %s\ntarget = %s\noutput = %s\n", source_fname, target_fname,out_fname);
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(mp.base_name, fname) ;
	mri_source = MRIread(source_fname) ;
	if (!mri_source)
		ErrorExit(ERROR_NOFILE, "%s: could not read source label volume %s",Progname, source_fname) ;

	mri_target = MRIread(target_fname) ;
	if (!mri_target)
		ErrorExit(ERROR_NOFILE, "%s: could not read target label volume %s", Progname, target_fname) ;

  switch (which)
  {
  case HIPPO:
    Gmri_orig_source = MRIcopy(mri_source, NULL) ;
    mri_source = HIPPOestimateIntensityImage(mri_source, Gmri_aseg, mri_target, NULL) ;
    MRIwrite(mri_source, "hint.mgz") ;
    break;
  case NONE:
    Gmri_orig_source = MRIcopy(mri_source, NULL) ;
    Gmri_orig_target = MRIcopy(mri_target, NULL) ;
    break ;
  }

	MRIboundingBox(mri_source, 0, &Gsource_box) ;
	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
		MRIwrite(mri_source, "s.mgz") ;

	MRIboundingBox(mri_target, 0, &Gtarget_box) ;
	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
		MRIwrite(mri_target, "t.mgz") ;

  if (Gdiag & DIAG_WRITE && mp.write_iterations > 0)
  {
    int req = snprintf(fname, STRLEN,
		       "%s_target", mp.base_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwriteImageViews(mri_target, fname, IMAGE_SIZE) ;	
    
    req = snprintf(fname, STRLEN, "%s_target.mgz", mp.base_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwrite(mri_target, fname) ;
  }
  
  if (transform == NULL)
    transform = TransformAlloc(LINEAR_RAS_TO_RAS, NULL) ;
  
  
  if (transform->type == MORPH_3D_TYPE)  // initializing m3d from a linear transform
    {
      printf("using previously create gcam...\n") ;
      ErrorExit(ERROR_UNSUPPORTED, "%s: can't initialize DCT with GCAM", Progname) ;
    }
  
  new_transform = 1 ;
  lta = ((LTA *)(transform->xform)) ;
  m_L = MRIrasXformToVoxelXform(mri_source, mri_target, lta->xforms[0].m_L, NULL) ;
  MatrixFree(&lta->xforms[0].m_L) ;
  
  lta->xforms[0].m_L = m_L ;
  printf("initializing DCT with vox->vox matrix:\n") ;
  MatrixPrint(stdout, m_L) ;
  mri_source = MRITransformedCenteredMatrix(mri_source, mri_target, m_L) ;
  if (Gdiag & DIAG_WRITE && mp.write_iterations > 0)
    MRIwrite(mri_source, "st.mgz") ;

  orig_skip = skip ;
  do
  {
    printf("------------------- sigma = %2.4f -----------------\n",
           sigma) ;
    mri_kernel = MRIgaussian1d(sigma, -1) ;
    mri_smooth_source = MRIconvolveGaussian(mri_source, NULL, mri_kernel) ;
    mri_smooth_target = MRIconvolveGaussian(mri_target, NULL, mri_kernel) ;
    for (i = 0 ; i < mp.npasses ; i++) 
    {
      printf("------------- outer loop iteration %d of %d - skip %d, %d ---------------\n",
             i+1, mp.npasses,skip, skip) ;
      
      dct = find_optimal_dct(dct, mri_source, mri_target, Gncoef, skip) ;
      
      if (use_powell)
        powell_minimize(dct, mri_smooth_source, mri_smooth_target) ;
      else
        quasi_newton_minimize(dct, mri_smooth_source, mri_smooth_target, &mp) ;

      if (apply_transform) 
      {
        MRI *mri_aligned ;
        char   fname[STRLEN] ;

        FileNameRemoveExtension(out_fname, fname) ;
        strcat(fname, ".mgz") ;
        mri_aligned = DCTapply(dct, mri_source, NULL, NULL, SAMPLE_NEAREST) ;
        if (mode_filters > 0)
        {
          MRI *mri_filtered ;
          mri_filtered = MRImodeFilter(mri_aligned, NULL, mode_filters) ;
          MRIfree(&mri_aligned) ;
          mri_aligned = mri_filtered ;
        }
        printf("writing transformed output volume to %s...\n", fname) ;
        MRIwrite(mri_aligned, fname) ;
        MRIfree(&mri_aligned) ;
      }
      skip /= 2 ;
    }

    skip = orig_skip ; sigma /= 2 ; 
    MRIfree(&mri_kernel) ; MRIfree(&mri_smooth_source) ; 
    MRIfree(&mri_smooth_target) ;
  } while (sigma > 0.25);
  {
    FILE *fp = fopen(out_fname, "w") ;
    printf("writing DCT to %s...\n", out_fname) ;
    DCTdump(dct, fp) ;
    fclose(fp) ;
  }
	if (apply_transform)
	{
		MRI *mri_aligned ;
		char   fname[STRLEN] ;
			
		FileNameRemoveExtension(out_fname, fname) ;
		strcat(fname, ".mgz") ;
    mri_aligned = DCTapply(dct, mri_source, NULL, NULL, SAMPLE_NEAREST) ;
    if (mode_filters > 0)
    {
      MRI *mri_filtered ;
      mri_filtered = MRImodeFilter(mri_aligned, NULL, mode_filters) ;
      MRIfree(&mri_aligned) ;
      mri_aligned = mri_filtered ;
    }
		printf("writing transformed output volume to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}
  DCTfree(&dct) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
	hours = seconds / (60*60) ;
  minutes = (seconds/60) % 60 ;
  seconds = seconds % 60 ;
  printf("registration took %d hours, %d minutes and %d seconds.\n", 
				 hours, minutes, seconds) ;
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
  else if (!stricmp(option, "CJ"))
  {
    mp.constrain_jacobian = 1;
		mp.l_jacobian = 0 ;
		mp.ratio_thresh = .25 ;
		mp.noneg = False ;
    printf("constraining jacobian to be in [%2.2f %2.2f]\n",
					 mp.ratio_thresh, 1/mp.ratio_thresh) ;
  }
  else if (!stricmp(option, "neg"))
  {
		mp.noneg = False ;
    printf("allowing negative vertices during morph\n") ;
  }
  else if (!stricmp(option, "morph_to"))
  {
		morph_to = 1 ;
    printf("morphing to atlas...\n") ;
  }
  else if (!stricmp(option, "find_label"))
  {
		find_label = atoi(argv[2]) ;
		x_ras = atof(argv[3]) ;
		y_ras = atof(argv[4]) ;
		z_ras = atof(argv[5]) ;
		nargs = 4 ;
    printf("finding label %s (%d) at (%2.1f, %2.1f, %2.1f)\n",
					 cma_label_to_name(find_label), find_label, x_ras, y_ras,z_ras) ;
  }
  else if (!stricmp(option, "scale_smoothness"))
  {
		mp.scale_smoothness = atoi(argv[2]) ;
		mp.npasses = 2 ;
    printf("%sscaling smooothness coefficient (default=1), and setting npasses=%d\n",
					 mp.scale_smoothness ? "" : "not ", mp.npasses) ;
		nargs = 1 ;
  }
  else if (!stricmp(option, "MOMENTUM") || !stricmp(option, "FIXED"))
  {
    mp.integration_type = GCAM_INTEGRATE_FIXED ;
    printf("using optimal time-step integration\n") ;
  }
  else if (!stricmp(option, "view"))
	{
		Gsx = atoi(argv[2]) ;
		Gsy = atoi(argv[3]) ;
		Gsz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("viewing voxel (%d, %d, %d)\n", Gsx, Gsy, Gsz) ;
	}
  else if (!stricmp(option, "LEVELS"))
  {
    mp.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", mp.levels) ;
  }
  else if (!stricmp(option, "area"))
  {
    mp.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_area=%2.3f\n", mp.l_area) ;
  }
  else if (!stricmp(option, "tol"))
  {
    mp.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("using tol=%2.3f\n", mp.tol) ;
  }
  else if (!stricmp(option, "si"))
  {
    smooth_intensities = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing gcam intensities with sigma=%2.2f\n",smooth_intensities);
  }
  else if (!stricmp(option, "sigma"))
  {
    mp.sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("using sigma=%2.3f\n", mp.sigma) ;
  }
	else if (!stricmp(option, "rthresh"))
  {
		mp.ratio_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using compression ratio threshold = %2.3f...\n", mp.ratio_thresh) ;
  }
	else if (!stricmp(option, "dt"))
	{
		mp.dt = atof(argv[2]) ;
		nargs = 1 ;
		printf("using dt = %2.3f\n", mp.dt) ;
	}
	else if (!stricmp(option, "passes"))
	{
		mp.npasses = atoi(argv[2]) ;
		nargs = 1 ;
		printf("integrating in %d passes (default=3)\n", mp.npasses) ;
	}
	else if (!stricmp(option, "skip"))
	{
		skip = atoi(argv[2]);
		printf("skipping %d voxels in source data...\n", skip) ;
		nargs = 1 ;
	}
	else if (!stricmp(option, "hippo"))
	{
		which = HIPPO ;
    flags = FLAGS_SOURCESPACE ;
    Gmri_aseg = MRIread(argv[2]) ;
    if (Gmri_aseg == NULL)
      ErrorExit(ERROR_BADFILE, "%s: could not read aseg from %s", argv[2]) ;
    nargs = 1 ;
		printf("assuming source is hires hippo and dst is aseg volume\n") ;
	}
	else if (!stricmp(option, "wm"))
	{
		which = WM ;
		printf("assuming source and target are wm volumes\n") ;
	}
	else if (!stricmp(option, "none"))
	{
		which = NONE ;
		printf("making no assumptions about labels (not angio or hippo\n") ;
	}
	else if (!stricmp(option, "upsample"))
	{
		upsample = atoi(argv[2]) ;
		printf("upsampling the GCAM %d times\n", upsample) ;
		nargs = 1 ;
	}
	else switch (*option)
	{
  case 'M':
    mp.momentum = atof(argv[2]) ;
    nargs = 1 ;
    printf("momentum = %2.2f\n", mp.momentum) ;
    break ;
	case 'N':
    Gncoef = atoi(argv[2]) ;
		nargs = 1 ;
		printf("using %d coefficients\n", Gncoef);
		break ;
	case 'T':
		printf("reading transform from %s...\n", argv[2]) ;
		transform = TransformRead(argv[2]) ;
		if (transform == NULL)
			ErrorExit(ERROR_NOFILE,"%s: could not read transform from %s\n",Progname,argv[2]);
		nargs = 1 ;
		break ;
	case 'I':
		source_intensity_fname = argv[2] ;
		nargs = 1 ;
		printf("reading intensity image from %s for debugging...\n", source_intensity_fname) ;
		break ;
	case 'F':
		mode_filters = atoi(argv[2]) ;
		nargs = 1 ;
		printf("applying %d mode filters before writing out transformed volume\n", mode_filters) ;
		break ;
	case 'B':
    mp.l_binary = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_binary=%2.3f\n", mp.l_binary) ;
    break ;
  case 'J':
    mp.l_jacobian = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_jacobian=%2.3f\n", mp.l_jacobian) ;
    break ;
  case 'A':
		mp.navgs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing gradient with %d averages...\n", mp.navgs) ;
    break ;
  case 'K':
    printf("setting exp_k to %2.2f (default=%2.2f)\n",
           atof(argv[2]), mp.exp_k) ;
    mp.exp_k = atof(argv[2]) ;
    nargs = 1 ;
    break ;
	case 'W':
		mp.write_iterations = atoi(argv[2]) ;
		Gdiag |= DIAG_WRITE ;
		nargs = 1 ;
		printf("setting write iterations = %d\n", mp.write_iterations) ;
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
				 "usage: %s <source> <destination> <output xform>\n",
				 Progname) ;
	exit(ecode) ;
}






#define NCORNERS 8
int
MRImapRegionToTargetMRI(MRI *mri_src, MRI *mri_dst, MRI_REGION *box)
{
	VECTOR *v1, *v2 = NULL ;
	MATRIX *m_vox2vox ;
	int    xs[NCORNERS], ys[NCORNERS], zs[NCORNERS] ;
	int    xd[NCORNERS], yd[NCORNERS], zd[NCORNERS], xmin, xmax, ymin, ymax, zmin, zmax, i ;

	xs[0] = box->x ;              ys[0] = box->y ;           zs[0] = box->z ; 
	xs[1] = box->x+box->dx-1 ;    ys[1] = box->y ;           zs[1] = box->z ; 
	xs[2] = box->x ;              ys[2] = box->y+box->dy-1 ; zs[2] = box->z ; 
	xs[3] = box->x+box->dx-1 ;    ys[3] = box->y+box->dy-1 ; zs[3] = box->z ; 
	xs[4] = box->x ;              ys[4] = box->y ;           zs[4] = box->z+box->dz-1 ; 
	xs[5] = box->x+box->dx-1 ;    ys[5] = box->y ;           zs[5] = box->z+box->dz-1 ; 
	xs[6] = box->x ;              ys[6] = box->y+box->dy-1 ; zs[6] = box->z+box->dz-1 ; 
	xs[7] = box->x+box->dx-1 ;    ys[7] = box->y+box->dy-1 ; zs[7] = box->z+box->dz-1 ; 

	m_vox2vox = MRIgetVoxelToVoxelXform(mri_src, mri_dst) ;
	v1 = VectorAlloc(4, MATRIX_REAL) ;
	VECTOR_ELT(v1, 4) = 1.0 ;

	xmax = ymax = zmax = 0 ;
	xmin = mri_dst->width ; ymin = mri_dst->height ; zmin = mri_dst->depth ;
	for (i = 0 ; i < NCORNERS ; i++)
	{
		V3_X(v1) = xs[i] ; V3_Y(v1) = ys[i] ; V3_Z(v1) = zs[i] ; 
		v2 = MatrixMultiply(m_vox2vox, v1, v2) ;
		xd[i] = nint(V3_X(v2)) ; yd[i] = nint(V3_Y(v2)) ; zd[i] = nint(V3_Z(v2)) ;
		if (xd[i] > xmax)
			xmax = xd[i] ;
		if (xd[i] < xmin)
			xmin = xd[i] ;

		if (yd[i] > ymax)
			ymax = yd[i] ;
		if (yd[i] < ymin)
			ymin = yd[i] ;

		if (zd[i] > zmax)
			zmax = zd[i] ;
		if (zd[i] < zmin)
			zmin = zd[i] ;
	}

	box->x = xmin ;         box->y = ymin ;         box->z = zmin ;
	box->dx = xmax-xmin+1 ; box->dy = ymax-ymin+1 ; box->dz = zmax-zmin+1 ; 
	MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ;
	return(NO_ERROR) ;
}

static DCT *
find_optimal_dct(DCT *dct, MRI *mri_source, MRI *mri_target, int ncoef, int skip)
{
  Gmri_source = mri_source ; Gmri_target = mri_target ;
  if (dct == NULL)
  {
    dct = DCTalloc(ncoef, mri_source) ;
    MatrixWrite(dct->m_x_basis, "xdct.mat", "xdct") ;
    MatrixWrite(dct->m_x_basis, "ydct.mat", "ydct") ;
    MatrixWrite(dct->m_x_basis, "zdct.mat", "zdct") ;
    //    write_snapshot(dct, mri_source, mri_target, "test.mgz") ;
  }

  return(dct) ;
}

static int
write_snapshot(DCT *dct, MRI *mri_source, MRI *mri_target, char *fname, int flags)
{
  MRI *mri_aligned ;
  char *cp, tmpstr[STRLEN] ;

  if (flags == FLAGS_SOURCESPACE)
    mri_aligned = DCTapply(dct, mri_source, NULL, NULL, SAMPLE_NEAREST) ;
  else
    mri_aligned = DCTapply(dct, mri_source, mri_target, NULL, SAMPLE_NEAREST) ;
  printf("writing snapshot to %s\n", fname) ;
  MRIwrite(mri_aligned, fname) ;

  strcpy(tmpstr, fname) ;
  cp = strrchr(tmpstr, '.') ;
  if (cp)
    *cp = 0 ;

  MRIwriteImageViews(mri_aligned, tmpstr, IMAGE_SIZE) ;	
  MRIfree(&mri_aligned) ;
  return(NO_ERROR) ;
}




static int
powell_minimize(DCT *dct, MRI *mri_source, MRI *mri_target) 
{
  float *p, **xi, fret, fstart;
  int   i, k, iter, nparms, r, c, start_t = 0 ;

  nparms = 3*dct->ncoef;
#if SCALE_INTENSITIES
  nparms++ ;
#endif
  p = vector(1, nparms) ;
  xi = matrix(1, nparms, 1, nparms) ;
  for (i = k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_xk, k) ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_yk, k) ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_zk, k) ;
#if SCALE_INTENSITIES
  p[i++] = dct->b ;
#endif

  Gncoef = dct->ncoef ;
  for (r = 1 ; r <= nparms ; r++)
    for (c = 1 ; c <= nparms ; c++)
      xi[r][c] = r == c ? 1 : 0 ;

  // TODO:  powell(p, xi, nparms, TOL, &iter, &fret, compute_sse);
  fstart = compute_sse(p) ;
  printf("starting sse = %2.3f\n", fstart) ;
  OpenPowell(p, xi, nparms, TOL, &iter, &fret, compute_sse);
  fret = compute_sse(p) ;
  printf("ending sse = %2.3f\n", fret) ;
  start_t += iter ;
  do {
    for (r = 1 ; r <= nparms ; r++)
      for (c = 1 ; c <= nparms ; c++)
        xi[r][c] = r == c ? 1 : 0 ;

    fstart = fret ;
    // TODO:    powell(p, xi, nparms, TOL, &iter, &fret, compute_sse);
    OpenPowell(p, xi, nparms, TOL, &iter, &fret, compute_sse);
    start_t += iter ;
    for (i = k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_xk, k) = p[i++] ;
    for (k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_yk, k) = p[i++] ;
    for (k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_zk, k) = p[i++] ;
    dct->b = p[i++] ;

    printf("%4.4d: best alignment after powell: %2.3f (%d steps)\n",
           start_t,fret, iter) ;
  } while (fret < fstart) ;

  free_matrix(xi, 1, nparms, 1, nparms) ;
  free_vector(p, 1, nparms) ;
  return(NO_ERROR) ;
}

static void
dfp_step_func(int itno, float sse, void *vparms, float *p)
{
  MP     *parms = (MP *)vparms ;
  float  rms ;
  int    ino = itno + parms->start_t+1 ;

  rms = compute_sse(p) ;
#if SCALE_INTENSITIES
  printf("%4.4d: rms = %2.3f, b = %2.4f\n", ino, rms, p[3*Gncoef+1]) ;
#else
  printf("%4.4d: rms = %2.3f\n", ino, rms) ;
#endif
  if (Gdiag & DIAG_WRITE && ((ino % parms->write_iterations) == 0))
  {
    static DCT *dct = NULL ;
    char fname[STRLEN] ;
    int  i, k ;

    if (dct && dct->ncoef != Gncoef)
      DCTfree(&dct) ;
    if (dct == NULL)
      dct = DCTalloc(Gncoef, Gmri_source);
    for (i = k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_xk, k) = p[i++] ;
    for (k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_yk, k) = p[i++] ;
    for (k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_zk, k) = p[i++] ;
    sprintf(fname, "%s_%4.4d.mgz", parms->base_name, ino) ;
    write_snapshot(dct, Gmri_orig_source, Gmri_orig_target, fname, flags) ;
  }
}

static int
quasi_newton_minimize(DCT *dct, MRI *mri_source, MRI *mri_target, GCA_MORPH_PARMS *gmp)
{
  float *p, fret, fstart;
  int   i, k, iter, nparms, start_t = gmp->start_t ;
  MP    _parms, *parms = &_parms ;

  memset(parms, 0, sizeof(*parms)) ;
  strcpy(parms->base_name, gmp->base_name) ;
  parms->write_iterations = gmp->write_iterations ;
  parms->start_t = gmp->start_t ;
  nparms = 3*dct->ncoef ;
#if SCALE_INTENSITIES
  nparms++ ;
#endif
  p = vector(1, nparms) ;
  for (i = k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_xk, k) ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_yk, k) ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_zk, k) ;
#if SCALE_INTENSITIES
  p[i++] = dct->b ;
#endif

  Gncoef = dct->ncoef ;
  fret = fstart = compute_sse(p) ;
  printf("%4.4d: rms = %2.3f, b = %2.4f\n", 0, fstart, dct->b) ;
  if (parms->write_iterations > 0 && parms->start_t == 0)
  {
    char fname[STRLEN] ;
    
    sprintf(fname, "%s_%4.4d.mgz", parms->base_name, 0) ;
    write_snapshot(dct, Gmri_orig_source, Gmri_orig_target, fname, flags) ;
  }
    
  do
  {
    fstart = fret ;
    OpenDFPMin(p, nparms, TOL, &iter, &fret,
               compute_sse,
               compute_gradient, dfp_step_func, parms, NULL) ;
    fret = compute_sse(p) ;
    start_t += iter ; parms->start_t += iter ;
  } while ((fstart-fret)/fstart > TOL) ;

  for (i = k = 1 ; k <= dct->ncoef ; k++)
    VECTOR_ELT(dct->v_xk, k) = p[i++] ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    VECTOR_ELT(dct->v_yk, k) = p[i++] ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    VECTOR_ELT(dct->v_zk, k) = p[i++] ;
#if SCALE_INTENSITIES
  dct->b = p[i++] ;
#endif
  
  free_vector(p, 1, nparms) ;
  gmp->start_t = start_t ;
  return(NO_ERROR) ;
}
static float
compute_sse(float *p) 
{
  static DCT *dct = NULL ;
  double  x, y, z, error, sse, source, target, xd, yd, zd, rms;
  int     i, k, x2, y2, z2, nvox ;
  MRI        *mri_source = Gmri_source, *mri_target = Gmri_target ;
  static MATRIX     *m_vox2vox = NULL ;
  static VECTOR     *v1, *v2 ;

  if (dct && dct->ncoef != Gncoef)
    DCTfree(&dct) ;
  if (dct == NULL)
    dct = DCTalloc(Gncoef, Gmri_source);
  for (i = k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_xk, k) = p[i++] ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    VECTOR_ELT(dct->v_yk, k) = p[i++] ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    VECTOR_ELT(dct->v_zk, k) = p[i++] ;
#if SCALE_INTENSITIES
  dct->b = p[i++] ;
#endif

  if (m_vox2vox == NULL)
  {
    m_vox2vox = MRIgetVoxelToVoxelXform(mri_source, mri_target) ;
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
  }

  // iterate over source voxels and use transform to map  to target
  x2 = Gsource_box.x+Gsource_box.dx ; y2 = Gsource_box.y+Gsource_box.dy ;
  z2 = Gsource_box.z+Gsource_box.dz ;
  nvox = Gsource_box.dx * Gsource_box.dy * Gsource_box.dz / (pow(skip+1, 3.0)) ;
  nvox = 0 ; sse = 0.0 ;
  for (x = Gsource_box.x ; x < x2 ; x += skip+1)
  {
    for (y = Gsource_box.y ; y < y2 ; y += skip+1)
    {
      for (z = Gsource_box.z ; z < z2 ; z += skip+1)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        source = MRIgetVoxVal(mri_source, x, y, z, 0) ;
#if 0
        if (FZERO(source))
          continue ;
#endif
        xd = V3_X(v2) ; yd = V3_Y(v2) ; zd = V3_Z(v2) ;
        DCTtransformPoint(dct, x, y, z, &xd, &yd, &zd) ;
        V3_X(v1) = xd ; V3_Y(v1) = yd ; V3_Z(v1) = zd ; 
        MatrixMultiply(m_vox2vox, v1, v2) ; // map to target voxel space
        xd = V3_X(v2) ; yd = V3_Y(v2) ; zd = V3_Z(v2) ;
        MRIsampleVolume(mri_target, xd, yd, zd, &target) ;
        error = (target - dct->b*source) ; sse += (error*error) ;
        nvox++ ;
      }
    }
  }
  
  rms = sqrt(sse/(double)nvox) ;
  return(rms) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static void
compute_gradient(float *p, float *g)
{
  static DCT *dct = NULL ;
  double  error, dx, dy, dz, target, source, xd, yd, zd ;
  int    i, k, x, y, z, x2, y2, z2, nparms ;
  MRI        *mri_source = Gmri_source, *mri_target = Gmri_target ;
  static MATRIX     *m_vox2vox = NULL ;
  static VECTOR     *v1, *v2 ;

  if (m_vox2vox == NULL)
  {
    m_vox2vox = MRIgetVoxelToVoxelXform(mri_source, mri_target) ;
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
  }

  if (dct && dct->ncoef != Gncoef)
    DCTfree(&dct) ;
  if (dct == NULL)
    dct = DCTalloc(Gncoef, Gmri_source);
  nparms = 3*dct->ncoef ;
#if SCALE_INTENSITIES
  nparms++ ;
#endif
  for (i = k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_xk, k) = p[i++] ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    VECTOR_ELT(dct->v_yk, k) = p[i++] ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    VECTOR_ELT(dct->v_zk, k) = p[i++] ;
#if SCALE_INTENSITIES
  dct->b = p[i++];
#endif

  // compute gradient and return it in g
  memset(g, 0, sizeof(g[0])*nparms) ;

  // iterate over target voxels and use inverse transform to map back to source
  x2 = Gsource_box.x+Gsource_box.dx ; y2 = Gsource_box.y+Gsource_box.dy ;
  z2 = Gsource_box.z+Gsource_box.dz ;
  for (x = Gsource_box.x ; x < x2 ; x += skip+1)
  {
    V3_X(v1) = x ;
    for (y = Gsource_box.y ; y < y2 ; y += skip+1)
    {
      V3_Y(v1) = y ;
      for (z = Gsource_box.z ; z < z2 ; z += skip+1)
      {
        if (Gx == x && Gy == y && Gz == z)
          DiagBreak() ;

        source = MRIgetVoxVal(mri_source, x, y, z, 0) ;
#if 0
        if (FZERO(source))
          continue ;
#endif
        DCTtransformPoint(dct, x, y, z, &xd, &yd, &zd) ;
        V3_X(v1) = xd ; V3_Y(v1) = yd ; V3_Z(v1) = zd ; 
        MatrixMultiply(m_vox2vox, v1, v2) ; // map to target voxel space
        xd = V3_X(v2) ; yd = V3_Y(v2) ; zd = V3_Z(v2) ;
        MRIsampleVolume(mri_target, xd, yd, zd, &target) ;
        MRIsampleVolumeGradient(mri_target, xd, yd, zd, &dx, &dy, &dz) ;

        error = (target - dct->b*source) ;
        for (i = k = 1 ; k <= dct->ncoef ; k++)
          g[i++] += error * dx * *MATRIX_RELT(dct->m_x_basis, x+1, k) ;
        for (k = 1 ; k <= dct->ncoef ; k++)
          g[i++] += error * dy * *MATRIX_RELT(dct->m_y_basis, y+1, k) ;
        for (k = 1 ; k <= dct->ncoef ; k++)
          g[i++] += error * dz * *MATRIX_RELT(dct->m_z_basis, z+1, k) ;
#if SCALE_INTENSITIES
        g[i++] += error * (-source) ;
#endif
      }
    }
  }
  DiagBreak() ;
}
