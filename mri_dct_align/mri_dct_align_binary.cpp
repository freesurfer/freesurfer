/**
 * @brief computes a nonlinear alignment using a discrete cosine transform
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
#include "voxlist.h"

#define NONMAX 0
#define PAD      10

static VOXEL_LIST *Gvl_source, *Gvl_target ;
static MRI *Gmri_source;
static int PADVOX = 1 ;

static float compute_powell_sse(float *p)  ;
static int    powell_minimize(VOXEL_LIST *vl_target,
                              VOXEL_LIST *vl_source,
                              DCT *dct,
                              MRI *mri_orig_source) ;
static int binary_label = 128 ;
static int check_angio_labels(MRI *mri_source, MRI *mri_target) ;
static double compute_recursive_optimum(DCT *dct, double scale, 
                                        VOXEL_LIST *vl_source, 
                                        VOXEL_LIST *vl_target, 
                                        int which_coord, int k,
                                        double best_overlap);

static float smooth_intensities = -1.0 ;
static int morph_to = 0 ;
static int find_label = -1 ;
static float x_ras = 0.0 ;
static float y_ras = 0.0 ;
static float z_ras = 0.0 ;
static int mode_filters = 0 ;

static int upsample = 0 ;
static int apply_transform = 1 ;

static double compute_global_dct_optimum(DCT *dct, MRI *mri_source, 
                                         MRI *mri_target, 
                                         VOXEL_LIST *vl_source, 
                                         VOXEL_LIST *vl_target) ;
static double compute_distance_transform_sse(VOXEL_LIST *vl_target,
                                             VOXEL_LIST *vl_source,DCT *dct) ;
static double compute_overlap(VOXEL_LIST *vl_target,
                              VOXEL_LIST *vl_source,
                              DCT *dct) ;
static double (*pf_overlap)(VOXEL_LIST *vl_target,
                            VOXEL_LIST *vl_source,
                            DCT *dct) = compute_overlap ;

static int write_snapshot(DCT *dct, MRI *mri_source, MRI *mri_target, const char *fname) ;
static DCT *find_optimal_dct(DCT *dct, MRI *mri_source, MRI *mri_target, 
                             VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, 
                             int ncoef,int skip) ;
int MRImapRegionToTargetMRI(MRI *mri_src, MRI *mri_dst, MRI_REGION *box) ;

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;


static char *source_intensity_fname = NULL ;
const char *Progname ;
static int target_label = 128 ;

static int skip = 2 ;
static double distance = 1.0 ;

static int Gncoef = 10 ;

static int non_artery_labels[] =
{
	Left_Common_IliacV,
	Right_Common_IliacV,
	Left_External_IliacV,
	Right_External_IliacV,
	Left_Internal_IliacV,
	Right_Internal_IliacV,
	Left_ObturatorV,
	Right_ObturatorV,
	Left_Internal_PudendalV,
	Right_Internal_PudendalV,
	Pos_Lymph,
	Neg_Lymph
} ;
#define NUM_NON_ARTERY_LABELS  (sizeof(non_artery_labels) / sizeof(non_artery_labels[0]))

#if 0
static int non_hippo_labels[] =
{
	entorhinal_cortex,
	Amygdala,
	Cerebral_White_Matter,
	Cerebral_Cortex,
	lateral_ventricle,
	Inf_Lat_Vent,
	Left_Cerebral_Cortex,
	Right_Cerebral_Cortex,
	Left_Cerebral_White_Matter,
	Right_Cerebral_White_Matter,
	Left_Inf_Lat_Vent,
	Right_Inf_Lat_Vent,
	Right_Lateral_Ventricle,
	Left_Lateral_Ventricle,
	Left_Thalamus_Proper,
	Right_Thalamus_Proper,
	Left_Thalamus,
	Right_Thalamus,
	Left_choroid_plexus,
	Right_choroid_plexus,
	Left_Amygdala,
	Right_Amygdala,
	Right_Pallidum,
	Left_Pallidum,
	Right_Putamen,
	Left_Putamen,
	Right_VentralDC,
	Left_VentralDC,
	Brain_Stem
} ;

#define NUM_NON_HIPPO_LABELS  (sizeof(non_hippo_labels) / sizeof(non_hippo_labels[0]))
#endif

static int target_aseg_label = Right_Hippocampus ;

static TRANSFORM  *transform = NULL ;
static GCA_MORPH_PARMS mp ;

#define NONE  0
#define ANGIO 1
#define HIPPO 2
#define WM    3
#define LABEL 4

static int which = ANGIO ;

int
main(int argc, char *argv[])
{
	char         **av, *source_fname, *target_fname, *out_fname, fname[STRLEN] ;
  int          ac, nargs, i, new_transform = 0, pad ;
	MRI          *mri_target, *mri_source, *mri_tmp, *mri_orig_source, *mri_orig_target ;
	MRI          *mri_dist_target = NULL, *mri_dist_source = NULL ;
	MRI_REGION   box ;
  Timer start ;
  int          msec, hours, minutes, seconds, label, j ;
	MATRIX       *m_L/*, *m_I*/ ;
	LTA          *lta ;
  DCT          *dct = NULL ;
  VOXEL_LIST    *vl_target, *vl_source ;


  mp.npasses = 2 ;

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
	mri_orig_target = MRIcopy(mri_target, NULL) ;

	// crop input volumes
	if (which == WM)
	{
		MRI *mri_tmp ;
		mri_tmp = MRIclone(mri_source, NULL) ;
		MRIcopyLabel(mri_source, mri_tmp, Left_Cerebral_White_Matter) ;
		MRIcopyLabel(mri_source, mri_tmp, Right_Cerebral_White_Matter) ;
		MRIcopyLabel(mri_source, mri_tmp, Left_Cerebellum_White_Matter) ;
		MRIcopyLabel(mri_source, mri_tmp, Right_Cerebellum_White_Matter) ;
		MRIcopyLabel(mri_source, mri_tmp, Brain_Stem) ;
		MRIfree(&mri_source) ; mri_source = mri_tmp ;
		MRIeraseBorders(mri_source, 1) ;
	}
	else if (which == LABEL)
	{
		MRI *mri_tmp ;
		mri_tmp = MRIclone(mri_source, NULL) ;
		MRIcopyLabel(mri_source, mri_tmp, target_aseg_label) ;
		MRIfree(&mri_source) ; mri_source = mri_tmp ;
	}
	MRIboundingBox(mri_source, 0, &box) ;
	pad = PADVOX ;
	printf("padding source with %d voxels...\n", pad) ;
	mri_tmp = MRIextractRegionAndPad(mri_source, NULL, &box, pad) ;
	MRIfree(&mri_source) ;
	mri_source = mri_tmp ;
	//	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
		MRIwrite(mri_source, "s.mgz") ;
	mri_orig_source = MRIcopy(mri_source, NULL) ;

	if (which == HIPPO || which == LABEL) // just copy label out of the target 
	{
		mri_tmp = MRIclone(mri_target, NULL) ; 
		MRIcopyLabel(mri_target, mri_tmp, target_aseg_label) ;
		MRIfree(&mri_target) ; mri_target = mri_tmp ;
		mp.target_label = target_label = target_aseg_label ;
	}
	else if (which == ANGIO)
	{
    MRIbinarize(mri_target, mri_target, 1, 0, binary_label) ;
    MRIbinarize(mri_source, mri_source, 1, 0, binary_label) ;
		mp.target_label = target_label ;
	}
	else if (which == WM)
	{
		MRI *mri_tmp ;
		mri_tmp = MRIclone(mri_target, NULL) ;
		MRIcopyLabel(mri_target, mri_tmp, Left_Cerebral_White_Matter) ;
		MRIcopyLabel(mri_target, mri_tmp, Right_Cerebral_White_Matter) ;
		MRIcopyLabel(mri_target, mri_tmp, Left_Cerebellum_White_Matter) ;
		MRIcopyLabel(mri_target, mri_tmp, Right_Cerebellum_White_Matter) ;
		MRIcopyLabel(mri_target, mri_tmp, Brain_Stem) ;
		MRIfree(&mri_target) ; mri_target = mri_tmp ;
		MRIeraseBorders(mri_target, 1) ;
	}
	MRIboundingBox(mri_target, 0, &box) ;
	pad = 0 ; // don't need to pad target volume
	if (!FZERO(mp.l_area_intensity) || !FZERO(mp.l_log_likelihood)  || !FZERO(mp.l_likelihood))
	{
		// need to provide some context for area-intensity case
#define EXTRA 5
		MRIcopy(mri_orig_target, mri_target) ;  // restore labels
		box.x -= EXTRA ; box.y -= EXTRA ; box.z -= EXTRA ; 
		box.dx += EXTRA*2 ; box.dy += EXTRA*2 ; box.dz += EXTRA*2 ; 
		MRIcropBoundingBox(mri_target, &box) ;
	}

  MRIbinarize(mri_target, mri_target, 1, 0, binary_label) ;
  MRIbinarize(mri_source, mri_source, 1, 0, binary_label) ;

	mri_tmp = MRIextractRegionAndPad(mri_target, NULL, &box, pad) ;
	MRIfree(&mri_target) ;mri_target = mri_tmp ;
	mri_tmp = MRIextractRegionAndPad(mri_orig_target, NULL, &box, pad) ;
	MRIfree(&mri_orig_target) ;
	mri_orig_target = mri_tmp ;
	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
		MRIwrite(mri_target, "t.mgz") ;

  if (pf_overlap == compute_distance_transform_sse)
  {
    printf("creating distance transforms...\n") ;
    
    mri_dist_source = MRIdistanceTransform(mri_source, NULL, binary_label, -1, DTRANS_MODE_SIGNED, NULL);
    mri_dist_target = MRIdistanceTransform(mri_target, NULL, binary_label, -1, DTRANS_MODE_SIGNED, NULL);
    MRIscalarMul(mri_dist_source, mri_dist_source, -1) ;
    MRIscalarMul(mri_dist_target, mri_dist_target, -1) ;
  }
	
	if (which == ANGIO)
	{
		for (i = 0 ; i < NUM_NON_ARTERY_LABELS ; i++)
		{
			label = non_artery_labels[i] ;
			MRIreplaceValues(mri_source, mri_source, label, 0) ;
			MRIreplaceValues(mri_target, mri_target, label, 0) ;
		}
		check_angio_labels(mri_source, mri_target) ;
	}
	else if (which == HIPPO)
	{
#if 0
		for (i = 0 ; i < NUM_NON_HIPPO_LABELS ; i++)
		{
			label = non_hippo_labels[i] ;
			MRIreplaceValues(mri_source, mri_source, label, 0) ;
		}
#endif
	}

	if (which == WM)  // debugging
	{
		mri_orig_target = mri_target ; mri_orig_source = mri_source ;
		mp.target_label = target_label ;
	}
  if (Gdiag & DIAG_WRITE && mp.write_iterations > 0)
  {
		sprintf(fname, "%s_target", mp.base_name) ;
    MRIwriteImageViews(mri_orig_target, fname, IMAGE_SIZE) ;	
		sprintf(fname, "%s_target.mgz", mp.base_name) ;
		MRIwrite(mri_orig_target, fname) ;
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
  mri_orig_source = MRITransformedCenteredMatrix(mri_orig_source, mri_target, m_L) ;
  MRIwrite(mri_orig_source, "st.mgz") ;

	if (mp.write_iterations != 0)
	{
		char fname[STRLEN] ;

		sprintf(fname, "%s_target.mgz", mp.base_name) ;
    printf("writing target volume to %s...\n", fname) ;
		if (mp.diag_morph_from_atlas)
			MRIwrite(mri_orig_target, fname) ;
		else
		{
      MRI    *mri_tmp ;
      mri_tmp = MRITransformedCenteredMatrix(mri_source, mri_target, m_L) ;
			MRIwrite(mri_tmp, fname) ;
      MRIfree(&mri_tmp) ;
		}
	}

  for (i = 0 ; i < mp.npasses ; i++) 
  {
    printf("------------- outer loop iteration %d of %d - skip %d, %d ---------------\n",
           i+1, mp.npasses,skip, skip) ;

    vl_target = VLSTcreate(mri_target, 1, MAX_LABEL+1, NULL, skip, 0);
    printf("%d voxels in target list\n", vl_target->nvox) ;
    vl_target->mri2 = mri_dist_target ;
      
    vl_source = VLSTcreate(mri_source, 1, MAX_LABEL+1, NULL, skip, 0) ;
    printf("%d voxels in source list\n", vl_source->nvox) ;
    vl_source->mri2 = mri_dist_source ;
    
    dct = find_optimal_dct(dct, mri_source, mri_target, vl_source, vl_target, Gncoef, skip) ;
    
    VLSTfree(&vl_source) ;
    vl_source = VLSTcreate(mri_source, 1, 255, NULL, skip/4, 0) ;
    vl_source->mri2 = mri_dist_source ;

    for (j = dct->ncoef ; j <= dct->ncoef ; j++)  // don't loop for now
    {
      DCT *dct_tmp ;

      printf("************ minimizing with %d coefficients *************\n", j) ;
      dct_tmp = DCTalloc(j, mri_source) ;
      DCTcopy(dct, dct_tmp) ;
      powell_minimize(vl_target, vl_source, dct_tmp, mri_orig_source) ;
      DCTcopy(dct_tmp, dct) ;
      DCTfree(&dct_tmp) ;
    }

    if (apply_transform) 
    {
      MRI *mri_aligned ;
      char   fname[STRLEN] ;
      double overlap ;

      FileNameRemoveExtension(out_fname, fname) ;
      strcat(fname, ".mgz") ;
      overlap = (*pf_overlap)(vl_target, vl_source, dct) ;
      mri_aligned = DCTapply(dct, mri_orig_source, NULL, NULL, SAMPLE_NEAREST) ;
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
    VLSTfree(&vl_target) ;
    VLSTfree(&vl_source) ;
  }
  {
    FILE *fp = fopen(fname, "w") ;
    printf("writing DCT to %s...\n", fname) ;
    DCTdump(dct, fp) ;
    fclose(fp) ;
  }
	if (apply_transform)
	{
		MRI *mri_aligned ;
		char   fname[STRLEN] ;
			
		FileNameRemoveExtension(out_fname, fname) ;
		strcat(fname, ".mgz") ;
    mri_aligned = DCTapply(dct, mri_orig_source, NULL, NULL, SAMPLE_NEAREST) ;
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
  else if (!stricmp(option, "distance")) 
  {
    pf_overlap = compute_distance_transform_sse ;
    printf("using distance transform for SSE\n") ;
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
    mp.navgs = 1024 ;
    mp.sigma = 0 ;
    mp.l_distance = 0 ;
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
	case 'P':
		PADVOX = atoi(argv[2]) ;
		nargs = 1 ;
		printf("padding gcam with %d voxels\n", PADVOX) ;
		break ;
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
	case 'L':
		which = LABEL ;
		target_aseg_label = atoi(argv[2]) ;
		printf("using %s %d as target label from source and destination\n", cma_label_to_name(target_aseg_label),target_aseg_label) ;
		nargs = 1 ;
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




static int
check_angio_labels(MRI *mri_source, MRI *mri_target)
{
	int   label, imin, imax, scount, tcount ;
	float smin, smax, tmin, tmax ;

	MRInonzeroValRange(mri_source, &smin, &smax) ;
	MRInonzeroValRange(mri_target, &tmin, &tmax) ;
	imin = nint(MIN(smin,tmin)) ; imax = nint(MAX(smax,tmax)) ;

	printf("checking labels %d (%s) --> %d (%s)\n", imin, 
				 cma_label_to_name(imin), imax, cma_label_to_name(imax)) ;
	for (label = imin ; label <= imax ; label++)
	{
		scount = MRIvoxelsInLabel(mri_source, label) ;
		tcount = MRIvoxelsInLabel(mri_target, label) ;
		if (scount == 0 && tcount > 0)
		{
			printf("removing label %s (%d) from target volume\n", cma_label_to_name(label), label) ;
			MRIreplaceValues(mri_target, mri_target, label, 0) ;
		}
		else if (tcount == 0 && scount > 0)
		{
			printf("removing label %s (%d) from source volume\n", cma_label_to_name(label), label) ;
			MRIreplaceValues(mri_source, mri_source, label, 0) ;
		}
	}
	return(NO_ERROR) ;
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
#define MAX_GLOBAL_COEFS 5 // 15 total - one in each dir

static DCT *
find_optimal_dct(DCT *dct, MRI *mri_source, MRI *mri_target, VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, int ncoef, int skip)
{
  DCT   *dct_tmp ;
  int   i, global_coefs=0 ;
  double sse ;

  if (dct == NULL)
  {
    dct = DCTalloc(ncoef,Gmri_source = mri_source) ;
    MatrixWrite(dct->m_x_basis, "xdct.mat", "xdct") ;
    MatrixWrite(dct->m_x_basis, "ydct.mat", "ydct") ;
    MatrixWrite(dct->m_x_basis, "zdct.mat", "zdct") ;
    write_snapshot(dct, mri_source, mri_target, "test.mgz") ;
    global_coefs = ncoef > MAX_GLOBAL_COEFS ? MAX_GLOBAL_COEFS : ncoef ;
    global_coefs = 0 ;
    for (i = 1 ; i <= global_coefs ; i++)
    {
      printf("***** searching DCT space with %d coefficients *******\n",i);
      dct_tmp = DCTalloc(i, mri_source) ;
      sse = compute_global_dct_optimum(dct_tmp, mri_source, mri_target, vl_source, vl_target) ;
      DCTcopy(dct_tmp, dct) ;
      DCTfree(&dct_tmp) ;
    }
  }
  else if (global_coefs >= dct->ncoef)
    sse = compute_global_dct_optimum(dct, mri_source, mri_target, vl_source, vl_target) ;

  return(dct) ;
}

static int
write_snapshot(DCT *dct, MRI *mri_source, MRI *mri_target, const char *fname)
{
  MRI *mri_aligned ;

  mri_aligned = DCTapply(dct, mri_source, NULL, NULL, SAMPLE_NEAREST) ;
  printf("writing snapshot to %s\n", fname) ;
  MRIwrite(mri_aligned, fname) ;
  MRIfree(&mri_aligned) ;
  return(NO_ERROR) ;
}

static double
compute_distance_transform_sse(VOXEL_LIST *vl_target,VOXEL_LIST *vl_source,DCT *dct)
{
  int     width, height, depth,
  hwidth, hheight, hdepth, i ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_target, *mri_source ;
  double  sse, error ;
  double  d1, d2, xd, yd, zd ;
  MATRIX  *m_L_inv, *m_L ;
  float   xf, yf, zf ;

  mri_target = vl_target->mri2 ;
  mri_source = vl_source->mri2 ;
  DCTtransformVoxlist(dct, vl_source) ;
  m_L = MRIgetVoxelToVoxelXform(mri_source, mri_target) ;
  m_L_inv = MatrixInverse(m_L, NULL) ;
  if (m_L_inv == NULL)
    ErrorExit
    (ERROR_BADPARM, "compute_distance_transform_sse: singular matrix.") ;


  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  *MATRIX_RELT(v2, 4, 1) = 1.0 ;

  width = mri_target->width ; height = mri_target->height; depth = mri_target->depth;
  hwidth = mri_source->width ; hheight = mri_source->height ; hdepth = mri_source->depth;

  /* go through both voxel lists and compute the sse
     map it to the source, and if the source hasn't been counted yet, count it.
  */

  sse = 0.0 ;
  for (i = 0 ; i < vl_source->nvox ; i++) 
  {
    xf = vl_source->xd[i] ; yf = vl_source->yd[i] ; zf = vl_source->zd[i] ;

    V3_X(v1) = xf ; V3_Y(v1) = yf ; V3_Z(v1) = zf ;
    MatrixMultiply(m_L, v1, v2) ;
    d1 = MRIgetVoxVal(vl_source->mri2, nint(xf), nint(yf), nint(zf), 0) ;
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
    error = d1-d2 ;
    sse += error*error ;
  }

#if 0
  /* now count target voxels that weren't mapped to in union */
  for (i = 0 ; i < vl_target->nvox ; i++) 
  {
    x = vl_target->xi[i] ;
    y = vl_target->yi[i] ;
    z = vl_target->zi[i] ;
    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    MatrixMultiply(m_L_inv, v1, v2) ;
    d1 = MRIgetVoxVal(vl_target->mri2, x, y, z, 0) ;

    xd = V3_X(v2) ;
    yd = V3_Y(v2) ;
    zd = V3_Z(v2) ;
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
    error = d1-d2 ;
    sse += error*error ;
  }
#endif

  VectorFree(&v1) ;
  VectorFree(&v2) ;
  MatrixFree(&m_L_inv) ;
  return(-sqrt(sse / (double)(vl_target->nvox + vl_source->nvox))) ;
}

static double
compute_overlap(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, DCT *dct) 
{
  int     x=0, y=0, z=0, i ;
  MRI     *mri_target, *mri_source ;
  double  val, source_vox, target_vox, xr, yr, zr, overlap ;
  float   xf, yf, zf ;
  static VECTOR  *v1 = NULL, *v2 ;
  MATRIX  *m_L ;

  DCTupdate(dct) ;
  mri_target = vl_target->mri ; mri_source = vl_source->mri ;
  m_L = MRIgetVoxelToVoxelXform(mri_source, mri_target) ;

  if (v1 == NULL) {
    v1 = VectorAlloc(4, MATRIX_REAL) ; v2 = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v1, 4, 1) = 1.0 ; *MATRIX_RELT(v2, 4, 1) = 1.0 ;
  }

  DCTtransformVoxlist(dct, vl_source) ;
  DCTinverseTransformVoxlist(dct, vl_target) ;

  /* first go through target volume and for every voxel that is on in it,
     map it to the source, and if the source is on, add one to the overlap
  */

  /* go  through source volume and for every voxel that is on in it,
     map it to the target, and if the target hasn't been counted yet, count it.
  */
  for (source_vox = 0.0, i = 0 ; i < vl_source->nvox ; i++) 
  {
    xf = vl_source->xd[i] ; yf = vl_source->yd[i] ; zf = vl_source->zd[i] ;

    V3_X(v1) = xf ; V3_Y(v1) = yf ; V3_Z(v1) = zf ;
    if (nint(x) == Gx && nint(y) == Gy && nint(z) == Gz)
      DiagBreak() ;
    MatrixMultiply(m_L, v1, v2) ;
    xr = V3_X(v2) ; yr = V3_Y(v2) ;zr = V3_Z(v2) ;
    MRIsampleVolume(mri_target, xr, yr, zr, &val) ;
    source_vox += (val/binary_label) ; // # of src vox mapped into target
  }

  // the DCT inversion already maps back to source space (to apply the DCT)
  for (target_vox = 0.0, i = 0 ; i < vl_target->nvox ; i++) 
  {
    xf = vl_target->xd[i] ; yf = vl_target->yd[i] ; zf = vl_target->zd[i] ;
    MRIsampleVolume(mri_source, xf, yf, zf, &val) ;
    target_vox += (val/binary_label) ; // # of target vox mapped into source
  }

  MatrixFree(&m_L) ;
  overlap = 0.5 * (source_vox / vl_source->nvox + target_vox / vl_target->nvox);
  return(overlap) ;
}

static double scales[] = { 100, 10.0, 5, 0.75, 0.5, 0.25} ;
#define NSCALES (sizeof(scales) / sizeof(scales[0]))

static double
compute_global_dct_optimum(DCT *dct, MRI *mri_source, 
                           MRI *mri_target, 
                           VOXEL_LIST *vl_source, 
                           VOXEL_LIST *vl_target)
{
  double overlap, scale ;
  int    i ;

  for (i = 0 ; i < NSCALES ; i++)
  {
    scale = scales[i] ;
    printf("***** searching DCT %dth scale %2.3f *******\n",i, scale);
    overlap = (*pf_overlap)(vl_target, vl_source, dct) ;
    overlap = 
      compute_recursive_optimum(dct, scale, vl_source, vl_target, 0, 0,
                                overlap) ;
  }

  return(overlap) ;
}
static double
compute_recursive_optimum(DCT *dct, double scale, 
                          VOXEL_LIST *vl_source, 
                          VOXEL_LIST *vl_target, 
                          int which_coord, int k,
                          double best_overlap)
{
  double overlap, best_coef, orig_coef=0.0,coef;
  float  *pk ;

  switch (which_coord)
  {
  case 0:  // x
    pk = &VECTOR_ELT(dct->v_xk, k+1) ;
    break; 
  case 1:  //y 
    pk = &VECTOR_ELT(dct->v_yk, k+1) ;
    break ;
  case 2:  // z
    pk = &VECTOR_ELT(dct->v_zk, k+1) ;
    break ;
  default:
    return(0) ;
    break;
  }

  best_coef = *pk ;
  
  for (coef = orig_coef - scale ; 
       coef <= orig_coef + scale*1.1;
       coef += scale )
  {
    *pk = coef ;
    if (k < dct->ncoef-1)
      overlap = compute_recursive_optimum(dct, scale, vl_source, vl_target,
                                which_coord, k+1, best_overlap);
    else if (which_coord < 2)
      overlap = compute_recursive_optimum(dct, scale, vl_source, vl_target,
                                          which_coord+1, 0, best_overlap);
    else
      overlap = (*pf_overlap)(vl_target, vl_source, dct) ;
    if (overlap > best_overlap)
    {
      best_overlap = overlap ;
      best_coef = coef ;
      printf("new optimum found at (%d, %d): %2.9f\n",
             which_coord, k, best_overlap) ;
    }
  }
  *pk = best_coef ;

  //  printf("scale %2.2f, changing %s: %d\n", scale, which_coord == 0 ? "x" : which_coord == 1 ? "y" : "z", k);
  return(best_overlap) ;
}

static int
powell_minimize(VOXEL_LIST *vl_target,
                VOXEL_LIST *vl_source,
                DCT *dct,
                MRI *mri_orig_source) 
{
  float *p, **xi, fret, fstart;
  int   i, k, iter, nparms, r, c, start_t = 0 ;

  p = vector(1, 3*dct->ncoef) ;
  nparms = 3*dct->ncoef ;
  xi = matrix(1, nparms, 1, nparms) ;
  for (i = k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_xk, k) ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_yk, k) ;
  for (k = 1 ; k <= dct->ncoef ; k++)
    p[i++] = VECTOR_ELT(dct->v_zk, k) ;

  Gncoef = dct->ncoef ;
  Gvl_target = vl_target ;
  Gvl_source = vl_source ;
  for (r = 1 ; r <= nparms ; r++)
    for (c = 1 ; c <= nparms ; c++)
      xi[r][c] = r == c ? 1 : 0 ;

  // TODO:  powell(p, xi, nparms, TOL, &iter, &fret, compute_powell_sse);
  fstart = compute_powell_sse(p) ;
  printf("starting sse = %2.3f\n", fstart) ;
  OpenPowell(p, xi, nparms, TOL, &iter, &fret, compute_powell_sse);
  printf("ending sse = %2.3f\n", fret) ;
  start_t += iter ;
  do {
    for (r = 1 ; r <= nparms ; r++)
      for (c = 1 ; c <= nparms ; c++)
        xi[r][c] = r == c ? 1 : 0 ;

    fstart = fret ;
    // TODO:    powell(p, xi, nparms, TOL, &iter, &fret, compute_powell_sse);
    OpenPowell(p, xi, nparms, TOL, &iter, &fret, compute_powell_sse);
    start_t += iter ;
    for (i = k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_xk, k) = p[i++] ;
    for (k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_yk, k) = p[i++] ;
    for (k = 1 ; k <= dct->ncoef ; k++)
      VECTOR_ELT(dct->v_zk, k) = p[i++] ;

    printf("%3.3d: best alignment after powell: %2.3f (%d steps)\n",
           start_t,fret, iter) ;
#if 0
    write_snapshot(vl_target->mri, mri_orig_source,
                   mat, &parms, parms.start_t++,conform,NULL);
#endif
  } while (fret < fstart) ;

  free_matrix(xi, 1, nparms, 1, nparms) ;
  free_vector(p, 1, nparms) ;
  return(NO_ERROR) ;
}

static float
compute_powell_sse(float *p) 
{
  static DCT *dct = NULL ;
  double  error ;
  int    i, k ;

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

  error = -(*pf_overlap)(Gvl_target, Gvl_source, dct) ;
  return(error) ;
}

