/**
 * @brief nonlinear alignment of binary images.
 *
 * Basically a binary implementation of the algorithm in:
 *
 * Fischl B, Salat DH, van der Kouwe AJW, Makris N, Ségonne F, Dale
 * AM. Sequence-Independent  Segmentation of Magnetic Resonance Images.
 * NeuroImage, 2004; 23 Suppl 1, S69-84.
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "fastmarching.h"
#include "voxlist.h"

#define NONMAX 0
#define PAD      10
#define MAX_DISTANCE 100

static int PADVOX = 1 ;

static int check_angio_labels(MRI *mri_source, MRI *mri_target) ;
static int find_gcam_node(GCA_MORPH *gcam, int label, 
			  float x_ras, float y_ras, float z_ras) ;

static int surf_flag = 0 ;
static float smooth_intensities = -1.0 ;
static int find_label = -1 ;
static float x_ras = 0.0 ;
static float y_ras = 0.0 ;
static float z_ras = 0.0 ;
static int mode_filters = 0 ;
static int aseg = 0 ;

static double min_sigma = 1.0 ;

static int upsample = 0 ;
static int apply_transform = 1 ;

int MRImapRegionToTargetMRI(MRI *mri_src, MRI *mri_dst, MRI_REGION *box) ;
#if 0
static MRI *estimate_densities(GCA_MORPH *gcam, MRI *mri_target, MRI *mri_intensity, MRI *mri_src) ;
#endif
static int write_snapshot(MRI *mri_target, MRI *mri_source, 
			  MATRIX *m_vox_xform, GCA_MORPH_PARMS *parms, 
			  int fno, int conform, const char *fname) ;

static int regrid = 0 ;

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;


static char *source_intensity_fname = NULL ;
const char *Progname ;
static int target_label = 128 ;

static int skip = 2 ;
static double distance = 1.0 ;

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
	Brain_Stem,
  left_fornix,
  right_fornix
} ;

#define NUM_NON_HIPPO_LABELS  (sizeof(non_hippo_labels) / sizeof(non_hippo_labels[0]))
static int target_aseg_label = Right_Hippocampus ;

static TRANSFORM  *transform = NULL ;
static GCA_MORPH_PARMS mp ;

#define NONE  0
#define ANGIO 1
#define HIPPO 2
#define WM    3
#define LABEL 4
#define ASEG  5

static int which = ANGIO ;

static MRI *mri_norm = NULL ;

int
main(int argc, char *argv[])
{
	char         **av, *source_fname, *target_fname, *out_fname, fname[STRLEN] ;
  int          ac, nargs, i, new_transform = 0, pad ;
	MRI          *mri_target = nullptr, *mri_source, *mri_tmp, *mri_orig_source, *mri_orig_target ;
#if NONMAX
	MRI          *mri_dist_target = nullptr, *mri_dist_source_sup, 
		*mri_dist_target_sup, *mri_dist_source = nullptr
     ;
#endif
	MRI_REGION   box ;
  Timer start ;
  int          msec, hours, minutes, seconds, label ;
	GCA_MORPH    *gcam ;
	MATRIX       *m_L/*, *m_I*/ ;
	LTA          *lta ;


	/* for nonlinear morph */
	mp.l_jacobian = 1 ;
	mp.l_distance = 1 ;
	mp.l_binary = .025 ;
	mp.dt = 0.005 ;
	mp.noneg = True ;
	mp.exp_k = 5 ;
	mp.momentum = 0.9 ;
#if 0
	if (FZERO(mp.l_smoothness))
		mp.l_smoothness = .001 ;
#endif
	mp.sigma = 8 ;
	mp.relabel_avgs = -1 ;
	mp.uncompress = 1 ;      // remove compression each time step
	mp.ratio_thresh = 0.25 ; // nodes with area/orig smaller than this are compressed
	mp.navgs = 256 ;
	mp.levels = 6 ;
	mp.integration_type = GCAM_INTEGRATE_BOTH ;
	mp.nsmall = 1 ;
	mp.reset_avgs = -1 ;
	mp.npasses = 3 ;
	mp.regrid = regrid? True : False ;
	mp.tol = 0.1 ;
	mp.niterations = 1000 ;
	mp.scale_smoothness = 1 ;
	
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
  mp.write_fname = out_fname ;
	printf("source = %s\ntarget = %s\noutput = %s\n", source_fname, target_fname,out_fname);
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(mp.base_name, fname) ;
	mri_source = MRIread(source_fname) ;
	if (!mri_source)
		ErrorExit(ERROR_NOFILE, "%s: could not read source label volume %s",Progname, source_fname) ;

  if (surf_flag)
  {
    MRI_SURFACE *mris ;
    mris = MRISread(target_fname) ;
    if (!mri_target)
      ErrorExit(ERROR_NOFILE, "%s: could not read target surface %s",
                Progname, target_fname) ;
    mri_target = MRISfillInterior(mris, .25, NULL) ;
    MRIreplaceValues(mri_target, mri_target, 1, target_aseg_label) ;
    MRISfree(&mris) ;
  }
  else
  {
    mri_target = MRIread(target_fname) ;
    if (!mri_target)
      ErrorExit(ERROR_NOFILE, "%s: could not read target label volume %s", Progname, target_fname) ;
  }
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
  if (aseg == 0)
  {
    MRIboundingBox(mri_source, 0, &box) ;
    pad = PADVOX ;
    printf("padding source with %d voxels...\n", pad) ;
    mri_tmp = MRIextractRegionAndPad(mri_source, NULL, &box, pad) ;
    MRIfree(&mri_source) ;
    mri_source = mri_tmp ;
    //	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
		MRIwrite(mri_source, "s.mgz") ;
  }
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
		MRIbinarize(mri_target, mri_target, 1, 0, target_label) ;
		MRIbinarize(mri_source, mri_source, 1, 0, target_label) ;
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
  if (aseg == 0)
  {
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
    
#if 1
    if (mri_norm)  // won't work if the norm and aseg aren't in the same voxel coords
    {
      MRI_REGION box_norm ;
      
      *(&box_norm) = *(&box) ;
      MRImapRegionToTargetMRI(mri_orig_target, mri_norm, &box_norm) ;
      mri_tmp = MRIextractRegionAndPad(mri_norm, NULL, &box_norm, pad) ;
      MRIfree(&mri_norm) ;
      mri_norm = mri_tmp ;
      //		if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
			MRIwrite(mri_norm, "n.mgz") ;
    }
#endif
    mri_tmp = MRIextractRegionAndPad(mri_target, NULL, &box, pad) ;
    MRIfree(&mri_target) ;mri_target = mri_tmp ;
    mri_tmp = MRIextractRegionAndPad(mri_orig_target, NULL, &box, pad) ;
    MRIfree(&mri_orig_target) ;
    mri_orig_target = mri_tmp ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_target, "t.mgz") ;
    
#if NONMAX
    printf("creating distance transforms...\n") ;
    
    mri_dist_source = MRIdistanceTransform(mri_source, NULL, target_label, -1, DTRANS_MODE_SIGNED, NULL);
    mri_dist_target = MRIdistanceTransform(mri_target, NULL, target_label, -1, DTRANS_MODE_SIGNED, NULL);
    MRIscalarMul(mri_dist_source, mri_dist_source, -1) ;
    MRIscalarMul(mri_dist_target, mri_dist_target, -1) ;
    MRIwrite(mri_dist_source, "dist_source.mgz") ;
    MRIwrite(mri_dist_target, "dist_target.mgz") ;
    mri_dist_source_sup = MRInonMaxSuppress(mri_dist_source, NULL, 0, 1) ;
    mri_dist_target_sup = MRInonMaxSuppress(mri_dist_target, NULL, 0, 1) ;
    MRIwrite(mri_dist_source_sup, "dist_source_sup.mgz") ;
#else
#endif
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
		MRIbinarize(mri_target, mri_target, 1, 0, target_label) ;
		MRIbinarize(mri_source, mri_source, 1, 0, target_label) ;
		mp.target_label = target_label ;
	}
  if (Gdiag & DIAG_WRITE && mp.write_iterations > 0)
  {
    int req = snprintf(fname, STRLEN, "%s_target", mp.base_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwriteImageViews(mri_orig_target, fname, IMAGE_SIZE) ;	
    req = snprintf(fname, STRLEN, "%s_target.mgz", mp.base_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    MRIwrite(mri_orig_target, fname) ;
  }
	
	mp.max_grad = 0.3*mri_source->xsize ;

	if (transform == NULL || transform->type != MORPH_3D_TYPE)  // initializing m3d from a linear transform
	{
    double det ;
		new_transform = 1 ;

    if (transform == NULL)
    {
      transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
      m_L = MRIgetVoxelToVoxelXform(mri_source, mri_target) ;
      lta = ((LTA *)(transform->xform)) ;
    }
    else
    {
      lta = ((LTA *)(transform->xform)) ;
      det = MatrixDeterminant( lta->xforms[0].m_L) ;
      m_L = MRIrasXformToVoxelXform(mri_source, mri_target, lta->xforms[0].m_L, NULL) ;
    }

    MatrixFree(&lta->xforms[0].m_L) ;
		lta->xforms[0].m_L = m_L ;
		if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
			write_snapshot(mri_target, mri_source, m_L, &mp, 0, 1,"linear_init");

    det = MatrixDeterminant( lta->xforms[0].m_L) ;
		printf("initializing GCAM with vox->vox matrix:\n") ;
		MatrixPrint(stdout, m_L) ;
		gcam = GCAMalloc(mri_target->width, mri_target->height, mri_target->depth);
		GCAMinit(gcam, mri_source, NULL, transform, 0) ;
		GCAMinitLabels(gcam, mri_orig_target) ;
    gcam->gca = gcaAllocMax(1, 1, 1, 
                            mri_target->width, mri_target->height, 
                            mri_target->depth,
                            0, 0) ;
    GCAMinitVolGeom(gcam, mri_source, mri_target) ;
	}
	else  /* use a previously create morph and integrate it some more */
	{
		printf("using previously create gcam...\n") ;
		gcam = (GCA_MORPH *)(transform->xform) ;
	}
	if (find_label >= 0)
	{
		GCAMvoxToRas(gcam) ;
		find_gcam_node(gcam, find_label, x_ras, y_ras, z_ras) ;
	}

	GCAMrasToVox(gcam, mri_source) ;
	if (gcam->width != mri_target->width ||
			gcam->height != mri_target->height ||
			gcam->depth != mri_target->depth)
		ErrorPrintf(ERROR_BADPARM, "%s: warning gcam (%d, %d, %d), doesn't match source vol (%d, %d, %d)",
							Progname, gcam->width, gcam->height, gcam->depth,
							mri_target->width, mri_target->height, mri_target->depth) ;

	if (mp.mri_binary == NULL)
		mp.mri_binary = MRIcopy(mri_target, NULL) ;
	mp.diag_mode_filter = mode_filters ;
	switch (which)
	{
	case ANGIO:
		mp.mri_diag = mri_target ;
		mp.diag_morph_from_atlas = 1 ;
		mp.diag_volume = GCAM_LABEL ;
		for (i = 0 ; i < NUM_NON_ARTERY_LABELS ; i++)
		{
			label = non_artery_labels[i] ;
			GCAMsetLabelStatus(gcam, label, GCAM_BINARY_ZERO) ;
		}
		break ;
	case HIPPO:
		for (i = 0 ; i < NUM_NON_HIPPO_LABELS ; i++)
		{
			label = non_hippo_labels[i] ;
			GCAMsetLabelStatus(gcam, label, GCAM_BINARY_ZERO) ;
		}
		// no break!!
	default:
    //		mp.mri_diag = mri_target ;
		mp.diag_volume = GCAM_LABEL ;
		break ;
	}

	if (mp.write_iterations != 0)
	{
		char fname[STRLEN] ;
		MRI  *mri_gca ;

		int req = snprintf(fname, STRLEN, "%s_target.mgz", mp.base_name) ;
		if( req >= STRLEN ) {
		  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
		}

		if (mp.diag_morph_from_atlas == 0)
		{
			printf("writing target volume to %s...\n", fname) ;
			MRIwrite(mri_orig_target, fname) ;
		}
		else
		{
			mri_gca = MRIclone(mri_source, NULL) ;
			GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
			printf("writing target volume to %s...\n", fname) ;
			MRIwrite(mri_gca, fname) ;
			MRIfree(&mri_gca) ;
		}
	}

	GCAMsetLabelStatus(gcam, 0, GCAM_NEVER_USE_LIKELIHOOD) ;
	mp.mri = mri_target ;
	if (mp.regrid == True && new_transform == 0)
		GCAMregrid(gcam, mri_target, PAD, &mp, &mri_source) ;
	if (upsample > 0)
	{
		MRI *mri_morphed ;
		char fname[STRLEN] ;

		mri_morphed = GCAMmorphFieldFromAtlas(gcam, mp.mri_diag, mp.diag_volume, 0, 
						      mp.diag_mode_filter) ;
		int req = snprintf(fname, STRLEN, "%s_orig.mgz", mp.base_name) ;
		if( req >= STRLEN ) {
		  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
		}
		MRIwrite(mri_morphed, fname) ;
	}
		
  if (aseg)
  {
    MRI *mri_target_dist, *mri_warp ;
    double old_sse, sse, pct_change ;
    int    done ;

    mri_warp = MRIallocSequence(mri_target->width, 
                                mri_target->height, mri_target->depth, 
                                MRI_FLOAT, 3) ;
    mp.ndtrans = NDTRANS_LABELS ;
    mp.dtrans_labels = dtrans_labels ;
    replace_labels(mri_source, mri_source, combine_labels, NCOMBINE_LABELS, &mp) ;
    replace_labels(mri_target, mri_target, combine_labels, NCOMBINE_LABELS, &mp) ;
    mri_target_dist = MRIcreateDistanceTransforms(mri_target, NULL,
                                                  MAX_DISTANCE,
                                                  mp.dtrans_labels,
                                                  mp.ndtrans) ;
    GCAMwriteWarpToMRI(gcam, mri_warp) ;
    sse = MRIlabelMorphSSE(mri_source, mri_target, mri_warp) ;
    do
    {
      int start_t = mp.start_t, steps ;
      old_sse = sse ;
      printf("calling demons registration with sigma = %2.1f\n", mp.sigma) ;
      GCAMdemonsRegister(gcam, mri_source, mri_target, &mp,
                         MAX_DISTANCE, NULL, mri_target_dist) ;
      sse = mp.last_sse ;
      pct_change = 100*(old_sse - sse) / (old_sse) ;
      done = pct_change < mp.tol ;
      printf("old sse %f, new sse %f, pct_change = %2.3f%%, done = %d\n",
             old_sse, sse, pct_change, done) ;
      steps = mp.start_t - start_t ;
      if (steps <= 1)  // couldn't make any progress at this scale - reduce it
        mp.sigma *= .75 ;

      done = (mp.sigma < min_sigma);
    } while (!done) ;

    MRIfree(&mri_warp) ;
  }
  else
  {
    for (i = 0 ; i < upsample ; i++)
      GCAMupsample2(gcam) ;
#if NONMAX
    GCAMregister(gcam, mri_dist_source, &mp) ; // atlas is source, morph target into register with it
#else
    GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, ORIGINAL_POSITIONS) ;
    GCAMstoreMetricProperties(gcam) ;
    GCAMregister(gcam, mri_source, &mp) ; // atlas is target, morph source into register with it
#endif
  }
	if (apply_transform)
	{
		MRI *mri_aligned ;
		char   fname[STRLEN] ;
			
		FileNameRemoveExtension(out_fname, fname) ;
		strcat(fname, ".mgz") ;
    mri_aligned = GCAMmorphToAtlas(mp.mri, gcam, NULL, -1, mp.diag_sample_type) ;
		printf("writing transformed output volume to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}
	GCAMvoxToRas(gcam) ;
	GCAMwrite(gcam, out_fname) ;
	GCAMrasToVox(gcam, mri_source) ;
	if (find_label >= 0)
		find_gcam_node(gcam, find_label, x_ras, y_ras, z_ras) ;

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

extern int gcam_write_grad ;
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
    mp.diag_morph_from_atlas = atoi(argv[2]) ;
    printf("morphing %s atlas...\n", mp.diag_morph_from_atlas ? "from" : "to") ;
  }
  else if (!stricmp(option, "spring"))
  {
		mp.l_spring = atof(argv[2]) ;
		nargs = 1;
    printf("setting l_spring = %2.2f\n", mp.l_spring) ;
  }
  else if (!stricmp(option, "aseg"))
  {
    aseg = 1 ;
    mp.tol = .25 ;
    mp.sigma = 14 ;
    mp.l_smoothness = 0.1 ;
    which = ASEG ;
    printf("assuming input volumes are aseg labelings (setting tol=%2.1f, l_smooth=%2.2f)\n",
           mp.tol, mp.l_smoothness) ;
  }
  else if (!stricmp(option, "uncompress"))
  {
		mp.uncompress = atoi(argv[2]) ;
		nargs = 1;
    printf("setting uncompress = %d\n", mp.uncompress) ;
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
  else if (!stricmp(option, "wg"))
  {
		gcam_write_grad = 1 ;
    printf("writing out gradient diagnostics\n") ;
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
  else if (!stricmp(option, "diag"))  {
    mp.mri_diag = MRIread(argv[2]) ;
    if (mp.mri_diag == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read diag volume from %s", Progname, argv[2]) ;
    nargs = 1 ;
    printf("writing diagnostics for input volume %s\n", argv[2]) ;
  }
  else if (!stricmp(option, "distance"))
  {
    distance = atof(argv[2]) ;
		nargs = 1 ;
    printf("expanding border by %2.1f mm every outer cycle\n", distance);
  }
  else if (!stricmp(option, "intensity") ||!stricmp(option, "ll"))
  {
    mp.l_log_likelihood = atof(argv[2]) ;
		nargs = 2 ;
    printf("setting l_log_likelihood = %2.1f\n", mp.l_log_likelihood );
		printf("reading intensity image from %s...\n", argv[3]) ;
		mri_norm = MRIread(argv[3]) ;
		if (mri_norm == NULL)
			ErrorExit(ERROR_BADPARM, "%s: could not read intensity image from %s\n", Progname, argv[3]) ;
  }
  else if (!stricmp(option, "likelihood"))
  {
    mp.l_likelihood = atof(argv[2]) ;
		nargs = 2 ;
    printf("setting l_likelihood = %2.1f\n", mp.l_likelihood );
		printf("reading intensity image from %s...\n", argv[3]) ;
		mri_norm = MRIread(argv[3]) ;
		if (mri_norm == NULL)
			ErrorExit(ERROR_BADPARM, "%s: could not read intensity image from %s\n", Progname, argv[3]) ;
  }
  else if (!stricmp(option, "noregrid"))
  {
		regrid = 0 ;
		mp.regrid = False ;
    printf("disabling regridding...\n") ;
  }
  else if (!stricmp(option, "regrid"))
  {
		regrid = atoi(argv[2]) ;
		mp.regrid = True ;
    printf("enabling regridding...\n") ;
    nargs = 1 ;
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
  else if (!stricmp(option, "area_intensity") || !stricmp(option, "aint"))
  {
    mp.l_area_intensity = atof(argv[2]) ;
    nargs = 2 ;
    printf("using l_area_intensity=%2.3f\n", mp.l_area_intensity) ;
		printf("reading intensity image from %s...\n", argv[3]) ;
		mri_norm = MRIread(argv[3]) ;
		if (mri_norm == NULL)
			ErrorExit(ERROR_BADPARM, "%s: could not read intensity image from %s\n", Progname, argv[3]) ;
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
  else if (!stricmp(option, "min_sigma"))
  {
    min_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("using min sigma=%2.3f\n", min_sigma) ;
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
	else if (!stricmp(option, "surf"))
	{
    surf_flag = 1 ;
    printf("interpreting target as a surface\n") ;
  }
	else if (!stricmp(option, "hippo"))
	{
		which = HIPPO ;
    mp.l_binary = 0.5 ;
    mp.l_smoothness = 0.1 ;
    mp.dt = 0.005 ;
    mp.levels = 7 ;
    mp.navgs = 1024 ;
    mp.sigma = 0 ;
    mp.l_distance = 0 ;
		printf("assuming source is hires hippo and dst is aseg volume\n") ;
    printf("setting l_binary=%2.4f, l_smooth=%2.4f\n", mp.l_binary, mp.l_smoothness);
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
	case 'D':
		mp.l_distance = atof(argv[2]) ;
		nargs = 1 ;
		printf("using l_distance = %2.3f\n", mp.l_distance) ;
		break ;
  case 'M':
    mp.momentum = atof(argv[2]) ;
    nargs = 1 ;
    printf("momentum = %2.2f\n", mp.momentum) ;
    break ;
	case 'N':
		mp.niterations = atoi(argv[2]) ;
		nargs = 1 ;
		printf("using niterations = %d\n", mp.niterations) ;
		break ;
	case 'S':
		mp.l_smoothness = atof(argv[2]) ;
		nargs = 1 ;
		printf("using l_smoothness = %2.3f\n", mp.l_smoothness) ;
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
				 "usage: %s <source> <target> <warp>\n",
				 Progname) ;
	exit(ecode) ;
}


static int
write_snapshot(MRI *mri_target, MRI *mri_source, MATRIX *m_vox_xform, 
	       GCA_MORPH_PARMS *parms, int fno, int conform, const char *in_fname)
{
	MRI *mri_aligned ;
	char fname[STRLEN] ;

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
		mri_aligned = MRITransformedCenteredMatrix(mri_source, mri_target, m_vox_xform) ;
	}
	if (in_fname) {
	  int req = snprintf(fname, STRLEN, "%s_%s", parms->base_name, in_fname) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
	} else {
	  int req = snprintf(fname, STRLEN, "%s_%03d", parms->base_name, fno) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
	}
	MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
	if (in_fname) {
	  int req = snprintf(fname, STRLEN, "%s_%s.mgz", parms->base_name, in_fname) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
	} else {
	  int req = snprintf(fname, STRLEN, "%s_%03d.mgz", parms->base_name, fno) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
	}
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
		if (in_fname) {
		  int req = snprintf(fname, STRLEN, "orig_%s_%s.mgz", parms->base_name, in_fname) ;
		  if( req >= STRLEN ) {
		    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
		  }
		} else {
		  int req = snprintf(fname, STRLEN, "orig_%s_%03d.mgz", parms->base_name, fno) ;
		  if( req >= STRLEN ) {
		    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
		  }
		}
		printf("writing snapshot to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}

#if 0
	if (mri_source_intensity)
	{
#if 0
		mri_aligned = MRIsrcTransformedCentered(mri_source_intensity, mri_target, m_vox_xform, SAMPLE_TRILINEAR) ;
#else
		mri_aligned = MRITransformedCenteredMatrix(mri_source_intensity, mri_target, m_vox_xform) ;
#endif
		if (in_fname) {
		  sprintf(fname, "intensity_%s_%s", parms->base_name, in_fname) ;
		} else {
		  sprintf(fname, "intensity_%s_%03d", parms->base_name, fno) ;
		}
		MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
		if (in_fname) {
		  sprintf(fname, "intensity_%s_%s.mgz", parms->base_name, in_fname) ;
		} else {
		  sprintf(fname, "intensity_%s_%03d.mgz", parms->base_name, fno) ;
		}
		printf("writing snapshot to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}
#endif

	return(NO_ERROR) ;
}

#if 0
static MRI *
estimate_densities(GCA_MORPH *gcam, MRI *mri_target, MRI *mri_intensities, MRI *mri_src)
{
	GCAM_LABEL_TRANSLATION_TABLE gcam_ltt ;
	int                          i ;

  memset(&gcam_ltt, 0, sizeof(gcam_ltt)) ;
	gcam_ltt.nlabels = 0 ;

	memset(gcam_ltt.means, 0, sizeof(gcam_ltt.means)) ;
	memset(gcam_ltt.scales, 0, sizeof(gcam_ltt.scales)) ;

	/* don't use inf_lat_vent label as it may be too small to
		 give reliable estimate of density */
  gcam_ltt.input_labels[gcam_ltt.nlabels] = alveus ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
    gcam_ltt.nlabels++ ;
  }
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_alveus ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
    gcam_ltt.nlabels++ ;
  }
  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_alveus ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
    gcam_ltt.nlabels++ ;
  }

  gcam_ltt.input_labels[gcam_ltt.nlabels] = perforant_pathway ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = parasubiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = presubiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = subiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_presubiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_presubiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_subiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_subiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA1 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_CA1 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_CA1 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;


  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_CA2_3 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_CA2_3 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = right_CA4_DG ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = left_CA4_DG ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA2 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA3 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA4 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = GC_DG ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = HATA ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = left_fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Cerebral_White_Matter ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }
	gcam_ltt.input_labels[gcam_ltt.nlabels] = right_fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }


	gcam_ltt.input_labels[gcam_ltt.nlabels] = right_fornix ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = left_fornix ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Cerebral_White_Matter ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;


	gcam_ltt.input_labels[gcam_ltt.nlabels] = fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }

	gcam_ltt.input_labels[gcam_ltt.nlabels] = lateral_ventricle ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Lateral_Ventricle ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = molecular_layer_HP ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = right_hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = left_hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = entorhinal_cortex ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_Cortex ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = molecular_layer_subiculum ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Amygdala ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Cerebral_White_Matter ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Cerebral_Cortex ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_Cortex ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Inf_Lat_Vent ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

#if 0
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_CADG_head ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_subiculum ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_CADG_head ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_subiculum ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	gcam_ltt.scales[gcam_ltt.nlabels] = 1.0 /* .95 */ ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }
#endif

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_choroid_plexus ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Right_Inf_Lat_Vent ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_choroid_plexus ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Hippocampus ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
  {
    gcam_ltt.second_labels[gcam_ltt.nlabels] = Left_Inf_Lat_Vent ;
    gcam_ltt.second_label_pct[gcam_ltt.nlabels] = .25 ;
		gcam_ltt.nlabels++ ;
  }

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Left_Thalamus ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Left_Thalamus_Proper ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Right_Thalamus ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Thalamus_Proper ;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;

	gcam_ltt.input_labels[gcam_ltt.nlabels] = Conn_Tissue ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Lateral_Ventricle;
	if (MRIlabelInVolume(mri_src, gcam_ltt.input_labels[gcam_ltt.nlabels]))
		gcam_ltt.nlabels++ ;


	for (i = 0 ; i < gcam_ltt.nlabels ; i++)
		if (FZERO(gcam_ltt.scales[i]))
			gcam_ltt.scales[i] = 1 ;
	return(GCAMinitDensities(gcam, mri_target, mri_intensities, &gcam_ltt)) ;
}
#endif

static int
find_gcam_node(GCA_MORPH *gcam, int label, float x0, float y0, float z0)
{
	int   x, y, z, nx, ny, nz ;
	double dist, min_dist, dx, dy, dz ;
	GCA_MORPH_NODE *gcamn ;

	printf("looking for label %s at (%2.1f, %2.1f, %2.1f)\n", cma_label_to_name(label),x0,y0,z0) ;
	min_dist = 1e10 ;

	nx = ny = nz = -1 ;
	for (x = 0 ; x < gcam->width ; x++)
	{
		for (y = 0 ; y < gcam->height ; y++)
		{
			for (z = 0 ; z < gcam->depth ; z++)
			{
				gcamn = &gcam->nodes[x][y][z] ;
				if (gcamn->label != label)
					continue ;
				dx = gcamn->x - x0 ;
				dy = gcamn->y - y0 ;
				dz = gcamn->z - z0 ;
				dist = sqrt(dx*dx + dy*dy + dz*dz);
				if (dist < min_dist)
				{
					min_dist = dist ;
					nx = x ; ny = y ; nz = z ; 
				}
			}
		}
	}

	if (nx >= 0)
	{
		printf("node found at (%d, %d, %d), setting G[xyz]\n", nx, ny, nz);
		Gx = nx ; Gy = ny ; Gz = nz ;
	}
	return(NO_ERROR) ;
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
