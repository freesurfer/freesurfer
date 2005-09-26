//
// mri_hires_register.c
//
// written by Bruce Fischl
// Nov. 9th ,2000
// 
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2005/09/26 17:34:25 $
// Revision       : $Revision: 1.1 $
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
#include "fastmarching.h"
#include "voxlist.h"

#define NONMAX 0
#define PAD      10
#define PADVOX   1

static int find_gcam_node(GCA_MORPH *gcam, int label, 
													float x_vox, float y_vox, float z_vox) ;

static int find_label = -1 ;
static float x_vox = 0.0 ;
static float y_vox = 0.0 ;
static float z_vox = 0.0 ;

static int apply_transform = 1 ;

//static MRI *estimate_densities(GCA_MORPH *gcam, MRI *mri_target, MRI *mri_intensity) ;
static int write_snapshot(MRI *mri_target, MRI *mri_source, 
													MATRIX *m_vox_xform, GCA_MORPH_PARMS *parms, 
													int fno, int conform, char *fname) ;

static int regrid = 0 ;

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;


static char *source_intensity_fname = NULL ;
char *Progname ;
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
	lateral_ventricle,	        
	entorhinal_cortex,         
	Amygdala,                  
	Cerebral_White_Matter,     
	Cerebral_Cortex,           
	Inf_Lat_Vent
} ;
#define NUM_NON_HIPPO_LABELS  (sizeof(non_hippo_labels) / sizeof(non_hippo_labels[0]))
static int target_hippo_label = Right_Hippocampus ;

static TRANSFORM  *transform = NULL ;
static GCA_MORPH_PARMS mp ;

#define NONE  0
#define ANGIO 1
#define HIPPO 2

static int which = ANGIO ;

int
main(int argc, char *argv[])
{
	char         **av, *source_fname, *target_fname, *out_fname, fname[STRLEN] ;
  int          ac, nargs, i, new_transform = 0, pad ;
	MRI          *mri_target, *mri_source, *mri_tmp, *mri_orig_source ;
#if NONMAX
	MRI          *mri_dist_target = NULL, *mri_dist_source_sup, 
		           *mri_dist_target_sup, *mri_dist_source = NULL ;
#endif
	MRI_REGION   box ;
  struct timeb start ;
  int          msec, minutes, seconds, label ;
	GCA_MORPH    *gcam ;
	MATRIX       *m_L, *m_I ;
	LTA          *lta ;


	/* for nonlinear morph */
	mp.l_jacobian = 1 ;
	mp.l_distance = 1 ;
	mp.l_binary = .025 ;
	mp.dt = 0.005 ;
	mp.noneg = True ;
	mp.exp_k = 20 ;
	mp.momentum = 0.9 ;
	if (FZERO(mp.l_smoothness))
		mp.l_smoothness = 1 ;
	mp.sigma = 8 ;
	mp.relabel_avgs = -1 ;
	mp.navgs = 256 ;
	mp.levels = 6 ;
	mp.integration_type = GCAM_INTEGRATE_BOTH ;
	mp.nsmall = 1 ;
	mp.reset_avgs = -1 ;
	mp.npasses = 3 ;
	mp.regrid = regrid? True : False ;
	mp.tol = 0.1 ;
	mp.niterations = 1000 ;
	
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
  strcpy(mp.base_name, fname) ;
	mri_source = MRIread(source_fname) ;
	if (!mri_source)
		ErrorExit(ERROR_NOFILE, "%s: could not read source label volume %s",
							Progname, source_fname) ;

	mri_target = MRIread(target_fname) ;
	if (!mri_target)
		ErrorExit(ERROR_NOFILE, "%s: could not read target label volume %s",
							Progname, target_fname) ;
	MRIboundingBox(mri_source, 0, &box) ;
	pad = (int)ceil(PADVOX * 
									MAX(mri_target->xsize,MAX(mri_target->ysize,mri_target->zsize)) / 
									MIN(mri_source->xsize,MIN(mri_source->ysize,mri_source->zsize))); 
	mri_tmp = MRIextractRegionAndPad(mri_source, NULL, &box, pad) ;
	printf("padding source with %d voxels...\n", pad) ;
	if (pad < 1)
		pad = 1 ;
	if (Gdiag & DIAG_WRITE && DIAG_VERBOSE)
		MRIwrite(mri_tmp, "t.mgz") ;
	MRIfree(&mri_source) ;
	mri_source = mri_tmp ;
	mri_orig_source = MRIcopy(mri_source, NULL) ;

	if (which == ANGIO)
	{
		for (i = 0 ; i < NUM_NON_ARTERY_LABELS ; i++)
		{
			label = non_artery_labels[i] ;
			MRIreplaceValues(mri_source, mri_source, label, 0) ;
			MRIreplaceValues(mri_target, mri_target, label, 0) ;
		}
	}
	else if (which == HIPPO)
	{
		MRI *mri_tmp ;

		for (i = 0 ; i < NUM_NON_HIPPO_LABELS ; i++)
		{
			label = non_hippo_labels[i] ;
			MRIreplaceValues(mri_source, mri_source, label, 0) ;
		}
		/* remove other labels from segmentation image */
		mri_tmp = MRIclone(mri_target, NULL) ;
		MRIcopyLabel(mri_target, mri_tmp, target_hippo_label) ;
		MRIfree(&mri_target) ; mri_target = mri_tmp ;
	}

	MRIbinarize(mri_target, mri_target, 1, 0, target_label) ;
	MRIbinarize(mri_source, mri_source, 1, 0, target_label) ;
#if 1
	if (mri_target->type != MRI_UCHAR)
	{
		mri_tmp = MRIchangeType(mri_target, MRI_UCHAR, 0, 255, 1) ;
		MRIfree(&mri_target) ;
		mri_target = mri_tmp ;
	}
	if (mri_source->type != MRI_UCHAR)
	{
		mri_tmp = MRIchangeType(mri_source, MRI_UCHAR, 0, 255, 1) ;
		MRIfree(&mri_source) ;
		mri_source = mri_tmp ;
	}
#endif
	MRIwrite(mri_target, "target_labels.mgz") ;
	MRIwrite(mri_source, "src_labels.mgz") ;

#if NONMAX
	printf("creating distance transforms...\n") ;
	
	mri_dist_source = MRIdistanceTransform(mri_source, NULL, target_label, -1, DTRANS_MODE_SIGNED);
	mri_dist_target = MRIdistanceTransform(mri_target, NULL, target_label, -1, DTRANS_MODE_SIGNED);
	MRIscalarMul(mri_dist_source, mri_dist_source, -1) ;
	MRIscalarMul(mri_dist_target, mri_dist_target, -1) ;
	MRIwrite(mri_dist_source, "dist_source.mgz") ;
	MRIwrite(mri_dist_target, "dist_target.mgz") ;
	mri_dist_source_sup = MRInonMaxSuppress(mri_dist_source, NULL, 0, 1) ;
	mri_dist_target_sup = MRInonMaxSuppress(mri_dist_target, NULL, 0, 1) ;
	MRIwrite(mri_dist_source_sup, "dist_source_sup.mgz") ;
#else
#endif
	
  if (Gdiag & DIAG_WRITE && mp.write_iterations > 0)
  {
		sprintf(fname, "%s_target", mp.base_name) ;
    MRIwriteImageViews(mri_target, fname, IMAGE_SIZE) ;	
		MRIwrite(mri_target, "target.mgz") ;
	}

	{
		mp.max_grad = 0.3*mri_source->xsize ;
		
		if (transform->type != MORPH_3D_TYPE)  // initializing m3d from a linear transform
		{
			new_transform = 1 ;
			lta = ((LTA *)(transform->xform)) ;
			m_L = MRIrasXformToVoxelXform(mri_source, mri_target, lta->xforms[0].m_L, NULL) ;
			MatrixFree(&lta->xforms[0].m_L) ;
			if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
				write_snapshot(mri_target, mri_source, m_L, &mp, 0, 1,"linear_init");

			// transform RAS xform of both source and non-max suppressed source to
			// incorporate xform
			printf("voxel xform:\n") ;
			MatrixPrint(stdout, m_L) ;
#if 0
			mri_tmp = MRITransformedCenteredMatrix(mri_dist_source_sup, mri_target, m_L) ;
			MRIfree(&mri_dist_source_sup) ;
			mri_dist_source_sup = mri_tmp ;
#endif
			mri_tmp = MRITransformedCenteredMatrix(mri_source, mri_target, m_L) ;
			MRIfree(&mri_source) ;
			mri_source = mri_tmp ;

			if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
			{
				MRIwrite(mri_source, "s1.mgz") ;
#if NONMAX
				MRIwrite(mri_dist_source_sup, "s1max.mgz") ;
#endif
			}
			m_I = MatrixIdentity(4, NULL) ;
			MRIrasXformToVoxelXform(mri_source, mri_target, m_I, m_L);
			MatrixFree(&m_I) ;

			/* make sure none of the labels are on the border */
			MRIboundingBox(mri_source, 0, &box) ;
			if (Gdiag & DIAG_WRITE && DIAG_VERBOSE)
				mri_tmp = MRIextractRegionAndPad(mri_source, NULL, &box, pad) ;
			MRIfree(&mri_source) ; mri_source = mri_tmp ;
			if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
				MRIwrite(mri_source, "source_xformed.mgz") ;

			lta->xforms[0].m_L = m_L ;
			printf("initializing GCAM with vox->vox matrix:\n") ;
			MatrixPrint(stdout, m_L) ;
#if NONMAX
			gcam = GCAMcreateFromIntensityImage(mri_dist_source_sup, mri_target, transform) ;
#else
			gcam = GCAMcreateFromIntensityImage(mri_source, mri_target, transform) ;
#endif
			GCAMinitLabels(gcam, mri_orig_source) ;
		}
		else  /* use a previously create morph and integrate it some more */
		{
			printf("using previously create gcam...\n") ;
			gcam = (GCA_MORPH *)(transform->xform) ;
			GCAMrasToVox(gcam, mri_target) ;
		}
		if (gcam->width != mri_source->width ||
				gcam->height != mri_source->height ||
				gcam->depth != mri_source->depth)
			ErrorExit(ERROR_BADPARM, "%s: warning gcam (%d, %d, %d), doesn't match source vol (%d, %d, %d)",
								Progname, gcam->width, gcam->height, gcam->depth,
								mri_source->width, mri_source->height, mri_source->depth) ;

		mp.mri_binary = MRIcopy(mri_target, NULL) ;
		switch (which)
		{
		case ANGIO:
			mp.mri_diag = mri_target ;
			mp.diag_morph_from_atlas = 1 ;
			mp.diag_mode_filter = 1 ;
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
			mp.mri_diag = mri_target ;
			mp.diag_morph_from_atlas = 1 ;
			mp.diag_mode_filter = 1 ;
			mp.diag_volume = GCAM_LABEL ;
			break ;
		}

    if (mp.write_iterations != 0)
    {
      char fname[STRLEN] ;
      MRI  *mri_gca ;

			sprintf(fname, "%s_target.mgz", mp.base_name) ;
			if (mp.diag_morph_from_atlas)
			{
				printf("writing target volume to %s...\n", fname) ;
				MRIwrite(mri_target, fname) ;
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
		if (find_label >= 0)
			find_gcam_node(gcam, find_label, x_vox, y_vox, z_vox) ;
#if NONMAX
		GCAMregister(gcam, mri_dist_target, &mp) ; // atlas is source, morph target into register with it
#else
		GCAMregister(gcam, mri_target, &mp) ; // atlas is source, morph target into register with it
#endif
		if (apply_transform)
		{
			MRI *mri_aligned ;
			char   fname[STRLEN] ;
			
			FileNameRemoveExtension(out_fname, fname) ;
			strcat(fname, ".mgz") ;
			mri_aligned = GCAMmorphFieldFromAtlas(gcam, mp.mri, GCAM_LABEL,0, 1) ;
			printf("writing transformed output volume to %s...\n", fname) ;
			MRIwrite(mri_aligned, fname) ;
			MRIfree(&mri_aligned) ;
		}
		GCAMvoxToRas(gcam) ;
		GCAMwrite(gcam, out_fname) ;
		GCAMrasToVox(gcam, mri_target) ;
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
  else if (!stricmp(option, "OPTIMAL"))
  {
    mp.integration_type = GCAM_INTEGRATE_OPTIMAL ;
    printf("using optimal time-step integration\n") ;
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
  else if (!stricmp(option, "intensity") ||!stricmp(option, "ll"))
  {
    mp.l_log_likelihood = atof(argv[2]) ;
		nargs = 1 ;
    printf("setting l_log_likelihood = %2.1f\n", mp.l_log_likelihood );
  }
  else if (!stricmp(option, "noregrid"))
  {
		regrid = 0 ;
		mp.regrid = False ;
    printf("disabling regridding...\n") ;
  }
  else if (!stricmp(option, "regrid"))
  {
		regrid = 1 ;
		mp.regrid = True ;
    printf("enabling regridding...\n") ;
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
  else if (!stricmp(option, "area_intensity"))
  {
    mp.l_area_intensity = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_area_intensity=%2.3f\n", mp.l_area_intensity) ;
  }
  else if (!stricmp(option, "tol"))
  {
    mp.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("using tol=%2.3f\n", mp.tol) ;
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
		printf("assuming source is hires hippo and dst is aseg volume\n") ;
	}
	else if (!stricmp(option, "none"))
	{
		which = NONE ;
		printf("making no assumptions about labels (not angio or hippo\n") ;
	}
	else switch (*option)
	{
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
    mp.exp_k = atof(argv[2]) ;
    printf("setting exp_k to %2.2f (default=%2.2f)\n",
           mp.exp_k, EXP_K) ;
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
write_snapshot(MRI *mri_target, MRI *mri_source, MATRIX *m_vox_xform, 
							 GCA_MORPH_PARMS *parms, int fno, int conform, char *in_fname)
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

#if 0
	if (mri_source_intensity)
	{
#if 0
		mri_aligned = MRIsrcTransformedCentered(mri_source_intensity, mri_target, m_vox_xform, SAMPLE_TRILINEAR) ;
#else
		mri_aligned = MRITransformedCenteredMatrix(mri_source_intensity, mri_target, m_vox_xform) ;
#endif
		if (in_fname)
			sprintf(fname, "intensity_%s_%s", parms->base_name, in_fname) ;
		else
			sprintf(fname, "intensity_%s_%03d", parms->base_name, fno) ;
		MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
		if (in_fname)
			sprintf(fname, "intensity_%s_%s.mgz", parms->base_name, in_fname) ;
		else
			sprintf(fname, "intensity_%s_%03d.mgz", parms->base_name, fno) ;
		printf("writing snapshot to %s...\n", fname) ;
		MRIwrite(mri_aligned, fname) ;
		MRIfree(&mri_aligned) ;
	}
#endif

	return(NO_ERROR) ;
}

#if 0
static MRI *
estimate_densities(GCA_MORPH *gcam, MRI *mri_target, MRI *mri_intensities)
{
	GCAM_LABEL_TRANSLATION_TABLE gcam_ltt ;

	gcam_ltt.nlabels = 0 ;

	memset(gcam_ltt.means, 0, sizeof(gcam_ltt.means)) ;

	/* don't use inf_lat_vent label as it may be too small to
		 give reliable estimate of density */
  gcam_ltt.input_labels[gcam_ltt.nlabels] = alveus ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 530 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = perforant_pathway ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 530 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = parasubiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = presubiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = subiculum ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA1 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA2 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA3 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
  gcam_ltt.input_labels[gcam_ltt.nlabels] = CA4 ;
  gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = GC_DG ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = HATA ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = fimbria ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 530 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = lateral_ventricle ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Lateral_Ventricle ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 325 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = molecular_layer_HP ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = hippocampal_fissure ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Lateral_Ventricle ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 325 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = entorhinal_cortex ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_Cortex ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = molecular_layer_subiculum ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Amygdala ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Hippocampus ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Cerebral_White_Matter ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_White_Matter ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 580 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Cerebral_Cortex ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Cerebral_Cortex ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 450 ;
	gcam_ltt.nlabels++ ;
	gcam_ltt.input_labels[gcam_ltt.nlabels] = Inf_Lat_Vent ;
	gcam_ltt.output_labels[gcam_ltt.nlabels] = Right_Lateral_Ventricle ;
	if (fix_intensity)
		gcam_ltt.means[gcam_ltt.nlabels] = 325 ;
	gcam_ltt.nlabels++ ;

	return(GCAMinitDensities(gcam, mri_target, mri_intensities, &gcam_ltt)) ;
}
#endif
static int
find_gcam_node(GCA_MORPH *gcam, int label, float x0, float y0, float z0)
{
	int   x, y, z, nx, ny, nz ;
	double dist, min_dist, dx, dy, dz ;
	GCA_MORPH_NODE *gcamn ;

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




