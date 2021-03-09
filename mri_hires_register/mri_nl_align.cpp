/**
 * @brief nonlinear volumetric alignment
 *
 * see (Fischl et al., Neuroimage, 2004)
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

#include "romp_support.h"


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
#include "mrisurf.h"

#define NONMAX 0
#define PAD      10
#define PADVOX   1


static int handle_expanded_ventricles = 0 ;
static char *label_ignore_name = NULL ;
static char *label_dist_name = NULL ;
static int apply_transform = 1 ;
static int erosions = 0;
static float scale_values = 1.0 ;

static int write_snapshot(MRI *mri_target, MRI *mri_source,
                          MATRIX *m_vox_xform, GCA_MORPH_PARMS *parms,
                          int fno, int conform, const char *fname) ;

static int regrid = 0 ;
static int nowmsa = 1 ; // remove all wmsa labels from the atlas

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;

static const char *source_surf = "";
static const char *target_surf = ".toM02100023.resample";
static char *mask_fname = NULL ;

const char *Progname ;

static int skip = 2 ;
static double distance = 1.0 ;
static int use_aseg = 0 ;
static int nozero = 1 ;
static int match_peak_intensity_ratio = 0 ;
static int match_mean_intensity = 1 ;

static char *ribbon_name = NULL ;

static TRANSFORM  *transform = NULL ;
static GCA_MORPH_PARMS mp ;

static int renormalize = 1 ;
static int rip = 0 ;
static MRI *mri_target_diag = NULL ;

static MRI *
replace_wmsa(MRI *mri_src, MRI *mri_dst) 
{
  int   x, y, z, label, lh, rh, nreplaced ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;
  for (nreplaced = x = 0 ; x < mri_dst->width; x++)
    for (y = 0 ; y < mri_dst->height; y++)
      for (z = 0 ; z < mri_dst->depth; z++)
      {
	label = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
	if (IS_WMSA(label))
	{
	  lh = MRIlabelsInNbhd(mri_dst,  x,  y,  z, 4, Left_Cerebral_White_Matter) ;
	  rh = MRIlabelsInNbhd(mri_dst,  x,  y,  z, 4, Right_Cerebral_White_Matter) ;
	  label = lh > rh ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
	  MRIsetVoxVal(mri_dst, x, y, z, 0, label) ;
	  nreplaced++ ;
	}
      }
  printf("%d WMSA labels replaced with WM\n", nreplaced) ;
  return(mri_dst) ;
}


int
main(int argc, char *argv[])
{
  char         **av, *source_fname, *target_fname, *out_fname, fname[STRLEN] ;
  int          ac, nargs, new_transform = 0, pad ;
  MRI          *mri_target, *mri_source, *mri_orig_source ;
  MRI_REGION   box ;
  Timer start ;
  int          msec, minutes, seconds ;
  GCA_MORPH    *gcam ;
  MATRIX       *m_L/*, *m_I*/ ;
  LTA          *lta ;
  int          n_omp_threads;

  /* initialize the morph params */
  memset(&mp, 0, sizeof(GCA_MORPH_PARMS));
  /* for nonlinear morph */
  mp.l_jacobian = 1 ;
  mp.min_sigma = 0.4 ;
  mp.l_distance = 0 ;
  mp.l_log_likelihood = .025 ;
  mp.dt = 0.005 ;
  mp.noneg = True ;
  mp.exp_k = 20 ;
  mp.diag_write_snapshots = 1 ;
  mp.momentum = 0.9 ;
  if (FZERO(mp.l_smoothness))
  {
    mp.l_smoothness = 2 ;
  }

  // elastic stuff
  mp.lame_mu = 0.38462 ;
  mp.lame_lambda = 0.57692 ;
#if 0
  mp.l_smoothness = 0 ;
  mp.l_elastic = 0 ;
#endif

  mp.uncompress = 0 ;
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

  start.reset() ;
  setRandomSeed(-1L) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  Progname = argv[0] ;
  ac = argc ;
  av = argv ;
  mp.mri_diag = NULL;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
  {
    usage_exit(1) ;
  }

#ifdef HAVE_OPENMP
  n_omp_threads = omp_get_max_threads();
  printf("%d avail.processors, using %d\n",omp_get_num_procs(),omp_get_max_threads());
  //printf("\n\n ======= NUMBER OF OPENMP THREADS = %d ======= \n", n_omp_threads);
#else
  n_omp_threads = 1;
#endif

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

  if (mri_source->type == MRI_INT)
  {
    MRI *mri_tmp = MRIchangeType(mri_source, MRI_FLOAT, 0, 1, 1) ;
    MRIfree(&mri_source);
    mri_source = mri_tmp ;
  }
  mri_target = MRIread(target_fname) ;
  if (!mri_target)
    ErrorExit(ERROR_NOFILE, "%s: could not read target label volume %s",
              Progname, target_fname) ;
  if (mri_target->type == MRI_INT)
  {
    MRI *mri_tmp = MRIchangeType(mri_target, MRI_FLOAT, 0, 1, 1) ;
    MRIfree(&mri_target); mri_target = mri_tmp ;
  }

  if (mri_target_diag == NULL)
    mri_target_diag = mri_target ;
  if (mask_fname)
  {
    MRI *mri_mask = MRIread(mask_fname) ;
    if (mri_mask == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load mask from %s\n", Progname,mask_fname) ;
    
    MRImask(mri_source, mri_mask, mri_source, 0, 0) ;
    MRImask(mri_target, mri_mask, mri_target, 0, 0) ;
    MRIfree(&mri_mask) ;
  }

  if (erosions > 0)
  {
    int n ;
    for (n = 0 ; n < erosions ; n++)
    {
      MRIerodeZero(mri_target, mri_target) ;
      MRIerodeZero(mri_source, mri_source) ;
    }
  }
  if (!FEQUAL(scale_values, 1.0))
  {
    MRIscalarMul(mri_source, mri_source, scale_values) ;
    MRIscalarMul(mri_target, mri_target, scale_values) ;
  }
  if (transform && transform->type == MORPH_3D_TYPE)
  {
    TransformRas2Vox(transform, mri_source,NULL) ;
  }
  if (use_aseg == 0)
  {
    if (match_peak_intensity_ratio)
      MRImatchIntensityRatio(mri_source, mri_target, mri_source, .8, 1.2,
                             100, 125) ;
    else if (match_mean_intensity)
    {
      MRImatchMeanIntensity(mri_source, mri_target, mri_source) ;
    }
    MRIboundingBox(mri_source, 0, &box) ;
    pad = (int)ceil(PADVOX *
                    MAX(mri_target->xsize,
                        MAX(mri_target->ysize,mri_target->zsize)) /
                    MIN(mri_source->xsize,
                        MIN(mri_source->ysize,mri_source->zsize)));
#if 0
    {
      MRI *mri_tmp ;
      if (pad < 1)
      {
        pad = 1 ;
      }
      printf("padding source with %d voxels...\n", pad) ;
      mri_tmp = MRIextractRegionAndPad(mri_source, NULL, &box, pad) ;
      if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON)
      {
        MRIwrite(mri_tmp, "t.mgz") ;
      }
      MRIfree(&mri_source) ;
      mri_source = mri_tmp ;
    }
#endif
  }
  if (label_ignore_name)
  {
    char path[STRLEN], fname[STRLEN] ;
    LABEL *area ;
    FileNamePath(mri_target->fname, path) ;
    sprintf(fname, "%s/%s", path, label_ignore_name) ;
    area = LabelRead(NULL, fname) ;
    if (area == NULL)
    {
      ErrorExit(ERROR_NOFILE,
                "%s: could not load label from %s", Progname, fname) ;
    }
    LabelFillVolume(mri_target, area, 0) ;
    LabelFree(&area) ;
#if 0
    FileNamePath(mri_source->fname, path) ;
    sprintf(fname, "%s/%s", path, label_ignore_name) ;
    area = LabelRead(NULL, fname) ;
    if (area == NULL)
    {
      ErrorExit(ERROR_NOFILE,
                "%s: could not load label from %s", Progname, fname) ;
    }
    LabelFillVolume(mri_source, area, 0) ;
    LabelFree(&area) ;
#endif
  }
  mri_orig_source = MRIcopy(mri_source, NULL) ;

  mp.max_grad = 0.3*mri_source->xsize ;
  mp.max_grad = .1 ;

  if (transform == NULL)
  {
    transform = TransformAlloc(LINEAR_VOXEL_TO_VOXEL, NULL) ;
  }

  if (transform->type != MORPH_3D_TYPE)  // initializing m3d from a linear transform
  {
    new_transform = 1 ;
    lta = ((LTA *)(transform->xform)) ;
    if (lta->type != LINEAR_VOX_TO_VOX)
    {
      printf("converting ras xform to voxel xform\n") ;
      m_L = MRIrasXformToVoxelXform(mri_source, 
                                    mri_target,
                                    lta->xforms[0].m_L,
                                    NULL) ;
      MatrixFree(&lta->xforms[0].m_L) ;
      lta->type = LINEAR_VOX_TO_VOX ;
    }
    else
    {
      printf("using voxel xform\n") ;
      m_L = lta->xforms[0].m_L ;
    }
#if 0
    if (Gsx >= 0)   // update debugging coords
    {
      VECTOR *v1, *v2 ;

      v1 = VectorAlloc(4, MATRIX_REAL) ;
      Gsx -= (box.x-pad) ;
      Gsy -= (box.y-pad) ;
      Gsz -= (box.z-pad) ;
      V3_X(v1) = Gsx ;
      V3_Y(v1) = Gsy ;
      V3_Z(v1) = Gsz ;
      VECTOR_ELT(v1,4) = 1.0 ;
      v2 = MatrixMultiply(m_L, v1, NULL) ;

      Gsx = nint(V3_X(v2)) ;
      Gsy = nint(V3_Y(v2)) ;
      Gsz = nint(V3_Z(v2)) ;
      MatrixFree(&v2) ;
      MatrixFree(&v1) ;
      printf("mapping by transform (%d, %d, %d) --> (%d, %d, %d) "
             "for rgb writing\n",
             Gx, Gy, Gz, Gsx, Gsy, Gsz) ;
    }
#endif
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      write_snapshot(mri_target, mri_source, m_L, &mp, 0, 1, "linear_init");
    }

    lta->xforms[0].m_L = m_L;
    printf("initializing GCAM with vox->vox matrix:\n") ;
    MatrixPrint(stdout, m_L) ;
    gcam = GCAMcreateFromIntensityImage(mri_source, mri_target, transform) ;
#if 0
    gcam->gca = gcaAllocMax(1, 1, 1,
                            mri_target->width, mri_target->height,
                            mri_target->depth,
                            0, 0) ;
#endif
    GCAMinitVolGeom(gcam, mri_source, mri_target) ;
    if (label_dist_name)
    {
      char path[STRLEN], fname[STRLEN] ;
      LABEL *area ;

      FileNamePath(mri_target->fname, path) ;
      sprintf(fname, "%s/%s", path, label_dist_name) ;
      area = LabelRead(NULL, fname) ;
      if (area == NULL)
      {
        ErrorExit(ERROR_NOFILE,
                  "%s: could not load label from %s", Progname, fname) ;
      }
      GCAMsetStatus(gcam, GCAM_IGNORE_DISTANCES) ;
      mp.l_distance = 1 ;
      GCAMpreserveLabelMetricProperties(gcam, area, mri_target) ;
      LabelFree(&area) ;
    }
    if (use_aseg)
    {
      if (nowmsa)
	replace_wmsa(mri_source, mri_source) ;
      if (ribbon_name)
      {
        char fname[STRLEN], path[STRLEN];
	const char *str;
	const char *hemi ;
        int  h, s, label ;
        MRI_SURFACE *mris_white, *mris_pial ;
        MRI         *mri ;

        for (s = 0 ; s <= 1 ; s++) // source and target
        {
          if (s == 0)
          {
            str = source_surf ;
            mri = mri_source ;
            FileNamePath(mri->fname, path) ;
            strcat(path, "/../surf") ;
          }
          else
          {
            mri = mri_target ;
            FileNamePath(mri->fname, path) ;
            strcat(path, "/../elastic") ;
            str = target_surf ;
          }
          // sorry - these values come from FreeSurferColorLUT.txt
          MRIreplaceValueRange(mri, mri, 1000, 1034, Left_Cerebral_Cortex) ;
          MRIreplaceValueRange(mri, mri, 1100, 1180, Left_Cerebral_Cortex) ;
          MRIreplaceValueRange(mri, mri, 2000, 2034, Right_Cerebral_Cortex) ;
          MRIreplaceValueRange(mri, mri, 2100, 2180, Right_Cerebral_Cortex) ;
          for (h = LEFT_HEMISPHERE ; h <= RIGHT_HEMISPHERE ; h++)
          {
            if (h == LEFT_HEMISPHERE)
            {
              hemi = "lh" ;
              label = Left_Cerebral_Cortex ;
            }
            else
            {
              label = Right_Cerebral_Cortex ;
              hemi = "rh" ;
            }
            sprintf(fname, "%s/%s%s.white", path, hemi, str) ;
            mris_white = MRISread(fname) ;
            if (mris_white == NULL)
            {
              ErrorExit(ERROR_NOFILE,
                        "%s: could not read surface %s", Progname, fname) ;
            }
            MRISsaveVertexPositions(mris_white, WHITE_VERTICES) ;
            sprintf(fname, "%s/%s%s.pial", path, hemi, str) ;
            mris_pial = MRISread(fname) ;
            if (mris_pial == NULL)
            {
              ErrorExit(ERROR_NOFILE,
                        "%s: could not read surface %s", Progname, fname) ;
            }
            MRISsaveVertexPositions(mris_pial, PIAL_VERTICES) ;
            if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
            {
              sprintf(fname, "sb.mgz") ;
              MRIwrite(mri_source, fname) ;
              sprintf(fname, "tb.mgz") ;
              MRIwrite(mri_target, fname) ;
            }

            insert_ribbon_into_aseg(mri, mri, mris_white, mris_pial, h) ;
            if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
            {
              sprintf(fname, "sa.mgz") ;
              MRIwrite(mri_source, fname) ;
              sprintf(fname, "ta.mgz") ;
              MRIwrite(mri_target, fname) ;
            }
            MRISfree(&mris_white) ;
            MRISfree(&mris_pial) ;
          }
        }
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        {
          sprintf(fname, "s.mgz") ;
          MRIwrite(mri_source, fname) ;
          sprintf(fname, "t.mgz") ;
          MRIwrite(mri_target, fname) ;
        }
      }
      GCAMinitLabels(gcam, mri_target) ;
      GCAMsetVariances(gcam, 1.0) ;
    }
  }
  else  /* use a previously create morph and integrate it some more */
  {
    printf("using previously create gcam...\n") ;
    gcam = (GCA_MORPH *)(transform->xform) ;
    GCAMrasToVox(gcam, mri_source) ;
    if (label_dist_name)
    {
      char path[STRLEN], fname[STRLEN] ;
      LABEL *area ;

      FileNamePath(mri_target->fname, path) ;
      sprintf(fname, "%s/%s", path, label_dist_name) ;
      area = LabelRead(NULL, fname) ;
      if (area == NULL)
      {
        ErrorExit(ERROR_NOFILE,
                  "%s: could not load label from %s", Progname, fname) ;
      }
      GCAMsetStatus(gcam, GCAM_IGNORE_DISTANCES) ;
      mp.l_distance = 100;
      GCAMpreserveLabelMetricProperties(gcam, area, mri_target) ;
      LabelFree(&area) ;
    }
    if (use_aseg)
    {
      GCAMinitLabels(gcam, mri_target) ;
      GCAMsetVariances(gcam, 1.0) ;
    }
    else
    {
      GCAMaddIntensitiesFromImage(gcam, mri_target) ;
    }
  }
  if (gcam->width != mri_source->width ||
      gcam->height != mri_source->height ||
      gcam->depth != mri_source->depth)
    ErrorExit(ERROR_BADPARM,
              "%s: warning gcam (%d, %d, %d), "
              "doesn't match source vol (%d, %d, %d)",
              Progname, gcam->width, gcam->height, gcam->depth,
              mri_source->width, mri_source->height, mri_source->depth) ;

  if (mp.mri_diag == NULL)
    mp.mri_diag = mri_source ;
  mp.diag_morph_from_atlas = 0 ;
  mp.diag_write_snapshots = 1 ;
  mp.diag_sample_type = use_aseg ? SAMPLE_NEAREST : SAMPLE_TRILINEAR ;
  mp.diag_volume = use_aseg ? GCAM_LABEL : GCAM_MEANS ;

  if (renormalize)
  {
    GCAMnormalizeIntensities(gcam, mri_source) ;
  }
  if (mp.write_iterations != 0)
  {
    char fname[STRLEN] ;
    MRI  *mri_gca ;

    if (getenv("DONT_COMPRESS"))
    {
      sprintf(fname, "%s_target.mgh", mp.base_name) ;
    }
    else
    {
      sprintf(fname, "%s_target.mgz", mp.base_name) ;
    }
    if (mp.diag_morph_from_atlas == 0)
    {
      printf("writing target volume to %s...\n", fname) ;
      MRIwrite(mri_target_diag, fname) ;
      sprintf(fname, "%s_target", mp.base_name) ;
      MRIwriteImageViews(mri_target_diag, fname, IMAGE_SIZE) ;
    }
    else
    {
      if (use_aseg)
      {
        mri_gca = GCAMwriteMRI(gcam, NULL, GCAM_LABEL) ;
      }
      else
      {
        mri_gca = MRIclone(mri_source, NULL) ;
        GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
      }
      printf("writing target volume to %s...\n", fname) ;
      MRIwrite(mri_gca, fname) ;
      sprintf(fname, "%s_target", mp.base_name) ;
      MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
      MRIfree(&mri_gca) ;
    }
  }

  if (nozero)  // if negative will run once without them and once with them
  {
    printf("disabling zero nodes\n") ;
    GCAMignoreZero(gcam, mri_target) ;
  }
  mp.mri = mri_target ;
  if (mp.regrid == True && new_transform == 0)
  {
    GCAMregrid(gcam, mri_target, PAD, &mp, &mri_source) ;
  }

  mp.write_fname = out_fname ;
  if (rip && Gx >= 0)
  {
    int x, y, z ;
    for (x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
	for (z = 0 ;z < gcam->depth ; z++)
	  if (abs(x-Gx) > 1 || abs(y-Gy) > 1 || abs(z-Gz) > 1)
	    gcam->nodes[x][y][z].invalid = GCAM_POSITION_INVALID ;
  }
  if (handle_expanded_ventricles)
  {
    GCA_MORPH_PARMS old_parms ;
    int               start_t ;

    memmove(&old_parms, (const void *)&mp, sizeof(old_parms)) ;
    mp.tol = .01 ;
    mp.l_label = 0 ;
    mp.l_smoothness = .1 ;   // defaults to 10 when renormalizing by alignment
    mp.uncompress = 1 ;
    mp.ratio_thresh = .25;
    mp.navgs = 1024*4 ;
    mp.integration_type = GCAM_INTEGRATE_OPTIMAL ;
    mp.noneg = 0 ;
    printf("registering ventricular system...\n") ;
    GCAMregisterVentricles(gcam, mri_source, &mp) ;
    start_t = mp.start_t ;
    memmove(&mp, (const void *)&old_parms, sizeof(old_parms)) ;
    mp.start_t = start_t ;
  }
  if (nozero < 0)
    mp.enable_zero_passes = mp.npasses-1 ;
  GCAMregister(gcam, mri_source, &mp) ; // atlas is target, morph target into register with it
  
  if (apply_transform)
  {
    MRI *mri_aligned ;
    char   fname[STRLEN] ;

    FileNameRemoveExtension(out_fname, fname) ;
    strcat(fname, ".mgz") ;
    mri_aligned =
      GCAMmorphToAtlas(mp.mri, gcam, NULL, -1, mp.diag_sample_type) ;
    printf("writing transformed output volume to %s...\n", fname) ;
    MRIwrite(mri_aligned, fname) ;
    MRIfree(&mri_aligned) ;
  }
  printf("writing warp vector field to %s\n", out_fname) ;
  //GCAMvoxToRas(gcam) ;
  GCAMwrite(gcam, out_fname) ;
  //GCAMrasToVox(gcam, mri_source) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("mri_nl_align registration took %d minutes and %d seconds.\n",  minutes, seconds) ;
  printf("#VMPC# mri_nl_align done VmPeak  %d\n",GetVmPeak());
  exit(0) ;
  return(0) ;
}

extern int gcam_write_neg, gcam_write_grad ;

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!stricmp(option, "debug_voxel"))
  {
    Gsx = Gx = atoi(argv[2]) ;
    Gsy = Gy = atoi(argv[3]) ;
    Gsz = Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "bigventricles"))
  {
    handle_expanded_ventricles = 1 ;
    printf("handling expanded ventricles...\n") ;
  }
  else if (!stricmp(option, "nobigventricles"))
  {
    handle_expanded_ventricles = 0 ;
    printf("not handling expanded ventricles...\n") ;
  }
  else if (!stricmp(option, "wmsa"))
  {
    nowmsa = atoi(argv[2]) ;
    printf("%sabling WMSA labels\n", nowmsa ? "dis" : "en") ;
  }
  else if (!stricmp(option, "uncompress"))
  {
    mp.uncompress = 1 ;
    printf("parms.uncompress=1...\n") ;
  }
  else if (!stricmp(option, "write_neg"))
  {
    gcam_write_neg = 1 ;
    printf("writing map of negative nodes\n") ;
  }
  else if (!stricmp(option, "write_grad"))
  {
    gcam_write_grad = 1 ;
    printf("writing gradient maps\n") ;
  }
  else if (!stricmp(option, "debug_node"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "OPTIMAL"))
  {
    mp.integration_type = GCAM_INTEGRATE_OPTIMAL ;
    printf("using optimal time-step integration\n") ;
  }
  else if (!stricmp(option, "RIP"))
  {
    rip = 1 ;
    printf("ripping all nodes except one being debugged\n") ;
  }
  else if (!stricmp(option, "mu"))
  {
    mp.lame_mu = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting Lame mu = %2.3f\n", mp.lame_mu) ;
  }
  else if (!stricmp(option, "lambda"))
  {
    mp.lame_lambda = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting Lame lambda = %2.3f\n", mp.lame_lambda) ;
  }
  else if (!stricmp(option, "NONEG"))
  {
    mp.noneg = atoi(argv[2]) ;
    nargs = 1 ;
    if (mp.noneg >= 0)
      printf("%s allowing temporary folds during numerical minimization\n",
	     mp.noneg ? "not" : "") ;
    else
      printf("allowing folds during numerical minimization and removing them before termination\n") ;
  }
  else if (!stricmp(option, "renormalize"))
  {
    renormalize = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%srenormalizing intensities\n", renormalize ? "" : "not ");
    if (renormalize == 0)
    {
      match_mean_intensity = match_peak_intensity_ratio = 0 ;
    }
  }
  else if (!stricmp(option, "label"))
  {
    label_ignore_name = argv[2] ;
    nargs = 1 ;
    printf("ignoring voxels in label %s\n", label_ignore_name) ;
  }
  else if (!stricmp(option, "label_dist"))
  {
    label_dist_name = argv[2] ;
    nargs = 1 ;
    printf("preserving metric properties in label %s\n", label_dist_name) ;
  }
  else if (!stricmp(option, "aseg"))
  {
    match_mean_intensity = match_peak_intensity_ratio = 0 ;
    use_aseg = 1 ;
    mp.l_dtrans = 1 ;
    mp.l_log_likelihood = 0 ;
    renormalize = 0 ;
    printf("treating inputs as segmentations\n") ;
    mp.dtrans_labels = dtrans_labels ;
    mp.ndtrans = NDTRANS_LABELS ;
  }
  else if (!stricmp(option, "diag2"))
  {
    mp.mri_diag2 = MRIread(argv[2]) ;
    if (mp.mri_diag2 == NULL)
    {
      ErrorExit(ERROR_NOFILE,
                "%s: could not read diag volume from %s", Progname, argv[2]) ;
    }
    nargs = 1 ;
    printf("writing d2 diagnostics for input volume %s\n", argv[2]) ;
  }
  else if (!stricmp(option, "diag_target") || !stricmp(option, "target_diag"))
  {
    mri_target_diag = MRIread(argv[2]) ;
    if (mri_target_diag == NULL)
    {
      ErrorExit(ERROR_NOFILE,
                "%s: could not read target diag volume from %s", Progname, argv[2]) ;
    }
    nargs = 1 ;
    printf("writing target image using input volume %s\n", argv[2]) ;
  }
  else if (!stricmp(option, "diag"))
  {
    mp.mri_diag = MRIread(argv[2]) ;
    if (mp.mri_diag == NULL)
    {
      ErrorExit(ERROR_NOFILE,
                "%s: could not read diag volume from %s", Progname, argv[2]) ;
    }
    nargs = 1 ;
    printf("writing d2 diagnostics for input volume %s\n", argv[2]) ;
  }
  else if (!stricmp(option, "MOMENTUM") || !stricmp(option, "FIXED"))
  {
    mp.integration_type = GCAM_INTEGRATE_FIXED ;
    printf("using fixed time-step integration\n") ;
  }
  else if (!stricmp(option, "distance"))
  {
    distance = atof(argv[2]) ;
    nargs = 1 ;
    printf("expanding border by %2.1f mm every outer cycle\n", distance);
  }
  else if (!stricmp(option, "dtrans"))
  {
    mp.l_dtrans = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting distance transform coefficient to %2.3f\n", mp.l_dtrans) ;
  }
  else if (!stricmp(option, "threads"))
  {
    int nthreads;
    sscanf(argv[2],"%d",&nthreads);
    setenv("OMP_NUM_THREADS",argv[2],1);
    #ifdef _OPENMP
    omp_set_num_threads(nthreads);
    #endif
    printf("setting threads to %d\n",nthreads);
    nargs = 1 ;
  }
  else if (!stricmp(option, "match_peak"))
  {
    match_peak_intensity_ratio = 1 ;
    match_mean_intensity = 0 ;
    printf("matching peak of intensity ratio histogram\n") ;
  }
  else if (!stricmp(option, "erode"))
  {
    erosions = atoi(argv[2]) ;
    nargs = 1 ;
    printf("eroding source and target image %d times before morphing\n",
           erosions) ;
  }
  else if (!stricmp(option, "match_mean"))
  {
    match_peak_intensity_ratio = 0 ;
    match_mean_intensity = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%smatching peak of intensity ratio histogram\n",
           match_mean_intensity ? "" : "not ") ;
  }
  else if (!stricmp(option, "intensity") ||!stricmp(option, "ll"))
  {
    mp.l_log_likelihood = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting l_log_likelihood = %2.3f\n", mp.l_log_likelihood );
  }
  else if (!stricmp(option, "likelihood"))
  {
    mp.l_likelihood = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting l_likelihood = %2.3f\n", mp.l_likelihood );
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
  else if (!stricmp(option, "area_smoothness") || !stricmp(option, "asmooth"))
  {
    mp.l_area_smoothness = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_area_smoothness=%2.3f\n", mp.l_area_smoothness) ;
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
  else if (!stricmp(option, "min_sigma"))
  {
    mp.min_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("using min sigma=%2.3f\n", mp.min_sigma) ;
  }
  else if (!stricmp(option, "ribbon"))
  {
    ribbon_name = argv[2] ;
    printf("reading ribbon from %s and inserting into aseg\n", ribbon_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "rthresh"))
  {
    mp.ratio_thresh = atof(argv[2]) ;
    mp.uncompress = 1 ;
    nargs = 1 ;
    printf("using compression ratio threshold = %2.3f...\n", mp.ratio_thresh) ;
  }
  else if (!stricmp(option, "scale"))
  {
    scale_values = atof(argv[2]) ;
    nargs = 1 ;
    printf("scaling input values by %2.3f\n", scale_values) ;
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
  else if (!stricmp(option, "mask"))
  {
    mask_fname = (argv[2]);
    printf("masking inputs with %s\n", mask_fname) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "debug_label"))
  {
    Gdiag_no = atoi(argv[2]);
    printf("debugging label %s (%d)\n", cma_label_to_name(Gdiag_no), Gdiag_no) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "apply"))
  {
    apply_transform = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%sapplying transform after registration\n",
           apply_transform ? "" : "not ") ;
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
    case 'E':
      mp.l_elastic = atof(argv[2]) ;
      nargs = 1 ;
      printf("using l_elastic = %2.3f\n", mp.l_elastic) ;
      break ;
    case 'T':
      printf("reading transform from %s...\n", argv[2]) ;
      transform = TransformRead(argv[2]) ;
      if (transform == NULL)
      {
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read transform from %s\n",
                  Progname,argv[2]);
      }
      nargs = 1 ;
      if (transform->type == LINEAR_VOX_TO_VOX)
      {
        printf("converting transform to ras....\n") ;
        LTAvoxelToRasXform((LTA *)(transform->xform), NULL, NULL) ;
      }
      break ;
    case 'I':
      printf("reading transform from %s...\n", argv[2]) ;
      transform = TransformRead(argv[2]) ;
      if (transform == NULL)
      {
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read transform from %s\n",
                  Progname,argv[2]);
      }
      TransformInvert(transform, NULL) ;
      TransformSwapInverse(transform) ;
      if (transform->type == LINEAR_VOX_TO_VOX)
      {
        printf("converting transform to ras....\n") ;
        LTAvoxelToRasXform((LTA *)(transform->xform), NULL, NULL) ;
      }
      nargs = 1 ;
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
    case 'Z':
    {
      int i ;
      i = atoi(argv[2]) ;
	
      if (i > 0)
	nozero = 0 ;
      else if (i < 0)
	nozero = -1 ;
      else
	nozero = 1 ;
		
      if (nozero >= 0)
	printf("%sdisabling zero image locations\n", nozero ? "" : "not ") ;
      else
	printf("disabling zero image locations for initial run, then enabling\n") ;

      nargs = 1 ;
      break ;
    }
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
      if (mp.write_iterations > 0)
      {
        Gdiag |= DIAG_WRITE ;
      }
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
  printf("usage: %s [options] <source> <target> <warp>\n",
         Progname) ;
  printf("\n") ;
  printf("Options:\n\n") ;
  printf("  -debug_voxel Gx Gy Gz\n");
  printf("  -debug_node Gx Gy Gz\n");
  printf("  -noneg val		: val = 0 or 1; <not> allowing temporary folds during numerical minimization (default = True)\n");
  printf("  -renormalize val	: val = 0 or 1; <not> renormalizing intensities\n");
  printf("  -aseg			: treating inputs as segmentations; no intenisty renormalization; log_likelihood = 0; distance transform coefficient (dtrans) = 1\n");
  printf("  -diag2 volname	: writing d2 diagnostics for input volume (*volname*)\n");
  printf("  -OPTIMAL		: using line search optimization (default = GCAM_INTEGRATE_BOTH)\n");
  printf("  -MOMENTUM		: using fixed time-step integration (default = GCAM_INTEGRATE_BOTH)\n");
  printf("  -FIXED		: using fixed time-step integration (default = GCAM_INTEGRATE_BOTH)\n");
  printf("  -distance val		: expanding border by *val* mm every outer cycle\n");
  printf("  -dtrans val		: setting distance transform coefficient to *val*\n");
  printf("  -match_peak		: matching peak of intensity ratio histogram\n");
  printf("  -erode val		: eroding source and target image *val* times before morphing\n");
  printf("  -match_mean val	: val = 0 or 1; <not> matching peak of intensity ratio histogram\n");
  printf("  -intensity val	: setting l_log_likelihood to *val* (default = 0.025)\n");
  printf("  -ll val		: setting l_log_likelihood to *val* (default = 0.025)\n");
  printf("  -noregrid		: disable regridding\n");
  printf("  -regrid		: enable regridding\n");
  printf("  -view Gx Gy Gz	: viewing voxel (Gx, Gy, Gz)\n");
  printf("  -levels val		: setting levels to *val* (default = 6)\n");
  printf("  -areasmoothness val	: setting l_area_smoothness to *val*\n");
  printf("  -asmooth val		: setting l_area_smoothness to *val*\n");
  printf("  -area val		: seeting l_area to *val*\n");
  printf("  -tol val		: setting tol to *val* (default = .1)\n");
  printf("  -sigma val		: seeting sigma to *val* (default = 8)\n");
  printf("  -min_sigma val	: setting minimum value for sigma to be *val* (default = .4)\n");
  printf("  -ribbon fname		: reading ribbon from *fname* and inserting into aseg\n");
  printf("  -rthresh val		: setting compression ratio threshold to *val*\n");
  printf("  -scale val		: scaling input values by *val*\n");
  printf("  -dt val		: setting dt to *val* (default = 0.005)\n");
  printf("  -passes	val	: integrating in *val* number of passes / go through all levels *val* number of times (default=3)\n");
  printf("  -skip val		: skipping *val* number of voxels in source data\n");
  printf("  -apply val		: val = 0 or 1; <not> applying transform after registration\n");
  printf("  -D val		: setting l_distance to *val* (default = 0)\n");
  printf("  -M val		: setting momentum to *val* (default = .9)\n");
  printf("  -N val		: setting number of iterations to *val* (default = 1000)\n");
  printf("  -s val		: setting l_smoothness to *val* (default = 2)\n");
  printf("  -T fname		: reading the forward transform / deformation from *fname*\n");
  printf("  -I fname		: reading the inverse of the transform / deformation in *fname* \n");
  printf("  -B val		: setting l_binary to *val*\n");
  printf("  -J val		: setting l_jacobian to *val* (default = 1)\n");
  printf("  -Z val		: val = 0 or 1; <not> disabling zero image locations\n");
  printf("  -a val		: smoothing gradient with *val* number of averages\n");
  printf("  -K val		: setting exp_k to val (default=20.0)\n");
  printf("  -W val		: write diagnostics at each *val* iteration\n");
  printf("  -?			: print usage\n");
  printf("  -U			: print usage \n");

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
  if (1 || conform)
  {
    mri_aligned = MRIalloc(mri_target->width, mri_target->height,
                           mri_target->depth,mri_source->type);
    MRIcopyHeader(mri_target, mri_aligned) ;
    MRIlinearTransformInterp(mri_source,
                             mri_aligned,
                             m_vox_xform,
                             SAMPLE_NEAREST);
  }
  else
  {
    mri_aligned = MRITransformedCenteredMatrix(mri_source,
                                               mri_target,
                                               m_vox_xform) ;
  }
  if (in_fname)
  {
    sprintf(fname, "%s_%s", parms->base_name, in_fname) ;
  }
  else
  {
    sprintf(fname, "%s_%03d", parms->base_name, fno) ;
  }
  MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
  if (in_fname)
  {
    sprintf(fname, "%s_%s.mgz", parms->base_name, in_fname) ;
  }
  else
  {
    sprintf(fname, "%s_%03d.mgz", parms->base_name, fno) ;
  }
  printf("writing snapshot to %s...\n", fname) ;
  MRIwrite(mri_aligned, fname) ;
  MRIfree(&mri_aligned) ;

  {
#if 0
    mri_aligned = MRIsrcTransformedCentered(mri_source,
                                            mri_target,
                                            m_vox_xform,
                                            SAMPLE_NEAREST) ;
#else
    mri_aligned = MRITransformedCenteredMatrix(mri_source,
                                               mri_target,
                                               m_vox_xform) ;
#endif
    if (in_fname)
    {
      sprintf(fname, "orig_%s_%s.mgz", parms->base_name, in_fname) ;
    }
    else
    {
      sprintf(fname, "orig_%s_%03d.mgz", parms->base_name, fno) ;
    }
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri_aligned, fname) ;
    MRIfree(&mri_aligned) ;
  }

  return(NO_ERROR) ;
}


