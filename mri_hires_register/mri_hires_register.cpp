/*
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
// mri_hires_register.c
//
// written by Bruce Fischl
// Nov. 9th ,2000
//
// Warning: Do not edit the following four lines.  CVS maintains them.
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
#include "numerics.h"
#include "fastmarching.h"
#include "voxlist.h"

#define DEFAULT_MAX_ANGLE       RADIANS(25)
static double MAX_ANGLE = DEFAULT_MAX_ANGLE ;
static double MAX_SCALE = 0.5 ;

#define HIRES_PAD       10
#define LOWRES_PAD      20

static int find_gcam_node(GCA_MORPH *gcam, int label,
                          float x_vox, float y_vox, float z_vox) ;

static int angio = 0 ;
static int find_label = -1 ;
static float x_vox = 0.0 ;
static float y_vox = 0.0 ;
static float z_vox = 0.0 ;

static int apply_transform = 1 ;

static int fix_intensity = 0 ;
static MRI *estimate_densities(GCA_MORPH *gcam,
                               MRI *mri_lowres,
                               MRI *mri_intensity) ;
static int write_snapshot(MRI *mri_lowres, MRI *mri_hires,
                          MATRIX *m_vox_xform, INTEGRATION_PARMS *parms,
                          int fno, int conform, const char *fname) ;

static double MAX_TRANS = 30 ;
static int regrid = 0 ;

static float compute_powell_sse(float *p) ;
static int powell_minimize(VOXEL_LIST *vl_lowres,
                           VOXEL_LIST *vl_hires,
                           MATRIX *mat);

static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;

static TRANSFORM *compute_optimal_transform(VOXEL_LIST *vl_lowres,
    VOXEL_LIST *vl_hires,
    INTEGRATION_PARMS *parms,
    TRANSFORM *transform) ;

static double compute_overlap(VOXEL_LIST *vl_lowres,
                              VOXEL_LIST *vl_hires,
                              MATRIX *m_L) ;
static double compute_distance_transform_sse(VOXEL_LIST *vl_lowres,
    VOXEL_LIST *vl_hires,
    MATRIX *m_L) ;

static double (*pf_overlap)(VOXEL_LIST *vl_lowres,
                            VOXEL_LIST *vl_hires,
                            MATRIX *m_L) = compute_overlap ;

static double find_optimal_translation(VOXEL_LIST *vl_lowres,
                                       VOXEL_LIST *vl_hires,
                                       MATRIX *m_L, float min_trans,
                                       float max_trans,
                                       float trans_steps, int nreductions);

static double find_optimal_linear_xform(VOXEL_LIST *vl_lowres,
                                        VOXEL_LIST *vl_hires,
                                        MATRIX *m_L,
                                        MATRIX *m_origin,
                                        float min_angle, float max_angle,
                                        float min_scale, float max_scale,
                                        float min_trans, float max_trans,
                                        float angle_steps, float scale_steps,
                                        float trans_steps,
                                        int nreductions);
const char *Progname ;
static int target_label = Right_Hippocampus ;

static int skip = 2 ;
static double distance = 1.0 ;

static MRI *mri_hires_intensity = NULL ;
static char *hires_intensity_fname = NULL ;


static int non_hippo_labels[] = {
                                  lateral_ventricle,
                                  entorhinal_cortex,
                                  Amygdala,
                                  Cerebral_White_Matter,
                                  Cerebral_Cortex,
                                  Inf_Lat_Vent
                                } ;
#define NUM_NON_HIPPO_LABELS  (sizeof(non_hippo_labels) / sizeof(non_hippo_labels[0]))

static int non_artery_labels[] = {
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


static INTEGRATION_PARMS parms ;
static TRANSFORM  *transform = NULL ;
static GCA_MORPH_PARMS mp ;

int
main(int argc, char *argv[]) {
  char       **av, *hires_fname, *aseg_fname,
  *intensity_fname, *out_fname, fname[STRLEN] ;
  int        ac, nargs, i, new_transform = 0 ;
  MRI        *mri_intensity, *mri_lowres, *mri_hires,
  *mri_tmp, *mri_target, *mri_dist_src = NULL,
                                         *mri_dist_dst = NULL ;
  VOXEL_LIST *vl_lowres, *vl_hires ;
  MRI_REGION  box ;
  Timer start ;
  int          msec, minutes, seconds, label ;

  parms.write_iterations = 0 ;
  parms.start_t = 0 ;

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
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
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

  if (angio)  // only use arterial labels for morph
  {
    pf_overlap = compute_distance_transform_sse ;
    //    if (transform == NULL)
    {
      for (i = 0 ; i < NUM_NON_ARTERY_LABELS ; i++) {
        label = non_artery_labels[i] ;
        MRIreplaceValues(mri_hires, mri_hires, label, 0) ;
      }
    }
  }
  else if (transform == NULL) /* remove non-hippo labels
                                                 if only doing linear morph */
  {
    for (i = 0 ; i < NUM_NON_HIPPO_LABELS ; i++) {
      label = non_hippo_labels[i] ;
      MRIreplaceValues(mri_hires, mri_hires, label, 0) ;
    }
  }

  if (mri_hires->type != MRI_UCHAR && angio == 0) {
    mri_tmp = MRIchangeType(mri_hires, MRI_UCHAR, 0, 255, 1) ;
    MRIfree(&mri_hires) ;
    mri_hires = mri_tmp ;
  }

  MRIboundingBox(mri_hires, 0, &box) ;
  box.x -= HIRES_PAD ;
  box.y -= HIRES_PAD ;
  box.z -= HIRES_PAD ;
  box.dx += 2*HIRES_PAD ;
  box.dy += 2*HIRES_PAD ;
  box.dz += 2*HIRES_PAD ;
  MRIcropBoundingBox(mri_hires, &box) ;
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
  if (angio) {
    MRIbinarize(mri_lowres, mri_lowres, 1, 0, 128) ;
    MRIbinarize(mri_hires, mri_hires, 1, 0, 128) ;
    MRIwrite(mri_lowres, "target_labels.mgz") ;
    MRIwrite(mri_hires, "src_labels.mgz") ;
#if 1
    if (mri_lowres->type != MRI_UCHAR) {
      mri_tmp = MRIchangeType(mri_lowres, MRI_UCHAR, 0, 255, 1) ;
      MRIfree(&mri_lowres) ;
      mri_lowres = mri_tmp ;
    }
    if (mri_hires->type != MRI_UCHAR) {
      mri_tmp = MRIchangeType(mri_hires, MRI_UCHAR, 0, 255, 1) ;
      MRIfree(&mri_hires) ;
      mri_hires = mri_tmp ;
    }
#endif
    target_label = 128 ;
    printf("creating distance transforms...\n") ;
#if 0
    mri_dist =
      MRIextractDistanceMap(mri_lowres,
                            NULL, 128,
                            2*MAX(MAX(mri_lowres->width,
                                      mri_lowres->height),
                                  mri_lowres->depth),1);
    MRIwrite(mri_dist, "dist.mgz") ;
#else
    if (transform == NULL) {
      mri_dist_src =
        MRIdistanceTransform(mri_hires, NULL, 128, -1, DTRANS_MODE_SIGNED, NULL);
      mri_dist_dst =
        MRIdistanceTransform(mri_lowres, NULL,128, -1, DTRANS_MODE_SIGNED, NULL);
      MRIwrite(mri_dist_src, "dist_src.mgz") ;
      MRIwrite(mri_dist_dst, "dist_dst.mgz") ;
    }
#endif
  }


  if (Gdiag & DIAG_WRITE && parms.write_iterations > 0) {
    int req = snprintf(fname, STRLEN, "%s_target", parms.base_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwriteImageViews(mri_lowres, fname, IMAGE_SIZE) ;
    req = snprintf(fname, STRLEN, "intensity_%s_target", parms.base_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwriteImageViews(mri_intensity, fname, IMAGE_SIZE) ;
    MRIwrite(mri_intensity, "intensity_target.mgz") ;
    MRIwrite(mri_lowres, "aseg_target.mgz") ;
  }

  if (hires_intensity_fname) {
    MRI *mri_tmp2 = MRIread(hires_intensity_fname) ;
    if (!mri_tmp2)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not read intensity volume %s...\n",
       hires_intensity_fname) ;
    if (mri_tmp2->type != MRI_UCHAR && mri_tmp2->type != MRI_FLOAT) {
      mri_hires_intensity =
        MRIchangeType(mri_tmp2, MRI_FLOAT, 0, 2000, 1) ;
      MRIfree(&mri_tmp2) ;
      mri_tmp2 = mri_hires_intensity ;
    }
  }

  if (transform == NULL)   /* compute optimal linear transform */
  {
    vl_lowres = VLSTcreate(mri_lowres,target_label,target_label,NULL,skip,0);
    vl_lowres->mri2 = mri_dist_dst ;
    for (i = 0 ; i < 3 ; i++) {
      printf("------------- outer loop iteration %d ---------------\n",i) ;
      vl_hires = VLSTcreate(mri_hires, 1, 255, NULL, skip, 0) ;
      vl_hires->mri2 = mri_dist_src ;

      transform = compute_optimal_transform(vl_lowres, vl_hires, &parms,
                                            transform) ;
      VLSTfree(&vl_hires) ;
      vl_hires = VLSTcreate(mri_hires, 1, 255, NULL, skip/4, 0) ;
      vl_hires->mri2 = mri_dist_src ;
      powell_minimize(vl_lowres, vl_hires,
                      ((LTA *)(transform->xform))->xforms[0].m_L) ;
      VLSTfree(&vl_hires) ;
      if (apply_transform) {
        MRI *mri_aligned ;
        MATRIX *m_vox_xform ;
        char   fname[STRLEN] ;

        FileNameRemoveExtension(out_fname, fname) ;
        strcat(fname, ".mgz") ;
        m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;
        mri_aligned = MRIclone(mri_lowres, NULL) ;
        MRIlinearTransformInterp
        (mri_hires, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
        printf("writing transformed output volume to %s...\n", fname) ;
        MRIwrite(mri_aligned, fname) ;
        MRIfree(&mri_aligned) ;
      }
      LTAvoxelToRasXform
      ((LTA *)(transform->xform), mri_hires, mri_lowres) ;
      TransformWrite(transform, out_fname) ;
      LTArasToVoxelXform
      ((LTA *)(transform->xform), mri_hires, mri_lowres) ;
      skip /= 2 ;
    }
    VLSTfree(&vl_lowres) ;

    printf("final vox2vox matrix:\n") ;
    MatrixPrint(stdout, ((LTA *)(transform->xform))->xforms[0].m_L) ;
    {
      MRI *mri_aligned, *mri_filtered ;
      char fname[STRLEN] ;
      int  i ;

#if 0
      mri_aligned =
        MRIsrcTransformedCentered
        (mri_hires, mri_lowres,
         ((LTA *)(transform->xform))->xforms[0].m_L, SAMPLE_NEAREST) ;
#else
      mri_aligned =
        MRITransformedCenteredMatrix
        (mri_hires, mri_lowres, ((LTA *)(transform->xform))->xforms[0].m_L) ;
#endif
      int req = snprintf(fname, STRLEN, "%sfinal.mgz", parms.base_name) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      MRIwrite(mri_aligned, fname) ;

      for (i = 1 ; i <= 10 ; i++) {
        int req = snprintf(fname, STRLEN, "%sfiltered%d.mgz", parms.base_name, i) ;
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        mri_filtered = MRImodeFilter(mri_aligned, NULL, i) ;
        printf("writing filtered image to %s\n", fname) ;
        MRIwrite(mri_filtered, fname) ;
        MRIfree(&mri_filtered) ;
      }
      MRIfree(&mri_aligned) ;
    }

    if (apply_transform) {
      MRI *mri_aligned ;
      MATRIX *m_vox_xform ;
      char   fname[STRLEN] ;

      FileNameRemoveExtension(out_fname, fname) ;
      strcat(fname, ".mgz") ;
      m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;
      mri_aligned = MRIclone(mri_lowres, NULL) ;
      MRIlinearTransformInterp
      (mri_hires, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
      printf("writing transformed output volume to %s...\n", fname) ;
      MRIwrite(mri_aligned, fname) ;
      MRIfree(&mri_aligned) ;
    }
    LTAvoxelToRasXform((LTA *)(transform->xform), mri_hires, mri_lowres) ;
    TransformWrite(transform, out_fname) ;
  } else   /* linear transform already computed - compute 3d morph */
  {
    GCA_MORPH *gcam ;
    MATRIX    *m_L, *m_I ;
    LTA       *lta ;

    strcpy(mp.base_name, parms.base_name) ;
    mp.max_grad = 0.3*mri_hires->xsize ;

    if (transform->type != MORPH_3D_TYPE)
      // initializing m3d from a linear transform
    {
      new_transform = 1 ;
      lta = ((LTA *)(transform->xform)) ;
      m_L =
        MRIrasXformToVoxelXform
        (mri_hires, mri_lowres, lta->xforms[0].m_L, NULL) ;
      MatrixFree(&lta->xforms[0].m_L) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        write_snapshot
        (mri_lowres, mri_hires, m_L, &parms, 0, 1,"linear_init");

      mri_tmp = MRITransformedCenteredMatrix(mri_hires, mri_lowres, m_L) ;
      MRIfree(&mri_hires) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_tmp, "b.mgz") ;
      mri_hires = MRImodeFilter(mri_tmp, NULL, 3) ;
      MRIfree(&mri_tmp) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_hires, "a.mgz") ;
      m_I = MatrixIdentity(4, NULL) ;
      MRIrasXformToVoxelXform(mri_hires, mri_lowres, m_I, m_L);
      MatrixFree(&m_I) ;

#if 0
      /* make sure none of the labels are on the border */
      MRIboundingBox(mri_hires, 0, &box) ;
      box.x -= HIRES_PAD ;
      box.y -= HIRES_PAD ;
      box.z -= HIRES_PAD ;
      box.dx += 2*HIRES_PAD ;
      box.dy += 2*HIRES_PAD ;
      box.dz += 2*HIRES_PAD ;
      mri_tmp = MRIextractRegion(mri_hires, NULL, &box) ;
      MRIfree(&mri_hires) ;
      mri_hires = mri_tmp ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_hires, "hires_xformed.mgz") ;
#endif

      lta->xforms[0].m_L = m_L ;
      printf("initializing GCAM with vox->vox matrix:\n") ;
      MatrixPrint(stdout, m_L) ;
      gcam =
        GCAMalloc(mri_hires->width, mri_hires->height, mri_hires->depth);

      GCAMinit(gcam, mri_lowres, NULL, transform, 0) ;
      GCAMinitVolGeom(gcam, mri_lowres, mri_hires) ;
    } else  /* use a previously create morph and integrate it some more */
    {
      printf("using previously create gcam...\n") ;
      gcam = (GCA_MORPH *)(transform->xform) ;
      if (FZERO(mp.l_area_intensity) && FZERO(mp.l_log_likelihood))
        GCAMrasToVox(gcam, mri_lowres) ;
      else
        GCAMrasToVox(gcam, mri_intensity) ;
      /*      GCAMremoveCompressedRegions(gcam, 0.1) ;*/
    }
    if (gcam->width != mri_hires->width ||
        gcam->height != mri_hires->height ||
        gcam->depth != mri_hires->depth)
      ErrorExit
      (ERROR_BADPARM,
       "%s: warning gcam (%d, %d, %d), "
       "doesn't match hires vol (%d, %d, %d)",
       Progname, gcam->width, gcam->height, gcam->depth,
       mri_hires->width, mri_hires->height, mri_hires->depth) ;
    mp.diag_morph_from_atlas = 1 ;
    mp.diag_volume = GCAM_LABEL ;
    mp.diag_mode_filter = 1 ;
    if (0 && regrid) {
      double pct_change ;
      int  niter = 0, nlevels = mp.levels,
                                level, npasses, navgs = mp.navgs,
                                                        level_steps, num_this_scale ;
      mp.levels = 1 ;

      for (npasses = 0 ; npasses < 3 ; npasses++) {
        mp.navgs = navgs ;
        for (level = nlevels ; level >= 0 ; level--) {
          num_this_scale = 0 ;
          do {
            GCAMinitLabels(gcam, mri_hires) ;
            if (!FZERO(mp.l_area_intensity) ||
                !FZERO(mp.l_log_likelihood))
              estimate_densities(gcam, mri_lowres, mri_intensity) ;

            if (distance > 0) {
              printf("expanding GCAM border by %2.3f mm\n",
                     distance) ;
              GCAMexpand(gcam, distance) ;
            }
            level_steps = mp.start_t ;
            GCAMregister(gcam, mri_lowres, &mp) ;
            level_steps = mp.start_t - level_steps ;
            GCAMvoxToRas(gcam) ;
            GCAMwrite(gcam, out_fname) ;
            GCAMrasToVox(gcam, mri_lowres) ;

            pct_change =
              100.0*(mp.start_rms - mp.end_rms) / (mp.start_rms) ;
            niter++ ;

            /* regrid */
            MRIfree(&mri_hires) ;
            gcam =
              GCAMregrid
              (gcam, mri_lowres, HIRES_PAD, &mp, &mri_hires) ;
            printf("outer loop iteration %d completed "
                   "%d steps (%d this scale), "
                   "rms = %2.3f (%2.1f%%) (was %2.1f)\n",
                   niter, level_steps, num_this_scale+1,
                   mp.end_rms, pct_change, mp.start_rms) ;
            if (level_steps == 0 || ++num_this_scale > 10)
              break ;
          } while (pct_change > level_steps*mp.tol) ;
          mp.navgs /= 4 ;
        }
      }
      mp.integration_type = GCAM_INTEGRATE_BOTH ;
      mp.navgs = navgs ;
      mp.tol *= .1 ;
      mp.levels = nlevels ;
      GCAMregister(gcam, mri_lowres, &mp) ;
    } else  /* don't regrid */
    {
      if (gcam->status == GCAM_UNLABELED)
        GCAMinitLabels(gcam, mri_hires) ;

      if (!DZERO(mp.l_area_intensity) || !DZERO(mp.l_log_likelihood)) {
        mp.mri_binary =
          estimate_densities(gcam, mri_lowres, mri_intensity) ;
        MRIfree(&mri_lowres) ;
        mri_lowres = mp.mri_binary ;
        mri_target = mri_intensity ;
      } else {
        mp.mri_binary = MRIcopy(mri_lowres, NULL) ;
        mri_target = mri_lowres ;
      }

      /* remove other labels from segmentation image */
      mri_tmp = MRIclone(mp.mri_binary, NULL) ;
      MRIcopyLabel(mp.mri_binary, mri_tmp, target_label) ;
      MRIfree(&mp.mri_binary) ;
      mp.mri_binary = mri_tmp ;

      /*      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)*/
      MRIwrite(mp.mri_binary, "aseg_target.mgz") ;
      if (angio)
        // doesn't do anything as the volume was binarized to 128 before
      {
        for (i = 0 ; i < NUM_NON_ARTERY_LABELS ; i++) {
          label = non_artery_labels[i] ;
          GCAMsetLabelStatus(gcam, label, GCAM_BINARY_ZERO) ;
        }
      }
      else {
        for (i = 0 ; i < NUM_NON_HIPPO_LABELS ; i++) {
          label = non_hippo_labels[i] ;
          GCAMsetLabelStatus(gcam, label, GCAM_BINARY_ZERO) ;
        }
      }

      GCAMsetLabelStatus(gcam, 0, GCAM_NEVER_USE_LIKELIHOOD) ;
      mp.mri = mri_target ;
      if (mp.regrid == True && new_transform == 0)
        GCAMregrid(gcam, mri_target, HIRES_PAD, &mp, &mri_hires) ;
      if (find_label >= 0)
        find_gcam_node(gcam, find_label, x_vox, y_vox, z_vox) ;
      GCAMregister(gcam, mri_target, &mp) ;
      if (apply_transform) {
        MRI *mri_aligned ;
        char   fname[STRLEN] ;

        FileNameRemoveExtension(out_fname, fname) ;
        strcat(fname, ".mgz") ;
        mri_aligned =
          GCAMmorphFieldFromAtlas(gcam, mp.mri, GCAM_LABEL,0, 1) ;
        printf("writing transformed output volume to %s...\n", fname) ;
        MRIwrite(mri_aligned, fname) ;
        MRIfree(&mri_aligned) ;
      }
      GCAMvoxToRas(gcam) ;
      GCAMwrite(gcam, out_fname) ;
      GCAMrasToVox(gcam, mri_target) ;
    }
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  printf("#VMPC# mri_hires_register VmPeak  %d\n",GetVmPeak());
  printf("registration took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  exit(0) ;
  return(0) ;
}

static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "FIX")) {
    fix_intensity = 1 ;
    printf("using predefined intensities for class means...\n") ;
  } else if (!stricmp(option, "OPTIMAL")) {
    mp.integration_type = GCAM_INTEGRATE_OPTIMAL ;
    printf("using optimal time-step integration\n") ;
  } else if (!stricmp(option, "find_label")) {
    find_label = atoi(argv[2]) ;
    x_vox = atof(argv[3]) ;
    y_vox = atof(argv[4]) ;
    z_vox = atof(argv[5]) ;
    nargs = 4 ;
    printf("finding label %s (%d) at (%2.1f, %2.1f, %2.1f)\n",
           cma_label_to_name(find_label), find_label, x_vox, y_vox,z_vox) ;
  } else if (!stricmp(option, "MOMENTUM") || !stricmp(option, "FIXED")) {
    mp.integration_type = GCAM_INTEGRATE_FIXED ;
    printf("using optimal time-step integration\n") ;
  } else if (!stricmp(option, "distance")) {
    distance = atof(argv[2]) ;
    nargs = 1 ;
    printf("expanding border by %2.1f mm every outer cycle\n", distance);
  } else if (!stricmp(option, "intensity") ||!stricmp(option, "ll")) {
    mp.l_log_likelihood = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting l_log_likelihood = %2.1f\n", mp.l_log_likelihood );
  } else if (!stricmp(option, "noregrid")) {
    regrid = 0 ;
    mp.regrid = False ;
    printf("disabling regridding...\n") ;
  } else if (!stricmp(option, "angio")) {
    angio = 1 ;
    printf("assuming inputs are vascular labelings...\n") ;
  } else if (!stricmp(option, "regrid")) {
    regrid = 1 ;
    mp.regrid = True ;
    printf("enabling regridding...\n") ;
  } else if (!stricmp(option, "view")) {
    Gsx = atoi(argv[2]) ;
    Gsy = atoi(argv[3]) ;
    Gsz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("viewing voxel (%d, %d, %d)\n", Gsx, Gsy, Gsz) ;
  } else if (!stricmp(option, "LEVELS")) {
    mp.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", mp.levels) ;
  } else if (!stricmp(option, "trans")) {
    MAX_TRANS = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting MAX_TRANS = %2.1f\n", MAX_TRANS) ;
  } else if (!stricmp(option, "area")) {
    mp.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_area=%2.3f\n", mp.l_area) ;
  } else if (!stricmp(option, "area_intensity")) {
    mp.l_area_intensity = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_area_intensity=%2.3f\n", mp.l_area_intensity) ;
  } else if (!stricmp(option, "max_angle")) {
    MAX_ANGLE = RADIANS(atof(argv[2])) ;
    nargs = 1 ;
    printf("setting max angle for search to %2.1f degrees\n",
           DEGREES(MAX_ANGLE)) ;
  } else if (!stricmp(option, "max_scale")) {
    MAX_SCALE = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting max scale for search to %2.1f --> %2.1f\n",
           1-MAX_SCALE, 1+MAX_SCALE) ;
  } else if (!stricmp(option, "tol")) {
    mp.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("using tol=%2.3f\n", mp.tol) ;
  } else if (!stricmp(option, "sigma")) {
    mp.sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("using sigma=%2.3f\n", mp.sigma) ;
  } else if (!stricmp(option, "rthresh")) {
    mp.ratio_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using compression ratio threshold = %2.3f...\n",
           mp.ratio_thresh) ;
  } else if (!stricmp(option, "dt")) {
    mp.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("using dt = %2.3f\n", mp.dt) ;
  } else if (!stricmp(option, "skip")) {
    skip = atoi(argv[2]);
    printf("skipping %d voxels in hires data...\n", skip) ;
    nargs = 1 ;
  } else switch (*option) {
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
        ErrorExit
        (ERROR_NOFILE,
         "%s: could not read transform from %s\n",Progname,argv[2]);
      nargs = 1 ;
      break ;
    case 'I':
      hires_intensity_fname = argv[2] ;
      nargs = 1 ;
      printf("reading intensity image from %s for debugging...\n",
             hires_intensity_fname) ;
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
usage_exit(int ecode) {
  printf(
    "usage: %s <hires labeling> <input intensity> <input aseg>"
    " <output xform>\n",
    Progname) ;
  exit(ecode) ;
}


static TRANSFORM *
compute_optimal_transform(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires,
                          INTEGRATION_PARMS *parms, TRANSFORM *transform) {
  MATRIX    *m_vox_xform, *m_origin, *m_inv_origin, *m_hires_vox2ras,
  *m_lowres_ras2vox, *m_trans, *m_tmp ;
  VECTOR    *v_cl, *v_ch ;
  double    hires_cent[3], lowres_cent[3], dx, dy, dz, scale,min_search_scale ;
  int       niter, nscales, good_step, done, trans ;
  double    old_max_overlap, max_overlap ;
  MRI       *mri_lowres, *mri_hires ;

  mri_lowres = vl_lowres->mri ;
  mri_hires = vl_hires->mri ;

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
  if (transform == NULL) {
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;

    m_lowres_ras2vox = MRIgetRasToVoxelXform(mri_lowres) ;
    m_hires_vox2ras = MRIgetVoxelToRasXform(mri_hires) ;
    MatrixMultiply(m_lowres_ras2vox, m_hires_vox2ras, m_vox_xform) ;
    printf("initial transform from direction cosines:\n") ;
    MatrixPrint(stdout, m_vox_xform) ;
    MatrixFree(&m_lowres_ras2vox) ;
    MatrixFree(&m_hires_vox2ras) ;

    dx = lowres_cent[0] - hires_cent[0] ;
    dy = lowres_cent[1] - hires_cent[1]  ;
    dz = lowres_cent[2] - hires_cent[2] ;

    v_cl = VectorAlloc(4, MATRIX_REAL) ;
    v_ch = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v_cl,4,1) = 1.0 ;
    *MATRIX_RELT(v_ch,4,1) = 1.0 ;
    V3_X(v_ch) = hires_cent[0] ;
    V3_Y(v_ch) = hires_cent[1] ;
    V3_Z(v_ch) = hires_cent[2] ;
    MatrixMultiply(m_vox_xform, v_ch, v_cl) ;
    dx = V3_X(v_cl) - lowres_cent[0] ;
    dy = V3_Y(v_cl) - lowres_cent[1] ;
    dz = V3_Z(v_cl) - lowres_cent[2] ;
    m_trans = MatrixIdentity(4, NULL) ;
    *MATRIX_RELT(m_trans, 1, 4) = -dx ;
    *MATRIX_RELT(m_trans, 2, 4) = -dy ;
    *MATRIX_RELT(m_trans, 3, 4) = -dz ;

    m_tmp = MatrixCopy(m_vox_xform, NULL) ;
    MatrixMultiply(m_trans, m_tmp, m_vox_xform) ;
    printf("after aligning centroids:\n") ;
    MatrixPrint(stdout, m_vox_xform) ;
    max_overlap = (*pf_overlap)(vl_lowres, vl_hires, m_vox_xform) ;
    printf("initial overlap = %2.4f...\n", max_overlap) ;

    MatrixFree(&m_trans) ;
    MatrixFree(&m_tmp) ;
    VectorFree(&v_cl) ;
    VectorFree(&v_ch) ;
    if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
      write_snapshot(mri_lowres,
                     mri_hires,
                     m_vox_xform,
                     parms,
                     parms->start_t,1,NULL);
    }
    parms->start_t++ ;
  } else
    m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;


  trans = MAX(MAX_TRANS,
              MAX(MAX(mri_hires->width,mri_hires->height),
                  mri_hires->depth)/8) ;
  max_overlap = find_optimal_translation(vl_lowres, vl_hires, m_vox_xform,
                                         -trans, trans, 5, 4) ;

  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
    write_snapshot(mri_lowres, mri_hires,
                   m_vox_xform, parms, parms->start_t,1,NULL);
  }
  parms->start_t++ ;
#define MIN_SCALES 3
  /////////////////// loop here ////////////////////////////////////////////
  niter = 0 ;
  nscales = 1 ;
  scale = 1.0 ;
  good_step = 0 ;
  done = 0 ;
  do {
    old_max_overlap = max_overlap ;
    printf("****************************************\n");
    printf("Nine parameter search.  iteration %d nscales = %d ...\n",
           niter+1, nscales);
    printf("****************************************\n");
    max_overlap = find_optimal_linear_xform(vl_lowres, vl_hires,
                                            m_vox_xform, m_origin,
                                            -MAX_ANGLE*scale,
                                            MAX_ANGLE*scale,
                                            1-MAX_SCALE*scale,
                                            1+MAX_SCALE*scale,
                                            -scale*MAX_TRANS,
                                            scale*MAX_TRANS,
                                            3, 3, 3, 2);

    if (parms->write_iterations != 0) {
      write_snapshot(mri_lowres, mri_hires,
                     ((LTA *)(transform->xform))->xforms[0].m_L,
                     parms, parms->start_t+niter, 1, NULL) ;

    }
    printf("Result so far: scale %2.3f: "
           "max overlap = %2.4f, old max overlap=%2.4f\n",
           scale,max_overlap, old_max_overlap) ;
    MatrixPrint(stderr, m_vox_xform);
    /* search a finer nbhd (if do-while continues) */
    if ((max_overlap <= old_max_overlap)) /* couldn't take a step */
    {
      scale *= 0.25 ;
      if (scale < min_search_scale)
        break ;
      good_step = 0 ;
      printf("reducing scale to %2.4f\n", scale) ;
      nscales++ ;
      done = (good_step == 0) ;
    } else
      good_step = 1 ; /* took at least one good step at this scale */

    niter++ ;
  } while (nscales < MIN_SCALES || (done == FALSE)) ;

  /*  m_vox_xform = compute_pca(mri_hires, mri_lowres) ;*/
  parms->start_t += niter ;
  MatrixFree(&m_origin) ;
  MatrixFree(&m_inv_origin) ;
  return(transform) ;
}


/* compute intersection of transformed hires with lowres divided by union */
static double
compute_overlap(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, MATRIX *m_L) {
  int     intersection, Union, x, y, z,
  width, height, depth, xd, yd, zd, label,
  hwidth, hheight, hdepth, i ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_lowres, *mri_hires ;
  static MRI *mri_intersection = NULL;

  mri_lowres = vl_lowres->mri ;
  mri_hires = vl_hires->mri ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  *MATRIX_RELT(v2, 4, 1) = 1.0 ;

  width = mri_lowres->width ;
  height = mri_lowres->height;
  depth = mri_lowres->depth;
  hwidth = mri_hires->width ;
  hheight = mri_hires->height ;
  hdepth = mri_hires->depth;

  mri_intersection = VLSTtoMri(vl_lowres, mri_intersection) ;

  /* first go through lowres volume and for every voxel that is on in it,
     map it to the hires, and if the hires is on, add one to the overlap
  */

  /* go  through lowres volume and for every voxel that is on in it,
     map it to the hires, and if the hires hasn't been counted yet, count it.
  */
#define IN_INTERSECTION 255
#define IN_UNION        254

  for (intersection = Union = i = 0 ; i < vl_hires->nvox ; i++) {
    x = vl_hires->xi[i] ;
    y = vl_hires->yi[i] ;
    z = vl_hires->zi[i] ;

    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    MatrixMultiply(m_L, v1, v2) ;
    xd = nint(V3_X(v2)) ;
    yd = nint(V3_Y(v2)) ;
    zd = nint(V3_Z(v2)) ;
    if (xd >= 0 && xd < width &&
        yd >= 0 && yd < height &&
        zd >= 0 && zd < depth) {
      label = MRIgetVoxVal(mri_intersection, xd, yd, zd, 0) ;
      switch (label) {
      case 0:  /* hires mapping outside of lowres label*/
        Union++ ;
        MRIsetVoxVal(mri_intersection, xd, yd, zd,0, IN_UNION) ;
        break ;
      case IN_UNION:             /* already processed
                                                                          one way or the other */
      case IN_INTERSECTION:
        break ;
      default:                   /* hires mapping into lowres label */
        intersection++ ;
        Union++ ;
        MRIsetVoxVal(mri_intersection, xd, yd, zd, 0, IN_INTERSECTION) ;
        break ;
      }
    } else
      Union++ ; /* penalize for mapping out of FOV */
  }

  /* now count lowres voxels that weren't mapped to in union */
  for (i = 0 ; i < vl_lowres->nvox ; i++) {
    x = vl_lowres->xi[i] ;
    y = vl_lowres->yi[i] ;
    z = vl_lowres->zi[i] ;
    label = MRIgetVoxVal(mri_intersection, x, y, z, 0) ;
    switch (label) {
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
  for (i = 0 ; i < vl_hires->nvox ; i++) {
    x = vl_hires->xi[i] ;
    y = vl_hires->yi[i] ;
    z = vl_hires->zi[i] ;

    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    MatrixMultiply(m_L, v1, v2) ;
    xd = nint(V3_X(v2)) ;
    yd = nint(V3_Y(v2)) ;
    zd = nint(V3_Z(v2)) ;
    if (xd >= 0 && xd < width &&
        yd >= 0 && yd < height &&
        zd >= 0 && zd < depth) {
      MRIsetVoxVal(mri_intersection, xd, yd, zd, 0, 0) ;
    }
  }

  VectorFree(&v1) ;
  VectorFree(&v2) ;
  return((double)intersection/(double)Union) ;
}

// compute mean squared error between distance transforms at voxel locations
static double
compute_distance_transform_sse(VOXEL_LIST *vl_lowres,
                               VOXEL_LIST *vl_hires,
                               MATRIX *m_L) {
  int     x, y, z, width, height, depth,
  hwidth, hheight, hdepth, i ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_lowres, *mri_hires ;
  double  sse, error ;
  double  d1, d2, xd, yd, zd ;
  MATRIX  *m_L_inv ;

  m_L_inv = MatrixInverse(m_L, NULL) ;
  if (m_L_inv == NULL)
    ErrorExit(ERROR_BADPARM,
              "compute_distance_transform_sse: singular matrix.") ;

  mri_lowres = vl_lowres->mri2 ;
  mri_hires = vl_hires->mri2 ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  *MATRIX_RELT(v2, 4, 1) = 1.0 ;

  width = mri_lowres->width ;
  height = mri_lowres->height;
  depth = mri_lowres->depth;
  hwidth = mri_hires->width ;
  hheight = mri_hires->height ;
  hdepth = mri_hires->depth;

  /* go through both voxel lists and compute the sse
     map it to the hires, and if the hires hasn't been counted yet, count it.
  */

  sse = 0.0 ;
  for (i = 0 ; i < vl_hires->nvox ; i++) {
    x = vl_hires->xi[i] ;
    y = vl_hires->yi[i] ;
    z = vl_hires->zi[i] ;

    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    MatrixMultiply(m_L, v1, v2) ;
    d1 = MRIgetVoxVal(vl_hires->mri2, x, y, z, 0) ;
    xd = V3_X(v2) ;
    yd = V3_Y(v2) ;
    zd = V3_Z(v2) ;
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
    MRIsampleVolume(vl_lowres->mri2, xd, yd, zd, &d2) ;
    error = d1-d2 ;
    sse += error*error ;
  }

  /* now count lowres voxels that weren't mapped to in union */
  for (i = 0 ; i < vl_lowres->nvox ; i++) {
    x = vl_lowres->xi[i] ;
    y = vl_lowres->yi[i] ;
    z = vl_lowres->zi[i] ;
    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    MatrixMultiply(m_L_inv, v1, v2) ;
    d1 = MRIgetVoxVal(vl_lowres->mri2, x, y, z, 0) ;

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
    MRIsampleVolume(vl_hires->mri2, xd, yd, zd, &d2) ;
    error = d1-d2 ;
    sse += error*error ;
  }

  VectorFree(&v1) ;
  VectorFree(&v2) ;
  MatrixFree(&m_L_inv) ;
  return(-sqrt(sse / (double)(vl_lowres->nvox + vl_hires->nvox))) ;
}
static double
find_optimal_translation(VOXEL_LIST *vl_lowres,  VOXEL_LIST *vl_hires,
                         MATRIX *m_L, float min_trans, float max_trans,
                         float trans_steps, int nreductions) {
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
  max_overlap = (*pf_overlap)(vl_lowres, vl_hires, m_L) ;

  for (i = 0 ; i <= nreductions ; i++) {
    delta = (max_trans-min_trans) / trans_steps ;
    if (FZERO(delta))
      return(max_overlap) ;
    if (Gdiag & DIAG_SHOW) {
      printf(
        "scanning translations %2.2f->%2.2f (step %2.1f) ",
        min_trans,max_trans, delta) ;
      fflush(stdout) ;
    }
    for (x_trans = min_trans ; x_trans <= max_trans ; x_trans += delta) {
      *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
      for (y_trans = min_trans ; y_trans <= max_trans ; y_trans += delta) {
        *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
        for (z_trans= min_trans ;
             z_trans <= max_trans ;
             z_trans += delta) {
          *MATRIX_RELT(m_trans, 3, 4) = z_trans ;
          if (nint((x_trans)) == -9 && nint((y_trans)) == -5 &&
              nint((z_trans)) == -7)
            DiagBreak() ;

          // get the transform
          m_L_tmp = MatrixMultiply(m_trans, m_L, m_L_tmp) ;
          // calculate the overlap
          overlap = (*pf_overlap)(vl_lowres, vl_hires, m_L_tmp) ;
          if (overlap > max_overlap) {
            max_overlap = overlap ;
            x_max = x_trans ;
            y_max = y_trans ;
            z_max = z_trans ;
#if 1
            printf("new max overlap %2.4f "
                   "found at (%2.1f, %2.1f, %2.1f)\n",
                   max_overlap, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }

    }

    if (Gdiag & DIAG_SHOW)
      printf(
        "max overlap = %2.4f @ (%2.1f, %2.1f, %2.1f)\n",
        max_overlap, x_max, y_max, z_max) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max ;
    // create a new transform by multiplying the previous one.
    MatrixMultiply(m_trans, m_L, m_L_tmp) ;
    MatrixCopy(m_L_tmp, m_L) ;
    x_max = y_max = z_max = 0.0 ;
    /* we've translated transform by old maxs */

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
                          int nreductions) {
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
  x_max_trans = y_max_trans = z_max_trans =
                                x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_overlap = (*pf_overlap)(vl_lowres, vl_hires, m_L) ;
  for (i = 0 ; i < nreductions ; i++) {
    delta_trans = (max_trans-min_trans) / (trans_steps-1) ;
    delta_scale = (max_scale-min_scale) / (scale_steps-1) ;
    delta_rot = (max_angle-min_angle) / (angle_steps-1) ;
    if (Gdiag & DIAG_SHOW) {
      printf("  scanning %2.2f degree nbhd (%2.1f)\n"
             "  scale %2.3f->%2.3f (step %2.3f), "
             "trans %2.2f->%2.2f (step %2.2f)\n",
             (float)DEGREES(max_angle), (float)DEGREES(delta_rot),
             min_scale,max_scale, delta_scale,
             min_trans, max_trans, delta_trans);
      fflush(stdout) ;
    }

    // scale /////////////////////////////////////////////////////////////
    for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta_scale) {
      /*      printf("x_scale = %2.3f\n", x_scale) ;*/
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ;
           y_scale <= max_scale ;
           y_scale += delta_scale) {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ;
             z_scale <= max_scale;
             z_scale += delta_scale) {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) =
            *MATRIX_RELT(m_scale, 2, 4) =
              *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

          // angle /////////////////////////////////
          for (x_angle = min_angle ;
               x_angle <= max_angle ;
               x_angle += delta_rot) {
            m_x_rot = MatrixReallocRotation
                      (4, x_angle, X_ROTATION, m_x_rot) ;
            for (y_angle = min_angle ;
                 y_angle <= max_angle ;
                 y_angle += delta_rot) {
              m_y_rot = MatrixReallocRotation
                        (4, y_angle, Y_ROTATION, m_y_rot);
              m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
              for (z_angle= min_angle;
                   z_angle <= max_angle;
                   z_angle += delta_rot) {
                m_z_rot = MatrixReallocRotation
                          (4, z_angle,Z_ROTATION,m_z_rot);
                m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
                m_tmp2 = MatrixMultiply
                         (m_rot, m_origin_inv, m_tmp2) ;
                MatrixMultiply(m_origin, m_tmp2, m_rot) ;

                m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
                m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3) ;

                // translation //////////////////////
                for (x_trans = min_trans ;
                     x_trans <= max_trans ;
                     x_trans += delta_trans) {
                  *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
                  for (y_trans = min_trans ;
                       y_trans <= max_trans ;
                       y_trans += delta_trans) {
                    *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
                    for (z_trans= min_trans ;
                         z_trans <= max_trans ;
                         z_trans += delta_trans) {
                      *MATRIX_RELT(m_trans, 3, 4) =
                        z_trans ;

                      m_L_tmp = MatrixMultiply
                                (m_trans, m_tmp3, m_L_tmp) ;
                      overlap =
                        (*pf_overlap)(vl_lowres,
                                      vl_hires,
                                      m_L_tmp) ;

                      if (overlap > max_overlap) {
                        max_overlap = overlap ;
                        x_max_scale = x_scale ;
                        y_max_scale = y_scale ;
                        z_max_scale = z_scale ;

                        x_max_rot = x_angle ;
                        y_max_rot = y_angle ;
                        z_max_rot = z_angle ;

                        x_max_trans = x_trans ;
                        y_max_trans = y_trans ;
                        z_max_trans = z_trans ;
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

    if (Gdiag & DIAG_SHOW) {
      printf("  max overlap = %2.4f @ R=(%2.3f,%2.3f,%2.3f),"
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


    x_max_trans = y_max_trans = z_max_trans = 0.0 ;  /* we've translated
                                                                          transform by
                                                                          old maxs */
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
  MatrixFree(&m_x_rot) ;
  MatrixFree(&m_y_rot) ;
  MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ;
  MatrixFree(&m_trans) ;
  MatrixFree(&m_tmp3) ;
  return(max_overlap) ;
}
#define NPARMS (4*4)
#ifdef TOL
#undef TOL
#endif
#define TOL 1e-12

static VOXEL_LIST *Gvl_lowres, *Gvl_hires ;

static int
powell_minimize(VOXEL_LIST *vl_lowres, VOXEL_LIST *vl_hires, MATRIX *mat) {
  float *p, **xi, fret, fstart;
  int   i, r, c, iter ;

  p = vector(1, NPARMS) ;
  xi = matrix(1, NPARMS, 1, NPARMS) ;
  for (i = r = 1 ; r <= 4 ; r++) {
    for (c = 1 ; c <= 4 ; c++) {
      p[i++] = *MATRIX_RELT(mat, r, c) ;
    }
  }

  Gvl_lowres = vl_lowres ;
  Gvl_hires = vl_hires ;
  for (r = 1 ; r <= NPARMS ; r++) {
    for (c = 1 ; c <= NPARMS ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
  do {
    for (r = 1 ; r <= NPARMS ; r++) {
      for (c = 1 ; c <= NPARMS ; c++) {
        xi[r][c] = r == c ? 1 : 0 ;
      }
    }

    fstart = fret ;
    OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
    for (i = r = 1 ; r <= 4 ; r++) {
      for (c = 1 ; c <= 4 ; c++) {
        *MATRIX_RELT(mat, r, c) = p[i++] ;
      }
    }
    *MATRIX_RELT(mat, 4, 1) = 0.0 ;
    *MATRIX_RELT(mat, 4, 2) = 0.0 ;
    *MATRIX_RELT(mat, 4, 3) = 0.0 ;
    *MATRIX_RELT(mat, 4, 4) = 1.0 ;
    printf("%3.3d: best alignment at after powell: %2.3f (%d steps)\n",
           parms.start_t,fret, iter) ;
    write_snapshot(vl_lowres->mri,
                   vl_hires->mri,
                   mat, &parms, parms.start_t++,1,NULL);
  } while (fret < fstart) ;

  free_matrix(xi, 1, NPARMS, 1, NPARMS) ;
  free_vector(p, 1, NPARMS) ;
  return(NO_ERROR) ;
}

static float
compute_powell_sse(float *p) {
  static MATRIX *mat = NULL ;
  float  error ;
  int    i, r, c ;

  if (mat == NULL)
    mat = MatrixAlloc(4, 4, MATRIX_REAL) ;
  for (i = r = 1 ; r <= 4 ; r++) {
    for (c = 1 ; c <= 4 ; c++) {
      *MATRIX_RELT(mat, r, c) = p[i++] ;
    }
  }
  *MATRIX_RELT(mat, 4, 1) = 0.0 ;
  *MATRIX_RELT(mat, 4, 2) = 0.0 ;
  *MATRIX_RELT(mat, 4, 3) = 0.0 ;
  *MATRIX_RELT(mat, 4, 4) = 1.0 ;
  error = -(*pf_overlap)(Gvl_lowres, Gvl_hires, mat) ;
  return(error) ;
}

static int
write_snapshot(MRI *mri_lowres, MRI *mri_hires, MATRIX *m_vox_xform,
               INTEGRATION_PARMS *parms, int fno, int conform, const char *in_fname) {
  MRI *mri_aligned ;
  char fname[STRLEN] ;
  LTA  *lta ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("hires->lowres vox->vox transform:\n") ;
    MatrixPrint(stdout, m_vox_xform) ;
  }
  if (conform) {
    mri_aligned = MRIclone(mri_lowres, NULL) ;
    MRIlinearTransformInterp
    (mri_hires, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
  } else {
#if 0
    mri_aligned = MRIsrcTransformedCentered
                  (mri_hires_intensity, mri_lowres, m_vox_xform, SAMPLE_TRILINEAR) ;
#else
    lta = LTAalloc(1, NULL) ;
    MatrixCopy(m_vox_xform, lta->xforms[0].m_L) ;
    mri_aligned = MRITransformedCenteredMatrix
                  (mri_hires_intensity, mri_lowres, m_vox_xform) ;
#endif
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
    mri_aligned = MRIsrcTransformedCentered
                  (mri_hires, mri_lowres, m_vox_xform, SAMPLE_NEAREST) ;
#else
    mri_aligned = MRITransformedCenteredMatrix
                  (mri_hires, mri_lowres, m_vox_xform) ;
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

  if (mri_hires_intensity) {
#if 0
    mri_aligned = MRIsrcTransformedCentered
                  (mri_hires_intensity, mri_lowres, m_vox_xform, SAMPLE_TRILINEAR) ;
#else
    mri_aligned = MRITransformedCenteredMatrix
                  (mri_hires_intensity, mri_lowres, m_vox_xform) ;
#endif
    if (in_fname) {
      int req = snprintf(fname, STRLEN, "intensity_%s_%s", parms->base_name, in_fname) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "intensity_%s_%03d", parms->base_name, fno) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    if (in_fname) {
      int req = snprintf(fname, STRLEN, "intensity_%s_%s.mgz", parms->base_name, in_fname) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "intensity_%s_%03d.mgz", parms->base_name, fno) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri_aligned, fname) ;
    MRIfree(&mri_aligned) ;
  }

  return(NO_ERROR) ;
}

static MRI *
estimate_densities(GCA_MORPH *gcam, MRI *mri_lowres, MRI *mri_intensities) {
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

  return(GCAMinitDensities(gcam, mri_lowres, mri_intensities, &gcam_ltt)) ;
}
static int
find_gcam_node(GCA_MORPH *gcam, int label, float x0, float y0, float z0) {
  int   x, y, z, nx, ny, nz ;
  double dist, min_dist, dx, dy, dz ;
  GCA_MORPH_NODE *gcamn ;

  min_dist = 1e10 ;

  nx = ny = nz = -1 ;
  for (x = 0 ; x < gcam->width ; x++) {
    for (y = 0 ; y < gcam->height ; y++) {
      for (z = 0 ; z < gcam->depth ; z++) {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label != label)
          continue ;
        dx = gcamn->x - x0 ;
        dy = gcamn->y - y0 ;
        dz = gcamn->z - z0 ;
        dist = sqrt(dx*dx + dy*dy + dz*dz);
        if (dist < min_dist) {
          min_dist = dist ;
          nx = x ;
          ny = y ;
          nz = z ;
        }
      }
    }
  }

  if (nx >= 0) {
    printf("node found at (%d, %d, %d), setting G[xyz]\n", nx, ny, nz);
    Gx = nx ;
    Gy = ny ;
    Gz = nz ;
  }
  return(NO_ERROR) ;
}




