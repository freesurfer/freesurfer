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
// mri_linear_align_binary.c
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
static double MAX_SCALE = 0.25 ;

static int hires_hippo = 0 ;
#define PAD       2

static int ncloses = 0 ;
static float binarize = 0 ;

static int conform = 0 ;

static int binary_label = 128 ;

static int nopowell = 0 ;
static int nfilter = 0 ;
static int apply_transform = 1 ;
static int use_target_label = 1 ;
static int surf_flag = 0 ;

static int write_snapshot(MRI *mri_target, MRI *mri_source,
                          MATRIX *m_vox_xform, MORPH_PARMS *parms,
                          int fno, int conform, char *fname, int nfilter) ;

static double MAX_TRANS = 30 ;

static float compute_powell_rigid_sse(float *p) ;
static float compute_powell_sse(float *p) ;
static int    powell_minimize(VOXEL_LIST *vl_target,
                              VOXEL_LIST *vl_source,
                              MATRIX *mat,
                              MRI *mri_orig_source) ;
static int powell_minimize_rigid(VOXEL_LIST *vl_target,
                                 VOXEL_LIST *vl_source,
                                 MATRIX *mat);


static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;

static TRANSFORM *compute_optimal_transform(VOXEL_LIST *vl_target,
    VOXEL_LIST *vl_source,
    MORPH_PARMS *parms,
    TRANSFORM *transform) ;

static double compute_overlap(VOXEL_LIST *vl_target,
                              VOXEL_LIST *vl_source,
                              MATRIX *m_L) ;
static double compute_distance_transform_sse(VOXEL_LIST *vl_target,
    VOXEL_LIST *vl_source,
    MATRIX *m_L) ;

static double (*pf_overlap)(VOXEL_LIST *vl_target,
                            VOXEL_LIST *vl_source,
                            MATRIX *m_L) = compute_overlap ;

static double find_optimal_translation(VOXEL_LIST *vl_target,
                                       VOXEL_LIST *vl_source,
                                       MATRIX *m_L, float min_trans,
                                       float max_trans,
                                       float trans_steps, int nreductions);

static double find_optimal_linear_xform(VOXEL_LIST *vl_target,
                                        VOXEL_LIST *vl_source,
                                        MATRIX *m_L,
                                        MATRIX *m_origin,
                                        float min_angle, float max_angle,
                                        float min_scale, float max_scale,
                                        float min_trans, float max_trans,
                                        float angle_steps, float scale_steps,
                                        float trans_steps,
                                        int nreductions,
                                        int rigid);
const char *Progname ;
static int target_label = Right_Hippocampus ;

static int wm = 0 ;
static int filled = 0 ;
static int target_skip = 2 ;
static int source_skip = 2 ;
static int npasses = 3 ;

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
static int non_hippo_labels[] = {
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
                                } ;
#define NUM_NON_ARTERY_LABELS  (sizeof(non_artery_labels) / sizeof(non_artery_labels[0]))
#define NUM_NON_HIPPO_LABELS  (sizeof(non_hippo_labels) / sizeof(non_hippo_labels[0]))


static MORPH_PARMS parms ;
static TRANSFORM  *transform = NULL ;

static MRI *mri_orig_source ;  // for debugging
int
main(int argc, char *argv[]) {
  char       **av, *source_fname, *target_fname, *out_fname, fname[STRLEN] ;
  int        ac, nargs, i ;
  MRI        *mri_target = nullptr, *mri_source, *mri_tmp, *mri_dist_src = nullptr,
      *mri_dist_target = nullptr, *mri_whole_source, *mri_whole_target ;
  VOXEL_LIST *vl_target, *vl_source ;
  MRI_REGION  box ;
  Timer start ;
  int          msec, minutes, seconds, label ;
  MATRIX       *m_L ;

  parms.write_iterations = 0 ;
  parms.start_t = 0 ;


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
  mri_whole_source = MRIcopy(mri_source, NULL) ; // to keep track of orig geom

  mri_orig_source = MRIcopy(mri_source, NULL) ;
  if (wm == 1)   // only doing white matter
  {
    MRI *mri_tmp ;
    mri_tmp = MRIclone(mri_source, NULL) ;
    MRIcopyLabel(mri_source, mri_tmp, Left_Cerebral_White_Matter) ;
    MRIcopyLabel(mri_source, mri_tmp, Right_Cerebral_White_Matter) ;
    MRIcopyLabel(mri_source, mri_tmp, Left_Cerebellum_White_Matter) ;
    MRIcopyLabel(mri_source, mri_tmp, Right_Cerebellum_White_Matter) ;
    MRIcopyLabel(mri_source, mri_tmp, Brain_Stem) ;
    MRIfree(&mri_source) ;
    mri_source = mri_tmp ;
    MRIeraseBorders(mri_source, 1) ;
  } else if (binarize >0)
    MRIbinarize(mri_source, mri_source, binarize, 0, binary_label) ;
  else if (!filled) {
    if (use_target_label)  // just copy this label out of volume
    {
      MRI *mri_tmp ;
      mri_tmp = MRIclone(mri_source, NULL) ;
      MRIcopyLabel(mri_source, mri_tmp, target_label) ;
      MRIfree(&mri_source) ;
      mri_source = mri_tmp ;
    } else {
      for (i = 0 ; i < NUM_NON_ARTERY_LABELS ; i++) {
        label = non_artery_labels[i] ;
        MRIreplaceValues(mri_source, mri_source, label, 0) ;
      }
      for (i = 0 ; i < NUM_NON_HIPPO_LABELS ; i++) {
        label = non_hippo_labels[i] ;
        MRIreplaceValues(mri_source, mri_source, label, 0) ;
      }
    }
  }

  for (i = 0 ; i < ncloses ; i++)
    MRIdilate(mri_source, mri_source) ;
  for (i = 0 ; i < ncloses ; i++)
    MRIerode(mri_source, mri_source) ;

  if (FZERO(binarize)/* && !use_target_label*/) {
    MRIboundingBox(mri_source, 0, &box) ;
    box.x -= PAD ;
    box.y -= PAD ;
    box.z -= PAD ;
    box.dx += 2*PAD ;
    box.dy += 2*PAD ;
    box.dz += 2*PAD ;
    MRIcropBoundingBox(mri_source, &box) ;
    mri_tmp = MRIextractRegion(mri_source, NULL, &box) ;
    MRIfree(&mri_source) ;
    mri_source = mri_tmp ;

    mri_tmp = MRIextractRegion(mri_orig_source, NULL, &box) ;
    MRIfree(&mri_orig_source) ;
    mri_orig_source = mri_tmp ;
  }

  if (surf_flag)
  {
    MRI_SURFACE *mris ;
    mris = MRISread(target_fname) ;
    if (!mri_target)
      ErrorExit(ERROR_NOFILE, "%s: could not read target surface %s",
                Progname, target_fname) ;
    mri_target = MRISfillInterior(mris, .25, NULL) ;
    MRIreplaceValues(mri_target, mri_target, 1, target_label) ;
    MRISfree(&mris) ;
  }
  else
  {
    mri_target = MRIread(target_fname) ;
    if (!mri_target)
      ErrorExit(ERROR_NOFILE, "%s: could not read target label volume %s",
                Progname, target_fname) ;
  }
  mri_whole_target = MRIcopy(mri_target, NULL) ; // to keep track of orig geom

  if (wm == 1)   // only doing white matter
  {
    MRI *mri_tmp ;
    mri_tmp = MRIclone(mri_target, NULL) ;
    MRIcopyLabel(mri_target, mri_tmp, Left_Cerebral_White_Matter) ;
    MRIcopyLabel(mri_target, mri_tmp, Right_Cerebral_White_Matter) ;
    MRIcopyLabel(mri_target, mri_tmp, Left_Cerebellum_White_Matter) ;
    MRIcopyLabel(mri_target, mri_tmp, Right_Cerebellum_White_Matter) ;
    MRIcopyLabel(mri_target, mri_tmp, Brain_Stem) ;
    MRIfree(&mri_target) ;
    mri_target = mri_tmp ;
    MRIeraseBorders(mri_target, 1) ;
  } else if (binarize > 0)
    MRIbinarize(mri_target, mri_target, binarize, 0, binary_label) ;
  if (target_label >= 0 && !binarize) {
    MRI *mri_tmp ;
    mri_tmp = MRIclone(mri_target, NULL) ;
    MRIcopyLabel(mri_target, mri_tmp, target_label) ;
    MRIfree(&mri_target) ;
    mri_target = mri_tmp ;
  }

  if (filled) {
    MRI *mri_tmp ;
    mri_tmp = MRIclone(mri_source, NULL) ;
    MRIcopyLabel(mri_source, mri_tmp, target_label) ;
    MRIfree(&mri_source) ;
    mri_source = mri_tmp ;
  }
  for (i = 0 ; i < ncloses ; i++)
    MRIdilate(mri_target, mri_target) ;
  for (i = 0 ; i < ncloses ; i++)
    MRIerode(mri_target, mri_target) ;
  MRIbinarize(mri_target, mri_target, 1, 0, binary_label) ;
  MRIbinarize(mri_source, mri_source, 1, 0, binary_label) ;
  MRIwrite(mri_target, "target_labels.mgz") ;
  MRIwrite(mri_source, "src_labels.mgz") ;
  if (mri_target->type != MRI_UCHAR) {
    mri_tmp = MRIchangeType(mri_target, MRI_UCHAR, 0, 255, 1) ;
    MRIfree(&mri_target) ;
    mri_target = mri_tmp ;
  }
  if (mri_source->type != MRI_UCHAR) {
    mri_tmp = MRIchangeType(mri_source, MRI_UCHAR, 0, 255, 1) ;
    MRIfree(&mri_source) ;
    mri_source = mri_tmp ;
  }

  if (wm || use_target_label)   // debugging
  {
    //    mri_orig_source = mri_source ;
  }
  if (pf_overlap == compute_distance_transform_sse) {
    printf("creating distance transforms...\n") ;
    mri_dist_src = MRIdistanceTransform(mri_source, NULL, binary_label, -1, DTRANS_MODE_SIGNED, NULL);
    mri_dist_target = MRIdistanceTransform(mri_target, NULL, binary_label, -1, DTRANS_MODE_SIGNED, NULL);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      MRIwrite(mri_dist_src, "dist_src.mgz") ;
      MRIwrite(mri_dist_target, "dist_target.mgz") ;
    }
  }

  if (Gdiag & DIAG_WRITE && parms.write_iterations > 0) {
    sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_target, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_target.mgz", parms.base_name) ;
    MRIwrite(mri_target, fname) ;
  }


  /* compute optimal linear transform */
  m_L = NULL ;  // to avoid compiler warning. npasses will never be 0
  for (i = 0 ; i < npasses ; i++) {
    printf("------------- outer loop iteration %d of %d - skip %d, %d ---------------\n",
           i+1, npasses,source_skip, target_skip) ;

    vl_target = VLSTcreate(mri_target,binary_label,binary_label,NULL,target_skip,0);
    // skip=0 for target!
    printf("%d voxels in target list\n", vl_target->nvox) ;
    vl_target->mri2 = mri_dist_target ;

    vl_source = VLSTcreate(mri_source, 1, 255, NULL, source_skip, 0) ;
    printf("%d voxels in source list\n", vl_source->nvox) ;
    vl_source->mri2 = mri_dist_src ;

    transform = compute_optimal_transform(vl_target, vl_source, &parms,transform) ;
    m_L = ((LTA *)(transform->xform))->xforms[0].m_L ;
    VLSTfree(&vl_source) ;
    vl_source = VLSTcreate(mri_source, 1, 255, NULL, source_skip/4, 0) ;
    vl_source->mri2 = mri_dist_src ;
    printf("calling powell minimization...\n") ;
    if ((nopowell == 0) && parms.rigid == 0)
      powell_minimize(vl_target, vl_source, m_L, mri_orig_source) ;
    else if (parms.rigid)
      powell_minimize_rigid(vl_target, vl_source, m_L) ;

    if (apply_transform) {
      MRI *mri_aligned ;
      char   fname[STRLEN] ;

      FileNameRemoveExtension(out_fname, fname) ;
      strcat(fname, ".mgz") ;
      if (binarize > 0) {
        mri_aligned = MRIclone(mri_target, NULL) ;
        MRIlinearTransformInterp(mri_source, mri_aligned, m_L, SAMPLE_NEAREST);
        MRIbinarize(mri_aligned, mri_aligned, 1, 0, 1) ;
      } else {
        mri_aligned = MRITransformedCenteredMatrix(mri_orig_source, mri_target, m_L) ;
      }
      printf("writing transformed output volume to %s...\n", fname) ;
      MRIwrite(mri_aligned, fname) ;
      MRIfree(&mri_aligned) ;

    }
    LTAvoxelToRasXform((LTA *)(transform->xform), mri_source, mri_target) ;
    TransformWrite(transform, out_fname) ;
    LTArasToVoxelXform((LTA *)(transform->xform), mri_source, mri_target) ;
    source_skip /= 2 ;
    target_skip /= 2 ;
    VLSTfree(&vl_target) ;
    VLSTfree(&vl_source) ;
  }

  m_L = ((LTA *)(transform->xform))->xforms[0].m_L ;
  printf("final vox2vox matrix:\n") ;
  MatrixPrint(stdout, m_L) ;
  {
    MRI *mri_aligned, *mri_filtered ;
    char fname[STRLEN] ;
    int  i ;

    mri_aligned = MRITransformedCenteredMatrix(mri_orig_source, mri_target, m_L) ;
    sprintf(fname, "%sfinal.mgz", parms.base_name) ;
    MRIwrite(mri_aligned, fname) ;

    for (i = 1 ; i <= nfilter ; i++) {
      sprintf(fname, "%sfiltered%d.mgz", parms.base_name, i) ;
      mri_filtered = MRImodeFilter(mri_aligned, NULL, i) ;
      printf("writing filtered image to %s\n", fname) ;
      MRIwrite(mri_filtered, fname) ;
      MRIfree(&mri_filtered) ;
    }
    MRIfree(&mri_aligned) ;
  }

  if (apply_transform) {
    MRI *mri_aligned ;
    char   fname[STRLEN] ;

    FileNameRemoveExtension(out_fname, fname) ;
    strcat(fname, ".mgz") ;
    mri_aligned = MRIclone(mri_target, NULL) ;
    MRIlinearTransformInterp(mri_orig_source, mri_aligned, m_L, SAMPLE_NEAREST);
    printf("writing transformed output volume to %s...\n", fname) ;
    MRIwrite(mri_aligned, fname) ;
    MRIfree(&mri_aligned) ;
  }
  LTAvoxelToRasXform((LTA *)(transform->xform), mri_source, mri_target) ;
  LTAsetVolGeom((LTA *)(transform->xform), mri_whole_source, mri_whole_target);
  TransformWrite(transform, out_fname) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
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
  } else if (!stricmp(option, "angio")) {
    pf_overlap = compute_distance_transform_sse ;
    printf("using distance transform to align angiograms\n") ;
  } else if (!stricmp(option, "nopowell")) {
    nopowell = 1 ;
    printf("not applying powell search\n") ;
  } else if (!stricmp(option, "view")) {
    Gsx = atoi(argv[2]) ;
    Gsy = atoi(argv[3]) ;
    Gsz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("viewing voxel (%d, %d, %d)\n", Gsx, Gsy, Gsz) ;
  } else if (!stricmp(option, "wm")) {
    wm = 1 ;
    pf_overlap = compute_distance_transform_sse ;
    target_label = -1 ;
    printf("aligning wm labels...\n") ;
  } else if (!stricmp(option, "target")) {
    target_label = atoi(argv[2]) ;
    pf_overlap = compute_distance_transform_sse ;
    printf("aligning label %d...\n", target_label) ;
    nargs = 1 ;
  } else if (!stricmp(option, "filled")) {
    filled = 1 ;
    target_label = atoi(argv[2]) ;
    pf_overlap = compute_distance_transform_sse ;
    printf("aligning label %d...\n", target_label) ;
    nargs = 1 ;
  } else if (!stricmp(option, "distance")) {
    pf_overlap = compute_distance_transform_sse ;
    printf("using distance transform for SSE\n") ;
  } else if (!stricmp(option, "trans")) {
    MAX_TRANS = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting MAX_TRANS = %2.1f\n", MAX_TRANS) ;
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
  } else if (!stricmp(option, "skip")) {
    source_skip = target_skip = atoi(argv[2]);
    printf("skipping %d voxels in both source and target...\n", source_skip) ;
    nargs = 1 ;
  } else if (!stricmp(option, "source_skip")) {
    source_skip = atoi(argv[2]);
    printf("skipping %d voxels in source...\n", source_skip) ;
    nargs = 1 ;
  } else if (!stricmp(option, "target_skip")) {
    target_skip = atoi(argv[2]);
    printf("skipping %d voxels in target...\n", target_skip) ;
    nargs = 1 ;
  } else switch (*option) {
  case 'H':
    hires_hippo = 1 ;
    //    pf_overlap = compute_distance_transform_sse ;
    use_target_label = 0 ;
    target_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("assuming source is hires hippo volume, and target is %s\n",
           cma_label_to_name(target_label)) ;
    break ;
  case 'S':
    surf_flag = 1 ;
    printf("interpreting target as a surface\n") ;
    break ;
  case 'L':
    target_label = atoi(argv[2]) ;
    printf("using %s %d as target label from source and destination\n", cma_label_to_name(target_label),target_label) ;
    use_target_label = 1 ;
    nargs = 1 ;
    break ;
  case 'B':
    binarize = atof(argv[2]) ;
    printf("binarizing image with thresh = %2.3f\n", binarize) ;
    nargs = 1 ;
    break ;
  case 'C':
    ncloses = atoi(argv[2]) ;
    printf("performing %dth order close on image\n", ncloses) ;
    nargs = 1 ;
    break ;
  case 'N':
    npasses = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using npasses=%d\n", npasses) ;
    break ;
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
  case 'F':
    nfilter = atoi(argv[2]) ;
    nargs = 1 ;
    printf("applying %d mode filters before writing out final volume\n",nfilter);
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
  printf("usage: %s <source> <target> <output xform>\n",Progname) ;
  printf("\t where options are:\n") ;
  printf("\t-h <target label>  -  assume source is hires hippocampal labeling\n");
  exit(ecode) ;
}


static TRANSFORM *
compute_optimal_transform(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source,
                          MORPH_PARMS *parms, TRANSFORM *transform) {
  MATRIX    *m_vox_xform, *m_origin, *m_inv_origin, *m_source_vox2ras,
  *m_target_ras2vox, *m_trans, *m_tmp ;
  VECTOR    *v_cl, *v_ch ;
  double    source_cent[3], target_cent[3],
  dx, dy, dz, scale,min_search_scale ;
  int       niter, nscales, good_step, done, trans ;
  double    old_max_overlap, max_overlap ;
  MRI       *mri_target, *mri_source ;

  done = FALSE;
  mri_target = vl_target->mri ;
  mri_source = vl_source->mri ;

#define MIN_LINEAR_SEARCH_SCALE 0.01
#define MIN_RIGID_SEARCH_SCALE  0.001
  min_search_scale =
    parms->rigid ? MIN_RIGID_SEARCH_SCALE : MIN_LINEAR_SEARCH_SCALE;
  m_origin = MatrixIdentity(4, NULL) ;
  MRIcenterOfMass(mri_source, source_cent, 0) ;
  MRIcenterOfMass(mri_target, target_cent, 0) ;
  *MATRIX_RELT(m_origin, 1, 4) = target_cent[0] ;
  *MATRIX_RELT(m_origin, 2, 4) = target_cent[1] ;
  *MATRIX_RELT(m_origin, 3, 4) = target_cent[2] ;
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;
  m_inv_origin = MatrixInverse(m_origin, NULL) ;
  nscales = 1 ;  // will get incremented if not 1st call
  if (transform == NULL) {
    scale = 1 ;
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;

    m_target_ras2vox = MRIgetRasToVoxelXform(mri_target) ;
    m_source_vox2ras = MRIgetVoxelToRasXform(mri_source) ;
    MatrixMultiply(m_target_ras2vox, m_source_vox2ras, m_vox_xform) ;
    printf("initial transform from direction cosines:\n") ;
    MatrixPrint(stdout, m_vox_xform) ;
    MatrixFree(&m_target_ras2vox) ;
    MatrixFree(&m_source_vox2ras) ;
    old_max_overlap = max_overlap =
                        (*pf_overlap)(vl_target, vl_source, m_vox_xform) ;
    printf("initial overlap = %2.4f...\n", max_overlap) ;

    dx = target_cent[0] - source_cent[0] ;
    dy = target_cent[1] - source_cent[1]  ;
    dz = target_cent[2] - source_cent[2] ;

    v_cl = VectorAlloc(4, MATRIX_REAL) ;
    v_ch = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v_cl,4,1) = 1.0 ;
    *MATRIX_RELT(v_ch,4,1) = 1.0 ;
    V3_X(v_ch) = source_cent[0] ;
    V3_Y(v_ch) = source_cent[1] ;
    V3_Z(v_ch) = source_cent[2] ;
    MatrixMultiply(m_vox_xform, v_ch, v_cl) ;
    dx = V3_X(v_cl) - target_cent[0] ;
    dy = V3_Y(v_cl) - target_cent[1] ;
    dz = V3_Z(v_cl) - target_cent[2] ;
    m_trans = MatrixIdentity(4, NULL) ;
    *MATRIX_RELT(m_trans, 1, 4) = -dx ;
    *MATRIX_RELT(m_trans, 2, 4) = -dy ;
    *MATRIX_RELT(m_trans, 3, 4) = -dz ;

    m_tmp = MatrixCopy(m_vox_xform, NULL) ;
    MatrixMultiply(m_trans, m_tmp, m_vox_xform) ;
    MatrixPrint(stdout, m_vox_xform) ;
    max_overlap = (*pf_overlap)(vl_target, vl_source, m_vox_xform) ;
    printf("after aligning centroids overlap = %2.4f...\n", max_overlap) ;
    if (max_overlap < old_max_overlap) {
      printf("resetting transform to initial one...\n") ;
      MatrixCopy(m_tmp, m_vox_xform) ;
      max_overlap = old_max_overlap ;
    }

    MatrixFree(&m_trans) ;
    MatrixFree(&m_tmp) ;
    VectorFree(&v_cl) ;
    VectorFree(&v_ch) ;
    if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
      write_snapshot
      (mri_target,
       mri_orig_source,
       m_vox_xform, parms, parms->start_t,conform,NULL, nfilter);
    }
    parms->start_t++ ;
  } else {
    scale = .25 ;   // only search fine scales the 2nd and subsequent times
    nscales++ ;
    m_vox_xform = ((LTA *)(transform->xform))->xforms[0].m_L ;
  }


  trans = MAX(MAX_TRANS,
              MAX(MAX(mri_source->width,mri_source->height),
                  mri_source->depth)/8) ;
  max_overlap = find_optimal_translation(vl_target, vl_source, m_vox_xform,
                                         -trans, trans, 9, 6) ;

  MatrixPrint(stdout, m_vox_xform) ;
  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
    write_snapshot(mri_target, mri_orig_source,
                   m_vox_xform, parms, parms->start_t,conform,NULL,nfilter);
  }
  parms->start_t++ ;
#define MIN_SCALES 3
  /////////////////////////// loop here //////////
  niter = 0 ;
  scale *= parms->rigid ? .25: 1.0 ;
  good_step = 0 ;
  do {
    old_max_overlap = max_overlap ;
    printf("****************************************\n");
    printf("%s parameter search.  iteration %d nscales = %d ...\n",
           parms->rigid ? "Six" : "Nine", niter+1, nscales);
    printf("****************************************\n");
    max_overlap = find_optimal_linear_xform(vl_target, vl_source,
                                            m_vox_xform, m_origin,
                                            -MAX_ANGLE*scale,
                                            MAX_ANGLE*scale,
                                            1-MAX_SCALE*scale,
                                            1+MAX_SCALE*scale,
                                            -scale*MAX_TRANS,
                                            scale*MAX_TRANS,
                                            3, 3, 3, 2, parms->rigid);
#if 0
    trans = MAX(MAX_TRANS,
                MAX(MAX(mri_source->width,mri_source->height),
                    mri_source->depth)/(16*scale)) ;
    max_overlap = find_optimal_translation
                  (vl_target, vl_source, m_vox_xform, -trans, trans, 5, 4) ;
#endif

    if (parms->write_iterations != 0) {
      write_snapshot(mri_target, mri_orig_source,
                     ((LTA *)(transform->xform))->xforms[0].m_L,
                     parms, parms->start_t+niter, conform, NULL, nfilter) ;

    }
    printf("Result so far: scale %2.3f: max overlap = %2.4f, "
           "old max overlap=%2.4f\n",
           scale,max_overlap, old_max_overlap) ;
    MatrixPrint(stderr, m_vox_xform);
    /* search a finer nbhd (if do-while continues) */
    if ((max_overlap <= old_max_overlap)) /* couldn't take a step */
    {
      scale *= 0.25 ;
      if (scale < min_search_scale)
        break ;
      printf("reducing scale to %2.4f\n", scale) ;
      nscales++ ;
      done = (good_step == 0) ;
      good_step = 0 ;
    } else
      good_step = 1 ; /* took at least one good step at this scale */

    niter++ ;
  } while (nscales < (MIN_SCALES+2*parms->rigid) || (done == FALSE)) ; // allow extra scales for rigid

  /*  m_vox_xform = compute_pca(mri_source, mri_target) ;*/
  parms->start_t += niter ;
  MatrixFree(&m_origin) ;
  MatrixFree(&m_inv_origin) ;
  return(transform) ;
}


static double best ;

/* compute intersection of transformed source with target divided by union */
#if 1
static double
compute_overlap(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *m_L) {
  int     x, y, z, width, height, depth, hwidth, hheight, hdepth, i ;
  MRI     *mri_target, *mri_source ;
  double  val, source_vox, target_vox, xr, yr, zr, overlap ;
  static VECTOR  *v1 = NULL, *v2 ;
  static MATRIX  *m_L_inv = NULL ;

  m_L_inv = MatrixInverse(m_L,m_L_inv) ;

  mri_target = vl_target->mri ;
  mri_source = vl_source->mri ;

  if (v1 == NULL) {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v1, 4, 1) = 1.0 ;
    *MATRIX_RELT(v2, 4, 1) = 1.0 ;
  }

  width = mri_target->width ;
  height = mri_target->height;
  depth = mri_target->depth;
  hdepth = mri_source->depth;
  hwidth = mri_source->width ;
  hheight = mri_source->height ;

  /* first go through target volume and for every voxel that is on in it,
     map it to the source, and if the source is on, add one to the overlap
  */

  /* go  through source volume and for every voxel that is on in it,
     map it to the target, and if the target hasn't been counted yet, count it.
  */
  for (source_vox = 0.0, i = 0 ; i < vl_source->nvox ; i++) {
    x = vl_source->xi[i] ;
    y = vl_source->yi[i] ;
    z = vl_source->zi[i] ;

    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    if (x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    MatrixMultiply(m_L, v1, v2) ;
    xr = V3_X(v2) ;
    yr = V3_Y(v2) ;
    zr = V3_Z(v2) ;
    MRIsampleVolume(mri_target, xr, yr, zr, &val) ;
    source_vox += (val/binary_label) ; // # of src vox mapped into target
  }

  for (target_vox = 0.0, i = 0 ; i < vl_target->nvox ; i++) {
    x = vl_target->xi[i] ;
    y = vl_target->yi[i] ;
    z = vl_target->zi[i] ;

    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    if (x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    MatrixMultiply(m_L_inv, v1, v2) ;
    xr = V3_X(v2) ;
    yr = V3_Y(v2) ;
    zr = V3_Z(v2) ;
    MRIsampleVolume(mri_source, xr, yr, zr, &val) ;
    target_vox += (val/binary_label) ; // # of target vox mapped into source
  }

  VectorFree(&v1) ;
  VectorFree(&v2) ;
  MatrixFree(&m_L_inv) ;
  overlap = 0.5 * (source_vox / vl_source->nvox + target_vox / vl_target->nvox);
  return(overlap) ;
}
#else
static double
compute_overlap(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *m_L) {
  int     x, y, z, width, height, depth, xd, yd, zd, label,
  hwidth, hheight, hdepth, i ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_target, *mri_source ;
  static MRI *mri_intersection = NULL;
  double  val, intersection, xr, yr, zr, Union, src_weight, target_weight ;

  mri_target = vl_target->mri ;
  mri_source = vl_source->mri ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  *MATRIX_RELT(v2, 4, 1) = 1.0 ;

  width = mri_target->width ;
  height = mri_target->height;
  depth = mri_target->depth;
  hdepth = mri_source->depth;
  hwidth = mri_source->width ;
  hheight = mri_source->height ;
  src_weight = 1.0 / (double)vl_source->nvox ;
  target_weight= 1.0 / (double)vl_target->nvox ;

  // write the target voxel values into the intersection volume
  mri_intersection = VLSTtoMri(vl_target, mri_intersection) ;

  /* first go through target volume and for every voxel that is on in it,
  map it to the source, and if the source is on, add one to the overlap
  */

  /* go  through source volume and for every voxel that is on in it,
  map it to the target, and if the target hasn't been counted yet, count it.
  */
#define IN_INTERSECTION 255
#define IN_UNION        254

  for (intersection = Union = 0.0, i = 0 ; i < vl_source->nvox ; i++) {
    x = vl_source->xi[i] ;
    y = vl_source->yi[i] ;
    z = vl_source->zi[i] ;

    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    if (x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    MatrixMultiply(m_L, v1, v2) ;
    xr = V3_X(v2) ;
    yr = V3_Y(v2) ;
    zr = V3_Z(v2) ;
    xd = nint(xr) ;
    yd = nint(yr) ;
    zd = nint(zr) ;
    if (xd == Gx && yd == Gy && zd == Gz)
      DiagBreak() ;
    if (xd >= 0 && xd < width &&
        yd >= 0 && yd < height &&
        zd >= 0 && zd < depth) {
      label = MRIgetVoxVal(mri_intersection, xd, yd, zd, 0) ;
      switch (label) {
      case 0:  // source mapping outside of target label
        // val = MRIgetVoxVal(mri_target, xd, yd, zd, 0) ;
        MRIsampleVolume(mri_target, xr, yr, zr, &val) ;
        Union += (1-(val/(float)binary_label)) ;
        MRIsetVoxVal(mri_intersection, xd, yd, zd,0, IN_UNION) ;
        break ;
      case IN_UNION:             /* already processed one way
        or the other */
      case IN_INTERSECTION:
        break ;
      default:                   /* not processed yet - target was on
        source mapping into target label */
        // val = MRIgetVoxVal(mri_target, xd, yd, zd, 0) ;
        MRIsampleVolume(mri_target, xr, yr, zr, &val) ;
        intersection += (val/(float)binary_label) ;
        Union++ ;
        MRIsetVoxVal(mri_intersection, xd, yd, zd, 0, IN_INTERSECTION) ;
        break ;
      }
    } else
      Union++ ; /* penalize for mapping out of FOV */
  }

  /* now count target voxels that weren't mapped to in union */
  for (i = 0 ; i < vl_target->nvox ; i++) {
    x = vl_target->xi[i] ;
    y = vl_target->yi[i] ;
    z = vl_target->zi[i] ;
    label = MRIgetVoxVal(mri_intersection, x, y, z, 0) ;
    switch (label) {
    case 0:  /* source mapping outside of target label*/
    case IN_UNION:             /* already processed one way or the other */
    case IN_INTERSECTION:
      break ;
    default:                   // target label that was never mapped
      Union++ ;
      break ;
    }
  }

  /* reset intersection volume */
  for (i = 0 ; i < vl_source->nvox ; i++) {
    x = vl_source->xi[i] ;
    y = vl_source->yi[i] ;
    z = vl_source->zi[i] ;

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
  for (i = 0 ; i < vl_target->nvox ; i++) {
    x = vl_target->xi[i] ;
    y = vl_target->yi[i] ;
    z = vl_target->zi[i] ;
    MRIsetVoxVal(mri_intersection, x, y, z, 0, 0) ;
  }

  VectorFree(&v1) ;
  VectorFree(&v2) ;
  if ((double)intersection/(double)Union > best)
    DiagBreak() ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_intersection, "i.mgz");
  return((double)intersection/(double)Union) ;
}
#endif
// compute mean squared error between distance transforms at voxel locations
static double
compute_distance_transform_sse(VOXEL_LIST *vl_target,
                               VOXEL_LIST *vl_source,
                               MATRIX *m_L) {
  int     x, y, z, width, height, depth,
  hwidth, hheight, hdepth, i ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_target, *mri_source ;
  double  sse, error ;
  double  d1, d2, xd, yd, zd ;
  MATRIX  *m_L_inv ;

  m_L_inv = MatrixInverse(m_L, NULL) ;
  if (m_L_inv == NULL)
    ErrorExit
    (ERROR_BADPARM, "compute_distance_transform_sse: singular matrix.") ;

  mri_target = vl_target->mri2 ;
  mri_source = vl_source->mri2 ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  *MATRIX_RELT(v2, 4, 1) = 1.0 ;

  width = mri_target->width ;
  height = mri_target->height;
  depth = mri_target->depth;
  hwidth = mri_source->width ;
  hheight = mri_source->height ;
  hdepth = mri_source->depth;

  /* go through both voxel lists and compute the sse
     map it to the source, and if the source hasn't been counted yet, count it.
  */

  sse = 0.0 ;
  for (i = 0 ; i < vl_source->nvox ; i++) {
    x = vl_source->xi[i] ;
    y = vl_source->yi[i] ;
    z = vl_source->zi[i] ;

    V3_X(v1) = x ;
    V3_Y(v1) = y ;
    V3_Z(v1) = z ;
    MatrixMultiply(m_L, v1, v2) ;
    d1 = MRIgetVoxVal(vl_source->mri2, x, y, z, 0) ;
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
    MRIsampleVolume(vl_target->mri2, xd, yd, zd, &d2) ;
    error = d1-d2 ;
    sse += error*error ;
  }

#if 1
  /* now count target voxels that weren't mapped to in union */
  for (i = 0 ; i < vl_target->nvox ; i++) {
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
find_optimal_translation(VOXEL_LIST *vl_target,  VOXEL_LIST *vl_source,
                         MATRIX *m_L, float min_trans, float max_trans,
                         float trans_steps, int nreductions) {
  MRI      *mri_target, *mri_source ;
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta,
  overlap, max_overlap, mean_trans ;
  int      i ;

  mri_target = vl_target->mri ;
  mri_source = vl_source->mri ;
  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_overlap = (*pf_overlap)(vl_target, vl_source, m_L) ;

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
          overlap = (*pf_overlap)(vl_target, vl_source, m_L_tmp) ;
          if (overlap > max_overlap) {
            max_overlap = overlap ;
            x_max = x_trans ;
            y_max = y_trans ;
            z_max = z_trans ;
#if 1
            printf("new max overlap %2.4f found at "
                   "(%2.2f, %2.2f, %2.2f)\n",
                   max_overlap, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }

    }

    if (Gdiag & DIAG_SHOW)
      printf(
        "max overlap = %2.4f @ (%2.2f, %2.2f, %2.2f)\n",
        max_overlap, x_max, y_max, z_max) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max ;
    // create a new transform by multiplying the previous one.
    MatrixMultiply(m_trans, m_L, m_L_tmp) ;
    MatrixCopy(m_L_tmp, m_L) ;
    x_max = y_max = z_max = 0.0 ;  /* we've translated
                                                      transform by old maxs */

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
find_optimal_linear_xform(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source,
                          MATRIX *m_L, MATRIX *m_origin,
                          float min_angle, float max_angle,
                          float min_scale, float max_scale,
                          float min_trans, float max_trans,
                          float angle_steps, float scale_steps,
                          float trans_steps,
                          int nreductions,
                          int rigid) {
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
  *m_tmp2, *m_scale, *m_trans, *m_tmp3 = NULL ;
  double   x_angle, y_angle, z_angle, x_max_rot, y_max_rot, z_max_rot,
  delta_rot, x_max_scale, y_max_scale, z_max_scale, delta_scale,
  x_trans, delta_trans, y_trans, z_trans,
  mean_angle, x_scale, y_scale, z_scale, mean_scale, x_max_trans,
  y_max_trans, overlap, max_overlap, z_max_trans, mean_trans ;
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
  x_max_trans = y_max_trans = z_max_trans =
                                x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_overlap = (*pf_overlap)(vl_target, vl_source, m_L) ;
  for (i = 0 ; i < nreductions ; i++) {
    delta_trans = (max_trans-min_trans) / (trans_steps-1) ;
    if (rigid) {
      max_scale = min_scale = 1 ;
      delta_scale = max_scale ;
    } else
      delta_scale = (max_scale-min_scale) / (scale_steps-1) ;
    delta_rot = (max_angle-min_angle) / (angle_steps-1) ;
    if (Gdiag & DIAG_SHOW) {
      if (rigid)
        printf("  scanning %2.4f degree nbhd (%2.4f)\n"
               "trans %2.2f->%2.2f (step %2.2f)\n",
               (float)DEGREES(max_angle), (float)DEGREES(delta_rot),
               min_trans, max_trans, delta_trans);
      else
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

          // angle ///////////////////////////
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

                // translation ////////
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
                      overlap = (*pf_overlap)(vl_target,
                                              vl_source,
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
      if (rigid)
        printf("  max overlap = %2.4f @ R=(%2.5f,%2.5f,%2.5f),"
               "T=(%2.3f,%2.3f,%2.3f)\n",
               max_overlap, DEGREES(x_max_rot), DEGREES(y_max_rot),
               DEGREES(z_max_rot), x_max_trans, y_max_trans,z_max_trans) ;
      else
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

    x_max_trans =
      y_max_trans =
        z_max_trans = 0.0 ;  /* we've translated transform by old maxs */
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
#define NPARMS (4*3)
#ifdef TOL
#undef TOL
#endif
#define TOL 1e-12

static VOXEL_LIST *Gvl_target, *Gvl_source ;

static int
powell_minimize(VOXEL_LIST *vl_target,
                VOXEL_LIST *vl_source,
                MATRIX *mat,
                MRI *mri_orig_source) {
  float *p, **xi, fret, fstart;
  int   i, r, c, iter ;

  p = vector(1, NPARMS) ;
  xi = matrix(1, NPARMS, 1, NPARMS) ;
  for (i = r = 1 ; r <= 3 ; r++)
    for (c = 1 ; c <= 4 ; c++)
      p[i++] = *MATRIX_RELT(mat, r, c) ;

  Gvl_target = vl_target ;
  Gvl_source = vl_source ;
  for (r = 1 ; r <= NPARMS ; r++)
    for (c = 1 ; c <= NPARMS ; c++)
      xi[r][c] = r == c ? 1 : 0 ;

  // TODO:  powell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
  OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
  do {
    for (r = 1 ; r <= NPARMS ; r++)
      for (c = 1 ; c <= NPARMS ; c++)
        xi[r][c] = r == c ? 1 : 0 ;

    fstart = fret ;
    // TODO:    powell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
    OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
    for (i = r = 1 ; r <= 3 ; r++)
      for (c = 1 ; c <= 4 ; c++)
        *MATRIX_RELT(mat, r, c) = p[i++] ;

    printf("%3.3d: best alignment after powell: %2.3f (%d steps)\n",
           parms.start_t,fret, iter) ;
    write_snapshot(vl_target->mri, mri_orig_source,
                   mat, &parms, parms.start_t++,conform,NULL, nfilter);
  } while (fret < fstart) ;

  free_matrix(xi, 1, NPARMS, 1, NPARMS) ;
  free_vector(p, 1, NPARMS) ;
  return(NO_ERROR) ;
}

static float
compute_powell_sse(float *p) {
  static MATRIX *mat = NULL ;
  double  error ;
  int    i, r, c ;

  if (mat == NULL)
    mat = MatrixIdentity(4, NULL) ;
  for (i = r = 1 ; r <= 3 ; r++)
    for (c = 1 ; c <= 4 ; c++)
      *MATRIX_RELT(mat, r, c) = p[i++] ;

  error = -(*pf_overlap)(Gvl_target, Gvl_source, mat) ;
  if (-error > best)
    DiagBreak() ;
  return(error) ;
}

#define NPARMS_RIGID (6)
static int
powell_minimize_rigid(VOXEL_LIST *vl_target, VOXEL_LIST *vl_source, MATRIX *mat) {
  float *p, **xi, fret, fstart ;
  int    r, c, iter, diag ;
  double xr, yr, zr, xt, yt, zt;

  // extract rigid body parameters from matrix
  MatrixToRigidParameters(mat, &xr, &yr, &zr, &xt, &yt, &zt) ;
  printf("initial rigid body parameters = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
         DEGREES(xr), DEGREES(yr), DEGREES(zr), xt, yt, zt) ;
  p = vector(1, NPARMS_RIGID) ;
  xi = matrix(1, NPARMS_RIGID, 1, NPARMS_RIGID) ;
  p[1] = xr ;
  p[2] = yr ;
  p[3] = zr ;
  p[4] = xt ;
  p[5] = yt ;
  p[6] = zt ;

  Gvl_target = vl_target ;
  Gvl_source = vl_source ;
  for (r = 1 ; r <= NPARMS_RIGID ; r++) {
    for (c = 1 ; c <= NPARMS_RIGID ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  diag = Gdiag ;
  Gdiag |= DIAG_VERBOSE ;
  OpenPowell(p, xi, NPARMS_RIGID, TOL, &iter, &fret, compute_powell_rigid_sse);
  MatrixFromRigidParameters(mat, p[1],p[2],p[3],p[4],p[5],p[6]) ;
  printf("%3.3d: best alignment at after powell: "
         "%2.5f (%d steps)\n\tparms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
         parms.start_t,fret, iter,
         DEGREES(p[1]), DEGREES(p[2]), DEGREES(p[3]), p[4], p[5], p[6]) ;
  write_snapshot(vl_target->mri, vl_source->mri,
                 mat, &parms, parms.start_t++,1,NULL, nfilter);
  Gdiag = diag ;
  do {
    for (r = 1 ; r <= NPARMS_RIGID ; r++) {
      for (c = 1 ; c <= NPARMS_RIGID ; c++) {
        xi[r][c] = r == c ? 1 : 0 ;
      }
    }

    fstart = fret ;
    OpenPowell(p, xi, NPARMS_RIGID, TOL, &iter, &fret, compute_powell_rigid_sse);
    MatrixFromRigidParameters(mat, p[1],p[2],p[3],p[4],p[5],p[6]) ;
#if 0
    *MATRIX_RELT(mat, 4, 1) = 0.0 ;
    *MATRIX_RELT(mat, 4, 2) = 0.0 ;
    *MATRIX_RELT(mat, 4, 3) = 0.0 ;
    *MATRIX_RELT(mat, 4, 4) = 1.0 ;
#endif
    printf("%3.3d: best alignment at after powell: "
           "%2.5f (%d steps)\n\tparms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
           parms.start_t,fret, iter,
           DEGREES(p[1]), DEGREES(p[2]), DEGREES(p[3]), p[4], p[5], p[6]) ;
    write_snapshot(vl_target->mri, vl_source->mri,
                   mat, &parms, parms.start_t++,1,NULL, nfilter);
  } while (fret < fstart) ;

  free_matrix(xi, 1, NPARMS_RIGID, 1, NPARMS_RIGID) ;
  free_vector(p, 1, NPARMS_RIGID) ;
  return(NO_ERROR) ;
}

static float
compute_powell_rigid_sse(float *p) {
  static MATRIX *mat = NULL ;
  float  error ;

  if (mat == NULL)
    mat = MatrixAlloc(4, 4, MATRIX_REAL) ;

  MatrixFromRigidParameters(mat, p[1],p[2],p[3],p[4],p[5],p[6]) ;
#if 0
  *MATRIX_RELT(mat, 4, 1) = 0.0 ;
  *MATRIX_RELT(mat, 4, 2) = 0.0 ;
  *MATRIX_RELT(mat, 4, 3) = 0.0 ;
  *MATRIX_RELT(mat, 4, 4) = 1.0 ;
#endif
  error = -(*pf_overlap)(Gvl_target, Gvl_source, mat) ;
  return(error) ;
}

static int
write_snapshot(MRI *mri_target, MRI *mri_source, MATRIX *m_vox_xform,
               MORPH_PARMS *parms, int fno, int conform, char *in_fname, int nfilter) {
  MRI *mri_aligned ;
  char fname[STRLEN] ;
  LTA  *lta ;
  int  i ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("source->target vox->vox transform:\n") ;
    MatrixPrint(stdout, m_vox_xform) ;
  }
  if (1 || conform) {
    mri_aligned = MRIalloc(mri_target->width, mri_target->height,
                           mri_target->depth,mri_source->type);
    MRIcopyHeader(mri_target, mri_aligned) ;
    MRIlinearTransformInterp
    (mri_source, mri_aligned, m_vox_xform, SAMPLE_NEAREST);
  } else {
    lta = LTAalloc(1, NULL) ;
    MatrixCopy(m_vox_xform, lta->xforms[0].m_L) ;
    //    mri_aligned = MRITransformedCenteredMatrix
    // (mri_source_intensity, mri_target, m_vox_xform) ;
    mri_aligned = MRITransformedCenteredMatrix
                  (mri_source, mri_target, m_vox_xform) ;
  }

  for (i = 1 ; i <= nfilter ; i++) 
  {
    MRI *mri_filtered ;

    mri_filtered = MRImodeFilter(mri_aligned, NULL, i) ;
    MRIfree(&mri_aligned) ; mri_aligned = mri_filtered ;
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
    mri_aligned = MRITransformedCenteredMatrix
                  (mri_source, mri_target, m_vox_xform) ;
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
