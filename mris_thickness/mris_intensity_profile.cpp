/**
 * @mris  file_intensity_profile.c
 * @brief program for computing intensity profiles across the cortical ribbon
 *
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "romp_support.h"

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "transform.h"
#include "fmriutils.h"
#include "mrishash.h"
#include "cma.h"
#include "fsinit.h"


int main(int argc, char *argv[]) ;

static int remove_bad_profiles(MRI_SURFACE *mris, MRI *mri_profiles, MRI *mri, float border_mm) ;
static float *mrisComputeWhiteMatterIntensities(MRI_SURFACE *mris, MRI *mri, float wsize_mm, float *norm, int navgs,
                                                MRI *mri_aseg, float mm_border) ;
MRI *MRISfindNearestVerticesAndMeasureCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, int nbhd_size,
    float max_thick, int normalize, int curv_thresh, float *norm);
static MRI *MRIfitPolynomial(MRI *mri_profiles, MRI *mri_poly, int order) ;
static MRI *MRIfitQuadratic(MRI *mri_profiles, MRI *mri_quad) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int pial_normal_avgs = 5 ;
static char *read_laplace_name = NULL ;
static char *write_surf_name = NULL ;
static const char *sphere_name = "sphere" ;
static int fmin_thick = 0 ;
static float laplace_res = 0.5 ;
static int laplace_thick = 0 ;
static char *write_thickness_fname = NULL ;
static INTEGRATION_PARMS parms ;

const char *Progname ;
static char *flat_name = NULL;
static int smooth_iters = 1 ;
static double flat_res = 0 ;

static int nbrs = 1 ;
static char pial_name[100] = "pial" ;
static char white_name[100] = WHITE_MATTER_NAME ;
#define MAX_LABELS 10000
static char *label_names[MAX_LABELS] ;
static int nlabels = 0 ;

static float wm_border_mm = 3.0 ;
static int zero_mean = 0 ;
static MRI *mri_aseg = NULL ;
double normal_sigma = 0 ;
static int inorm = 0 ;
static int nbhd_size = 2 ;
static float max_thick = 5.0 ;
static int remove_bad = 0 ;
static char *sdir = NULL ;
static int num_erode = 0 ;
static float thresh ;
static int normalize = 0 ;
static int use_normal = 0 ;
static int use_pial = 1 ;
static char *curv_fname = NULL ;
static int  curv_thresh = 0 ;
static int navgs = 0 ;
#define MIN_BORDER_DIST 5.0  // mm from border
#define MAX_SAMPLES 20

static char *wm_norm_fname = NULL ;
static int max_samples = MAX_SAMPLES ;
static int quadfit = 0 ;
static int polyfit = 0 ;
static int extra = 0 ;

static float norm_gw = 0 ;
static float norm_white ;
static float norm_csf = 0 ;
static float norm_pial = 0 ;
static float norm_mid = 0 ;
static float norm_mean = 0 ;
static float norm_median = 0 ;

static char *overlay_fname = NULL ;
static int overlay_t0 ;
static int overlay_t1 ;

static int ratio_offsets[4] ;
static int ratio = 0 ;

static int use_normals = 0 ;
static double dist_in, dist_out ;

/* The following specifies the src and dst volumes of the input FSL/LTA transform */
MRI          *lta_src = 0;
MRI          *lta_dst = 0;
static int invert = 0 ;
static char *xform_fname = NULL;
static char *mean_outname = NULL;
static MRI *MRISmeasureCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, int nsamples) ;
static MRI *MRISmeasureNormalCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, double dist_in, double dist_out, int nsamples)  ;

int
main(int argc, char *argv[]) {
  char          *out_fname, *sname, buf[STRLEN], *cp, fname[STRLEN], *hemi ;
  int           nargs ;
  MRI_SURFACE   *mris ;
  LABEL         *label = NULL ;
  MRI           *mri, *mri_profiles ;
  float         *norm = NULL ;


  nargs = handleVersionOption(argc, argv, "mris_intensity_profile");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.dt = 0.01 ;
  parms.remove_neg = 1 ;
  parms.momentum = .1; 
  parms.niterations = 1000 ;
  parms.l_nlarea = 0 ;
  parms.l_thick_min = 0 ;
  parms.l_tsmooth = 0 ;
  parms.l_thick_parallel = 1 ;
  parms.l_thick_spring = 1 ;
  parms.l_ashburner_triangle = 1 ;
  parms.l_ashburner_lambda = .1 ;
  parms.l_tspring = 0;
  parms.l_thick_normal = .01;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.tol = 1e-3 ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit() ;

  sname = argv[1] ;
  hemi = argv[2] ;
  out_fname = argv[4] ;
  if (!sdir) {
    sdir = buf ;
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, pial_name) ;
  fprintf(stderr, "reading pial surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
  if (nbrs > 1)
    MRISresetNeighborhoodSize(mris, nbrs) ;


  if (curv_fname) {
    if (MRISreadCurvatureFile(mris, curv_fname)  != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s", curv_fname) ;
  }

  fprintf(stderr, "reading intensity volume from %s...\n", argv[3]) ;
  mri = MRIread(argv[3]) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s",
              Progname, argv[3]) ;
  if (inorm) {
    int   x, y, z, n ;
    float mean, val ;

    for (mean = 0.0, n = 0, x = mri->width/2-5 ; x <= mri->width/2+5 ; x++)
      for (y = 0 ; y <= 5 ; y++)
        for (z = mri->depth/2-5 ; z <= mri->depth/2+5 ; z++) {
          val = MRIgetVoxVal(mri, x, y, z, 0) ;
          mean += val ;
          n++ ;
        }
    mean /= (float)n ;
    printf("mean noise level %2.1f - scaling by %2.2f (DOF=%d)\n", mean, 10.0/mean, n) ;
    MRIscalarMul(mri, mri, 10.0/mean) ;
  }


  if (wm_norm_fname) {
    MRI *mri_wm ;
    double mean ;

    mri_wm = MRIread(wm_norm_fname) ;
    if (mri_wm == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read wm volume from %s...", wm_norm_fname) ;
    MRIbinarize(mri_wm, mri_wm, MIN_WM_VAL, 0, 1) ;
    mean = MRImeanInLabel(mri, mri_wm, 1) ;
    printf("mean wm = %2.0f, normalizing by %2.2f\n", mean, 100/mean) ;
    MRIscalarMul(mri, mri, 100.0/mean) ;
    MRIfree(&mri_wm) ;
  }
  if (MRISreadOriginalProperties(mris, white_name) != NO_ERROR)
    ErrorExit(Gerror, "%s: could not read white matter surface", Progname) ;

  if (use_normals)
  {
    MRI *mri_profiles ;

    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    mri_profiles = MRISmeasureNormalCorticalIntensityProfiles(mris, mri, dist_in, dist_out, max_samples)  ;
    printf("writing cortical intensity profiles to %s...\n", out_fname) ;
    MRIwrite(mri_profiles, out_fname) ;
    exit(0) ;
  }

  
  fprintf(stderr, "measuring gray matter intensity profile...\n") ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
  if (use_pial)
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  else
    MRIScomputeMetricProperties(mris) ; // use surface normals from white surface

    

  //use the specified xform if it's given
  if (xform_fname != NULL) {
    LTA *lta = 0;
    int transform_type;

    printf("INFO: Applying transformation from file %s...\n",
           xform_fname);
    transform_type =  TransformFileNameType(xform_fname);

    if (transform_type == MNI_TRANSFORM_TYPE ||
        transform_type == TRANSFORM_ARRAY_TYPE ||
        transform_type == REGISTER_DAT ||
        transform_type == FSLREG_TYPE
       ) {

      printf("Reading transform ...\n");
      lta = LTAreadEx(xform_fname) ;

      if (!lta)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
                  Progname, xform_fname) ;

      if (transform_type == FSLREG_TYPE) {
        if (lta_src == 0 || lta_dst == 0) {
          fprintf(stderr, "ERROR: fslmat does not have information on the src and dst volumes\n");
          fprintf(stderr, "ERROR: you must give options '-src' and '-dst' to specify the src and dst volume infos for the registration\n");
        }


        LTAmodifySrcDstGeom(lta, lta_src, lta_dst); // add src and dst information
        //The following is necessary to interpret FSLMAT correctly!!!
        LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      }
      if (lta->xforms[0].src.valid == 0) {
        if (lta_src == 0) {
          fprintf(stderr, "The transform does not have the valid src volume info.\n");
          fprintf(stderr, "Either you give src volume info by option -src or\n");
          fprintf(stderr, "make the transform to have the valid src info.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        } else {
          LTAmodifySrcDstGeom(lta, lta_src, NULL); // add src information
        }
      }
      if (lta->xforms[0].dst.valid == 0) {
        if (lta_dst == 0) {
          fprintf(stderr, "The transform does not have the valid dst volume info.\n");
          fprintf(stderr, "Either you give src volume info by option -dst or\n");
          fprintf(stderr, "make the transform to have the valid dst info.\n");
          fprintf(stderr, "If the dst was average_305, then you can set\n");
          fprintf(stderr, "environmental variable USE_AVERAGE305 true\n");
          fprintf(stderr, "instead.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        } else {
          LTAmodifySrcDstGeom(lta, NULL, lta_dst); // add  dst information
        }
      }

      if (invert) {
        VOL_GEOM vgtmp;
        LT *lt;
        MATRIX *m_tmp = lta->xforms[0].m_L ;
        lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
        MatrixFree(&m_tmp) ;
        lt = &lta->xforms[0];
        if (lt->dst.valid == 0 || lt->src.valid == 0) {
          fprintf(stderr, "WARNING:***************************************************************\n");
          fprintf(stderr, "WARNING:dst volume infor is invalid.  Most likely produce wrong inverse.\n");
          fprintf(stderr, "WARNING:***************************************************************\n");
        }
        copyVolGeom(&lt->dst, &vgtmp);
        copyVolGeom(&lt->src, &lt->dst);
        copyVolGeom(&vgtmp, &lt->src);
      }


      LTAchangeType(lta, LINEAR_RAS_TO_RAS);

      printf("---------------------------------\n");
      printf("INFO: Transform Matrix \n");
      MatrixPrint(stdout,lta->xforms[0].m_L);
      printf("---------------------------------\n");

      // convert surface into hires volume surface if needed
      //////////////////////////////////////////////////////////////////////
      MRISsurf2surfAll(mris, mri, lta);

      if (lta_src)
        MRIfree(&lta_src);
      if (lta_dst)
        MRIfree(&lta_dst);

    } else {
      fprintf(stderr, "unknown transform type in file %s\n",
              xform_fname);
      exit(1);
    }
  }
  // now convert the surface to be in the volume ras coords
  else if (mriConformed(mri) == 0 && mris->useRealRAS == 0) {
#if 0
    LINEAR_TRANSFORM *lt = 0;
    MATRIX           *m_L = 0;
    TRANSFORM        *xform ;
    MRI              *mri_cor ;
    LTA              *lta ;

    printf("volume is not conformed - transforming surfaces....\n") ;
    sprintf(fname, "%s/%s/mri/T1/COR-256", sdir, sname) ;
    if (!FileExists(fname))
      sprintf(fname, "%s/%s/mri/T1.mgz", sdir, sname) ;
    fprintf(stderr, "conformed volume from %s...\n", fname) ;
    mri_cor = MRIread(fname) ;
    if (!mri_cor)
      ErrorExit(ERROR_NOFILE, "%s: could not read conformed volume from %s",
                Progname, fname) ;

    fprintf(stderr, "allocating identity RAS-to-RAS xform...\n") ;
    // allocate xform->xform
    xform = TransformAlloc(MNI_TRANSFORM_TYPE, NULL);
    if (!xform)
      ErrorExit(ERROR_NOFILE, "%s: could not allocate hires xform %s", Progname, argv[3]) ;
    lta = (LTA *) (xform->xform);
    lt = &lta->xforms[0];
    lt->sigma = 1.0f ;
    lt->x0 = lt->y0 = lt->z0 = 0 ;
    lta->type = LINEAR_RAS_TO_RAS;
    m_L = lt->m_L;
    MatrixIdentity(4, m_L);
    getVolGeom(mri, &lt->dst);
    getVolGeom(mri_cor, &lt->src);
    if (mris->useRealRAS == 0)
      lt->src.c_r = lt->src.c_a = lt->src.c_s = 0 ;
    if (Gdiag & DIAG_WRITE) {
      MRI *mri_tmp;
      mri_tmp = MRIclone(mri, NULL);
      MRIapplyRASlinearTransform(mri_cor, mri_tmp, lt->m_L);
      MRIwrite(mri_tmp, "T1_xformed.mgz") ;
      MRIfree(&mri_tmp) ;
    }
    MRIfree(&mri_cor) ;
    ///////////////////////////////////////////////////////////////////////
    // convert surface into hires volume surface if needed
    //////////////////////////////////////////////////////////////////////
    MRISsurf2surfAll(mris, mri, lta);
#endif
  }

  if (laplace_thick)
  {
    MRI *mri_laplace ;
    MRI *mri_dist_white, *mri_dist_pial ;
    int nwmissing, npmissing, nmissing, vno ;
    double xv, yv, zv, white_val, pial_val ;
    VERTEX *v ;
    FILE   *fp ;
    
    if (read_laplace_name)
    {
      mri_laplace = MRIread(read_laplace_name) ;
      if (mri_laplace == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read laplace volume from %s",
		  read_laplace_name) ;
    }
    else
    {
      MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
      mri_dist_pial = MRIScomputeDistanceToSurface(mris, NULL, laplace_res) ;
      
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;  // actually white
      MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
      mri_dist_white = MRIScomputeDistanceToSurface(mris, NULL, laplace_res) ;
      for (nmissing = nwmissing = npmissing = vno = 0 ; vno < mris->nvertices ; vno++)
      {
	v = &mris->vertices[vno] ;
	MRISsurfaceRASToVoxel(mris, mri_dist_pial, v->pialx, v->pialy, v->pialz, &xv, &yv, &zv);
	MRIsampleVolumeFrameType(mri_dist_pial, xv, yv, zv, 0, SAMPLE_NEAREST, &pial_val) ;
	MRISsurfaceRASToVoxel(mris, mri_dist_white, v->whitex, v->whitey, v->whitez, &xv, &yv, &zv);
	MRIsampleVolumeFrameType(mri_dist_white, xv, yv, zv, 0, SAMPLE_NEAREST, &white_val) ;
	if (fabs(white_val) > laplace_res)
	  nwmissing++ ;
	if (fabs(pial_val) > laplace_res)
	  npmissing++ ;
	if ((fabs(pial_val) > laplace_res) || (fabs(white_val) > laplace_res))
	  nmissing++ ;
      }
      printf("%d of %d pial surface nodes not resolved - %2.3f %%\n",
	     npmissing, mris->nvertices, 100.0*npmissing/mris->nvertices) ;
      printf("%d of %d gray/white surface nodes not resolved - %2.3f %%\n",
	     nwmissing, mris->nvertices, 100.0*nwmissing/mris->nvertices) ;
      MRIfree(&mri_dist_white) ; MRIfree(&mri_dist_pial) ;
      MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
      fp = fopen("laplace_missing.txt", "a") ;
      fprintf(fp, "%f %d %d %f %d %f %d %f\n", 
	      laplace_res,  mris->nvertices,
	      nwmissing, 100.0*nwmissing/mris->nvertices,
	      npmissing, 100.0*npmissing/mris->nvertices,
	      nmissing, 100.0*nmissing/mris->nvertices) ;
      fclose(fp) ;
      MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
      mri_laplace = MRISsolveLaplaceEquation(mris, NULL, laplace_res) ;
    }
    mri_profiles = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, max_samples) ;
    MRISmeasureLaplaceStreamlines(mris, mri_laplace, mri, mri_profiles) ;
  }
  else if (fmin_thick)
  {
    char              *cp, surf_fname[STRLEN], fname[STRLEN] ;

    printf("computing variational vector field...\n") ;
    
    if (parms.base_name[0] == 0) {
      
      FileNameOnly(out_fname, fname) ;
      cp = strchr(fname, '.') ;
      if (cp)
        strcpy(parms.base_name, cp+1) ;
      else
        strcpy(parms.base_name, sphere_name) ;
      cp = strrchr(parms.base_name, '.') ;
      if (cp)
        *cp = 0 ;
    }

    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    if (MRISreadVertexPositions(mris, sphere_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, sphere_name) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    MRIScomputeSurfaceNormals(mris, PIAL_VERTICES, pial_normal_avgs) ;

    // read in icosahedron data (highly tessellated one)
    cp = getenv("FREESURFER_HOME");
    if (cp == NULL)
      ErrorExit(ERROR_BADPARM, "%s: FREESURFER_HOME not defined in environment", cp) ;
    sprintf(surf_fname,"%s/lib/bem/ic7.tri",cp);

    MRISripZeroThicknessRegions(mris) ;
    MRISstoreRipFlags(mris) ;
    MRISminimizeThicknessFunctional(mris, &parms, max_thick) ;
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    if (write_surf_name)
    {
      int vno ;
      MRI    *mri_vector ;
      mri_vector = MRIallocSequence(mris->nvertices, 1,1,MRI_FLOAT,3) ;

      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
	VERTEX *v ;
	
	v = &mris->vertices[vno] ;
	if (vno == Gdiag_no)
	  DiagBreak() ;
	
	if (v->ripflag)
	  continue ;
	
	MRIsetVoxVal(mri_vector, vno, 0, 0, 0, v->pialx-v->whitex) ;
	MRIsetVoxVal(mri_vector, vno, 0, 0, 1, v->pialy-v->whitey) ;
	MRIsetVoxVal(mri_vector, vno, 0, 0, 2, v->pialz-v->whitez) ;
      }
      printf("writing vector field to %s\n", write_surf_name) ;
      MRIwrite(mri_vector, write_surf_name) ;
      MRIfree(&mri_vector) ;
    }
    mri_profiles = MRISmeasureCorticalIntensityProfiles(mris, mri, MAX_SAMPLES) ;
  }

  if (nlabels > 0) {
    int l ;
    char label_name[STRLEN] ;
    LABEL *ltotal = NULL ;

    for (l = 0 ; l < nlabels ; l++) {
      sprintf(label_name, "%s/%s/label/%s.%s.label", sdir, sname, hemi,label_names[l]) ;

      label = LabelRead(NULL, label_name) ;
      if (!label)
        ErrorExit(ERROR_NOFILE, "%s: could not read label file %s...\n", Progname,
                  label_name) ;
      if (num_erode > 0) {
        printf("eroding label %d times, npoints went from %d ", num_erode,label->n_points) ;
        LabelErode(label, mris, num_erode) ;
        printf("to %d ", label->n_points) ;
      }
      ltotal = LabelCombine(label, ltotal) ;
    }
    if (nlabels == 0)
      ltotal = LabelInFOV(mris, mri, MIN_BORDER_DIST) ;

    LabelRipRestOfSurfaceWithThreshold(ltotal, mris, thresh) ;
  }

  if (norm_white > 0)
    norm = mrisComputeWhiteMatterIntensities(mris, mri, 15, NULL, 100, mri_aseg, wm_border_mm) ;
  if (flat_name)
  {
    if (MRISreadFlattenedCoordinates(mris, flat_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not load flat patch %d\n", Progname, flat_name) ;
    mri_profiles = MRIScomputeFlattenedVolume(mris, mri, flat_res, max_samples, normalize,NULL, smooth_iters, 0, 0) ;
  }
  else if (!laplace_thick && !fmin_thick)
    mri_profiles =
      MRISfindNearestVerticesAndMeasureCorticalIntensityProfiles(mris, mri, nbhd_size, max_thick, normalize, curv_thresh, norm) ;

  mri_profiles->tr = 1 ;
  if (remove_bad)
    remove_bad_profiles(mris, mri_profiles, mri, 20) ;
  if (norm_csf > 0)
  {
    MRI    *mri_ven, *mri_tmp ;
    double csf_mean ;
    char   fname[STRLEN], tmp[STRLEN] ;
    FILE   *fp ;
    int    nvox, n ;

    if (mri_aseg == NULL)
      ErrorExit(ERROR_BADPARM, "%s: must specify aseg volume with -aseg", Progname) ;

    mri_ven = MRIclone(mri_aseg, NULL) ;
    MRIcopyLabel(mri_aseg, mri_ven, Left_Lateral_Ventricle) ;
    MRIcopyLabel(mri_aseg, mri_ven, Right_Lateral_Ventricle) ;
    mri_tmp = MRIcopy(mri_ven, NULL) ;
    // erode as many times as we can while still having enough for stable estimate of CSF
    n = 0 ;
    do
    {
      MRIcopy(mri_tmp, mri_ven) ;
      MRIerode(mri_ven, mri_tmp) ;
      nvox = MRItotalVoxelsOn(mri_tmp, 1) ;
      n++ ;
      if (nvox < 100)
        break ;
    } while ((nvox > 250) || (n < 2)) ;  
    csf_mean = MRImeanInLabel(mri,  mri_ven, Left_Lateral_Ventricle) ;
    csf_mean += MRImeanInLabel(mri,  mri_ven, Right_Lateral_Ventricle) ;
    csf_mean /= 2 ;
    MRIfree(&mri_tmp) ;

    sprintf(fname, "%s.csf.ven.mgz", FileNameRemoveExtension(out_fname, tmp)) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_ven, fname) ;
    sprintf(fname, "%s.csf.dat", FileNameRemoveExtension(out_fname, tmp)) ;
    printf("normalizing by %d-eroded csf mean %2.0f, storing in %s\n", n, csf_mean, fname) ;
    fp = fopen(fname, "w") ; 
    if (fp)
    {
      fprintf(fp, "%f\n", csf_mean) ; 
      fclose(fp) ;
    }
    else
      printf("*************** couldn't open file %s ********************\n", fname) ;
    MRIfree(&mri_ven) ;
    MRIscalarMul(mri_profiles, mri_profiles, norm_csf/csf_mean) ;
  }
  if (navgs > 0) {
    MRIwriteFrame(mri_profiles, "lh.f0.mgz", 0) ;
    printf("smoothing profiles %d times\n", navgs) ;
    MRISsmoothFrames(mris, mri_profiles, navgs) ;
    MRIwriteFrame(mri_profiles, "lh.f0.smooth.mgz", 0) ;
  }

  if (mean_outname != NULL) {
    MRI *mri_mean;
    int vtx, index;
    double sum_profile;

    printf("Compute the mean of the intensity profile ...\n");

    mri_mean = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 1);

    for (vtx = 0; vtx < mris->nvertices; vtx++) {
      sum_profile = 0;

      if (mris->vertices[vtx].ripflag) {
        MRIsetVoxVal(mri_mean, vtx, 0, 0, 0, sum_profile);
        continue ;
      }

      for (index= 0; index < mri_profiles->nframes; index++)
        sum_profile += MRIFseq_vox(mri_profiles, vtx, 0, 0, index);

      sum_profile /= (float)mri_profiles->nframes;

      MRIsetVoxVal(mri_mean, vtx, 0, 0, 0, sum_profile);
    }

    MRIScopyMRI(mris, mri_mean, 0, "curv");

    printf("writing average intensity profiles to curvature file %s...\n", mean_outname) ;
    MRISwriteCurvature(mris, mean_outname);

    MRIfree(&mri_mean);
  }


  if (normal_sigma > 0)
  {
    MRI *mri_tmp ;

    if (Gdiag_no >= 0)
    {
#if 0
      float v1, v2, v3 ;
      v1 = MRIgetVoxVal(mri_profiles, Gdiag_no, 0, 0, 0) ;
      v2 = MRIgetVoxVal(mri_profiles, Gdiag_no, 0, 0, 1) ;
      v3 = MRIgetVoxVal(mri_profiles, Gdiag_no, 0, 0, 2) ;
#endif
    }
    printf("smoothing profiles in normal direction\n") ;
    mri_tmp = fMRItemporalGaussian(mri_profiles, normal_sigma, NULL);
    MRIfree(&mri_profiles) ; mri_profiles = mri_tmp ; 
    if (Gdiag_no >= 0)
    {
#if 0
      float v1, v2, v3 ;
      v1 = MRIgetVoxVal(mri_profiles, Gdiag_no, 0, 0, 0) ;
      v2 = MRIgetVoxVal(mri_profiles, Gdiag_no, 0, 0, 1) ;
      v3 = MRIgetVoxVal(mri_profiles, Gdiag_no, 0, 0, 2) ;
      DiagBreak() ;
#endif
    }
  }
  if (quadfit) {
    MRI *mri_quad ;

    mri_quad = MRIfitQuadratic(mri_profiles, NULL) ;
    printf("writing best fitting quadratic fits to %s...\n", out_fname) ;
    MRIwrite(mri_quad, out_fname) ;
  } else if (polyfit > 0) {
    MRI *mri_poly ;
    char fname[STRLEN], ext[STRLEN] ;

    mri_poly = MRIfitPolynomial(mri_profiles, NULL, polyfit) ;
    if (zero_mean)
    {
      printf("removing mean from profiles...\n") ;
      MRIzeroMeanTimecourse(mri_profiles, mri_profiles) ;
    }

    printf("writing cortical intensity profiles to %s...\n", out_fname) ;
    MRIwrite(mri_profiles, out_fname) ;
    FileNameExtension(out_fname, ext) ;
    FileNameRemoveExtension(out_fname, fname) ;
    strcat(fname, "_poly.mgz") ;
    printf("writing best fitting %dth order polynomial fits to %s...\n", polyfit, fname) ;
    MRIwrite(mri_poly, fname) ;
  } else {
    if (zero_mean)
      MRIzeroMeanTimecourse(mri_profiles, mri_profiles) ;
    printf("writing cortical intensity profiles to %s...\n", out_fname) ;
    MRIwrite(mri_profiles, out_fname) ;
    if (write_thickness_fname)
      MRISwriteCurvature(mris, write_thickness_fname) ;
  }
  if (overlay_fname)
  {
    MRI *mri_overlay ;
    int vno, t ;
    double avg, avg2 ;

    mri_overlay = MRIalloc(mri_profiles->width, 1, 1, MRI_FLOAT) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (ratio)
      {
        for (avg = 0.0, t = ratio_offsets[0] ; t <= ratio_offsets[1] ; t++)
        {
          avg += MRIgetVoxVal(mri_profiles, vno, 0, 0, t) ;
        }
        avg /= (1.0+ratio_offsets[1]-ratio_offsets[0]);
        for (avg2 = 0.0, t = ratio_offsets[2] ; t <= ratio_offsets[3] ; t++)
        {
          avg2 += MRIgetVoxVal(mri_profiles, vno, 0, 0, t) ;
        }
        avg2 /= (1.0+ratio_offsets[3]-ratio_offsets[2]);
        if (!FZERO(avg2))
          avg /= avg2 ;
      }
      else
      {
        for (avg = 0.0, t = overlay_t0 ; t <= overlay_t1 ; t++)
        {
          avg += MRIgetVoxVal(mri_profiles, vno, 0, 0, t) ;
        }
        avg /= (1.0+overlay_t1-overlay_t0);
      }
      MRIsetVoxVal(mri_overlay, vno, 0, 0, 0, avg);
    }
    fprintf(stderr, "writing overlay to %s...\n", overlay_fname) ;
    MRIwrite(mri_overlay, overlay_fname) ;
    MRIfree(&mri_overlay) ;
  }
  MRIfree(&mri_profiles) ;

  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "wm_border")) {
    wm_border_mm = atof(argv[2]) ;
    fprintf(stderr,  "using %2.1f mm border region in wm normalization\n",
            wm_border_mm) ;
    nargs = 1 ;
  } else if (!stricmp(option, "write_thickness") || 
	     (!stricmp(option, "thickness"))) {
    write_thickness_fname = argv[2] ;
    fprintf(stderr,  "writing thickness file to %s\n", write_thickness_fname);
    nargs = 1 ;
  } else if (!stricmp(option, "write_surf")) {
    write_surf_name = argv[2] ;
    printf("writing variational pial surface locations to %s\n", write_surf_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tol")) {
    parms.tol = atof(argv[2]) ;
    printf("using tol = %f\n", parms.tol) ;
    nargs = 1 ;
  } else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    printf("setting neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  } else if (!stricmp(option, "rl")) {
    read_laplace_name = argv[2] ;
    printf("reading Laplace solution volume from %s\n", read_laplace_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "triangle")) {
    parms.l_ashburner_triangle = atof(argv[2]) ;
    printf("using l_ashburner_triangle = %2.3f\n", parms.l_ashburner_triangle) ;
    nargs = 1 ;
  } else if (!stricmp(option, "new") || !stricmp(option, "fmin") || !stricmp(option, "variational")) {
    fmin_thick = 1 ;
    fprintf(stderr,  "using variational thickness measurement\n") ;
  } else if (!stricmp(option, "flatten")) {
    flat_name = argv[2] ;
    flat_res = atof(argv[3]) ;
    fprintf(stderr,  "reading flattened coordinates from %s, and resolution %2.2f\n", flat_name,
            flat_res) ;
    nargs = 2 ;
  } else if (!stricmp(option, "pial")) {
    strcpy(pial_name, argv[2]) ;
    fprintf(stderr,  "reading pial surface from file named %s\n", pial_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "white")) {
    strcpy(white_name, argv[2]) ;
    fprintf(stderr,  "reading white surface from file named %s\n", white_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "dt")) {
    parms.dt = atof(argv[2]) ;
    fprintf(stderr,  "setting dt=%2.3f\n", parms.dt) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tsmooth")) {
    parms.l_tsmooth = atof(argv[2]) ;
    fprintf(stderr,  "setting l_tsmooth=%2.3f\n", parms.l_tsmooth) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tnormal")) {
    parms.l_thick_normal = atof(argv[2]) ;
    fprintf(stderr,  "setting l_thick_normal=%2.3f\n", parms.l_thick_normal) ;
    nargs = 1 ;
  } else if (!stricmp(option, "normal") || !stricmp(option, "normals")) {
    use_normals = 1 ;
    dist_in = atof(argv[2]) ;
    dist_out = atof(argv[3]) ;
    fprintf(stderr,  "computing profiles along surface normal in interval [-%2.2f, %2.2f] with %d samples (del=%2.3f)\n",
	    -dist_in, dist_out, max_samples, (dist_in+dist_out)/(float)(max_samples-1)) ;
    nargs = 2 ;
  } else if (!stricmp(option, "tspring")) {
    parms.l_thick_spring = atof(argv[2]) ;
    fprintf(stderr,  "setting l_thick_spring=%2.3f\n", parms.l_thick_spring) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tparallel")) {
    parms.l_thick_parallel = atof(argv[2]) ;
    fprintf(stderr,  "setting l_thick_parallel=%2.3f\n", parms.l_thick_parallel) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tmin")) {
    parms.l_thick_min = atof(argv[2]) ;
    fprintf(stderr,  "setting l_thick_min=%2.3f\n", parms.l_thick_min) ;
    nargs = 1 ;
  } else if (!stricmp(option, "aseg")) {
    mri_aseg = MRIread(argv[2]) ;
    if (mri_aseg == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read %s", Progname, argv[2]) ;
    nargs = 1 ;
    printf("using aseg volume %s to constrain wm voxels\n", argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "laplace") || !stricmp(option, "laplacian")) {
    laplace_thick = 1 ;
    laplace_res = atof(argv[2]) ;
    fprintf(stderr,  "using Laplacian thickness measurement with PDE resolution = %2.3fmm\n",laplace_res) ;
    nargs = 1 ;
  } else if (!stricmp(option, "overlay")) {
    overlay_fname = argv[2] ;
    overlay_t0 = atoi(argv[3]) ;
    overlay_t1 = atoi(argv[4]) ;
    fprintf(stderr,  "computing overlay as average in interval [%d %d]\n",
            overlay_t0, overlay_t1) ;
    nargs = 3 ;
  } else if (!stricmp(option, "ratio")) {
    int i ;
    overlay_fname = argv[2] ;
    ratio = 1 ;

    for (i = 0 ; i < 4 ; i++)
      ratio_offsets[i] = atoi(argv[3+i]) ;
    fprintf(stderr,  "computing overlay as ratio of intervals [%d %d] / [%d %d]\n",
            ratio_offsets[0], ratio_offsets[1], ratio_offsets[2], ratio_offsets[3]);
    nargs = 5 ;
  } else if (!stricmp(option, "nsigma")) {
    normal_sigma = atof(argv[2]) ;
    fprintf(stderr,  "applying surface normal smoothing with sigma=%2.2f\n",
            normal_sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "norm_gw")) {
    norm_gw = atof(argv[2]) ;
    fprintf(stderr,  "normalizing intensities to have be %2.2f at gray/white boundary\n", norm_gw) ;
    nargs = 1 ;
  } else if (!stricmp(option, "norm_white")) {
    norm_white = atof(argv[2]) ;
    fprintf(stderr,  "normalizing intensities to have be %2.2f relative to white matter\n", norm_white) ;
    nargs = 1 ;
  } else if (!stricmp(option, "norm_csf")) {
    norm_csf = atof(argv[2]) ;
    fprintf(stderr,  "normalizing intensities to have be %2.2f relative to csf\n", norm_csf) ;
    nargs = 1 ;
  } else if (!stricmp(option, "norm_pial")) {
    norm_pial = atof(argv[2]) ;
    fprintf(stderr,  "normalizing intensities to have be %2.2f at gray/csf boundary\n", norm_pial) ;
    nargs = 1 ;
  } else if (!stricmp(option, "norm_mid")) {
    norm_mid = atof(argv[2]) ;
    fprintf(stderr,  "normalizing intensities to have be %2.2f at mid point of ribbon\n", norm_mid) ;
    nargs = 1 ;
  } else if (!stricmp(option, "norm_median")) {
    norm_median = atof(argv[2]) ;
    fprintf(stderr,  "normalizing intensities to have median value %2.3f\n", norm_median) ;
    nargs = 1 ;
  } else if (!stricmp(option, "norm_mean")) {
    norm_mean = atof(argv[2]) ;
    fprintf(stderr,  "normalizing intensities to have mean value %2.3f\n", norm_mean) ;
    nargs = 1 ;
  } else if (!stricmp(option, "wmnorm")) {
    wm_norm_fname = argv[2] ;
    fprintf(stderr,  "reading wm normalization volume from %s\n", wm_norm_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "erode")) {
    num_erode = atoi(argv[2]) ;
    fprintf(stderr,  "eroding label %d times\n", num_erode) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sdir")) {
    sdir = argv[2] ;
    nargs = 1 ;
  } else if (!stricmp(option, "use_normal")) {
    use_normal = atoi(argv[2]) ;
    printf("sampling along surface normal\n") ;
    nargs = 1 ;
  } else if (!stricmp(option, "use_pial")) {
    use_pial = atoi(argv[2]) ;
    if (use_pial == 0)
      use_normal = 1 ;  // must use surface normal if no pial surface
    printf("%susing pial surface to compute profiles\n", use_pial ? "" : "not ") ;
    nargs = 1 ;
  } else if (!stricmp(option, "normalize")) {
    normalize = 1 ;
    printf("normalizing profiles to be same length\n") ;
  } else if (!stricmp(option, "inorm")) {
    inorm = 1 ;
    printf("normalizing intensities to mean background level\n") ;
  } else if (!stricmp(option, "nsamples") || !stricmp(option, "samples")) {
    max_samples = atoi(argv[2]) ;
    normalize = 1 ;
    nargs = 1 ;
    printf("normalizing profiles to have %d samples\n", max_samples) ;
  } else if (!stricmp(option, "max")) {
    max_thick = atof(argv[2]) ;
    fprintf(stderr,  "limiting maximum cortical thickness to %2.2f mm.\n",
            max_thick) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sphere")) {
    sphere_name = argv[2] ;
    fprintf(stderr,  "using %s as spherical surface\n", sphere_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "mean")) {
    mean_outname = argv[2];
    fprintf(stderr,  "output the intensity profile mean to file %s.\n",
            mean_outname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "xform") ||
             !stricmp(option, "at")
            ) {
    xform_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "transform file name is %s\n", xform_fname);
  } else if (!stricmp(option, "ait")
            ) {
    xform_fname = argv[2];
    invert = 1;
    nargs = 1;
    fprintf(stderr, "Inversely apply the transform given by %s\n", xform_fname);
  } else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "invert")) {
    invert = 1;
    fprintf(stderr, "apply the inverse of the given registration transform\n");
  } else if (!stricmp(option, "openmp")) {
    char str[STRLEN] ;
    sprintf(str, "OMP_NUM_THREADS=%d", atoi(argv[2]));
    putenv(str) ;
#ifdef HAVE_OPENMP
    omp_set_num_threads(atoi(argv[2]));
#else
    fprintf(stderr, "Warning - built without openmp support\n");
#endif
    nargs = 1 ;
    fprintf(stderr, "Setting %s\n", str) ;
  } else if (!stricmp(option, "lta_src") ||
             !stricmp(option, "src")
            ) {
    fprintf(stderr, "src volume for the given transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the src volume...\n");
    lta_src = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_src) {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else if (!stricmp(option, "lta_dst") ||
             !stricmp(option, "dst")
            ) {
    fprintf(stderr, "dst volume for the transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the dst volume...\n");
    lta_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_dst) {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else switch (toupper(*option)) {
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("using write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'M':
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      printf("using momentum = %2.1f\n", parms.momentum) ;
      break ;
    case 'Z':
      zero_mean = 1 ;
      printf("making profiles zero mean\n") ;
      break ;
    case 'B':
      remove_bad = atoi(argv[2]) ;
      nargs = 1 ;
      printf("%sremoving bad intensity profiles\n", remove_bad ? "" : "not ") ;
      break ;
    case 'E':
      extra = atoi(argv[2]) ;
      nargs = 1 ;
      printf("adding %d extra points to profile\n", extra) ;
      max_samples += 2*extra ;
      break ;
    case 'C':
      curv_fname = argv[2] ;
      curv_thresh = atoi(argv[3]) ;
      printf("reading curvature from %s, and limiting calculations to %s regions\n",
             curv_fname, curv_thresh > 0 ? "sulcal" : "gyral") ;
      nargs = 2 ;
      break ;
    case 'P':
      polyfit = atoi(argv[2]) ;
      printf("using %dth order polynomial fit\n", polyfit) ;
      nargs = 1 ;
      break ;
    case 'Q':
      quadfit = 1 ;
      printf("writing out best-fitting quadratic curvature\n") ;
      break ;
    case 'L':
      label_names[nlabels] = argv[2] ;
      nargs = 1 ;
      printf("limiting profile calculation to label %s\n", label_names[nlabels]) ;
      nlabels++ ;
      break ;
    case 'N':
      nbhd_size = atoi(argv[2]) ;
      fprintf(stderr, "using neighborhood size=%d\n", nbhd_size) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1;
      printf("smoothing profiles %d times across space\n", navgs) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'T':
      thresh = atof(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <subject name> <hemi> <volume> <output file>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program computes the intensity profile of the cortical ribbon\n"
          "and writes the resulting measurement into a 'curvature' file\n"
          "<output file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-write_surf %%s write the variational pial surface target locations \n") ;
  fprintf(stderr, "-sdir %%s specifies the SUBJECTS_DIR \n") ;
  fprintf(stderr, "-white %%s specifies WHITE surface filename.\n") ;
  fprintf(stderr, "-pial %%s specifies PIAL surface filename.\n") ;
  fprintf(stderr, "-normalize normalized sampling w.r.t. thickness  \n") ;
  fprintf(stderr, "-mean %%s outputs the mean profile-integral to the specified file (output is in curv format).\n") ;
  fprintf(stderr, "-xform %%s specify the registration xform that maps the T1 volume (from which the surfaces were computed) to the input volume to be sampled.\n") ;
  fprintf(stderr, "-src %%s source volume used when computing the registration xform\n");
  fprintf(stderr, "-dst %%s destination volume used when computing the registration xform\n");
  fprintf(stderr, "-invert this flag asks to apply the registration xform inversely \n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}


MRI *
MRISfindNearestVerticesAndMeasureCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, int nbhd_size,
							   float max_thick,
							   int normalize, int curv_thresh, float *norm) {
  int     vno, n, vlist[100000], vtotal, ns, i, vnum, min_n,
  pial_vno, nsamples, min_vox ;
  float   d, dx, dy, dz, dist, min_dist, nx, ny, nz, dot, sample_dist, thick,
  white_mode, gray_mode, min_gray, max_gray, csf_mode, min_gray_at_wm ;
  double  x, y, z, xv, yv, zv, val ;
  MRI     *mri_profiles ;
  float   scale = 1 ;

  MRIScomputeClassModes(mris, mri, &white_mode, &gray_mode, &csf_mode, NULL, NULL, NULL) ;
  min_gray = (gray_mode + 2*csf_mode)/3 ;
  min_gray_at_wm = (1.5*gray_mode + csf_mode)/2.5 ;
  max_gray = (2*white_mode+gray_mode)/3 ;
  printf("bounding GM intensities in [%2.0f,%2.0f]\n", min_gray, max_gray) ;
  mri_profiles = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, max_samples) ;

  /* current vertex positions are gray matter, orig are white matter */
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX * const v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv) ;
    if (MRIindexNotInVolume(mri, xv, yv, zv))
      v->ripflag = 1 ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    if (!(vno % 25000))
      fprintf(stdout, "%d of %d vertices processed\n", vno,mris->nvertices) ;
    VERTEX * const v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (curv_thresh != 0 && curv_thresh*v->curv < 0)
      continue ;
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->x - v->origx ;
    dy = v->y - v->origy ;
    dz = v->z - v->origz ;
    pial_vno = vno ;
    min_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    v->marked = 1 ;
    vtotal = 1 ;
    vlist[0] = vno ;
    min_n = 0 ;
    for (ns = 1 ; ns <= nbhd_size ; ns++) {
      vnum = 0 ;  /* will be # of new neighbors added to list */
      for (i = 0 ; i < vtotal ; i++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[i]];
        VERTEX                * const vn  = &mris->vertices         [vlist[i]] ;
        if (vn->ripflag)
          continue ;
        if (vn->marked && vn->marked < ns-1)
          continue ;
        for (n = 0 ; n < vnt->vnum ; n++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n]] ;
          if (vn2->ripflag || vn2->marked)  /* already processed */
            continue ;
          vlist[vtotal+vnum++] = vnt->v[n] ;
          vn2->marked = ns ;
          dx = vn2->x-v->origx ;
          dy = vn2->y-v->origy ;
          dz = vn2->z-v->origz ;
          dot = dx*nx + dy*ny + dz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dot = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist) {
            min_n = ns ;
            min_dist = dist ;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON)
              fprintf(stdout, "%d --> %d = %2.3f\n",
                      vno,vnt->v[n], dist) ;
            pial_vno = vnt->v[n] ;
          }
        }
      }
      vtotal += vnum ;
    }

    // unmark stuff for next time
    for (n = 0 ; n < vtotal ; n++) {
      VERTEX * const vn = &mris->vertices[vlist[n]] ;
      if (vn->ripflag)
        continue ;
      vn->marked = 0 ;
    }

    pial_vno = vno ;  // disable shortest distance!!!
    {
      VERTEX const * const vn2 = &mris->vertices[pial_vno] ;
      if (use_normal)
      {
        dx = vn2->x-v->origx ;
        dy = vn2->y-v->origy ;
        dz = vn2->z-v->origz ;
        thick = sqrt(dx*dx + dy*dy + dz*dz) ;
        dx = v->nx ; dy = v->ny ; dz = v->nz ;
      }
      else  // use vector pointing from white to pial
      {
        dx = vn2->x-v->origx ;
        dy = vn2->y-v->origy ;
        dz = vn2->z-v->origz ;
        thick = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (FZERO(thick) == 0)
        {
          dx /= thick ;
          dy /= thick ;
          dz /= thick ;
        }
      }
    }
    
    if (use_pial == 0)
      thick = max_thick ;  // uniform length for samples
#define SAMPLE_DIST 0.5
    if (Gdiag_no == vno)
      DiagBreak() ;
    if (normalize) {
      sample_dist = thick / (max_samples-1);
      if (norm_white > 0 && norm)
        scale = norm_white / norm[vno] ;
      if (norm_mean > 0)
      {
        for (nsamples = 0, d = 0.0, scale = 0.0 ; nsamples < max_samples ; d += sample_dist, nsamples++) 
        {
          x = v->origx + d*dx ; y = v->origy + d*dy ; z = v->origz + d*dz ;
          MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
          MRIsampleVolume(mri, xv, yv, zv, &val) ;
          scale += val ;
        }
        scale /= (double) nsamples ;
        scale = norm_mean / scale ;
      }
      if (norm_pial > 0 || norm_mid > 0 || norm_gw > 0)
      {
        if (norm_pial > 0)
        {
          scale = norm_pial ;
          d = thick ;
        }
        else if (norm_mid > 0)
        {
          scale = norm_mid ;
          d = thick/2 ;
        }
        else
        {
          d = 0 ;
          scale = norm_gw ;
        }
        x = v->origx + d*dx ; y = v->origy + d*dy ; z = v->origz + d*dz ;
        MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        if (!FZERO(val))
          scale = scale / val ;
        else
          scale = 1 ;
      }
      for (nsamples = 0, d = 0.0 ; nsamples < max_samples ; d += sample_dist, nsamples++) {
        x = v->origx + d*dx ; y = v->origy + d*dy ; z = v->origz + d*dz ;
        MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        val *= scale ;
        MRIFseq_vox(mri_profiles, vno, 0, 0, nsamples) = val ;
      }
    } else {
      sample_dist = SAMPLE_DIST * mri->xsize ; 
      //   for (nsamples = 0, d = 0.0 ; d <= max_thick ; d += sample_dist, nsamples++)
      if (norm_mean > 0)
      {
        for (nsamples = 0, d = 0.0, scale = 0.0 ; nsamples < max_samples ; d += sample_dist, nsamples++) 
        {
          x = v->origx + d*dx ; y = v->origy + d*dy ; z = v->origz + d*dz ;
          MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
          MRIsampleVolume(mri, xv, yv, zv, &val) ;
          scale += val ;
        }
        scale /= (double) nsamples ;
        scale = norm_mean / scale ;
      }
      if (norm_pial > 0 || norm_mid > 0 || norm_gw > 0)
      {
        if (norm_pial > 0)
        {
          scale = norm_pial ;
          d = thick ;
        }
        else if (norm_mid > 0)
        {
          scale = norm_mid ;
          d = thick/2 ;
        }
        else
        {
          d = 0 ;
          scale = norm_gw ;
        }
        x = v->origx + d*dx ; y = v->origy + d*dy ; z = v->origz + d*dz ;
        MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        if (!FZERO(val))
          scale = scale / val ;
        else
          scale = 1 ;
      }
      for (nsamples = 0, d = -extra*sample_dist ;nsamples < max_samples; d += sample_dist, nsamples++) {
#if 0
        if (d > thick)  // so each row has the same # of cols
        {
          MRIFseq_vox(mri_profiles, vno, 0, 0, nsamples) = -1 ;
        } else
#endif
          // otherwise...
        {
          x = v->origx + d*dx ;
          y = v->origy + d*dy ;
          z = v->origz + d*dz ;
          MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
          MRIsampleVolume(mri, xv, yv, zv, &val) ;
          val *= scale ;
          MRIFseq_vox(mri_profiles, vno, 0, 0, nsamples) = val ;
        }
      }
    }
  }


  // now check profiles and correct them where they are obviously wrong
  min_vox = (int)ceil(MIN_BORDER_DIST/mri->xsize) ;  // within 2cm of edge
  if (FZERO(norm_gw) && 0)
  {
    int  bad ;

    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      if (!(vno % 25000))
        fprintf(stdout, "%d of %d vertices processed\n", vno,mris->nvertices) ;
      VERTEX * const v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag)
        continue ;
      if (curv_thresh != 0 && curv_thresh*v->curv < 0)
        continue ;
      bad = 0 ;
      for (nsamples = 0 ; nsamples < 4 ; nsamples++) {
        val = MRIgetVoxVal(mri_profiles, vno, 0, 0, nsamples) ;
        if (val < min_gray_at_wm)
          bad = 1 ;
      }
      MRISsurfaceRASToVoxelCached(mris, mri, v->whitex, v->whitey, v->whitez, &xv, &yv, &zv) ;
      if (xv < min_vox || yv < min_vox || zv < min_vox ||
          xv > mri->width-min_vox || yv > mri->height-min_vox || zv > mri->depth-min_vox)
        continue ;   // near edge of FOV
      if (bad == 0)  // check to make sure it's not all too dark
      {
        for (nsamples = 0 ; nsamples < 2*mri_profiles->nframes/3 ; nsamples++) {
          val = MRIgetVoxVal(mri_profiles, vno, 0, 0, nsamples) ;
          if (val < min_gray)
            bad = 1 ;
        }
      }
      if (bad)
        v->ripflag = 1 ;
    }
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag == 1) {
      int n ;
      for (n = 0 ; n < vt->vnum ; n++)
        mris->vertices[vt->v[n]].ripflag = 2 ;
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX * const v = &mris->vertices[vno] ;
    if (v->ripflag == 2)
      v->ripflag = 1 ;
  }
  //  MRISunrip(mris) ;
  return(mri_profiles) ;
}
static MRI *
MRIfitQuadratic(MRI *mri_profiles, MRI *mri_quad) {
  int     width, depth, height, x, y, z, wsize, whalf, t, tk, ti, nsamples ;
  float   max_curv, a /*, b, c */ ;
  MATRIX  *mX, *mXinv = NULL ;
  VECTOR  *vY, *vP = NULL ;
  max_curv=0;

  wsize = 5 ;
  whalf = (wsize-1)/2 ;
  mX = MatrixAlloc(wsize, 3, MATRIX_REAL) ;
  vY = VectorAlloc(wsize, MATRIX_REAL) ;

  width = mri_profiles->width ;
  height = mri_profiles->height ;
  depth = mri_profiles->depth ;
  nsamples = mri_profiles->nframes ;
  if (mri_quad == NULL) {
    mri_quad = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_profiles, mri_quad) ;
  }

  for (x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      for (z = 0 ; z < depth ; z++) {
        if (x == Gdiag_no)
          DiagBreak() ;
        max_curv = 0 ;
        for (t = whalf ; t < nsamples-whalf ; t++) {
          for (tk = -whalf ; tk <= whalf ; tk++) {
            ti = t+tk ;
            *MATRIX_RELT(mX, tk+whalf+1, 1) = tk*tk ;
            *MATRIX_RELT(mX, tk+whalf+1, 2) = tk ;
            *MATRIX_RELT(mX, tk+whalf+1, 3) = 1.0f ;

            VECTOR_ELT(vY, tk+whalf+1) = MRIgetVoxVal(mri_profiles, x, y, z, ti) ;
          }
          mXinv = MatrixPseudoInverse(mX, mXinv) ;
          if (mXinv == NULL)
            continue ;
          vP = MatrixMultiply(mXinv, vY, vP) ;
          a = RVECTOR_ELT(vP, 1) ;
//          b = RVECTOR_ELT(vP, 2) ;
//          c = RVECTOR_ELT(vP, 3);
          if (a > max_curv)
            max_curv = a ;
          MRIsetVoxVal(mri_quad, x, y, z, 0, max_curv) ;
        }
      }
    }
  }

  return(mri_quad) ;
}
static MRI *
MRIfitPolynomial(MRI *mri_profiles, MRI *mri_poly, int order) {
  int     width, depth, height, x, y, z, i, p ;
  MATRIX  *mX, *mXinv = NULL ;
  VECTOR  *vY, *vP = NULL ;
  double  d ;

  mX = MatrixAlloc(mri_profiles->nframes, order+1, MATRIX_REAL) ;
  vY = VectorAlloc(mri_profiles->nframes, MATRIX_REAL) ;

  width = mri_profiles->width ;
  height = mri_profiles->height ;
  depth = mri_profiles->depth ;
  if (mri_poly == NULL) {
    mri_poly = MRIallocSequence(width, height, depth, MRI_FLOAT, order+1) ;
    MRIcopyHeader(mri_profiles, mri_poly) ;
  }

  for (x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      for (z = 0 ; z < depth ; z++) {
        if (x == Gdiag_no)
          DiagBreak() ;
        for (i = 0 ; i < mri_profiles->nframes ; i++) {
          *MATRIX_RELT(mX, i+1, 1) = 1 ;
          for (p = 1 ; p <= order ; p++)
            *MATRIX_RELT(mX, i+1, p+1) = pow(i,p) ;
          VECTOR_ELT(vY, i+1) = MRIgetVoxVal(mri_profiles, x, y, z, i) ;
        }

        mXinv = MatrixPseudoInverse(mX, mXinv) ;
        if (mXinv == NULL)
          continue ;
        vP = MatrixMultiply(mXinv, vY, vP) ;
        for (p = 0 ; p <= order ; p++)
          MRIsetVoxVal(mri_poly, x, y, z, p, VECTOR_ELT(vP,p+1)) ;

        // replace intensity profile with fit version
        for (i = 0 ; i < mri_profiles->nframes ; i++) {
          d = VECTOR_ELT(vP,1) ;
          for (p = 1 ; p <= order ; p++)
            d += VECTOR_ELT(vP,p+1)*pow(i,p) ;
          MRIsetVoxVal(mri_profiles, x, y, z, i, d) ;
        }
      }
    }
  }

  return(mri_poly) ;
}
static int
remove_bad_profiles(MRI_SURFACE *mris, MRI *mri_profiles, MRI *mri, float border_mm) {
  int     nvox, vno, i, good, unknown_index, med_index  ;
  VERTEX  *v ;
  double  xv, yv, zv ;
  int     annot_index;

  if (MRISreadAnnotation(mris, "ad_aparc") != NO_ERROR) {
    if (MRISreadAnnotation(mris, "aparc") != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read annotation file ad_aparc", Progname) ;
  }

  CTABfindName(mris->ct, "Unknown", &unknown_index) ;
  CTABfindName(mris->ct, "Medial_wall", &med_index) ;
  printf("unknown index = %d, med = %d\n", unknown_index, med_index) ;

  nvox = (int)ceil(border_mm / mri->xsize) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv) ;
    CTABfindAnnotation(mris->ct, v->annotation, &annot_index);
    if (xv < nvox || yv < nvox || zv < nvox ||
        xv > mri->width-nvox || yv > mri->height-nvox || zv > mri->depth-nvox ||
        (annot_index == unknown_index) ||
        (annot_index == med_index))
      v->ripflag = 1 ;  // to prevent it from being included in the smoothing
    else {
      for (good = i = 0; i < mri_profiles->nframes ; i++) {
        if (!FZERO(MRIgetVoxVal(mri_profiles, vno, 0, 0, i)))
          good = 1 ;
      }
      if (good == 0)
        v->ripflag = 1 ;
    }
  }
  // blank all the bad profiles
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag == 0)
      continue ;
    for (i = 0; i < mri_profiles->nframes ; i++)
      MRIsetVoxVal(mri_profiles, vno, 0, 0, i, 0) ;
  }
  return(NO_ERROR) ;
}
static float *
mrisComputeWhiteMatterIntensities(MRI_SURFACE *mris, MRI *mri, float wsize_mm, float *norm, int navgs, MRI *mri_aseg, float mm_border)
{
  int     vno, whalf, num, xk, yk, zk, i, whalf_orig, label, erosions ;
  VERTEX  *v ;
  MRI     *mri_interior ;
  VECTOR  *v1, *v2 ;
  MATRIX  *m_vox2vox_aseg, *m_vox2vox ;
  float   avg_wm, res ;
  double  xv, yv, zv, xi, yi, zi, val ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  res = MAX(0.5, mri->xsize/2) ;
  mri_interior = MRISfillInterior(mris, res, NULL) ;
  erosions = MAX(nint(1.5/mri_interior->xsize), nint(mm_border/mri_interior->xsize)) ;
  for (i = 0 ; i < erosions ; i++)
    MRIerode(mri_interior, mri_interior) ;   // avoid partial volume voxels.
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_interior, "i.mgz") ;
  whalf_orig = (nint(wsize_mm / mri_interior->xsize)-1)/2 ;
  if (norm == NULL)
    norm = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (norm == NULL)
    ErrorExit(ERROR_NOMEMORY, "mrisComputeWhiteMatterIntensities: couldn't allocate array") ;


  m_vox2vox = MRIgetVoxelToVoxelXform(mri_interior, mri) ;
  if (mri_aseg)
    m_vox2vox_aseg = MRIgetVoxelToVoxelXform(mri_interior, mri_aseg) ;
  else
    m_vox2vox_aseg = NULL ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    MRISsurfaceRASToVoxel(mris, mri_interior, v->x, v->y, v->z, &xv, &yv, &zv) ;
    whalf = whalf_orig ;
    do
    {
      for (avg_wm = 0.0, num = 0, xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = xv+xk ;
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          yi = yv+yk ;
          for (zk = -whalf ; zk <= whalf ; zk++)
          {
            zi = zv+zk ;
            MRIsampleVolume(mri_interior, xi, yi, zi, &val) ;
            if (val > 0.5)  // in wm
            {
              V3_X(v1) = xi ; V3_Y(v1) = yi ; V3_Z(v1) = zi ;
              if (mri_aseg)
              {
                MatrixMultiply(m_vox2vox_aseg, v1, v2) ;
                MRIsampleVolumeType(mri_aseg, V3_X(v2), V3_Y(v2), V3_Z(v2), &val, SAMPLE_NEAREST) ;
              }
              else
                val = (float)Left_Cerebral_White_Matter ;

              label = (int)val ;
              if (IS_WM(label) == 0)
                continue ;
              MatrixMultiply(m_vox2vox, v1, v2) ;
              MRIsampleVolume(mri, V3_X(v2), V3_Y(v2), V3_Z(v2), &val) ;

              avg_wm += val ;
              num++ ;
            }
          }
        }
      }
      whalf += 3 ;
      if (whalf > 3*whalf_orig)
        break ;
    } while (num < 3) ;
    if (num > 0)
      norm[vno] = avg_wm / num;
    else
      norm[vno] = 1 ;  // couldn't find any wm
  }

  MRISimportValVector(mris, norm) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRISwriteValues(mris, "lh.norm.mgz") ;
  MRISaverageVals(mris, navgs) ;
  MRISexportValVector(mris, norm) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRISwriteValues(mris, "lh.norm.smooth.mgz") ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ;
  if (m_vox2vox_aseg)
    MatrixFree(&m_vox2vox_aseg) ;
  MRIfree(&mri_interior) ;
  return(norm) ;
}


int
MRISsampleFaceNormals(MHT *mht, MRI_SURFACE *mris, double x, double y, double z, int which, 
                      float *px, float *py, float *pz)
{
  float    d, dtotal, dx, dy, dz, xc, yc, zc ;
  int      n, fno ;
  double   fdist ;
  FACE     *face ;
  VERTEX   *v ;

  MHTfindClosestFaceGeneric(mht, mris, x, y, z, 1000, -1, -1, &face, &fno, &fdist) ;
  face = &mris->faces[fno] ;

  xc = yc = zc = 0.0 ; // to get rid of mac warnings
  for (dtotal = 0.0, n = 0  ; n < VERTICES_PER_FACE ; n++)
  {
    v = &mris->vertices[face->v[n]] ;
    MRISvertexCoord2XYZ_float(v, CANONICAL_VERTICES, &dx, &dy, &dz) ;
    dx -= x ; dy -= y ; dz -= z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    dtotal += d ;
  }

  *px = *py = *pz = 0 ;
  for (n = 0  ; n < VERTICES_PER_FACE ; n++)
  {
    v = &mris->vertices[face->v[n]] ;
    MRISvertexCoord2XYZ_float(v, CANONICAL_VERTICES, &dx, &dy, &dz) ;
    dx -= x ; dy -= y ; dz -= z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    d = 1-(2*d)/dtotal;
    if (d < 0)
    {
      DiagBreak() ;
      continue ;
    }

    MRISvertexCoord2XYZ_float(v, which, &xc, &yc, &zc) ;
    *px += d*xc ; *py += d*yc ; *pz += d*zc ;
  }


  return(NO_ERROR) ;
}
static MRI *
MRISmeasureCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, int nsamples) 
{
  MRI     *mri_profiles ;
  int     vno, n ;
  double  del, dist, dx, dy, dz, len, xv, yv, zv, x, y, z, val ;
  VERTEX  *v ;

  mri_profiles = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, nsamples);

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    dx = v->pialx - v->whitex ; dy = v->pialy - v->whitey ; dz = v->pialz - v->whitez ;
    len = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (len < 0.01)
      continue ;
    dx /= len ; dy /= len ; dz /= len ;
    del = len / (nsamples-1) ;
    for (dist = 0.0, n = 0 ; n < nsamples ; n++, dist += del) 
    {
      x = v->whitex + dist*dx ;
      y = v->whitey + dist*dy ;
      z = v->whitez + dist*dz ;
      MRISsurfaceRASToVoxel(mris, mri, x, y, z, &xv, &yv, &zv);
      MRIsampleVolumeFrameType(mri, xv, yv, zv, 0, SAMPLE_TRILINEAR, &val) ;
      MRIsetVoxVal(mri_profiles, vno, 0, 0, n, val) ;
    }
  }

  return(mri_profiles) ;
}

static MRI *
MRISmeasureNormalCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, double dist_in, double dist_out, int nsamples) 
{
  MRI     *mri_profiles ;
  int     vno, n ;
  double  del, dist, dx, dy, dz, len, xv, yv, zv, x, y, z, val ;
  VERTEX  *v ;

  mri_profiles = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, nsamples);

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    dx = v->nx ; dy = v->ny ; dz = v->nz ;
    len = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (len < 0.01)
      continue ;
    dx /= len ; dy /= len ; dz /= len ;
    len = dist_in + dist_out ;
    del = len / (nsamples-1) ;
    for (dist = -dist_in, n = 0 ; n < nsamples ; n++, dist += del) 
    {
      x = v->x + dist*dx ;
      y = v->y + dist*dy ;
      z = v->z + dist*dz ;
      MRISsurfaceRASToVoxel(mris, mri, x, y, z, &xv, &yv, &zv);
      MRIsampleVolumeFrameType(mri, xv, yv, zv, 0, SAMPLE_TRILINEAR, &val) ;
      MRIsetVoxVal(mri_profiles, vno, 0, 0, n, val) ;
      if (vno == Gdiag_no)
	printf("d %2.2f: (%2.1f, %2.1f, %2.1f) --> (%d, %d, %d): %2.1f\n",
	       dist, x, y, z, (int)(nint(xv)),(int)(nint(yv)),(int)(nint(zv)),val);
    }
  }

  return(mri_profiles) ;
}
