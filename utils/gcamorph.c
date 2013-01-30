/**
 * @file  gcamorph.c
 * @brief Utilities to morph the Gaussian Classifier Atlas (gca) data
 *
 * Reference:
  * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
  * Structures in the Human Brain", Fischl et al.
  * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2013/01/30 21:45:04 $
 *    $Revision: 1.281 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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

// Control which portions are done on the GPU
#define SHOW_EXEC_LOC 0

#ifdef FS_CUDA
// GCAmorphGPU requires Fermi
#ifdef GCAMORPH_ON_GPU
#define GCAM_CMP_GPU

#define GCAM_LABEL_ENERGY_GPU
#define GCAM_LLENERGY_GPU
#define GCAM_JACOBENERGY_GPU
#define GCAM_COMPUTE_RMS_GPU
#define GCAM_SMOOTHNESS_ENERGY_GPU
#define GCAM_FIND_OPTIMAL_TIMESTEP_GPU

#define GCAM_SMOOTH_TERM_GPU
#define GCAM_JACOB_TERM_GPU
#define GCAM_LL_TERM_GPU

#define GCAM_LABEL_TERM_MAINLOOP_GPU
#define GCAM_LABEL_TERM_REMOVE_OUTLIERS_GPU
#define GCAM_LABEL_TERM_COPYDELTAS_GPU
#define GCAM_LABEL_TERM_POSTANT_GPU
#define GCAM_LABEL_TERM_FINAL_GPU
#define GCAM_LABEL_TERM_GPU

#define GCAM_SMOOTH_GRADIENT_GPU

#define GCAM_COMPUTE_GRADIENT_GPU

#define GCAM_COPY_NODE_POSITIONS_GPU

#define GCAM_REGISTER_LEVEL_GPU
#define GCAM_REGISTER_PIPELINE_GPU

#else
// Have to turn everything off
#endif
#endif


#define MALLOC_CHECK_ 2

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "error.h"
#include "fio.h"
#include "diag.h"
#include "gca.h"
#include "gcamorph.h"
#include "transform.h"
#include "macros.h"
#include "proto.h"
#include "mrimorph.h"
#include "mrinorm.h"
#include "mriBSpline.h"
#include "matrix.h"
#include "cma.h"
#include "mri.h"
#include "tags.h"
#include "utils.h"
#include "fio.h"
#include "mri_circulars.h"

#include "chronometer.h"

#ifdef FS_CUDA
#include "gcamfots_cuda.h"
#endif

#include "gcamorphtestutils.h"

#if WITH_DMALLOC
#include <dmalloc.h>
#endif

#define GCAM_VERSION   1.0
#define MIN_STD (2.0)
#define MIN_VAR (MIN_STD*MIN_STD)

#ifndef FSIGN
#define FSIGN(f)  (((f) < 0) ? -1 : 1)
#endif

extern const char* Progname;

int dtrans_label_to_frame(GCA_MORPH_PARMS *mp, int label) ;
MRI *MRIcomposeWarps(MRI *mri_warp1, MRI *mri_warp2, MRI *mri_dst) ;
int fix_borders(GCA_MORPH *gcam)  ;

int gcam_write_grad = 0 ;
int gcam_write_neg = 0 ;

#if 1
int dtrans_labels[] =
{
  Left_Thalamus_Proper,
  Right_Thalamus_Proper,
  Left_Putamen,
  Right_Putamen,
  Left_Pallidum,
  Right_Pallidum,
  Left_Lateral_Ventricle,
  Right_Lateral_Ventricle,
  Left_Caudate,
  Right_Caudate,
  Left_Cerebral_White_Matter,
  Right_Cerebral_White_Matter,
  Left_Hippocampus,
  Right_Hippocampus,
  Left_Amygdala,
  Right_Amygdala,
  Left_VentralDC,
  Right_VentralDC,
  Brain_Stem,
  Left_Inf_Lat_Vent,
  Right_Inf_Lat_Vent,
  Left_Cerebral_Cortex,
  Right_Cerebral_Cortex,
  Left_Cerebellum_White_Matter,
  Right_Cerebellum_White_Matter,
  Left_Cerebellum_Cortex,
  Right_Cerebellum_Cortex,
  Unknown,
} ;
int combine_labels[][2] =
{
  { Right_Pallidum, Left_Pallidum},
  { Right_Putamen, Left_Putamen},
  { Right_Amygdala, Left_Amygdala},
  { Right_Hippocampus, Left_Hippocampus},
  { Left_Inf_Lat_Vent, Left_Lateral_Ventricle},
  { Right_Inf_Lat_Vent, Right_Lateral_Ventricle},
  { Right_Caudate, Left_Caudate},
  { Right_VentralDC, Left_VentralDC},
} ;

#else
int combine_labels[][2] =
{
  { Right_Pallidum, Left_Pallidum},
  { Right_Putamen, Left_Putamen},
  { Right_Amygdala, Left_Amygdala},
  { Right_Hippocampus, Left_Hippocampus},
  { Left_Inf_Lat_Vent, Left_Lateral_Ventricle},
  { Right_Inf_Lat_Vent, Right_Lateral_Ventricle},
  { Right_Caudate, Left_Caudate},
  { Right_VentralDC, Left_VentralDC},
} ;
int dtrans_labels[] =
{
  Left_Lateral_Ventricle,
  Right_Lateral_Ventricle
} ;
#endif

#define _NDTRANS_LABELS (sizeof(dtrans_labels) / sizeof(dtrans_labels[0]))
#define _NCOMBINE_LABELS (sizeof(combine_labels) / (2*sizeof(combine_labels[0][0])))
int NCOMBINE_LABELS = _NCOMBINE_LABELS ;
int NDTRANS_LABELS = _NDTRANS_LABELS ;

//static int GCAMsetNegativeNodeStatus(GCA_MORPH *gcam, int status) ;
HISTOGRAM *gcamJacobianHistogram(GCA_MORPH *gcam, HISTOGRAM *h);
int gcamComputePeriventricularWMDeformation(GCA_MORPH *gcam, MRI *mri) ;
double gcamMaxGradient(GCA_MORPH *gcam) ;
int gcamCheck(GCA_MORPH *gcam, MRI *mri) ;
int gcamWriteDiagnostics(GCA_MORPH *gcam) ;
int gcamShowCompressed(GCA_MORPH *gcam, FILE *fp) ;
MATRIX *gcamComputeOptimalTargetLinearTransform(GCA_MORPH *gcam,
    MATRIX *m_L,
    double reg) ;
int    gcamApplyLinearTransform(GCA_MORPH *gcam, MATRIX *m_L) ;
int gcamComputeTargetGradient(GCA_MORPH *gcam) ;
int gcamExpansionTerm(GCA_MORPH *gcam, MRI *mri, double l_expansion) ;
double gcamExpansionEnergy(GCA_MORPH *gcam, MRI *mri) ;
int gcamSetTotalMovement(GCA_MORPH *gcam, float alpha) ;
int gcamCheckJacobian(GCA_MORPH *gcam, float min_j, float max_j) ;
int gcamConstrainJacobian(GCA_MORPH *gcam,
                          MRI *mri,
                          GCA_MORPH_PARMS *parms) ;
int gcamComputeMostLikelyDirection(GCA_MORPH *gcam,
                                   MRI *mri,
                                   double x0, double y0, double z0,
                                   double target, MRI *mri_kernel,
                                   MRI *mri_nbhd,
                                   double *pdx, double *pdy, double *pdz) ;
NODE_LOOKUP_TABLE *gcamCreateNodeLookupTable(GCA_MORPH *gcam,
    MRI *mri,
    NODE_LOOKUP_TABLE *nlt) ;
//static int gcamSetGradientToNbrAverage(GCA_MORPH *gcam, int x, int y, int z);
int gcamFreeNodeLookupTable(NODE_LOOKUP_TABLE **pnlt) ;
MRI *gcamCreateJacobianImage(GCA_MORPH *gcam) ;

#if 0
static int gcamMLElabelAtLocation(GCA_MORPH *gcam,
                                  int x, int y, int z, float *vals) ;
#endif

int gcamRemoveNegativeNodes(GCA_MORPH *gcam,
                            MRI *mri,
                            GCA_MORPH_PARMS *parms) ;
int gcamRemoveCompressedNodes(GCA_MORPH *gcam,
                              MRI *mri,
                              GCA_MORPH_PARMS *parms,
                              float compression_ratio) ;
double gcamFindOptimalTimeStep(GCA_MORPH *gcam,
                               GCA_MORPH_PARMS *parms,
                               MRI *mri) ;
int  gcamAreaTermAtNode(GCA_MORPH *gcam, double l_area,
                        int i, int j, int k, double *pdx, double *pdy,
                        double *pdz) ;
int  gcamVolumeChangeTermAtNode(GCA_MORPH *gcam,
                                MRI *mri, double l_area,
                                int i, int j, int k,
                                double *pdx, double *pdy, double *pdz) ;

int   finitep(float f) ;

int write_snapshot(GCA_MORPH *gcam,
                   MRI *mri,
                   GCA_MORPH_PARMS *parms,
                   int iter) ;
int  log_integration_parms(FILE *fp, GCA_MORPH_PARMS *parms) ;
int gcamLimitGradientMagnitude(GCA_MORPH *gcam,
                               GCA_MORPH_PARMS *parms,
                               MRI *mri) ;


int gcamMultiscaleTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth,
                       double l_multiscale) ;
int gcamDistanceTransformTerm( GCA_MORPH *gcam,
                               MRI *mri_source,
                               MRI *mri,
                               double l_dtrans, GCA_MORPH_PARMS *mp);
double gcamDistanceTransformEnergy( GCA_MORPH *gcam,
				    MRI *mri_source,
                                    MRI *mri,
                                    GCA_MORPH_PARMS *mp) ;
int gcamLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth,
                       double l_likelihood, GCA_MORPH_PARMS *parms) ;
double gcamMapEnergy(GCA_MORPH *gcam, MRI *mri) ;

double gcamBinaryEnergy(GCA_MORPH *gcam, MRI *mri) ;
double gcamAreaIntensityEnergy(GCA_MORPH *gcam,
                               MRI *mri,
                               NODE_LOOKUP_TABLE *nlt) ;
int gcamMapTerm(GCA_MORPH *gcam,
                MRI *mri,
                MRI *mri_smooth,
                double l_map) ;

int gcamBinaryTerm(GCA_MORPH *gcam,
                   MRI *mri,
                   MRI *mri_smooth,
                   MRI *mri_dist,
                   double l_binary) ;
int gcamAreaIntensityTerm(GCA_MORPH *gcam,
                          MRI *mri, MRI *mri_smooth,
                          double l_area_intensity,
                          NODE_LOOKUP_TABLE *nlt, float sigma) ;



double gcamMultiscaleEnergy(GCA_MORPH *gcam, MRI *mri) ;

int gcamDistanceTerm(GCA_MORPH *gcam, MRI *mri, double l_distance) ;
double gcamDistanceEnergy(GCA_MORPH *gcam, MRI *mri) ;



int gcamLSmoothnessTerm(GCA_MORPH *gcam,
                        MRI *mri,
                        double l_smoothness) ;
int gcamElasticTerm( const GCA_MORPH *gcam, GCA_MORPH_PARMS *parms ) ;
double gcamElasticEnergy( const GCA_MORPH *gcam, GCA_MORPH_PARMS *parms ) ;
double gcamLSmoothnessEnergy(GCA_MORPH *gcam, MRI *mri) ;

int gcamSpringTerm(GCA_MORPH *gcam,
                   double l_spring,
                   double ratio_thresh) ;
double gcamSpringEnergy(GCA_MORPH *gcam, double ratio_thresh) ;

//static int gcamInvalidSpringTerm(GCA_MORPH *gcam, double l_spring) ;

int gcamAreaSmoothnessTerm(GCA_MORPH *gcam,
                           MRI *mri,
                           double l_jacobian) ;
int gcamAreaTerm(GCA_MORPH *gcam, double l_jacobian) ;
double gcamAreaEnergy(GCA_MORPH *gcam) ;



int check_gcam( const GCAM *gcam ) ;
MATRIX *gcamComputeOptimalLinearTransformInRegion
(GCA_MORPH *gcam, MRI *mri_mask, MATRIX *m_L,
 int x0, int y0, int z0, int whalf);
#define DEFAULT_PYRAMID_LEVELS 3
#define MAX_PYRAMID_LEVELS     20
#define MAX_EXP          200
#define GCAMN_SUB(mns1, mns2, v)          \
  V3_LOAD(v, mns1->x - mns2->x, mns1->y - mns2->y, mns1->z - mns2->z)

static int Ginvalid = 0 ;
static int Galigned = 0 ;

#define NODE_SAMPLE_VAR  (.25*.25)
#define MIN_NODE_DIST    (1)
#define MIN_NODE_DIST_SQ (MIN_NODE_DIST*MIN_NODE_DIST)

void GCAMwriteGeom( const GCA_MORPH *gcam, znzFile file )
{
  char buf[512];

  // src volume info
  znzwriteInt(gcam->image.valid, file);
  znzwriteInt(gcam->image.width, file);
  znzwriteInt(gcam->image.height, file);
  znzwriteInt(gcam->image.depth, file);
  znzwriteFloat(gcam->image.xsize, file) ;
  znzwriteFloat(gcam->image.ysize, file) ;
  znzwriteFloat(gcam->image.zsize, file) ;
  znzwriteFloat(gcam->image.x_r, file) ;
  znzwriteFloat(gcam->image.x_a, file) ;
  znzwriteFloat(gcam->image.x_s, file) ;
  znzwriteFloat(gcam->image.y_r, file) ;
  znzwriteFloat(gcam->image.y_a, file) ;
  znzwriteFloat(gcam->image.y_s, file) ;
  znzwriteFloat(gcam->image.z_r, file) ;
  znzwriteFloat(gcam->image.z_a, file) ;
  znzwriteFloat(gcam->image.z_s, file) ;
  znzwriteFloat(gcam->image.c_r, file) ;
  znzwriteFloat(gcam->image.c_a, file) ;
  znzwriteFloat(gcam->image.c_s, file) ;
  memset(buf, 0, 512*sizeof(char));
  strcpy(buf, gcam->image.fname);
  znzwrite(buf, sizeof(char), 512, file);
  // dst volume info
  znzwriteInt(gcam->atlas.valid, file);
  znzwriteInt(gcam->atlas.width, file);
  znzwriteInt(gcam->atlas.height, file);
  znzwriteInt(gcam->atlas.depth, file);
  znzwriteFloat(gcam->atlas.xsize, file) ;
  znzwriteFloat(gcam->atlas.ysize, file) ;
  znzwriteFloat(gcam->atlas.zsize, file) ;
  znzwriteFloat(gcam->atlas.x_r, file) ;
  znzwriteFloat(gcam->atlas.x_a, file) ;
  znzwriteFloat(gcam->atlas.x_s, file) ;
  znzwriteFloat(gcam->atlas.y_r, file) ;
  znzwriteFloat(gcam->atlas.y_a, file) ;
  znzwriteFloat(gcam->atlas.y_s, file) ;
  znzwriteFloat(gcam->atlas.z_r, file) ;
  znzwriteFloat(gcam->atlas.z_a, file) ;
  znzwriteFloat(gcam->atlas.z_s, file) ;
  znzwriteFloat(gcam->atlas.c_r, file) ;
  znzwriteFloat(gcam->atlas.c_a, file) ;
  znzwriteFloat(gcam->atlas.c_s, file) ;
  memset(buf, 0, 512*sizeof(char));
  strcpy(buf, gcam->atlas.fname);
  znzwrite(buf, sizeof(char), 512, file);
}


void GCAMreadGeom(GCA_MORPH *gcam, znzFile file)
{
  // src volume info
  gcam->image.valid = znzreadInt(file);
  gcam->image.width = znzreadInt(file);
  gcam->image.height = znzreadInt(file);
  gcam->image.depth = znzreadInt(file);
  gcam->image.xsize = znzreadFloat(file) ;
  gcam->image.ysize = znzreadFloat(file) ;
  gcam->image.zsize = znzreadFloat(file) ;
  gcam->image.x_r = znzreadFloat(file) ;
  gcam->image.x_a = znzreadFloat(file) ;
  gcam->image.x_s = znzreadFloat(file) ;
  gcam->image.y_r = znzreadFloat(file) ;
  gcam->image.y_a = znzreadFloat(file) ;
  gcam->image.y_s = znzreadFloat(file) ;
  gcam->image.z_r = znzreadFloat(file) ;
  gcam->image.z_a = znzreadFloat(file) ;
  gcam->image.z_s = znzreadFloat(file) ;
  gcam->image.c_r = znzreadFloat(file) ;
  gcam->image.c_a = znzreadFloat(file) ;
  gcam->image.c_s = znzreadFloat(file) ;
  memset(gcam->image.fname, 0, 512*sizeof(char));
  znzread(gcam->image.fname, sizeof(char), 512, file);
  // dst volume info
  gcam->atlas.valid = znzreadInt(file);
  gcam->atlas.width = znzreadInt(file);
  gcam->atlas.height = znzreadInt(file);
  gcam->atlas.depth = znzreadInt(file);
  gcam->atlas.xsize = znzreadFloat(file) ;
  gcam->atlas.ysize = znzreadFloat(file) ;
  gcam->atlas.zsize = znzreadFloat(file) ;
  gcam->atlas.x_r = znzreadFloat(file) ;
  gcam->atlas.x_a = znzreadFloat(file) ;
  gcam->atlas.x_s = znzreadFloat(file) ;
  gcam->atlas.y_r = znzreadFloat(file) ;
  gcam->atlas.y_a = znzreadFloat(file) ;
  gcam->atlas.y_s = znzreadFloat(file) ;
  gcam->atlas.z_r = znzreadFloat(file) ;
  gcam->atlas.z_a = znzreadFloat(file) ;
  gcam->atlas.z_s = znzreadFloat(file) ;
  gcam->atlas.c_r = znzreadFloat(file) ;
  gcam->atlas.c_a = znzreadFloat(file) ;
  gcam->atlas.c_s = znzreadFloat(file) ;
  memset(gcam->atlas.fname, 0, 512*sizeof(char));
  znzread(gcam->atlas.fname, sizeof(char), 512, file);
}


int
GCAMwrite( const GCA_MORPH *gcam, const char *fname )
{
  znzFile         file;
  //FILE            *fp=0 ;
  int             x, y, z ;
  GCA_MORPH_NODE  *gcamn ;
  int gzipped = 0;

  if (strstr(fname, ".m3z"))
  {
    //    printf("GCAMwrite:: m3z loop\n");
    gzipped = 1;
  }

  file = znzopen(fname, "wb", gzipped) ;
  if (znz_isnull(file))
  {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "GCAMwrite(%s): could not open file",fname)) ;
  }

  znzwriteFloat(GCAM_VERSION, file) ;
  znzwriteInt(gcam->width, file) ;
  znzwriteInt(gcam->height, file) ;
  znzwriteInt(gcam->depth, file) ;
  znzwriteInt(gcam->spacing, file) ;
  znzwriteFloat(gcam->exp_k, file) ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        znzwriteFloat(gcamn->origx, file) ;
        znzwriteFloat(gcamn->origy, file) ;
        znzwriteFloat(gcamn->origz, file) ;

        znzwriteFloat(gcamn->x, file) ;
        znzwriteFloat(gcamn->y, file) ;
        znzwriteFloat(gcamn->z, file) ;

        znzwriteInt(gcamn->xn, file) ;
        znzwriteInt(gcamn->yn, file) ;
        znzwriteInt(gcamn->zn, file) ;
      }
    }
  }
  znzwriteInt(TAG_GCAMORPH_GEOM, file);
  GCAMwriteGeom(gcam, file);

  znzwriteInt(TAG_GCAMORPH_TYPE, file) ;
  znzwriteInt(gcam->type, file) ;

  znzwriteInt(TAG_GCAMORPH_LABELS, file) ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        znzwriteInt(gcamn->label, file) ;
      }
    }
  }
  if (gcam->m_affine)
  {
    znzwriteInt(TAG_MGH_XFORM, file) ;
    //MatrixAsciiWriteInto(file, gcam->m_affine) ;
    znzWriteMatrix(file, gcam->m_affine);
  }

  znzclose(file);

  return(NO_ERROR) ;
}

/*-------------------------------------------------------------------------
  GCAMwriteInverse() - saves GCAM inverse to
  talairach.m3z.inv.{x,y,z}.mgh in the same directory as the GCAM (ie,
  mri/transforms). It does not write the GCAM. If the inverse already
  exists, it is NOT read in and will be overwritten.

  This function can be run in  one of two ways:
    1. GCAMwriteInverse(gcamfile,NULL) - reads in gcam, inverts,
       saves inverse, then frees gcam.
    2. GCAMwriteInverse(gcamfile,gcam) - pass gcam as arg. If
       it has not been inverted, then it computes the inverse.
       It then saves the inverse (does not free gcam).

  If the inverse must be computed, then it reads in the header for
  mri/orig.mgz. See also GCAMreadAndInvert().
  -----------------------------------------------------------------*/
int GCAMwriteInverse(const char *gcamfname, GCA_MORPH *gcam)
{
  char *gcamdir, *mridir, tmpstr[2000];
  MRI *mri;
  int freegcam;

  // Read in gcam if not passed
  freegcam=0;
  if (gcam == NULL)
  {
    printf("Reading %s \n",gcamfname);
    gcam = GCAMread(gcamfname);
    if (gcam == NULL)
    {
      return(1);
    }
    freegcam=1;
  }
  gcamdir  = fio_dirname(gcamfname);

  // Check whether inverse has been computed, if not, do so now
  if (gcam->mri_xind == NULL)
  {
    // Need a template MRI in order to compute inverse
    mridir   = fio_dirname(gcamdir);
    //sprintf(tmpstr,"%s/orig.mgz",mridir);
    sprintf(tmpstr, "%s", (gcam->image).fname);
    mri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
    if (mri==NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
      GCAMfree(&gcam);
      return(1);
    }
    printf("Inverting GCAM\n");
    GCAMinvert(gcam, mri);
    free(mridir);
  }

  printf("Saving inverse \n");
  sprintf(tmpstr,"%s/talairach.m3z.inv.x.mgz",gcamdir);
  MRIwrite(gcam->mri_xind,tmpstr);
  sprintf(tmpstr,"%s/talairach.m3z.inv.y.mgz",gcamdir);
  MRIwrite(gcam->mri_yind,tmpstr);
  sprintf(tmpstr,"%s/talairach.m3z.inv.z.mgz",gcamdir);
  MRIwrite(gcam->mri_zind,tmpstr);

  if (freegcam)
  {
    GCAMfree(&gcam);
  }
  free(gcamdir);

  return(0);
}

int GCAMwriteInverseNonTal(const char *gcamfname, GCA_MORPH *gcam)
{
  char tmpstr[2000];
  MRI *mri;
  int freegcam;

  // Read in gcam if not passed
  freegcam=0;
  if (gcam == NULL)
  {
    printf("Reading %s \n",gcamfname);
    gcam = GCAMread(gcamfname);
    if (gcam == NULL)
    {
      return(1);
    }
    freegcam=1;
  }

  // Check whether inverse has been computed, if not, do so now
  if (gcam->mri_xind == NULL)
  {
    // Need a template MRI in order to compute inverse
    sprintf(tmpstr, "%s", (gcam->image).fname);
    mri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
    if (mri==NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
      GCAMfree(&gcam);
      return(1);
    }
    printf("Inverting GCAM\n");
    GCAMinvert(gcam, mri);
  }

  printf("Saving inverse \n");
  sprintf(tmpstr,"%s.inv.x.mgz",gcamfname);
  MRIwrite(gcam->mri_xind,tmpstr);
  sprintf(tmpstr,"%s.inv.y.mgz",gcamfname);
  MRIwrite(gcam->mri_yind,tmpstr);
  sprintf(tmpstr,"%s.inv.z.mgz",gcamfname);
  MRIwrite(gcam->mri_zind,tmpstr);

  if (freegcam)
  {
    GCAMfree(&gcam);
  }

  return(0);
}
/*----------------------------------------------------------------------
  GCAMreadAndInvert() - reads morph and inverts. If the inverse
  exists on disk, then it is read in instead of being computed. The
  files that compose the inverse morph are talairach.m3z.inv.{x,y,z}.mgh,
  which are assumed to live in the same directory as the forward morph,
  which is assumed to be in mri/transforms. This is important because,
  if the morph must be explicitly inverted, it will read in the header
  for mri/orig.mgz (or die trying).
  ----------------------------------------------------------------------*/
GCA_MORPH *GCAMreadAndInvert(const char *gcamfname)
{
  GCA_MORPH *gcam;
  char tmpstr[2000], *gcamdir, *mridir;
  MRI *mri;
  int xexists, yexists, zexists;

  // Read GCAM
  gcam = GCAMread(gcamfname);
  if (gcam == NULL)
  {
    return(NULL);
  }

  gcamdir  = fio_dirname(gcamfname);

  // Check whether inverse morph (talairach.m3z.inv.{x,y,z}.mgh)  exists
  xexists=0;
  sprintf(tmpstr,"%s/talairach.m3z.inv.x.mgz",gcamdir);
  if (fio_FileExistsReadable(tmpstr))
  {
    xexists = 1;
  }
  yexists=0;
  sprintf(tmpstr,"%s/talairach.m3z.inv.y.mgz",gcamdir);
  if (fio_FileExistsReadable(tmpstr))
  {
    yexists = 1;
  }
  zexists=0;
  sprintf(tmpstr,"%s/talairach.m3z.inv.z.mgz",gcamdir);
  if (fio_FileExistsReadable(tmpstr))
  {
    zexists = 1;
  }

  if (xexists && yexists && zexists)
  {
    // Inverse found on disk, just read it in
    printf("Reading morph inverse\n");
    sprintf(tmpstr,"%s/talairach.m3z.inv.x.mgz",gcamdir);
    printf("Reading %s\n",tmpstr);
    gcam->mri_xind = MRIread(tmpstr);
    if (gcam->mri_xind == NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
    }

    sprintf(tmpstr,"%s/talairach.m3z.inv.y.mgz",gcamdir);
    printf("Reading %s\n",tmpstr);
    gcam->mri_yind = MRIread(tmpstr);
    if (gcam->mri_yind == NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
    }

    sprintf(tmpstr,"%s/talairach.m3z.inv.z.mgz",gcamdir);
    printf("Reading %s\n",tmpstr);
    gcam->mri_zind = MRIread(tmpstr);
    if (gcam->mri_zind == NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
    }
  }
  else
  {
    // Must invert explicitly
    printf("Inverting Morph\n");
    mridir   = fio_dirname(gcamdir);
    //  sprintf(tmpstr,"%s/orig.mgz",mridir);
    // Need template mri
    sprintf(tmpstr, "%s", (gcam->image).fname);
    mri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
    if (mri==NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
      return(NULL);
    }
    GCAMinvert(gcam, mri);
    free(mridir);
  }

  free(gcamdir);
  return(gcam) ;
}

GCA_MORPH *GCAMreadAndInvertNonTal(const char *gcamfname)
{
  GCA_MORPH *gcam;
  char tmpstr[2000], *gcamdir;
  MRI *mri;
  int xexists, yexists, zexists;

  // Read GCAM
  gcam = GCAMread(gcamfname);
  if (gcam == NULL)
  {
    return(NULL);
  }

  gcamdir  = fio_dirname(gcamfname);

  // Check whether inverse morph ($gcamfname.inv.{x,y,z}.mgh)  exists
  xexists=0;
  sprintf(tmpstr,"%s.inv.x.mgz",gcamfname);
  if (fio_FileExistsReadable(tmpstr))
  {
    xexists = 1;
  }
  yexists=0;
  sprintf(tmpstr,"%s.inv.y.mgz",gcamfname);
  if (fio_FileExistsReadable(tmpstr))
  {
    yexists = 1;
  }
  zexists=0;
  sprintf(tmpstr,"%s.inv.z.mgz",gcamfname);
  if (fio_FileExistsReadable(tmpstr))
  {
    zexists = 1;
  }

  if (xexists && yexists && zexists)
  {
    // Inverse found on disk, just read it in
    printf("Reading morph inverse\n");
    sprintf(tmpstr,"%s.inv.x.mgz",gcamfname);
    printf("Reading %s\n",tmpstr);
    gcam->mri_xind = MRIread(tmpstr);
    if (gcam->mri_xind == NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
    }

    sprintf(tmpstr,"%s.inv.y.mgz",gcamfname);
    printf("Reading %s\n",tmpstr);
    gcam->mri_yind = MRIread(tmpstr);
    if (gcam->mri_yind == NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
    }

    sprintf(tmpstr,"%s.inv.z.mgz",gcamfname);
    printf("Reading %s\n",tmpstr);
    gcam->mri_zind = MRIread(tmpstr);
    if (gcam->mri_zind == NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
    }
  }
  else
  {
    // Must invert explicitly
    printf("Inverting Morph\n");
    // Need template mri
    sprintf(tmpstr, "%s", (gcam->image).fname);
    if (!fio_FileExistsReadable(tmpstr)) {	// Look in parent directory
      char *mridir, altname[2000];

      mridir = fio_dirname(gcamdir);
      sprintf(altname, "%s/mri/norm.mgz", mridir);

      if (fio_FileExistsReadable(altname)) {
        printf("WARN: cannot acccess %s\n", tmpstr);
        printf("WARN: using %s\n", altname);
        strcpy(tmpstr, altname);
      }
      else {
        mridir = fio_dirname(mridir);
        sprintf(altname, "%s/mri/norm.mgz", mridir);

        if (fio_FileExistsReadable(altname)) {
          printf("WARN: cannot acccess %s\n", tmpstr);
          printf("WARN: using %s\n", altname);
          strcpy(tmpstr, altname);
        }
      }
    }
    mri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
    if (mri==NULL)
    {
      printf("ERROR: reading %s\n",tmpstr);
      return(NULL);
    }
    GCAMinvert(gcam, mri);
    GCAMwriteInverseNonTal(gcamfname, gcam);
  }

  free(gcamdir);
  return(gcam) ;
}

//! A thin pipeline segment to save a GPU copy
void GCAMregisterPipeline( GCA_MORPH *gcam,
                           MRI *mri,
                           MRI *mri_smooth,
                           GCA_MORPH_PARMS *parms,
                           double *last_rms,
                           int *level_steps,
                           int i )
{

#ifdef GCAM_REGISTER_PIPELINE_GPU
  printf( "%s: On GPU\n", __FUNCTION__ );

  GCAMregisterPipelineGPU( gcam, mri, mri_smooth, parms,
                           last_rms, level_steps, i );
#else
  *last_rms = GCAMcomputeRMS( gcam, mri, parms );
  if( i==0 )
  {
    parms->start_rms = *last_rms;
  }
  *level_steps = parms->start_t;
  GCAMregisterLevel( gcam, mri, mri_smooth, parms);
#endif
}


int
GCAMregister(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  char   fname[STRLEN] ;
  int    level, i, level_steps, navgs, l2, relabel,
         orig_relabel, start_t = 0, passno ;
  MRI    *mri_smooth = NULL, *mri_kernel ;
  double base_sigma, pct_change, rms, last_rms = 0.0, label_dist,
    orig_dt,l_smooth, start_rms=0.0, l_orig_smooth, l_elastic, l_orig_elastic ;

  if (FZERO(parms->min_sigma))
  {
    parms->min_sigma = 0.4 ;
  }

  // debugging
  if (mri)
  {
    MATRIX *m_vox2ras ;
    VECTOR  *v1, *v2 ;
    double  rtx, rty, rtz, rx, ry, rz, d, dmin, rxmin, rymin, rzmin ;
    int     xmin, ymin, zmin, whalf, xk, yk, zk, xi, yi, zi, x, y, z ;
    GCA_MORPH_NODE *gcamn_nbr, *gcamn ;

    xmin = ymin = zmin = 0.0 ;

    m_vox2ras = MRIgetVoxelToRasXform(mri) ;
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v1,4,1) = 1.0 ;
    *MATRIX_RELT(v2,4,1) = 1.0 ;

    dmin = 10000000 ;
    rtx = 29.5 ;
    rty = 38.75 ;
    rtz = -48.44 ;
    for (x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          if (gcamn->label == 0)
          {
            V3_X(v1) = gcamn->x ;
            V3_Y(v1) = gcamn->y ;
            V3_Z(v1) = gcamn->z ;
            MatrixMultiply(m_vox2ras, v1, v2) ;
            rx = V3_X(v2) ;
            ry = V3_Y(v2) ;
            rz = V3_Z(v2) ;
            d = sqrt(SQR(rtx-rx)+SQR(rty-ry)+SQR(rtz-rz)) ;
            if (d < dmin)
            {
              rxmin = rx ;
              rymin = ry ;
              rzmin = rz ;
              xmin = x ;
              ymin = y ;
              zmin = z ;
              dmin = d ;
            }
          }
        }
    whalf = 1 ;
    gcamn = &gcam->nodes[xmin][ymin][zmin] ;
    for (xk = -whalf ; xk <= whalf  ; xk++)
      for (yk = -whalf ; yk <= whalf ; yk++)
        for (zk = -whalf ; zk <= whalf ; zk++)
        {
          xi = xk+xmin ;
          yi = yk+ymin ;
          zi = zk+zmin ;
          if ((xi < 0 || xi >= gcam->width))
          {
            continue ;
          }
          if ((yi < 0 || yi >= gcam->height))
          {
            continue ;
          }
          if ((zi < 0 || zi >= gcam->depth))
          {
            continue ;
          }
          gcamn_nbr = &gcam->nodes[xi][yi][zi] ;
          if (gcamn_nbr->label > 0)
          {
            DiagBreak() ;
          }
        }


    MatrixFree(&m_vox2ras) ;
    VectorFree(&v1) ;
    VectorFree(&v2) ;
  }

  /*
     sorry, I know this is a hack, but I'm trying not to break any of the
     old morph stuff used for standard aseg generation. Here if both binary and
     intensity images are used, the binary *must* be passed
     in the parms structure,
     and it must already have been transformed to be in the same voxel space as
     the intensity image, otherwise things are too much of a mess.
  */
  if (!DZERO(parms->l_dtrans))
    parms->mri_dist_map =
        GCAMcreateDistanceTransforms(mri, parms->mri, NULL, 40.0, gcam, &parms->mri_atlas_dist_map);

  if (!DZERO(parms->l_binary))
  {
    MRI *mri_tmp ;

    parms->mri_dist_map =
      MRIdistanceTransform
      (parms->mri_binary, NULL, parms->target_label, -1,
       DTRANS_MODE_SIGNED, NULL);
#if 0
    if (DZERO(parms->l_area_intensity))
    {
      mri_tmp = MRIbinarize(mri, NULL, 1, 0, 1) ;
      parms->mri_binary = MRIchangeType(mri_tmp, MRI_FLOAT, 0, 1, 1) ;
      MRIfree(&mri_tmp) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        MRIwrite(parms->mri_binary, "bin.mgz") ;
      }
    }
    else  /* using both binary and intensity images for registration */
#endif
    {
      if (parms->mri_binary == NULL)
        ErrorExit
        (ERROR_BADPARM, \
         "GCAMregister: must fill parms->mri_binary "
         "pointer when area_intensity also specified") ;
      mri_tmp = MRIbinarize(parms->mri_binary, NULL, 1, 0, 1) ;
      parms->mri_binary = MRIchangeType(mri_tmp, MRI_FLOAT, 0, 1, 1) ;
      MRIfree(&mri_tmp) ;
    }
  }
#if 0
  else
  {
    parms->mri_dist_map = NULL ;
  }
#endif


  navgs = parms->navgs ;
  orig_dt = parms->dt ;
  l_orig_smooth = l_smooth = parms->l_smoothness;
  l_orig_elastic = l_elastic = parms->l_elastic;
  if (DZERO(parms->exp_k))
  {
    parms->exp_k = EXP_K ;
  }
  if (parms->levels < 0)
  {
    parms->levels = DEFAULT_PYRAMID_LEVELS ;
  }
  else if (parms->levels >= MAX_PYRAMID_LEVELS)
  {
    parms->levels = MAX_PYRAMID_LEVELS ;
  }

  gcam->exp_k = parms->exp_k ;
  parms->mri = mri ;
  //////////////////////////////////////////////////////////
  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s.log", parms->base_name) ;
    if (parms->log_fp == NULL)
    {
      if (parms->start_t == 0)
      {
        parms->log_fp = fopen(fname, "w") ;
      }
      else
      {
        parms->log_fp = fopen(fname, "a") ;
      }
    }
  }
  else
  {
    parms->log_fp = NULL ;
  }

  base_sigma = parms->sigma ;

  // GCAMinit() did this at the end
  // make node to have the max_prior label values
  if (parms->relabel_avgs >= parms->navgs && parms->relabel)
  {
    GCAMcomputeLabels(mri, gcam) ;
  }
  else
  {
    GCAMcomputeMaxPriorLabels(gcam) ;
  }

  orig_relabel = relabel = parms->relabel ;
  if (parms->npasses <= 0)
  {
    parms->npasses = 1 ;
  }
  for (passno = 0 ; passno < parms->npasses ; passno++)
  {
    printf("**************** pass %d of %d ************************\n",
           passno+1, parms->npasses);
    if (parms->log_fp)
      fprintf(parms->log_fp,
              "**************** pass %d of %d ************************\n",
              passno+1, parms->npasses);
    if (passno == parms->enable_zero_passes)
    {
      printf("enabling zero nodes\n") ;
      if (parms->log_fp)
	fprintf(parms->log_fp, "enabling zero nodes\n") ;
      GCAMremoveIgnoreZero(gcam, parms->mri) ;
    }
      
    parms->navgs = navgs ;
    label_dist = parms->label_dist ;
    for (level = parms->levels-1 ; level >= 0 ; level--)
    {
      rms = parms->start_rms ;
      if (parms->reset_avgs == parms->navgs)
      {
        printf("resetting metric properties...\n") ;
        GCAMcopyNodePositions
        (gcam, ORIGINAL_POSITIONS, SAVED_ORIGINAL_POSITIONS) ;
        GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, ORIGINAL_POSITIONS) ;
        gcamComputeMetricProperties(gcam) ;
        GCAMstoreMetricProperties(gcam) ;
      }
      parms->relabel = (parms->relabel_avgs >= parms->navgs) ;
      parms->sigma = base_sigma ;
      if (FZERO(parms->l_smoothness) && !FZERO(parms->l_elastic))
      {
        parms->l_elastic = l_elastic ;
	parms->l_elastic = l_elastic / ((parms->navgs)+1) ;
        printf("setting elastic coefficient to %2.3f\n", parms->l_elastic);
      }
      else
      {
#if 0
        parms->l_smoothness = l_smooth / (sqrt(parms->navgs)+1) ;
#else
        if (DZERO(parms->l_binary ) &&
            DZERO(parms->l_area_intensity) && parms->scale_smoothness)
        {
          parms->l_smoothness = l_smooth / ((parms->navgs)+1) ;
        }
        else
        {
          parms->l_smoothness = l_smooth ;
        }
#endif
        printf("setting smoothness coefficient to %2.3f\n",
               parms->l_smoothness);
      }
      for (l2 = 0 ; l2 < parms->levels ; l2++)  // different sigma levels
      {
        if (mri)
        {
          if (DZERO(parms->sigma))
          {
            mri_smooth = MRIcopy(mri, mri_smooth) ;
            if (parms->mri_binary)
              parms->mri_binary_smooth =
                MRIcopy(parms->mri_binary, parms->mri_binary_smooth) ;
          }
          else
          {
            //      if (Gdiag & DIAG_SHOW)
            printf("blurring input image with Gaussian with sigma=%2.3f...\n",
                   parms->sigma) ;
            mri_kernel = MRIgaussian1d(parms->sigma, 100) ;
            if (parms->mri_binary)
              parms->mri_binary_smooth = MRIconvolveGaussian(parms->mri_binary,
                                         parms->mri_binary_smooth,
                                         mri_kernel) ;
            mri_smooth = MRIconvolveGaussian(mri, mri_smooth, mri_kernel) ;
            MRIfree(&mri_kernel) ;
          }
          if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
          {
            MRIwrite(mri_smooth, "bin_smooth.mgz") ;
          }
        }

        if (DZERO(start_rms) && parms->regrid)
        {
          start_rms = GCAMcomputeRMS(gcam, mri, parms) ;
          start_t = parms->start_t ;
        }

        i = 0 ;
        do
        {
          if (((level != (parms->levels-1)) || (i > 0)) && parms->relabel)
          {
            if (mri)
            {
              GCAMcomputeLabels(mri, gcam) ;
            }
            if (parms->write_iterations != 0)
            {
              char fname[STRLEN] ;
              MRI  *mri_gca, *mri_tmp ;
              mri_gca = MRIclone(mri, NULL) ;
              GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
              if (mri_gca->nframes > 1)
              {
                printf("gcamorph: extracting %dth frame\n",
                       mri_gca->nframes-1) ;
                mri_tmp = MRIcopyFrame(mri_gca, NULL, mri_gca->nframes-1, 0) ;
                MRIfree(&mri_gca) ;
                mri_gca = mri_tmp ;
              }
              sprintf(fname, "%s_target%d", parms->base_name, level) ;
              MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
              sprintf(fname, "%s_target%d.mgz", parms->base_name, level) ;
              printf("writing target volume to %s...\n", fname) ;
              MRIwrite(mri_gca, fname) ;
              MRIfree(&mri_gca) ;
            }
          }

          GCAMregisterPipeline( gcam, mri, mri_smooth, parms,
                                &last_rms, &level_steps, i );


          if (parms->write_fname &&
              (((Gdiag & DIAG_WRITE)) || (getenv("WRITE_GCAM") != NULL)))
          {
            char path[STRLEN], fname_only[STRLEN] ;
            FileNameOnly(parms->write_fname, fname_only) ;
            FileNameRemoveExtension(fname_only, fname_only) ;
            FileNamePath(parms->write_fname, path) ;
            sprintf(fname, "%s/%s_level%d.m3z", path,fname_only,level);
            // strcpy(fname, parms->write_fname) ;
            printf("writing results of level to %s...\n", fname) ;
            //            GCAMvoxToRas(gcam) ;
            GCAMwrite(gcam, fname) ;
            //            GCAMrasToVox(gcam, mri) ;
          }
          parms->end_rms = rms = GCAMcomputeRMS(gcam, mri, parms) ;
          level_steps =
            parms->start_t - level_steps ;   /* # of steps taken in
                                                GCAMregisterLevel */
          if (level_steps == 0)
          {
            level_steps = 1 ;
          }
          pct_change = 100.0*(last_rms-rms)/(level_steps*last_rms) ;
#if 0
          printf("iter %d: last rms %2.3f, rms %2.3f, pct_change %2.3f/iter\n",
                 i+1, last_rms, rms, pct_change) ;
#endif
          i++ ;
#if 0
          if (parms->regrid == True)
          {
            GCAMregrid(gcam, mri, 0, parms, NULL) ;
          }
          else if (i >= 0)  /* don't bother iterating if not regridding */
#else
          if (i >= 0)  /* don't bother iterating if not regridding */
          {
            break ;
          }
#endif
          }
        while (pct_change > parms->tol) ;
        parms->sigma /= 4 ;
        if (parms->sigma < parms->min_sigma)
        {
          break ;
        }
        /*      GCAMremoveCompressedRegions(gcam, 0.1) ;*/
      }
      parms->navgs /= 4 ;
#if 0
      if (Gdiag & DIAG_WRITE &&
          parms->write_iterations > 0 && level > 0 && gcam_write_grad)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s_level%d.m3z",
                parms->base_name, parms->levels-level) ;
        GCAMvoxToRas(gcam) ;
        GCAMwrite(gcam, fname) ;
        GCAMrasToVox(gcam, mri) ;
      }
#endif
      if (parms->regrid > 1 && level == 0 /*&& passno >= parms->npasses-1*/)
      {
        int steps = parms->start_t - start_t ;
        pct_change = 100.0*(start_rms-rms)/(steps*start_rms) ;
        printf
        ("checking regridding - %d steps and rms %2.3f --> %2.3f (%2.3f)\n",
         steps, start_rms, rms, pct_change) ;
        if (pct_change > (parms->tol/10))
        {
          GCAMregrid(gcam, mri, 0, parms, NULL) ;
          parms->navgs = navgs ;
          level = parms->levels ;
          last_rms = GCAMcomputeRMS(gcam, mri, parms) ;
        }
      }
      if (parms->navgs < parms->min_avgs)
      {
        break ;
      }
    }
    parms->label_dist = label_dist ;
    // if not scaling smoothness and making multiple passes,
    // relax smoothness contstraint
    if ((passno < parms->npasses-1) && (parms->scale_smoothness == 0))
    {
      parms->l_smoothness /= 4 ;
      parms->l_distance /= 4 ;
      l_smooth /= 4 ;
      parms->l_distance /= 4 ;
      parms->l_elastic /= 4 ;  // relax elasticity slower as we want to find the solution with maximal elasticity
      l_elastic /= 4 ;
    }
    if (parms->navgs < parms->min_avgs)
    {
      break ;
    }
  }

  if (parms->noneg < 0 && gcam->neg >0)
    gcamRemoveNegativeNodes(gcam, mri, parms) ;
  if (mri_smooth)
  {
    MRIfree(&mri_smooth) ;
  }

  parms->sigma = base_sigma ;
  if (parms->log_fp)
  {
    fclose(parms->log_fp) ;
    parms->log_fp = NULL ;
  }

  parms->relabel = orig_relabel ;
  parms->l_smoothness = l_orig_smooth ;
  parms->l_elastic = l_orig_elastic ;
  parms->navgs = navgs ;
  parms->dt = orig_dt ;
  if (parms->mri_dist_map)
    MRIfree(&parms->mri_dist_map) ;
  if (parms->mri_atlas_dist_map)
    MRIfree(&parms->mri_atlas_dist_map) ;

  return(NO_ERROR) ;
}

GCA_MORPH *
GCAMread(const char *fname)
{
  GCA_MORPH       *gcam ;
  znzFile file;
  int             x, y, z, width, height, depth ;
  GCA_MORPH_NODE  *gcamn ;
  float           version ;
  int             tag;
  int gzipped = 0;

  if (!fio_FileExistsReadable(fname))
  {
    printf("ERROR: cannot find or read %s\n",fname);
    return(NULL);
  }

  if (strstr(fname, ".m3z"))
  {
    gzipped = 1;
  }

  file = znzopen(fname, "rb", gzipped);

  if (znz_isnull(file))
  {
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "GCAMread(%s): could not open file", fname)) ;
  }

  version = znzreadFloat(file) ;
  if (version != GCAM_VERSION)
  {
    znzclose(file);
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "GCAMread(%s): invalid version # %2.3f\n", fname, version)) ;
  }
  width  = znzreadInt(file) ;
  height = znzreadInt(file) ;
  depth  = znzreadInt(file) ;
  gcam   = GCAMalloc(width, height, depth) ;

  gcam->spacing = znzreadInt(file) ;
  gcam->exp_k   = znzreadFloat(file) ;

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        gcamn->origx = znzreadFloat(file) ;
        gcamn->origy = znzreadFloat(file) ;
        gcamn->origz = znzreadFloat(file) ;

        gcamn->x = znzreadFloat(file) ;
        gcamn->y = znzreadFloat(file) ;
        gcamn->z = znzreadFloat(file) ;

        gcamn->xn = znzreadInt(file) ;
        gcamn->yn = znzreadInt(file) ;
        gcamn->zn = znzreadInt(file) ;

        // if all the positions are zero, then this is not a valid point
        // mark invalid = 1
        if (FZERO(gcamn->origx) && FZERO(gcamn->origy) && FZERO(gcamn->origz)
            && FZERO(gcamn->x) && FZERO(gcamn->y) && FZERO(gcamn->z))
        {
          gcamn->invalid = GCAM_POSITION_INVALID ;
        }
        else
        {
	  if (x == 0 || x == width-1 ||
	      y == 0 || y == height-1 ||
	      z == 0 || z == depth-1)
	    gcamn->invalid = GCAM_AREA_INVALID ;
	  else
	    gcamn->invalid = GCAM_VALID ;
        }
      }
    }
  }
  gcam->det = 1 ;
  gcam->image.valid = 0; // make src invalid
  gcam->atlas.valid = 0; // makd dst invalid
  while (znzreadIntEx(&tag, file))
  {
    switch (tag)
    {
    case TAG_GCAMORPH_LABELS:
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      {
        printf("reading labels out of gcam file...\n") ;
      }
      gcam->status = GCAM_LABELED ;
      for (x = 0 ; x < width ; x++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            gcamn = &gcam->nodes[x][y][z] ;
            gcamn->label = znzreadInt(file) ;
            if (gcamn->label != 0)
            {
              DiagBreak() ;
            }
          }
        }
      }
      break ;
    case TAG_GCAMORPH_GEOM:
      GCAMreadGeom(gcam, file);
      if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
      {
        fprintf(stderr,"GCAMORPH_GEOM tag found.  "
                "Reading src and dst information.\n");
        fprintf(stderr, "src geometry:\n");
        writeVolGeom(stderr, &gcam->image);
        fprintf(stderr, "dst geometry:\n");
        writeVolGeom(stderr, &gcam->atlas);
      }
      break ;
    case TAG_GCAMORPH_TYPE:
      gcam->type = znzreadInt(file) ;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      {
        printf("gcam->type = %s\n", gcam->type == GCAM_VOX ? "vox" : "ras") ;
      }
      break ;
    case TAG_MGH_XFORM:
      //gcam->m_affine = MatrixAsciiReadFrom(fp, NULL) ;
      gcam->m_affine = znzReadMatrix(file) ;
      gcam->det = MatrixDeterminant(gcam->m_affine) ;
      break ;
    }
  }

  znzclose(file);

  GCAMcomputeOriginalProperties(gcam) ;
  gcamComputeMetricProperties(gcam) ;
  GCAMcopyNodePositions(gcam, ORIGINAL_POSITIONS, SAVED_ORIGINAL_POSITIONS) ;
  return(gcam) ;
}


GCA_MORPH *
GCAMalloc( const int width, const int height, const int depth )
{
  GCA_MORPH  *gcam ;
  int        x, y, z ;

  gcam = (GCA_MORPH *)calloc(1, sizeof(GCA_MORPH)) ;
  if (!gcam)
    ErrorExit
    (ERROR_NOMEMORY, "GCAMalloc: could not allocate GCA_MORPH struct") ;

  gcam->width =  width  ;
  gcam->height = height ;
  gcam->depth =  depth  ;
  gcam->spacing = 1 ;  // may be changed by the user later

  gcam->nodes = (GCA_MORPH_NODE ***)calloc(width, sizeof(GCA_MORPH_NODE **)) ;
  if (!gcam->nodes)
  {
    ErrorExit(ERROR_NOMEMORY, "GCAMalloc: could not allocate nodes") ;
  }

  for (x = 0 ; x < gcam->width ; x++)
  {
    gcam->nodes[x] =
      (GCA_MORPH_NODE **)calloc(gcam->height, sizeof(GCA_MORPH_NODE *)) ;
    if (!gcam->nodes[x])
    {
      ErrorExit(ERROR_NOMEMORY, "GCAMalloc: could not allocate %dth **",x) ;
    }

#define ELECTRIC_FENCE 0
#if ELECTRIC_FENCE
    {
      GCA_MORPH_NODE *buf ;

      buf =
        (GCA_MORPH_NODE *)calloc(gcam->depth*gcam->height,
                                 sizeof(GCA_MORPH_NODE)) ;
      if (buf == NULL)
        ErrorExit
        (ERROR_NO_MEMORY,
         "GCAMalloc(%d, %d, %d): could not allocate "
         "%d bytes for %dth bslice\n",
         height, width, depth,
         (gcam->depth*gcam->height*sizeof(GCA_MORPH_NODE)), x);
      for (y = 0 ; y < gcam->height ; y++)
      {
        gcam->nodes[x][y] = buf+(y*gcam->depth) ;
      }
    }
#else
    for (y = 0 ; y < gcam->height ; y++)
    {
      gcam->nodes[x][y] = (GCA_MORPH_NODE *)calloc(gcam->depth,
                          sizeof(GCA_MORPH_NODE)) ;
      if (!gcam->nodes[x][y])
        ErrorExit(ERROR_NOMEMORY,
                  "GCAMalloc: could not allocate %d,%dth *",x,y);
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcam->nodes[x][y][z].origx = x ;
        gcam->nodes[x][y][z].origy = y ;
        gcam->nodes[x][y][z].origz = z ;
        gcam->nodes[x][y][z].x = x ;
        gcam->nodes[x][y][z].y = y ;
        gcam->nodes[x][y][z].z = z ;
      }
    }
#endif
  }
  initVolGeom(&gcam->image);
  initVolGeom(&gcam->atlas);
  return(gcam) ;
}

int
GCAMinitVolGeom(GCAM *gcam, MRI *mri_image, MRI *mri_atlas)
{
  if (gcam->type == GCAM_VOX) /* in some other voxel coords */
  {
    GCAMrasToVox(gcam, mri_image) ;
  }

  getVolGeom(mri_image, &gcam->image) ;
  getVolGeom(mri_atlas, &gcam->atlas) ;
  if (gcam->gca)
  {
    GCAreinit(mri_atlas, gcam->gca) ;
  }
  return(NO_ERROR) ;
}

int
GCAMpreserveLabelMetricProperties(GCA_MORPH *gcam, LABEL *area, MRI *mri)
{
  int              n, xv, yv, zv ;
  GCA_MORPH_NODE   *gcamn ;
  double           xvd, yvd, zvd ;
  LABEL_VERTEX     *lv ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    MRIsurfaceRASToVoxel(mri, lv->x, lv->y, lv->z, &xvd, &yvd, &zvd) ;
    xv = nint(xvd) ;
    yv = nint(yvd) ;
    zv = nint(zvd) ;
    if (xv == Gx && yv == Gy && zv == Gz)
    {
      DiagBreak() ;
    }
    gcamn = &gcam->nodes[xv][yv][zv] ;
    gcamn->status &= ~(GCAM_IGNORE_DISTANCES) ;
    gcamn->status |= GCAM_DISCOUNT_LIKELIHOOD ;
  }
  return(NO_ERROR) ;
}
int
GCAMinit(GCA_MORPH *gcam,
         MRI *mri_image, GCA *gca,
         TRANSFORM *transform,
         int relabel)
{
  GCA_MORPH_NODE  *gcamn ;
  GC1D            *gc ;
  GCA_PRIOR       *gcap ;
  int             invalid = 0, x, y, z, width, height, depth, n, max_n,
                  max_label, label, xform_allocated = 0 ;
  float           max_p, ox, oy, oz ;

  if (!mri_image)
  {
    ErrorExit(ERROR_BADPARM,
              "GCAMinit() must be called with valid mri_image.\n");
  }

  if (transform == NULL)
  {
    xform_allocated = 1 ;
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, mri_image) ;
  }

  gcam->ninputs = mri_image->nframes ;

  // save geometry information
  getVolGeom(mri_image, &gcam->image);
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;

  TransformInvert(transform, mri_image) ;
  if (gca)
  {
    // check consistency
    GCAsetVolGeom(gca, &gcam->atlas); // fname is still "unknown"
    if (width != gca->prior_width
        || height != gca->prior_height
        || depth != gca->prior_depth)
    {
      fprintf
      (stderr,
       "GCA_MORPH (%d, %d, %d) must agree with GCA prior (%d, %d, %d)\n",
       gcam->width, gcam->height, gcam->depth, gca->prior_width, \
       gca->prior_height, gca->prior_depth);
      ErrorExit(ERROR_BADPARM, "Exiting ....\n");
    }
    gcam->gca = gca ;
    gcam->spacing = gca->prior_spacing ;


    // use gca information
    for (x = 0 ; x < width ; x++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (z = 0 ; z < depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          gcamn = &gcam->nodes[x][y][z] ;
          gcap = &gca->priors[x][y][z] ;
          max_p = 0 ;
          max_n = -1 ;
          max_label = 0 ;

          // find the label which has the max p
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            label = gcap->labels[n] ;   // get prior label
            if (label == Gdiag_no)
            {
              DiagBreak() ;
            }
            if (label >= MAX_CMA_LABEL)
            {
              printf("invalid label %d at (%d, %d, %d) in prior volume\n",
                     label, x, y, z);
            }
            if (gcap->priors[n] >= max_p) // update the max_p and max_label
            {
              max_n = n ;
              max_p = gcap->priors[n] ;
              max_label = gcap->labels[n] ;
            }
          }
          gcamn->xn = x ;
          gcamn->yn = y ;
          gcamn->zn = z ;
          // here mri info is used
          if (!GCApriorToSourceVoxelFloat(gca, mri_image, transform, x, y, z,
                                          &ox, &oy, &oz) || 1)  // ALWAYS!!!
          {
	    if (x > 0 && x < width-1 &&
		y > 0 && y < height-1 &&
		z > 0 && z < depth-1)
	      gcamn->invalid = GCAM_VALID;
            // if inside the volume, then
            gcamn->x = gcamn->origx = ox ;
            gcamn->y = gcamn->origy = oy ;
            gcamn->z = gcamn->origz = oz ;
#if 1
            gcamn->label = max_label ;
            gcamn->n = max_n ;
            gcamn->prior = max_p ;
            gc = GCAfindPriorGC(gca, x, y, z, max_label) ;
            if (gc == NULL)
            {
              int xn, yn, zn ;
              GCApriorToNode(gca, x, y, z, &xn, &yn, &zn);
              gc = GCAfindClosestValidGC(gca, xn, yn, zn, max_label,0);
            }
            // gc can be NULL
            if (gc == NULL)
            {
              DiagBreak() ;
            }
            gcamn->gc = gc ;
            gcamn->log_p = 0 ;
#endif
            if (x == Gx && y == Gy && z == Gz)
            {
              printf("node(%d,%d,%d) --> MRI (%2.1f, %2.1f, %2.1f),"
                     " label %s (%d) = %2.1f\n",
                     x, y, z, ox, oy, oz, cma_label_to_name(max_label),
                     max_label,
                     gcamn->gc ? gcamn->gc->means[0] : -1) ;
              DiagBreak() ;
#if 0
              Gvx = nint(ox) ;
              Gvy = nint(oy) ;
              Gvz = nint(oz) ;
              /* for writing out image views */
#endif
            }
          }// !GCA
          else // went outside the volume but we must initialize
          {
            gcamn->invalid = GCAM_POSITION_INVALID; // mark invalid
            invalid++ ;
            gcamn->x = gcamn->origx = 0 ;
            gcamn->y = gcamn->origy = 0 ;
            gcamn->z = gcamn->origz = 0 ;
            gcamn->label = 0 ; // unknown
            gcamn->n = 0;
            gcamn->prior = 1.;
            gcamn->gc = 0 ;
            gcamn->log_p = 0 ;
          }
        }
      }
    }

    if (relabel)
    {
      GCAMcomputeLabels(mri_image, gcam) ;
    }
    else
    {
      GCAMcomputeMaxPriorLabels(gcam) ;
    }
  }
  else   /* no gca specified */
  {
    VECTOR *v1, *v2 ;
    MATRIX *m_L ;

    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v1, 4, 1) = 1.0 ;
    *MATRIX_RELT(v2, 4, 1) = 1.0 ;
    m_L = ((LTA *)(transform->xform))->xforms[0].m_L;
    for (x = 0 ; x < width ; x++)
    {
      V3_X(v1) = x ;
      for (y = 0 ; y < height ; y++)
      {
        V3_Y(v1) = y ;
        for (z = 0 ; z < depth ; z++)
        {
          V3_Z(v1) = z ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          gcamn = &gcam->nodes[x][y][z] ;

          gcamn->xn = x ;
          gcamn->yn = y ;
          gcamn->zn = z ;

          MatrixMultiply(m_L, v1, v2) ;
          ox = V3_X(v2) ;
          oy = V3_Y(v2) ;
          oz = V3_Z(v2) ;
          // here mri info is used
          if ((MRIindexNotInVolume(mri_image, ox, oy, oz) == 0) || 1) // ALWAYS!!
          {
            gcamn->invalid = GCAM_VALID;
            // if inside the volume, then
            gcamn->x = gcamn->origx = ox ;
            gcamn->y = gcamn->origy = oy ;
            gcamn->z = gcamn->origz = oz ;
            gcamn->gc = NULL ;
            gcamn->log_p = 0 ;

            if (x == Gx && y == Gy && z == Gz)
            {
              printf("node(%d,%d,%d) --> MRI (%2.1f, %2.1f, %2.1f)\n",
                     x, y, z, ox, oy, oz) ;
              DiagBreak() ;
            }
          }// !GCA
          else // went outside the volume but we must initialize NEVER!!
          {
            if (gcamn->label > 0)
            {
              DiagBreak() ;
            }
            gcamn->invalid = GCAM_POSITION_INVALID; // mark invalid
            gcamn->x = gcamn->origx = 0 ;
            gcamn->y = gcamn->origy = 0 ;
            gcamn->z = gcamn->origz = 0 ;
            gcamn->label = 0 ; // unknown
            gcamn->n = 0;
            gcamn->prior = 1.;
            gcamn->gc = 0 ;
            invalid++ ;
            gcamn->log_p = 0 ;
          }
        }
      }
    }

    MatrixFree(&v1) ;
    MatrixFree(&v2) ;
    gcam->spacing = 1 ;
  }

  if (transform->type != MORPH_3D_TYPE)
  {
    LTA *lta = (LTA *)transform->xform  ;
    gcam->m_affine = MatrixCopy(lta->xforms[0].m_L, NULL) ;
    gcam->det = MatrixDeterminant(gcam->m_affine) ;
    fprintf(stderr,"det(m_affine) = %2.2f (predicted orig area = %2.1f)\n",
	    gcam->det, gcam->spacing*gcam->spacing*gcam->spacing/gcam->det) ;
  }
  gcamComputeMetricProperties(gcam) ;
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
      {
        gcam->nodes[x][y][z].orig_area = gcam->nodes[x][y][z].area ;
        gcam->nodes[x][y][z].orig_area1 = gcam->nodes[x][y][z].area1 ;
        gcam->nodes[x][y][z].orig_area2 = gcam->nodes[x][y][z].area2 ;
#if 0
        if (((gcam->nodes[x][y][z].orig_area <= 0) ||
             (gcam->nodes[x][y][z].orig_area1 <= 0) ||
             (gcam->nodes[x][y][z].orig_area2 <= 0) ||
             (gcam->nodes[x][y][z].area <= 0)) &&
            (gcam->nodes[x][y][z].invalid == GCAM_VALID))
#endif
          if (((gcam->nodes[x][y][z].orig_area1 <= 0) &&
               (gcam->nodes[x][y][z].orig_area2 <= 0)) &&
              (gcam->nodes[x][y][z].invalid == GCAM_VALID))
          {
            gcam->nodes[x][y][z].invalid = GCAM_AREA_INVALID ;
            invalid++ ;
          }
      }

#if 0
  gcamComputeMetricProperties(gcam) ;
  while (gcam->neg > 0)
  {
    MRI *mri_warp ;
    mri_warp = GCAMwriteWarpToMRI(gcam, NULL) ;
    GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri_warp) ;
    MRIfree(&mri_warp) ;
  }
  GCAMstoreMetricProperties(gcam) ;
#endif

  if (gcam->mri_xind)
  {
    MRIfree(&gcam->mri_xind) ;
    MRIfree(&gcam->mri_yind) ;
    MRIfree(&gcam->mri_zind) ;
  }
  gcam->mri_xind = MRIalloc(mri_image->width,
                            mri_image->height,
                            mri_image->depth,
                            MRI_FLOAT) ;
  gcam->mri_yind = MRIalloc(mri_image->width,
                            mri_image->height,
                            mri_image->depth,
                            MRI_FLOAT) ;
  gcam->mri_zind = MRIalloc(mri_image->width,
                            mri_image->height,
                            mri_image->depth,
                            MRI_FLOAT) ;
  for (x = 0 ; x < mri_image->width ; x++)
    for (y = 0 ; y < mri_image->height ; y++)
      for (z = 0 ; z < mri_image->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        TransformSample(transform, x, y, z, &ox, &oy, &oz) ;
        MRIsetVoxVal(gcam->mri_xind, x, y,z, 0, ox/gcam->spacing) ;
        MRIsetVoxVal(gcam->mri_yind, x, y,z, 0, oy/gcam->spacing) ;
        MRIsetVoxVal(gcam->mri_zind, x, y,z, 0, oz/gcam->spacing) ;
      }

  GCAMcopyNodePositions(gcam, ORIGINAL_POSITIONS, SAVED_ORIGINAL_POSITIONS) ;
  GCAMcomputeOriginalProperties(gcam) ;
  gcamComputeMetricProperties(gcam) ;
  gcam->type = GCAM_VOX ;
  if (xform_allocated)
  {
    TransformFree(&transform) ;
  }
  return(NO_ERROR) ;
}

////////////////////////////////////////////////
// user is responsible for freeing gca inside
// the gcam.
int
GCAMfree(GCA_MORPH **pgcam)
{
  GCA_MORPH *gcam ;

  gcam = *pgcam ;
  *pgcam = NULL ;

  GCAMfreeContents(gcam) ;
  free(gcam) ;
  return(NO_ERROR) ;
}
////////////////////////////////////////////////
// user is responsible for freeing gca inside
// the gcam.
int
GCAMfreeContents(GCA_MORPH *gcam)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  GCAMfreeInverse(gcam) ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->gc &&
            gcam->gca == NULL)  // don't free the gcs if
          // they are part of the gca
        {
          free_gcs(gcamn->gc, 1, gcam->ninputs) ;
        }
      }
      free(gcam->nodes[x][y]) ;
    }
    free(gcam->nodes[x]) ;
  }
  free(gcam->nodes) ;
  return(NO_ERROR) ;
}


/*
  Note:  d [ I(r)' C I(r), r] = delI * C * I(r)
  where delI(r) = 3xn, C=nxn, and I(r)=nx1
*/

#define GCAM_LLT_OUTPUT 0

int
gcamLogLikelihoodTerm( GCA_MORPH *gcam,
                       const MRI *mri,
                       const MRI *mri_smooth,
                       double l_log_likelihood )
{
#ifdef GCAM_LL_TERM_GPU
  if( DZERO(l_log_likelihood) )
  {
    return( NO_ERROR );
  }

  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamLogLikelihoodTermGPU( gcam, mri, mri_smooth, l_log_likelihood );
#else
  int             x=0, y=0, z=0, n=0 /*,label*/ ;
  int             i;
  int             nthreads=1, tid=0;
  double          dx=0.0, dy=0.0, dz=0.0, norm=0.0;
  float           vals[_MAX_FS_THREADS][MAX_GCA_INPUTS] ;
  GCA_MORPH_NODE  *gcamn=NULL ;
  MATRIX          *m_delI[_MAX_FS_THREADS], *m_inv_cov[_MAX_FS_THREADS] ;
  VECTOR          *v_means[_MAX_FS_THREADS], *v_grad[_MAX_FS_THREADS] ;

  if (DZERO(l_log_likelihood))
  {
    return(NO_ERROR) ;
  }

#if GCAM_LLT_OUTPUT
  const unsigned int gcamLLToutputFreq = 10;
  static unsigned int nCalls = 0;

  if( (nCalls%gcamLLToutputFreq)==0 )
  {
    char fname[STRLEN];
    const unsigned int nOut = nCalls/gcamLLToutputFreq;

    snprintf( fname, STRLEN-1, "gcamLLTinput%04u", nOut );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mriLLTinput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri, fname );

    snprintf( fname, STRLEN-1, "mrismoothLLTinput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri_smooth, fname );
  }

  nCalls++;
#endif

#ifdef HAVE_OPENMP
  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
#else
  nthreads = 1;
#endif

  for (i=0; i<nthreads; i++)
  {
    m_delI[i] = MatrixAlloc(3, gcam->ninputs, MATRIX_REAL) ;
    m_inv_cov[i] = MatrixAlloc(gcam->ninputs, gcam->ninputs, MATRIX_REAL) ;
    v_means[i] = VectorAlloc(gcam->ninputs, 1) ;
    v_grad[i] = VectorAlloc(3, MATRIX_REAL) ;
  }

#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate(tid,y,z,gcamn,n,norm,dx,dy,dz,vals,m_delI,m_inv_cov,v_means,v_grad) shared(gcam,mri,Gx,Gy,Gz,Gvx,Gvy,Gvz) schedule(static,1)
#endif

  for (x = 0 ; x < gcam->width ; x++)
  {
#ifdef HAVE_OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {

        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        if (fabs(gcamn->x-Gvx)<1 &&
            fabs(gcamn->y-Gvy)<1  &&
            fabs(gcamn->z-Gvz)<1)
        {
          DiagBreak() ;
        }

        if( gcamn->status &
            (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD) )
        {
          continue;
        }

        /* don't use unkown nodes unless they border
           something that's not unknown */
        if( IS_UNKNOWN(gcamn->label) &&
            different_neighbor_labels(gcam, x,y,z,1) == 0 )
        {
          continue;
        }

        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals[tid], gcam->ninputs);

        if (!gcamn->gc)
        {
          MatrixClear(v_means[tid]) ;
          MatrixIdentity(gcam->ninputs, m_inv_cov[tid]) ;
          MatrixScalarMul
          (m_inv_cov[tid], 1.0/(MIN_VAR), m_inv_cov[tid]) ; /* variance=4 is min */
        }
        else
        {
#if 0
          if (parms->relabel)
          {
            label =
              GCAcomputeMAPlabelAtLocation
              (gcam->gca, x,y,z,vals,&n,&gcamn->log_p);
            if (label == gcamn->label)  /* already correct label -
                                           don't move anywhere */
            {
              continue ;
            }
          }
#endif
          load_mean_vector(gcamn->gc, v_means[tid], gcam->ninputs) ;
          load_inverse_covariance_matrix(gcamn->gc, m_inv_cov[tid], gcam->ninputs) ;
        }

        for (n = 0 ; n < gcam->ninputs ; n++)
        {
          MRIsampleVolumeGradientFrame(mri_smooth,
                                       gcamn->x, gcamn->y, gcamn->z,
                                       &dx, &dy, &dz, n) ;
          norm = sqrt(dx*dx+dy*dy+dz*dz) ;
          if (!FZERO(norm))  /* don't worry about magnitude of gradient */
          {
            dx /= norm ;
            dy /= norm ;
            dz /= norm ;
          }
          *MATRIX_RELT(m_delI[tid], 1, n+1) = dx ;
          *MATRIX_RELT(m_delI[tid], 2, n+1) = dy ;
          *MATRIX_RELT(m_delI[tid], 3, n+1) = dz ;
          VECTOR_ELT(v_means[tid], n+1) -= vals[tid][n] ;
#define MAX_ERROR 1000
          if (fabs(VECTOR_ELT(v_means[tid], n+1)) > MAX_ERROR)
            VECTOR_ELT(v_means[tid], n+1) =
              MAX_ERROR * FSIGN(VECTOR_ELT(v_means[tid], n+1)) ;
        }

        MatrixMultiply(m_inv_cov[tid], v_means[tid], v_means[tid]) ;

        if (IS_UNKNOWN(gcamn->label))
        {
          if (zero_vals(vals[tid], gcam->ninputs))
          {
            if (Gx == x && Gy == y && Gz == z)
              printf("discounting unknown label at (%d, %d, %d) "
                     "due to skull strip difference\n",
                     x, y, z) ;
            /* probably difference in skull stripping (vessels present or
               absent) - don't let it dominate */
            if (VECTOR_ELT(v_means[tid],1) > .5)  /* don't let it be more
                                                than 1/2 stds away */
            {
              VECTOR_ELT(v_means[tid],1) = .5 ;
            }
          }
#if 0
          else if (VECTOR_ELT(v_means,1) > 2)  /* don't let it be more
                                               than 2 stds away */
          {
            VECTOR_ELT(v_means,1) = 2 ;
          }
#endif
        }
        MatrixMultiply(m_delI[tid], v_means[tid], v_grad[tid]) ;

        gcamn->dx += l_log_likelihood*V3_X(v_grad[tid]) ;
        gcamn->dy += l_log_likelihood*V3_Y(v_grad[tid]) ;
        gcamn->dz += l_log_likelihood*V3_Z(v_grad[tid]) ;

        if (x == Gx && y == Gy && z == Gz)
        {
          printf
          ("ll_like: node(%d,%d,%d): dI=(%2.1f,%2.1f,%2.1f), "
           "grad=(%2.2f,%2.2f,%2.2f), "
           "node %2.2f+-%2.2f, MRI=%2.1f\n",
           x, y, z, dx, dy, dz, gcamn->dx, gcamn->dy, gcamn->dz,
           gcamn->gc ? gcamn->gc->means[0] : 0.0,
           gcamn->gc ? sqrt(covariance_determinant
                            (gcamn->gc,
                             gcam->ninputs)) : 0.0, vals[tid][0]) ;
        }
      }
    }
  }

  for (i=0; i<nthreads; i++)
  {
    MatrixFree(&m_delI[i]) ;
    MatrixFree(&m_inv_cov[i]) ;
    VectorFree(&v_means[i]) ;
    VectorFree(&v_grad[i]) ;
  }
#endif
  return(NO_ERROR) ;
}

/*
 */
int
gcamLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth,
                   double l_likelihood, \
                   GCA_MORPH_PARMS *parms)
{
  int             x, y, z, len, xi, yi, zi, x0, y0, z0, half_len ;
  double            dx, dy, dz ;
  float           val, mean, var ;
  GCA_MORPH_NODE  *gcamn ;
  MRI             *mri_nbhd, *mri_kernel ;

  if (DZERO(l_likelihood))
  {
    return(NO_ERROR) ;
  }

  len = (int)nint(4.0f * parms->sigma)+1 ;
  if (ISEVEN(len))   /* ensure it's odd */
  {
    len++ ;
  }
  half_len = (len-1)/2 ;

  mri = MRISeqchangeType(mri, MRI_FLOAT, 0, 0, 1) ;
  mri_nbhd = MRIallocSequence(len,len,len, MRI_FLOAT, gcam->ninputs) ;
  // must copy header info.  we are missing here...........
  mri_kernel = MRIgaussian1d(parms->sigma, 0) ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        if (fabs(gcamn->x-Gvx)<1 &&
            fabs(gcamn->y-Gvy)<1  &&
            fabs(gcamn->z-Gvz)<1)
        {
          DiagBreak() ;
        }
        if (gcamn->status &
            (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD) ||
            (gcamn->gc == NULL))
        {
          continue ;
        }

        x0 = nint(gcamn->x) ;
        y0 = nint(gcamn->y) ;
        z0 = nint(gcamn->z) ;
        if ((x0+half_len >= mri->width) ||
            (y0+half_len >= mri->height) ||
            (z0+half_len >= mri->depth))
        {
          continue ;
        }
        MRIextractInto(mri, mri_nbhd,
                       x0-half_len, y0-half_len, z0-half_len,
                       len, len, len, 0, 0, 0) ;
        if (gcam->ninputs == 1)
        {
          var = gcamn->gc->covars[0] ;
          mean = gcamn->gc->means[0] ;

          for (xi = 0 ;  xi < mri_nbhd->width ; xi++)
          {
            for (yi = 0 ;  yi < mri_nbhd->height ; yi++)
            {
              for (zi = 0 ;  zi < mri_nbhd->depth ; zi++)
              {
                val = MRIFvox(mri_nbhd, xi, yi, zi) ;
                val = (val - mean) * (val-mean) / var ;
                MRIFvox(mri_nbhd, xi, yi, zi) = sqrt(val) ;
              }
            }
          }
        }
        else  /* haven't written the multi-input case yet */
        {
          ErrorExit
          (ERROR_UNSUPPORTED, "vector-based likelihood not written yet") ;
        }
        MRIconvolveGaussian(mri_nbhd, mri_nbhd, mri_kernel) ;
        dx = -MRIgetVoxDx(mri_nbhd, half_len, half_len, half_len, 0) ;
        dy = -MRIgetVoxDy(mri_nbhd, half_len, half_len, half_len, 0) ;
        dz = -MRIgetVoxDz(mri_nbhd, half_len, half_len, half_len, 0) ;
        if (gcamn->status & GCAM_DISCOUNT_LIKELIHOOD)
        {
          dx *= .1 ;
          dy *= .1 ;
          dz *= .1 ;
        }

        gcamn->dx += l_likelihood*dx ;
        gcamn->dy += l_likelihood*dy ;
        gcamn->dz += l_likelihood*dz ;

        if (x == Gx && y == Gy && z == Gz)
          printf("l_like: node(%d,%d,%d): dp=(%2.2f,%2.2f,%2.2f), "
                 "node %2.2f+-%2.2f\n",
                 x, y, z, gcamn->dx, gcamn->dy, gcamn->dz,
                 gcamn->gc ? gcamn->gc->means[0] : 0.0,
                 gcamn->gc ? sqrt(covariance_determinant
                                  (gcamn->gc, gcam->ninputs)) : 0.0) ;
      }

  MRIfree(&mri_kernel) ;
  MRIfree(&mri_nbhd) ;
  MRIfree(&mri) ;
  return(NO_ERROR) ;
}

#define DEBUG_LL_SSE 0
#if DEBUG_LL_SSE
static float ***last_sse = NULL;
#endif

#define GCAM_LLENERGY_OUTPUT 0

double
gcamLogLikelihoodEnergy( const GCA_MORPH *gcam, MRI *mri)
{
  /*
    NOTE: The only reason *mri is not declared 'const' is because
    I can't figure out how to make MRIwrite accept a const MRI
    structure.
    One might be perturbed that writing an MRI out to disc
    changes it....
  */
#if GCAM_LLENERGY_OUTPUT
  const unsigned int gcamLLEoutputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%gcamLLEoutputFreq)==0 )
  {
    char fname[STRLEN];
    snprintf( fname, STRLEN-1, "gcamLLEinput%04u", nCalls/gcamLLEoutputFreq );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );
    snprintf( fname, STRLEN-1, "mriLLEinput%04u.mgz", nCalls/gcamLLEoutputFreq );
    MRIwrite( mri, fname );
  }
#endif

  double sse = 0.0;

#ifndef GCAM_LLENERGY_GPU
  double error=0.0 ;
  int x=0, y=0, z=0;
#endif

#if DEBUG_LL_SSE
  int max_x, max_y, max_z ;
  max_x = max_y = max_z = 0 ;
#endif


#if DEBUG_LL_SSE
  double max_increase = 0, increase ;
  if (last_sse == NULL)
  {
    last_sse = (float ***)calloc(gcam->width, sizeof(float*)) ;
    for (x = 0 ; x < gcam->width ; x++)
    {
      last_sse[x] = (float **)calloc(gcam->height, sizeof(float*)) ;
      for (y = 0 ; y < gcam->height ; y++)
      {
        last_sse[x][y] = (float *)calloc(gcam->depth, sizeof(float)) ;
      }
    }
  }
#endif

#ifdef GCAM_LLENERGY_GPU
  if( gcam->ninputs != 0 )
  {
    printf( "%s: ninputs = %i\n", __FUNCTION__, gcam->ninputs );
  }
#endif

#ifdef GCAM_LLENERGY_GPU
#if SHOW_EXEC_LOC
  printf( "%s: CUDA call\n", __FUNCTION__ );
#endif
  sse = gcamLogLikelihoodEnergyGPU( gcam, mri );
#else
#if SHOW_EXEC_LOC
  printf( "%s: CPU call\n", __FUNCTION__ );
#endif
#ifdef HAVE_OPENMP
  #pragma omp parallel for reduction (+:sse) firstprivate(y,z,error) shared(gcam,mri) schedule(static,1)
#endif
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      float vals[MAX_GCA_INPUTS];
      for (z = 0 ; z < gcam->depth ; z++)
      {

        // Debugging breakpoint
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak();
        }

        // Shorthand way of accessing current node
        const GCA_MORPH_NODE* /* const */ gcamn = &gcam->nodes[x][y][z] ;

        // Don't operate on invalid nodes
        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }
        check_gcam(gcam) ;

        // Check for ignore
        if ( gcamn->status &
             (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD) )
        {
          continue ;
        }

        /* don't use unkown nodes unless they border
           something that's not unknown */
        if (IS_UNKNOWN(gcamn->label) &&
            (different_neighbor_labels(gcam, x,y,z,1) == 0) )
        {

          continue ;
        }

        check_gcam(gcam) ;

        // Load up the MRI values (which will do trilinear interpolation)
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->ninputs);
        check_gcam(gcam);

        // Compute 'error' for this node
        if( gcamn->gc )
        {
          error = GCAmahDist(gcamn->gc, vals, gcam->ninputs)
                  + log(covariance_determinant(gcamn->gc, gcam->ninputs));
        }
        else
        {
          int n ;
          // Note that the for loop sets error=0 on the first iteration
          for (n = 0, error = 0.0 ; n < gcam->ninputs ; n++)
          {
            error += (vals[n]*vals[n]/MIN_VAR) ;
          }
        }

        check_gcam(gcam) ;

        // Random output
        if (x == Gx && y == Gy && z == Gz)
          printf("E_like: node(%d,%d,%d) -> "
                 "(%2.1f,%2.1f,%2.1f), target=%2.1f+-%2.1f, val=%2.1f\n",
                 x, y, z, gcamn->x, gcamn->y, gcamn->z,
                 gcamn->gc ? gcamn->gc->means[0] : 0.0,
                 gcamn->gc ? sqrt(covariance_determinant
                                  (gcamn->gc, gcam->ninputs)) : 0.0, vals[0]) ;

        check_gcam(gcam) ;
#if DEBUG_LL_SSE
        if (last_sse[x][y][z] < (.9*error) && !FZERO(last_sse[x][y][z]))
        {
          DiagBreak() ;
          increase = error - last_sse[x][y][z] ;
          if (increase > max_increase)
          {
            max_increase = increase ;
            max_x = x ;
            max_y = y ;
            max_z = z ;
          }
        }
        last_sse[x][y][z] = (error) ;
#endif

        // Accumulate onto main variable
        sse += (error) ;

        // Enable debugging breakpoint if not finite
        if (!finitep(sse))
        {
          DiagBreak() ;
        }
      }
    }
  }
#endif

#if DEBUG_LL_SSE
  if (Gdiag & DIAG_SHOW)
    printf("max increase %2.2f at (%d, %d, %d)\n",
           max_increase, max_x, max_y, max_z) ;
#endif


#if GCAM_LLENERGY_OUTPUT
  nCalls++;
#endif

  return(sse) ;
}


int
gcamDistanceTerm(GCA_MORPH *gcam, MRI *mri, double l_distance)
{
  double          dx, dy, dz, error, d0, d, xdelta, ydelta, zdelta ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  if (DZERO(l_distance))
  {
    return(NO_ERROR) ;
  }
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if ((gcamn->invalid == GCAM_POSITION_INVALID) ||
            (gcamn->status & GCAM_IGNORE_DISTANCES))
        {
          continue;
        }

        num = 0 ;
        xdelta = ydelta = zdelta = 0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          xn = MAX(0,xn) ;
          xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            yn = MAX(0,yn) ;
            yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!xk && !yk && !zk)
              {
                continue ;
              }
              zn = z+zk ;
              zn = MAX(0,zn) ;
              zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;

              if (gcamn_nbr->invalid/* == GCAM_POSITION_INVALID*/)
              {
                continue;
              }
#if 0
              if (gcamn_nbr->label != gcamn->label)
              {
                continue ;
              }
#endif
              dx = gcamn->origx - gcamn_nbr->origx ;
              dy = gcamn->origy - gcamn_nbr->origy ;
              dz = gcamn->origz - gcamn_nbr->origz ;
              d0 = sqrt(dx*dx + dy*dy + dz*dz) ;
              dx = gcamn->x - gcamn_nbr->x ;
              dy = gcamn->y - gcamn_nbr->y ;
              dz = gcamn->z - gcamn_nbr->z ;
              d = sqrt(dx*dx + dy*dy + dz*dz) ;
              if (FZERO(d))
              {
                continue ;
              }
              error = (d-d0)/d ;
              xdelta -= error*dx ;
              ydelta -= error*dy;
              zdelta -= error*dz ;
              num++ ;
            }
          }
        }
        num = 1 ;
        if (num > 0)
        {
          xdelta /= num ;
          ydelta /= num ;
          zdelta /= num ;
        }
        xdelta *= l_distance ;
        ydelta *= l_distance ;
        zdelta *= l_distance ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_dist: node(%d,%d,%d) distance term (%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, xdelta, ydelta, zdelta) ;
        gcamn->dx += xdelta ;
        gcamn->dy += ydelta ;
        gcamn->dz += zdelta ;
      }

  return(NO_ERROR) ;
}

double
gcamDistanceEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double          sse = 0.0, dx, dy, dz, error, node_sse, d0, d ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if ((gcamn->invalid == GCAM_POSITION_INVALID) ||
            (gcamn->status & GCAM_IGNORE_DISTANCES))
        {
          continue;
        }

        num = 0 ;
        node_sse = 0.0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          xn = MAX(0,xn) ;
          xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            yn = MAX(0,yn) ;
            yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!xk && !yk && !zk)
              {
                continue ;
              }
              zn = z+zk ;
              zn = MAX(0,zn) ;
              zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;

              if (gcamn_nbr->invalid/* == GCAM_POSITION_INVALID*/)
              {
                continue;
              }

#if 0
              if (gcamn_nbr->label != gcamn->label)
              {
                continue ;
              }
#endif
              dx = gcamn->origx - gcamn_nbr->origx ;
              dy = gcamn->origy - gcamn_nbr->origy ;
              dz = gcamn->origz - gcamn_nbr->origz ;
              d0 = sqrt(dx*dx + dy*dy + dz*dz) ;
              dx = gcamn->x - gcamn_nbr->x ;
              dy = gcamn->y - gcamn_nbr->y ;
              dz = gcamn->z - gcamn_nbr->z ;
              d = sqrt(dx*dx + dy*dy + dz*dz) ;
              error = d0-d ;
              num++ ;
              node_sse += error*error ;
            }
          }
        }
        num = 1 ;
        if (num > 0)
        {
          sse += node_sse/num ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("E_dist: node(%d,%d,%d) distance sse %2.3f\n",
                 x, y, z, node_sse/num) ;
      }

  return(sse) ;
}

#define AREA_NEIGHBORS 8
const float jac_scale = 10 ;
#if 1
int
gcamJacobianTerm(GCA_MORPH *gcam, const MRI *mri,
                 double l_jacobian, double ratio_thresh)
{

#ifdef GCAM_JACOB_TERM_GPU
  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamJacobianTermGPU( gcam, l_jacobian, jac_scale );
#else
  int            i=0, j=0, k=0, num /*, xi, yi, zi, xk, yk, zk = 0*/;
  int            n_omp_threads=1;
  double         dx=0.0, dy=0.0, dz=0.0, norm=0.0;
  double         orig_area=0.0, ratio=0.0, max_norm ;
  double         mn[_MAX_FS_THREADS]; /* _MAX_FS_THREADS is in utils.h */
  GCA_MORPH_NODE *gcamn=NULL ;


  if (DZERO(l_jacobian))
  {
    return(NO_ERROR) ;
  }

  num = 0;
#ifdef HAVE_OPENMP
  #pragma omp parallel for reduction (+:num) firstprivate (j,k,gcamn,ratio,orig_area,ratio_thresh) shared(gcam) schedule(static,1)
#endif
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        orig_area = gcamn->orig_area ;
        if (FZERO(orig_area))
        {
          continue ;
        }

        ratio = gcamn->area / orig_area ;
        if (ratio < ratio_thresh)
        {
          num++ ;
#if 0
          for (xk = -1 ; xk <= 1 ; xk++)
          {
            xi = i+xk ;
            if (xi < 0)
            {
              xi = 0 ;
            }
            if (xi >= gcam->width)
            {
              xi = gcam->width-1 ;
            }
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = i+yk ;
              if (yi < 0)
              {
                yi = 0 ;
              }
              if (yi >= gcam->height)
              {
                yi = gcam->height-1 ;
              }
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = i+zk ;
                if (zi < 0)
                {
                  zi = 0 ;
                }
                if (zi >= gcam->depth)
                {
                  zi = gcam->depth-1 ;
                }
                gcamn = &gcam->nodes[xi][yi][zi] ;
                /*
                     gcamn->dx = gcamn->dy = gcamn->dz = 0 ;
                */
              }
            }
          }
#endif
        }
      }
    }
  }

  if (DIAG_VERBOSE_ON)
  {
    printf("  %d nodes compressed more than %2.2f\n", num, ratio_thresh) ;
  }

  max_norm = 0.0 ;

#ifdef HAVE_OPENMP
  #pragma omp parallel
  {
    n_omp_threads = omp_get_num_threads();
  }
#else
  n_omp_threads = 1;
#endif

  int tid=0;
  for (i=0; i<_MAX_FS_THREADS; i++)
  {
    mn[i] = 0.0;
  }

#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate (j,k,gcamn,dx,dy,dz,norm,tid,mn) shared(gcam) schedule(static,1)
#endif
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        dx = gcamn->dx ;
        dy = gcamn->dy ;
        dz = gcamn->dz ;
        norm = sqrt(dx*dx + dy*dy + dz*dz) ;
#ifdef HAVE_OPENMP
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif
        if (norm > mn[tid])
        {
          mn[tid] = norm ;
        }
      }
    }
  }

  for (i = 0; i<n_omp_threads; i++)
    if (mn[tid] > max_norm)
    {
      max_norm = mn[tid];
    }

#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate(j,k,gcamn,dx,dy,dz,norm) shared(gcam,mri,l_jacobian,Gx,Gy,Gz,max_norm)  schedule(static,1)
#endif
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        gcamJacobianTermAtNode(gcam, mri, l_jacobian, i, j, k, &dx, &dy, &dz) ;
#if 1
        norm = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (norm > max_norm*jac_scale && max_norm > 0 && norm > 1)
          /* don't let it get too big, otherwise it's the
             only thing that happens */
        {
          dx *= max_norm/norm ;
          dy *= max_norm/norm ;
          dz *= max_norm/norm ;
        }
#endif
        gcam->nodes[i][j][k].dx += dx ;
        gcam->nodes[i][j][k].dy += dy ;
        gcam->nodes[i][j][k].dz += dz ;
      }
    }
  }

#endif

  return(NO_ERROR) ;
}

#else

int
gcamJacobianTerm(GCA_MORPH *gcam, MRI *mri,
                 double l_jacobian, double ratio_thresh)
{
  int            i, j, k, num, xi, yi, zi, xk, yk, zk = 0 ;
  double         dx, dy, dz, norm, max_norm ;
  GCA_MORPH_NODE *gcamn ;


  // don't let jacobian term dominate
  if (DIAG_VERBOSE_ON)
  {
    printf("%d nodes compressed more than %2.2f\n", num, ratio_thresh) ;
  }

  max_norm = 0.0 ;
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        dx = gcamn->dx ;
        dy = gcamn->dy ;
        dz = gcamn->dz ;
        norm = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (norm > max_norm)
        {
          max_norm = norm ;
        }
      }
    }
  }

  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        if (FZERO(gcamn->orig_area))
        {
          continue ;
        }
        gcamJacobianTermAtNode(gcam, mri, l_jacobian, i, j, k, &dx, &dy, &dz) ;
#if 0
        norm = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (norm > max_norm*jac_scale && norm > 1)
          /* don't let it get too big, otherwise it's
             the only thing that happens */
        {
          dx *= max_norm/norm ;
          dy *= max_norm/norm ;
          dz *= max_norm/norm ;
        }
#endif
        gcam->nodes[i][j][k].jx = dx ;
        gcam->nodes[i][j][k].jy = dy ;
        gcam->nodes[i][j][k].jz = dz ;
      }
    }
  }
#if 0
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        num = 0 ;
        max_norm = 0.0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = i+xk ;
          if (xi < 0 || xi >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = i+yk ;
            if (yi < 0 ||  (yi >= gcam->height))
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = i+zk ;
              if (zi < 0 || (zi >= gcam->depth))
              {
                continue ;
              }
              num++ ;
              gcamn = &gcam->nodes[xi][yi][zi] ;
              norm = sqrt(gcamn->dx*gcamn->dx +
                          gcamn->dy*gcamn->dy +
                          gcamn->dz*gcamn->dz) ;
              if (norm > max_norm)
              {
                max_norm = norm ;
              }
            }
          }
        }
        gcamn = &gcam->nodes[i][j][k] ;
        norm = sqrt(gcamn->jx*gcamn->jx +
                    gcamn->jy*gcamn->jy +
                    gcamn->jz*gcamn->jz) ;
        if (norm > 2*max_norm && max_norm > 0)
        {
          gcamn->jx = gcamn->jx*(2*max_norm/norm) ;
          gcamn->jy = gcamn->jy*(2*max_norm/norm) ;
          gcamn->jz = gcamn->jz*(2*max_norm/norm) ;
        }
      }
    }
  }
#endif

  // now copy jacobian derivatives into total gradient
  for (num = i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        if (i == Gx && j == Gy && k == Gz)
        {
          printf("node(%d,%d,%d) --> Jacobian term: (%2.1f, %2.1f, %2.1f)\n",
                 i, j, k, gcamn->jx, gcamn->jy, gcamn->jz) ;
          DiagBreak() ;
        }
        gcamn->dx += gcamn->jx ;
        gcamn->dy += gcamn->jy ;
        gcamn->dz += gcamn->jz ;
      }
    }
  }

  return(NO_ERROR) ;
}
#endif

int
gcamAreaTerm(GCA_MORPH *gcam, double l_area)
{
  int    i, j, k ;
  double dx, dy, dz ;
  GCA_MORPH_NODE *gcamn ;

  if (DZERO(l_area))
  {
    return(NO_ERROR) ;
  }
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }
        if (FZERO(gcamn->orig_area))
        {
          continue ;
        }
        gcamAreaTermAtNode(gcam, l_area, i, j, k, &dx, &dy, &dz) ;
        gcam->nodes[i][j][k].dx += dx ;
        gcam->nodes[i][j][k].dy += dy ;
        gcam->nodes[i][j][k].dz += dz ;
      }
    }
  }
  return(NO_ERROR) ;
}

int
gcamAreaSmoothnessTerm(GCA_MORPH *gcam, MRI *mri, double l_area_smoothness)
{
  int            i, j, k ;
  double         dx, dy, dz, Ix, Iy, Iz, norm, nc ;
  GCA_MORPH_NODE *gcamn ;

  if (DZERO(l_area_smoothness))
  {
    return(NO_ERROR) ;
  }
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[i][j][k] ;
        if (gcamn->invalid == GCAM_AREA_INVALID ||
            gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue ;
        }
        MRIsampleVolumeGradientFrame(mri, gcamn->x, gcamn->y, gcamn->z,
                                     &Ix, &Iy, &Iz, 0) ;
        gcamAreaTermAtNode(gcam, l_area_smoothness, i, j, k, &dx, &dy, &dz) ;

        // take out component in intensity gradient direction
        norm = sqrt(Ix*Ix+Iy*Iy+Iz*Iz) ;
        if (FZERO(norm))
        {
          norm = 1 ;
        }
        Ix /= norm ;
        Iy /= norm ;
        Iz /= norm ;
        nc = Ix*dx + Iy*dy + Iz*dz ;
        dx -= Ix*nc ;
        dy -= Iy*nc ;
        dz -= Iz*nc ;
        if (i == Gx && j == Gy && k == Gz)
          printf("l_AS(%d, %d, %d): (%2.2f, %2.2f, %2.2f)\n",
                 i, j, k, dx, dy, dz) ;
        gcam->nodes[i][j][k].dx += dx ;
        gcam->nodes[i][j][k].dy += dy ;
        gcam->nodes[i][j][k].dz += dz ;
      }
    }
  }
  return(NO_ERROR) ;
}

int
gcamJacobianTermAtNode( GCA_MORPH *gcam, const MRI *mri,
                        double l_jacobian,
                        int i, int j, int k,
                        double *pdx, double *pdy, double *pdz )
{
  GCA_MORPH_NODE *gcamn, *gcamni, *gcamnj, *gcamnk ;
  float          delta, ratio ;
  int            n, width, height, depth, invert ;
  /* these are all thread local variables */
  /* _MAX_FS_THREADS is defined in utils.h */
  static VECTOR  *v_i[_MAX_FS_THREADS] = {NULL}, *v_j[_MAX_FS_THREADS], *v_k[_MAX_FS_THREADS], *v_j_x_k[_MAX_FS_THREADS], *v_i_x_j[_MAX_FS_THREADS],*v_k_x_i[_MAX_FS_THREADS],*v_grad[_MAX_FS_THREADS],
                                         *v_tmp[_MAX_FS_THREADS] ;
  double         exponent, orig_area, area ;

  *pdx = *pdy = *pdz = 0 ;
  if (DZERO(l_jacobian))
  {
    return(NO_ERROR) ;
  }

#ifdef HAVE_OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;

  if (!v_i[tid])   /* initialize */
  {
    /*printf("tid  doing VectorAlloc = %d\n", tid);*/
    /*
    Note that, while the positions are all stored as
    double in the GCA_MORPH_NODE, matrices always
    hold float (or complex float) data.
    */
    v_i[tid] = VectorAlloc(3, MATRIX_REAL) ;
    v_j[tid] = VectorAlloc(3, MATRIX_REAL) ;
    v_k[tid] = VectorAlloc(3, MATRIX_REAL) ;
    v_grad[tid] = VectorAlloc(3, MATRIX_REAL) ;
    v_j_x_k[tid] = VectorAlloc(3, MATRIX_REAL) ;
    v_i_x_j[tid] = VectorAlloc(3, MATRIX_REAL) ;
    v_k_x_i[tid] = VectorAlloc(3, MATRIX_REAL) ;
    v_tmp[tid] = VectorAlloc(3, MATRIX_REAL) ;
  }
  else
  {
    V3_CLEAR(v_grad[tid]) ;
  }

  for (n = 0 ; n < AREA_NEIGHBORS ; n++)
  {
    /* assign gcamn pointers to appropriate nodes */
    invert = 1 ;
    switch (n)
    {
    default:
    case 0:    /* first do central node */
      if ((i+1 >= width) || (j+1 >= height) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i+1][j][k] ;
      gcamnj = &gcam->nodes[i][j+1][k] ;
      gcamnk = &gcam->nodes[i][j][k+1] ;
      break ;
    case 1:       /*  i-1 */
      if ((i == 0) || (j+1 >= height) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i-1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i-1][j+1][k] ;
      gcamnk = &gcam->nodes[i-1][j][k+1] ;
      break ;
    case 2:       /* j-1 */
      if ((i+1 >= width) || (j == 0) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j-1][k] ;
      gcamni = &gcam->nodes[i+1][j-1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j-1][k+1] ;
      break ;
    case 3:      /* k-1 */
      if ((i+1 >= width) || (j+1 >= height) || (k == 0))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j][k-1] ;
      gcamni = &gcam->nodes[i+1][j][k-1] ;
      gcamnj = &gcam->nodes[i][j+1][k-1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    case 4:  // left-handed central node
      if ((i == 0) || (j == 0) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i-1][j][k] ;
      gcamnj = &gcam->nodes[i][j-1][k] ;
      gcamnk = &gcam->nodes[i][j][k-1] ;
      break ;
    case 5:       /*  i+1 */
      if ((i+1 >= width) || (j == 0) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i+1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i+1][j-1][k] ;
      gcamnk = &gcam->nodes[i+1][j][k-1] ;
      break ;
    case 6:       /* j+1 */
      if ((i == 0) || (j+1 >= height) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j+1][k] ;
      gcamni = &gcam->nodes[i-1][j+1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j+1][k-1] ;
      break ;
    case 7:      /* k+1 */
      if ((i == 0) || (j == 0) || (k+1 >= depth))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k+1] ;
      gcamni = &gcam->nodes[i-1][j][k+1] ;
      gcamnj = &gcam->nodes[i][j-1][k+1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    }
    if (invert > 0)  // right-handed coordinate system
    {
      orig_area = gcamn->orig_area1 ;
      area = gcamn->area1 ;
    }
    else  // left-handed coordinate system
    {
      orig_area = gcamn->orig_area2 ;
      area = gcamn->area2 ;
    }
    if (FZERO(orig_area))
    {
      continue ;
    }

    if (gcamn->invalid  == GCAM_POSITION_INVALID ||
        gcamni->invalid == GCAM_POSITION_INVALID ||
        gcamnj->invalid == GCAM_POSITION_INVALID ||
        gcamnk->invalid == GCAM_POSITION_INVALID)
    {
      continue;
    }
    if (gcamn->invalid ||
        gcamni->invalid ||
        gcamnj->invalid ||
        gcamnk->invalid)
    {
      DiagBreak() ;
    }

    /* num++ ; why is this here?? chc */

    /* compute cross products and area delta */
    /*
      Continuing the note about the matrices above,
      I believe that the subtraction (of double data)
      will be performed in double, but the result
      then stored in single precision
    */
    GCAMN_SUB(gcamni, gcamn, v_i[tid]) ;
    GCAMN_SUB(gcamnj, gcamn, v_j[tid]) ;
    GCAMN_SUB(gcamnk, gcamn, v_k[tid]) ;

    ratio = area / orig_area ;
    if (ratio < 0.1 && ratio > 0)
    {
      DiagBreak() ;
    }
    if (area < 0)
    {
      DiagBreak() ;
    }

    exponent = gcam->exp_k*ratio ;
    if (exponent > MAX_EXP)
    {
      exponent = MAX_EXP ;
    }

    /* don't use -k, since we are moving in the negative gradient direction */
    delta = (invert * gcam->exp_k / orig_area) * (1.0 / (1.0+exp(exponent))) ;

    if (fabs(delta) > 10000)
    {
      DiagBreak() ;
    }

    /* compute cross-products and add the appropriate
      (i.e. scaled by area difference) cross-products to the gradient */
    /*
      Note that by this point, all the data have been downgraded to
      single precision.
      But there are still a lot of subtractions going on.
    */
    switch (n)
    {
    default:
    case 4:
    case 0:    /* first do central node */
      V3_CROSS_PRODUCT(v_j[tid], v_k[tid], v_j_x_k[tid]) ;
      V3_CROSS_PRODUCT(v_k[tid], v_i[tid], v_k_x_i[tid]) ;
      V3_CROSS_PRODUCT(v_i[tid], v_j[tid], v_i_x_j[tid]) ;
      V3_ADD(v_i_x_j[tid], v_j_x_k[tid], v_tmp[tid]) ;
      V3_ADD(v_k_x_i[tid], v_tmp[tid], v_tmp[tid]) ;
      V3_SCALAR_MUL(v_tmp[tid], -delta, v_tmp[tid]) ;
      break ;
    case 5:       /*  i+1 */
    case 1:       /*  i-1 */
      V3_CROSS_PRODUCT(v_j[tid], v_k[tid], v_j_x_k[tid]) ;
      V3_SCALAR_MUL(v_j_x_k[tid], delta, v_tmp[tid]) ;
      break ;
    case 6:      /* j+1 */
    case 2:      /* j-1 */
      V3_CROSS_PRODUCT(v_k[tid], v_i[tid], v_k_x_i[tid]) ;
      V3_SCALAR_MUL(v_k_x_i[tid], delta, v_tmp[tid]) ;
      break ;
    case 7:      /* k+1 */
    case 3:      /* k-1 */
      V3_CROSS_PRODUCT(v_i[tid], v_j[tid], v_i_x_j[tid]) ;
      V3_SCALAR_MUL(v_i_x_j[tid], delta, v_tmp[tid]) ;
      break ;
    }
    V3_ADD(v_tmp[tid], v_grad[tid], v_grad[tid]) ;
  }

  *pdx = l_jacobian*V3_X(v_grad[tid]) ;
  *pdy = l_jacobian*V3_Y(v_grad[tid]) ;
  *pdz = l_jacobian*V3_Z(v_grad[tid]) ;

  if (fabs(*pdx) > 0.02 || fabs(*pdy) > 0.02 || fabs(*pdz) > 0.02)
  {
    DiagBreak() ;
  }
  if (i == Gx && j == Gy && k == Gz)
  {
    gcamn = &gcam->nodes[i][j][k] ;
    printf("l_jaco: node(%d,%d,%d): area=%2.4f, orig_area=%2.4f, "
           "grad=(%2.3f,%2.3f,%2.3f)\n",
           i, j, k, gcamn->area,gcamn->orig_area, *pdx, *pdy, *pdz) ;
    if (fabs(*pdx) > 0.02 || fabs(*pdy) > 0.02 || fabs(*pdz) > 0.02)
    {
      DiagBreak() ;
    }
  }
  return(NO_ERROR) ;
}
#define NLT_PAD 10

int
gcamVolumeChangeTermAtNode(GCA_MORPH *gcam, MRI *mri, double l_area,
                           int i, int j, int k,
                           double *pdx, double *pdy, double *pdz)
{
  GCA_MORPH_NODE *gcamn, *gcamni, *gcamnj, *gcamnk ;
  float          delta, total_delta = 0 ;
  int            n, width = 0, height = 0, depth = 0, num, invert ;
  static VECTOR  *v_i = NULL, *v_j, *v_k,
                  *v_j_x_k, *v_i_x_j,*v_k_x_i,*v_grad, *v_tmp ;
  double          area, orig_area, del_v_scale, uk, image_val, error ;

  del_v_scale=0;

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  if (!v_i)   /* initialize */
  {
    v_i = VectorAlloc(3, MATRIX_REAL) ;
    v_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k = VectorAlloc(3, MATRIX_REAL) ;
    v_grad = VectorAlloc(3, MATRIX_REAL) ;
    v_j_x_k = VectorAlloc(3, MATRIX_REAL) ;
    v_i_x_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k_x_i = VectorAlloc(3, MATRIX_REAL) ;
    v_tmp = VectorAlloc(3, MATRIX_REAL) ;
  }
  else
  {
    V3_CLEAR(v_grad) ;
  }

  if (i == Gx && j == Gy && k == Gz)
  {
    DiagBreak() ;
  }
  for (num = n = 0 ; n < AREA_NEIGHBORS ; n++)
  {
    /* assign gcamn pointers to appropriate nodes */
    invert = 1 ;
    switch (n)
    {
    default:
    case 0:    /* first do central node */
      if ((i+1 >= width) || (j+1 >= height) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i+1][j][k] ;
      gcamnj = &gcam->nodes[i][j+1][k] ;
      gcamnk = &gcam->nodes[i][j][k+1] ;
      break ;
    case 1:       /*  i-1 */
      if ((i == 0) || (j+1 >= height) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i-1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i-1][j+1][k] ;
      gcamnk = &gcam->nodes[i-1][j][k+1] ;
      break ;
    case 2:       /* j-1 */
      if ((i+1 >= width) || (j == 0) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j-1][k] ;
      gcamni = &gcam->nodes[i+1][j-1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j-1][k+1] ;
      break ;
    case 3:      /* k-1 */
      if ((i+1 >= width) || (j+1 >= height) || (k == 0))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j][k-1] ;
      gcamni = &gcam->nodes[i+1][j][k-1] ;
      gcamnj = &gcam->nodes[i][j+1][k-1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    case 4:
      if ((i == 0) || (j == 0) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i-1][j][k] ;
      gcamnj = &gcam->nodes[i][j-1][k] ;
      gcamnk = &gcam->nodes[i][j][k-1] ;
      break ;
    case 5:       /*  i+1 */
      if ((i+1 >= width) || (j == 0) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i+1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i+1][j-1][k] ;
      gcamnk = &gcam->nodes[i+1][j][k-1] ;
      break ;
    case 6:       /* j+1 */
      if ((i == 0) || (j+1 >= height) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j+1][k] ;
      gcamni = &gcam->nodes[i-1][j+1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j+1][k-1] ;
      break ;
    case 7:      /* k+1 */
      if ((i == 0) || (j == 0) || (k+1 >= depth))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k+1] ;
      gcamni = &gcam->nodes[i-1][j][k+1] ;
      gcamnj = &gcam->nodes[i][j-1][k+1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    }

    //
    if (gcamn->invalid  == GCAM_POSITION_INVALID ||
        gcamni->invalid == GCAM_POSITION_INVALID ||
        gcamnj->invalid == GCAM_POSITION_INVALID ||
        gcamnk->invalid == GCAM_POSITION_INVALID ||
        gcamn->gc == NULL ||
        DZERO(gcamn->sum_ci_vi_ui))
    {
      continue;
    }

    num++ ;

    uk = gcamn->gc->means[0] ;
    /*
      this works because the intensity gradients are scaled to remove their
      norm, so this piece of the whole gradient needs to be scaled down as
      well, otherwise it dominates the gradient.
    */
#if 1
    del_v_scale =
      (uk*(gcamn->sum_ci_vi - gcamn->area) -
       gcamn->sum_ci_vi_ui) / (gcamn->sum_ci_vi_ui*gcamn->sum_ci_vi_ui) ;
#else
    del_v_scale =
      (uk*(gcamn->sum_ci_vi - gcamn->area) -
       gcamn->sum_ci_vi_ui) / (gcamn->sum_ci_vi*gcamn->sum_ci_vi) ;
    del_v_scale /= (gcamn->sum_ci_vi_ui) ;
#endif
    MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z,
                             0, SAMPLE_TRILINEAR,&image_val) ;
    error = image_val - gcamn->predicted_val ;
    del_v_scale *= error ;

    /* compute cross products and area delta */
    GCAMN_SUB(gcamni, gcamn, v_i) ;
    GCAMN_SUB(gcamnj, gcamn, v_j) ;
    GCAMN_SUB(gcamnk, gcamn, v_k) ;

    if (invert > 0)  // right-handed coordinate system
    {
      orig_area = gcamn->orig_area1 ;
      area = gcamn->area1 ;
    }
    else  // left-handed coordinate system
    {
      orig_area = gcamn->orig_area2 ;
      area = gcamn->area2 ;
    }

    delta = invert * del_v_scale ;
    total_delta += delta ;  // for diagnostic purposes

    if (fabs(delta) > 10000)
    {
      DiagBreak() ;
    }

    /* compute cross-products and add the appropriate
      (i.e. scaled by area difference) cross-products to the gradient */
    switch (n)
    {
    default:
    case 4:
    case 0:    /* first do central node */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_ADD(v_i_x_j, v_j_x_k, v_tmp) ;
      V3_ADD(v_k_x_i, v_tmp, v_tmp) ;
      V3_SCALAR_MUL(v_tmp, -delta, v_tmp) ;
      break ;
    case 5:       /*  i+1 */
    case 1:       /*  i-1 */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_SCALAR_MUL(v_j_x_k, delta, v_tmp) ;
      break ;
    case 6:      /* j+1 */
    case 2:      /* j-1 */
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_SCALAR_MUL(v_k_x_i, delta, v_tmp) ;
      break ;
    case 7:      /* k+1 */
    case 3:      /* k-1 */
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_SCALAR_MUL(v_i_x_j, delta, v_tmp) ;
      break ;
    }
    V3_ADD(v_tmp, v_grad, v_grad) ;
  }

  *pdx = l_area*V3_X(v_grad) ;
  *pdy = l_area*V3_Y(v_grad) ;
  *pdz = l_area*V3_Z(v_grad) ;

  if (i == Gx && j == Gy && k == Gz)
  {
    gcamn = &gcam->nodes[i][j][k] ;
    printf("l_volume: node(%d,%d,%d): uk=%2.4f, "
           "total_delta = %2.3f, grad=(%2.3f,%2.3f,%2.3f)\n",
           i, j, k, gcamn->gc->means[0], del_v_scale, *pdx, *pdy, *pdz) ;
    if (fabs(*pdx) > 0.02 || fabs(*pdy) > 0.02 || fabs(*pdz) > 0.02)
    {
      DiagBreak() ;
    }
  }
  return(NO_ERROR) ;
}

int
gcamAreaTermAtNode(GCA_MORPH *gcam, double l_area,
                   int i, int j, int k, double *pdx, double *pdy, double *pdz)
{
  GCA_MORPH_NODE *gcamn, *gcamni, *gcamnj, *gcamnk ;
  float          delta ;
  int            n, width = 0, height = 0, depth = 0, num, invert ;
  static VECTOR  *v_i = NULL, *v_j, *v_k, *v_j_x_k, *v_i_x_j,*v_k_x_i,*v_grad,
                  *v_tmp ;
  double         orig_area ;

  *pdx = *pdy = *pdz = 0.0 ;
  if (DZERO(l_area))
  {
    return(NO_ERROR) ;
  }

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  if (!v_i)   /* initialize */
  {
    v_i = VectorAlloc(3, MATRIX_REAL) ;
    v_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k = VectorAlloc(3, MATRIX_REAL) ;
    v_grad = VectorAlloc(3, MATRIX_REAL) ;
    v_j_x_k = VectorAlloc(3, MATRIX_REAL) ;
    v_i_x_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k_x_i = VectorAlloc(3, MATRIX_REAL) ;
    v_tmp = VectorAlloc(3, MATRIX_REAL) ;
  }
  else
  {
    V3_CLEAR(v_grad) ;
  }

  for (num = n = 0 ; n < AREA_NEIGHBORS ; n++)
  {
    /* assign gcamn pointers to appropriate nodes */
    invert = 1 ;
    switch (n)
    {
    default:
    case 0:    /* first do central node */
      if ((i+1 >= width) || (j+1 >= height) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i+1][j][k] ;
      gcamnj = &gcam->nodes[i][j+1][k] ;
      gcamnk = &gcam->nodes[i][j][k+1] ;
      break ;
    case 1:       /*  i-1 */
      if ((i == 0) || (j+1 >= height) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i-1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i-1][j+1][k] ;
      gcamnk = &gcam->nodes[i-1][j][k+1] ;
      break ;
    case 2:       /* j-1 */
      if ((i+1 >= width) || (j == 0) || (k+1 >= depth))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j-1][k] ;
      gcamni = &gcam->nodes[i+1][j-1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j-1][k+1] ;
      break ;
    case 3:      /* k-1 */
      if ((i+1 >= width) || (j+1 >= height) || (k == 0))
      {
        continue ;
      }
      gcamn = &gcam->nodes[i][j][k-1] ;
      gcamni = &gcam->nodes[i+1][j][k-1] ;
      gcamnj = &gcam->nodes[i][j+1][k-1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    case 4:
      if ((i == 0) || (j == 0) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i-1][j][k] ;
      gcamnj = &gcam->nodes[i][j-1][k] ;
      gcamnk = &gcam->nodes[i][j][k-1] ;
      break ;
    case 5:       /*  i+1 */
      if ((i+1 >= width) || (j == 0) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i+1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i+1][j-1][k] ;
      gcamnk = &gcam->nodes[i+1][j][k-1] ;
      break ;
    case 6:       /* j+1 */
      if ((i == 0) || (j+1 >= height) || (k == 0))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j+1][k] ;
      gcamni = &gcam->nodes[i-1][j+1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j+1][k-1] ;
      break ;
    case 7:      /* k+1 */
      if ((i == 0) || (j == 0) || (k+1 >= depth))
      {
        continue ;
      }
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k+1] ;
      gcamni = &gcam->nodes[i-1][j][k+1] ;
      gcamnj = &gcam->nodes[i][j-1][k+1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    }
    orig_area = gcamn->orig_area ;
    if (FZERO(orig_area))
    {
      continue ;
    }

    //
    if (gcamn->invalid  == GCAM_POSITION_INVALID ||
        gcamni->invalid == GCAM_POSITION_INVALID ||
        gcamnj->invalid == GCAM_POSITION_INVALID ||
        gcamnk->invalid == GCAM_POSITION_INVALID)
    {
      continue;
    }

    num++ ;

    /* compute cross products and area delta */
    GCAMN_SUB(gcamni, gcamn, v_i) ;
    GCAMN_SUB(gcamnj, gcamn, v_j) ;
    GCAMN_SUB(gcamnk, gcamn, v_k) ;

    delta = invert * (orig_area - gcamn->area) ;

    if (fabs(delta) > 10000)
    {
      DiagBreak() ;
    }

    /* compute cross-products and add the appropriate
      (i.e. scaled by area difference) cross-products to the gradient */
    switch (n)
    {
    default:
    case 4:
    case 0:    /* first do central node */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_ADD(v_i_x_j, v_j_x_k, v_tmp) ;
      V3_ADD(v_k_x_i, v_tmp, v_tmp) ;
      V3_NORM(v_tmp) ;
      V3_SCALAR_MUL(v_tmp, -delta, v_tmp) ;
      break ;
    case 5:       /*  i+1 */
    case 1:       /*  i-1 */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_NORM(v_j_x_k) ;
      V3_SCALAR_MUL(v_j_x_k, delta, v_tmp) ;
      break ;
    case 6:      /* j+1 */
    case 2:      /* j-1 */
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_NORM(v_k_x_i) ;
      V3_SCALAR_MUL(v_k_x_i, delta, v_tmp) ;
      break ;
    case 7:      /* k+1 */
    case 3:      /* k-1 */
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_NORM(v_i_x_j) ;
      V3_SCALAR_MUL(v_i_x_j, delta, v_tmp) ;
      break ;
    }
    V3_ADD(v_tmp, v_grad, v_grad) ;
  }

  *pdx = l_area*V3_X(v_grad) ;
  *pdy = l_area*V3_Y(v_grad) ;
  *pdz = l_area*V3_Z(v_grad) ;

  if (i == Gx && j == Gy && k == Gz)
  {
    gcamn = &gcam->nodes[i][j][k] ;
    printf("l_area: node(%d,%d,%d): area=%2.4f, "
           "orig_area=%2.4f, grad=(%2.3f,%2.3f,%2.3f)\n",
           i, j, k, gcamn->area,gcamn->orig_area, *pdx, *pdy, *pdz) ;
    if (fabs(*pdx) > 0.02 || fabs(*pdy) > 0.02 || fabs(*pdz) > 0.02)
    {
      DiagBreak() ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/

void SetGinvalid( const int val )
{
  Ginvalid = val;
}

#define GCAM_CMP_OUTPUT 0
#if 1
int
gcamComputeMetricProperties(GCA_MORPH *gcam)
{
#if GCAM_CMP_OUTPUT
  static unsigned int nCalls = 0;
  const unsigned int gcamCMPoutputFreq=10;
  if( (nCalls%gcamCMPoutputFreq)==0 )
  {
    char fname[STRLEN];
    snprintf( fname, STRLEN-1, "gcamCMPinput%04u", nCalls/gcamCMPoutputFreq );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );
  }
#endif

#ifdef GCAM_CMP_GPU
  gcamComputeMetricPropertiesGPU( gcam, &Ginvalid );
#else
#if SHOW_EXEC_LOC
  printf( "%s: CPU call\n", __FUNCTION__ );
#endif
  double         area1=0.0, area2=0.0 ;
  int            i=0, j=0, k=0, width, height, depth, num=0, neg=0 ;
  int            nthreads=1,tid=0;
  int      gcam_neg_counter[_MAX_FS_THREADS], Ginvalid_counter[_MAX_FS_THREADS];
  GCA_MORPH_NODE *gcamn=NULL, *gcamni=NULL, *gcamnj=NULL, *gcamnk=NULL ;
  VECTOR         *v_i[_MAX_FS_THREADS], *v_j[_MAX_FS_THREADS], *v_k[_MAX_FS_THREADS] ;

  // Ginvalid has file scope and static storage.....
  Ginvalid = 0 ;
#ifdef HAVE_OPENMP
  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
#else
  nthreads = 1;
#endif

  for (i=0; i<nthreads; i++)
  {
    v_i[i] = VectorAlloc(3, MATRIX_REAL) ;
    v_j[i] = VectorAlloc(3, MATRIX_REAL) ;
    v_k[i]  = VectorAlloc(3, MATRIX_REAL) ;
    gcam_neg_counter[i] = 0;
    Ginvalid_counter[i] = 0;
  }
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  gcam->neg = 0 ;

#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate (tid,j,k,gcamn,neg,num,gcamni,gcamnj,gcamnk,area1,area2) shared(gcam,Gx,Gy,Gz,v_i,v_j,v_k,gcam_neg_counter,Ginvalid_counter) schedule(static,1)
#endif
  for (i = 0 ; i < width ; i++)
  {
#ifdef HAVE_OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for (j = 0 ; j < height ; j++)
    {
      for (k = 0 ; k < depth ; k++)
      {
        // get node at this point
        gcamn = &gcam->nodes[i][j][k] ;
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }

        // Test to see if current location is valid
        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          /* Ginvalid++ ; */
          Ginvalid_counter[tid] ++;
          continue;
        }

        neg = num = 0 ;
        gcamn->area = 0.0 ;


        // Compute Jacobean determinants on the 'right'
        if ((i < width-1) && (j < height-1) && (k < depth-1))
        {
          gcamni = &gcam->nodes[i+1][j][k] ;
          gcamnj = &gcam->nodes[i][j+1][k] ;
          gcamnk = &gcam->nodes[i][j][k+1] ;

          if (gcamni->invalid != GCAM_POSITION_INVALID &&
              gcamnj->invalid != GCAM_POSITION_INVALID &&
              gcamnk->invalid != GCAM_POSITION_INVALID)
          {

            num++ ;
            GCAMN_SUB(gcamni, gcamn, v_i[tid]) ;
            GCAMN_SUB(gcamnj, gcamn, v_j[tid]) ;
            GCAMN_SUB(gcamnk, gcamn, v_k[tid]) ;
            // (v_j (x) v_k) (.) v_i (volume)
            area1 = VectorTripleProduct(v_j[tid], v_k[tid], v_i[tid]) ;
            if (area1 <= 0)
            {
              neg = 1 ;
              DiagBreak() ;
            }

            // Store the 'right' Jacobean determinant
            gcamn->area1 = area1 ;

            // Accumulate onto common determinant
            gcamn->area += area1 ;
          }
        }
        else
        {
          // Going to the 'right' would fall out of the volume
          gcamn->area1 = 0 ;
        }


        // Compute Jacobean determinants on the 'left'
        if ((i > 0) && (j > 0) && (k > 0))  /* left-hand coordinate system */
        {
          gcamni = &gcam->nodes[i-1][j][k] ;
          gcamnj = &gcam->nodes[i][j-1][k] ;
          gcamnk = &gcam->nodes[i][j][k-1] ;

          if (gcamni->invalid != GCAM_POSITION_INVALID &&
              gcamnj->invalid != GCAM_POSITION_INVALID &&
              gcamnk->invalid != GCAM_POSITION_INVALID)
          {
            /* invert v_i so that coordinate system is right-handed */
            num++ ;
            GCAMN_SUB(gcamn, gcamni, v_i[tid]) ; // Note args swapped compared to above
            GCAMN_SUB(gcamnj, gcamn, v_j[tid]) ;
            GCAMN_SUB(gcamnk, gcamn, v_k[tid]) ;
            // add two volume
            area2 = VectorTripleProduct(v_j[tid], v_k[tid], v_i[tid]);

            // Store the 'left' Jacobean determinant
            gcamn->area2 = area2 ;

            if (area2 <= 0)
            {
              neg = 1 ;
              DiagBreak() ;
            }

            // Accumulate onto common determinant
            gcamn->area += area2 ;
          }
        }
        else
        {
          // Going to the 'left' would fall out of the volume
          gcamn->area2 = 0 ;
        }

        // Check if at least one Jacobean determinant was computed
        if (num > 0)
        {
          // Store the average of computed determinants in the common determinant
          gcamn->area = gcamn->area / (float)num ; // average volume
        }
        else
        {
          // If no determinants computed, this node becomes invalid
          if (i == Gx && j == Gy && k == Gz)
          {
            DiagBreak() ;
          }
          gcamn->invalid = GCAM_AREA_INVALID ;
          gcamn->area = 0 ;
        }

        // Keep track of determinants which have become negative
        if ((gcamn->invalid == GCAM_VALID) && neg && (gcamn->orig_area > 0))
        {

          if (i > 0 && j > 0 && k > 0&&
              i < gcam->width-1 && j < gcam->height-1 && k < gcam->depth-1)
          {
            DiagBreak() ;
          }

          if (gcam->neg == 0 && getenv("SHOW_NEG"))
          {
            printf("node (%d, %d, %d), label %s (%d) - NEGATIVE!\n",
                   i, j, k, cma_label_to_name(gcamn->label), gcamn->label) ;
          }
          /* gcam->neg++ ; */
          gcam_neg_counter[tid] ++;
        }

        // Add to count of invalid locations
        if (gcamn->invalid)
        {
          /* Ginvalid++ ; */
          Ginvalid_counter[tid] ++;
        }
      }
    }
  }

  for (i=0; i<nthreads; i++)
  {
    VectorFree(&v_i[i]) ;
    VectorFree(&v_j[i]) ;
    VectorFree(&v_k[i]) ;
  }

  for (i = 0; i<nthreads; i++)
  {
    gcam->neg += gcam_neg_counter[i];
    Ginvalid += Ginvalid_counter[i];
  }

#endif

#if GCAM_CMP_OUTPUT
#if 0
  if( (nCalls%gcamCMPoutputFreq)==0 )
  {
    snprintf( fname, STRLEN-1, "after%04u", nCalls/gcamCMPoutputFreq );
    fname[STRLEN-1] = '\0';
    WriteGCAMforMetricProperties( gcam, fname );
  }
#endif
  nCalls++;
#endif

  return(NO_ERROR) ;
}
#else
int
gcamComputeMetricProperties(GCA_MORPH *gcam)
{
  double         area ;
  int            i, j, k, width, height, depth, num, n, invert ;
  GCA_MORPH_NODE *gcamn, *gcamni, *gcamnj, *gcamnk ;
  VECTOR         *v_i, *v_j, *v_k ;

  Ginvalid = 0 ;
  v_i = VectorAlloc(3, MATRIX_REAL) ;
  v_j = VectorAlloc(3, MATRIX_REAL) ;
  v_k = VectorAlloc(3, MATRIX_REAL) ;
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  gcam->neg = 0 ;
  for (i = 0 ; i < width ; i++)
  {
    for (j = 0 ; j < height ; j++)
    {
      for (k = 0 ; k < depth ; k++)
      {
        // get node at this point
        gcamn = &gcam->nodes[i][j][k] ;
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }
        area = 0.0 ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          Ginvalid++ ;
          continue;
        }

        for (area = 0.0, num = n = 0 ; n < AREA_NEIGHBORS ; n++)
        {
          /* assign gcamn pointers to appropriate nodes */
          invert = 1 ;
          switch (n)
          {
          default:
          case 0:    /* first do central node */
            if ((i+1 >= width) || (j+1 >= height) || (k+1 >= depth))
            {
              continue ;
            }
            gcamn = &gcam->nodes[i][j][k] ;
            gcamni = &gcam->nodes[i+1][j][k] ;
            gcamnj = &gcam->nodes[i][j+1][k] ;
            gcamnk = &gcam->nodes[i][j][k+1] ;
            break ;
          case 1:       /*  i-1 */
            if ((i == 0) || (j+1 >= height) || (k+1 >= depth))
            {
              continue ;
            }
            gcamn = &gcam->nodes[i-1][j][k] ;
            gcamni = &gcam->nodes[i][j][k] ;
            gcamnj = &gcam->nodes[i-1][j+1][k] ;
            gcamnk = &gcam->nodes[i-1][j][k+1] ;
            break ;
          case 2:       /* j-1 */
            if ((i+1 >= width) || (j == 0) || (k+1 >= depth))
            {
              continue ;
            }
            gcamn = &gcam->nodes[i][j-1][k] ;
            gcamni = &gcam->nodes[i+1][j-1][k] ;
            gcamnj = &gcam->nodes[i][j][k] ;
            gcamnk = &gcam->nodes[i][j-1][k+1] ;
            break ;
          case 3:      /* k-1 */
            if ((i+1 >= width) || (j+1 >= height) || (k == 0))
            {
              continue ;
            }
            gcamn = &gcam->nodes[i][j][k-1] ;
            gcamni = &gcam->nodes[i+1][j][k-1] ;
            gcamnj = &gcam->nodes[i][j+1][k-1] ;
            gcamnk = &gcam->nodes[i][j][k] ;
            break ;
          case 4:
            if ((i == 0) || (j == 0) || (k == 0))
            {
              continue ;
            }
            invert = -1 ;
            gcamn = &gcam->nodes[i][j][k] ;
            gcamni = &gcam->nodes[i-1][j][k] ;
            gcamnj = &gcam->nodes[i][j-1][k] ;
            gcamnk = &gcam->nodes[i][j][k-1] ;
            break ;
          case 5:       /*  i+1 */
            if ((i+1 >= width) || (j == 0) || (k == 0))
            {
              continue ;
            }
            invert = -1 ;
            gcamn = &gcam->nodes[i+1][j][k] ;
            gcamni = &gcam->nodes[i][j][k] ;
            gcamnj = &gcam->nodes[i+1][j-1][k] ;
            gcamnk = &gcam->nodes[i+1][j][k-1] ;
            break ;
          case 6:       /* j+1 */
            if ((i == 0) || (j+1 >= height) || (k == 0))
            {
              continue ;
            }
            invert = -1 ;
            gcamn = &gcam->nodes[i][j+1][k] ;
            gcamni = &gcam->nodes[i-1][j+1][k] ;
            gcamnj = &gcam->nodes[i][j][k] ;
            gcamnk = &gcam->nodes[i][j+1][k-1] ;
            break ;
          case 7:      /* k+1 */
            if ((i == 0) || (j == 0) || (k+1 >= depth))
            {
              continue ;
            }
            invert = -1 ;
            gcamn = &gcam->nodes[i][j][k+1] ;
            gcamni = &gcam->nodes[i-1][j][k+1] ;
            gcamnj = &gcam->nodes[i][j-1][k+1] ;
            gcamnk = &gcam->nodes[i][j][k] ;
            break ;
          }

          if (gcamn->invalid  == GCAM_POSITION_INVALID ||
              gcamni->invalid == GCAM_POSITION_INVALID ||
              gcamnj->invalid == GCAM_POSITION_INVALID ||
              gcamnk->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }

          num++ ;
          GCAMN_SUB(gcamni, gcamn, v_i) ;
          GCAMN_SUB(gcamnj, gcamn, v_j) ;
          GCAMN_SUB(gcamnk, gcamn, v_k) ;
          // (v_j (x) v_k) (.) v_i (volume)
          area += invert * VectorTripleProduct(v_j, v_k, v_i) ;
        }
        gcamn = &gcam->nodes[i][j][k] ;
        if (num > 0)
        {
          gcamn->area = area / (float)num ;  // average volume
        }
        else
        {
          if (i == Gx && j == Gy && k == Gz)
          {
            DiagBreak() ;
          }
          gcamn->invalid = GCAM_AREA_INVALID ;
          gcamn->area = 0 ;
        }
        if (gcamn->area < 0.4 || gcamn->area > 0.5)
        {
          DiagBreak() ;
        }
        if ((gcamn->invalid == 0) && (area <= 0) && (gcamn->orig_area > 0))
        {
          if (i > 0 && j > 0 && k > 0&&
              i < gcam->width-1 && j < gcam->height-1 && k < gcam->depth-1)
          {
            DiagBreak() ;
          }
          gcam->neg++ ;
        }
        if (gcamn->invalid)
        {
          Ginvalid++ ;
        }
      }
    }
  }

  VectorFree(&v_i) ;
  VectorFree(&v_j) ;
  VectorFree(&v_k) ;
  return(NO_ERROR) ;
}
#endif


#define GCAM_JACOBENERGY_OUTPUT 0

double
gcamJacobianEnergy( const GCA_MORPH *gcam, MRI *mri)
{

  double sse = 0;

#if GCAM_JACOBENERGY_OUTPUT
  const unsigned int gcamJacobEnergyOutputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%gcamJacobEnergyOutputFreq)==0 )
  {
    char fname[STRLEN];
    snprintf( fname, STRLEN-1,
              "gcamJacobEnergyInput%04u", nCalls/gcamJacobEnergyOutputFreq );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );
    snprintf( fname, STRLEN-1,
              "mriJacobEnergyInput%04u.mgz", nCalls/gcamJacobEnergyOutputFreq );
    MRIwrite( mri, fname ); // 'mri' is not const because of this call
  }
#endif

#ifdef GCAM_JACOBENERGY_GPU
#if SHOW_EXEC_LOC
  printf( "%s: CUDA call\n", __FUNCTION__ );
#endif
  if( gcam->ninputs != 1 )
  {
    printf( "%s: ninputs = %i\n", __FUNCTION__, gcam->ninputs );
  }

  sse = gcamJacobianEnergyGPU( gcam, mri );
#else
#if SHOW_EXEC_LOC
  printf( "%s: CPU call\n", __FUNCTION__ );
#endif
  double          delta=0.0, ratio=0.0, exponent=0.0, thick ;
  int             i=0, j=0, k=0, width, height, depth ;
  GCA_MORPH_NODE *gcamn=NULL ;

  thick = mri ? mri->thick : 1.0 ;
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;

  // Note sse initialised to zero here
  sse = 0.0f;
#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate(j,k,gcamn,ratio,exponent,delta) shared(width,height,depth,gcam) reduction(+:sse) schedule(static,1)
#endif
  for (i = 0 ; i < width ; i++)
  {
    for (j = 0 ; j < height ; j++)
    {
      for (k = 0 ; k < depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid)
        {
          continue;
        }

        /* scale up the area coefficient if the area of the current node is
          close to 0 or already negative */
        if (!FZERO(gcamn->orig_area1))
        {
          ratio = gcamn->area1 / gcamn->orig_area1 ;
          exponent = -gcam->exp_k*ratio ;
          if (exponent > MAX_EXP)
          {
            delta = 0.0 ;
          }
          else
          {
            delta = log(1+exp(exponent)) /*   / gcam->exp_k */ ;
          }

          sse += delta * thick ;

          if (!finitep(delta) || !finitep(sse))
          {
            DiagBreak() ;
          }

          if (i == Gx && j == Gy && k == Gz)
          {
            printf("E_jaco: node(%d,%d,%d): area1=%2.4f, error=%2.3f\n",
                   i, j, k, gcamn->area1,delta);
          }

          if (!FZERO(delta))
          {
            DiagBreak() ;
          }
        }

        if (!FZERO(gcamn->orig_area2))
        {
          ratio = gcamn->area2 / gcamn->orig_area2 ;
          exponent = -gcam->exp_k*ratio ;

          if (exponent > MAX_EXP)
          {
            delta = MAX_EXP ;
          }
          else
          {
            delta = log(1+exp(exponent)) /*   / gcam->exp_k */ ;
          }

          sse += delta * thick ;

          if (!finitep(delta) || !finitep(sse))
          {
            DiagBreak() ;
          }

          if (i == Gx && j == Gy && k == Gz)
          {
            printf("E_jaco: node(%d,%d,%d): area2=%2.4f, error=%2.3f\n",
                   i, j, k, gcamn->area2,delta);
          }

          if (!FZERO(delta))
          {
            DiagBreak() ;
          }
        }
      }
    }
  }
#endif


#if GCAM_JACOBENERGY_OUTPUT
  nCalls++;
#endif

  return(sse) ;
}

double
gcamAreaEnergy(GCA_MORPH *gcam)
{
  double          sse = 0.0, error ;
  int             i, j, k, width, height, depth ;
  GCA_MORPH_NODE *gcamn ;

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (sse = 0.0, i = 0 ; i < width ; i++)
  {
    for (j = 0 ; j < height ; j++)
    {
      for (k = 0 ; k < depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid)
        {
          continue;
        }

        error = gcamn->area - gcamn->orig_area ;

        sse += (error*error) ;
        if (!finitep(error) || !finitep(sse))
        {
          DiagBreak() ;
        }
        if (i == Gx && j == Gy && k == Gz)
          printf("E_area: node(%d,%d,%d): area=%2.4f, error=%2.3f\n",
                 i, j, k, gcamn->area,error);
      }
    }
  }

  return(sse) ;
}

/*
  GCAMmorphFromAtlas:
  Applied inverse gcam morph to input. Currently NN and non-NN interpolation
  cases are treated differently.
*/
MRI *
GCAMmorphFromAtlas(MRI *mri_in, GCA_MORPH *gcam, MRI *mri_morphed, int sample_type)
{
  if (!sample_type)   // NN interpolation
  {
    TRANSFORM   _transform, *transform = &_transform ;
    int         x, y, z, f ;
    double      xr, yr, zr, scale, val ;
    float       xf, yf, zf ;

    if (mri_morphed == NULL)
      mri_morphed = MRIallocSequence(gcam->image.width,
                                     gcam->image.height,
                                     gcam->image.depth,
                                     mri_in->type,
                                     mri_in->nframes);

    useVolGeomToMRI(&gcam->image, mri_morphed);
    transform->type = MORPH_3D_TYPE ;
    transform->xform = (void *)gcam ;
    TransformInvert(transform, mri_morphed) ; //NOTE: if inverse exists, it does nothing
    scale = 1;
    //scale = gcam->spacing / mri_in->xsize ;

    for (x = 0; x < mri_morphed->width; x++)
      for (y = 0; y < mri_morphed->height; y++)
        for (z = 0; z < mri_morphed->depth; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          if (gcam->gca)
          {
            if (GCAsourceVoxelToPriorReal(gcam->gca, mri_morphed, transform,
                                          x, y, z, &xr, &yr, &zr)
                != NO_ERROR)
            {
              continue ;
            }
          }
          else
          {
            if (GCAMsampleInverseMorph(gcam,
                                       (float)x, (float)y, (float)z,
                                       &xf, &yf, &zf) != NO_ERROR)
            {
              continue ;
            }
            xr = (double)xf ;
            yr = (double)yf ;
            zr = (double)zf ;
          }

          xr *= scale ;
          yr *= scale ;
          zr *= scale ;

          for (f = 0; f < mri_morphed->nframes; f++)
          {
            MRIsampleVolumeFrameType(mri_in, xr, yr, zr, f, sample_type, &val) ;
            MRIsetVoxVal(mri_morphed, x, y, z, f, val) ;
          }
        }

    return(mri_morphed) ;
  }
  else
  {
#if 0
    MRI *mri_warp ;

    mri_warp = MRIallocSequence(mri_in->width,
                                mri_in->height,
                                mri_in->depth,
                                MRI_FLOAT, 3) ;
    GCAMwriteInverseWarpToMRI(gcam, mri_warp) ;

    mri_morphed = MRIapplyMorph(mri_in, mri_warp, NULL, sample_type) ;
    useVolGeomToMRI(&gcam->image, mri_morphed);

    MRIfree(&mri_warp) ;
    return(mri_morphed) ;
#else
    int        width, height, depth, frames, x, y, z, f,
               xm1, ym1, zm1, xp1, yp1, zp1 ;
    float      xd, yd, zd, dx, dy, dz, thick ;
    double     weight;
    MRI        *mri_weights, *mri_ctrl, *mri_s_morphed ;

    // GCAM is a non-linear voxel-to-voxel transform
    // it also assumes that the uniform voxel size
    if ( (mri_in->xsize != mri_in->ysize)
         || (mri_in->xsize != mri_in->zsize)
         || (mri_in->ysize != mri_in->zsize))
    {
      ErrorExit
      (ERROR_BADPARM,
       "non-uniform volumes cannot be used for GCAMmorphFromAtlas()\n");
    }

    width = mri_in->width ;
    height = mri_in->height ;
    depth = mri_in->depth ;
    frames = mri_in->nframes ;
    thick = mri_in->thick ;

    float orig_vals[frames];
    float vals[frames];

    // uses the input volume size
    mri_weights = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_in, mri_weights);
    //mri_s_morphed = MRIalloc(width, height, depth, MRI_FLOAT) ;
    mri_s_morphed = MRIallocSequence(width, height, depth, MRI_FLOAT, frames) ;
    MRIcopyHeader(mri_in, mri_s_morphed);

    MRI_BSPLINE * bspline = NULL;
    if (sample_type == SAMPLE_CUBIC_BSPLINE)
    {
      bspline = MRItoBSpline(mri_in,NULL,3);
    }


    // loop over input volume indices
    for (x = 0 ; x < width ; x++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (z = 0 ; z < depth ; z++)
        {

          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          /* compute voxel coordinates of this morph point */
          if (!GCAMsampleMorph(gcam,
                               (float)x*thick, (float)y*thick, (float)z*thick,
                               &xd, &yd, &zd))
          {
            xd /= thick ;
            yd /= thick ;
            zd /= thick ;  /* voxel coords */

            /* now use trilinear interpolation */
            xm1 = (int)floor(xd) ;
            ym1 = (int)floor(yd) ;
            zm1 = (int)floor(zd) ;
            xp1 = xm1 + 1 ;
            yp1 = ym1 + 1 ;
            zp1 = zm1 + 1 ;

            /* make sure they are within bounds */
            xm1 = mri_in->xi[xm1] ;
            ym1 = mri_in->yi[ym1] ;
            zm1 = mri_in->zi[zm1] ;
            xp1 = mri_in->xi[xp1] ;
            yp1 = mri_in->yi[yp1] ;
            zp1 = mri_in->zi[zp1] ;


            dx = xd - xm1 ;
            dy = yd - ym1 ;
            dz = zd - zm1 ;

#if 1
            if (sample_type == SAMPLE_CUBIC_BSPLINE)
              // recommended to externally call this and keep mri_coeff
              // if image is resampled often (e.g. in registration algo)
            {
              MRIsampleSeqBSpline(bspline, x, y, z, orig_vals, 0,frames-1);
            }
            else
            {
              MRIsampleSeqVolumeType(mri_in, x, y, z, 
                                     orig_vals, 0, frames-1, sample_type);
            }
#else
            //orig_val = MRIgetVoxVal(mri_in, x, y, z, 0) ;
            MRIsampleSeqVolume(mri_in, x, y, z, orig_vals, 0, frames-1) ;
#endif

            /* now compute contribution to each of 8 nearest voxels */
            if (xm1 == Gx && ym1 == Gy && zm1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (1-dx) * (1-dy) * (1-dz) ;
            MRIFvox(mri_weights,xm1,ym1,zm1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xm1,ym1,zm1,f) 
                += weight * orig_vals[f] ;
            }


            if (xp1 == Gx && ym1 == Gy && zm1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (dx) * (1-dy) * (1-dz) ;
            MRIFvox(mri_weights,xp1,ym1,zm1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xp1,ym1,zm1,f) 
                += weight * orig_vals[f] ;
            }


            if (xm1 == Gx && yp1 == Gy && zm1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (1-dx) * (dy) * (1-dz) ;
            MRIFvox(mri_weights,xm1,yp1,zm1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xm1,yp1,zm1,f) 
                += weight * orig_vals[f] ;
            }


            if (xm1 == Gx && ym1 == Gy && zp1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (1-dx) * (1-dy) * (dz) ;
            MRIFvox(mri_weights,xm1,ym1,zp1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xm1,ym1,zp1,f) 
                += weight * orig_vals[f] ;
            }


            if (xp1 == Gx && yp1 == Gy && zm1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (dx) * (dy) * (1-dz) ;
            MRIFvox(mri_weights,xp1,yp1,zm1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xp1,yp1,zm1,f) 
                += weight * orig_vals[f] ;
            }

            if (xp1 == Gx && ym1 == Gy && zp1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (dx) * (1-dy) * (dz) ;
            MRIFvox(mri_weights,xp1,ym1,zp1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xp1,ym1,zp1,f) 
                += weight * orig_vals[f] ;
            }

            if (xm1 == Gx && yp1 == Gy && zp1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (1-dx) * (dy) * (dz) ;
            MRIFvox(mri_weights,xm1,yp1,zp1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xm1,yp1,zp1,f) 
                += weight * orig_vals[f] ;
            }

            if (xp1 == Gx && yp1 == Gy && zp1 == Gz)
            {
              DiagBreak() ;
            }
            weight = (dx) * (dy) * (dz) ;
            MRIFvox(mri_weights,xp1,yp1,zp1) += weight ;
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,xp1,yp1,zp1,f) 
                += weight * orig_vals[f] ;
            }
          }
        }
      }
    }
    if (bspline)
    {
      MRIfreeBSpline(&bspline);
    }
    if (sample_type == SAMPLE_CUBIC_BSPLINE)
    {
      bspline = MRItoBSpline(mri_s_morphed,NULL,3);
    }

    /* now normalize weights and values */
    for (x = 0 ; x < width ; x++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (z = 0 ; z < depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          weight = (float)MRIFvox(mri_weights,x,y,z) ;
          if (!FZERO(weight))
          {
            if (sample_type == SAMPLE_CUBIC_BSPLINE)
              // recommended to externally call this and keep mri_coeff
              // if image is resampled often (e.g. in registration algo)
            {
              MRIsampleSeqBSpline(bspline, x, y, z, vals, 0, frames-1);
            }
            else
            {
              MRIsampleSeqVolume(mri_s_morphed, x, y, z, vals, 0, frames-1);
            }
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_s_morphed,x,y,z,f) = 
                (float)((double)vals[f]/weight) ;
            }
          }
        }
      }
    }
    if (bspline)
    {
      MRIfreeBSpline(&bspline);
    }

    /* copy from short image to BUFTYPE one */
    if (!mri_morphed)
    {
      mri_morphed = MRIclone(mri_in, NULL) ;
    }
    else
    {
      MRIclear(mri_morphed) ;
    }

    for (x = 0 ; x < width ; x++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (z = 0 ; z < depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          switch (mri_morphed->type)
          {
          case MRI_UCHAR:
            for (f = 0 ; f < frames ; f++)
            {
              MRIseq_vox(mri_morphed,x,y,z,f) = 
                (uchar)MRIFseq_vox(mri_s_morphed,x,y,z,f) ;
            }
            break ;
          case MRI_SHORT:
            for (f = 0 ; f < frames ; f++)
            {
              MRISseq_vox(mri_morphed,x,y,z,f) = 
                (short)MRIFseq_vox(mri_s_morphed,x,y,z,f) ;
            }
            break ;
          case MRI_FLOAT:
            for (f = 0 ; f < frames ; f++)
            {
              MRIFseq_vox(mri_morphed,x,y,z,f) = 
                MRIFseq_vox(mri_s_morphed,x,y,z,f) ;
            }
            break ;
          default:
            ErrorReturn(NULL,
                        (ERROR_UNSUPPORTED,
                         "GCAMmorphFromAtlas: unsupported volume type %d",
                         mri_morphed->type)) ;
            break ;

          }
        }
      }
    }

    MRIfree(&mri_s_morphed) ;

    /* run soap bubble to fill in remaining holes */
#if 1
    //mri_ctrl = MRIcloneDifferentType(mri_in, MRI_UCHAR) ;
    mri_ctrl = MRIcloneDifferentType(mri_weights, MRI_UCHAR) ;
#else
    mri_ctrl = mri_weights ;
#endif

    for (x = 0 ; x < width ; x++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (z = 0 ; z < depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          weight = (float)MRIFvox(mri_weights,x,y,z) ;
          if (weight > .1)
          {
            MRIvox(mri_ctrl, x, y, z) = 1 ;
          }
          else
          {
            MRIvox(mri_ctrl, x, y, z) = 0 ;
          }
        }
      }
    }

#if 1
    MRIfree(&mri_weights) ;
#endif
    /*
      disable Voronoi stuff - it will extrapolate morph in regions
      that won't be accurate anyway.
      MRIbuildVoronoiDiagram(mri_morphed, mri_ctrl, mri_morphed) ;
    */

    if( sample_type != 0) // non-NN interpolation
    {
      MRIsoapBubble(mri_morphed, mri_ctrl, mri_morphed, 3*gcam->spacing, -1) ;
    }
    MRIfree(&mri_ctrl) ;

    // use gcam src information to the morphed image
    useVolGeomToMRI(&gcam->image, mri_morphed);
    return(mri_morphed) ;
#endif
  }
}

/*
  GCAMmorphPlistFromAtlas:
  Applies inverse gcam morph to input point list.
*/
int
GCAMmorphPlistFromAtlas(int N,
                        float *points_in,
                        GCA_MORPH *gcam,
                        float *points_out)
{
  int sample_type = 0; // more options for sampling types: on the to do list
  int counter = 0;

  if (!sample_type)   // NN interpolation
  {
    TRANSFORM   _transform, *transform = &_transform ;
    double      x, y, z;
    double      xr, yr, zr, scale;

    transform->type = MORPH_3D_TYPE ;
    transform->xform = (void *)gcam ;
    TransformInvert(transform, NULL) ; //NOTE: if inverse exists / mri = NULL, it does nothing
    scale = 1;
    //scale = gcam->spacing / mri_in->xsize ;

    for (counter = 0; counter < N; counter++)
    {
      x = points_in[counter*3];
      y = points_in[counter*3+1];
      z = points_in[counter*3+2];
      if (x == Gx && y == Gy && z == Gz)
      {
        DiagBreak() ;
      }
      //if (GCAsourceVoxelToPriorReal(gcam->gca, NULL, transform,
      //                              x, y, z, &xr, &yr, &zr) != NO_ERROR)
      if (GCAsourceFloatVoxelToPriorReal(gcam->gca, NULL, transform, 
					 x, y, z, &xr, &yr, &zr) != NO_ERROR)
        continue ;
      
      xr *= scale ;
      yr *= scale ;
      zr *= scale ;
      points_out[counter*3] = xr;
      points_out[counter*3+1] = yr;
      points_out[counter*3+2] = zr;
    }
  }
  return(1);
}

/*
  GCAMmorphPlistToSource:
  Transforms a list of points from the target space back to the source space.
*/
int
GCAMmorphPlistToSource(int N,
                       float *points_in, 
                       GCA_MORPH *gcam,
                       float *points_out)
{
  int sample_type = 0; // more options for sampling types: on the to do list
  int counter = 0;

  if (!sample_type)   // NN interpolation
  {
    TRANSFORM   _transform, *transform = &_transform ;
    float      x, y, z;
    float      xr, yr, zr, scale;

    transform->type = MORPH_3D_TYPE ;
    transform->xform = (void *)gcam ;
    TransformInvert(transform, NULL) ; //NOTE: if inverse exists / mri = NULL, it does nothing
    scale = 1;

    for (counter = 0; counter < N; counter++)
    {
      x = points_in[counter*3];
      y = points_in[counter*3+1];
      z = points_in[counter*3+2];
      if (x == Gx && y == Gy && z == Gz)
      {
        DiagBreak() ;
      }
      if (GCApriorToSourceVoxelFloat(gcam->gca, NULL, transform,
                                     x, y, z, &xr, &yr, &zr) != NO_ERROR)
      {
        continue ;
      }
      xr *= scale ;
      yr *= scale ;
      zr *= scale ;
      points_out[counter*3] = xr;
      points_out[counter*3+1] = yr;
      points_out[counter*3+2] = zr;
    }
  }
  return(1);
}

MRI *
GCAMmorphToAtlas(MRI *mri_src,
                 GCA_MORPH *gcam,
                 MRI *mri_morphed,
                 int frame,
                 int sample_type)
{
  int        width, height, depth, x, y, z, start_frame, end_frame ;
  int        out_of_gcam;
  float      xd, yd, zd ;
  double       val, xoff, yoff, zoff ;

  if (frame >= 0 && frame < mri_src->nframes)
  {
    start_frame = end_frame = frame ;
  }
  else
  {
    start_frame = 0 ;
    end_frame = mri_src->nframes-1 ;
  }

#if 0
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
#else
  // should be 256^3
  width  = gcam->width*gcam->spacing ;
  height = gcam->height*gcam->spacing ;
  depth  = gcam->depth*gcam->spacing ;
#endif


  // GCAM is a non-linear voxel-to-voxel transform
  // it also assumes that the uniform voxel size
  if (mri_morphed)
  {
    if ( (mri_src->xsize != mri_src->ysize)
         || (mri_src->xsize != mri_src->zsize)
         || (mri_src->ysize != mri_src->zsize))
    {
      ErrorExit
      (ERROR_BADPARM,
       "non-uniform volumes cannot be used for GCAMmorphToAtlas()\n");
    }
  }
  if (!mri_morphed)
  {
    // alloc with FOV same as gcam
    mri_morphed = MRIallocSequence
                  //(width, height, depth, MRI_FLOAT, frame < 0 ? mri_src->nframes : 1) ;
                  (width, height, depth, mri_src->type, frame < 0 ? mri_src->nframes : 1) ;
    MRIcopyHeader(mri_src, mri_morphed) ;
  }

  if (getenv("MGH_TAL"))
  {
    xoff = -7.42 ;
    yoff = 24.88 ;
    zoff = -18.85 ;
    printf("INFO: adding MGH tal offset (%2.1f, %2.1f, %2.1f) to xform\n",
           xoff, yoff, zoff) ;
  }
  else
  {
    xoff = yoff = zoff = 0 ;
  }

  MRI_BSPLINE * bspline = NULL;
  if (sample_type == SAMPLE_CUBIC_BSPLINE)
  {
    bspline = MRItoBSpline(mri_src,NULL,3);
  }

  // x, y, z are the col, row, and slice (and xyz) in the gcam/target volume
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        // Should not divide by src thick
        //out_of_gcam = GCAMsampleMorph(gcam, (float)x*mri_src->thick,
        //   (float)y*mri_src->thick,
        //   (float)z*mri_src->thick,
        //   &xd, &yd, &zd);

        // Convert target-crs to input-crs
        out_of_gcam = GCAMsampleMorph(gcam,
                                      (float)x, (float)y, (float)z,
                                      &xd, &yd, &zd);

        if (!out_of_gcam)
        {
          // Should not divide by src thick. If anything,
          // divide by target thick,
          // but its always 1 here anyway
          //xd /= mri_src->thick ;
          // yd /= mri_src->thick ; zd /= mri_src->thick ;
          xd += xoff ;
          yd += yoff ;
          zd += zoff ;
          for (frame = start_frame ; frame <= end_frame ; frame++)
          {
            if (nint(xd) == Gx && nint(yd) == Gy && nint(zd) == Gz)
            {
              DiagBreak() ;
            }

            if (xd > -1 && yd > -1 && zd > 0 &&
                xd < mri_src->width &&
                yd < mri_src->height &&
                zd < mri_src->depth)
            {
              if (sample_type == SAMPLE_CUBIC_BSPLINE)
              {
                MRIsampleBSpline(bspline, xd, yd, zd, frame, &val);
              }
              else
                MRIsampleVolumeFrameType(mri_src, xd, yd, zd,
                                         frame, sample_type, &val) ;
              // printf("Within GCAMmorphToAtlas: (%d, %d, %d): (%f, %f, %f): %f \n", x, y, z, xd, yd, zd, val) ;
            }
            else
            {
              val = 0.0 ;
            }
            MRIsetVoxVal(mri_morphed, x, y, z, frame-start_frame, val) ;
          }
        }
      }
    }
  }
  if (bspline)
  {
    MRIfreeBSpline(&bspline);
  }

  // copy the gcam dst information to the morphed volume
  if (getenv("USE_AVERAGE305"))
  {
    fprintf(stderr,
            "INFO: Environmental variable USE_AVERAGE305 set\n");
    fprintf(stderr,
            "INFO: Modifying dst c_(r,a,s), using average_305 values\n");
    mri_morphed->c_r = -0.0950;
    mri_morphed->c_a = -16.5100;
    mri_morphed->c_s = 9.7500;
    mri_morphed->ras_good_flag = 1;
    // now we cache transform and thus we have to do the following whenever
    // we change direction cosines
    MRIreInitCache(mri_morphed);
  }
  else
  {
    useVolGeomToMRI(&gcam->atlas, mri_morphed);
  }

  return(mri_morphed) ;
}
MRI *
GCAMmorphToAtlasWithDensityCorrection(MRI *mri_src, GCA_MORPH *gcam,
                                      MRI *mri_morphed, int frame)
{
  int        width, height, depth, x, y, z, start_frame, end_frame ;
  float      xd, yd, zd ;
  double       val, jacobian ;
  MRI        *mri_jacobian ;

  mri_morphed = GCAMmorphToAtlas(mri_src, gcam, mri_morphed, frame, SAMPLE_TRILINEAR) ;
  if (!mri_morphed)
  {
    return(NULL) ;
  }

  width = mri_morphed->width ;
  height = mri_morphed->height ;
  depth = mri_morphed->depth ;
  if (frame >= 0 && frame < mri_src->nframes)
  {
    start_frame = end_frame = frame ;
  }
  else
  {
    start_frame = 0 ;
    end_frame = mri_src->nframes-1 ;
  }
  mri_jacobian = gcamCreateJacobianImage(gcam) ;
  if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_jacobian, "jac.mgz") ;
  }
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        xd = (float)x/gcam->spacing ;
        yd = (float)y/gcam->spacing ;
        zd = (float)z/gcam->spacing ;
        if (MRIindexNotInVolume(mri_jacobian, xd, yd, zd))
        {
          continue ;
        }
        if (MRIsampleVolume(mri_jacobian, xd, yd, zd, &jacobian) == NO_ERROR)
        {
          for (frame = start_frame ; frame <= end_frame ; frame++)
          {
            val = MRIgetVoxVal(mri_morphed, x, y, z, frame) ;
            val *= jacobian ;
            MRIsetVoxVal(mri_morphed, x, y, z, frame-start_frame, val) ;
          }
        }
      }
    }
  }

  MRIfree(&mri_jacobian) ;

  return(mri_morphed) ;
}
MRI *
GCAMmorphToAtlasType(MRI *mri_src, GCA_MORPH *gcam, MRI *mri_morphed,
                     int frame, int interp_type)
{
  int        width, height, depth, x, y, z, start_frame, end_frame ;
  float      xd, yd, zd ;
  double       val ;

  if (frame >= 0 && frame < mri_src->nframes)
  {
    start_frame = end_frame = frame ;
  }
  else
  {
    start_frame = 0 ;
    end_frame = mri_src->nframes-1 ;
  }

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  // GCAM is a non-linear voxel-to-voxel transform
  // it also assumes that the uniform voxel size
  if (mri_morphed)
  {
    if ( (mri_src->xsize != mri_src->ysize)
         || (mri_src->xsize != mri_src->zsize)
         || (mri_src->ysize != mri_src->zsize))
    {
      ErrorExit
      (ERROR_BADPARM,
       "non-uniform volumes cannot be used for GCAMmorphToAtlas()\n");
    }
  }
  if (!mri_morphed)
  {
    mri_morphed = MRIallocSequence
                  (width, height, depth, mri_src->type, frame < 0 ? mri_src->nframes : 1) ;
    MRIcopyHeader(mri_src, mri_morphed) ;
  }

  MRI_BSPLINE * bspline = NULL;
  if (interp_type == SAMPLE_CUBIC_BSPLINE)
  {
    bspline = MRItoBSpline(mri_src,NULL,3);
  }

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        if (!GCAMsampleMorph(gcam, (float)x*mri_src->thick,
                             (float)y*mri_src->thick, (float)z*mri_src->thick,
                             &xd, &yd, &zd))
        {
          xd /= mri_src->thick ;
          yd /= mri_src->thick ;
          zd /= mri_src->thick ;
          for (frame = start_frame ; frame <= end_frame ; frame++)
          {
            if (xd > -1 && yd > -1 && zd > 0 &&
                xd < width && yd < height && zd < depth)
            {
              if (interp_type == SAMPLE_CUBIC_BSPLINE)
                // recommended to externally call this and keep mri_coeff
                // if image is resampled often (e.g. in registration algo)
              {
                MRIsampleBSpline(bspline, xd, yd, zd, frame, &val);
              }
              else
                MRIsampleVolumeFrameType
                (mri_src, xd, yd, zd, frame, interp_type, &val) ;
            }
            else
            {
              val = 0.0 ;
            }
            MRIsetVoxVal(mri_morphed, x, y, z, frame-start_frame, val) ;
          }
        }
      }
    }
  }

  // copy the gcam dst information to the morphed volume
  useVolGeomToMRI(&gcam->atlas, mri_morphed);
  if (bspline)
  {
    MRIfreeBSpline(&bspline);
  }

  return(mri_morphed) ;
}


int
log_integration_parms(FILE *fp, GCA_MORPH_PARMS *parms)
{
  char  *cp, host_name[STRLEN] ;

  cp = getenv("HOST") ;
  if (cp)
  {
    strcpy(host_name, cp) ;
  }
  else
  {
    strcpy(host_name, "unknown") ;
  }

  if (!DZERO(parms->l_elastic))
  {
    fprintf(fp,"l_elastic=%2.2f, lambda=%2.2f, mu=%2.2f ", parms->l_elastic, parms->lame_lambda, parms->lame_mu) ;
  }
  if (!DZERO(parms->l_binary))
  {
    fprintf(fp,"l_binary=%2.2f ", parms->l_binary) ;
  }
  if (!DZERO(parms->l_dtrans))
  {
    fprintf(fp,"l_dtrans=%2.2f ", parms->l_dtrans) ;
  }
  if (!DZERO(parms->l_area_intensity))
  {
    fprintf(fp,"l_area_intensity=%2.2f ", parms->l_area_intensity) ;
  }
  if (!DZERO(parms->l_map))
  {
    fprintf(fp,"l_map=%2.2f ", parms->l_map) ;
  }
  if (!DZERO(parms->l_area_smoothness))
  {
    fprintf(fp,"l_area_smoothness=%2.2f ", parms->l_area_smoothness) ;
  }
  if (!DZERO(parms->l_expansion))
  {
    fprintf(fp,"l_expansion=%2.2f ", parms->l_expansion) ;
  }
  if (!DZERO(parms->l_jacobian))
  {
    fprintf(fp,"l_jacobian=%2.2f ", parms->l_jacobian) ;
  }
  if (!DZERO(parms->l_label))
  {
    fprintf(fp,"l_label=%2.2f ", parms->l_label) ;
  }
  if (!DZERO(parms->l_distance))
  {
    fprintf(fp,"l_dist=%2.2f ", parms->l_distance) ;
  }
  if (!DZERO(parms->l_log_likelihood))
  {
    fprintf(fp,"l_likelihood=%2.2f ", parms->l_log_likelihood) ;
  }
  if (!DZERO(parms->l_multiscale))
  {
    fprintf(fp,"l_multiscale=%2.2f ", parms->l_multiscale) ;
  }
  if (!DZERO(parms->l_smoothness))
  {
    fprintf(fp,"l_smoothness=%2.2f ", parms->l_smoothness) ;
  }
  if (!DZERO(parms->l_area))
  {
    fprintf(fp,"l_area=%2.2f ", parms->l_area) ;
  }

  fprintf
  (fp,
   "\ntol=%2.2e, dt=%2.2e, exp_k=%2.1f, momentum=%2.2f, "
   "levels=%d, niter=%d, "
   "lbl_dist=%2.2f, avgs=%d, sigma=%2.1f,type=%d, relabel=%d, neg=%s\n",
   parms->tol, parms->dt, parms->exp_k,parms->momentum, parms->levels,
   parms->niterations, parms->label_dist,parms->navgs,
   parms->sigma,parms->integration_type,
   parms->relabel, parms->noneg > 0 ? "no" : "yes") ;
  return(NO_ERROR) ;
}

#if 0
int
gcamCropNegativeAreaNodeGradient(GCA_MORPH *gcam, float crop)
{
  int            i, j, k /*, num, xi, yi, zi, xk, yk, zk*/ ;
  /*  double         dx, dy, dz, norm, orig_area, ratio ;*/
  GCA_MORPH_NODE *gcamn ;

  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid)
        {
          continue;
        }

        if (FZERO(gcamn->orig_area) || gcamn->area > 0.0)
        {
          continue ;
        }
        gcamn->dx *= crop ;
        gcamn->dy *= crop ;
        gcamn->dz *= crop ;
#if 0
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = i+xk ;
          if (xi < 0)
          {
            xi = 0 ;
          }
          if (xi >= gcam->width)
          {
            xi = gcam->width-1 ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = i+yk ;
            if (yi < 0)
            {
              yi = 0 ;
            }
            if (yi >= gcam->height)
            {
              yi = gcam->height-1 ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = i+zk ;
              if (zi < 0)
              {
                zi = 0 ;
              }
              if (zi >= gcam->depth)
              {
                zi = gcam->depth-1 ;
              }
              gcamn = &gcam->nodes[xi][yi][zi] ;
              gcamn->dx = gcamn->dy = gcamn->dz = 0 ;
            }
          }
        }
#endif
      }
    }
  }

  return(NO_ERROR) ;
}
#endif

// globals to track conditions by which loop can terminate early
#ifndef GCAM_REGISTER_LEVEL_GPU
static int nodesCompressed1, nodesCompressed2, nodesCompressed3;
#endif
static int inconsistentLabelNodes;
static float maxGradient, gradientArea;

int
GCAMregisterLevel(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth,
                  GCA_MORPH_PARMS *parms)
{
#ifdef GCAM_REGISTER_LEVEL_GPU
  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamRegisterLevelGPU( gcam, mri, mri_smooth, parms );
#else
  int  n, nsmall, done = 0, which = GCAM_INTEGRATE_OPTIMAL,
                  max_small, increasing, good_step,reduced, good_step_ever ;
  double rms, last_rms, pct_change, orig_dt, min_dt, orig_j,
         tol, last_pct_change ;
  GCA_MORPH_PARMS jacobian_parms ;

  if (parms->integration_type == GCAM_INTEGRATE_FIXED)
  {
    which = GCAM_INTEGRATE_FIXED ;
  }
  max_small = parms->nsmall ;
  gcamClearMomentum(gcam) ;
  jacobian_parms = *parms ;
  jacobian_parms.l_likelihood =
    jacobian_parms.l_label = \
                             jacobian_parms.l_distance = \
                                 jacobian_parms.l_smoothness = 0.0 ;
  orig_dt = parms->dt ;
  nsmall = 0 ;
  pct_change = 0.0 ;
  if (parms->integration_type == GCAM_INTEGRATE_FIXED)
  {
    increasing = 0/*1*/ ;
  }  /* will be set to 0 if new step is
                              smaller than last */
  else
  {
    increasing = 0 ;  /* don't use it for optimal time-step stuff */
  }
  if (parms->log_fp)
  {
    if (mri)
      fprintf
      (parms->log_fp,
       "GCAMregisterLevel: using voxel thickness %2.1f, "
       "navgs=%d, relabel=%d\n",
       mri->xsize, parms->navgs, parms->relabel) ;
    log_integration_parms(parms->log_fp, parms) ;
    fflush(parms->log_fp) ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    log_integration_parms(stdout, parms) ;
  }

  if (parms->start_t == 0 &&
      parms->l_area_smoothness > 0 && (Gdiag & DIAG_WRITE))
  {
    char fname[STRLEN] ;
    HISTOGRAM *h ;
    sprintf(fname, "%s_jac%4.4d.plt", parms->base_name, parms->start_t) ;
    h = gcamJacobianHistogram(gcam, NULL) ;
    printf("writing jacobian histogram to %s\n", fname) ;
    HISTOplot(h, fname) ;
    HISTOfree(&h) ;
  }
  if (parms->write_iterations && (Gdiag & DIAG_WRITE) && parms->start_t == 0)
  {
    write_snapshot(gcam, mri, parms, 0) ;
  }
  if (parms->uncompress)
  {
    gcamRemoveCompressedNodes(gcam, mri, parms, parms->ratio_thresh) ;
  }

  GCAMremoveStatus(gcam, GCAM_LABEL_NODE) ;
  GCAMremoveStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ;
  last_rms = GCAMcomputeRMS(gcam, mri, parms) ;
  if (parms->log_fp)
  {
    fprintf(parms->log_fp, "%04d: dt=%2.3f, rms=%2.3f, neg=%d, invalid=%d",
            0, 0.0f, last_rms, gcam->neg, Ginvalid) ;
    if (parms->l_binary > 0)
      fprintf
      (parms->log_fp, ", aligned = %d (%2.3f%%)\n",
       Galigned,
       100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
    else
    {
      fprintf(parms->log_fp, "\n") ;
    }
    fflush(parms->log_fp) ;
    fflush(parms->log_fp) ;
  }
  //  if (Gdiag & DIAG_SHOW)
  {
    printf("%04d: dt=%2.3f, rms=%2.3f, neg=%d, invalid=%d",
           0, 0.0f, last_rms, gcam->neg, Ginvalid) ;
    if (parms->l_binary > 0)
      printf(", aligned = %d (%2.3f%%)\n",
             Galigned,
             100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
    else
    {
      printf("\n") ;
    }
    fflush(stdout);
  }

  orig_j = parms->l_jacobian ;
  tol = parms->tol ;
  good_step = good_step_ever = 0 ;
  reduced = 0 ;

  // globals to track conditions by which loop can terminate early
  nodesCompressed1 = nodesCompressed2 = nodesCompressed3 = 1000;
  inconsistentLabelNodes = 1000;
  maxGradient = gradientArea = 1000;
  if (Gdiag & DIAG_SHOW)
  {
    gcamShowCompressed(gcam, stdout) ;
  }

  gcamCheck(gcam, mri) ;
  for (n = parms->start_t ; n < parms->start_t+parms->niterations ; n++)
  {
#if 0
    printf("Loop %d of %d...\n", n,parms->start_t+parms->niterations);
    // check for early termination, rather than looping needlessly
    if ( (n>200) \
         && (nodesCompressed1==0) &&
         (nodesCompressed2==0) &&
         (nodesCompressed3==0) \
         && (inconsistentLabelNodes==0) \
         && (maxGradient==0) && (gradientArea==0))
    {
      printf("GCAMRegisterLevel appears complete.");
      break;
    }
#endif

    /*              parms->l_jacobian = 1 ;*/
    GCAMremoveStatus(gcam, GCAM_LABEL_NODE) ;
    GCAMremoveStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ;
    gcamComputeGradient(gcam, mri, mri_smooth, parms) ;
    parms->l_jacobian = orig_j ;
    gcamWriteDiagnostics(gcam) ;
    switch (parms->integration_type)
    {
    case GCAM_INTEGRATE_OPTIMAL:
      parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ; /* will search around this value */
      min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
      parms->dt = min_dt ;
#if 0
      if (gcam->neg > 0 && parms->noneg == True && 0)  // disabled
      {
	printf("negative nodes detected - removing likelihood from those lattice locations and recomputing\n") ;
	GCAMsetNegativeNodeStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ;
	GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
	gcamClearMomentum(gcam) ;
	gcamComputeMetricProperties(gcam) ;
	gcamComputeGradient(gcam, mri, mri_smooth, parms) ;
	parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ; /* will search around this value */
	min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
	parms->dt = min_dt ;
      }
#endif
      break ;
    case GCAM_INTEGRATE_FIXED:
      min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
      break ;
    case GCAM_INTEGRATE_BOTH:
      if (which == GCAM_INTEGRATE_OPTIMAL)
      {
        check_gcam(gcam) ;
        parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ; /* will search around  this value */
        min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
        check_gcam(gcam) ;
        parms->dt = min_dt ;
        max_small = parms->nsmall ;
        tol = parms->tol ;
      }
      else  /* take some momentum steps */
      {
        max_small = 2*parms->nsmall ;
        tol = parms->tol/2 ;
      }
      break ;
    default:
      min_dt = parms->dt ;
      ErrorExit(ERROR_BADPARM,
                "GCAMregisterLevel: unknown integration type %d",
                parms->integration_type) ;
    }

    GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, SAVED2_POSITIONS) ;
    gcamApplyGradient(gcam, parms) ;
    gcamComputeMetricProperties(gcam) ;
    gcamCheck(gcam, mri) ;
    if (gcam->neg > 0)
    {
      if (gcam_write_neg)  // diagnostic
      {
	MRI *mri ;
	char fname[STRLEN] ;

	mri = GCAMwriteMRI(gcam, NULL, GCAM_NEG) ;
	sprintf(fname, "%s_neg_%4.4d.mgz", parms->base_name, n+1) ;
	printf("writing neg nodes to %s...\n", fname) ;
	MRIwrite(mri, fname) ;
	MRIfree(&mri) ;
	mri = GCAMwriteMRI(gcam, NULL, GCAM_MIN_AREA) ;
	sprintf(fname, "%s_minarea_%4.4d.mgz", parms->base_name, n+1) ;
	printf("writing minarea to %s...\n", fname) ;
	MRIwrite(mri, fname) ;
	MRIfree(&mri) ;
	mri = GCAMwriteMRI(gcam, NULL, GCAM_LOG_MIN_AREA) ;
	sprintf(fname, "%s_log_minarea_%4.4d.mgz", parms->base_name, n+1) ;
	printf("writing log minarea to %s...\n", fname) ;
	MRIwrite(mri, fname) ;
	MRIfree(&mri) ;
      }


#if 0      
      printf("%d negative nodes detected - removing likelihood from those lattice locations and recomputing\n",gcam->neg) ;
      GCAMsetNegativeNodStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ;
      GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
      gcamClearMomentum(gcam) ;
      gcamComputeMetricProperties(gcam) ;
      gcamComputeGradient(gcam, mri, mri_smooth, parms) ;
      GCAMremoveStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ;
      parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ; /* will search around this value */
      min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
      parms->dt = min_dt ;
#endif
    }
    if (parms->constrain_jacobian)
    {
      gcamConstrainJacobian(gcam, mri, parms) ;
    }
    if (Gdiag & DIAG_SHOW)
    {
      gcamShowCompressed(gcam, stdout) ;
    }

    //    GCAMcomputeOriginalProperties(gcam) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      GCAMwrite(gcam, "test.m3z") ;
    }
    if (parms->uncompress/* && which == GCAM_INTEGRATE_OPTIMAL*/)
    {
      gcamRemoveCompressedNodes(gcam, mri, parms, parms->ratio_thresh) ;
    }

    if (gcam->neg > 0 && parms->noneg == True)
    {
      int i = 0 ;

      if (Gdiag & DIAG_SHOW)
        printf("%3.3d: %d negative nodes - clearing momentum...\n",
               i+1, gcam->neg) ;
      GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
      gcamClearMomentum(gcam) ;
      gcamComputeMetricProperties(gcam) ;
      gcamApplyGradient(gcam, parms) ;
      gcamComputeMetricProperties(gcam) ;
      while (gcam->neg > 0)
      {
#if 0
        gcamCropNegativeAreaNodeGradient(gcam, 0.9) ;
#else
        parms->dt *= 0.5 ;
        if (Gdiag & DIAG_SHOW)
          printf("%3.3d: %d negative nodes - reducing timestep to %2.6f...\n",
                 i+1, gcam->neg, parms->dt) ;
        reduced = 1 ;  /* don't reset parms->dt until we
                          switch integration types */
#endif
        gcamUndoGradient(gcam) ;
        gcamComputeMetricProperties(gcam) ;
        gcamApplyGradient(gcam, parms) ;
        gcamComputeMetricProperties(gcam) ;
        if (++i > 50)
        {
          if (gcam->neg > 0)  /* couldn't find a  step without folds */
          {
            GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
            gcamComputeMetricProperties(gcam) ;
          }
          break ;
        }
      }
    }
    min_dt = parms->dt ;

    if (parms->noneg >= 0)
    {
      gcamRemoveNegativeNodes(gcam, mri, parms) ;
      if (gcam->neg > 0)   // couldn't unfold everything
      {
	printf("---------- unfolding failed - restoring original position --------------------\n");
	GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
	gcamComputeMetricProperties(gcam) ;
      }
    }
    if (parms->l_area_smoothness > 0 && (Gdiag & DIAG_WRITE))
    {
      char fname[STRLEN] ;
      HISTOGRAM *h ;
      sprintf(fname, "%s_jac%4.4d.plt", parms->base_name, n+1) ;
      h = gcamJacobianHistogram(gcam, NULL) ;
      printf("writing jacobian histogram to %s\n", fname) ;
      HISTOplot(h, fname) ;
      HISTOfree(&h) ;
    }

    if (parms->write_iterations > 0 &&
        !((n+1) % parms->write_iterations) && (Gdiag & DIAG_WRITE))
    {
      write_snapshot(gcam, mri, parms, n+1) ;
    }
    rms = GCAMcomputeRMS(gcam, mri, parms) ;
    last_pct_change = pct_change ;
    if (FZERO(last_rms))
    {
      pct_change = 0.0 ;
    }
    else
    {
      pct_change = 100.0*(last_rms-rms)/last_rms ;
    }
    if ((pct_change < last_pct_change) || FZERO(pct_change))
    {
      if (increasing && (Gdiag & DIAG_SHOW))
      {
        printf("pct change decreased\n") ;
      }
      increasing = 0 ;
    }
    else /* could check for last_pct_change == 0 here */
    {
      increasing = 1 ;
    }
    if (pct_change <= 0)
    {
      if (Gdiag & DIAG_SHOW)
      {
        printf("rms increased - undoing step...\n") ;
      }
      increasing = 0 ;
      if (parms->constrain_jacobian == 0)
      {
        if (parms->uncompress == False)  // otherwise it could be the
          // uncompressing that caused the sse to dec.
        {
          done = 1 ;
        }
        GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
        gcamComputeMetricProperties(gcam) ;
        rms = GCAMcomputeRMS(gcam, mri, parms) ;
      }
    }

    if (parms->log_fp)
    {
      fprintf(parms->log_fp,
              "%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d",
              n+1, min_dt, rms, pct_change, gcam->neg, Ginvalid) ;
      if (parms->l_binary > 0)
        fprintf
        (parms->log_fp, ", aligned = %d (%2.3f%%)\n",
         Galigned,
         100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
      else
      {
        fprintf(parms->log_fp, "\n") ;
      }
      fflush(parms->log_fp) ;
    }

    //  if (Gdiag & DIAG_SHOW)
    {
      printf("%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d",
             n+1, min_dt, rms, pct_change, gcam->neg, Ginvalid) ;
      if (parms->l_binary > 0)
        printf
        (", aligned = %d (%2.3f%%)\n",
         Galigned,
         100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
      else
      {
        printf("\n") ;
      }
      fflush(stdout);
    }

    if ((pct_change < tol) && !increasing)
    {
      int compressed ;

      if ((max_small > 1) && (Gdiag & DIAG_SHOW))
        printf("\tpct change < tol %2.3f, nsmall = %d of %d\n",
               tol, nsmall+1, max_small) ;
      // if we ever took a good step since the last regridding,
      // and we tried to uncomress and failed, regrid
      gcamCheck(gcam, mri) ;
      if (good_step_ever && parms->regrid > 1 &&
          ((compressed =
              GCAMcountCompressedNodes(gcam, parms->ratio_thresh)) > 0) && 0)
      {
        // lattice is messed up - regrid if regridding
        printf("could not remove %d compressed nodes - regridding....\n",
               compressed) ;
        GCAMregrid(gcam, mri, 0, parms, NULL) ;
        nsmall = 0 ;
        which = GCAM_INTEGRATE_OPTIMAL ;
        rms = GCAMcomputeRMS(gcam, mri, parms) ;
        good_step = good_step_ever = done = 0 ;
      }
      else if ((++nsmall >= max_small) || (pct_change <= 0))
      {
        if (parms->integration_type == GCAM_INTEGRATE_BOTH)
        {
          if (!good_step)
          {
            done++ ;
          }
          if (done >= 2)  /* couldn't take a step with either technique */
          {
            if (good_step_ever && parms->regrid > 1 && 0)
            {
              GCAMregrid(gcam, mri, 0, parms, NULL) ;
              rms = GCAMcomputeRMS(gcam, mri, parms) ;
              nsmall = 0 ;
              which = GCAM_INTEGRATE_OPTIMAL ;
              rms = GCAMcomputeRMS(gcam, mri, parms) ;
              good_step = good_step_ever = done = 0 ;
            }
            else
            {
              n++ ;
              break ;
            }
          }
          else  /* switch integration types */
          {
            reduced = 0 ;
            if (which == GCAM_INTEGRATE_FIXED)
            {
              increasing = 0 ;
              which = GCAM_INTEGRATE_OPTIMAL ;
              parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
              if (parms->uncompress)
              {
                gcamRemoveCompressedNodes(gcam, mri,parms,parms->ratio_thresh);
              }
              /* will search around this value */
            }
            else
            {
#if 0
              if (good_step && parms->regrid == True)
              {
                GCAMregrid(gcam, mri, 0, parms, NULL) ;
                rms = GCAMcomputeRMS(gcam, mri, parms) ;
              }
              else
#endif
              {
                increasing = /*1 */0;
                pct_change = 0.0 ;
                which = GCAM_INTEGRATE_FIXED ;
#if 0
                min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
#else
                if (!DZERO(min_dt))
                {
                  parms->dt = min_dt ;
                }
                else
                {
                  min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
                }
#endif
              }
              if (Gdiag & DIAG_SHOW)
                printf("\tswitching integration type to %s (done=%d)\n",
                       which == GCAM_INTEGRATE_FIXED ? "fixed" : "optimal",
                       done);
            }
            good_step = nsmall = 0 ;
            gcamClearMomentum(gcam) ;
          }
        }
        else
        {
          n++ ;
          break ;
        }
      }
    }
    else if (pct_change >= tol)    /* took at least one good step */
    {
      good_step = good_step_ever = 1 ;
      done = 0 ;    /* for integration type == BOTH, apply both types again */
      nsmall = 0 ;  /* start counting small steps again */
    }
#if 0
    if (parms->relabel)
    {
      GCAMcomputeLabels(mri, gcam) ;
      last_rms = GCAMcomputeRMS(gcam, mri, parms) ;
    }
    else
#endif
      last_rms = rms ;
    /*          parms->label_dist -= 0.5 ;*/
    if (parms->label_dist < 0)
    {
      parms->label_dist = 0 ;
    }
    if (Gprofile > 0 && n >= Gprofile)
    {
      printf("exiting to generate profiling results\n") ;
      exit(0) ;
    }
  }

  parms->start_t = n ;
  parms->dt = orig_dt ;

#endif
  return(NO_ERROR) ;
}


int
mriFillRegion(MRI *mri, int x, int y, int z, int fill_val, int whalf)
{
  int   xi, xk, yi, yk, zi, zk ;

  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri->zi[z+zk] ;
        MRIsetVoxVal(mri, xi, yi, zi, 0, fill_val) ;
      }
    }
  }
  return(NO_ERROR) ;
}


int
write_snapshot(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms, int iter)
{
  char           fname[STRLEN], base_name[STRLEN] ;
  MRI            *mri_morphed, *mri_samples ;
  GCA_MORPH_NODE *gcamn ;
  int            x, y, z, xv, yv, zv ;
  static         int write_samples = -1, write_labels = -1 ;

  /* hack to do things differently for hires->lowres registration */
  if (parms->diag_morph_from_atlas /*!DZERO(parms->l_binary)*/)
    mri_morphed = GCAMmorphFieldFromAtlas(gcam,
                                          parms->mri_diag ? parms->mri_diag :
                                          parms->mri, parms->diag_volume,0,
                                          parms->diag_mode_filter) ;
  else
    mri_morphed =
      GCAMmorphToAtlas(parms->mri_diag ? parms->mri_diag : parms->mri,
                       gcam, NULL, -1, parms->diag_sample_type) ;
  sprintf(base_name, "%s_%4.4d", parms->base_name, iter) ;
  if (getenv("DONT_COMPRESS"))
  {
    sprintf(fname, "%s.mgh", base_name) ;
  }
  else
  {
    sprintf(fname, "%s.mgz", base_name) ;
  }
  printf("writing snapshot to %s\n", fname) ;
  MRIwrite(mri_morphed, fname) ;

  if (getenv("GCAM_WRITE_MEAN") != NULL)
  {
    /* hack to do things differently for hires->lowres registration */
    if (parms->diag_morph_from_atlas)
      mri_morphed =
        GCAMmorphFieldFromAtlas(gcam, parms->mri_diag ? parms->mri_diag :
                                parms->mri, GCAM_MEANS,0, 0) ;
    else
      mri_morphed =
        GCAMmorphToAtlas(parms->mri_diag ? parms->mri_diag : parms->mri,
                         gcam, NULL, -1, SAMPLE_TRILINEAR) ;
    sprintf(base_name, "%s_%4.4d", parms->base_name, iter) ;
    if (getenv("DONT_COMPRESS"))
    {
      sprintf(fname, "%s_intensity.mgh", base_name) ;
    }
    else
    {
      sprintf(fname, "%s_intensity.mgz", base_name) ;
    }
    printf("writing snapshot to %s\n", fname) ;
    MRIwrite(mri_morphed, fname) ;
  }

  if ((parms->diag_morph_from_atlas == 0) || (parms->diag_write_snapshots))
  {
    MRIwriteImageViews(mri_morphed, base_name, IMAGE_SIZE) ;
  }
  MRIfree(&mri_morphed) ;

  if (parms->mri_diag2)
  {
    /* hack to do things differently for hires->lowres registration */
    if (parms->diag_morph_from_atlas /*!DZERO(parms->l_binary)*/)
      mri_morphed = GCAMmorphFieldFromAtlas(gcam,
                                            parms->mri_diag2,
                                            parms->diag_volume,0,
                                            parms->diag_mode_filter) ;
    else
      mri_morphed =
        GCAMmorphToAtlas(parms->mri_diag2,
                         gcam, NULL, -1, parms->diag_sample_type) ;
    sprintf(base_name, "d2.%s_%4.4d", parms->base_name, iter) ;
    if (getenv("DONT_COMPRESS"))
    {
      sprintf(fname, "%s.mgh", base_name) ;
    }
    else
    {
      sprintf(fname, "%s.mgz", base_name) ;
    }
    printf("writing snapshot to %s\n", fname) ;
    MRIwrite(mri_morphed, fname) ;

    if (getenv("GCAM_WRITE_MEAN") != NULL)
    {
      /* hack to do things differently for hires->lowres registration */
      if (parms->diag_morph_from_atlas)
        mri_morphed =
          GCAMmorphFieldFromAtlas(gcam, parms->mri_diag ? parms->mri_diag :
                                  parms->mri, GCAM_MEANS,0, 0) ;
      else
        mri_morphed =
          GCAMmorphToAtlas(parms->mri_diag ? parms->mri_diag : parms->mri,
                           gcam, NULL, -1, SAMPLE_TRILINEAR) ;
      sprintf(base_name, "d2.%s_%4.4d", parms->base_name, iter) ;
      if (getenv("DONT_COMPRESS"))
      {
        sprintf(fname, "%s_intensity.mgh", base_name) ;
      }
      else
      {
        sprintf(fname, "%s_intensity.mgz", base_name) ;
      }
      printf("writing snapshot to %s\n", fname) ;
      MRIwrite(mri_morphed, fname) ;
    }

    if ((parms->diag_morph_from_atlas == 0) || (parms->diag_write_snapshots))
    {
      MRIwriteImageViews(mri_morphed, base_name, IMAGE_SIZE) ;
    }
    MRIfree(&mri_morphed) ;

  }
  mri_samples = MRIclone(parms->mri, NULL) ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        xv = mri_samples->xi[nint(gcamn->x)] ;
        yv = mri_samples->yi[nint(gcamn->y)] ;
        zv = mri_samples->zi[nint(gcamn->z)] ;
        MRIvox(mri_samples, xv, yv, zv) = gcamn->label ;
      }
  if (write_samples < 0)
  {
    write_samples = getenv("GCAM_WRITE_SAMPLES") != NULL ;
  }
  if (write_labels < 0)
  {
    write_labels = getenv("GCAM_WRITE_LABELS") != NULL ;
  }

  if (write_samples > 0)
  {
    sprintf(fname, "%s_fsamples_%4.4d.mgz", parms->base_name, iter) ;
    printf("writing samples to %s....\n", fname) ;
    MRIwrite(mri_samples, fname) ;
  }

  if (write_labels > 0)
  {
    MRIclear(mri_samples) ;
    for (x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;

          if (gcamn->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }
#if 0
          GCApriorToVoxel(gcam->gca, mri_samples, x, y, z, &xv, &yv, &zv) ;
#else
          xv = mri_samples->xi[nint(gcamn->x)] ;
          yv = mri_samples->yi[nint(gcamn->y)] ;
          zv = mri_samples->zi[nint(gcamn->z)] ;
#endif
          if (gcamn->status &
              (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD))
          {
            mriFillRegion(mri_samples, xv, yv, zv, gcamn->label, 1) ;
          }
        }
    sprintf(fname, "%s_labels_%4.4d.mgz", parms->base_name, iter) ;
    printf("writing label samples to %s....\n", fname) ;
    MRIwrite(mri_samples, fname) ;
  }
  MRIfree(&mri_samples) ;

  return(NO_ERROR) ;
}

int boundsCheckf(float x, float y, float z, int width, int height, int depth)
{
  if (x >= width)
  {
    return ERROR_BADPARM;
  }
  else if (y >= height)
  {
    return ERROR_BADPARM;
  }
  else if (z >= depth)
  {
    return ERROR_BADPARM;
  }
  else if (x < 0)
  {
    return ERROR_BADPARM;
  }
  else if (y < 0)
  {
    return ERROR_BADPARM;
  }
  else if (z < 0)
  {
    return ERROR_BADPARM;
  }
  else
  {
    return NO_ERROR;
  }
}

/*----------------------------------------------------------------------------
  The GCAM is an isotropic FOV with voxel size gcam->spacing and dimension
  gcam->{width,height,depth} with 3 frames. The values stored at a voxel
  in this volume are the CRS of this voxel in the Unmorphed space. The morphed
  and unmorphed spaces are assumed to have the same FOV and geometry, but
  possibly different resolution.
  ---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
  GCAMsampleMorph() - given the CRS in the Morphed-anat-space volume
  (but not really), computes the CRS in the UnMorphed-anatomical-space
  volume. This should produce the inverse of GCAMsampleInverseMorph().
  ---------------------------------------------------------------------------*/
int
GCAMsampleMorph( const GCA_MORPH *gcam, float x, float y, float z,
                 float *pxd, float *pyd, float *pzd)
{
  int            xm, xp, ym, yp, zm, zp, width, height, depth ;
  float          xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */
  int            errCode = NO_ERROR;

  /* x, y, z are in MRI voxel coords */
  x /= gcam->spacing ;
  y /= gcam->spacing ;
  z /= gcam->spacing ;
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;

  if ((errCode = boundsCheckf(x, y, z, width, height, depth)) != NO_ERROR)
  {
    return errCode;
  }

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;
  if (
    (gcam->nodes[xm][ym][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xm][ym][zp].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xm][yp][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xm][yp][zp].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][ym][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][ym][zp].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][yp][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][yp][zp].invalid == GCAM_POSITION_INVALID))
  {
    return(ERROR_BADPARM) ;
  }

  *pxd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].x +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].x +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].x +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].x +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].x +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].x +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].x +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].x ;
  *pyd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].y +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].y +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].y +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].y +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].y +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].y +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].y +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].y ;
  *pzd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].z +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].z +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].z +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].z +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].z +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].z +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].z +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].z ;

  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------------
  GCAMsampleInverseMorph() - given the CRS in the anatomical volume, computes
  the CRS in the Morph volume. GCAM must have been inverted. This should
  produce the inverse of GCAMsampleMorph().
  ---------------------------------------------------------------------------*/
int GCAMsampleInverseMorph(GCA_MORPH *gcam,
                           float  cAnat,  float  rAnat,  float  sAnat,
                           float *cMorph, float *rMorph, float *sMorph)
{
  int err;
  double v;

  if (gcam->mri_xind == NULL)
  {
    printf("ERROR: GCAMsampleInverseMorph(): GCAM not inverted\n");
    return(1);
  }

  err = MRIsampleVolume(gcam->mri_xind, cAnat,rAnat,sAnat, &v);
  if (err)
  {
    return(1);
  }
  *cMorph = v * gcam->spacing;

  err = MRIsampleVolume(gcam->mri_yind, cAnat,rAnat,sAnat, &v);
  if (err)
  {
    return(1);
  }
  *rMorph = v * gcam->spacing;

  err = MRIsampleVolume(gcam->mri_zind, cAnat,rAnat,sAnat, &v);
  if (err)
  {
    return(1);
  }
  *sMorph = v * gcam->spacing;

  return(0);
}

/*----------------------------------------------------------------------
  GCAMsampleMorphCheck() - checks the morph by going forwards and backwards.
  Returns 0 if error is less than thresh. Returns 1 if greater. Returns -1
  if there was some other problem. Typical distances are less than 0.5, which
  does not seem all that great to me. Also, not all input CRSs will be valid.
  -------------------------------------------------------------------------*/
int GCAMsampleMorphCheck(GCA_MORPH *gcam, float thresh,
                         float cMorph, float rMorph, float sMorph)
{
  int err;
  float cAnat, rAnat, sAnat;
  float cMorph2, rMorph2, sMorph2;
  float d;

  // Compute the CRS in the Anat space from the CRS in the Morph space
  err = GCAMsampleMorph(gcam, cMorph, rMorph, sMorph, &cAnat, &rAnat, &sAnat);
  if (err)
  {
    return(-1);
  }

  // Feedback the CRS in the Anat space to compute the CRS in the Morph space
  err = GCAMsampleInverseMorph(gcam, cAnat, rAnat, sAnat,
                               &cMorph2, &rMorph2, &sMorph2);
  if (err)
  {
    return(-1);
  }

  // Compare CRS in Morph space
  d = sqrt((cMorph-cMorph2)*(cMorph-cMorph2) +
           (rMorph-rMorph2)*(rMorph-rMorph2) +
           (sMorph-sMorph2)*(sMorph-sMorph2));

  if (Gdiag_no > 0)
  {
    printf("Mcrs=(%5.1f,%5.1f,%5.1f) Mcrs2=(%5.1f,%5.1f,%5.1f), "
           "d=%5.2f th=%g Acrs=(%5.1f,%5.1f,%5.1f)\n",
           cMorph, rMorph, sMorph,
           cMorph2, rMorph2, sMorph2,  d, thresh, cAnat, rAnat, sAnat);

  }

  if (d > thresh)
  {
    return(1);
  }
  return(0);
}
/*----------------------------------------------------------------------------
  GCAMsampleMorphRAS() - given an RAS coord in the Morph space
  (xyzMorph), compute the RAS in the Anat space (xyzAnat).  The Anat
  space RAS is in "tkregister" or "surface" space (ie, cras=0). See
  also GCAMsampleInverseMorphRAS().
  --------------------------------------------------------------------------*/
int GCAMsampleMorphRAS(GCA_MORPH *gcam,
                       float xMorph, float yMorph, float zMorph,
                       float  *xAnat,  float  *yAnat,  float  *zAnat)
{
  static int init_needed = 1, err;
  static MATRIX *ras=NULL, *crs=NULL, *vox2ras=NULL, *ras2vox=NULL;
  float  cMorph, rMorph, sMorph;
  float  cAnat, rAnat, sAnat;

  if (init_needed)
  {
    // Use static so that all this matrix overhead does not have to be
    // done on each call.  Init tkreg vox2ras with the inverse volume
    // as this will be the same as the anatomical volume. The only
    // time this will fail is if this function is called twice with
    // two different volume sizes (eg, 1mm^3 res and a high
    // res). Multiple calls with different subjects is not a problem
    // as long as they are the same voxel size as the tkreg vox2ras
    // does not use individual information.
    if (gcam->mri_xind == NULL)
    {
      printf("ERROR: GCAsampleMorphRAS(): gcam not inverted\n");
      return(1);
    }
    vox2ras = MRIxfmCRS2XYZtkreg(gcam->mri_xind);
    ras2vox = MatrixInverse(vox2ras,ras2vox);
    ras = MatrixAlloc(4,1,MATRIX_REAL);
    ras->rptr[4][1] = 1;
    crs= MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
    init_needed = 0;
  }

  ras->rptr[1][1] = xMorph;
  ras->rptr[2][1] = yMorph;
  ras->rptr[3][1] = zMorph;
  crs = MatrixMultiply(ras2vox,ras,crs);
  cMorph = crs->rptr[1][1];
  rMorph = crs->rptr[2][1];
  sMorph = crs->rptr[3][1];

  err = GCAMsampleMorph(gcam, cMorph, rMorph, sMorph, &cAnat, &rAnat, &sAnat);
  if (err)
  {
    return(1);
  }

  crs->rptr[1][1] = cAnat;
  crs->rptr[2][1] = rAnat;
  crs->rptr[3][1] = sAnat;
  ras = MatrixMultiply(vox2ras,crs,ras);

  *xAnat = ras->rptr[1][1];
  *yAnat = ras->rptr[2][1];
  *zAnat = ras->rptr[3][1];

  return(0);
}
/*----------------------------------------------------------------------------
  GCAMsampleInverseMorphRAS() - given an RAS coord in the Anat space
  (xyzAnat), compute the RAS in the Morph space (xyzMorph).  The Anat
  space RAS is in "tkregister" or "surface" space (ie, cras=0). See
  also GCAMsampleMorphRAS().
  --------------------------------------------------------------------------*/
int GCAMsampleInverseMorphRAS(GCA_MORPH *gcam,
                              float xAnat, float yAnat, float zAnat,
                              float *xMorph, float *yMorph, float *zMorph)
{
  static int init_needed = 1, err;
  static MATRIX *ras=NULL, *crs=NULL, *vox2ras=NULL, *ras2vox=NULL;
  float  cMorph=0., rMorph=0., sMorph=0.;
  float  cAnat, rAnat, sAnat;

  if (init_needed)
  {
    // Use static so that all this matrix overhead does not have to be
    // done on each call.  Init tkreg vox2ras with the inverse volume
    // as this will be the same as the anatomical volume. The only
    // time this will fail is if this function is called twice with
    // two different volume sizes (eg, 1mm^3 res and a high
    // res). Multiple calls with different subjects is not a problem
    // as long as they are the same voxel size as the tkreg vox2ras
    // does not use individual information.
    if (gcam->mri_xind == NULL)
    {
      printf("ERROR: GCAsampleMorphRAS(): gcam not inverted\n");
      return(1);
    }
    vox2ras = MRIxfmCRS2XYZtkreg(gcam->mri_xind);
    ras2vox = MatrixInverse(vox2ras,ras2vox);
    ras = MatrixAlloc(4,1,MATRIX_REAL);
    ras->rptr[4][1] = 1;
    crs= MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
    init_needed = 0;
  }

  ras->rptr[1][1] = xAnat;
  ras->rptr[2][1] = yAnat;
  ras->rptr[3][1] = zAnat;
  crs = MatrixMultiply(ras2vox,ras,crs);
  cAnat = crs->rptr[1][1];
  rAnat = crs->rptr[2][1];
  sAnat = crs->rptr[3][1];

  err = GCAMsampleInverseMorph(gcam, cAnat, rAnat, sAnat,
                               &cMorph, &rMorph, &sMorph);
  if (err)
  {
    return(1);
  }

  crs->rptr[1][1] = cMorph;
  crs->rptr[2][1] = rMorph;
  crs->rptr[3][1] = sMorph;
  ras = MatrixMultiply(vox2ras,crs,ras);

  *xMorph = ras->rptr[1][1];
  *yMorph = ras->rptr[2][1];
  *zMorph = ras->rptr[3][1];

  return(0);
}
/*-----------------------------------------------------------------------
  GCAMmorphSurf() - compute the vertex xyz in the morph space. Replaces
  vertex xyz. Does not recompute the metric properties of the surface.
  ---------------------------------------------------------------------*/
int GCAMmorphSurf(MRIS *mris, GCA_MORPH *gcam)
{
  int vtxno,err;
  VERTEX *v;
  float Mx, My, Mz;

  //printf("Appling Inverse Morph \n");
  for (vtxno = 0; vtxno < mris->nvertices; vtxno++)
  {
    v = &(mris->vertices[vtxno]);
    err = GCAMsampleInverseMorphRAS(gcam, v->x, v->y, v->z, &Mx, &My, &Mz);
    if (err)
    {
      printf("WARNING: GCAMmorphSurf(): error converting vertex %d\n",vtxno);
      printf("  Avxyz = (%g,%g,%g), Mvxyz = (%g,%g,%g), \n",
             v->x, v->y, v->z, Mx, My, Mz);
      printf(" ... Continuing\n");
    }
    // pack it back into the vertex
    v->x = Mx;
    v->y = My;
    v->z = Mz;
  }
  return(0);
}

double
GCAMcomputeRMS(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  double rms;

#ifdef GCAM_COMPUTE_RMS_GPU
  rms = gcamComputeRMSonGPU( gcam, mri, parms );
#else
  float   nvoxels ;
  double  sse ;

  check_gcam(gcam) ;
  sse = gcamComputeSSE(gcam, mri, parms) ;
  check_gcam(gcam) ;
  nvoxels = gcam->width*gcam->height*gcam->depth ;
  rms = sqrt(sse/nvoxels) ;
#endif

  return(rms) ;
}

#define BIN_SCALE 1
double
gcamComputeSSE(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  double sse;

  double elastic_sse=0 ;
  double ms_sse, l_sse, s_sse, ls_sse, j_sse, d_sse, a_sse;
  double nvox, label_sse, map_sse, dtrans_sse;
  double binary_sse, area_intensity_sse, spring_sse, exp_sse;

  if (!DZERO(parms->l_area_intensity))
  {
    parms->nlt = gcamCreateNodeLookupTable(gcam, parms->mri, parms->nlt) ;
  }

  nvox = gcam->width*gcam->height*gcam->depth ;
  area_intensity_sse =
    binary_sse = spring_sse = ls_sse = exp_sse = dtrans_sse =
                                         label_sse = map_sse = a_sse = sse = ms_sse = l_sse = s_sse = j_sse = d_sse = 0.0;

  check_gcam(gcam) ;
  gcamComputeMetricProperties(gcam) ;
  check_gcam(gcam) ;
  if (!DZERO(parms->l_log_likelihood) || !DZERO(parms->l_likelihood))
    l_sse = MAX(parms->l_log_likelihood, parms->l_likelihood) *
            gcamLogLikelihoodEnergy(gcam, mri) ;
  if (!DZERO(parms->l_multiscale))
  {
    ms_sse = parms->l_multiscale * gcamMultiscaleEnergy(gcam, mri) ;
  }

  if (!DZERO(parms->l_dtrans))
  {
    dtrans_sse = gcamDistanceTransformEnergy(gcam, mri, parms->mri_dist_map, parms) ;
  }
  check_gcam(gcam) ;
  if (!DZERO(parms->l_label))
    label_sse = parms->l_label *
                gcamLabelEnergy(gcam, mri, parms->label_dist) ;
  if (!DZERO(parms->l_binary))
  {
    binary_sse = parms->l_binary * gcamBinaryEnergy(gcam, parms->mri_binary) ;
  }
  if (!DZERO(parms->l_area_intensity))
    area_intensity_sse =
      parms->l_area_intensity *
      gcamAreaIntensityEnergy(gcam, parms->mri, parms->nlt) ;

  check_gcam(gcam) ;
  if (!DZERO(parms->l_map))
  {
    map_sse = parms->l_map * gcamMapEnergy(gcam, mri) ;
  }
  if (!DZERO(parms->l_expansion))
  {
    exp_sse = parms->l_expansion * gcamExpansionEnergy(gcam, mri) ;
  }
  if (!DZERO(parms->l_distance))
  {
    d_sse = parms->l_distance * gcamDistanceEnergy(gcam, mri) ;
  }
  if (!DZERO(parms->l_jacobian))
  {
    j_sse = parms->l_jacobian * gcamJacobianEnergy(gcam, mri) ;
  }
  if (!DZERO(parms->l_area))
  {
    a_sse = parms->l_area * gcamAreaEnergy(gcam) ;
  }
  if (!DZERO(parms->l_area_smoothness))
  {
    a_sse = parms->l_area_smoothness * gcamAreaEnergy(gcam) ;
  }
  if (!DZERO(parms->l_smoothness))
  {
    s_sse = parms->l_smoothness * gcamSmoothnessEnergy(gcam, mri) ;
  }
  if (!DZERO(parms->l_lsmoothness))
  {
    ls_sse = parms->l_lsmoothness * gcamLSmoothnessEnergy(gcam, mri) ;
  }
  if (!DZERO(parms->l_spring))
  {
    spring_sse = parms->l_spring * gcamSpringEnergy(gcam, parms->ratio_thresh);
  }
  if (!DZERO(parms->l_elastic))
  {
    elastic_sse = parms->l_elastic * gcamElasticEnergy(gcam, parms);
  }

  if (Gdiag & DIAG_SHOW)
  {
    printf("\t");
    if (!DZERO(parms->l_elastic))
    {
      printf("elastic_sse = %2.2f ", elastic_sse/nvox) ;
    }
    if (!DZERO(parms->l_map))
    {
      printf("map_sse = %2.2f ", map_sse/nvox) ;
    }
    if (!DZERO(parms->l_spring))
    {
      printf("spring_sse = %2.2f ", spring_sse/nvox) ;
    }
    if (!DZERO(parms->l_dtrans))
    {
      printf("dtrans_rms = %2.4f ", sqrt(dtrans_sse/nvox)) ;
    }
    if (!DZERO(parms->l_log_likelihood))
    {
      printf("l_rms = %2.2f ", sqrt(l_sse/nvox)) ;
    }
    if (!DZERO(parms->l_multiscale))
    {
      printf("ms_rms = %2.2f ", sqrt(ms_sse/nvox)) ;
    }
    if (!DZERO(parms->l_expansion))
    {
      printf("l_exp = %2.2f ", sqrt(exp_sse/nvox)) ;
    }
    if (!DZERO(parms->l_binary))
    {
      printf("b_rms = %2.3f ", sqrt(binary_sse/(BIN_SCALE*nvox))) ;
    }
    if (!DZERO(parms->l_area_intensity))
    {
      printf("ai_rms = %2.3f ", sqrt(area_intensity_sse/nvox)) ;
    }
    if (!DZERO(parms->l_distance))
    {
      printf("d_sse = %2.2f ", d_sse/nvox) ;
    }
    if (!DZERO(parms->l_jacobian))
    {
      printf("j_rms = %2.4f ", sqrt(j_sse/nvox)) ;
    }
    if (!DZERO(parms->l_smoothness))
    {
      printf("s_rms = %2.4f ", sqrt(s_sse/nvox)) ;
    }
    if (!DZERO(parms->l_area) || !DZERO(parms->l_area_smoothness))
    {
      printf("a_sse = %2.2f ", a_sse/nvox) ;
    }
    printf("\n") ;
  }
  sse =
    spring_sse+area_intensity_sse+binary_sse+l_sse + ms_sse +
    s_sse + j_sse + d_sse + a_sse + label_sse + map_sse + exp_sse + dtrans_sse + elastic_sse ;
  return(sse) ;
}

int
gcamComputeGradient
(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, GCA_MORPH_PARMS *parms)
{

#ifdef GCAM_COMPUTE_GRADIENT_GPU
  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamComputeGradientGPU( gcam, mri, mri_smooth, parms );
#else
  static int i = 0 ;

  // make dx = dy = 0
  gcamClearGradient(gcam) ;
  gcamComputeMetricProperties(gcam) ;
  gcamMapTerm(gcam, mri, mri_smooth, parms->l_map)  ;
  gcamLabelTerm(gcam, mri, parms->l_label, parms->label_dist, parms->mri_twm)  ;
  gcamAreaIntensityTerm(gcam, mri, mri_smooth,
                        parms->l_area_intensity, parms->nlt, parms->sigma) ;
  gcamBinaryTerm(gcam, parms->mri_binary, parms->mri_binary_smooth,
                 parms->mri_dist_map, parms->l_binary) ;
  gcamExpansionTerm(gcam, mri, parms->l_expansion)  ;
  gcamLikelihoodTerm(gcam, mri, mri_smooth, parms->l_likelihood, parms)  ;
  gcamDistanceTransformTerm(gcam, mri, parms->mri_dist_map, parms->l_dtrans, parms)  ;
  gcamLogLikelihoodTerm(gcam, mri, mri_smooth, parms->l_log_likelihood)  ;
  gcamMultiscaleTerm(gcam, mri, mri_smooth, parms->l_multiscale)  ;
  gcamDistanceTerm(gcam, mri, parms->l_distance)  ;
  gcamElasticTerm(gcam, parms)  ;
  gcamAreaSmoothnessTerm(gcam, mri_smooth, parms->l_area_smoothness)  ;
  gcamAreaTerm(gcam, parms->l_area)  ;
  gcamSmoothnessTerm(gcam, mri, parms->l_smoothness)  ;
  gcamLSmoothnessTerm(gcam, mri, parms->l_lsmoothness)  ;
  gcamSpringTerm(gcam, parms->l_spring, parms->ratio_thresh)  ;
  //  gcamInvalidSpringTerm(gcam, 1.0)  ;
  //
  gcamJacobianTerm(gcam, mri, parms->l_jacobian,parms->ratio_thresh)  ;
  // The following appears to be a null operation, based on current #ifdefs
  gcamLimitGradientMagnitude(gcam, parms, mri) ;


  if (i == Gdiag_no)
  {
    DiagBreak() ;
  }
  if (parms->write_iterations > 0 && \
      (Gdiag & DIAG_WRITE) && \
      getenv("GCAM_YGRAD") != NULL && \
      (parms->l_label > 0))
  {
    char fname[STRLEN] ;
    MRI *mri_grad ;

#if 0
    mri_grad = GCAMwriteMRI(gcam, NULL, GCAM_Y_GRAD) ;
#else
    mri_grad = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_Y_GRAD, 0, 0);
#endif
    sprintf(fname, "%s_ygrad_before_%4.4d.mgz", parms->base_name, i) ;
    printf("writing y gradient to %s...\n", fname) ;
    MRIwrite(mri_grad, fname) ;
    MRIfree(&mri_grad) ;

  }

  if ((Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) || (gcam_write_grad && i<10) || (gcam_write_grad>1))
  {
    char fname[STRLEN] ;
    MRI *mri ;
    printf("writing gradients to ...%s_d[xyz]b_%4.4d.mgz\n", parms->base_name, i) ;
    mri = GCAMwriteMRI(gcam, NULL, GCAM_X_GRAD) ;
    sprintf(fname, "%s_dxb_%4.4d.mgz", parms->base_name, i) ;
    MRIwrite(mri, fname) ;
    sprintf(fname, "%s_dyb_%4.4d.mgz", parms->base_name, i) ;
    GCAMwriteMRI(gcam, mri, GCAM_Y_GRAD) ;
    MRIwrite(mri, fname);
    sprintf(fname, "%s_dzb_%4.4d.mgz", parms->base_name, i) ;
    GCAMwriteMRI(gcam, mri, GCAM_Z_GRAD) ;
    MRIwrite(mri, fname) ;
    MRIfree(&mri) ;
    if (i == 0)
    {
      MRIwrite(mri_smooth, "s.mgz") ;
    }
  }

  gcamSmoothGradient(gcam, parms->navgs) ;
  fix_borders(gcam) ;

  if (parms->write_iterations > 0 &&
      (Gdiag & DIAG_WRITE) &&
      getenv("GCAM_YGRAD_AFTER") != NULL)
  {
    char fname[STRLEN] ;
    MRI *mri_grad ;

#if 0
    mri_grad = GCAMwriteMRI(gcam, NULL, GCAM_Y_GRAD) ;
#else
    mri_grad = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_Y_GRAD, 0, 0);
#endif
    sprintf(fname, "%s_ygrad_after_%4.4d.mgz", parms->base_name, i) ;
    printf("writing smoothed y gradient to %s...\n", fname) ;
    MRIwrite(mri_grad, fname) ;
    MRIfree(&mri_grad) ;

  }


  if ((Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) || (gcam_write_grad && i<10) || (gcam_write_grad>1))
  {
    MRI *mri ;
    char fname[STRLEN] ;
    printf("writing gradients to ...%s_d[xyz]_%4.4d.mgz\n", parms->base_name, i) ;
    mri = GCAMwriteMRI(gcam, NULL, GCAM_X_GRAD) ;
    sprintf(fname, "%s_dx_%4.4d.mgz", parms->base_name, i) ;
    MRIwrite(mri, fname) ;
    sprintf(fname, "%s_dy_%4.4d.mgz", parms->base_name, i) ;
    GCAMwriteMRI(gcam, mri, GCAM_Y_GRAD) ;
    MRIwrite(mri, fname) ;
    sprintf(fname, "%s_dz_%4.4d.mgz", parms->base_name, i) ;
    GCAMwriteMRI(gcam, mri, GCAM_Z_GRAD) ;
    MRIwrite(mri, fname) ;
    MRIfree(&mri) ;
  }
  i++ ;  /* for debugging */
#endif

  return(NO_ERROR) ;
}

double
gcamMaxGradient(GCA_MORPH *gcam)
{
  int            x, y, z ;
  double         dx, dy, dz, norm, max_grad ;
  GCA_MORPH_NODE *gcamn ;

  max_grad = 0.0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        dx = (gcamn->dx) ;
        dy = (gcamn->dy) ;
        dz = (gcamn->dz) ;
        norm = sqrt(dx*dx+dy*dy+dz*dz) ;
        if (norm > max_grad)
        {
          max_grad = norm ;
        }
      }
  return(max_grad) ;
}

int
gcamLimitGradientMagnitude(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms, MRI *mri)
{
  /*!
    This appears to be a null operation, according to the #ifdefs
    both here and in gcamJacobianTerm
  */
  int            x, y, z, xmax, ymax, zmax ;
  double         norm, max_norm, dx, dy, dz /*, scale*/ ;
  GCA_MORPH_NODE *gcamn ;
  float          dt ;

  dt = parms->dt ;
  max_norm = 0.0 ;
  xmax = ymax = zmax = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;

#if 0
        if (FZERO(gcamn->orig_area))
        {
          gcamn->dx = gcamn->dy = gcamn->dz = 0 ;
        }
#endif
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        // dt*length of gradient
        dx = (gcamn->dx + gcamn->jx) ;
        dy = (gcamn->dy + gcamn->jy) ;
        dz = (gcamn->dz + gcamn->jz) ;
        norm = sqrt(dx*dx+dy*dy+dz*dz) ;
#if 0
        if (dt*norm > 3*parms->max_grad)
        {
          scale = 3*parms->max_grad / (dt*norm) ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          gcamn->dx *= scale ;
          gcamn->dy *= scale ;
          gcamn->dz *= scale ;
          norm =
            dt*sqrt(gcamn->dx*gcamn->dx+gcamn->dy*gcamn->dy+gcamn->dz*gcamn->dz) ;
        }
#endif
        // get max norm and its position
        if (norm > max_norm)
        {
          max_norm = norm ;
          xmax = x ;
          ymax = y ;
          zmax = z ;
        }
	if (norm > parms->max_grad)
	{
	  double scale = parms->max_grad / norm ;
          gcamn->dx *= scale ;
          gcamn->dy *= scale ;
          gcamn->dz *= scale ;
	  if (x == Gx && y == Gy && z == Gz)
	    printf("scaling gradient @ (%d, %d, %d) --> (%2.2f, %2.2f, %2.2f)\n", x, y, z, gcamn->dx, gcamn->dy, gcamn->dz) ;
	}
      }
  {
    float vals[MAX_GCA_INPUTS] ;
    int   r ;
#if 0
    int memoryUsed = 0;
#endif
    gcamn = &gcam->nodes[xmax][ymax][zmax] ;
    // print the info at this position
    maxGradient=max_norm;
    gradientArea=gcam->nodes[xmax][ymax][zmax].area;
    if (Gdiag & DIAG_SHOW)
    {
      printf("   max grad %2.3f mm @ (%d, %d, %d), "
             "Area=%2.3f, new/orig=%2.3f, l=%s ",
             max_norm, xmax, ymax, zmax,
             gcam->nodes[xmax][ymax][zmax].area,
             gcam->nodes[xmax][ymax][zmax].area/gcam->nodes[xmax][ymax][zmax].orig_area,
             cma_label_to_name(gcam->nodes[xmax][ymax][zmax].label));
      fflush(stdout);
      if (parms->l_dtrans > 0)
      {
	double dist = -10000; ;
        int frame = dtrans_label_to_frame(parms, gcamn->label) ;
	if (frame >= 0)
	  MRIsampleVolumeFrameType(parms->mri_dist_map, gcamn->x, gcamn->y, gcamn->z, frame, SAMPLE_TRILINEAR,&dist) ;

        printf("dist:target = %2.2f:%2.2f ", dist, gcamn->target_dist) ;
      }
      else
      {
        printf("vals(means) = ") ;
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->ninputs) ;
        for (r = 0 ; r < gcam->ninputs ; r++)
          printf("%2.1f (%2.1f)  ",
                 vals[r], gcamn->gc ? gcamn->gc->means[r] :-0.0);
      }
      fflush(stdout);
#if 0
      if (memoryUsed=getMemoryUsed() != -1)
      {
        printf("memory used: %d Kbytes\n", getMemoryUsed());
      }
#endif
      printf("\n") ;
    }
  }

#if 0
  if (max_norm > parms->max_grad)
  {
    scale = parms->max_grad / max_norm ;

    printf("scaling by %2.3f based on max gradient %2.3f mm @ (%d, %d, %d)\n",
           scale, max_norm, xmax, ymax, zmax) ;
    for (x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;

          if (gcamn->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }

          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          gcamn->dx *= scale ;
          gcamn->dy *= scale ;
          gcamn->dz *= scale ;
        }
  }
#endif
  return(NO_ERROR) ;
}


int
gcamSmoothnessTerm( GCA_MORPH *gcam,
                    const MRI *mri,
                    const double l_smoothness)
{

#ifdef GCAM_SMOOTH_TERM_GPU
  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamSmoothnessTermGPU( gcam, l_smoothness );
#else
  double          vx=0.0, vy=0.0, vz=0.0, vnx=0.0, vny=0.0, vnz=0.0;
  double          dx=0.0, dy=0.0, dz=0.0 ;
  int             x=0, y=0, z=0, xk=0, yk=0, zk=0, xn=0, yn=0, zn=0 ;
  int             width, height, depth, num=0 ;
  GCA_MORPH_NODE  *gcamn=NULL, *gcamn_nbr=NULL ;

  if (DZERO(l_smoothness))
  {
    return(NO_ERROR) ;
  }
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate (y,z,gcamn,vx,vy,vz,dx,dy,dz,num,xk,xn,yk,yn,zk,zn,gcamn_nbr,vnx,vny,vnz) shared(gcam,Gx,Gy,Gz) schedule(static,1)
#endif
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        vx = gcamn->x - gcamn->origx ;
        vy = gcamn->y - gcamn->origy ;
        vz = gcamn->z - gcamn->origz ;
        dx = dy = dz = 0.0f ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_smoo: node(%d,%d,%d): V=(%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, vx, vy, vz) ;
        num = 0 ;

        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          xn = MAX(0,xn) ;
          xn = MIN(width-1,xn) ;

          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            yn = MAX(0,yn) ;
            yn = MIN(height-1,yn) ;

            for (zk = -1 ; zk <= 1 ; zk++)
            {

              if (!zk && !yk && !xk)
              {
                continue ;
              }

              zn = z+zk ;
              zn = MAX(0,zn) ;
              zn = MIN(depth-1,zn) ;

              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;

              if (gcamn_nbr->invalid == GCAM_POSITION_INVALID)
              {
                continue;
              }

#if 0
              if (gcamn_nbr->label != gcamn->label)
              {
                continue ;
              }
#endif

              vnx = gcamn_nbr->x - gcamn_nbr->origx ;
              vny = gcamn_nbr->y - gcamn_nbr->origy ;
              vnz = gcamn_nbr->z - gcamn_nbr->origz ;

              dx += (vnx-vx) ;
              dy += (vny-vy) ;
              dz += (vnz-vz) ;

              if ((x == Gx && y == Gy && z == Gz) &&
                  (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
              {
                printf("\tnode(%d,%d,%d): V=(%2.2f,%2.2f,%2.2f), "
                       "DX=(%2.2f,%2.2f,%2.2f)\n",
                       xn, yn, zn, vnx, vny, vnz, vnx-vx, vny-vy, vnz-vz) ;
              }

              num++ ;
            }
          }
        }
        /*        num = 1 ;*/
        if (num)
        {
          dx = dx * l_smoothness / num ;
          dy = dy * l_smoothness / num ;
          dz = dz * l_smoothness / num ;
        }

        if (x == Gx && y == Gy && z == Gz)
        {
          printf("l_smoo: node(%d,%d,%d): DX=(%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, dx, dy, dz) ;
        }

        gcamn->dx += dx ;
        gcamn->dy += dy ;
        gcamn->dz += dz ;

      }
    }
  }
#endif
  return(NO_ERROR) ;
}



#define GCAM_SMOOTHNESS_OUTPUT 0

double
gcamElasticEnergy( const GCA_MORPH *gcam, GCA_MORPH_PARMS *parms )
{
  MRI    *mri_warp, *mri_sobel[3] ;
  int    x, y, z, dim1, dim2 ;
  double energy, val1, val2, lambda, mu, rigid_energy, volume_change, svals[3][3] ;

  if (FZERO(parms->l_elastic))
  {
    return(0.0) ;
  }
  lambda = 0.57692 ;
  mu = 0.38462 ; // matches CVS
  mri_warp = GCAMwriteWarpToMRI(gcam, NULL) ;
  mri_sobel[0] = MRIsobelFrame(mri_warp, NULL, NULL, 0) ;
  mri_sobel[1] = MRIsobelFrame(mri_warp, NULL, NULL, 1) ;
  mri_sobel[2] = MRIsobelFrame(mri_warp, NULL, NULL, 2) ;

  for (energy = 0.0, x = 0 ; x < mri_warp->width ; x++)
    for (y = 0 ; y < mri_warp->height ; y++)
      for (z = 0 ; z < mri_warp->depth ; z++)
      {
        for (dim1 = 0 ; dim1 < 3 ; dim1++)
          for (dim2 = 0 ; dim2 < 3 ; dim2++)
          {
            svals[dim1][dim2] = MRIgetVoxVal(mri_sobel[dim1], x, y, z, dim2) ;
          }

        rigid_energy = volume_change = 0 ;
        for (dim1 = 0 ; dim1 < 3 ; dim1++)
          for (dim2 = 0 ; dim2 < 3 ; dim2++)
          {
            val1 = svals[dim1][dim1] ;
            val2 = svals[dim2][dim2] ;
            rigid_energy += val1*val2 ;

            val1 = svals[dim1][dim2] ;
            val2 = svals[dim2][dim1] ;
            volume_change += SQR(val1 + val2);
          }
        energy += rigid_energy * mu/4 + volume_change * lambda/2 ;
      }

  for (dim1 = 0 ; dim1 < 3 ; dim1++)
  {
    MRIfree(&mri_sobel[dim1]) ;
  }
  MRIfree(&mri_warp) ;
  return(energy) ;
}

int
gcamElasticTerm( const GCA_MORPH *gcam, GCA_MORPH_PARMS *parms )
{
  MRI    *mri_warp, *mri_tmp, *mri_divergence, *mri_laplacian, *mri_grad, *mri_kernel ;
  int    x, y, z, dim ;
  double val1, val2, grad[3], lambda, mu ;
  GCA_MORPH_NODE  *gcamn ;

  if (DZERO(parms->l_elastic))
  {
    return(NO_ERROR) ;
  }

  mri_kernel = MRIgaussian1d(parms->sigma, 100) ;
  lambda = parms->lame_lambda ;
  mu = parms->lame_mu ; // matches CVS
  mri_warp = GCAMwriteWarpToMRI(gcam, NULL) ;

  mri_tmp = MRIdivergence(mri_warp, NULL) ;
  mri_divergence = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
  MRIfree(&mri_tmp) ;
  mri_grad = MRIsobel(mri_divergence, NULL, NULL) ;      // gradient of the divergence
  MRIfree(&mri_divergence) ;
  mri_tmp = MRIlaplacian(mri_warp, NULL) ;
  mri_laplacian = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
  MRIfree(&mri_kernel) ;
  MRIfree(&mri_tmp) ;

  for (x = 0 ; x < mri_warp->width ; x++)
    for (y = 0 ; y < mri_warp->height ; y++)
      for (z = 0 ; z < mri_warp->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        for (dim = 0 ; dim < 3 ; dim++)
        {
          val1 = MRIgetVoxVal(mri_grad, x, y, z, dim) ;
          val2 = MRIgetVoxVal(mri_laplacian, x, y, z, dim) ;
          grad[dim] = val1 * (lambda + mu) + val2 * mu ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("node(%d, %d, %d), elastic term = (%2.2f, %2.2f, %2.2f)\n",
                 x, y, z, parms->l_elastic * grad[0], parms->l_elastic * grad[1], parms->l_elastic * grad[2]) ;
        gcamn->dx += parms->l_elastic * grad[0] ;
        gcamn->dy += parms->l_elastic * grad[1] ;
        gcamn->dz += parms->l_elastic * grad[2] ;
      }

  MRIfree(&mri_warp) ;
  MRIfree(&mri_laplacian) ;
  MRIfree(&mri_grad) ;
  return(NO_ERROR) ;
}

double
gcamSmoothnessEnergy( const GCA_MORPH *gcam, const MRI *mri )
{
  /*!
    Computes a load of derivatives (of some description).
    I don't know the purpose of the MRI argument, since it's
    not actually used by anything
  */
  double sse = 0.0;

#if GCAM_SMOOTHNESS_OUTPUT
  const unsigned int gcamSmoothnessOutputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%gcamSmoothnessOutputFreq)==0 )
  {
    char fname[STRLEN];
    snprintf( fname, STRLEN-1,
              "gcamSmoothInput%04u", nCalls/gcamSmoothnessOutputFreq );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );
  }
  nCalls++;
#endif

#ifdef GCAM_SMOOTHNESS_ENERGY_GPU
#if SHOW_EXEC_LOC
  printf( "%s: CUDA call\n", __FUNCTION__ );
#endif

  sse = gcamSmoothnessEnergyGPU( gcam );
#else
#if SHOW_EXEC_LOC
  printf( "%s: CPU call\n", __FUNCTION__ );
#endif
  double vx=0, vy=0, vz=0, vnx=0, vny=0, vnz=0, error=0, node_sse=0;
  double dx=0, dy=0, dz=0 ;
  int x=0, y=0, z=0, xk=0, yk=0, zk=0, xn=0, yn=0, zn=0 ;
  int width=0, height=0, depth=0, num=0 ;
  const GCA_MORPH_NODE  *gcamn=NULL, *gcamn_nbr=NULL ;

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;

#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate(y,z,vx,vy,vz,vnx,vny,vnz,error,node_sse,dx,dy,dz,xk,yk,zk,xn,yn,zn,num,gcamn,gcamn_nbr) reduction(+:sse) shared(gcam,Gx,Gy,Gz) schedule(static,1)
#endif
  // Loop over all voxels
  for (x = 0 ; x < gcam->width ; x++ )
  {
    for (y = 0 ; y < gcam->height ; y++ )
    {
      for (z = 0 ; z < gcam->depth ; z++ )
      {

        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        // Compute differences from original
        vx = gcamn->x - gcamn->origx ;
        vy = gcamn->y - gcamn->origy ;
        vz = gcamn->z - gcamn->origz ;
        num = 0 ;
        node_sse = 0.0 ;

        // Loop over 3^3 voxels centred on current
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          xn = MAX(0,xn) ;
          xn = MIN(width-1,xn) ;

          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            yn = MAX(0,yn) ;
            yn = MIN(height-1,yn) ;

            for (zk = -1 ; zk <= 1 ; zk++)
            {

              // Don't use self
              if (!xk && !yk && !zk)
              {
                continue ;
              }

              zn = z+zk ;
              zn = MAX(0,zn) ;
              zn = MIN(depth-1,zn) ;

              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;

              if (gcamn_nbr->invalid == GCAM_POSITION_INVALID)
              {
                continue;
              }
#if 0
              if (gcamn_nbr->label != gcamn->label)
              {
                continue ;
              }
#endif
              vnx = gcamn_nbr->x - gcamn_nbr->origx ;
              vny = gcamn_nbr->y - gcamn_nbr->origy ;
              vnz = gcamn_nbr->z - gcamn_nbr->origz ;

              dx = vnx-vx ;
              dy = vny-vy ;
              dz = vnz - vz ;

              error = dx*dx + dy*dy + dz*dz ;

              num++ ;
              node_sse += error ;
            }
          }
        }

        /*        num = 1 ;*/
        if (num > 0)
        {
          sse += node_sse/num ;
        }

        if (x == Gx && y == Gy && z == Gz)
        {
          printf("E_smoo: node(%d,%d,%d) smoothness sse %2.3f (%d nbrs)\n",
                 x, y, z, node_sse/num, num) ;
        }
      }
    }
  }
#endif

  return(sse) ;
}


int
gcamLSmoothnessTerm(GCA_MORPH *gcam, MRI *mri, double l_smoothness)
{
  double          vx, vy, vz, vnx, vny, vnz, dx, dy, dz ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  if (DZERO(l_smoothness))
  {
    return(NO_ERROR) ;
  }
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        vx = gcamn->x - gcamn->origx ;
        vy = gcamn->y - gcamn->origy ;
        vz = gcamn->z - gcamn->origz ;
        dx = dy = dz = 0.0f ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_smoo: node(%d,%d,%d): V=(%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, vx, vy, vz) ;
        num = 0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          xn = MAX(0,xn) ;
          xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            yn = MAX(0,yn) ;
            yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!zk && !yk && !xk)
              {
                continue ;
              }
              zn = z+zk ;
              zn = MAX(0,zn) ;
              zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;

              if (gcamn_nbr->invalid/* == GCAM_POSITION_INVALID*/)
              {
                continue;
              }

              if (gcamn_nbr->label != gcamn->label)
              {
                continue ;
              }
              vnx = gcamn_nbr->x - gcamn_nbr->origx ;
              vny = gcamn_nbr->y - gcamn_nbr->origy ;
              vnz = gcamn_nbr->z - gcamn_nbr->origz ;
              dx += (vnx-vx) ;
              dy += (vny-vy) ;
              dz += (vnz-vz) ;
              if ((x == Gx && y == Gy && z == Gz) &&
                  (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
                printf("\tnode(%d,%d,%d): V=(%2.2f,%2.2f,%2.2f), "
                       "DX=(%2.2f,%2.2f,%2.2f)\n",
                       xn, yn, zn, vnx, vny, vnz, vnx-vx, vny-vy, vnz-vz) ;
              num++ ;
            }
          }
        }
        /*        num = 1 ;*/
        if (num)
        {
          dx = dx * l_smoothness / num ;
          dy = dy * l_smoothness / num ;
          dz = dz * l_smoothness / num ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("l_smoo: node(%d,%d,%d): DX=(%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, dx, dy, dz) ;
        gcamn->dx += dx ;
        gcamn->dy += dy ;
        gcamn->dz += dz ;
      }
  return(NO_ERROR) ;
}


double
gcamLSmoothnessEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double sse = 0.0, vx, vy, vz, vnx, vny, vnz, error, node_sse, dx, dy, dz ;
  int x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        vx = gcamn->x - gcamn->origx ;
        vy = gcamn->y - gcamn->origy ;
        vz = gcamn->z - gcamn->origz ;
        num = 0 ;
        node_sse = 0.0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          xn = MAX(0,xn) ;
          xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            yn = MAX(0,yn) ;
            yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!xk && !yk && !zk)
              {
                continue ;
              }
              zn = z+zk ;
              zn = MAX(0,zn) ;
              zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;

              if (gcamn_nbr->invalid/* == GCAM_POSITION_INVALID*/)
              {
                continue;
              }
              if (gcamn_nbr->label != gcamn->label)
              {
                continue ;
              }
              vnx = gcamn_nbr->x - gcamn_nbr->origx ;
              vny = gcamn_nbr->y - gcamn_nbr->origy ;
              vnz = gcamn_nbr->z - gcamn_nbr->origz ;
              dx = vnx-vx ;
              dy = vny-vy ;
              dz = vnz - vz ;
              error = dx*dx + dy*dy + dz*dz ;
              num++ ;
              node_sse += error ;
            }
          }
        }
        /*        num = 1 ;*/
        if (num > 0)
        {
          sse += node_sse/num ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("E_smoo: node(%d,%d,%d) smoothness sse %2.3f (%d nbrs)\n",
                 x, y, z, node_sse/num, num) ;
      }
  return(sse) ;
}


int
gcamSpringTerm(GCA_MORPH *gcam, double l_spring, double ratio_thresh)
{
  double          xc, yc, zc, dx, dy, dz ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  if (DZERO(l_spring))
  {
    return(NO_ERROR) ;
  }
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
#if 0
        if (gcamn->area / gcamn->orig_area >= ratio_thresh)
        {
          continue ;
        }
#endif

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        // compute centroid
        num = 0 ;
        xc = yc = zc = 0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          if (xn < 0 || xn >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            if (yn < 0 || yn >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zn = z+zk ;
              if (zn < 0 || zn >= gcam->depth)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;
              if (gcamn_nbr->invalid == GCAM_POSITION_INVALID)
              {
                continue;
              }
              xc += gcamn_nbr->x ;
              yc += gcamn_nbr->y ;
              zc += gcamn_nbr->z ;

              num++ ;
            }
          }
        }
        if (num == 0)
        {
          continue ;
        }
        xc /= (float)num ;
        yc /= (float)num ;
        zc /= (float)num ;

        dx = xc - gcamn->x ;
        dy = yc - gcamn->y ;
        dz = zc - gcamn->z ;
        gcamn->dx += l_spring * dx ;
        gcamn->dy += l_spring * dy ;
        gcamn->dz += l_spring * dz ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_spring: node(%d,%d,%d) "
                 "spring term (%2.3f, %2.3f, %2.3f))\n",
                 x, y, z, l_spring*dx, l_spring*dy, l_spring*dz) ;
      }
  return(NO_ERROR) ;
}
#if 0
int
gcamInvalidSpringTerm(GCA_MORPH *gcam, double l_spring)
{
  double          xc, yc, zc, dx, dy, dz ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  if (DZERO(l_spring))
  {
    return(NO_ERROR) ;
  }
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == 0)
        {
          continue;
        }

        // compute centroid
        num = 0 ;
        xc = yc = zc = 0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          if (xn < 0 || xn >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            if (yn < 0 || yn >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zn = z+zk ;
              if (zn < 0 || zn >= gcam->depth)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;
              if (gcamn_nbr->invalid == GCAM_POSITION_INVALID)
              {
                continue;
              }
              xc += gcamn_nbr->x ;
              yc += gcamn_nbr->y ;
              zc += gcamn_nbr->z ;

              num++ ;
            }
          }
        }
        if (num == 0)
        {
          continue ;
        }
        xc /= (float)num ;
        yc /= (float)num ;
        zc /= (float)num ;

        dx = xc - gcamn->x ;
        dy = yc - gcamn->y ;
        dz = zc - gcamn->z ;
        gcamn->dx += l_spring * dx ;
        gcamn->dy += l_spring * dy ;
        gcamn->dz += l_spring * dz ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_spring: node(%d,%d,%d) spring term "
                 "(%2.3f, %2.3f, %2.3f))\n",
                 x, y, z, l_spring*dx, l_spring*dy, l_spring*dz) ;
      }
  return(NO_ERROR) ;
}
#endif


double
gcamSpringEnergy(GCA_MORPH *gcam, double ratio_thresh)
{
  double          sse = 0.0, xc, yc, zc, node_sse, dx, dy, dz ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->area / gcamn->orig_area >= ratio_thresh)
        {
          continue ;
        }

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        // compute centroid
        num = 0 ;
        xc = yc = zc = 0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ;
          if (xn < 0 || xn >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ;
            if (yn < 0 || yn >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zn = z+zk ;
              if (zn < 0 || zn >= gcam->depth)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;
              if (gcamn_nbr->invalid == GCAM_POSITION_INVALID)
              {
                continue;
              }
              xc += gcamn_nbr->x ;
              yc += gcamn_nbr->y ;
              zc += gcamn_nbr->z ;

              num++ ;
            }
          }
        }
        xc /= (float)num ;
        yc /= (float)num ;
        zc /= (float)num ;
        dx = gcamn->x - xc ;
        dy = gcamn->y - yc ;
        dz = gcamn->z - zc ;
        node_sse = dx*dx + dy*dy + dz*dz ;
        if (num > 0)
        {
          sse += node_sse/num ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("E_spring: node(%d,%d,%d) spring sse %2.3f (%d nbrs)\n",
                 x, y, z, node_sse/num, num) ;
      }
  return(sse) ;
}

int
gcamClearGradient(GCA_MORPH *gcam)
{
  int   x, y, z ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
        gcam->nodes[x][y][z].dx =
          gcam->nodes[x][y][z].dy = gcam->nodes[x][y][z].dz = 0.0;
  return(NO_ERROR) ;
}


int
gcamApplyGradient(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms)
{
  int            x, y, z ;
  float          dx, dy, dz, dt, momentum ;
  GCA_MORPH_NODE *gcamn ;

  dt = parms->dt ;
  momentum = parms->momentum ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        dx = gcamn->dx*dt + gcamn->odx*momentum ;
        dy = gcamn->dy*dt + gcamn->ody*momentum ;
        dz = gcamn->dz*dt + gcamn->odz*momentum ;
        gcamn->odx = dx ;
        gcamn->ody = dy ;
        gcamn->odz = dz ;

        if (x == Gx && y == Gy && z == Gz)
          printf("GRAD: node(%d,%d,%d): moving by (%2.3f, %2.3f, %2.3f) from "
                 "(%2.2f,%2.2f,%2.2f) to ",
                 x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
        gcamn->x += dx;
        gcamn->y += dy;
        gcamn->z += dz;
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("(%2.2f,%2.2f,%2.2f)\n", gcamn->x, gcamn->y, gcamn->z) ;
        }
      }
  if (!DZERO(parms->l_area_intensity))
  {
    parms->nlt = gcamCreateNodeLookupTable(gcam, parms->mri, parms->nlt) ;
  }
  return(NO_ERROR) ;
}


int
gcamClearMomentum(GCA_MORPH *gcam)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        gcamn->odx = gcamn->ody = gcamn->odz = 0;

      }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
finitep(float f)
{
  if (!finite(f))
  {
    return(0) ;
  }
  if (fabs(f) > 1e5)
  {
    return(1) ;
  }
  return(1) ;
}



int
GCAMcomputeLabels(MRI *mri, GCA_MORPH *gcam)
{
  int            x, y, z, width, height, depth, label, n, nchanged = 0 ;
  float          vals[MAX_GCA_INPUTS] ;
  GCA_MORPH_NODE *gcamn ;
  GCA_PRIOR      *gcap ;
  GC1D           *gc ;

  if (gcam->gca == NULL)
  {
    return(NO_ERROR) ;
  }

  // gcam usually has prior_width, prior_height, prior_depth
  width = gcam->width  ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z];

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        // get the grey scale values from mri at floating point voxel position
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->ninputs) ;
        label =
          GCAcomputeMAPlabelAtLocation(gcam->gca, x,y,z,vals,&n,&gcamn->log_p);
        // n is the most probable index in labels array
        gcap = &gcam->gca->priors[x][y][z] ;

        if (n >= 0)
        {
          if (label != gcamn->label)
          {
            nchanged++ ;
          }
          gcamn->label = label ;
          gcamn->n = n ;
          gcamn->prior = gcap->priors[n] ;
          gcamn->gc = gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;

          if (x == Gx && y == Gy && z == Gz)
            printf
            ("RELABEL: node(%d, %d, %d): label %s (%d), "
             "mean %2.1f+-%2.1f, prior %2.1f, MRI=%2.0f\n",
             x, y, z, cma_label_to_name(label), label,
             gcamn->gc ? gcamn->gc->means[0] : 0.0,
             gcamn->gc ? sqrt(covariance_determinant(gcamn->gc,
                              gcam->ninputs)) : 0.0,
             gcamn->prior,vals[0]) ;
        }
        else  /* out of FOV probably */
        {
          gcamn->label = label ;
          gcamn->n = 0 ;
          gcamn->prior = 1.0 ;
          gcamn->gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;
          if (x == Gx && y == Gy && z == Gz)
            printf("RELABEL: node(%d, %d, %d): label %s (%d), "
                   "mean %2.1f+-%2.1f, prior %2.1f\n",
                   x, y, z, cma_label_to_name(label), label,
                   0.0, 0.0, gcamn->prior) ;
        }
        if (gcamn->gc == NULL)
        {
          int xn, yn, zn ;
          GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn);
          gcamn->gc = GCAfindClosestValidGC(gcam->gca, xn, yn, zn, label,0);
        }
      }

  fprintf(stderr, "label assignment complete, %d changed (%2.2f%%)\n",
         nchanged, 100.0*(float)nchanged/(width*height*depth)) ;

  return(nchanged) ;
}
int
GCAMcomputeMaxPriorLabels(GCA_MORPH *gcam)
{
  int            x, y, z, width, height, depth, label, n, nchanged = 0,max_n ;
  GCA_MORPH_NODE *gcamn ;
  GC1D           *gc ;
  GCA_PRIOR      *gcap ;
  double         max_prior ;

  if (gcam->gca == NULL || gcam->gca->priors == NULL)
  {
    return(NO_ERROR) ;
  }

  width = gcam->width  ;
  height = gcam->height ;
  depth = gcam->depth ;
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        gcap = &gcam->gca->priors[x][y][z] ;
        for (max_prior = -1, max_n = -1, n = 0 ; n < gcap->nlabels ; n++)
          if  (gcap->priors[n] >= max_prior)
          {
            max_prior = gcap->priors[n] ;
            max_n = n ;
          }
        n = max_n ;
        if (n >= 0)
        {
          label = gcap->labels[n] ;
          if (label != gcamn->label)
          {
            nchanged++ ;
          }
          gcamn->label = label ;
          gcamn->n = n ;
          gcamn->prior = gcap->priors[n] ;
          gcamn->gc = gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;
          if (gc == NULL && gcam->gca != NULL) // find a nearby one
          {
            int xn, yn, zn ;
            GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn);
            gcamn->gc = gc = GCAfindClosestValidGC(gcam->gca, xn, yn, zn, label,0);
          }
        }
        else  /* out of FOV probably */
        {
          int xn, yn, zn ;
          GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn);
          gcamn->invalid = GCAM_POSITION_INVALID ;
          gcamn->label = label = 0 ;
          gcamn->n = 0 ;
          gcamn->prior = 1.0 ;
          gcamn->gc = gc = GCAfindClosestValidGC(gcam->gca, xn, yn, zn, label,0); // was NULL
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("MPRIOR: node(%d, %d, %d): label %s (%d), "
                 "mean %2.1f+-%2.1f, prior %2.1f\n",
                 x, y, z, cma_label_to_name(label), label,
                 gcamn->gc ? gcamn->gc->means[0] : 0.0,
                 gcamn->gc ? sqrt(covariance_determinant(gcamn->gc,
                                  gcam->ninputs)) : 0.0,
                 gcamn->prior) ;
      }

  fprintf(stderr, "label assignment complete, %d changed (%2.2f%%)\n",
         nchanged, 100.0*(float)nchanged/(width*height*depth)) ;

  return(NO_ERROR) ;
}
MRI *
GCAMbuildMostLikelyVolume(GCA_MORPH *gcam, MRI *mri)
{
  int            x,  y, z, xn, yn, zn, width, depth, height, n ;
  GCA_MORPH_NODE *gcamn ;
  float          val ;

  // error check
  if (!mri)
    ErrorExit
    (ERROR_BADPARM, "GCAbuildMostLikelyVolume called with null MRI.\n");
  if (mri->width != gcam->atlas.width
      || mri->height != gcam->atlas.height
      || mri->depth != gcam->atlas.depth)
    ErrorExit(ERROR_BADPARM, \
              "GCAbuildMostLikelyVolume called with mri dimension "
              "being different from M3D.\n");

  // set direction cosines etc.
  useVolGeomToMRI(&gcam->atlas, mri);

  width = mri->width ;
  depth = mri->depth ;
  height = mri->height ;

  if (gcam->gca)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          if (!GCAvoxelToPrior(gcam->gca, mri, x, y, z, &xn, &yn, &zn))
          {
            if (xn == Gx && yn == Gy && zn == Gz)
            {
              DiagBreak() ;
            }
            gcamn = &gcam->nodes[xn][yn][zn] ;

            if (gcamn->invalid == GCAM_POSITION_INVALID)
            {
              continue;
            }

            for (n = 0 ; n < gcam->ninputs ; n++)
            {
              if (gcamn->gc)
              {
                val = gcamn->gc->means[n] ;
              }
              else
              {
                val = 0 ;
              }

              switch (mri->type)
              {
              default:
                ErrorReturn
                (NULL,
                 (ERROR_UNSUPPORTED,
                  "GCAMbuildMostLikelyVolume: unsupported image type %d",
                  mri->type)) ;
                break ;
              case MRI_SHORT:
                MRISseq_vox(mri, x, y, z, n) = nint(val) ;
                break ;
              case MRI_UCHAR:
                MRIseq_vox(mri, x, y, z, n) = nint(val) ;
                break ;
              case MRI_FLOAT:
                MRIFseq_vox(mri, x, y, z, n) = val ;
                break ;
              }
            }
          }// !GCA
        }
      }
    }
  }
  else   // no gca
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          gcamn = &gcam->nodes[x][y][z] ;

          if (gcamn->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }

          for (n = 0 ; n < gcam->ninputs ; n++)
          {
            if (gcamn->gc)
            {
              val = gcamn->gc->means[n] ;
            }
            else
            {
              val = 0 ;
            }

            switch (mri->type)
            {
            default:
              ErrorReturn
              (NULL,
               (ERROR_UNSUPPORTED,
                "GCAMbuildMostLikelyVolume: unsupported image type %d",
                mri->type)) ;
              break ;
            case MRI_SHORT:
              MRISseq_vox(mri, x, y, z, n) = nint(val) ;
              break ;
            case MRI_UCHAR:
              MRIseq_vox(mri, x, y, z, n) = nint(val) ;
              break ;
            case MRI_FLOAT:
              MRIFseq_vox(mri, x, y, z, n) = val ;
              break ;
            }
          }
        }
      }
    }
  }


  return(mri) ;
}

MRI *
GCAMbuildVolume(GCA_MORPH *gcam, MRI *mri)
{
  int            x,  y, z, xn, yn, zn, width, depth, height, n ;
  GCA_MORPH_NODE *gcamn ;
  float          val ;

  // error check
  if (!mri)
    ErrorExit
    (ERROR_BADPARM, "GCAbuildMostLikelyVolume called with null MRI.\n");
  if (mri->width != gcam->atlas.width
      || mri->height != gcam->atlas.height
      || mri->depth != gcam->atlas.depth)
    ErrorExit(ERROR_BADPARM, \
              "GCAbuildMostLikelyVolume called with mri dimension "
              "being different from M3D.\n");

  // set direction cosines etc.
  useVolGeomToMRI(&gcam->atlas, mri);

  width = mri->width ;
  depth = mri->depth ;
  height = mri->height ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (!GCAvoxelToPrior(gcam->gca, mri, x, y, z, &xn, &yn, &zn))
        {
          gcamn = &gcam->nodes[xn][yn][zn] ;

          if (gcamn->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }

          for (n = 0 ; n < gcam->ninputs ; n++)
          {
            if (gcamn->gc)
            {
              val = gcamn->gc->means[n] ;
            }
            else
            {
              val = 0 ;
            }

            switch (mri->type)
            {
            default:
              ErrorReturn
              (NULL,
               (ERROR_UNSUPPORTED,
                "GCAMbuildMostLikelyVolume: unsupported image type %d",
                mri->type)) ;
              break ;
            case MRI_SHORT:
              MRISseq_vox(mri, x, y, z, n) = nint(val) ;
              break ;
            case MRI_UCHAR:
              MRIseq_vox(mri, x, y, z, n) = nint(val) ;
              break ;
            case MRI_FLOAT:
              MRIFseq_vox(mri, x, y, z, n) = val ;
              break ;
            }
          }
        } // !GCA
      }
    }
  }

  return(mri) ;
}

int
GCAMinvert(GCA_MORPH *gcam, MRI *mri)
{
  int            x, y, z, width, height, depth, xv, yv, zv ;
  MRI            *mri_ctrl, *mri_counts ;
  GCA_MORPH_NODE *gcamn ;
  double           xf, yf, zf ;
  float          num ;

#if 1
  if (gcam->mri_xind)   /* already inverted */
  {
    return(NO_ERROR) ;
  }
#else
  if (gcam->mri_xind)
  {
    MRIfree(&gcam->mri_xind) ;
  }
  if (gcam->mri_yind)
  {
    MRIfree(&gcam->mri_yind) ;
  }
  if (gcam->mri_zind)
  {
    MRIfree(&gcam->mri_zind) ;
  }
#endif

  // verify the volume size ////////////////////////////////////////////
  if (mri->width != gcam->image.width
      || mri->height != gcam->image.height
      || mri->depth != gcam->image.depth)
    ErrorExit(ERROR_BADPARM, \
              "mri passed volume size ( %d %d %d ) is different from "
              "the one used to create M3D data ( %d %d %d )\n",
              mri->width, mri->height,mri->depth,gcam->image.width,
              gcam->image.height,gcam->image.depth);

  // use mri
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  // mri_xind, yind, zind
  gcam->mri_xind = MRIalloc(width, height, depth, MRI_FLOAT) ;
  MRIcopyHeader(mri, gcam->mri_xind);
  gcam->mri_yind = MRIalloc(width, height, depth, MRI_FLOAT) ;
  MRIcopyHeader(mri, gcam->mri_yind);
  gcam->mri_zind = MRIalloc(width, height, depth, MRI_FLOAT) ;
  MRIcopyHeader(mri, gcam->mri_zind);
  mri_counts = MRIalloc(width, height, depth, MRI_FLOAT) ;
  // mri_ctrl
  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR) ;
  MRIcopyHeader(mri, mri_ctrl);

  if (!gcam->mri_xind || !gcam->mri_yind || !gcam->mri_zind || !mri_ctrl)
    ErrorExit(ERROR_NOMEMORY,
              "GCAMinvert: could not allocated %dx%dx%d index volumes",
              width, height, depth) ;

  // going through gcam volume (x,y,z)
  // gcam volume points could be mapped to many points in xind, yind, and zind
  for (z = 0 ; z < gcam->depth ; z++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (x = 0 ; x < gcam->width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        // find nodes
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        // get the source volume position
        xf = gcamn->x ;
        yf = gcamn->y ;
        zf = gcamn->z ;
        // make them within the range of index /////////////////////////////
        if (xf < 0)
        {
          xf = 0 ;
        }
        if (yf < 0)
        {
          yf = 0 ;
        }
        if (zf < 0)
        {
          zf = 0 ;
        }
        if (xf >= width)
        {
          xf = width-1 ;
        }
        if (yf >= height)
        {
          yf = height-1 ;
        }
        if (zf >= depth)
        {
          zf = depth-1 ;
        }
        xv = nint(xf) ;
        yv = nint(yf) ;
        zv = nint(zf) ;
        ////////////////////////////////////////////////////////////////////
        // cache position
        // xind[xv][yv][zv] += x

        if (abs(xf-Gx)<1 && abs(yf-Gy) <1 && abs(zf-Gz) < 1)
        {
          DiagBreak() ;
        }
        MRIinterpolateIntoVolume
        (gcam->mri_xind, (double)xf, (double)yf, (double)zf, (double)x) ;
        // src -> gcam volume position
        MRIinterpolateIntoVolume
        (gcam->mri_yind, (double)xf, (double)yf, (double)zf, (double)y) ;
        MRIinterpolateIntoVolume
        (gcam->mri_zind, (double)xf, (double)yf, (double)zf, (double)z) ;

        // mark counts (how many went in)
        MRIinterpolateIntoVolume
        (mri_counts, (double)xf, (double)yf, (double)zf, (double)1.0) ;
#ifndef __OPTIMIZE__
        if (xv == Ggca_x && yv == Ggca_y && zv == Ggca_z)
        {
          if (xv >=0 && yv >= 0 && zv >= 0)
            fprintf(stderr,
                    "src (%d, %d, %d) corresponds to gcam (%d, %d, %d)\n",
                    xv, yv, zv, x, y, z);
        }
#endif
      }
    }
  }

  if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
  {
    MRIwrite(gcam->mri_xind, "xi.mgz") ;
    MRIwrite(mri_counts, "c.mgz") ;
  }

  // xind, yind, zind is of size (width, height, depth)
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        // get count
        num = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
        if (num == 0)
        {
          continue ;  /* nothing there */
        }
        // give average gcam position for this points
        MRIFvox(gcam->mri_xind, x, y, z) =
          MRIFvox(gcam->mri_xind, x, y, z)/(float)num ;
        MRIFvox(gcam->mri_yind, x, y, z) =
          MRIFvox(gcam->mri_yind, x, y, z)/(float)num ;
        MRIFvox(gcam->mri_zind, x, y, z) =
          MRIFvox(gcam->mri_zind, x, y, z)/(float)num ;
        MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED ;
      }
    }
  }

  MRIfree(&mri_counts) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("performing soap bubble of x indices...\n") ;
  }
  MRIbuildVoronoiDiagram(gcam->mri_xind, mri_ctrl, gcam->mri_xind) ;
  if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
  {
    MRIwrite(gcam->mri_xind, "xi.mgz") ;
  }
  MRIsoapBubble(gcam->mri_xind, mri_ctrl, gcam->mri_xind, 50, 1) ;
  if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
  {
    MRIwrite(gcam->mri_xind, "xis.mgz") ;
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("performing soap bubble of y indices...\n") ;
  }
  MRIbuildVoronoiDiagram(gcam->mri_yind, mri_ctrl, gcam->mri_yind) ;
  MRIsoapBubble(gcam->mri_yind, mri_ctrl, gcam->mri_yind, 50, 1) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("performing soap bubble of z indices...\n") ;
  }
  MRIbuildVoronoiDiagram(gcam->mri_zind, mri_ctrl, gcam->mri_zind) ;
  MRIsoapBubble(gcam->mri_zind, mri_ctrl, gcam->mri_zind, 50, 1) ;
  MRIfree(&mri_ctrl) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(gcam->mri_xind, "xi.mgz") ;
    MRIwrite(gcam->mri_yind, "yi.mgz") ;
    MRIwrite(gcam->mri_zind, "zi.mgz") ;
  }
#ifndef __OPTIMIZE__
  xv = Ggca_x;
  yv = Ggca_y;
  zv = Ggca_z;
  if (xv >= 0 && yv >= 0 && zv >= 0)
    fprintf(stderr,
            "src (%d, %d, %d) corresponds to gcam (%2.1f, %2.1f, %2.1f)\n",
            xv, yv, zv,
            MRIFvox(gcam->mri_xind, xv, yv, zv),
            MRIFvox(gcam->mri_yind, xv, yv, zv),
            MRIFvox(gcam->mri_zind, xv, yv, zv));
#endif
  return(NO_ERROR) ;
}

int
GCAMfreeInverse(GCA_MORPH *gcam)
{
  if (gcam->mri_xind)
  {
    MRIfree(&gcam->mri_xind) ;
  }
  if (gcam->mri_yind)
  {
    MRIfree(&gcam->mri_yind) ;
  }
  if (gcam->mri_zind)
  {
    MRIfree(&gcam->mri_zind) ;
  }
  return(NO_ERROR) ;
}


int
gcamUndoGradient(GCA_MORPH *gcam)
{
  int            x, y, z ;
  float          dx, dy, dz ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        dx = gcamn->odx ;
        dy  = gcamn->ody ;
        dz = gcamn->odz ;
        gcamn->odx = gcamn->ody = gcamn->odz = 0 ;   /* turn off momentum */

        if (x == Gx && y == Gy && z == Gz)
          printf("UNGRAD: node(%d,%d,%d): "
                 "moving by (%2.3f, %2.3f, %2.3f) from "
                 "(%2.1f,%2.1f,%2.1f) to ",
                 x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
        gcamn->x -= dx;
        gcamn->y -= dy;
        gcamn->z -= dz;
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("(%2.1f,%2.1f,%2.1f)\n", gcamn->x, gcamn->y, gcamn->z) ;
        }
      }
  return(NO_ERROR) ;
}

int
gcamReadMRI(GCA_MORPH *gcam, MRI *mri, int which)
{
  int            x, y, z, i ;
  float          d ;
  GCA_MORPH_NODE *gcamn ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri, "m.mgz") ;
  }
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        d = MRIFvox(mri, x, y, z)  ;
        switch (which)
        {
        default:
          ErrorExit(ERROR_BADPARM,
                    "gcamReadMRI(%d): unknown which parameter", which) ;
          break ;
        case GCAM_X:
          gcamn->x = d ;
          break ;
        case GCAM_Y:
          gcamn->y = d;
          break ;
        case GCAM_Z:
          gcamn->z = d;
          break ;
        case GCAM_X_GRAD:
          gcamn->dx = d ;
          break ;
        case GCAM_Y_GRAD:
          gcamn->dy = d;
          break ;
        case GCAM_Z_GRAD:
          gcamn->dz = d;
          break ;
        case GCAM_LABEL:
          gcamn->label = nint(d) ;
          break ;
        case GCAM_ORIGX:
          gcamn->origx = d ;
          break ;
        case GCAM_ORIGY:
          gcamn->origy = d ;
          break ;
        case GCAM_ORIGZ:
          gcamn->origz = d ;
          break ;
        case GCAM_NODEX:
          gcamn->xn = d ;
          break ;
        case GCAM_NODEY:
          gcamn->yn = d ;
          break ;
        case GCAM_NODEZ:
          gcamn->zn = d ;
          break ;
        case GCAM_MEANS:
          if (gcamn->gc != NULL)
          {
            for (i = 0 ; i < mri->nframes ; i++)
            {
              d = MRIFseq_vox(mri, x, y, z,i)  ;
              gcamn->gc->means[i] = d ;
            }
          }
          break ;
        case GCAM_COVARS:
          if (gcamn->gc != NULL)
          {
            for (i = 0 ; i < mri->nframes ; i++)
            {
              d = MRIFseq_vox(mri, x, y, z,i)  ;
              gcamn->gc->covars[i] = d ;
              if (FZERO(d))
              {
                DiagBreak() ;
              }
            }
          }
          break ;
        }
      }
  return(NO_ERROR) ;
}

MRI *
GCAMwriteMRI(GCA_MORPH *gcam, MRI *mri, int which)
{
  int            x, y, z ;
  float          d ;
  GCA_MORPH_NODE *gcamn ;

  if (!mri)
  {
    mri = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
    MRIcopyVolGeomToMRI(mri, &gcam->atlas) ;
    MRIsetResolution(mri, gcam->atlas.xsize*gcam->spacing, gcam->atlas.ysize*gcam->spacing, gcam->atlas.zsize*gcam->spacing) ;
    MRIreInitCache(mri) ;
  }
  if (!mri)
    ErrorExit(ERROR_NOMEMORY, "gcamWrite: could not allocate %dx%dx%d MRI\n",
              gcam->width, gcam->height, gcam->depth) ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        switch (which)
        {
        default:
          ErrorExit
          (ERROR_BADPARM,
           "GCAMwriteMRI(%d): unknown which parameter", which) ;
          d = 0 ;
          break ;
	case GCAM_NEG:
	  if (gcamn->invalid == 0 && gcamn->orig_area1 > 0 && gcamn->orig_area2 > 0)
	    d = (gcamn->area1 <= 0 || gcamn->area2 <= 0) ;
	  else
	    d = 0 ;
	  if (d  > 0)
	    DiagBreak() ;
	  break ;
	case GCAM_MIN_AREA:
	  if (gcamn->invalid == 0)
	    d = MIN(gcamn->area1, gcamn->area2) ;
	  else
	    d = 0 ;
	  break ;
	case GCAM_LOG_MIN_AREA:
	  if (gcamn->invalid == 0)
	  {
	    d = MIN(gcamn->area1, gcamn->area2) ;
	    if (d > 1)
	      d = 0 ;
	    else if (d > 0)
	      d = -log10(fabs(d)) ;
	    else
	      d = log10(fabs(d)) ;
	  }
	  else
	    d = 0 ;
	  break ;
        case GCAM_X_GRAD:
          d = gcamn->dx ;
          break ;
        case GCAM_Y_GRAD:
          d = gcamn->dy ;
          break ;
        case GCAM_Z_GRAD:
          d = gcamn->dz ;
          break ;
        case GCAM_JACOBIAN:
          d = FZERO(gcamn->orig_area) ? 1.0 : gcamn->area / gcamn->orig_area ;
          break ;
        case GCAM_LOG_JACOBIAN:
          d = FZERO(gcamn->orig_area) ? 1.0 : gcamn->area / gcamn->orig_area ;
	  if (d < 0)
	    d = -100 ;
	  else
	    d = log(d) ;
          break ;
        case GCAM_AREA:
          d = gcamn->area ;
          break ;
        case GCAM_ORIG_AREA:
          d = gcamn->orig_area ;
          break ;
        case GCAM_X:
          d = gcamn->x ;
          break ;
        case GCAM_Y:
          d = gcamn->y ;
          break ;
        case GCAM_Z:
          d = gcamn->z ;
          break ;
        case GCAM_ORIGX:
          d = gcamn->origx ;
          break ;
        case GCAM_ORIGY:
          d = gcamn->origy ;
          break ;
        case GCAM_ORIGZ:
          d = gcamn->origz ;
          break ;
        case GCAM_LABEL:
          d = gcamn->label ;
          break ;
        case GCAM_MEANS:
          d = gcamn->gc ? gcamn->gc->means[0] : 0;
	  break ;
        }
        MRIFvox(mri, x, y, z) = d ;
      }
  return(mri) ;
}

int
gcamSmoothGradient(GCA_MORPH *gcam, int navgs)
{
#ifdef GCAM_SMOOTH_GRADIENT_GPU
  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamSmoothGradientGPU( gcam, navgs );
#else

#if 1
  MRI   *mri_tmp = NULL, *mri_kernel ;
  int   i ;

  if (navgs <= 0)     /* no smoothing */
  {
    return(NO_ERROR) ;
  }
#if 0
  mri_kernel = MRIgaussian1d(sqrt((float)navgs)*M_PI/2.0, 0) ;
#else
  mri_kernel = MRIgaussian1d(sqrt((float)navgs*2/M_PI), 0) ;
#endif
  if ((Gx >= 0 && Gy >= 0 && Gz >= 0) && (Gdiag & DIAG_SHOW))
    printf("before smoothing %d times: grad = (%2.3f, %2.3f, %2.3f)\n",
           navgs,
           gcam->nodes[Gx][Gy][Gz].dx,
           gcam->nodes[Gx][Gy][Gz].dy,
           gcam->nodes[Gx][Gy][Gz].dz) ;

  for (i = GCAM_X_GRAD ; i <= GCAM_Z_GRAD ; i++)
  {
    char fname[STRLEN] ;
    mri_tmp = GCAMwriteMRI(gcam, mri_tmp, i) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "grad%d_before.mgz", i) ;
      MRIwrite(mri_tmp, fname) ;
    }
    MRIconvolveGaussian(mri_tmp, mri_tmp, mri_kernel) ;
    gcamReadMRI(gcam, mri_tmp, i) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "grad%d_after.mgz", i) ;
      MRIwrite(mri_tmp, fname) ;
    }
  }

  MRIfree(&mri_tmp) ;
  MRIfree(&mri_kernel) ;

#else
  static double *dx, *dy, *dz = NULL ;
  int            x, y, z, xk, yk, zk, xi, yi, zi, index, i, num ;
  GCA_MORPH_NODE *gcamn, *gcamn_nbr ;

  if (navgs <= 0)
  {
    return(NO_ERROR) ;
  }

  if (!dz)
  {
    dx = (double *)calloc
         (gcam->width*gcam->height*gcam->depth, sizeof(double)) ;
    dy = (double *)calloc
         (gcam->width*gcam->height*gcam->depth, sizeof(double)) ;
    dz = (double *)calloc
         (gcam->width*gcam->height*gcam->depth, sizeof(double)) ;
    if (!dx || !dy || !dz)
      ErrorExit(ERROR_NOMEMORY,
                "gcamSmoothGradient: could not allocate grad buffers") ;
  }

  if ((Gx >= 0 && Gy >= 0 && Gz >= 0) && (Gdiag & DIAG_SHOW))
    printf("before smoothing %d times: grad = (%2.3f, %2.3f, %2.3f)\n",
           navgs,
           gcam->nodes[Gx][Gy][Gz].dx,
           gcam->nodes[Gx][Gy][Gz].dy,
           gcam->nodes[Gx][Gy][Gz].dz) ;

  for (i = 0 ; i < navgs ; i++)
  {
    for (index = x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
        for (z = 0 ; z < gcam->depth ; z++, index++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          gcamn = &gcam->nodes[x][y][z] ;

          if (gcamn->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }

          dx[index] = dy[index] = dz[index] = 0 ;
#if 0
          if (x == Gx && y == Gy && z == Gz)
            printf("GRAD: node(%d,%d,%d): moving by "
                   "(%2.3f, %2.3f, %2.3f) from "
                   "(%2.1f,%2.1f,%2.1f) to ",
                   x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
          if (x == Gx && y == Gy && z == Gz)
          {
            printf("(%2.1f,%2.1f,%2.1f)\n", gcamn->x, gcamn->y, gcamn->z) ;
          }
#endif
          for (num = 0, xk = -1 ; xk <= 1 ; xk++)
          {
            xi = x+xk ;
            if (xi < 0)
            {
              xi = 0 ;
            }
            if (xi >= gcam->width)
            {
              xi = gcam->width-1 ;
            }
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = y+yk ;
              if (yi < 0)
              {
                yi = 0 ;
              }
              if (yi >= gcam->height)
              {
                yi = gcam->height-1 ;
              }
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = z+zk ;
                if (zi < 0)
                {
                  zi = 0 ;
                }
                if (zi >= gcam->depth)
                {
                  zi = gcam->depth-1 ;
                }
                gcamn_nbr = &gcam->nodes[xi][yi][zi] ;

                if (gcamn_nbr->invalid/* == GCAM_POSITION_INVALID*/)
                {
                  continue;
                }
#if 1
                if (gcamn_nbr->label != gcamn->label)
                {
                  continue ;
                }
#endif
                num++ ;
                dx[index] += gcamn_nbr->dx ;
                dy[index] += gcamn_nbr->dy ;
                dz[index] += gcamn_nbr->dz ;
              }
            }
          }
          dx[index] /= (double)num ;
          dy[index] /= (double)num ;
          dz[index] /= (double)num ;
        }
    for (index = x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
        for (z = 0 ; z < gcam->depth ; z++, index++)
        {
          gcamn = &gcam->nodes[x][y][z] ;

          if (gcamn->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }

          gcamn->dx = dx[index] ;
          gcamn->dy = dy[index] ;
          gcamn->dz = dz[index] ;

        }
  }

#endif

  if ((Gx >= 0 && Gy >= 0 && Gz >= 0) && (Gdiag & DIAG_SHOW))
    printf("after smoothing %d times: grad = (%2.3f, %2.3f, %2.3f)\n",
           navgs,
           gcam->nodes[Gx][Gy][Gz].dx,
           gcam->nodes[Gx][Gy][Gz].dy,
           gcam->nodes[Gx][Gy][Gz].dz) ;

#endif

  return(NO_ERROR) ;
}


int
GCAMsetStatus(GCA_MORPH *gcam, int status)
{

  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        gcamn->status = status ;
      }
  return(NO_ERROR) ;
}


int
GCAMsetLabelStatus(GCA_MORPH *gcam, int label, int status)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->label == label)
        {
          gcamn->status = status ;
        }
      }
  return(NO_ERROR) ;
}
#define MAX_SAMPLES 100

#define GCAM_FOTS_OUTPUT 0
#define GCAM_FOTS_TIMERS 0

double
gcamFindOptimalTimeStep(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms, MRI *mri)
{
#if GCAM_FOTS_TIMERS
  Chronometer tFOTS;

  InitChronometer( &tFOTS );
  StartChronometer( &tFOTS );
#endif

  double min_dt;

#ifndef GCAM_FIND_OPTIMAL_TIMESTEP_GPU
  MATRIX   *mX, *m_xTx, *m_xTx_inv, *m_xTy, *mP, *m_xT ;
  double   rms, min_rms, dt_in[MAX_SAMPLES],
           rms_out[MAX_SAMPLES],orig_dt,
           a, b, c, max_dt, start_dt, start_rms, pct_change ;
  VECTOR   *vY ;
  int      N, i, Gxs, Gys, Gzs,bad, suppressed=0, prev_neg ;
  long     diag ;
#endif

#if GCAM_FOTS_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;

  if( (nCalls%outputFreq) == 0 )
  {
    char fname[STRLEN];

    snprintf( fname, STRLEN-1, "gcam%05u", nCalls/outputFreq );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mri%05u.mgz", nCalls/outputFreq );
    MRIwrite( mri, fname );
  }
#endif

#ifdef GCAM_FIND_OPTIMAL_TIMESTEP_GPU
#if SHOW_EXEC_LOC
  printf( "%s: CUDA call\n", __FUNCTION__ );
#endif
  min_dt = gcamFindOptimalTimestepGPU( gcam, parms, mri );
#else
#if SHOW_EXEC_LOC
  printf( "%s: CPU call\n", __FUNCTION__ );
#endif
  gcamClearMomentum(gcam) ;
  // disable diagnostics so we don't see every time step sampled
  Gxs = Gx ;
  Gys = Gy ;
  Gzs = Gz ;
  diag = Gdiag ;
  Gx = Gy = Gz = -1 ;
  Gdiag = 0 ;

  start_rms = min_rms = GCAMcomputeRMS(gcam, mri, parms) ;
  min_dt = 0 ;

  /* first find right order of magnitude for time step */
  orig_dt = parms->dt ;
  if (DZERO(orig_dt))
  {
    orig_dt = 1e-6 ;
    ErrorPrintf(ERROR_BADPARM, "orig dt = 0 in gcamFindOptimalTimeStep");
  }

  /* define a pretty broad search range initially */
  start_dt = (sqrt(parms->navgs)+1)*orig_dt / (10*16.0*16.0) ;
  max_dt =   (sqrt(parms->navgs)+1)*orig_dt * (10*16.0*16.0) ;

  i = 0 ;
  do
  {
    if (i++ > 0)
    {
      start_dt *= 0.01 ;
    }
    if (DEQUAL(start_dt, 0))
    {
      break ;
    }
    bad = 0 ;
    prev_neg = 0 ;
    for (parms->dt = start_dt ; parms->dt <= max_dt ; parms->dt*=4)
    {
      gcamApplyGradient(gcam, parms) ;
      rms = GCAMcomputeRMS(gcam, mri, parms) ;
      gcamUndoGradient(gcam) ;
      if ((gcam->neg > 0) && 0)
      {
        pct_change = 100.0*(start_rms-min_rms)/start_rms ;
        //    if (pct_change < parms->tol)
        {
          if (getenv("SHOW_NEG"))
          {
            printf("!!!!!!!!!!reducing gradient magnitude around %d negative nodes\n",gcam->neg);
          }
          gcamSuppressNegativeGradients(gcam, 0.5) ;
          if (suppressed++ < 1000 &&
              ((gcam->neg < prev_neg) || (prev_neg <= 0)))
          {
            prev_neg = gcam->neg ;
            parms->dt /= 4 ;  // try again at this scale
            continue ;
          }
        }
      }
      else
      {
        prev_neg = suppressed = 0 ;
      }

      if (gcam->neg > 0 && gcam_write_grad)  // diagnostic
      {
        int x, y, z, k ;
        GCA_MORPH_NODE *gcamn ;
        for (k = x = 0 ; x < gcam->width ; x++)
        {
          for (y = 0 ; y < gcam->height ; y++)
          {
            for (z = 0 ; z < gcam->depth ; z++)
            {
              gcamn = &gcam->nodes[x][y][z] ;
              if (gcamn->invalid == GCAM_VALID && gcamn->area < 0)
              {
                if (k++ <  10)
                {
                  printf("gcam(%d, %d, %d) neg, areas = %2.2f, %2.2f\n", x, y, z,gcamn->area1,gcamn->area2) ;
                }
              }
            }
          }
        }
      }

      if (gcam->neg > 0 && getenv("SHOW_NEG"))
      {
        printf("%d negative nodes at dt=%2.6f\n", gcam->neg, parms->dt) ;
      }

      if (gcam->neg > 0 && DEQUAL(parms->dt, start_dt))
      {
        start_dt *= 0.1 ;
        parms->dt /= 4 ;
        continue ;
      }
#if 0
      if (gcam->neg <100)
      {
        GCA_MORPH_NODE *gcamn ;
        int x,y,z;
        if (bad++ > 10)
        {
          break ;
        }
        for (x=0; x<gcam->width; x++)
        {
          for (y=0; y<gcam->height; y++)
          {
            for (z=0; z<gcam->depth; z++)
            {
              gcamn = &gcam->nodes[x][y][z] ;
              if (gcamn->area < 0)
              {
                gcamSetGradientToNbrAverage(gcam, x, y, z) ;
              }
            }
          }
        }
        parms->dt /= 4 ;  /* do it again */
        continue ;
      }
      if (gcam->neg == 0)
      {
        bad = 0 ;
      }
#endif

      if (gcam->neg && parms->noneg == True)
      {
        break ;
      }
      if (rms < min_rms)
      {
        min_rms = rms ;
        min_dt = parms->dt ;
      }
      if (DEQUAL(min_dt, max_dt))
      {
        max_dt *= 10 ;
      }
    }
    if (i > 10)
    {
      break ;
    }
    max_dt = start_dt*10 ; /* will only iterate if min was at start_dt */
  }
  while (DEQUAL(min_dt,start_dt)) ;

  dt_in[0] = min_dt ;
  rms_out[0] = min_rms ;
  parms->dt = dt_in[1] = dt_in[0] - dt_in[0]*.2 ;
  gcamApplyGradient(gcam, parms) ;
  rms = rms_out[1] = GCAMcomputeRMS(gcam, mri, parms) ;
  gcamUndoGradient(gcam) ;
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg <= 0))
  {
    min_rms = rms ;
    min_dt = parms->dt ;
  }

  parms->dt = dt_in[2] = dt_in[0] - dt_in[0]*.4 ;
  gcamApplyGradient(gcam, parms) ;
  rms = rms_out[2] = GCAMcomputeRMS(gcam, mri, parms) ;
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg <= 0))
  {
    min_rms = rms ;
    min_dt = parms->dt ;
  }
  gcamUndoGradient(gcam) ;

  parms->dt = dt_in[3] = dt_in[0] + dt_in[0]*.2 ;
  gcamApplyGradient(gcam, parms) ;
  rms = rms_out[3] = GCAMcomputeRMS(gcam, mri, parms) ;
  gcamUndoGradient(gcam) ;
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg <= 0))
  {
    min_rms = rms ;
    min_dt = parms->dt ;
  }


  parms->dt = dt_in[4] = dt_in[0] + dt_in[0]*.4 ;
  gcamApplyGradient(gcam, parms) ;
  rms = rms_out[4] = GCAMcomputeRMS(gcam, mri, parms) ;
  gcamUndoGradient(gcam) ;
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg <= 0))
  {
    min_rms = rms ;
    min_dt = parms->dt ;
  }


  /* now compute location of minimum of best quadratic fit */
  N = 5 ;  /* min_dt +- .1*min_dt and .2*min_dt */
  mX = MatrixAlloc(N, 3, MATRIX_REAL) ;
  vY = VectorAlloc(N, MATRIX_REAL) ;

  for (i = 1 ; i <= N ; i++)
  {
    *MATRIX_RELT(mX, i, 1) = dt_in[i-1] * dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 2) = 2*dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 3) = 1.0f ;

    VECTOR_ELT(vY, i) = rms_out[i-1] ;
  }

  m_xT = MatrixTranspose(mX, NULL) ;
  m_xTx = MatrixMultiply(m_xT, mX, NULL) ;
  m_xTx_inv = MatrixInverse(m_xTx, NULL) ;
  if (m_xTx_inv)
  {
    m_xTy = MatrixMultiply(m_xT, vY, NULL) ;
    mP = MatrixMultiply(m_xTx_inv, m_xTy, NULL) ;
    a = RVECTOR_ELT(mP, 1) ;
    b = RVECTOR_ELT(mP, 2) ;
    c = RVECTOR_ELT(mP, 3);
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,
              "(a,b,c) = (%2.3f, %2.3f, %2.3f), predicted min at %2.3f\n",
              a, b, c, -b/a) ;
    if (!finitep(a))
    {
      DiagBreak() ;
    }
    MatrixFree(&mP) ;
    MatrixFree(&m_xTx_inv) ;
    MatrixFree(&m_xTy) ;
    if (finitep(a) && !FZERO(a))
    {
      parms->dt = -b/a ;
      gcamApplyGradient(gcam, parms) ;
      rms = GCAMcomputeRMS(gcam, mri, parms) ;
      gcamUndoGradient(gcam) ;
      if (rms < min_rms && (gcam->neg == 0 || parms->noneg <= 0))
      {
        min_rms = rms ;
        min_dt = parms->dt ;
      }
    }
  }
  MatrixFree(&m_xT) ;
  MatrixFree(&m_xTx) ;
  MatrixFree(&mX) ;
  VectorFree(&vY) ;

  gcamComputeMetricProperties(gcam) ;
  parms->dt = orig_dt ;
  Gx = Gxs ;
  Gy = Gys ;
  Gz = Gzs ;
  Gdiag = diag ;

#endif

#if GCAM_FOTS_OUTPUT
  nCalls++;
#endif

#if GCAM_FOTS_TIMERS
  StopChronometer( &tFOTS );
  printf( "%s: Complete in %9.3f ms\n",
          __FUNCTION__, GetChronometerValue( &tFOTS ) );
#endif

  return(min_dt) ;
}

#define GCAM_LABELENERGY_OUTPUT 0

double
gcamLabelEnergy( const GCA_MORPH *gcam,
                 const MRI *mri,
                 const double label_dist)
{

#if GCAM_LABELENERGY_OUTPUT
  const unsigned int gcamLabelEnergyOutputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%gcamLabelEnergyOutputFreq) == 0 )
  {
    char fname[STRLEN];
    snprintf( fname, STRLEN-1,
              "gcamLabelEnergyinput%04u", nCalls / gcamLabelEnergyOutputFreq );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );
  }
  nCalls++;
#endif

  float           sse ;

#ifdef GCAM_LABEL_ENERGY_GPU
#if SHOW_EXEC_LOC
  printf( "%s: CUDA call\n", __FUNCTION__ );
#endif
  sse = gcamLabelEnergyGPU( gcam );
#else
#if SHOW_EXEC_LOC
  printf( "%s: CPU call\n", __FUNCTION__ );
#endif
  int             x, y, z ;
  GCA_MORPH_NODE  *gcamn ;

  sse = 0.0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }


        gcamn = &gcam->nodes[x][y][z] ;

        if ((gcamn->invalid == GCAM_POSITION_INVALID) ||
            ((gcamn->status & GCAM_LABEL_NODE) == 0))
        {
          continue;
        }

        sse += fabs(gcamn->label_dist) ;
      }
    }
  }

#endif


  return(sse) ;
}

// ------------------------------------------------------------


#include "voxlist.h"




#define GCAM_REMOVE_LABEL_OUTLIERS_OUTPUT 0

int
remove_label_outliers( GCA_MORPH *gcam,
                       MRI *mri_dist,
                       const int whalf,
                       const double thresh )
{

  int nremoved;

#ifdef GCAM_LABEL_TERM_REMOVE_OUTLIERS_GPU

  printf( "%s: On GPU\n", __FUNCTION__ );
  nremoved = gcamRemoveLabelOutliersGPU( gcam, mri_dist, whalf, thresh );

#else
  int         nremoved_total, n, i, vox_to_examine, niters ;
  MRI         *mri_std, *mri_ctrl, *mri_tmp, *mri_ctrl_tmp ;
  VOXEL_LIST  *vl ;
  float       diff, val0, oval ;
  int         delete, xv, yv, zv, xo, yo, zo, x, y, z ;
  GCA_MORPH_NODE  *gcamn, *gcamn_sup, *gcamn_inf ;
  double      max_change ;


#if GCAM_REMOVE_LABEL_OUTLIERS_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%outputFreq) == 0 )
  {
    unsigned int nOut = nCalls/outputFreq;

    char fname[STRLEN];

    snprintf( fname, STRLEN-1, "gcamRemoveLabelOutliersInput%04u", nOut );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mriRemoveLabelOutliersInput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri_dist, fname );
  }
  nCalls++;

#endif

  mri_ctrl = MRIcloneDifferentType(mri_dist, MRI_UCHAR) ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->status & GCAM_LABEL_NODE)
        {
          MRIsetVoxVal(mri_ctrl, x, y, z, 0, CONTROL_MARKED);
        }
      }
    }
  }

  niters = 100 ;
  for (nremoved_total = i = 0 ; i < niters ; i++)
  {
    mri_std = MRIstdNonzero(mri_dist, NULL, mri_dist, whalf*2+1) ;
    vl = VLSTcreate(mri_std, .001, 1e10, NULL, 0, 0) ;
    VLSTsort(vl, vl) ;

    vox_to_examine = ceil(vl->nvox*.05) ;
    for (nremoved = n = 0 ; n < vox_to_examine ; n++)
    {
      x = vl->xi[n] ;
      y = vl->yi[n] ;
      z = vl->zi[n] ;
      gcamn = &gcam->nodes[x][y][z] ;

      if (x == Gx && y == Gy && z == Gz)
      {
        DiagBreak() ;
      }

      if ((gcamn->status & GCAM_LABEL_NODE) == 0 || (gcamn->status & GCAM_MANUAL_LABEL))
      {
        continue ;
      }

      val0 = MRIgetVoxVal(mri_dist, x, y, z, 0) ;
      delete = 0 ;
      for (xo = -whalf ; xo <= whalf && !delete ; xo++)
      {
        xv = mri_dist->xi[x+xo] ;
        for (yo = -whalf ; yo <= whalf && !delete ; yo++)
        {
          yv = mri_dist->yi[y+yo] ;
          for (zo = -whalf ; zo <= whalf && !delete ; zo++)
          {
            zv = mri_dist->zi[z+zo] ;
            oval = MRIgetVoxVal(mri_dist, xv, yv, zv, 0) ;
            if (!FZERO(oval))
            {
              diff = fabs(oval - val0) ;
              if (diff > thresh)   /* shouldn't ever be this big */
              {
                if (fabs(val0) > fabs(oval) && (val0*oval < 0))
                  /*if their signs are different and difference is big */
                {
                  delete = 1 ;
                }
              }
            }
          }
        }
      }
      if (delete)
      {
        nremoved++ ;
        MRIFvox(mri_dist, x, y, z) = 0 ;
        MRIsetVoxVal(mri_ctrl, x, y, z, 0, CONTROL_NONE) ;
        gcamn->dy = 0 ;
        //        gcamn->status = GCAM_USE_LIKELIHOOD ;
        gcamn_inf = &gcam->nodes[x][y+1][z] ;
        if ((gcamn_inf->status & GCAM_LABEL_NODE) && ((gcamn_inf->status & GCAM_MANUAL_LABEL) == 0))
        {
          nremoved++ ;
          //          gcamn_inf->status = GCAM_USE_LIKELIHOOD ;
          MRIsetVoxVal(mri_dist, x, y+1, z, 0, 0) ;
          MRIsetVoxVal(mri_ctrl, x, y+1, z, 0, CONTROL_NONE) ;
          gcamn_inf->dy = 0 ;
        }
        gcamn_sup = &gcam->nodes[x][y-1][z] ;
        if ((gcamn_sup->status & GCAM_LABEL_NODE) && ((gcamn_sup->status & GCAM_MANUAL_LABEL) == 0))
        {
          nremoved++ ;
          //          gcamn_sup->status = GCAM_USE_LIKELIHOOD ;
          gcamn_sup->dy = 0 ;
          MRIsetVoxVal(mri_dist, x, y-1, z, 0, 0) ;
          MRIsetVoxVal(mri_ctrl, x, y-1, z, 0, CONTROL_NONE) ;
        }
      }
    }

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      char fname[STRLEN] ;
      sprintf(fname, "dist_after%d.mgz",i) ;
      MRIwrite(mri_dist, fname) ;
    }
    nremoved_total += nremoved ;
    MRIfree(&mri_std) ;
    if (nremoved == 0)   // nothing happened
    {
      break ;
    }

    VLSTfree(&vl) ; // keep it for last iteration
  }




  // now use soap bubble smoothing to estimate label offset of deleted locations
  mri_tmp = NULL ;
  mri_ctrl_tmp = NULL;
  for (i = 0 ; i < 100 ; i++)
  {
    max_change = 0.0 ;
    mri_tmp = MRIcopy(mri_dist, mri_tmp) ;
    mri_ctrl_tmp = MRIcopy( mri_ctrl, mri_ctrl_tmp );

    for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
      {
        for (z = 0 ; z < gcam->depth ; z++)
        {
          int    xi, yi, zi, xk, yk, zk, num ;
          double mean ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          gcamn = &gcam->nodes[x][y][z] ;

          if( MRIgetVoxVal(mri_ctrl, x, y, z, 0) == CONTROL_MARKED ||
              (gcamn->status & GCAM_LABEL_NODE) == 0 )
          {
            continue ;
          }

          for (xk = -1, num = 0, mean = 0.0 ; xk <= 1 ; xk++)
          {
            xi = mri_ctrl->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_ctrl->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = mri_ctrl->zi[z+zk] ;
                if (MRIgetVoxVal(mri_ctrl, xi, yi, zi, 0) == CONTROL_MARKED ||
                    MRIgetVoxVal(mri_ctrl, xi, yi, zi, 0) == CONTROL_NBR)
                {
                  mean += MRIgetVoxVal(mri_dist, xi, yi, zi, 0) ;
                  num++ ;
                }
              }
            }
          }
          if (num > 0)
          {
            float val ;
            val = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;
            mean /= num ;
            if (fabs(val-mean) > max_change)
            {
              max_change = fabs(val-mean) ;
            }
            MRIsetVoxVal(mri_tmp, x, y, z, 0, mean) ;
            MRIsetVoxVal(mri_ctrl_tmp, x, y, z, 0, CONTROL_TMP) ;
            gcamn->status = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ;
          }
        }
      }
    }

    MRIcopy(mri_tmp, mri_dist) ;
    MRIcopy( mri_ctrl_tmp, mri_ctrl );

#if 0
    {
      char fname[STRLEN] ;
      sprintf(fname, "dist_after%d.mgz",i) ;
      MRIwrite(mri_dist, fname) ;
    }
#endif

    MRIreplaceValuesOnly(mri_ctrl, mri_ctrl, CONTROL_TMP, CONTROL_NBR) ;
    if (max_change < 0.05)
    {
      break ;
    }
  }

  MRIfree(&mri_ctrl);
  MRIfree(&mri_tmp) ;
  MRIfree(&mri_ctrl_tmp);

#if GCAM_REMOVE_LABEL_OUTLIERS_OUTPUT
  if( nremoved != 0 )
  {
    printf( "%s: nremoved = %i on nCall %i\n",
            __FUNCTION__, nremoved, nCalls-1 );
  }
#endif

#endif

  return(nremoved) ;
}



// ====================================================
// Separate out some operations from gcamLabelTerm


#define MAX_MLE_DIST 1

#define GCAM_LABEL_COPYDELTAS_OUTPUT 0

void gcamLabelTermCopyDeltas( GCA_MORPH *gcam,
                              const MRI* mri_dist,
                              const double l_label )
{

#ifdef GCAM_LABEL_TERM_COPYDELTAS_GPU

  printf( "%s: On GPU\n", __FUNCTION__ );

  gcamLabelTermCopyDeltasGPU( gcam, mri_dist, l_label );
#else
  int x, y, z;
  GCA_MORPH_NODE *gcamn;

#if GCAM_LABEL_COPYDELTAS_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%outputFreq) == 0 )
  {
    unsigned int nOut = nCalls/outputFreq;

    char fname[STRLEN];

    snprintf( fname, STRLEN-1, "gcamLabelCopyDeltasInput%04u", nOut );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mriLabelCopyDeltasInput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri_dist, fname );
  }
  nCalls++;
#endif

  // copy deltas from mri_dist into gcam struct
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        double dy ;

        gcamn = &gcam->nodes[x][y][z] ;

        if ((gcamn->invalid/* == GCAM_POSITION_INVALID*/) ||
            ((gcamn->status & GCAM_LABEL_NODE) == 0))
        {
          continue;
        }

        dy = MRIgetVoxVal(mri_dist, x, y,z, 0) ;
        gcamn->label_dist = dy ;   /* for use in label energy */

        if (fabs(dy) > MAX_MLE_DIST)
        {
          dy = dy * MAX_MLE_DIST / fabs(dy) ;
        }

        gcamn->dy += l_label * dy ;
      }
    }
  }
#endif

}


// ----------------------

#define GCAM_LABEL_POSTANTCONSIST_OUTPUT 0

int gcamLabelTermPostAntConsistency( GCA_MORPH *gcam,
                                     MRI* mri_dist )
{
  /*!
    The main loop within this routine is very nasty.
    Each voxel depends on those on either side in x,
    and in z.
  */
  int nremoved;
#ifdef GCAM_LABEL_TERM_POSTANT_GPU

  printf( "%s: On GPU\n", __FUNCTION__ );

  nremoved = gcamLabelTermPostAntConsistencyGPU( gcam, mri_dist );
#else
  int x, y, z;

  GCA_MORPH_NODE *gcamn;
  GCA_MORPH_NODE *gcamn_medial, *gcamn_lateral;
  GCA_MORPH_NODE *gcamn_post, *gcamn_ant;

#if GCAM_LABEL_POSTANTCONSIST_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%outputFreq) == 0 )
  {
    unsigned int nOut = nCalls/outputFreq;

    char fname[STRLEN];

    snprintf( fname, STRLEN-1, "gcamLabelPostAntConsistInput%04u", nOut );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mriLabelPostAntConsistInput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri_dist, fname );
  }
  nCalls++;
#endif


  nremoved = 0;

  /* do posterior/anterior consistency check */
  for (x = gcam->width-1 ; x >=0 ; x--)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {

        gcamn = &gcam->nodes[x][y][z] ;

        if ((gcamn->invalid/* == GCAM_POSITION_INVALID*/) ||
            ((gcamn->status & GCAM_LABEL_NODE) == 0) ||
            ((gcam->status & GCAM_MANUAL_LABEL)))
        {
          continue;
        }

        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if ((fabs(gcamn->x-Gvx)<=gcam->spacing) &&
            (fabs(gcamn->y-Gvy)<=gcam->spacing) &&
            (fabs(gcamn->z-Gvz)<=gcam->spacing))
        {
          DiagBreak() ;
        }

        /*
           if sign of nodes to left and right (medial and lateral)
        are both different, then don't trust this one.
        */
        if ((x < gcam->width-1) && (x > 0))
        {
          gcamn_medial = &gcam->nodes[x-1][y][z] ;
          gcamn_lateral = &gcam->nodes[x+1][y][z] ;
          if (((gcamn_medial->status & GCAM_LABEL_NODE) == 0) ||
              ((gcamn_lateral->status & GCAM_LABEL_NODE) == 0))
            /* only if they are both label nodes */
          {
            continue ;
          }

          if ((gcamn_medial->dy*gcamn_lateral->dy > 0) &&
              (gcamn_medial->dy*gcamn->dy < 0))
          {
            gcamn->status = GCAM_USE_LIKELIHOOD ;
            MRIsetVoxVal(mri_dist, x, y, z, 0, 0) ;
            nremoved++ ;
            continue ;
          }
        }
        if ((z < gcam->depth-1) && (z > 0))
        {
          gcamn_post = &gcam->nodes[x][y][z-1] ;
          gcamn_ant = &gcam->nodes[x][y][z+1] ;
          if (((gcamn_ant->status & GCAM_LABEL_NODE) == 0) ||
              ((gcamn_post->status & GCAM_LABEL_NODE) == 0))
            /* only if they are both label nodes */
          {
            continue ;
          }

          if ((gcamn_ant->dy*gcamn_post->dy > 0) &&
              (gcamn_ant->dy*gcamn->dy < 0))
          {
            gcamn->status = GCAM_USE_LIKELIHOOD ;
            MRIsetVoxVal(mri_dist, x, y, z, 0, 0) ;
            nremoved++ ;
            continue ;
          }
        }

      }
    }
  }
#endif

  return( nremoved );
}


// ----------------------


#define GCAM_LABELFINAL_OUTPUT 0

int gcamLabelTermFinalUpdate( GCA_MORPH *gcam,
                              const MRI* mri_dist,
                              const double l_label )
{
#ifdef GCAM_LABEL_TERM_FINAL_GPU
  int num;

  printf( "%s: On GPU\n", __FUNCTION__ );

  num = gcamLabelTermFinalUpdateGPU( gcam, mri_dist, l_label );

  return( num );
#else
  int num, x, y, z;
  GCA_MORPH_NODE *gcamn;

#if GCAM_LABELFINAL_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%outputFreq) == 0 )
  {
    unsigned int nOut = nCalls/outputFreq;

    char fname[STRLEN];

    snprintf( fname, STRLEN-1, "gcamLabelFinalInput%04u", nOut );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mriLabelFinalInput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri_dist, fname );
  }
  nCalls++;
#endif

  num = 0;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {

        gcamn = &gcam->nodes[x][y][z] ;

        if ((gcamn->invalid/* == GCAM_POSITION_INVALID*/) ||
            ((gcamn->status & GCAM_LABEL_NODE) == 0))
        {
          continue;
        }

        gcamn->dy = l_label * MRIgetVoxVal(mri_dist, x, y,z, 0) ;

        if (fabs(gcamn->dy)/l_label >= 1)
        {
          num++ ;
        }

        gcamn->label_dist = gcamn->dy ;   /* for use in label energy */

      }
    }
  }

  return( num );
#endif
}

// ----------------------

#define GCAM_LABEL_MAINLOOP_OUTPUT 0

//! Sliced out of gcamLabelTerm
void gcamLabelTermMainLoop( GCA_MORPH *gcam, const MRI *mri,
                            MRI *mri_dist,
                            const double l_label, const double label_dist, MRI *mri_twm )
{
#ifdef GCAM_LABEL_TERM_MAINLOOP_GPU

  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamLabelTermMainLoopGPU( gcam, mri, mri_dist, l_label, label_dist );

#else
  int x, y, z;
  int xn, yn, zn, ctrl_point_found ;
  int best_label, sup_wm, sup_ven, wm_label;
  int n;
  double dy;
  GCA_MORPH_NODE *gcamn, *gcamn_inf, *gcamn_sup;
  GCA_NODE *gcan;
  GC1D *wm_gc;
  float vals[MAX_GCA_INPUTS], yi, yk, min_dist;

#if GCAM_LABEL_MAINLOOP_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%outputFreq) == 0 )
  {
    unsigned int nOut = nCalls/outputFreq;

    char fname[STRLEN];

    snprintf( fname, STRLEN-1, "gcamLabelMainLoopInput%04u", nOut );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mri_distLabelMainLoopInput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri_dist, fname );

    snprintf( fname, STRLEN-1, "mriLabelMainLoopInput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri, fname );

    snprintf( fname, STRLEN-1, "gcaLabelMainLoop%04u.gca", nOut );
    GCAwrite( gcam->gca, fname );
  }
  nCalls++;
#endif


  if (mri_twm != NULL)
  {
    MRI *mri_morphed =
      GCAMmorphToAtlas(mri_twm, gcam, NULL, -1, SAMPLE_NEAREST) ;
    MRIwrite(mri_morphed, "twm.morphed.mgz") ;
    MRIfree(&mri_morphed) ;
  }

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        ctrl_point_found = 0 ;
        gcamn->status &= ~GCAM_MANUAL_LABEL ;

        // Skip invalid nodes
        if( gcamn->invalid == GCAM_POSITION_INVALID )
        {
          continue;
        }

        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (x == Gx && y == (Gy-1) && z == Gz)
        {
          DiagBreak() ;
        }
        if (x == Gx && y == (Gy+1) && z == Gz)
        {
          DiagBreak() ;
        }

        // Can't do top or bottom layers
        if( (y == gcam->height-1) || (y == 0) )
        {
          continue;
        }

        if ((fabs(gcamn->x-Gvx)<=gcam->spacing) &&
            (fabs(gcamn->y-Gvy)<=gcam->spacing) &&
            (fabs(gcamn->z-Gvz)<=gcam->spacing))
        {
          DiagBreak() ;
        }

        /* only process nodes which are hippocampus superior to something else,
           or white matter inferior to hippocampus
        */
        if (y == 0 || y == gcam->height-1 ||
            gcamn->y == 0 || gcamn->y == mri->height-1)
        {
          continue;
        }

        if( !IS_HIPPO(gcamn->label) && !IS_WM(gcamn->label) )
        {
          continue;
        }

        if ( !IS_WM(gcamn->label) )     /* only do white matter for now */
        {
          continue;
        }

        if (fabs(2*x-107) <= 2 && fabs(2*y-162)<=2 && fabs(2*z-133)<=2)
        {
          DiagBreak() ;
        }

        gcamn_inf = &gcam->nodes[x][y+1][z] ;
        gcamn_sup = &gcam->nodes[x][y-1][z] ;
        if (
          ((IS_HIPPO(gcamn->label) && IS_WM(gcamn_inf->label)) ||
           (IS_WM(gcamn->label) && IS_HIPPO(gcamn_sup->label))) == 0)
        {
          continue ;  /* only hippo above wm, or wm below hippo */
        }

        if (nint(gcamn->x) == Gsx && nint(gcamn->z) == Gsz)
        {
          DiagBreak() ;
        }
        if (IS_HIPPO(gcamn->label))
        {
          load_vals(mri, gcamn->x, gcamn->y+1, gcamn->z, vals, gcam->ninputs);
        }
        else
        {
          load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->ninputs);
        }

#if 0
        label = gcamMLElabelAtLocation(gcam, x, y, z, vals) ;
        if (IS_WM(label))  /* already has wm immediately inferior */
        {
          continue ;
        }
#endif
        if (GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn) != NO_ERROR)
        {
          continue ;
        }

        if ((IS_HIPPO(gcamn->label) && gcamn->label == Left_Hippocampus) ||
            (IS_WM(gcamn->label) &&
             gcamn->label == Left_Cerebral_White_Matter))
        {
          wm_label = Left_Cerebral_White_Matter ;
        }
        else
        {
          wm_label = Right_Cerebral_White_Matter ;
        }

        wm_gc = GCAfindPriorGC(gcam->gca, x, y, z, wm_label) ;

        if (wm_gc == NULL)
        {
          continue;
        }

        gcan = GCAbuildRegionalGCAN(gcam->gca, xn, yn, zn, 3) ;

        // ventral DC is indistinguishible from temporal wm pretty much
        for (n = 0 ; n < gcan->nlabels; n++)
        {
          if ((gcan->labels[n] == Left_VentralDC ||
               gcan->labels[n] == Right_VentralDC ||
               gcan->labels[n] == Brain_Stem) &&
              gcan->gcs[n].means[0] > 90)
          {
            gcan->labels[n] = wm_label ;
          }
        }

        dy = 0 ;
        min_dist = label_dist+1 ;
        sup_ven = sup_wm = 0 ;  /* if can't find any wm superior, then
                                   must be partial volume and don't trust */
#define SAMPLE_DIST 0.1
        for (yk = -label_dist ; yk <= label_dist ; yk += SAMPLE_DIST)
        {
          yi = gcamn->y+yk ;   /* sample inferiorly */
          if ((yi >= (mri->height-1)) || (yi <= 0))
          {
            break ;
          }

          if (mri_twm && MRIgetVoxVal(mri_twm, nint(gcamn->x), nint(yi), nint(gcamn->z), 0) >0)
          {
            min_dist = yk ;
            ctrl_point_found = 1 ;
            break ;
          }
          load_vals(mri, gcamn->x, yi, gcamn->z, vals, gcam->ninputs) ;
          best_label = GCAmaxLikelihoodLabel(gcan, vals, gcam->ninputs, NULL) ;
          if (yk < 0)
          {
            if (IS_CSF(best_label))
            {
              sup_ven++ ;
            }
            else if (sup_ven < 3/SAMPLE_DIST)
            {
              sup_ven = 0 ;
            }
          }

          if (IS_CSF(best_label) &&
              sup_ven > 2/SAMPLE_DIST && yk < 0)
            // shouldn't have to go through CSF to get to wm superiorly
          {
            min_dist = label_dist+1 ;
          }

          if (best_label != gcamn->label)
          {
            continue ;
          }

          if (yk < 0 && IS_WM(best_label))
          {
            sup_wm = 1 ;
          }

          if (fabs(yk) < fabs(min_dist))
          {
            if (is_temporal_wm(gcam, mri, gcan, gcamn->x, yi, gcamn->z, gcam->ninputs))
            {
              min_dist = yk ;
            }
          }
        }

        /* if inferior to lateral ventricle (as opposed to
           temporal horn) and can't find any
           wm above then it means the wm is partial-volumed
           and don't trust estimate */
        if (sup_ven && sup_wm == 0)
        {
          min_dist = label_dist+1 ;
        }

        if (min_dist > label_dist)  /* couldn't find any labels that match */
        {
          double log_p, max_log_p ;

          /* wm may be partial volumed - look in smaller
             nbhd for most likely location of wm */
          min_dist = 0 ;
          max_log_p = -1e20 ;
          for (yk = -label_dist/3 ; yk <= label_dist/3 ; yk += SAMPLE_DIST)
          {
            yi = gcamn->y+yk ;
            if ((yi >= (mri->height-1)) || (yi <= 0))
            {
              break ;
            }

            load_vals(mri, gcamn->x, yi, gcamn->z, vals, gcam->ninputs) ;
            log_p = GCAcomputeConditionalLogDensity
                    (wm_gc, vals, gcam->ninputs, wm_label);
            if (log_p > max_log_p)
            {
              max_log_p = log_p ;
              min_dist = yk ;
            }
          }
        }
        else if (!ctrl_point_found)   // adjust estimated position to be at most likely value
        {
          double log_p, max_log_p ;
          double  ykmin, ykmax ;

          /* wm may be partial volumed - look in smaller
             nbhd for most likely location of wm */
          max_log_p = -1e20 ;
          ykmin = min_dist-1 ;
          ykmax = min_dist+1 ;
          for (yk = ykmin ; yk <= ykmax ; yk += SAMPLE_DIST)
          {
            yi = gcamn->y+yk ;   /* sample inferiorly */
            if ((yi >= (mri->height-1)) || (yi <= 0))
            {
              break ;
            }
            load_vals(mri, gcamn->x, yi, gcamn->z, vals, gcam->ninputs) ;
            log_p =
              GCAcomputeConditionalLogDensity
              (wm_gc, vals, gcam->ninputs, wm_label);
            if (mri_twm && MRIgetVoxVal(mri_twm, nint(gcamn->x), nint(yi), nint(gcamn->z), 0) > 0)
            {
              ctrl_point_found = 1 ;
              log_p = 0 ;   // manually specified temporal lobe wm point
            }
            if (log_p > max_log_p)
            {
              max_log_p = log_p ;
              min_dist = yk ;
            }
          }
        }

        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (x == Gx && y == (Gy-1) && z == Gz)
        {
          DiagBreak() ;
        }
        if (x == Gx && y == (Gy+1) && z == Gz)
        {
          DiagBreak() ;
        }

        gcamn->label_dist = MRIFvox(mri_dist, x, y, z) = min_dist ;
#if 1
        if (fabs(min_dist) > MAX_MLE_DIST && !ctrl_point_found)
        {
          min_dist = MAX_MLE_DIST * min_dist / fabs(min_dist) ;
        }
#endif
        dy = min_dist ;
        if (!FZERO(min_dist))
        {
          gcamn->status = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ;
          if (ctrl_point_found)
          {
            gcamn->status |= GCAM_MANUAL_LABEL ;
          }
          if (IS_WM(gcamn_inf->label) &&
              ((gcamn_inf->status & GCAM_LABEL_NODE)==0))
          {
            gcamn_inf->status = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ;
            gcamn_inf->dy += (l_label)*dy ;
            if (ctrl_point_found)
            {
              gcamn_inf->status |= GCAM_MANUAL_LABEL ;
            }
            gcamn_inf->label_dist =
              MRIFvox(mri_dist, x, y+1, z) = gcamn->label_dist ;
            if (x == Gx && (y+1) == Gy && z == Gz)
              printf("l_label: node(%d,%d,%d): dy = %2.2f\n",
                     x, y+1, z, gcamn_inf->dy)  ;
          }
          if (IS_HIPPO(gcamn_sup->label) &&
              ((gcamn_sup->status & GCAM_LABEL_NODE)==0))
          {
            gcamn_sup->status = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ;
            if (ctrl_point_found)
            {
              gcamn_sup->status |= GCAM_MANUAL_LABEL ;
            }
            gcamn_sup->dy += (l_label)*dy ;
            gcamn_sup->label_dist =
              MRIFvox(mri_dist, x, y-1, z) = gcamn->label_dist ;
            if (x == Gx && (y-1) == Gy && z == Gz)
            {
              printf("l_label: node(%d,%d,%d): dy = %2.2f\n", x, y-1, z, gcamn_sup->dy)  ;
            }
          }
        }
        GCAfreeRegionalGCAN(&gcan) ;
      }
    }
  }
#endif
  DiagBreak() ;
}


// ----------------------

#define GCAM_LABEL_TERM_TIMERS 0

void SetInconsistentLabelNodes( const int val )
{
  /*!
    This is a call back from the GPU version of Label Term.
    I'm not sure if the inconsistentLabelNodes global is
    ever actually used.
  */
  inconsistentLabelNodes = val;
}



#define GCAM_LABEL_TERM_OUTPUT 0

int
gcamLabelTerm( GCA_MORPH *gcam, const MRI *mri,
               double l_label, double label_dist , MRI *mri_twm)
{

#ifdef GCAM_LABEL_TERM_MAINLOOP_GPU

  printf( "%s: On GPU\n", __FUNCTION__ );
  gcamLabelTermGPU( gcam, mri, l_label, label_dist );


#else

  int num;
  MRI *mri_dist ;
  int nremoved ;

  if( DZERO(l_label) )
  {
    return(NO_ERROR) ;
  }

#if GCAM_LABEL_TERM_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( (nCalls%outputFreq) == 0 )
  {
    unsigned int nOut = nCalls/outputFreq;

    char fname[STRLEN];

    snprintf( fname, STRLEN-1, "gcamLabelTermInput%04u", nOut );
    fname[STRLEN-1] = '\0';
    WriteGCAMoneInput( gcam, fname );

    snprintf( fname, STRLEN-1, "mriLabelTermInput%04u.mgz", nOut );
    MRIwrite( (MRI*)mri, fname );

    snprintf( fname, STRLEN-1, "gcaLabelTerm%04u.gca", nOut );
    GCAwrite( gcam->gca, fname );
  }
  nCalls++;
#endif

  mri_dist = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
  MRIsetResolution(mri_dist, gcam->spacing, gcam->spacing, gcam->spacing) ;

  GCAMresetLabelNodeStatus( gcam ) ;

  gcamLabelTermMainLoop( gcam, mri, mri_dist, l_label, label_dist , mri_twm);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_dist, "dist_before.mgz") ;
  }

  /* do neighborhood consistency check */
  nremoved = remove_label_outliers(gcam, mri_dist, 2, 3) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("%d inconsistent label nodes removed...\n", nremoved) ;
  }

  // Copy to some random global (which I think is actually unused)
  inconsistentLabelNodes=nremoved;

  gcamLabelTermCopyDeltas( gcam, mri_dist, l_label );

  /*
    The following gives different results on the CPU and GPU.
    It is built around a very nasty, order dependent loop.
  */
  nremoved += gcamLabelTermPostAntConsistency( gcam, mri_dist );

  num = gcamLabelTermFinalUpdate( gcam, mri_dist, l_label );

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_dist, "dist_after.mgz") ;
  }

  MRIfree(&mri_dist) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("\t%d nodes for which label term applies\n", num) ;
  }


#endif

  return(NO_ERROR) ;
}


int
gcamMapTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_map)
{
  int             x=0, y=0, z=0, n=0, i=0 ;
  int     xn[_MAX_FS_THREADS], yn[_MAX_FS_THREADS], zn[_MAX_FS_THREADS];
  int             tid;
  double          node_prob=0.0, prob=0.0, dx=0.0, dy=0.0, dz=0.0, norm=0.0 ;
  float           vals[_MAX_FS_THREADS][MAX_GCA_INPUTS] ;
  GCA_MORPH_NODE  *gcamn=NULL ;
  GCA_PRIOR       *gcap=NULL ;
  GCA_NODE        *gcan=NULL ;
  GC1D            *gc=NULL ;
  MATRIX          *m_delI[_MAX_FS_THREADS], *m_inv_cov[_MAX_FS_THREADS] ;
  VECTOR          *v_means[_MAX_FS_THREADS], *v_grad[_MAX_FS_THREADS] ;

  if (DZERO(l_map))
  {
    return(0) ;
  }
  // 3 x ninputs
#ifdef HAVE_OPENMP
  #pragma omp parallel
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  m_delI[tid] = MatrixAlloc(3, gcam->ninputs, MATRIX_REAL) ;
  // ninputs x ninputs
  m_inv_cov[tid] = MatrixAlloc(gcam->ninputs, gcam->ninputs, MATRIX_REAL) ;
  // ninputs x 1
  v_means[tid] = VectorAlloc(gcam->ninputs, 1) ;
  // 3 x 1
  v_grad[tid] = VectorAlloc(3, MATRIX_REAL) ;

#ifdef HAVE_OPENMP
  #pragma omp parallel for firstprivate(tid,i,y,z,gcamn,gcan,gcap,n,norm,dx,dy,dz,gc,node_prob,prob) shared (gcam,Gx,Gy,Gz,mri_smooth,l_map) schedule(static,1)
#endif
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        if (gcamn->status & (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD))
        {
          continue ;
        }
        /////////////////////
        if (!GCApriorToNode(gcam->gca, x, y, z, &xn[tid], &yn[tid], &zn[tid]))
        {
          gcan = &gcam->gca->nodes[xn[tid]][yn[tid]][zn[tid]] ;
          gcap = &gcam->gca->priors[x][y][z] ;
          // get the values from mri
          load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals[tid], gcam->ninputs) ;

          for (n = 0 ; n < gcam->ninputs ; n++)
          {
            // get dx, dy, dz
            MRIsampleVolumeGradientFrame(mri_smooth, gcamn->x, gcamn->y, gcamn->z, \
                                         &dx, &dy, &dz, n) ;
            // magnitude
            norm = sqrt(dx*dx+dy*dy+dz*dz) ;
            // non-zero, then normalize
            if (!FZERO(norm))  /* don't worry about magnitude of gradient */
            {
              dx /= norm ;
              dy /= norm ;
              dz /= norm ;
            }
            // store   3 x ninputs
            *MATRIX_RELT(m_delI[tid], 1, n+1) = dx ;
            *MATRIX_RELT(m_delI[tid], 2, n+1) = dy ;
            *MATRIX_RELT(m_delI[tid], 3, n+1) = dz ;
          }

          if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
          {
            printf("l_map: node(%d,%d,%d), delI=(%2.2f, %2.2f, %2.2f)\n", x, y, z, dx, dy, dz) ;
          }

          dx = dy = dz = 0.0f ;
          for (node_prob = 0.0, n = 0 ; n < gcap->nlabels ; n++)
          {
            gc = GCAfindGC(gcam->gca, xn[tid], yn[tid], zn[tid], gcap->labels[n]) ;
            if (!gc)
            {
              continue ;
            }
            // mean
            load_mean_vector(gc, v_means[tid], gcam->ninputs) ;
            // inv_cov
            load_inverse_covariance_matrix(gc, m_inv_cov[tid], gcam->ninputs) ;
            // get prob
            prob = GCAcomputeConditionalDensity(gc, vals[tid], gcam->ninputs, gcap->labels[n]) ;
            // v_mean = mean - vals
            for (i = 0 ; i < gcam->ninputs ; i++)
            {
              VECTOR_ELT(v_means[tid], i+1) -= vals[tid][i] ;
            }
            // v_mean = inv_cov * (mean - vals)
            MatrixMultiply(m_inv_cov[tid], v_means[tid], v_means[tid]) ;
            // v_grad = delI * inv_cov * (mean - value)
            MatrixMultiply(m_delI[tid], v_means[tid], v_grad[tid]) ;

            if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
              printf("l_map: node(%d,%d,%d), label %s: p=%2.3f (%2.3f), D=(%2.1f,%2.1f,%2.1f)\n",
                     x, y, z, cma_label_to_name(gcap->labels[n]), prob, gcap->priors[n],
                     prob*V3_X(v_grad[tid]), prob*V3_Y(v_grad[tid]), prob*V3_Z(v_grad[tid])) ;

            dx += prob*V3_X(v_grad[tid]) ;
            dy += prob*V3_Y(v_grad[tid]) ;
            dz += prob*V3_Z(v_grad[tid]) ;

            node_prob += prob ;
          }

          if (DZERO(node_prob))
          {
            node_prob = 1e-6 ;
          }

          dx *= (1.0/node_prob) ;
          dy *= (1.0/node_prob) ;
          dz *= (1.0/node_prob) ;
          if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
            printf("l_map: node(%d,%d,%d), 1/p=%2.3f, grad=(%2.2f,%2.2f,%2.2f)\n",
                   x, y, z, 1.0/node_prob, dx, dy, dz) ;

          gcamn->dx += l_map * dx ;
          gcamn->dy += l_map * dy ;
          gcamn->dz += l_map * dz ;
        } //!GCA
      }

  MatrixFree(&m_delI[tid]) ;
  MatrixFree(&m_inv_cov[tid]) ;
  VectorFree(&v_means[tid]) ;
  VectorFree(&v_grad[tid]) ;
  return(NO_ERROR) ;
}


double
gcamMapEnergy(GCA_MORPH *gcam, MRI *mri)
{
  int             x, y, z, n, xn, yn, zn ;
  double          sse, node_prob, prob ;
  GCA_MORPH_NODE  *gcamn ;
  GCA_PRIOR       *gcap ;
  GCA_NODE        *gcan ;
  GC1D            *gc ;
  float            vals[MAX_GCA_INPUTS] ;

  sse = 0.0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        if (gcamn->status & (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD))
        {
          continue ;
        }
        if (!GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn))
        {
          gcan = &gcam->gca->nodes[xn][yn][zn] ;
          gcap = &gcam->gca->priors[x][y][z] ;
          load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->ninputs) ;
          for (node_prob = 0.0, n = 0 ; n < gcap->nlabels ; n++)
          {
            gc = GCAfindGC(gcam->gca, xn, yn, zn, gcap->labels[n]) ;
            if (!gc)
            {
              continue ;
            }

            prob = GCAcomputeConditionalDensity(gc, vals, gcam->ninputs, gcap->labels[n]) ;
            if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
              printf("E_map: node(%d,%d,%d), label %s: p=%2.3f (%2.3f)\n", x, y, z,
                     cma_label_to_name(gcap->labels[n]), prob, gcap->priors[n]) ;
            node_prob += prob ;
          }

          if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
          {
            printf("E_map: node(%d,%d,%d), -log(p)=%2.3f\n", x, y, z, -log(node_prob)) ;
          }
          if (FZERO(node_prob))
          {
            node_prob = 1e-6 ;  /* something big?? */
          }
          sse -= log(node_prob) ;
        } // !GCA
      }

  return(sse) ;
}

#if 0
int
gcamMLElabelAtLocation(GCA_MORPH *gcam, int x, int y, int z, float *vals)
{
  int   max_n ;
  float log_p ;

  return(GCAcomputeMLElabelAtLocation(gcam->gca, x, y, z, vals, &max_n, &log_p)) ;
}
#endif

int
GCAMstoreMetricProperties(GCA_MORPH *gcam)
{
  int             x, y, z ;

  gcamComputeMetricProperties(gcam) ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcam->nodes[x][y][z].orig_area = gcam->nodes[x][y][z].area ;
        gcam->nodes[x][y][z].orig_area1 = gcam->nodes[x][y][z].area1 ;
        gcam->nodes[x][y][z].orig_area2 = gcam->nodes[x][y][z].area2 ;
        if ((FZERO(gcam->nodes[x][y][z].orig_area) ||
             (gcam->nodes[x][y][z].orig_area < 0)) &&
            (gcam->nodes[x][y][z].invalid == 0))
        {
          gcam->nodes[x][y][z].invalid = GCAM_AREA_INVALID ;
        }
      }
  return(NO_ERROR) ;
}

int
GCAMcomputeOriginalProperties(GCA_MORPH *gcam)
{
  GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, SAVED_POSITIONS) ;
  GCAMcopyNodePositions(gcam, ORIGINAL_POSITIONS, CURRENT_POSITIONS) ;
  gcamComputeMetricProperties(gcam) ;
  GCAMstoreMetricProperties(gcam) ;
  GCAMcopyNodePositions(gcam, SAVED_POSITIONS, CURRENT_POSITIONS) ;
  gcamComputeMetricProperties(gcam) ;
  return(NO_ERROR) ;
}

int
GCAMcopyNodePositions(GCA_MORPH *gcam, int from, int to)
{
#ifdef GCAM_COPY_NODE_POSITIONS_GPU
  printf( "%s: On GPU\n", __FUNCTION__ );
  GCAMcopyNodePositionsGPU( gcam, from, to );
#else
  int             x, y, z ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        switch (from)
        {
        default:
          ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, \
                                      "GCAMcopyNodePositions: unsupported from %d",from)) ;
        case ORIGINAL_POSITIONS:
          switch (to)
          {
          default:
            ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, \
                                        "GCAMcopyNodePositions: unsupported to %d",to)) ;
          case SAVED_ORIGINAL_POSITIONS:
            gcamn->saved_origx = gcamn->origx ;
            gcamn->saved_origy = gcamn->origy ;
            gcamn->saved_origz = gcamn->origz ;
            break ;
          case SAVED_POSITIONS:
            gcamn->xs = gcamn->origx ;
            gcamn->ys = gcamn->origy ;
            gcamn->zs = gcamn->origz ;
            break ;
          case CURRENT_POSITIONS:
            gcamn->x = gcamn->origx ;
            gcamn->y = gcamn->origy ;
            gcamn->z = gcamn->origz ;
            break ;

          }
          break ;

        case SAVED_ORIGINAL_POSITIONS:
          switch (to)
          {
          default:
            ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, \
                                        "GCAMcopyNodePositions: unsupported to %d",to)) ;
          case SAVED_POSITIONS:
            gcamn->xs = gcamn->saved_origx ;
            gcamn->ys = gcamn->saved_origy ;
            gcamn->zs = gcamn->saved_origz ;
            break ;
          case CURRENT_POSITIONS:
            gcamn->x = gcamn->saved_origx ;
            gcamn->y = gcamn->saved_origy ;
            gcamn->z = gcamn->saved_origz ;
            break ;
          case ORIGINAL_POSITIONS:
            gcamn->origx = gcamn->saved_origx ;
            gcamn->origy = gcamn->saved_origy ;
            gcamn->origz = gcamn->saved_origz ;
            break ;
          }
          break ;

        case SAVED_POSITIONS:
          switch (to)
          {
          default:
            ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, \
                                        "GCAMcopyNodePositions: unsupported to %d",to)) ;
          case ORIGINAL_POSITIONS:
            gcamn->origx = gcamn->xs ;
            gcamn->origy = gcamn->ys ;
            gcamn->origz = gcamn->zs ;
            break ;
          case CURRENT_POSITIONS:
            gcamn->x = gcamn->xs ;
            gcamn->y = gcamn->ys ;
            gcamn->z = gcamn->zs ;
            break ;
          case SAVED_ORIGINAL_POSITIONS:
            gcamn->saved_origx = gcamn->xs ;
            gcamn->saved_origy = gcamn->ys ;
            gcamn->saved_origz = gcamn->zs ;
            break ;
          }
          break ;

        case SAVED2_POSITIONS:
          switch (to)
          {
          default:
            ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, \
                                        "GCAMcopyNodePositions: unsupported to %d",to)) ;
          case ORIGINAL_POSITIONS:
            gcamn->origx = gcamn->xs2 ;
            gcamn->origy = gcamn->ys2 ;
            gcamn->origz = gcamn->zs2 ;
            break ;
          case CURRENT_POSITIONS:
            gcamn->x = gcamn->xs2 ;
            gcamn->y = gcamn->ys2 ;
            gcamn->z = gcamn->zs2 ;
            break ;
          case SAVED_ORIGINAL_POSITIONS:
            gcamn->saved_origx = gcamn->xs2 ;
            gcamn->saved_origy = gcamn->ys2 ;
            gcamn->saved_origz = gcamn->zs2 ;
            break ;
          }
          break ;

        case CURRENT_POSITIONS:
          switch (to)
          {
          default:
            ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, \
                                        "GCAMcopyNodePositions: unsupported to %d",to)) ;
          case ORIGINAL_POSITIONS:
            gcamn->origx = gcamn->x ;
            gcamn->origy = gcamn->y ;
            gcamn->origz = gcamn->z ;
            break ;
          case SAVED_POSITIONS:
            gcamn->xs = gcamn->x ;
            gcamn->ys = gcamn->y ;
            gcamn->zs = gcamn->z ;
            break ;
          case SAVED2_POSITIONS:
            gcamn->xs2 = gcamn->x ;
            gcamn->ys2 = gcamn->y ;
            gcamn->zs2 = gcamn->z ;
            break ;
          case SAVED_ORIGINAL_POSITIONS:
            gcamn->saved_origx = gcamn->x ;
            gcamn->saved_origy = gcamn->y ;
            gcamn->saved_origz = gcamn->z ;
            break ;
          }

          break ;

        }
      }

#endif
  return(NO_ERROR) ;
}

int
gcamRemoveCompressedNodes(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms,
                          float compression_ratio)
{
  GCA_MORPH_PARMS saved_parms = *parms ;
  double          min_dt, orig_dt, rms, last_rms, pct_change, max_grad, lattice_spacing ;
  int             old_compressed, new_compressed, i, delta, integration_type, done, nsmall, good_step, nfixed = 0, start_neg ;

  new_compressed = GCAMcountCompressedNodes(gcam, compression_ratio) ;
  if (new_compressed <= 0)
  {
    return(NO_ERROR) ;
  }

  parms->l_distance = parms->l_log_likelihood = parms->l_spring = parms->l_multiscale =
                        parms->l_area_intensity = parms->l_expansion = parms->l_area_smoothness =
                              parms->l_binary = parms->l_area = parms->l_elastic = parms->l_smoothness = parms->l_label = 0.0 ;
  parms->navgs = 0 ;
  parms->l_area = 0.0 ;
  parms->l_jacobian = 1 ;
  parms->dt = 0.1 ;
  parms->ratio_thresh = compression_ratio ;
  parms->exp_k = 5 ;
  parms->integration_type = GCAM_INTEGRATE_BOTH ;  // iterate between fixed and optimal
  integration_type = GCAM_INTEGRATE_OPTIMAL ;      // start with optimal
  parms->navgs = 0 ;
  parms->dt = orig_dt = 1e-6 ;
  start_neg = gcam->neg ;

  last_rms = rms = GCAMcomputeRMS(gcam, mri, parms) ;
  printf(" starting rms=%2.5f, compressed=%d, removing compressions in lattice....\n", rms, new_compressed) ;
  if (Gdiag & DIAG_SHOW)
  {
    gcamShowCompressed(gcam, stdout) ;
  }
  if (parms->log_fp)
  {
    fprintf(parms->log_fp, "starting: ") ;
    gcamShowCompressed(gcam, parms->log_fp) ;
  }
  done = good_step = nsmall = i = 0 ;
  lattice_spacing = MIN(gcam->atlas.xsize, MIN(gcam->atlas.ysize,gcam->atlas.zsize)) ;
  do
  {
    old_compressed = new_compressed ;
    gcamComputeGradient(gcam, mri, mri, parms) ;
    gcamSmoothGradient(gcam, parms->navgs) ;
    max_grad = gcamMaxGradient(gcam) ;
    if (max_grad < 0.01 && new_compressed > 10)
    {
      DiagBreak() ;
    }
    if (DZERO(max_grad))
    {
      max_grad = 1 ;
    }
    switch (integration_type)
    {
    case GCAM_INTEGRATE_OPTIMAL:
      parms->dt = .1*lattice_spacing;
      parms->dt = MIN(parms->dt, 0.01*lattice_spacing/max_grad) ;
      min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
      break ;
    default:
    case GCAM_INTEGRATE_FIXED:
      nfixed++ ;
      min_dt = .1*lattice_spacing/(2*max_grad) ;
      min_dt = .1*lattice_spacing ;
      if (min_dt < orig_dt)
      {
        min_dt = orig_dt ;
      }
      break ;
    }
    parms->dt = min_dt ;
    gcamApplyGradient(gcam, parms) ;
    last_rms = rms ;
    rms = GCAMcomputeRMS(gcam, mri, parms) ;
    new_compressed = GCAMcountCompressedNodes(gcam, compression_ratio) ;
    if (gcam->neg > 0 && start_neg == 0 && parms->noneg == True)
    {
      gcamUndoGradient(gcam) ;
      rms = last_rms ;
      new_compressed = old_compressed ;
    }
    parms->dt = orig_dt ;
    pct_change = 100.0*(last_rms-rms)/(last_rms) ;
    delta = old_compressed-new_compressed ;

    printf("  iter %d, %c dt=%2.6f: new compressed %d, old_compressed %d, delta %d, rms=%2.5f\n", \
           ++i, integration_type == GCAM_INTEGRATE_OPTIMAL ? 'O' : 'F',
           min_dt, new_compressed, old_compressed, old_compressed-new_compressed, rms) ;
    if (Gdiag & DIAG_SHOW)
    {
      gcamShowCompressed(gcam, stdout) ;
    }
    if (new_compressed == 0)
    {
      break ;
    }
    // iterate until compressed nodes don't decrease, but always for at least 3 iters
    if (pct_change < parms->tol && delta <= 0)
    {
      if (integration_type == GCAM_INTEGRATE_FIXED)
      {
        if (nsmall++ > 1)
        {
          if (good_step == 0) // couldn't take any steps
          {
            done = 1 ;
          }
          integration_type = GCAM_INTEGRATE_OPTIMAL ;
          good_step = 0 ;
          printf("switching integration type to optimal\n") ;
        }
      }
      else
      {
        printf("switching integration type to fixed\n") ;
        nfixed = 0 ;
        integration_type = GCAM_INTEGRATE_FIXED ;
      }
    }
    else
    {
      nsmall = 0 ;
      good_step = 1 ;
    }
    if (nfixed > 10)
    {
      printf("switching integration type to optimal\n") ;
      good_step = 0 ;
      nfixed = 0 ;
      integration_type = GCAM_INTEGRATE_OPTIMAL ;
    }
  }
  while (done == 0 && i < 150) ;

  if (parms->log_fp)
  {
    fprintf(parms->log_fp, "ending %d: ", i) ;
    if (gcamShowCompressed(gcam, parms->log_fp) == 0)
    {
      fprintf(parms->log_fp, "\n") ;
    }
    fflush(parms->log_fp) ;
  }
  *parms = *(&saved_parms) ;
  return(NO_ERROR) ;
}


#define MAX_NEG_ITER 100

int
gcamRemoveNegativeNodes(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  GCA_MORPH_PARMS saved_parms = *parms ;
  double          min_dt, orig_dt = parms->orig_dt, rms, last_rms, pct_change ;
  int             old_neg, new_neg, i ;
  MRI             *mri_warp = NULL ;

  if (gcam->neg <= 0)
  {
    return(NO_ERROR) ;
  }

  mri_warp = GCAMwriteWarpToMRI(gcam, mri_warp) ;
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri_warp) ;
  MRIfree(&mri_warp) ;

  if (gcam->neg <= 0)
  {
    return(NO_ERROR) ;
  }
  parms->noneg = 0 ;
  parms->l_distance = parms->l_log_likelihood = parms->l_binary = parms->l_multiscale =
                        parms->l_spring = parms->l_area = parms->l_elastic = parms->l_smoothness = parms->l_label = 0.0 ;
  parms->navgs = 0 ;
  parms->l_area = 0.0 ;
  parms->l_jacobian = 1 ;
  parms->dt = 0.1 ;
  parms->tol = 0.01 ;

  last_rms = rms = GCAMcomputeRMS(gcam, mri, parms) ;
  printf("starting rms=%2.3f, neg=%d, removing folds in lattice....\n", rms, gcam->neg) ;
  new_neg = gcam->neg ;
  i = 0 ;
  do
  {
    old_neg = new_neg ;
    gcamComputeGradient(gcam, mri, mri, parms) ;
    min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
    parms->dt = min_dt ;
    gcamApplyGradient(gcam, parms) ;
    last_rms = rms ;
    rms = GCAMcomputeRMS(gcam, mri, parms) ;
    parms->dt = orig_dt ;
    new_neg = gcam->neg ;
    pct_change = 100.0*(last_rms-rms)/(last_rms) ;
    printf("iter %d, dt=%2.6f: new neg %d, old_neg %d, delta %d, rms=%2.3f (%2.3f%%)\n",
           ++i, min_dt, new_neg, old_neg, old_neg-new_neg, rms, pct_change) ;
  }
  while ((new_neg > 0) && (pct_change > parms->tol) && (i < MAX_NEG_ITER)) ;

  *parms = *(&saved_parms) ;
  return(NO_ERROR) ;
}



int
check_gcam( const GCAM *gcam )
{
  if ((abs(gcam->width) > 10000) ||
      (abs(gcam->height) > 10000) ||
      (abs(gcam->depth) > 10000))
  {
    DiagBreak() ;
    printf("gcam dimensions invalid - (%d, %d, %d)\n", gcam->width, gcam->height, gcam->depth) ;
    return(ERROR_BADPARM) ;
  }
#if 0
  int  x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (z = 0 ; z < gcam->depth ; z++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (x = 0 ; x < gcam->width ; x++)
      {
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }

        if (gcamn->gc)
        {
          gcamn->gc->means[0] = gcamn->gc->means[0]*1;
        }
      }
    }
  }
#endif

  return(0) ;
}

#define MAX_TEMPORAL_WM 3

int
is_temporal_wm( const GCA_MORPH *gcam, const MRI *mri,
                const GCA_NODE *gcan,
                float xf, float yf, float zf, int ninputs )
{
  int  yk, label, nwhite ;
  float vals[MAX_GCA_INPUTS], yi ;


  nwhite = 0 ;
  for (yk = 1 ; yk < 2*MAX_TEMPORAL_WM ; yk++)
  {
    yi = yf-yk ;
    if (yi < 0)
    {
      break ;
    }
    load_vals(mri, xf, yi, zf, vals, gcam->ninputs) ;
    label = GCAmaxLikelihoodLabel(gcan, vals, gcam->ninputs, NULL) ;
    if (IS_WM(label) || IS_THALAMUS(label))
    {
      nwhite++ ;
    }
  }
  for (yk = -1 ; yk > -2*MAX_TEMPORAL_WM ; yk--)
  {
    yi = yf-yk ;
    if (yi < 0 || yi >= mri->height)
    {
      break ;
    }
    load_vals(mri, xf, yi, zf, vals, gcam->ninputs) ;
    label = GCAmaxLikelihoodLabel(gcan, vals, gcam->ninputs, NULL) ;
    if (IS_WM(label) || IS_THALAMUS(label))
    {
      nwhite++ ;
    }
    else
    {
      break ;  /* if moving inferiorly and find non-wm voxel, then stop counting - should be hippo */
    }
  }

  return(nwhite <= MAX_TEMPORAL_WM) ;  /* must be main body of white matter - too much of it */
}

int
GCAMapplyTransform(GCA_MORPH *gcam, TRANSFORM *transform)
{
  int            x, y, z ;
  float          xf, yf, zf ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        TransformSample(transform, (float)gcamn->x, (float)gcamn->y, (float)gcamn->z,
                        &xf, &yf, &zf) ;
        gcamn->x = xf ;
        gcamn->y = yf ;
        gcamn->z = zf ;

        TransformSample(transform, (float)gcamn->origx, (float)gcamn->origy, (float)gcamn->origz,
                        &xf, &yf, &zf) ;
        gcamn->origx = xf ;
        gcamn->origy = yf ;
        gcamn->origz = zf ;
      }
    }
  }
  if (transform->type == MORPH_3D_TYPE)
  {
    GCA_MORPH_PARMS parms ;
    MRI             *mri_tmp ;

    GCAMcomputeOriginalProperties(gcam) ;
    gcamComputeMetricProperties(gcam) ;
    GCAMmarkNegativeNodesInvalid(gcam) ;
    gcamComputeMetricProperties(gcam) ;
    mri_tmp = MRIalloc(gcam->image.width, gcam->image.height, gcam->image.depth, MRI_UCHAR);
    memset(&parms, 0, sizeof(parms)) ;
    gcamRemoveNegativeNodes(gcam, mri_tmp, &parms) ;
    MRIfree(&mri_tmp) ;
  }
  return(NO_ERROR) ;
}
int
GCAMapplyInverseTransform(GCA_MORPH *gcam, TRANSFORM *transform)
{
  int            x, y, z ;
  float          xf, yf, zf ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        TransformSampleInverseFloat(transform, (float)gcamn->x, (float)gcamn->y, (float)gcamn->z,
                                    &xf, &yf, &zf) ;
        gcamn->x = xf ;
        gcamn->y = yf ;
        gcamn->z = zf ;

        TransformSampleInverseFloat(transform, (float)gcamn->origx, (float)gcamn->origy, (float)gcamn->origz,
                                    &xf, &yf, &zf) ;
        gcamn->origx = xf ;
        gcamn->origy = yf ;
        gcamn->origz = zf ;
      }
    }
  }
  if (transform->type == MORPH_3D_TYPE)
  {
    GCA_MORPH_PARMS parms ;
    MRI             *mri_tmp ;

    GCAMcomputeOriginalProperties(gcam) ;
    gcamComputeMetricProperties(gcam) ;
    GCAMmarkNegativeNodesInvalid(gcam) ;
    gcamComputeMetricProperties(gcam) ;
    mri_tmp = MRIalloc(gcam->image.width, gcam->image.height, gcam->image.depth, MRI_UCHAR);
    memset(&parms, 0, sizeof(parms)) ;
    gcamRemoveNegativeNodes(gcam, mri_tmp, &parms) ;
    MRIfree(&mri_tmp) ;
  }
  return(NO_ERROR) ;
}
int
GCAMmarkNegativeNodesInvalid(GCA_MORPH *gcam)
{
  int  x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (z = 0 ; z < gcam->depth ; z++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (x = 0 ; x < gcam->width ; x++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if ((gcamn->area <= 0 || gcamn->area1 <= 0 || gcamn->area2 <= 0 || gcamn->orig_area <= 0) &&
            (gcamn->invalid == 0))
        {
          gcamn->invalid = GCAM_AREA_INVALID ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}
int
GCAMmarkNegativeBorderNodesInvalid(GCA_MORPH *gcam)
{
  int  x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (z = 0 ; z < gcam->depth ; z++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (x = 0 ; x < gcam->width ; x++)
      {
        if (x != 0 && y != 0 && z != 0 &&
            x != gcam->width-1 &&
            y != gcam->height-1 &&
            z != gcam->depth-1)
        {
          continue ;  // must be on at least one border plane
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if ((gcamn->area <= 0 || gcamn->area1 <= 0 || gcamn->area2 <= 0 || gcamn->orig_area <= 0) &&
            (gcamn->invalid == 0))
        {
          gcamn->invalid = GCAM_AREA_INVALID ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

int
zero_vals(float *vals, int nvals)
{
  int n, z ;

  for (z = 1, n = 0 ; n < nvals ; n++)
    if (!FZERO(vals[n]))
    {
      z = 0 ;
      break ;
    }

  return(z) ;
}

int
different_neighbor_labels( const GCA_MORPH *gcam,
                           const int x, const int y, const int z,
                           const int whalf )
{
  int        label, num, i, j, k ;

  label = gcam->nodes[x][y][z].label ;
  for (num = 0, i = x-whalf ; i <= x+whalf ; i++)
  {
    if (i < 0 || i >= gcam->width)
    {
      continue ;
    }
    for (j = y-whalf ; j <= y+whalf ; j++)
    {
      if (j < 0 || j >= gcam->height)
      {
        continue ;
      }
      for (k = z-whalf ; k <= z+whalf ; k++)
      {
        if (k < 0 || k >= gcam->height)
        {
          continue ;
        }
        if (i == 0 && j == 0 && k == 0)
        {
          continue ;
        }
        if (label != gcam->nodes[i][j][k].label)
        {
          num++ ;
        }
      }
    }
  }
  return(num) ;
}


MRI *
gcamCreateJacobianImage(GCA_MORPH *gcam)
{
  GCA_MORPH_NODE *gcamn ;
  MRI            *mri ;
  int            x, y, z ;
  float          jacobian ;

  gcamComputeMetricProperties(gcam) ;
  mri = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
  MRIsetResolution(mri, gcam->spacing, gcam->spacing, gcam->spacing) ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (!FZERO(gcamn->area))
        {
          jacobian = gcamn->orig_area / gcamn->area ;
        }
        else
        {
          if (FZERO(gcamn->orig_area))
          {
            jacobian = 1 ;
          }
          else
          {
            jacobian = .001 ;
          }
        }

        if (x == Gx && y == Gy && z == Gz)
        {
          printf("jacobian(%d, %d, %d) = %2.3f\n", x, y, z, jacobian) ;
        }
        MRIsetVoxVal(mri, x, y, z, 0, jacobian) ;
      }
    }
  }

  return(mri) ;
}

MRI *
GCAMextract_density_map(MRI *mri_seg, MRI *mri_intensity, GCA_MORPH *gcam, int target_label, MRI *mri_out)
{
  int    whalf, x, y, z /*, xv, yv, zv*/ ;
  //  float  xa, ya, za ;
  MRI    *mri_density ;
  double   volume, det ;
  GCA_MORPH_NODE *gcamn ;
  MATRIX *m_L = NULL ;

  whalf = 1 ;
  if (mri_out == NULL)
  {
    mri_out = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
    MRIcopyVolGeomToMRI(mri_out, &gcam->atlas) ;
    MRIsetResolution(mri_out, gcam->spacing, gcam->spacing, gcam->spacing) ;
    MRIreInitCache(mri_out) ;
  }
  if (mri_out->type != MRI_FLOAT)
    ErrorReturn(NULL, (ERROR_BADPARM, \
                       "GCAMextract_density_map: output volume type must be MRI_FLOAT (%d)",
                       mri_out->type)) ;


  /* go through each source voxel.
     If it is the right label, partial volume correct it, then map it into the
     destination space, and trilinearly interpolate it into that volume
  */


  mri_density = MRImakeDensityMap(mri_seg, mri_intensity, target_label, NULL,
                                  mri_seg->xsize) ;
  if (Gx >= 0)
  {
    printf("density at (%d, %d, %d) = %2.3f\n", Gx, Gy, Gz, MRIgetVoxVal(mri_density, Gx, Gy, Gz, 0)) ;
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_density, "density.mgz") ;
  }

#if 0
  GCAMinvert(gcam, mri_seg) ;
  for (x = 0  ; x < mri_seg->width ; x++)
  {
    for (y = 0  ; y < mri_seg->height ; y++)
    {
      for (z = 0  ; z < mri_seg->width ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        volume = MRIgetVoxVal(mri_density, x, y, z, 0) ;
        GCAMsample(gcam, x, y, z, &xa, &ya, &za) ; /* convert to atlas coordinates */
        if (!finitep(volume))
        {
          DiagBreak() ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("voxel(%d, %d, %d) maps to atlas (%2.1f, %2.1f, %2.1f)\n",
                 x, y, z, xa, ya, za) ;
#if 0
        MRIinterpolateIntoVolume(mri_out, (double)xa, (double)ya, (double)za, (double)volume) ;
#else
        xv = nint(xa) ;
        yv = nint(ya) ;
        zv = nint(za) ;
        if (xv == Gx && yv == Gy && zv == Gz)
        {
          DiagBreak() ;
        }
        MRIsetVoxVal(mri_out, xv, yv, zv,0, MRIgetVoxVal(mri_out,xv,yv,zv,0)+volume) ;
#endif
      }
    }
  }
#else
  for (x = 0  ; x < mri_out->width ; x++)
  {
    for (y = 0  ; y < mri_out->height ; y++)
    {
      for (z = 0  ; z < mri_out->width ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (!finitep(MRIFvox(mri_out, x, y, z)))
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolume(mri_density, gcamn->x, gcamn->y, gcamn->z, &volume) ;
        if (volume > 0)
        {
          m_L = gcamComputeOptimalLinearTransformInRegion(gcam, mri_density, m_L, x, y, z, whalf);
          det = MatrixDeterminant(m_L) ;
          volume *= det ;
        }
        MRIsetVoxVal(mri_out, x, y, z, 0, volume) ;
      }
    }
  }
#endif
  for (x = 0  ; x < mri_out->width ; x++)
  {
    for (y = 0  ; y < mri_out->height ; y++)
    {
      for (z = 0  ; z < mri_out->width ; z++)
      {
        if (!finitep(MRIFvox(mri_out, x, y, z)))
        {
          DiagBreak() ;
        }
      }
    }
  }

  MRIfree(&mri_density) ;
  return(mri_out) ;
}

int
GCAMsample(GCA_MORPH *gcam, float xv, float yv, float zv, float *px, float *py, float *pz)
{
  double            xt, yt, zt ;
#if 0
  int             xi, yi, zi ;
#endif
  if (!gcam->mri_xind)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, "TransformSample: gcam has not been inverted!")) ;

  // the following should not happen /////////////////
  if (xv < 0)
  {
    xv = 0 ;
  }
  if (xv >= gcam->mri_xind->width)
  {
    xv = gcam->mri_xind->width-1 ;
  }
  if (yv < 0)
  {
    yv = 0 ;
  }
  if (yv >= gcam->mri_yind->height)
  {
    yv = gcam->mri_yind->height-1 ;
  }
  if (zv < 0)
  {
    zv = 0 ;
  }
  if (zv >= gcam->mri_zind->depth)
  {
    zv = gcam->mri_zind->depth-1 ;
  }

#if 0
  xi = nint(xv) ;
  yi = nint(yv) ;
  zi = nint(zv) ;
  xt = MRIFvox(gcam->mri_xind, xi, yi, zi)*gcam->spacing ;
  yt = MRIFvox(gcam->mri_yind, xi, yi, zi)*gcam->spacing ;
  zt = MRIFvox(gcam->mri_zind, xi, yi, zi)*gcam->spacing ;
#else
  MRIsampleVolume(gcam->mri_xind, xv, yv, zv, &xt) ;
  MRIsampleVolume(gcam->mri_yind, xv, yv, zv, &yt) ;
  MRIsampleVolume(gcam->mri_zind, xv, yv, zv, &zt) ;
#endif
  *px = xt ;
  *py = yt ;
  *pz = zt ;
  return(NO_ERROR) ;
}
int
GCAMinitLabels(GCA_MORPH *gcam, MRI *mri_labeled)
{
  int    x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (mri_labeled == NULL)
        {
          gcamn->label = gcamn->gc ? 128 : 0 ;
        }
        else
        {
          gcamn->label = MRIgetVoxVal(mri_labeled, x, y, z, 0) ;
        }
        if (gcamn->label > 0)
        {
          DiagBreak() ;
        }
        else
        {
          DiagBreak() ;
        }
      }


  gcam->status = GCAM_LABELED ;
  return(NO_ERROR) ;
}
#define MAX_NODES 1000
#define AINT_NBHD_SIZE 0

double
gcamAreaIntensityEnergy(GCA_MORPH *gcam, MRI *mri, NODE_LOOKUP_TABLE *nlt)
{
  double         sse, val, image_val, sum_v, sum_uv, error, c, distsq, x0, y0, z0, dx, dy, dz, max_increase ;
  int            x, y, z, xn, yn, zn, i, debug = 0, nbrs=0,
                                         x1, y1, z1, xk, yk, zk, xi, yi, zi, xmax, ymax, zmax;
  GCA_MORPH_NODE *gcamn, *gcamn_nbr ;
  NODE_BUCKET    *nb ;

  sse = 0.0 ;
  max_increase = 0.0 ;
  zmax = ymax = xmax = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0 || gcamn->invalid)
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z, 0, SAMPLE_TRILINEAR,&image_val) ;
        xn = nint(gcamn->x) ;
        yn = nint(gcamn->y) ;
        zn = nint(gcamn->z) ;
        x0 = gcamn->x ;
        y0 = gcamn->y ;
        z0 = gcamn->z ;
        nbrs = 0 ;
        sum_uv = sum_v = 0.0 ;
        debug = 0 ;  /* for diagnostics */

        gcamn->sum_ci_vi_ui = gcamn->sum_ci_vi = 0 ;
        for (nbrs = 0, xk = -AINT_NBHD_SIZE ; xk <= AINT_NBHD_SIZE ; xk++)
        {
          xi = xn+xk+NLT_PAD ;
          if (xi < 0 || xi >= nlt->width)
          {
            continue ;
          }
          for (yk = -AINT_NBHD_SIZE ; yk <= AINT_NBHD_SIZE ; yk++)
          {
            yi = yn+yk+NLT_PAD ;
            if (yi < 0 || yi >= nlt->height)
            {
              continue ;
            }
            for (zk = -AINT_NBHD_SIZE ; zk <= AINT_NBHD_SIZE ; zk++)
            {
              zi = zn+zk+NLT_PAD ;
              if (zi < 0 || zi >= nlt->depth)
              {
                continue ;
              }

              nb = &nlt->nodes[xi][yi][zi] ;
              for (i = 0 ; i < nb->nnodes ; i++)
              {
                x1 = nb->node_bins[i].x ;
                y1 = nb->node_bins[i].y ;
                z1 = nb->node_bins[i].z ;   /* gcamn indices */
                gcamn_nbr = &gcam->nodes[x1][y1][z1] ;
                dx = gcamn_nbr->x-xn ;
                dy = gcamn_nbr->y-yn ;
                dz = gcamn_nbr->z-zn ;
                distsq = dx*dx + dy*dy + dz*dz ;

                if (distsq < MIN_NODE_DIST_SQ)  /* less than 1 voxel away */
                {
                  nbrs++ ;

#if 0
                  c = exp(-distsq/NODE_SAMPLE_VAR) ; /* interpolation weight */
#else
                  c = 1.0 ;
#endif
                  val = gcamn_nbr->gc->means[0] ;
                  sum_v += c*gcamn_nbr->area ;
                  sum_uv += c*gcamn_nbr->area*val ;
                  gcamn->sum_ci_vi += gcamn_nbr->area ;
                  gcamn->sum_ci_vi_ui += gcamn_nbr->area*val ;
                }
              }
            }
          }
        }

        if (nbrs == 0)   // out of bounds - shouldn't happen I don't think
        {
          sum_v = gcamn->area ;
          sum_uv = gcamn->area*gcamn->gc->means[0] ;
        }
        if (FZERO(sum_v))
        {
          DiagBreak() ;
          error = image_val ;
          continue ;
        }
        else
        {
          error = image_val - sum_uv/sum_v ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
          printf("E_area_int: node(%d,%d,%d): image %2.1f, node=%2.3f (v=%2.4f)\n",
                 x, y, z, image_val, gcamn->gc->means[0], gcamn->area) ;
          printf("            partial volume intensity = %2.1f (error = %2.1f, nbrs=%d)\n", \
                 sum_uv/sum_v, error, nbrs) ;
        }
        if (!finitep(error))
        {
          DiagBreak() ;
        }
        if (fabs(error) > 1 && (y < 30 && y > 31))
        {
          DiagBreak() ;
        }
        if (error*error-gcamn->last_se > max_increase)
        {
          max_increase = error*error-gcamn->last_se ;
          xmax = x ;
          ymax = y ;
          zmax = z ;
          DiagBreak() ;
        }
        gcamn->last_se = error*error ;
        sse += (error*error) ;
        if (!finitep(sse))
        {
          DiagBreak() ;
        }
      }


  return(sse) ;
}

double
gcamBinaryEnergy(GCA_MORPH *gcam, MRI *mri)
{
  int            x, y, z ;
  double         sse, error ;
  double           val ;
  GCA_MORPH_NODE *gcamn ;

  sse = 0.0 ;
  Galigned = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolumeType(mri, gcamn->x, gcamn->y, gcamn->z,
                            &val, SAMPLE_TRILINEAR) ;
        if ((gcamn->label == 0) || (gcamn->status & GCAM_BINARY_ZERO))
        {
          error = 0-val ;
        }
        else
        {
          error = 1-val ;
        }
        if ((abs(error)<.1) && gcamn->invalid == GCAM_VALID)
        {
          Galigned++ ;
        }
        else
        {
          DiagBreak() ;
          if (gcamn->label > 0)
          {
            DiagBreak() ;
          }
        }

        sse += (error*error) ;
        if (x == Gx && y == Gy && z == Gz)
          printf("E_bnry: node(%d,%d,%d) binary se %2.3f (label=%d, val=%2.3f)\n", \
                 x, y, z, error*error, gcamn->label, val) ;

        if (error*error > gcamn->last_se)
        {
          DiagBreak() ;
        }
        gcamn->last_se = error*error ;
      }

  return(sse*BIN_SCALE) ;
}


#define MAX_NBR_NODES   30000

int
gcamAreaIntensityTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_area_intensity,
                      NODE_LOOKUP_TABLE *nlt, float sigma)
{
  double         val, image_val, error, dx, dy, dz, sum_cv, Ix, Iy, Iz,
                 c,distsq, sum_cuv, dxI, dyI, dzI, predicted_val, norm ;
  int            x, y, z, xn, yn, zn, i, x1, y1, z1, wsize, xi, yi, zi, xk, yk, zk, nbrs ;
  GCA_MORPH_NODE *gcamn, *gcamn_nbr ;
  MRI            *mri_dx, *mri_dy, *mri_dz, *mri_ctrl, *mri_kernel, *mri_nbhd ;
  NODE_BUCKET    *nb ;

  if (DZERO(l_area_intensity))
  {
    return(NO_ERROR) ;
  }

  mri_kernel = MRIgaussian1d(sigma, -1) ;
  wsize = MAX(MIN((2*nint(4*sigma)/2)+1, 7),21) ;
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
  {
    printf("allocating %d x %d x %d nbhd (sigma=%2.1f)...\n",wsize,wsize,wsize,sigma) ;
  }
  mri_nbhd = MRIalloc(wsize, wsize, wsize, MRI_FLOAT) ;

  mri_dx = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
  mri_dy = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
  mri_dz = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
  mri_ctrl = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_UCHAR) ;

  // now form the two terms of the gradient at each node
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0 || gcamn->invalid)
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z, 0, SAMPLE_TRILINEAR,&image_val) ;
        xn = nint(gcamn->x) ;
        yn = nint(gcamn->y) ;
        zn = nint(gcamn->z) ;
        sum_cv = sum_cuv = 0 ;
        for (nbrs = 0, xk = -AINT_NBHD_SIZE ; xk <= AINT_NBHD_SIZE ; xk++)
        {
          xi = xn+xk+NLT_PAD ;
          if (xi < 0 || xi >= nlt->width)
          {
            continue ;
          }
          for (yk = -AINT_NBHD_SIZE ; yk <= AINT_NBHD_SIZE ; yk++)
          {
            yi = yn+yk+NLT_PAD ;
            if (yi < 0 || yi >= nlt->height)
            {
              continue ;
            }
            for (zk = -AINT_NBHD_SIZE ; zk <= AINT_NBHD_SIZE ; zk++)
            {
              zi = zn+zk+NLT_PAD ;
              if (zi < 0 || zi >= nlt->depth)
              {
                continue ;
              }
              nb = &nlt->nodes[xi][yi][zi] ;
              for (i = 0 ; i < nb->nnodes ; i++)
              {
                x1 = nb->node_bins[i].x ;
                y1 = nb->node_bins[i].y ;
                z1 = nb->node_bins[i].z ;   /* gcamn indices */
                gcamn_nbr = &gcam->nodes[x1][y1][z1] ;
                dx = gcamn_nbr->x-xn ;
                dy = gcamn_nbr->y-yn ;
                dz = gcamn_nbr->z-zn ;
                distsq = dx*dx + dy*dy + dz*dz ;

                if (distsq < MIN_NODE_DIST_SQ)  /* less than 1 voxel away */
                {
                  nbrs++ ;

#if 0
                  c = exp(-distsq/NODE_SAMPLE_VAR) ; /* interpolation weight */
#else
                  c = 1.0 ;
#endif
                  val = gcamn_nbr->gc->means[0] ;
                  sum_cv += c*gcamn_nbr->area ;
                  sum_cuv += c*gcamn_nbr->area*val ;
                }
              }
            }
          }
        }

        if (FZERO(sum_cv) || nbrs == 0)   // out of bounds - shouldn't happen I don't think
        {
          sum_cv = gcamn->area ;
          sum_cuv = gcamn->area*gcamn->gc->means[0] ;
          error = image_val ;
          DiagBreak() ;
        }
        gcamn->sum_ci_vi = sum_cv ;
        gcamn->sum_ci_vi_ui = sum_cuv ;

        predicted_val = sum_cuv / sum_cv ;
        gcamn->predicted_val = predicted_val ;
      }

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0 || gcamn->invalid)
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z, 0, SAMPLE_TRILINEAR,&image_val) ;
        xn = nint(gcamn->x) ;
        yn = nint(gcamn->y) ;
        zn = nint(gcamn->z) ;

        /* gradient has two compononents,
          one that pushes the node towards the image location with the predicted intensity
          one that squeezes or expands the node if it is brighter than the predicted val.
        */
        error = image_val - gcamn->predicted_val ;
        gcamComputeMostLikelyDirection(gcam, mri, gcamn->x, gcamn->y, gcamn->z, gcamn->predicted_val,
                                       mri_kernel, mri_nbhd, &Ix, &Iy, &Iz) ;

        norm = sqrt(Ix*Ix+Iy*Iy+Iz*Iz) ;
        if (!FZERO(norm))  /* don't worry about magnitude of gradient */
        {
          Ix /= norm ;
          Iy /= norm ;
          Iz /= norm ;
        }
        dxI = fabs(error)*l_area_intensity * Ix / sqrt(gcamn->gc->covars[0]) ;
        dyI = fabs(error)*l_area_intensity * Iy / sqrt(gcamn->gc->covars[0]) ;
        dzI = fabs(error)*l_area_intensity * Iz / sqrt(gcamn->gc->covars[0]);

        gcamVolumeChangeTermAtNode(gcam, mri, l_area_intensity, x, y, z, &dx, &dy, &dz) ;
#define DISABLE_VOLUME_CHANGE_TERM 0
#if DISABLE_VOLUME_CHANGE_TERM
        dx = dy = dz = 0 ;
#endif
        if (Gx == x && y == Gy && z == Gz)
        {
          printf("            volume    grad=(%2.3f, %2.3f, %2.3f)\n", dx, dy, dz) ;
          printf("            intensity grad=(%2.3f, %2.3f, %2.3f)\n", dxI, dyI, dzI) ;
        }

        if (!finitep(dx) || !finitep(dy) || !finitep(dz) ||
            !finitep(dxI) || !finitep(dyI) || !finitep(dzI))
        {
          DiagBreak() ;
        }
#if !DISABLE_VOLUME_CHANGE_TERM
        MRIFvox(mri_dx, x, y, z) = (dx+dxI) ;
        MRIFvox(mri_dy, x, y, z) = (dy+dyI) ;
        MRIFvox(mri_dz, x, y, z) = (dz+dzI) ;
        MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED ;
        gcamn->dx += (dx+dxI) ;
        gcamn->dy += (dy+dyI) ;
        gcamn->dz += (dz+dzI) ;
#else
        MRIFvox(mri_dx, x, y, z) = (dxI) ;
        MRIFvox(mri_dy, x, y, z) = (dyI) ;
        MRIFvox(mri_dz, x, y, z) = (dzI) ;
        MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED ;
        gcamn->dx += (dxI) ;
        gcamn->dy += (dyI) ;
        gcamn->dz += (dzI) ;
#endif


        if (Gx == x && y == Gy && z == Gz)
        {
          DiagBreak() ;
          printf("l_area_int: node(%d,%d,%d): image %2.1f, uk=%2.1f (v=%2.4f)\n",
                 x, y, z, image_val, gcamn->gc->means[0], gcamn->area) ;
          printf("            partial volume intensity = %2.1f (error = %2.1f), gradI(%2.3f, %2.3f, %2.3f)\n",
                 gcamn->predicted_val, error, dxI, dyI, dzI) ;
          printf("            volume change gradient = (%2.3f, %2.3f, %2.3f)\n",dx, dy, dz) ;
        }
      }


  /* propagate gradients outwards to non-labeled vertices so that they
     move as well. Otherwise they will cause folds at the borders and
     nothing will go anywhere
  */
#if 0
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_dx, "dx.mgz") ;
    MRIwrite(mri_dy, "dy.mgz") ;
    MRIwrite(mri_dz, "dz.mgz") ;
  }
  MRIbuildVoronoiDiagram(mri_dx, mri_ctrl, mri_dx) ;
  MRIbuildVoronoiDiagram(mri_dy, mri_ctrl, mri_dy) ;
  MRIbuildVoronoiDiagram(mri_dz, mri_ctrl, mri_dz) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_dx, "dxv.mgz") ;
    MRIwrite(mri_dy, "dyv.mgz") ;
    MRIwrite(mri_dz, "dzv.mgz") ;
  }
  MRIsoapBubble(mri_dx, mri_ctrl, mri_dx, 10, -1) ;
  MRIsoapBubble(mri_dy, mri_ctrl, mri_dy, 10, -1) ;
  MRIsoapBubble(mri_dz, mri_ctrl, mri_dz, 10, -1) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_dx, "dxs.mgz") ;
    MRIwrite(mri_dy, "dys.mgz") ;
    MRIwrite(mri_dz, "dzs.mgz") ;
  }
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label != 0 || gcamn->invalid)
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        dx = MRIFvox(mri_dx, x, y, z) ;
        dy = MRIFvox(mri_dy, x, y, z) ;
        dz = MRIFvox(mri_dz, x, y, z) ;
        gcamn->dx += dx ;
        gcamn->dy += dy ;
        gcamn->dz += dz ;
      }
#endif

  MRIfree(&mri_dx) ;
  MRIfree(&mri_dy) ;
  MRIfree(&mri_dz) ;
  MRIfree(&mri_ctrl) ;
  MRIfree(&mri_nbhd) ;
  MRIfree(&mri_kernel) ;
  return(NO_ERROR) ;
}


int
gcamBinaryTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, MRI *mri_dist, double l_binary)
{
  int            x, y, z ;
  double         sse, error, dx, dy, dz, odx, ody, odz, norm ;
  double           val ;
  GCA_MORPH_NODE *gcamn ;

  if (DZERO(l_binary))
  {
    return(NO_ERROR) ;
  }

  l_binary *= BIN_SCALE ;

  sse = 0.0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolumeType(mri, gcamn->x, gcamn->y, gcamn->z, &val, SAMPLE_TRILINEAR) ;
#if 0
        MRIsampleVolumeGradientFrame(mri_smooth, gcamn->x, gcamn->y, gcamn->z, &odx, &ody, &odz, 0) ;
#else
        MRIsampleVolumeGradientFrame(mri_dist, gcamn->x, gcamn->y, gcamn->z, &odx, &ody, &odz, 0) ;
        odx *= -1 ;
        ody *= -1 ;
        odz *= -1 ;  // because distance map gradient is neg of what we want
#endif
        norm = sqrt(odx*odx + ody*ody + odz*odz) ;
        if (DZERO(norm))
        {
          continue ;
        }
        odx /= norm ;
        ody /= norm ;
        odz /= norm ;
        if ((gcamn->label == 0) || (gcamn->status & GCAM_BINARY_ZERO))
        {
          error = 0-val ;
        }
        else
        {
          error = 1-val ;
        }

        if (gcamn->label > 0 && gcamn->invalid == GCAM_VALID)
        {
          DiagBreak() ;
        }
        error *= l_binary ;
        if (!FZERO(odx) || !FZERO(ody) || !FZERO(odz))
        {
          DiagBreak() ;
        }


        dx = error * odx ;
        dy = error * ody ;
        dz = error * odz ;
        gcamn->dx += dx ;
        gcamn->dy += dy ;
        gcamn->dz += dz ;
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("node(%d,%d,%d) --> (%2.2f, %2.2f, %2.2f), dI=(%2.3f, %2.3f, %2.3f), label %d, binary %2.3f, grad=(%2.4f,%2.4f,%2.4f)\n",\
                 x, y, z, gcamn->x, gcamn->y, gcamn->z, odx, \
                 ody, odz, gcamn->label, val, dx, dy, dz);
          DiagBreak() ;
        }
      }

  return(sse) ;
}

#if 0
MRI *
GCAMmorphFieldFromAtlas(GCA_MORPH *gcam, MRI *mri, int which, int save_inversions, int filter)
{
  MRI            *mri_morphed, *mri_ctrl, *mri_xind, *mri_yind, *mri_zind,
                 *mri_counts, *mri_mapped ;
  int            x, y, z, xd, yd, zd, nfound = 0, xmin, ymin, zmin,
                                      xmax, ymax, zmax, nmissed, nlabels, inverted = 0, i ;
  GCA_MORPH_NODE *gcamn ;
  VECTOR         *v1, *v2 ;
  MATRIX         *m_lowres_vox2ras, *m_hires_ras2vox, *m_vox2vox,
                 *m_hires_vox2ras ;
  float          cx, cy, cz, rcx, rcy, rcz, rxmin, rxmax, rymin, rymax,
                 rzmin, rzmax ;
  double           xr, yr, zr, num ;
  int            type, nframes = 1, label ;

  switch (which)
  {
  case GCAM_LABEL:
    type = MRI_SHORT /*mri->type*/ ;
    break ;
  case GCAM_DIAG_VOL:
  case GCAM_NODEX:
  case GCAM_NODEY:
  case GCAM_NODEZ:
  case GCAM_ORIGX:
  case GCAM_ORIGY:
  case GCAM_ORIGZ:
    type = MRI_FLOAT ;
    break ;
  case GCAM_MEANS:
    nframes = gcam->ninputs ;
    type = MRI_FLOAT ;
    break ;
  case GCAM_COVARS:
    nframes = (gcam->ninputs * (gcam->ninputs+1)) / 2;
    type = MRI_FLOAT ;
    break ;
  default:
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "GCAMmorphFieldFromAtlas: unsupported field %d",which));
  }
#if 1
  mri_morphed = MRIalloc(gcam->width*1.25, gcam->height*1.25, gcam->depth*3, type) ;
  mri_morphed->x_r = mri->x_r;
  mri_morphed->x_a = mri->x_a;
  mri_morphed->x_s = mri->x_s;
  mri_morphed->y_r = mri->y_r;
  mri_morphed->y_a = mri->y_a;
  mri_morphed->y_s = mri->y_s;
  mri_morphed->z_r = mri->z_r;
  mri_morphed->z_a = mri->z_a;
  mri_morphed->z_s = mri->z_s;
  mri_morphed->xsize = gcam->atlas.xsize ;
  mri_morphed->ysize = gcam->atlas.ysize ;
  mri_morphed->zsize = gcam->atlas.zsize ;
#else
  mri_morphed = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_UCHAR) ;
  mri_morphed->x_r = gcam->atlas.x_r;
  mri_morphed->x_a = gcam->atlas.x_a;
  mri_morphed->x_s = gcam->atlas.x_s;
  mri_morphed->y_r = gcam->atlas.y_r;
  mri_morphed->y_a = gcam->atlas.y_a;
  mri_morphed->y_s = gcam->atlas.y_s;
  mri_morphed->z_r = gcam->atlas.z_r;
  mri_morphed->z_a = gcam->atlas.z_a;
  mri_morphed->z_s = gcam->atlas.z_s;
  mri_morphed->xsize = gcam->atlas.xsize ;
  mri_morphed->ysize = gcam->atlas.ysize ;
  mri_morphed->zsize = gcam->atlas.zsize ;
#endif

  /* compute centroid of labels */
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1,4,1) = 1.0 ;
  *MATRIX_RELT(v2,4,1) = 1.0 ;
  nlabels = 0 ;
  cx = cy = cz = rcx = rcy = rcz = 0 ;
  m_lowres_vox2ras = MRIgetVoxelToRasXform(mri) ;
  m_hires_vox2ras = MRIgetVoxelToRasXform(mri_morphed) ;
  rxmin = rymin = rzmin = 1000000.0 ;
  rxmax = rymax = rzmax = -1000000.0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (x >= mri->width || y >= mri->height || z >= mri->depth)
        {
          continue ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (which == GCAM_DIAG_VOL)
        {
          label = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
        }
        else
        {
          label = gcamn->label ;
        }
        //    if (label > 0)
        {
          V3_X(v1) = gcamn->x ;
          V3_Y(v1) = gcamn->y ;
          V3_Z(v1) = gcamn->z ;
          MatrixMultiply(m_lowres_vox2ras, v1, v2) ;
          xr = V3_X(v2) ;
          yr = V3_Y(v2) ;
          zr = V3_Z(v2) ;
          rcx += xr ;
          rcy += yr ;
          rcz += zr ;
          if (xr < rxmin)
          {
            rxmin = xr ;
          }
          if (yr < rymin)
          {
            rymin = yr ;
          }
          if (zr < rzmin)
          {
            rzmin = zr ;
          }
          if (xr > rxmax)
          {
            rxmax = xr ;
          }
          if (yr > rymax)
          {
            rymax = yr ;
          }
          if (zr > rzmax)
          {
            rzmax = zr ;
          }
          cx += x ;
          cy += y ;
          cz += z ;
          nlabels++ ;
        }
      }
    }
  }
  cx /= nlabels ;
  cy /= nlabels ;
  cz /= nlabels ;
  rcx /= nlabels ;
  rcy /= nlabels ;
  rcz /= nlabels ;

  rcx = (rxmin + rxmax) / 2 ;
  rcy = (rymin + rymax) / 2 ;
  rcz = (rzmin + rzmax) / 2 ;
#if 0
  x = nint(cx) ;
  y = nint(cy) ;
  z = nint(cz) ;
  gcamn = &gcam->nodes[x][y][z] ;
  V3_X(v1) = gcamn->x ;
  V3_Y(v1) = gcamn->y ;
  V3_Z(v1) = gcamn->z ;
  MatrixMultiply(m_lowres_vox2ras, v1, v2) ;
  mri_morphed->c_r = V3_X(v2) ;
  mri_morphed->c_a = V3_Y(v2) ;
  mri_morphed->c_s = V3_Z(v2) ;
#else
  mri_morphed->c_r = rcx ;
  mri_morphed->c_a = rcy ;
  mri_morphed->c_s = rcz ;
#endif
  m_hires_ras2vox = MRIgetRasToVoxelXform(mri_morphed) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("using hires morphed ras2vox xform:\n") ;
    MatrixPrint(stdout, m_hires_ras2vox) ;

    printf("setting c_ras of morphed volume to (%2.1f, %2.1f, %2.1f)\n",\
           mri_morphed->c_r, mri_morphed->c_a,mri_morphed->c_s);
  }
  MatrixMultiply(m_hires_ras2vox, v2, v1) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("maps to (%2.1f, %2.1f, %2.1f)\n", V3_X(v1), V3_Y(v1), V3_Z(v1)) ;
  }

  m_vox2vox = MatrixMultiply(m_hires_ras2vox, m_lowres_vox2ras, NULL) ;


  xmin = ymin = zmin = 1000000;
  xmax = ymax = zmax = -1 ;
  nmissed = 0 ;  /* compute size of volume we will need */
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x >= mri->width || y >= mri->height || z >= mri->depth)
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        /* xform from gcam to lowres volume */
        V3_X(v1) = gcamn->x ;
        V3_Y(v1) = gcamn->y ;
        V3_Z(v1) = gcamn->z ;

        /* xform from lowres volume to hires volume */
        MatrixMultiply(m_vox2vox, v1, v2) ;   /* v2 is now morphed voxel coord */
        xr = V3_X(v2) ;
        yr = V3_Y(v2) ;
        zr = V3_Z(v2) ;
        xd = nint(xr) ;
        yd = nint(yr) ;
        zd = nint(zr) ;

        if (which == GCAM_DIAG_VOL)
        {
          label = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
        }
        else
        {
          label = gcamn->label ;
        }
        //    if (label > 0)
        {
          if (xd < xmin)
          {
            xmin = xd ;
          }
          if (yd < ymin)
          {
            ymin = yd ;
          }
          if (zd < zmin)
          {
            zmin = zd ;
          }
          if (xd > xmax)
          {
            xmax = xd ;
          }
          if (yd > ymax)
          {
            ymax = yd ;
          }
          if (zd > zmax)
          {
            zmax = zd ;
          }
        }
      }
    }
  }

#define PAD 20
  MRIfree(&mri_morphed) ;
  mri_morphed = MRIallocSequence(xmax-xmin+PAD, ymax-ymin+PAD, zmax-zmin+PAD,
                                 type, nframes) ;
  mri_morphed->x_r = mri->x_r;
  mri_morphed->x_a = mri->x_a;
  mri_morphed->x_s = mri->x_s;
  mri_morphed->y_r = mri->y_r;
  mri_morphed->y_a = mri->y_a;
  mri_morphed->y_s = mri->y_s;
  mri_morphed->z_r = mri->z_r;
  mri_morphed->z_a = mri->z_a;
  mri_morphed->z_s = mri->z_s;
  mri_morphed->xsize = gcam->atlas.xsize ;
  mri_morphed->ysize = gcam->atlas.ysize ;
  mri_morphed->zsize = gcam->atlas.zsize ;

  /* map center of bounding box to same RAS point as before*/
  V3_X(v1) = (xmax+xmin)/2 ;
  V3_Y(v1) = (ymax+ymin)/2 ;
  V3_Z(v1) = (zmax+zmin)/2 ;
  MatrixMultiply(m_hires_vox2ras, v1, v2) ;
  mri_morphed->c_r = V3_X(v2) ;
  mri_morphed->c_a = V3_Y(v2) ;
  mri_morphed->c_s = V3_Z(v2) ;

  mri_morphed->c_r = rcx ;
  mri_morphed->c_a = rcy ;
  mri_morphed->c_s = rcz ;
  MatrixFree(&m_hires_vox2ras) ;
  MatrixFree(&m_hires_ras2vox) ;
  m_hires_ras2vox = MRIgetRasToVoxelXform(mri_morphed) ;
  MatrixMultiply(m_hires_ras2vox, m_lowres_vox2ras, m_vox2vox) ;

  /* to get a point in the morphed hires space, use gcam to map from hires->lowres,
     then to RAS using lowres_vox2ras, then the morphed hires_ras2vox */
  mri_mapped = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_UCHAR) ;
  mri_counts = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_morphed, mri_mapped) ;
  MRIcopyHeader(mri_morphed, mri_counts) ;
  if (save_inversions)  /* use previously computed ones (if they exist) */
  {
    mri_xind = gcam->mri_xind ;
    mri_yind = gcam->mri_yind ;
    mri_zind = gcam->mri_zind ;
  }
  else
  {
    mri_xind = mri_yind = mri_zind = NULL ;
  }
  if (mri_xind == NULL)
  {
    mri_xind = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
    mri_yind = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
    mri_zind = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_morphed, mri_xind) ;
    MRIcopyHeader(mri_morphed, mri_yind) ;
    MRIcopyHeader(mri_morphed, mri_zind) ;
    mri_ctrl = MRIclone(mri_mapped, NULL) ;
    MRIcopyHeader(mri_morphed, mri_ctrl) ;
    inverted = 0 ;
    if (save_inversions)  /* use previously computed ones (if they exist) */
    {
      gcam->mri_xind = mri_xind ;
      gcam->mri_yind = mri_yind ;
      gcam->mri_zind = mri_zind ;
    }
  }
  else
  {
    inverted = 1 ;
  }
  xmin = ymin = zmin = 1000000;
  xmax = ymax = zmax = -1 ;
  nmissed = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        /* xform from gcam to target volume */
        V3_X(v1) = gcamn->x ;
        V3_Y(v1) = gcamn->y ;
        V3_Z(v1) = gcamn->z ;

        /* xform from lowres volume to hires volume */
        MatrixMultiply(m_vox2vox, v1, v2) ;   /* v2 is now morphed voxel coord */
        xr = V3_X(v2) ;
        yr = V3_Y(v2) ;
        zr = V3_Z(v2) ;
        xd = nint(xr) ;
        yd = nint(yr) ;
        zd = nint(zr) ;

        if (inverted == 0)
        {
          MRIinterpolateIntoVolume(mri_counts, xr, yr, zr, (double)1.0) ;
          MRIinterpolateIntoVolume(mri_xind, xr, yr, zr, (double)x) ;
          MRIinterpolateIntoVolume(mri_yind, xr, yr, zr, (double)y) ;
          MRIinterpolateIntoVolume(mri_zind, xr, yr, zr, (double)z) ;
        }
        if (xd == Gx && yd == Gy && zd == Gz)
        {
          DiagBreak() ;
        }

        if (which == GCAM_DIAG_VOL)
        {
          if (x >= mri->width || y >= mri->height || z >= mri->depth)
          {
            continue ;
          }
          label = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
        }
        else
        {
          label = gcamn->label ;
        }
        //    if (label > 0)
        {
          if (xd < xmin)
          {
            xmin = xd ;
          }
          if (yd < ymin)
          {
            ymin = yd ;
          }
          if (zd < zmin)
          {
            zmin = zd ;
          }
          if (xd > xmax)
          {
            xmax = xd ;
          }
          if (yd > ymax)
          {
            ymax = yd ;
          }
          if (zd > zmax)
          {
            zmax = zd ;
          }
#if 0
          if (nfound < 5 || (abs(nfound -100000)<5))
            printf("gcamn(%d, %d, %d), label %d --> (%d, %d, %d)\n",
                   x, y, z, gcamn->label, xd, yd, zd) ;
#endif
          nfound++ ;
        }
        if (MRIindexNotInVolume(mri_morphed, xd, yd, zd) == 0)
        {
          if (gcamn->label > 0)
          {
            DiagBreak() ;
          }
          if (MRIgetVoxVal(mri_morphed, xd, yd, zd,0) > 0)
          {
            if (x == Gx && y == Gz && z == Gz)
            {
              DiagBreak() ;
            }
          }
          switch (which)
          {
          case GCAM_DIAG_VOL:
            MRIsetVoxVal(mri_morphed, xd, yd,zd, 0, MRIgetVoxVal(mri,x,y,z,0)) ;
            break ;
          case GCAM_LABEL:
            MRIsetVoxVal(mri_morphed, xd, yd,zd, 0, gcamn->label) ;
            break ;
          case GCAM_NODEX:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->xn ;
            break ;
          case GCAM_NODEY:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->yn ;
            break ;
          case GCAM_NODEZ:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->zn ;
            break ;
          case GCAM_ORIGX:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->origx ;
            break ;
          case GCAM_ORIGY:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->origy ;
            break ;
          case GCAM_ORIGZ:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->origz ;
            break ;
          case GCAM_COVARS:
            if (gcamn->gc != NULL)
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, xd, yd,zd,i) = gcamn->gc->covars[i] ;
              }
            break ;
          case GCAM_MEANS:
            if (gcamn->gc != NULL)
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, xd, yd,zd,i) = gcamn->gc->means[i] ;
              }
            break ;
          }
          if (xd == Gx && yd == Gy && zd == Gz)
          {
            DiagBreak() ;
          }
          if (which == GCAM_MEANS || which == GCAM_COVARS)
          {
            if (gcamn->label > 0)
            {
              MRIvox(mri_mapped, xd, yd, zd) = 1 ;
            }
          }
          else
          {
            MRIvox(mri_mapped, xd, yd, zd) = 1 ;
          }
        }
        else if (gcamn->label > 0)
        {
          nmissed++ ;
          DiagBreak() ;
        }
      }
    }
  }

  if (inverted == 0)
  {
    for (z = 0 ; z < mri_morphed->depth ; z++)
    {
      for (y = 0 ; y < mri_morphed->height ; y++)
      {
        for (x = 0 ; x < mri_morphed->width ; x++)
        {
          //     if (x >= mri->width || y >= mri->height || z >= mri->depth)
          //      continue ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          // get count
          num = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
          if (num == 0)
          {
            continue ;  /* nothing there */
          }
          // give average gcam position for this points
          MRIFvox(mri_xind, x, y, z) = MRIFvox(mri_xind, x, y, z)/(float)num ;
          MRIFvox(mri_yind, x, y, z) = MRIFvox(mri_yind, x, y, z)/(float)num ;
          MRIFvox(mri_zind, x, y, z) = MRIFvox(mri_zind, x, y, z)/(float)num ;
          MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED ;
        }
      }
    }

    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_xind, "xi.mgz") ;
      MRIwrite(mri_yind, "yi.mgz") ;
      MRIwrite(mri_zind, "zi.mgz") ;
    }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      printf("performing soap bubble of x indices...\n") ;
    }
    MRIbuildVoronoiDiagram(mri_xind, mri_ctrl, mri_xind) ;
    MRIsoapBubble(mri_xind, mri_ctrl, mri_xind, 5, -1) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      printf("performing soap bubble of y indices...\n") ;
    }
    MRIbuildVoronoiDiagram(mri_yind, mri_ctrl, mri_yind) ;
    MRIsoapBubble(mri_yind, mri_ctrl, mri_yind, 5, -1) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      printf("performing soap bubble of z indices...\n") ;
    }
    MRIbuildVoronoiDiagram(mri_zind, mri_ctrl, mri_zind) ;
    MRIsoapBubble(mri_zind, mri_ctrl, mri_zind, 5, -1) ;
    MRIfree(&mri_ctrl) ;

    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_xind, "xis.mgz") ;
      MRIwrite(mri_yind, "yis.mgz") ;
      MRIwrite(mri_zind, "zis.mgz") ;
    }
  }

  /* now go through morphed volume and sample back into gcam */
  for (x = 0 ; x < mri_morphed->width ; x++)
  {
    for (y = 0 ; y < mri_morphed->height ; y++)
    {
      for (z = 0 ; z < mri_morphed->depth ; z++)
      {
        //    if (x >= mri->width || y >= mri->height || z >= mri->depth)
        //     continue ;
        if (MRIvox(mri_mapped, x, y, z) > 0)  /* already mapped */
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        xd = nint(MRIgetVoxVal(mri_xind, x, y, z, 0)) ;
        yd = nint(MRIgetVoxVal(mri_yind, x, y, z, 0)) ;
        zd = nint(MRIgetVoxVal(mri_zind, x, y, z, 0)) ;
        if (xd >= 0 && xd < gcam->width &&
            yd >= 0 && yd < gcam->height &&
            zd >= 0 && zd < gcam->depth)
        {
          gcamn = &gcam->nodes[xd][yd][zd] ;
          switch (which)
          {
          case GCAM_DIAG_VOL:
            MRIsetVoxVal(mri_morphed, x, y,z, 0, MRIgetVoxVal(mri,xd,yd,zd,0)) ;
            break ;
          case GCAM_LABEL:
            MRIsetVoxVal(mri_morphed, x, y,z, 0, gcamn->label) ;
            break ;
          case GCAM_NODEX:
            MRIFvox(mri_morphed, x, y,z) = gcamn->xn ;
            break ;
          case GCAM_NODEY:
            MRIFvox(mri_morphed, x, y,z) = gcamn->yn ;
            break ;
          case GCAM_NODEZ:
            MRIFvox(mri_morphed, x, y,z) = gcamn->zn ;
            break ;
          case GCAM_ORIGX:
            MRIFvox(mri_morphed, x, y,z) = gcamn->origx ;
            break ;
          case GCAM_ORIGY:
            MRIFvox(mri_morphed, x, y,z) = gcamn->origy ;
            break ;
          case GCAM_ORIGZ:
            MRIFvox(mri_morphed, x, y,z) = gcamn->origz ;
            break ;
          case GCAM_MEANS:
            if (gcamn->gc != NULL)
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, x, y,z,i) = gcamn->gc->means[i] ;
              }
            break ;
          case GCAM_COVARS:
            if (gcamn->gc != NULL)
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, x, y,z,i) = gcamn->gc->covars[i] ;
              }
            break ;
          }
        }
        else
        {
          MRIsetVoxVal(mri_morphed, x, y, z, 0,0) ;
        }
      }
    }
  }

  if ((Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) || nmissed > 0)
  {
    printf("%d of %d labels mapped outside volume (%2.1f%%), c=(%2.1f,%2.1f,%2.1f)\n",
           nmissed, nlabels, 100.0*nmissed/nlabels, cx, cy, cz) ;

    printf("bounding box = (%d, %d, %d) --> (%d, %d, %d)\n",
           xmin, ymin, zmin, xmax, ymax, zmax) ;
    MRIwrite(mri_morphed, "m.mgz") ;
  }
  switch (which)
  {
  case GCAM_LABEL:
  {
    MRI *mri_filtered ;
    if (filter > 0)
    {
      mri_filtered = MRImodeFilter(mri_morphed, NULL, filter) ;
      MRIfree(&mri_morphed) ;
      mri_morphed = mri_filtered ;
    }
    break ;
  }
  default:
    break ;
  }

  MRIreInitCache(mri_morphed); // recalculate the stored transforms
  MatrixFree(&m_lowres_vox2ras) ;
  MatrixFree(&m_hires_ras2vox);
  MatrixFree(&m_vox2vox) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  MRIfree(&mri_counts) ;
  if (save_inversions == 0)
  {
    MRIfree(&mri_xind) ;
    MRIfree(&mri_yind) ;
    MRIfree(&mri_zind) ;
  }
  MRIfree(&mri_mapped) ;
  MRIremoveNaNs(mri_morphed, mri_morphed) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    static int mno = 0 ;
    char fname[STRLEN] ;
    sprintf(fname, "m%d.mgz", mno++) ;
    MRIwrite(mri_morphed, fname) ;
  }
  return(mri_morphed) ;
}
#else
MRI *
GCAMmorphFieldFromAtlas(GCA_MORPH *gcam, MRI *mri, int which, int save_inversions, int filter)
{
  MRI            *mri_morphed, *mri_ctrl, *mri_xind, *mri_yind, *mri_zind,
                 *mri_counts, *mri_mapped, *mri_tmp, *mri_tmp2 ;
  int            x, y, z, xd, yd, zd, nfound = 0, xmin, ymin, zmin,
                                      xmax, ymax, zmax, nmissed, nlabels, inverted = 0, i ;
  GCA_MORPH_NODE *gcamn ;
  VECTOR         *v1, *v2 ;
  MATRIX         *m_target_vox2ras, *m_morphed_ras2vox, *m_vox2vox ;
  double           xr, yr, zr, num, val ;
  int            type, nframes = 1, label, pad ;
  MRI_REGION     box ;

  gcamn = NULL;
  switch (which)
  {
  case GCAM_LABEL:
    type = MRI_SHORT /*mri->type*/ ;
    break ;
  case GCAM_DIAG_VOL:
  case GCAM_NODEX:
  case GCAM_NODEY:
  case GCAM_NODEZ:
  case GCAM_ORIGX:
  case GCAM_ORIGY:
  case GCAM_ORIGZ:
  case GCAM_X_GRAD:
  case GCAM_Y_GRAD:
  case GCAM_Z_GRAD:
  case GCAM_JACOBIAN:
  case GCAM_AREA:
  case GCAM_ORIG_AREA:
    type = MRI_FLOAT ;
    break ;
  case GCAM_MEANS:
    nframes = gcam->ninputs ;
    type = MRI_FLOAT ;
    break ;
  case GCAM_COVARS:
    nframes = (gcam->ninputs * (gcam->ninputs+1)) / 2;
    type = MRI_FLOAT ;
    break ;
  default:
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "GCAMmorphFieldFromAtlas: unsupported field %d",which));
  }

  // allocate an initial morphed volume and see if morphed image will fit
  mri_tmp = MRIalloc(gcam->width, gcam->height, gcam->depth, type) ;
#if 0
  mri_tmp->x_r = mri->x_r;
  mri_tmp->x_a = mri->x_a;
  mri_tmp->x_s = mri->x_s;
  mri_tmp->y_r = mri->y_r;
  mri_tmp->y_a = mri->y_a;
  mri_tmp->y_s = mri->y_s;
  mri_tmp->z_r = mri->z_r;
  mri_tmp->z_a = mri->z_a;
  mri_tmp->z_s = mri->z_s;
  mri_tmp->c_r = mri->c_r;
  mri_tmp->c_a = mri->c_a;
  mri_tmp->c_s = mri->c_s ;
#else
  mri_tmp->x_r = gcam->atlas.x_r;
  mri_tmp->x_a = gcam->atlas.x_a;
  mri_tmp->x_s = gcam->atlas.x_s;
  mri_tmp->y_r = gcam->atlas.y_r;
  mri_tmp->y_a = gcam->atlas.y_a;
  mri_tmp->y_s = gcam->atlas.y_s;
  mri_tmp->z_r = gcam->atlas.z_r;
  mri_tmp->z_a = gcam->atlas.z_a;
  mri_tmp->z_s = gcam->atlas.z_s;
  mri_tmp->c_r = gcam->atlas.c_r;
  mri_tmp->c_a = gcam->atlas.c_a;
  mri_tmp->c_s = gcam->atlas.c_s ;
#endif
  mri_tmp->xsize = gcam->atlas.xsize ;
  mri_tmp->ysize = gcam->atlas.ysize ;
  mri_tmp->zsize = gcam->atlas.zsize ;
  MRIreInitCache(mri_tmp) ;

  // compute matrix that takes target voxels and maps them to morphed ones
  m_target_vox2ras = MRIgetVoxelToRasXform(mri) ;
  m_morphed_ras2vox = MRIgetRasToVoxelXform(mri_tmp) ;
  m_vox2vox = MatrixMultiply(m_morphed_ras2vox, m_target_vox2ras, NULL) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1,4,1) = 1.0 ;
  *MATRIX_RELT(v2,4,1) = 1.0 ;

  // go through volume and see how much bigger it has to be to fit
  pad = 0 ;
  xmin = mri_tmp->width;
  ymin = mri_tmp->height ;
  zmin = mri_tmp->depth ;
  xmax = ymax = zmax = 0 ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (which == GCAM_DIAG_VOL)
        {
          if (x >= mri->width || y >= mri->height || z >= mri->depth)
          {
            continue ;
          }
          label = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
        }
        else
        {
          label = gcamn->label ;
        }
        if (label == 0 && which ==GCAM_LABEL)
        {
          continue ;
        }
        V3_X(v1) = gcamn->x ;
        V3_Y(v1) = gcamn->y ;
        V3_Z(v1) = gcamn->z ;

        /* xform from target volume to morphed volume */
        MatrixMultiply(m_vox2vox, v1, v2) ;   /* v2 is now morphed voxel coord */
        xr = V3_X(v2) ;
        yr = V3_Y(v2) ;
        zr = V3_Z(v2) ;
        xd = nint(xr) ;
        yd = nint(yr) ;
        zd = nint(zr) ;

        if (x == Gx && y == Gy && z == Gz)
        {
          printf("(xd, yd, zd) = (%d, %d, %d)\n", xd, yd, zd) ;
        }

        if (xd < 0 || yd < 0 || zd < 0)
        {
          DiagBreak() ;
        }
        xmin = MIN(xmin, xd) ;
        ymin = MIN(ymin, yd) ;
        zmin = MIN(zmin, zd) ;
        xmax = MAX(xmax, xd) ;
        ymax = MAX(ymax, yd) ;
        zmax = MAX(zmax, zd) ;
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("(xmin, ymin, zmin) = (%d, %d, %d); (xmax, ymax, zmax) = (%d, %d, %d);\n", xmin, ymin, zmin, xmax, ymax, zmax) ;
        }
        if (xd >= 0 && yd >= 0 && zd >= 0 &&
            xd < mri_tmp->width && yd < mri_tmp->height && zd < mri_tmp->depth)
        {
          MRIsetVoxVal(mri_tmp, xd, yd, zd, 0, label) ;  // for debugging
        }
      }
    }
  }

  //if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  MRIwrite(mri_tmp, "tmp1.mgz") ;

  pad = MAX(-xmin, MAX(-ymin, -zmin)) ;
  pad = MAX(pad, (xmax-mri_tmp->width)+1) ;
  pad = MAX(pad, (ymax-mri_tmp->height)+1) ;
  pad = MAX(pad, (zmax-mri_tmp->depth)+1) ;
  pad = MAX(1,pad) ;
  printf("padding by %d voxels\n", pad) ;
  mri_tmp2 = MRIextractRegionAndPad(mri_tmp, NULL, NULL, pad) ;
  printf("Left MRIextractRegionAndPad\n") ;
  MRIfree(&mri_tmp) ;
  MatrixFree(&m_morphed_ras2vox) ;
  m_morphed_ras2vox = MRIgetRasToVoxelXform(mri_tmp2) ;
  MatrixMultiply(m_morphed_ras2vox, m_target_vox2ras, m_vox2vox) ;

  // go through volume and see how much we can crop it
  xmin = mri_tmp2->width;
  ymin = mri_tmp2->height ;
  zmin = mri_tmp2->depth ;
  xmax = ymax = zmax = 0 ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (which == GCAM_DIAG_VOL)
        {
          if (x >= mri->width || y >= mri->height || z >= mri->depth)
          {
            continue ;
          }
          label = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
        }
        else
        {
          label = gcamn->label ;
        }
        if (label == 0)
          //if(0)
        {
          continue ;
        }
        V3_X(v1) = gcamn->x ;
        V3_Y(v1) = gcamn->y ;
        V3_Z(v1) = gcamn->z ;

        /* xform from target volume to morphed volume */
        MatrixMultiply(m_vox2vox, v1, v2) ;   /* v2 is now morphed voxel coord */
        xr = V3_X(v2) ;
        yr = V3_Y(v2) ;
        zr = V3_Z(v2) ;
        xd = nint(xr) ;
        yd = nint(yr) ;
        zd = nint(zr) ;

        xmin = MIN(xmin, xd) ;
        ymin = MIN(ymin, yd) ;
        zmin = MIN(zmin, zd) ;
        xmax = MAX(xmax, xd) ;
        ymax = MAX(ymax, yd) ;
        zmax = MAX(zmax, zd) ;

        if (xmin < -100 || ymin < -100 || zmin < -100)
        {
          DiagBreak() ;
        }
        if (xd >= 0 && yd >= 0 && zd >= 0 &&
            xd < mri_tmp2->width && yd < mri_tmp2->height && zd < mri_tmp2->depth)
        {
          MRIsetVoxVal(mri_tmp2, xd, yd, zd, 0, label) ;  // for debugging
        }
      }
    }
  }

  //if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  MRIwrite(mri_tmp2, "tmp2.mgz") ;
  if (xmin < -100)
  {
    DiagBreak() ;
  }
  printf("bounding box: (%d, %d, %d) --> (%d, %d, %d)\n", xmin, ymin, zmin, xmax, ymax, zmax) ;
  box.x = xmin ;
  box.dx = xmax-xmin+1 ;
  box.y = ymin ;
  box.dy = ymax-ymin+1 ;
  box.z = zmin ;
  box.dz = zmax-zmin+1 ;
  mri_morphed = MRIextractRegion(mri_tmp2, NULL, &box) ;
  MRIfree(&mri_tmp2) ;

  MatrixFree(&m_morphed_ras2vox) ;
  m_morphed_ras2vox = MRIgetRasToVoxelXform(mri_morphed) ;
  MatrixMultiply(m_morphed_ras2vox, m_target_vox2ras, m_vox2vox) ;

  nlabels = 0 ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("using morphed morphed ras2vox xform:\n") ;
    MatrixPrint(stdout, m_morphed_ras2vox) ;

    printf("setting c_ras of morphed volume to (%2.1f, %2.1f, %2.1f)\n",\
           mri_morphed->c_r, mri_morphed->c_a,mri_morphed->c_s);
  }

  m_vox2vox = MatrixMultiply(m_morphed_ras2vox, m_target_vox2ras, NULL) ;

  /* to get a point in the morphed morphed space, use gcam to map from morphed->target,
  then to RAS using target_vox2ras, then the morphed morphed_ras2vox */
  mri_mapped = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_UCHAR) ;
  mri_counts = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_morphed, mri_mapped) ;
  MRIcopyHeader(mri_morphed, mri_counts) ;
  /* use previously computed ones (if they exist) */
  mri_xind = gcam->mri_xind ;
  mri_yind = gcam->mri_yind ;
  mri_zind = gcam->mri_zind ;
  if (mri_xind == NULL)
  {
    mri_xind = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
    mri_yind = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
    mri_zind = MRIalloc(mri_morphed->width, mri_morphed->height, mri_morphed->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_morphed, mri_xind) ;
    MRIcopyHeader(mri_morphed, mri_yind) ;
    MRIcopyHeader(mri_morphed, mri_zind) ;
    mri_ctrl = MRIclone(mri_mapped, NULL) ;
    MRIcopyHeader(mri_morphed, mri_ctrl) ;
    inverted = 0 ;
    gcam->mri_xind = mri_xind ;
    gcam->mri_yind = mri_yind ;
    gcam->mri_zind = mri_zind ;
  }
  else
  {
    inverted = 1 ;
  }
  xmin = ymin = zmin = 1000000;
  xmax = ymax = zmax = -1 ;
  nmissed = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        /* xform from gcam to target volume */
        V3_X(v1) = gcamn->x ;
        V3_Y(v1) = gcamn->y ;
        V3_Z(v1) = gcamn->z ;

        /* xform from target volume to morphed volume */
        MatrixMultiply(m_vox2vox, v1, v2) ;   /* v2 is now morphed voxel coord */
        xr = V3_X(v2) ;
        yr = V3_Y(v2) ;
        zr = V3_Z(v2) ;
        xd = nint(xr) ;
        yd = nint(yr) ;
        zd = nint(zr) ;

        if (inverted == 0)
        {
          MRIinterpolateIntoVolume(mri_counts, xr, yr, zr, (double)1.0) ;
          MRIinterpolateIntoVolume(mri_xind, xr, yr, zr, (double)x) ;
          MRIinterpolateIntoVolume(mri_yind, xr, yr, zr, (double)y) ;
          MRIinterpolateIntoVolume(mri_zind, xr, yr, zr, (double)z) ;
        }
        if (xd == Gx && yd == Gy && zd == Gz)
        {
          DiagBreak() ;
        }

        if (which == GCAM_DIAG_VOL)
        {
          if (x >= mri->width || y >= mri->height || z >= mri->depth)
          {
            continue ;
          }
          label = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
        }
        else
        {
          label = gcamn->label ;
        }
        //    if (label > 0)
        {
          if (xd < xmin)
          {
            xmin = xd ;
          }
          if (yd < ymin)
          {
            ymin = yd ;
          }
          if (zd < zmin)
          {
            zmin = zd ;
          }
          if (xd > xmax)
          {
            xmax = xd ;
          }
          if (yd > ymax)
          {
            ymax = yd ;
          }
          if (zd > zmax)
          {
            zmax = zd ;
          }
#if 0
          if (nfound < 5 || (abs(nfound -100000)<5))
            printf("gcamn(%d, %d, %d), label %d --> (%d, %d, %d)\n",
                   x, y, z, gcamn->label, xd, yd, zd) ;
#endif
          nfound++ ;
        }
        if (MRIindexNotInVolume(mri_morphed, xd, yd, zd) == 0)
        {
          if (gcamn->label > 0)
          {
            DiagBreak() ;
          }
          if (MRIgetVoxVal(mri_morphed, xd, yd, zd,0) > 0)
          {
            if (x == Gx && y == Gz && z == Gz)
            {
              DiagBreak() ;
            }
          }
          switch (which)
          {
          case GCAM_DIAG_VOL:
            MRIsetVoxVal(mri_morphed, xd, yd,zd, 0, MRIgetVoxVal(mri,x,y,z,0)) ;
            break ;
          case GCAM_LABEL:
            MRIsetVoxVal(mri_morphed, xd, yd,zd, 0, gcamn->label) ;
            break ;
          case GCAM_NODEX:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->xn ;
            break ;
          case GCAM_NODEY:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->yn ;
            break ;
          case GCAM_NODEZ:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->zn ;
            break ;
          case GCAM_ORIGX:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->origx ;
            break ;
          case GCAM_ORIGY:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->origy ;
            break ;
          case GCAM_ORIGZ:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->origz ;
            break ;
          case GCAM_X_GRAD:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->dx ;
            break ;
          case GCAM_Y_GRAD:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->dy ;
            break ;
          case GCAM_Z_GRAD:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->dz ;
            break ;
          case GCAM_JACOBIAN:
            if (!FZERO(gcamn->orig_area))
            {
              MRIFvox(mri_morphed, xd, yd,zd) = gcamn->area / gcamn->orig_area ;
            }
            else
            {
              MRIFvox(mri_morphed, xd, yd,zd) = 1.0 ;
            }
            break ;
          case GCAM_AREA:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->area ;
            break ;
          case GCAM_ORIG_AREA:
            MRIFvox(mri_morphed, xd, yd,zd) = gcamn->orig_area ;
            break ;
          case GCAM_COVARS:
            if (gcamn->gc != NULL)
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, xd, yd,zd,i) = gcamn->gc->covars[i] ;
              }
            break ;
          case GCAM_MEANS:
            if (gcamn->gc != NULL)
            {
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, xd, yd,zd,i) = gcamn->gc->means[i] ;
              }
            }
            else
            {
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, xd, yd,zd,i) = 0 ;
              }
            }
            break ;
          }
          if (xd == Gx && yd == Gy && zd == Gz)
          {
            DiagBreak() ;
          }
          if (which == GCAM_MEANS || which == GCAM_COVARS)
          {
            //      if (gcamn->label > 0)
            MRIvox(mri_mapped, xd, yd, zd) = 1 ;
          }
          else
          {
            MRIvox(mri_mapped, xd, yd, zd) = 1 ;
          }
        }
        else if (gcamn->label > 0)
        {
          nmissed++ ;
          DiagBreak() ;
        }
      }
    }
  }

  if (inverted == 0)
  {
    for (z = 0 ; z < mri_morphed->depth ; z++)
    {
      for (y = 0 ; y < mri_morphed->height ; y++)
      {
        for (x = 0 ; x < mri_morphed->width ; x++)
        {
          //     if (x >= mri->width || y >= mri->height || z >= mri->depth)
          //      continue ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          // get count
          num = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
          if (num == 0)
          {
            continue ;  /* nothing there */
          }
          // give average gcam position for this points
          MRIFvox(mri_xind, x, y, z) = MRIFvox(mri_xind, x, y, z)/(float)num ;
          MRIFvox(mri_yind, x, y, z) = MRIFvox(mri_yind, x, y, z)/(float)num ;
          MRIFvox(mri_zind, x, y, z) = MRIFvox(mri_zind, x, y, z)/(float)num ;
          MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED ;
        }
      }
    }

    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_xind, "xi.mgz") ;
      MRIwrite(mri_yind, "yi.mgz") ;
      MRIwrite(mri_zind, "zi.mgz") ;
    }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      printf("performing soap bubble of x indices...\n") ;
    }
    MRIbuildVoronoiDiagram(mri_xind, mri_ctrl, mri_xind) ;
    MRIsoapBubble(mri_xind, mri_ctrl, mri_xind, 5, -1) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      printf("performing soap bubble of y indices...\n") ;
    }
    MRIbuildVoronoiDiagram(mri_yind, mri_ctrl, mri_yind) ;
    MRIsoapBubble(mri_yind, mri_ctrl, mri_yind, 5, -1) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      printf("performing soap bubble of z indices...\n") ;
    }
    MRIbuildVoronoiDiagram(mri_zind, mri_ctrl, mri_zind) ;
    MRIsoapBubble(mri_zind, mri_ctrl, mri_zind, 5, -1) ;
    MRIfree(&mri_ctrl) ;

    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_xind, "xis.mgz") ;
      MRIwrite(mri_yind, "yis.mgz") ;
      MRIwrite(mri_zind, "zis.mgz") ;
    }
  }

  /* now go through morphed volume and sample back into gcam */
  if (which == GCAM_ORIGX || which == GCAM_ORIGY || which == GCAM_ORIGZ)
  {
    mri_tmp = GCAMwriteMRI(gcam, NULL, which) ;
  }
  else
  {
    mri_tmp = NULL ;
  }
  for (x = 0 ; x < mri_morphed->width ; x++)
  {
    for (y = 0 ; y < mri_morphed->height ; y++)
    {
      for (z = 0 ; z < mri_morphed->depth ; z++)
      {
        //    if (x >= mri->width || y >= mri->height || z >= mri->depth)
        //     continue ;
        if (MRIvox(mri_mapped, x, y, z) > 0)  /* already mapped */
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        xd = MRIgetVoxVal(mri_xind, x, y, z, 0) ;
        yd = MRIgetVoxVal(mri_yind, x, y, z, 0) ;
        zd = MRIgetVoxVal(mri_zind, x, y, z, 0) ;
        if (mri_tmp == NULL)   // not sampling from volume - use nearest nbr in gcam
        {
          xd = nint(xd) ;
          yd = nint(yd) ;
          zd = nint(zd) ;
        }
        if (xd >= 0 && xd < gcam->width &&
            yd >= 0 && yd < gcam->height &&
            zd >= 0 && zd < gcam->depth)
        {
          if (mri_tmp == NULL)
          {
            gcamn = &gcam->nodes[nint(xd)][nint(yd)][nint(zd)] ;
          }
          else
          {
            MRIsampleVolume(mri_tmp, xd, yd, zd, &val) ;
          }
          switch (which)
          {
          case GCAM_DIAG_VOL:
            MRIsetVoxVal(mri_morphed, x, y,z, 0, MRIgetVoxVal(mri,xd,nint(yd),nint(zd),0)) ;
            break ;
          case GCAM_LABEL:
            MRIsetVoxVal(mri_morphed, x, y,z, 0, gcamn->label) ;
            break ;
          case GCAM_NODEX:
            MRIFvox(mri_morphed, x, y,z) = gcamn->xn ;
            break ;
          case GCAM_NODEY:
            MRIFvox(mri_morphed, x, y,z) = gcamn->yn ;
            break ;
          case GCAM_NODEZ:
            MRIFvox(mri_morphed, x, y,z) = gcamn->zn ;
            break ;
          case GCAM_X_GRAD:
            MRIFvox(mri_morphed, x, y,z) = gcamn->dx ;
            break ;
          case GCAM_Y_GRAD:
            MRIFvox(mri_morphed, x, y,z) = gcamn->dy ;
            break ;
          case GCAM_Z_GRAD:
            MRIFvox(mri_morphed, x, y,z) = gcamn->dz ;
            break ;
          case GCAM_JACOBIAN:
            if (!FZERO(gcamn->orig_area))
            {
              MRIFvox(mri_morphed, x, y,z) = gcamn->area / gcamn->orig_area ;
            }
            else
            {
              MRIFvox(mri_morphed, x, y,z) = 1.0 ;
            }
            break ;
          case GCAM_AREA:
            MRIFvox(mri_morphed, x, y,z) = gcamn->area ;
            break ;
          case GCAM_ORIG_AREA:
            MRIFvox(mri_morphed, x, y,z) = gcamn->orig_area ;
            break ;
          case GCAM_ORIGX:
            MRIFvox(mri_morphed, x, y,z) = val ;
            break ;
          case GCAM_ORIGY:
            MRIFvox(mri_morphed, x, y,z) = val ;
            break ;
          case GCAM_ORIGZ:
            MRIFvox(mri_morphed, x, y,z) = val ;
            break ;
          case GCAM_MEANS:
            if (gcamn->gc != NULL)
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, x, y,z,i) = gcamn->gc->means[i] ;
              }
            break ;
          case GCAM_COVARS:
            if (gcamn->gc != NULL)
              for (i = 0 ; i < mri_morphed->nframes ; i++)
              {
                MRIFseq_vox(mri_morphed, x, y,z,i) = gcamn->gc->covars[i] ;
              }
            break ;
          }
        }
        else
        {
          MRIsetVoxVal(mri_morphed, x, y, z, 0,0) ;
        }
      }
    }
  }

  if ((Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) || nmissed > 0)
  {
    printf("%d of %d labels mapped outside volume (%2.1f%%)\n",
           nmissed, nlabels, 100.0*nmissed/nlabels) ;

    printf("bounding box = (%d, %d, %d) --> (%d, %d, %d)\n",
           xmin, ymin, zmin, xmax, ymax, zmax) ;
    MRIwrite(mri_morphed, "m.mgz") ;
  }
  switch (which)
  {
  case GCAM_LABEL:
  {
    MRI *mri_filtered ;
    if (filter > 0)
    {
      mri_filtered = MRImodeFilter(mri_morphed, NULL, filter) ;
      MRIfree(&mri_morphed) ;
      mri_morphed = mri_filtered ;
    }
    break ;
  }
  default:
    break ;
  }

  MRIreInitCache(mri_morphed); // recalculate the stored transforms
  MatrixFree(&m_target_vox2ras) ;
  MatrixFree(&m_morphed_ras2vox);
  MatrixFree(&m_vox2vox) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  MRIfree(&mri_counts) ;
  if (save_inversions == 0)
  {
    MRIfree(&gcam->mri_xind) ;
    MRIfree(&gcam->mri_yind) ;
    MRIfree(&gcam->mri_zind) ;
  }
  MRIfree(&mri_mapped) ;
  MRIremoveNaNs(mri_morphed, mri_morphed) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    static int mno = 0 ;
    char fname[STRLEN] ;
    sprintf(fname, "m%d.mgz", mno++) ;
    MRIwrite(mri_morphed, fname) ;
  }
  return(mri_morphed) ;
}
#endif
int
GCAMexpand(GCA_MORPH *gcam, float distance)
{
  int            i, j, k, num ;
  double         cx, cy, cz, dx, dy, dz, norm ;
  GCA_MORPH_NODE *gcamn ;

  cx = cy = cz = 0 ;
  for (num = i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        if (gcamn->label == 0)
        {
          continue ;
        }
        num++ ;
        cx += gcamn->x ;
        cy += gcamn->y ;
        cz += gcamn->z ;
      }
    }
  }
  cx /= num ;
  cy /= num ;
  cz /= num ;

  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        if (i > 0 && j > 0 && k > 0&&
            i < gcam->width-1 && j < gcam->height-1 && k < gcam->depth-1)
        {
          continue ;  /* not a border node */
        }
        gcamn = &gcam->nodes[i][j][k] ;
        dx = gcamn->x - cx ;
        dy = gcamn->y - cy ;
        dz = gcamn->z - cz ;
        norm = sqrt(dx*dx + dy*dy+dz*dz) ;
        if (FZERO(norm))
        {
          continue ;
        }
        dx /= norm ;
        dy /= norm ;
        dz /= norm ;
        gcamn->x += dx*distance ;
        gcamn->y += dy*distance ;
        gcamn->z += dz*distance ;
        gcamn->invalid |= GCAM_BORDER_NODE ;
      }
    }
  }
  return(NO_ERROR) ;
}
GCA_MORPH *
GCAMregrid(GCA_MORPH *gcam, MRI *mri_dst, int pad, GCA_MORPH_PARMS *parms, MRI **pmri)
{
  MRI        *mri_tmp, *mri_labeled, *mri_means, *mri_covars,
             *mri_origx, *mri_origy, *mri_origz, *mri_tmp2,
             *mri_xn, *mri_yn, *mri_zn ;
  MRI_REGION box ;
  MATRIX     *m_I, *m_vox2vox ;
  TRANSFORM  *transform ;
  GCA_MORPH  *gcam_new ;
  int        x,y,z,type, use_means = 0, ninputs, spacing ;
  GCA_MORPH_NODE *gcamn ;
  float     means[10000], vars[10000];
  double    exp_k ;
  VOL_GEOM  image, atlas ;
  GCA       *gca ;

  memset(means, 0, sizeof(means)) ;
  memset(vars, 0, sizeof(vars)) ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0 || gcamn->gc == NULL)
        {
          continue ;
        }
        use_means = 1 ;
        means[gcamn->label] = gcamn->gc->means[0] ;
        vars[gcamn->label] = gcamn->gc->covars[0] ;
      }
    }
  }

  if (DIAG_VERBOSE_ON)
  {
    write_snapshot(gcam, mri_dst, parms, 100) ;
  }
  printf("regridding GCAM...\n") ;
  if (parms->log_fp)
  {
    fprintf(parms->log_fp, "regridding GCAM...\n") ;
  }
  transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;

  if (gcam->mri_xind)
  {
    MRIfree(&gcam->mri_xind) ;
  }
  if (gcam->mri_yind)
  {
    MRIfree(&gcam->mri_yind) ;
  }
  if (gcam->mri_zind)
  {
    MRIfree(&gcam->mri_zind) ;
  }
  mri_tmp = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_LABEL,1, 0) ;
  if (Gx >= 0 && Gx < gcam->width && Gy < gcam->height && Gz < gcam->depth)
  {
    int            ox = Gx, oy = Gy, oz = Gz;
    VECTOR         *v1, *v2 ;
    GCA_MORPH_NODE *gcamn ;
    MATRIX         *m_lowres_vox2ras, *m_morphed_ras2vox, *m_vox2vox ;

    m_lowres_vox2ras = MRIgetVoxelToRasXform(parms->mri) ;
    m_morphed_ras2vox = MRIgetRasToVoxelXform(mri_tmp) ;
    m_vox2vox = MatrixMultiply(m_morphed_ras2vox, m_lowres_vox2ras, NULL) ;
    gcamn = &gcam->nodes[Gx][Gy][Gz] ;
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v1,4,1) = 1.0 ;
    V3_X(v1) = gcamn->x ;
    V3_Y(v1) = gcamn->y ;
    V3_Z(v1) = gcamn->z ;
    v2 = MatrixMultiply(m_vox2vox, v1, NULL) ;
    gcamn = &gcam->nodes[Gx][Gy][Gz] ;
    Gx = nint(V3_X(v2)) ;
    Gy = nint(V3_Y(v2)) ;
    Gz = nint(V3_Z(v2)) ;
    printf("G[xyz]  (%d, %d, %d) --> (%d, %d, %d)\n", ox,oy,oz, Gx,Gy,Gz) ;

    MatrixFree(&m_lowres_vox2ras) ;
    MatrixFree(&m_morphed_ras2vox) ;
    VectorFree(&v1) ;
    VectorFree(&v2) ;
    MatrixFree(&m_vox2vox) ;
  }
  MRIboundingBox(mri_tmp, 0, &box) ;
  box.x -= pad ;
  box.y -= pad ;
  box.z -= pad ;
  box.dx += 2*pad ;
  box.dy += 2*pad ;
  box.dz += 2*pad ;
  MRIcropBoundingBox(mri_tmp, &box) ;

  if (Gx >= 0)
  {
    int ox = Gx, oy = Gy, oz = Gz;
    Gx -= box.x ;
    Gy -= box.y ;
    Gz -= box.z ;
    printf("G[xyz]  (%d, %d, %d) --> (%d, %d, %d)\n", ox,oy,oz, Gx,Gy,Gz) ;
  }
  mri_labeled = MRIextractRegion(mri_tmp, NULL, &box) ;
  MRIfree(&mri_tmp) ;

  if (use_means)
  {
    mri_covars = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_COVARS,1, 0) ;
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_covars, "c.mgz") ;
    }
    mri_means = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_MEANS,1, 0) ;
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_means, "m.mgz") ;
    }
  }
  else
  {
    mri_means = mri_covars = NULL ;
  }

  mri_origx = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_ORIGX,1, 0) ;
  mri_origy = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_ORIGY,1, 0) ;
  mri_origz = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_ORIGZ,1, 0) ;
  if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_origz, "z.mgz") ;
  }
  mri_xn = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_NODEX,1, 0) ;
  mri_yn = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_NODEY,1, 0) ;
  mri_zn = GCAMmorphFieldFromAtlas(gcam, parms->mri, GCAM_NODEZ,0, 0) ;

  mri_tmp2 = MRIextractRegion(mri_xn, NULL, &box) ;
  MRIfree(&mri_xn) ;
  mri_xn = mri_tmp2 ;
  mri_tmp2 = MRIextractRegion(mri_yn, NULL, &box) ;
  MRIfree(&mri_yn) ;
  mri_yn = mri_tmp2 ;
  mri_tmp2 = MRIextractRegion(mri_zn, NULL, &box) ;
  MRIfree(&mri_zn) ;
  mri_zn = mri_tmp2 ;
  if (use_means)
  {
    mri_tmp2 = MRIextractRegion(mri_means, NULL, &box) ;
    MRIfree(&mri_means) ;
    mri_means = mri_tmp2 ;
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_means, "m2.mgz") ;
    }
    mri_tmp2 = MRIextractRegion(mri_covars, NULL, &box) ;
    MRIfree(&mri_covars) ;
    mri_covars = mri_tmp2 ;
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      MRIwrite(mri_covars, "c2.mgz") ;
    }
  }

  mri_tmp2 = MRIextractRegion(mri_origx, NULL, &box) ;
  MRIfree(&mri_origx) ;
  mri_origx = mri_tmp2 ;
  mri_tmp2 = MRIextractRegion(mri_origy, NULL, &box) ;
  MRIfree(&mri_origy) ;
  mri_origy = mri_tmp2 ;
  mri_tmp2 = MRIextractRegion(mri_origz, NULL, &box) ;
  MRIfree(&mri_origz) ;
  mri_origz = mri_tmp2 ;


  m_I = MatrixIdentity(4, NULL) ;
  m_vox2vox = ((LTA *)(transform->xform))->xforms[0].m_L;
  MRIrasXformToVoxelXform(mri_labeled, mri_dst, m_I, m_vox2vox);
  MatrixFree(&m_I) ;

  // store things we want to preserve in new gcam
  *(&atlas) = *(&gcam->atlas) ;
  *(&image) = *(&gcam->image) ;
  gca = gcam->gca ;
  type = gcam->type ;
  exp_k = gcam->exp_k ;
  ninputs = gcam->ninputs ;
  spacing = gcam->spacing ;

  GCAMfreeContents(gcam) ; /*MRIfree(&parms->mri_binary) ;*/
  gcam_new = GCAMalloc(mri_labeled->width, mri_labeled->height,mri_labeled->depth);


  *gcam = *gcam_new ;  // copy new one into current one

  // restore saved settings
  gcam->exp_k = exp_k ;
  *(&gcam->atlas) = *(&atlas) ;
  *(&gcam->image) = *(&image) ;
  gcam->ninputs = ninputs ;
  gcam->spacing = spacing ;
  gcam->gca = gca ;
  gcam->type = type ;

  free(gcam_new) ;

  if (Gx < 0 || Gy < 0 || Gz < 0 || Gx >= gcam->width ||
      Gy >= gcam->height || Gz >= gcam->depth)
  {
    Gx = Gy = Gz = -1 ;
  }
  GCAMinit(gcam, mri_labeled, NULL, transform, 0) ;
  GCAMinitLabels(gcam, mri_labeled) ;

  GCAMinitVolGeom(gcam, mri_labeled, mri_dst) ;
  gcamReadMRI(gcam, mri_xn, GCAM_NODEX) ;
  gcamReadMRI(gcam, mri_yn, GCAM_NODEY) ;
  gcamReadMRI(gcam, mri_zn, GCAM_NODEZ) ;
  gcamReadMRI(gcam, mri_origx, GCAM_ORIGX) ;
  gcamReadMRI(gcam, mri_origy, GCAM_ORIGY) ;
  gcamReadMRI(gcam, mri_origz, GCAM_ORIGZ) ;
  if (use_means)
  {
#if 1
    for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
      {
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          if (gcamn->label == 0)
          {
            continue ;
          }
          gcamn->gc = alloc_gcs(1, GCA_NO_MRF, gcam->ninputs) ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
        }
      }
    }
    gcamReadMRI(gcam, mri_means, GCAM_MEANS) ;
    gcamReadMRI(gcam, mri_covars, GCAM_COVARS) ;
#else
    for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
      {
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          if (gcamn->label == 0 || gcamn->gc == NULL)
          {
            continue ;
          }
          gcamn->gc = alloc_gcs(1, GCA_NO_MRF, gcam->ninputs) ;
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          gcamn->gc->means[0] = means[gcamn->label] ;
          gcamn->gc->covars[0] = vars[gcamn->label] ;
        }
      }
    }
#endif
    MRIfree(&mri_means) ;
    MRIfree(&mri_covars) ;
  }
  MRIfree(&mri_origx) ;
  MRIfree(&mri_origy) ;
  MRIfree(&mri_origz) ;
  MRIfree(&mri_xn) ;
  MRIfree(&mri_yn) ;
  MRIfree(&mri_zn) ;
  TransformFree(&transform) ;
  if (pmri)
  {
    *pmri = mri_labeled ;
  }
  else
  {
    MRIfree(&mri_labeled);
  }
  if (DIAG_VERBOSE_ON)
  {
    write_snapshot(gcam, mri_dst, parms, 101) ;
  }


  {
    int x, y, z, xi, yi, zi, xk, yk, zk,num ;
    GCA_MORPH_NODE *gcamn ;

    for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
      {
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          if (gcamn->label == 0)
          {
            continue ;
          }
          for (num = 0, xk = -1 ; xk <= 1 ; xk++)
          {
            xi = x+xk ;
            if (xi < 0 || xi >= gcam->width)
            {
              continue ;
            }
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = y+yk ;
              if (yi < 0 || yi >= gcam->height)
              {
                continue ;
              }
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = z+zk ;
                if (zi < 0 || zi >= gcam->depth)
                {
                  continue ;
                }
                if (gcam->nodes[xi][yi][zi].label >0)
                {
                  num++ ;
                }
              }
            }
          }
          if (num <= 1)
          {
            DiagBreak() ;
          }
        }
      }
    }
  }

  gcamComputeMetricProperties(gcam) ;
  return(gcam) ;
}

MRI *
GCAMinitDensities(GCA_MORPH *gcam, MRI *mri_lowres_seg, MRI *mri_intensities,
                  GCAM_LABEL_TRANSLATION_TABLE *gcam_ltt)
{
  GCA_MORPH_NODE *gcamn ;
  int            i, in_label, out_label, bin, x, y, z, x1, y1, z1, x2, y2, z2, n ;
  MATRIX         *m_vox2vox ;
  MRI            *mri_xformed_aseg ;
  HISTOGRAM      *h, *hsmooth;
  float          val, means[1000], mean ;
  MRI_REGION     bbox ;

  gcam->ninputs = 1 ;
  memset(means, 0, sizeof(means)) ;

  /* transform gcam to be a node->vox map of this intensity image */
  GCAMrasToVox(gcam, mri_intensities) ;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_lowres_seg, mri_intensities) ;
  mri_xformed_aseg = MRIclone(mri_intensities, NULL) ;
  MRIlinearTransformInterp(mri_lowres_seg, mri_xformed_aseg, m_vox2vox, SAMPLE_NEAREST) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_intensities, "int.mgz") ;
    MRIwrite(mri_xformed_aseg, "aseg_xformed.mgz") ;
  }

  /* compute bounding box for hires labels on current image */
  x1 = mri_intensities->width ;
  y1 = mri_intensities->height ;
  z1 = mri_intensities->depth ;
  x2 = y2 = z2 = -1 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0)
        {
          continue ;
        }
        if (gcamn->y < 53)
        {
          DiagBreak() ;
        }
        if (gcamn->x < x1)
        {
          x1 = MAX(0, floor(gcamn->x)) ;
        }
        if (gcamn->y < y1)
        {
          y1 = MAX(0, floor(gcamn->y)) ;
        }
        if (gcamn->z < z1)
        {
          z1 = MAX(0, floor(gcamn->z)) ;
        }

        if (gcamn->x > x2)
        {
          x2 = MIN(mri_intensities->width-1, ceil(gcamn->x)) ;
        }
        if (gcamn->y > y2)
        {
          y2 = MIN(mri_intensities->height-1, ceil(gcamn->y)) ;
        }
        if (gcamn->z > z2)
        {
          z2 = MIN(mri_intensities->depth-1, ceil(gcamn->z)) ;
        }
        if (gcamn->label > 0)
        {
          DiagBreak() ;
        }
      }
    }
  }
  bbox.x = x1 ;
  bbox.y = y1 ;
  bbox.z = z1 ;
  bbox.dx = (x2-x1+1) ;
  bbox.dy = (y2-y1+1) ;
  bbox.dz = (z2-z1+1) ;

  /* expand it by 1 cm in all directions */
  bbox.x -= nint(10*mri_intensities->xsize) ;
  bbox.dx += nint(10*mri_intensities->xsize) ;
  bbox.y -= nint(10*mri_intensities->ysize) ;
  bbox.dy += nint(10*mri_intensities->ysize) ;
  bbox.z -= nint(10*mri_intensities->zsize) ;
  bbox.dz += nint(10*mri_intensities->zsize) ;
  MRIcropBoundingBox(mri_intensities, &bbox) ;

  for (i = 0 ; i < gcam_ltt->nlabels ; i++)
  {
    in_label = gcam_ltt->input_labels[i] ;
    out_label = gcam_ltt->output_labels[i] ;
    if (in_label == Gdiag_no)
    {
      DiagBreak() ;
    }
    if (FZERO(gcam_ltt->means[i])) /* estimate it from image */
    {
      MRI *mri_tmp, *mri_eroded ;
      mri_tmp = MRIclone(mri_xformed_aseg, NULL) ;
      MRIcopyLabel(mri_xformed_aseg, mri_tmp, out_label) ;
      mri_eroded = MRIerode6(mri_tmp, NULL) ;  // account for residual misregistration/distortion
      MRIfree(&mri_tmp) ;
      if (MRIvoxelsInLabel(mri_eroded, out_label) >= 5000)
      {
        mri_tmp = MRIerode6(mri_eroded, NULL) ;
        MRIfree(&mri_eroded) ;
        mri_eroded = mri_tmp ;
      }
      if (MRIvoxelsInLabel(mri_eroded, out_label) >= 25)
      {
        h = MRIhistogramLabelRegion(mri_intensities, mri_eroded, &bbox, out_label, 0) ;
        HISTOfillHoles(h) ;
      }
      else
      {
        h = MRIhistogramLabelRegion(mri_intensities, mri_xformed_aseg, &bbox, out_label, 0) ;
      }
      MRIfree(&mri_eroded) ;
      HISTOclearZeroBin(h) ;
      hsmooth = HISTOsmooth(h, NULL, 2) ;
      HISTOplot(h, "h.plt") ;
      HISTOplot(hsmooth, "hsmooth.plt") ;
      if (IS_LAT_VENT(out_label))
      {
        bin = HISTOfindFirstPeak(hsmooth, 3, .1) ;
      }
      else
      {
        bin = HISTOfindHighestPeakInRegion(hsmooth, 0, hsmooth->nbins) ;
      }
      val = hsmooth->bins[bin] ;
      HISTOfree(&h) ;
      HISTOfree(&hsmooth) ;
      if (out_label != in_label && FZERO(means[out_label]) && (gcam_ltt->second_labels[i] != 0))
        printf("setting intensity of label %d (%s) to %2.1f\n", out_label,
               cma_label_to_name(out_label), val) ;

      means[out_label] = val ;
      means[in_label] = gcam_ltt->scales[i]*val ;
    }
    else
    {
      means[in_label] = gcam_ltt->means[i] ;
    }
    if (gcam_ltt->second_labels[i] == 0) // otherwise will reestimate later
      printf("setting intensity of label %d (%s-->%s) to %2.1f\n", in_label,
             cma_label_to_name(in_label), cma_label_to_name(out_label), means[in_label]) ;
  }

  // now estimate densities for classes that are linear combinations of others (e.g. choroid)
  for (i = 0 ; i < gcam_ltt->nlabels ; i++)
  {
    in_label = gcam_ltt->input_labels[i] ;
    out_label = gcam_ltt->output_labels[i] ;
    if (in_label == Gdiag_no)
    {
      DiagBreak() ;
    }
    if (gcam_ltt->second_labels[i] > 0)
    {
      double pct ;
      pct = gcam_ltt->second_label_pct[i] ;
      means[in_label] = (1-pct)*means[in_label] + pct*means[gcam_ltt->second_labels[i]] ;
      printf("setting linear combo intensity of label %d (%s-->%s) to %2.1f\n", in_label,
             cma_label_to_name(in_label), cma_label_to_name(out_label), means[in_label]) ;
    }
  }

  /*  MRIfree(&mri_xformed_aseg) ;*/
  MatrixFree(&m_vox2vox) ;

#if 0
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (gcamn->label == left_CA2_3)
        {
          gcamn->label = right_CA2_3 ;  /// only for now!!!!!!!!!!!!!!!
        }
      }
    }
  }
#endif

  for (mean = 0.0, n = x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (gcamn->label == 0)
        {
          continue ;
        }
        gcamn->gc = alloc_gcs(1, GCA_NO_MRF, 1) ;
        if (FZERO(means[gcamn->label])) // hasn't been estimated yet
        {
          h = MRIhistogramLabelRegion(mri_intensities, mri_xformed_aseg, &bbox, gcamn->label, 0) ;
          HISTOclearZeroBin(h) ;
          hsmooth = HISTOsmooth(h, NULL, 2) ;
          HISTOplot(h, "h.plt") ;
          HISTOplot(hsmooth, "hsmooth.plt") ;
          bin = HISTOfindHighestPeakInRegion(hsmooth, 0, hsmooth->nbins) ;
          val = hsmooth->bins[bin] ;
          means[gcamn->label] = val ;
          printf("setting intensity of label %d (%s) to %2.1f\n", gcamn->label,
                 cma_label_to_name(gcamn->label), means[gcamn->label]) ;
        }

        gcamn->gc->means[0] = means[gcamn->label] ;
        gcamn->gc->covars[0] = SQR(0.05*gcamn->gc->means[0]) ;
        mean += gcamn->gc->means[0] ;
        n++ ;
        if (FZERO(gcamn->gc->covars[0]))
        {
          DiagBreak() ;
          gcamn->gc->covars[0] = 5 ;
        }
      }
    }
  }

  if (n > 0)
  {
    mean /= n ;
  }
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (gcamn->label == 0)
        {
          continue ;
        }
        if (FZERO(gcamn->gc->means[0]))
        {
          DiagBreak() ;
        }
        gcamn->gc->covars[0] = SQR(0.05 * mean) ;
      }
    }
  }
  return(mri_xformed_aseg);
}
GCA_MORPH *
GCAMlinearTransform(GCA_MORPH *gcam_src, MATRIX *m_vox2vox, GCA_MORPH *gcam_dst)
{
  int             x, y, z ;
  VECTOR          *v1, *v2 ;
  GCA_MORPH_NODE  *gcamn_src, *gcamn_dst ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  v2 = MatrixCopy(v1, NULL) ;

  for (x = 0 ; x < gcam_src->width ; x++)
  {
    for (y = 0 ; y < gcam_src->height ; y++)
    {
      for (z = 0 ; z < gcam_src->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn_src= &gcam_src->nodes[x][y][z] ;
        gcamn_dst= &gcam_dst->nodes[x][y][z] ;
        V3_X(v1) = gcamn_src->x ;
        V3_Y(v1) = gcamn_src->y ;
        V3_Z(v1) = gcamn_src->z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        gcamn_dst->x = V3_X(v2) ;
        gcamn_dst->y = V3_Y(v2) ;
        gcamn_dst->z = V3_Z(v2) ;
      }
    }
  }

  VectorFree(&v1) ;
  VectorFree(&v2) ;
  return(gcam_dst) ;
}

int
GCAMrasToVox(GCA_MORPH *gcam, MRI *mri)
{
  MATRIX   *m ;
  int      x, y, z ;
  GCA_MORPH_NODE  *gcamn ;
  VECTOR          *v1, *v2 ;

  if (gcam->type == GCAM_VOX)
  {
    GCAMvoxToRas(gcam) ;  /* convert it to RAS coords */
  }

  if (mri == NULL)
  {
    m = VGgetRasToVoxelXform(&gcam->image, NULL, 0) ;
  }
  else
  {
    m = MRIgetRasToVoxelXform(mri) ;
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("ras 2 vox matrix:\n") ;
    MatrixPrint(stdout, m) ;
  }

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  v2 = MatrixCopy(v1, NULL) ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label > 0)
        {
          DiagBreak() ;
        }

        V3_X(v1) = gcamn->x ;
        V3_Y(v1) = gcamn->y ;
        V3_Z(v1) = gcamn->z ;
        MatrixMultiply(m, v1, v2) ;
        gcamn->x = V3_X(v2) ;
        gcamn->y = V3_Y(v2) ;
        gcamn->z = V3_Z(v2) ;
        if (nint(gcamn->x) == Gx && nint(gcamn->y) == Gy && nint(gcamn->z) == Gz)
        {
          DiagBreak() ;
        }

        V3_X(v1) = gcamn->origx ;
        V3_Y(v1) = gcamn->origy ;
        V3_Z(v1) = gcamn->origz ;
        MatrixMultiply(m, v1, v2) ;
        gcamn->origx = V3_X(v2) ;
        gcamn->origy = V3_Y(v2) ;
        gcamn->origz = V3_Z(v2) ;
        if (gcamn->invalid == GCAM_AREA_INVALID)
        {
          gcamn->invalid = GCAM_VALID ;
        }
      }
    }
  }

  GCAMcomputeOriginalProperties(gcam) ;  // orig node positions have changed to vox coords
  gcam->type = GCAM_VOX ;
  MatrixFree(&m) ;
  if (mri)
  {
    getVolGeom(mri, &gcam->image) ;
  }

  return(NO_ERROR)  ;
}

int
GCAMvoxToRas(GCA_MORPH *gcam)
{
  MATRIX   *m ;
  int      x, y, z ;
  GCA_MORPH_NODE  *gcamn ;
  VECTOR          *v1, *v2 ;

  if (gcam->type == GCAM_RAS)
  {
    return(NO_ERROR) ;
  }
  m = vg_getVoxelToRasXform(&gcam->image) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("vox 2 ras matrix:\n") ;
    MatrixPrint(stdout, m) ;
  }

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v1, 4, 1) = 1.0 ;
  v2 = MatrixCopy(v1, NULL) ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label > 0)
        {
          DiagBreak() ;
        }
        V3_X(v1) = gcamn->x ;
        V3_Y(v1) = gcamn->y ;
        V3_Z(v1) = gcamn->z ;
        MatrixMultiply(m, v1, v2) ;
        gcamn->x = V3_X(v2) ;
        gcamn->y = V3_Y(v2) ;
        gcamn->z = V3_Z(v2) ;

        V3_X(v1) = gcamn->origx ;
        V3_Y(v1) = gcamn->origy ;
        V3_Z(v1) = gcamn->origz ;
        MatrixMultiply(m, v1, v2) ;
        gcamn->origx = V3_X(v2) ;
        gcamn->origy = V3_Y(v2) ;
        gcamn->origz = V3_Z(v2) ;
      }
    }
  }

  gcam->type = GCAM_RAS ;
  MatrixFree(&m) ;
  return(NO_ERROR) ;
}


NODE_LOOKUP_TABLE *
gcamCreateNodeLookupTable(GCA_MORPH *gcam, MRI *mri, NODE_LOOKUP_TABLE *nlt)
{
  int                 x, y, z, xi, yi, zi, index ;
  GCA_MORPH_NODE      *gcamn ;
  MRI                 *mri_tmp ;
  NODE_BUCKET         *nb ;

  if (nlt)
  {
    gcamFreeNodeLookupTable(&nlt) ;
  }
  nlt = (NODE_LOOKUP_TABLE *)calloc(1, sizeof(NODE_LOOKUP_TABLE)) ;
  if (!nlt)
  {
    ErrorExit(ERROR_NOMEMORY, "gcamCreateNodeLookupTable: could not allocate NLT");
  }
  nlt->width = 2*NLT_PAD+mri->width ;
  nlt->height = 2*NLT_PAD+mri->height;
  nlt->depth = 2*NLT_PAD+mri->depth;
  nlt->nodes = (NODE_BUCKET ***)calloc(nlt->width, sizeof(NODE_BUCKET **)) ;
  if (!nlt->nodes)
  {
    ErrorExit(ERROR_NOMEMORY, "gcamCreateNodeLookupTable: could not allocate ***") ;
  }

  for (x = 0 ; x < nlt->width ; x++)
  {
    nlt->nodes[x] = (NODE_BUCKET **)calloc(nlt->height, sizeof(NODE_BUCKET *)) ;
    if (!nlt->nodes[x])
    {
      ErrorExit(ERROR_NOMEMORY, "gcamCreateNodeLookupTable: could not allocate %dth **",x) ;
    }

    for (y = 0 ; y < nlt->height ; y++)
    {
      nlt->nodes[x][y] = (NODE_BUCKET *)calloc(nlt->depth, sizeof(NODE_BUCKET)) ;
      if (!nlt->nodes[x][y])
      {
        ErrorExit(ERROR_NOMEMORY,"gcamCreateNodeLookupTable: could not allocate %d,%dth *",x,y);
      }
    }
  }

  for (x = 0 ; x < gcam->width; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0)
        {
          continue ;  /* only labeled nodes */
        }
        xi = nint(gcamn->x)+NLT_PAD ;
        yi = nint(gcamn->y)+NLT_PAD ;
        zi = nint(gcamn->z)+NLT_PAD ;
        xi = MAX(0, xi) ;
        xi = MIN(nlt->width-1, xi) ;
        yi = MAX(0, yi) ;
        yi = MIN(nlt->height-1, yi) ;
        zi = MAX(0, zi) ;
        zi = MIN(nlt->depth-1, zi) ;
        if (xi == Gx && yi == Gy && zi == Gz)
        {
          DiagBreak() ;
        }
        nlt->nodes[xi][yi][zi].nnodes++ ;
      }

  /* now allocate the node buckets */
  mri_tmp = MRIalloc(nlt->width, nlt->height, nlt->depth, MRI_SHORT) ;
  for (x = 0 ; x < nlt->width; x++)
    for (y = 0 ; y < nlt->height ; y++)
      for (z = 0 ; z < nlt->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (nlt->nodes[x][y][z].nnodes == 0)
        {
          nlt->nodes[x][y][z].node_bins = NULL ;
          continue ;
        }
        MRISvox(mri_tmp, x, y, z) = nlt->nodes[x][y][z].nnodes ;
        nlt->nodes[x][y][z].node_bins =
          (NODE_BIN *)calloc(nlt->nodes[x][y][z].nnodes,sizeof(NODE_BIN)) ;
        if (NULL == nlt->nodes[x][y][z].node_bins)
          ErrorExit
          (ERROR_NOMEMORY, \
           "gcamCreateNodeLookupTable: "
           "could not allocate %d bins at (%d, %d, %d)", \
           nlt->nodes[x][y][z].nnodes,x,y,z) ;
        nlt->nodes[x][y][z].nnodes = 0 ;  /* reset counts to 0 */
      }



  /* go through each gcam node and find
     where it points in the image, then add the node
     to the lookup table */
  for (x = 0 ; x < gcam->width; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        int max ;
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0)
        {
          continue ;  /* only labeled nodes */
        }
        xi = nint(gcamn->x)+NLT_PAD ;
        yi = nint(gcamn->y)+NLT_PAD ;
        zi = nint(gcamn->z)+NLT_PAD ;
        xi = MAX(0, xi) ;
        xi = MIN(nlt->width-1, xi) ;
        yi = MAX(0, yi) ;
        yi = MIN(nlt->height-1, yi) ;
        zi = MAX(0, zi) ;
        zi = MIN(nlt->depth-1, zi) ;
        max = MRISvox(mri_tmp, xi, yi, zi) ;
        index = nlt->nodes[xi][yi][zi].nnodes++ ;
        if (index >= max)
        {
          DiagBreak() ;
        }
        if (xi < 0 || xi >= nlt->width ||
            yi < 0 || yi >= nlt->height ||
            zi < 0 || zi >= nlt->depth)
        {
          DiagBreak() ;
        }
        nb = &nlt->nodes[xi][yi][zi] ;
        nb->node_bins[index].x = x ;
        // store for mapping image voxel -> node later
        nb->node_bins[index].y = y ;
        nb->node_bins[index].z = z ;
      }
  MRIfree(&mri_tmp) ;
  return(nlt) ;
}


int
gcamFreeNodeLookupTable(NODE_LOOKUP_TABLE **pnlt)
{
  NODE_LOOKUP_TABLE  *nlt ;
  int                x, y, z ;

  nlt = *pnlt ;
  *pnlt = NULL ;

  for (x = 0 ; x < nlt->width ; x++)
  {
    for (y = 0 ; y < nlt->height ; y++)
    {
      for (z = 0 ; z < nlt->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (nlt->nodes[x][y][z].node_bins == NULL)
        {
          continue ;
        }
        free(nlt->nodes[x][y][z].node_bins) ;
      }
      free(nlt->nodes[x][y]) ;
    }
    free(nlt->nodes[x]) ;
  }
  free(nlt->nodes) ;
  free(nlt) ;
  return(NO_ERROR) ;
}
int
GCAMresetLabelNodeStatus(GCA_MORPH *gcam)
{
  GCAMremoveStatus(gcam, GCAM_LABEL_NODE) ;
  GCAMremoveStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ;
  return(NO_ERROR) ;
}

int
GCAMaddStatus(GCA_MORPH *gcam, int status)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        gcamn->status |= status ;
      }
  return(NO_ERROR) ;
}

int
GCAMremoveStatus(GCA_MORPH *gcam, int status)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  status = ~status ;  /* remove the original status bits */
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        gcamn->status &= status ;
      }
  return(NO_ERROR) ;
}


#if 0
int
gcamSetGradientToNbrAverage(GCA_MORPH *gcam, int x, int y, int z)
{
  GCA_MORPH_NODE *gcamn, *gcamn_nbr ;
  double         dx, dy, dz ;
  int            num, xi, yi, zi, xk, yk, zk ;

  gcamn = &gcam->nodes[x][y][z];
  num = 0 ;
  dx = dy = dz = 0.0 ;
  for (zk = -1 ; zk <= 1 ; zk++)
  {
    zi = z+zk ;
    if ((zi <0) || (zi >= gcam->depth))
    {
      continue ;
    }
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = y+yk ;
      if ((yi <0) || (yi >= gcam->height))
      {
        continue ;
      }
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = x+xk ;
        if ((xi <0) || (xi >= gcam->width))
        {
          continue ;
        }
        gcamn_nbr = &gcam->nodes[xi][yi][zi];
        if (gcamn_nbr->area < 0)
        {
          continue ;
        }
        dx += gcamn_nbr->dx ;
        dy += gcamn_nbr->dy ;
        dz += gcamn_nbr->dz ;
        num++ ;
      }
    }
  }
  if (num == 0)
  {
    return(NO_ERROR) ;
  }

  gcamn->dx  = dx / num ;
  gcamn->dy  = dy / num ;
  gcamn->dz  = dz / num ;

  return(NO_ERROR) ;
}
#endif
int
GCAMremoveCompressedRegions(GCA_MORPH *gcam, float min_ratio)
{
  int             ncompressed_before, ncompressed_after, niter ;
  GCA_MORPH_PARMS parms ;

  ncompressed_before =
    ncompressed_after = GCAMcountCompressedNodes(gcam, min_ratio) ;
  if (ncompressed_before <= 0)
  {
    return(NO_ERROR) ;
  }

  if (Gdiag & DIAG_SHOW)
    printf("%d nodes compressed more than %2.2f - "
           "using Jacobian term to decompress\n",
           ncompressed_before, min_ratio) ;
  memset(&parms, 0, sizeof(parms)) ;

  parms.l_jacobian = 1.0 ;
  parms.dt = 0.005 ;
  parms.noneg = True ;
  parms.exp_k = 20 ;
  parms.levels = 1 ;
  parms.sigma = 1.0 ;
  parms.relabel_avgs = -1 ;
  parms.integration_type = GCAM_INTEGRATE_OPTIMAL ;
  parms.nsmall = 1 ;
  parms.reset_avgs = -1 ;
  parms.tol = 1e-3 ;
  parms.niterations = 20 ;
  strcpy(parms.base_name, "jacobian") ;
  niter = 0 ;
  do
  {
    ncompressed_before = ncompressed_after ;
    GCAMregister(gcam, NULL, &parms) ;
    ncompressed_after = GCAMcountCompressedNodes(gcam, min_ratio) ;
    if (Gdiag & DIAG_SHOW)
      printf("after decompression - %d nodes compressed more than %2.2f\n", \
             ncompressed_after, min_ratio) ;
    if (niter++ > 5)
    {
      break ;
    }
  }
  while (ncompressed_after < ncompressed_before) ;
  return(NO_ERROR) ;
}

int
GCAMcountCompressedNodes(GCA_MORPH *gcam, float min_ratio)
{
  int             i, j, k, ncompressed, compressed ;
  GCA_MORPH_NODE  *gcamn ;
  float           ratio1, ratio2 ;

  gcamComputeMetricProperties(gcam) ;
  for (ncompressed = i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }
        if (gcamn->invalid == GCAM_AREA_INVALID ||
            gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue ;
        }

        compressed = 0 ;
        if (!DZERO(gcamn->orig_area1))
        {
          ratio1 = gcamn->area1 / gcamn->orig_area1 ;
          if (ratio1 < min_ratio)
          {
            compressed = 1 ;
          }
        }
        if (!DZERO(gcamn->orig_area2))
        {
          ratio2 = gcamn->area2 / gcamn->orig_area2 ;
          if (ratio2 < min_ratio)
          {
            compressed = 1 ;
          }
        }
        if (compressed)
        {
          ncompressed++ ;
        }
      }
    }
  }
  return(ncompressed) ;
}

int
gcamComputeMostLikelyDirection(GCA_MORPH *gcam, MRI *mri,
                               double x0, double y0, double z0,
                               double target, MRI *mri_kernel, MRI *mri_nbhd,
                               double *pdx, double *pdy, double *pdz)
{
  int   xi, yi, zi, xk, yk, zk, x, y, z, whalf, x1, y1, z1 ;
  double val, max_val, max_dist, dist, max_x, max_y, max_z ;
  static MRI *mri_int = NULL ;

  whalf = (mri_nbhd->width-1)/2 ;
  x = nint(x0) ;
  y = nint(y0) ;
  z = nint(z0) ;
  MRIclear(mri_nbhd) ;
  MRIsampleVolume(mri, x0, y0, z0, &val) ;
  if (mri_int == NULL)
  {
    mri_int =
      MRIextractInto
      (mri, mri_int, x0-whalf, y0-whalf, z0-whalf,
       mri_nbhd->width,mri_nbhd->height, mri_nbhd->depth, 0, 0, 0) ;
    MRIcopyHeader(mri_int, mri_nbhd) ;
  }
  // transform origin to local coordinates
  x0 = x0 - x + whalf ;
  y0 = y0 - y + whalf ;
  z0 = z0 - z + whalf ;

  max_val = -SQR(val - target) ;
  max_dist = 0 ;
  max_x = x0 ;
  max_y = y0 ;
  max_z = z0 ;

  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = x+xk ;
    if (xi < 0 || xi >= mri->width)
    {
      continue ;
    }
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi =y+yk ;
      if (yi < 0 || yi >= mri->height)
      {
        continue ;
      }
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi =z+zk ;
        if (zi < 0 || zi >= mri->depth)
        {
          continue ;
        }
        val = MRIgetVoxVal(mri, xi, yi, zi, 0) ;
        x1 = xk+whalf ;
        y1 = yk+whalf ;
        z1 = zk+whalf ;  // mri_nbhd voxel coords
        MRIsetVoxVal(mri_int, x1, y1, z1, 0, val) ;
        val = -SQR(val - target) ;  // sort of log(likelihood)
        MRIsetVoxVal(mri_nbhd, x1, y1, z1, 0, val) ;
        if (FEQUAL(val, max_val) ||
            (val > max_val))  // if same likelihood, pick closer one
        {
          dist = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0)) ;
          if ((val > max_val) || (dist < max_dist))  // new maximum
          {
            max_x = (double)x1 ;
            max_y = (double)y1 ;
            max_z = (double)z1 ;
            max_val = val ;
            max_dist = dist ;
          }
        }
      }
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_nbhd, "n.mgz") ;
  }
  // exact position of origin
#if 1
  MRIconvolveGaussian(mri_nbhd, mri_nbhd, mri_kernel) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_int, "ni.mgz") ;
    MRIwrite(mri_nbhd, "ns.mgz") ;
  }

  MRIsampleVolumeGradientFrame(mri_nbhd, x0, y0, z0, pdx, pdy, pdz, 0) ;
#else
  *pdx = max_x - x0 ;
  *pdy = max_y - y0 ;
  *pdz = max_z - z0 ;
#endif
  return(NO_ERROR) ;
}

int
GCAMaddIntensitiesFromImage(GCA_MORPH *gcam, MRI *mri_target)
{
  int             x, y, z, num ;
  double            val, mean ;
  GCA_MORPH_NODE  *gcamn ;

  if (gcam->ninputs == 0)
  {
    gcam->ninputs = 1 ;
  }
  num = 0 ;
  mean = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        val = MRIgetVoxVal(mri_target, x, y, z, 0) ;
#if 0
        if (FZERO(val)/* && MRIneighborsOn3x3(mri_target, x, y, z, 1) == 0*/)
        {
          continue ;
        }
#endif
        gcamn = &gcam->nodes[x][y][z] ;
	if (x == 0 || x == gcam->width-1 ||
	    y == 0 || y == gcam->height-1 ||
	    z == 0 || z == gcam->depth-1)
	  gcamn->invalid = GCAM_AREA_INVALID ;
	  
        gcamn->gc = alloc_gcs(1, GCA_NO_MRF, gcam->ninputs) ;
        gcamn->gc->means[0] = val ;
        if (!DZERO(val))
        {
          num++ ;
          mean += val ;
          gcamn->label = 128 ;   // something real here
        }
        else
        {
          DiagBreak() ;
        }
        gcamn->gc->covars[0] = SQR(0.05*val) ;  // assume SNR=20 or so
      }
    }
  }

  if (num > 0)  // fill in 0 variances
  {
    mean /= num ;
    for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
      {
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          if (gcamn->gc == NULL)
          {
            continue ;
          }
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          //     if (gcamn->gc->covars[0] < 0.05*mean)
          gcamn->gc->covars[0] = SQR(0.05*mean) ;
        }
      }
    }
  }
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->gc == NULL)
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        //     if (gcamn->gc->covars[0] < 0.05*mean)
        //        if (gcamn->gc->covars[0] < MIN_VAR)
        gcamn->gc->covars[0] = MIN_VAR ;
      }
    }
  }

  return(NO_ERROR) ;
}

GCA_MORPH *
GCAMcreateFromIntensityImage
(MRI *mri_source, MRI *mri_target, TRANSFORM *transform)
{
  GCA_MORPH       *gcam ;
  int             x, y, z, num ;
  double            val, mean ;
  GCA_MORPH_NODE  *gcamn ;

  gcam = GCAMalloc(mri_target->width, mri_target->height, mri_target->depth) ;
  GCAMinit(gcam, mri_source, NULL, transform, 0) ;

  num = 0 ;
  mean = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        val = MRIgetVoxVal(mri_target, x, y, z, 0) ;
#if 0
        if (FZERO(val)/* && MRIneighborsOn3x3(mri_target, x, y, z, 1) == 0*/)
        {
          continue ;
        }
#endif
        gcamn = &gcam->nodes[x][y][z] ;
        gcamn->gc = alloc_gcs(1, GCA_NO_MRF, gcam->ninputs) ;
        gcamn->gc->means[0] = val ;
        if (!DZERO(val))
        {
          num++ ;
          mean += val ;
          gcamn->label = 128 ;   // something real here
        }
        else
        {
          DiagBreak() ;
        }
        gcamn->gc->covars[0] = SQR(0.05*val) ;  // assume SNR=20 or so
      }
    }
  }

  if (num > 0)  // fill in 0 variances
  {
    mean /= num ;
    for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
      {
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          if (gcamn->gc == NULL)
          {
            continue ;
          }
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          //     if (gcamn->gc->covars[0] < 0.05*mean)
          gcamn->gc->covars[0] = SQR(0.05*mean) ;
        }
      }
    }
  }
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->gc == NULL)
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        //     if (gcamn->gc->covars[0] < 0.05*mean)
        //        if (gcamn->gc->covars[0] < MIN_VAR)
        gcamn->gc->covars[0] = MIN_VAR ;
      }
    }
  }

  GCAMinitVolGeom(gcam, mri_source, mri_target) ;
  // GCAMinitLabels(gcam, NULL) ;
  return(gcam) ;
}
int
GCAMsetVariances(GCA_MORPH *gcam, float var)
{
  int             x, y, z, r ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->gc == NULL)
        {
          continue ;
        }
        for (r = 0 ; r < gcam->ninputs ; r++)
        {
          gcamn->gc->covars[r] = var ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

int
GCAMthresholdLikelihoodStatus(GCAM *gcam, MRI *mri, float thresh)
{
  int             x, y, z ;
  double            val ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        if (val <= thresh)
        {
          gcamn->status |= (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD) ;
        }
        else
        {
          DiagBreak() ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

int
gcamConstrainJacobian(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  double   alpha, delta, lower_alpha, upper_alpha, rms,
           min_jacobian, max_jacobian, start_rms ;
  int      done = 0 ;

  start_rms = GCAMcomputeRMS(gcam, mri, parms) ;

  min_jacobian = parms->ratio_thresh ;
  max_jacobian = 1/min_jacobian ;

  // find the max alpha that allows an admissable jacobian
  if (gcamCheckJacobian(gcam, min_jacobian, max_jacobian))
  {
    return(NO_ERROR) ;
  }

  // not currently valid, reduce total movement
  GCAMcopyNodePositions
  (gcam, CURRENT_POSITIONS, SAVED_POSITIONS) ;  // store current

  upper_alpha = 1 ;

  // find lower alpha
  for (lower_alpha = upper_alpha * .75 ; lower_alpha > 0 ; lower_alpha *= .75)
  {
    gcamSetTotalMovement(gcam, lower_alpha) ;
    if (gcamCheckJacobian(gcam, min_jacobian, max_jacobian))
    {
      break ;  // found a valid movement
    }
  }
  rms = GCAMcomputeRMS(gcam, mri, parms) ;

  do
  {
    // binary search
    alpha = (upper_alpha + lower_alpha) / 2 ;
    gcamSetTotalMovement(gcam, alpha) ;
    if (gcamCheckJacobian(gcam, min_jacobian, max_jacobian))
    {
      delta = alpha - lower_alpha ;
      lower_alpha = alpha ;
      rms = GCAMcomputeRMS(gcam, mri, parms) ;
    }
    else  // invalid movement - make upper alpha lower
    {
      delta = upper_alpha - alpha ;
      upper_alpha = alpha ;
    }
    done = 100*delta < parms->tol ;
  }
  while (!done) ;
  gcamSetTotalMovement(gcam, lower_alpha) ;
  printf("optimal alpha = %2.2f\n", lower_alpha) ;
  return(NO_ERROR) ;
}


int
gcamCheckJacobian(GCA_MORPH *gcam, float min_jacobian, float max_jacobian)
{
  int             i, j, k ;
  GCA_MORPH_NODE  *gcamn ;
  float           ratio ;

  //  gcamComputeMetricProperties(gcam) ;
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }
        if (gcamn->invalid == GCAM_AREA_INVALID ||
            gcamn->invalid == GCAM_POSITION_INVALID
            || DZERO(gcamn->orig_area))
        {
          continue ;
        }

        ratio = gcamn->area / gcamn->orig_area ;
        if (ratio < min_jacobian || ratio > max_jacobian)
        {
          return(0) ;
        }
      }
    }
  }
  return(1) ;
}


int
gcamSetTotalMovement(GCA_MORPH *gcam, float alpha)
{
  int             i, j, k ;
  double          dx, dy, dz ;
  GCA_MORPH_NODE  *gcamn ;

  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        dx = gcamn->xs - gcamn->origx ;
        dy = gcamn->ys - gcamn->origy ;
        dz = gcamn->zs - gcamn->origz ;
        gcamn->x = gcamn->origx + alpha*dx ;
        gcamn->y = gcamn->origy + alpha*dy ;
        gcamn->z = gcamn->origz + alpha*dz ;
      }
    }
  }

  gcamComputeMetricProperties(gcam) ;
  return(NO_ERROR) ;
}

int
GCAMreinitWithLTA(GCA_MORPH *gcam, LTA *lta, MRI *mri, GCA_MORPH_PARMS *parms)
{
  GCA_MORPH_NODE  *gcamn ;
  GCA             *gca ;
  int             x, y, z, width, height, depth, i, level, navgs, first=1 ;
  float           min_dt, orig_dt, last_rms, rms, pct_change, lin_rms ;
  VECTOR          *v_src, *v_dst ;
  MATRIX          *m_avg, *m_L ;
  GCA_MORPH_PARMS lp, *lin_parms = &lp ;

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.log", parms->base_name) ;
    if (parms->log_fp == NULL)
    {
      if (parms->start_t == 0)
      {
        parms->log_fp = fopen(fname, "w") ;
      }
      else
      {
        parms->log_fp = fopen(fname, "a") ;
      }
    }
  }
  else
  {
    parms->log_fp = NULL ;
  }

  *lin_parms = *parms ;
  lin_parms->l_smoothness = lin_parms->l_elastic =
                              lin_parms->l_expansion = lin_parms->l_distance =
                                    lin_parms->l_lsmoothness = lin_parms->l_spring = 0.0 ;
  m_avg = MatrixCopy(lta->xforms[0].m_L, NULL) ;
  for (i = 1 ; i < lta->num_xforms ; i++)
  {
    MatrixAdd(lta->xforms[i].m_L, m_avg, m_avg) ;
  }
  MatrixScalarMul(m_avg, 1.0/(float)lta->num_xforms, m_avg) ;

  gca = gcam->gca ;
  v_src = VectorAlloc(4, MATRIX_REAL) ;
  v_dst = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_src,4) = 1.0 ;
  VECTOR_ELT(v_dst,4) = 1.0 ;

  // save geometry information
  width = gcam->width ;
  height = gcam->height ;
  depth = gcam->depth ;

  // check consistency
  GCAsetVolGeom(gca, &gcam->atlas); // fname is still "unknown"
  if (width != gca->prior_width
      || height != gca->prior_height
      || depth != gca->prior_depth)
  {
    fprintf
    (stderr,
     "GCA_MORPH (%d, %d, %d) must agree with GCA prior (%d, %d, %d)\n",
     gcam->width, gcam->height, gcam->depth, gca->prior_width, \
     gca->prior_height, gca->prior_depth);
    ErrorExit(ERROR_BADPARM, "Exiting ....\n");
  }
  gcam->gca = gca ;
  gcam->spacing = gca->prior_spacing ;
  gcam->exp_k = parms->exp_k ;
  parms->mri = mri ;

  last_rms = GCAMcomputeRMS(gcam, mri, parms) ;
  if (parms->log_fp)
  {
    fprintf(parms->log_fp, "%04d: dt=%2.6f, rms=%2.3f, neg=%d, invalid=%d\n",
            parms->start_t, 0.0, last_rms, gcam->neg, Ginvalid) ;
    fflush(parms->log_fp) ;
  }
  if ((parms->write_iterations > 0) && (Gdiag & DIAG_WRITE))
  {
    write_snapshot(gcam, mri, parms, parms->start_t++) ;  // snap before changing initialization
  }

  // compute target positions from polyaffine transform
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        // see if this label is in the LTA
        for (i = 0 ; i < lta->num_xforms ; i++)
          if (lta->xforms[i].label == gcamn->label)
          {
            break ;
          }

        // now apply this xform to the orig point
        // here mri info is used
        V3_X(v_src) = gcamn->x ;
        V3_Y(v_src) = gcamn->y ;
        V3_Z(v_src) = gcamn->z ;
        MatrixMultiply(m_avg, v_src, v_dst) ;
        if (i < lta->num_xforms)  // measure distance between label specific xform and avg one
        {
          MatrixMultiply(lta->xforms[i].m_L, v_src, v_dst) ;
          gcamn->dx = V3_X(v_dst) - gcamn->x ;
          gcamn->dy = V3_Y(v_dst) - gcamn->y ;
          gcamn->dz = V3_Z(v_dst) - gcamn->z ;
          gcamn->xs = V3_X(v_dst) ;
          gcamn->ys = V3_Y(v_dst) ;
          gcamn->zs = V3_Z(v_dst) ;
          gcamn->status |= GCAM_TARGET_DEFINED ;
        }
        else  // no label specific xform - don't know where it should go
        {
          gcamn->xs = gcamn->x ;
          gcamn->ys = gcamn->y ;
          gcamn->zs = gcamn->z ;
        }

        gcamn->log_p = 0 ;
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
      }
    }
  }

  // recompute original properties
  GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, ORIGINAL_POSITIONS) ;
  gcamComputeMetricProperties(gcam) ;
  GCAMstoreMetricProperties(gcam) ;

  orig_dt = parms->dt ;


  m_L = NULL ;
  for (navgs = parms->navgs, level = parms->levels ; level >= 0 ; level--)
  {
    do
    {
      do
      {
        gcamClearGradient(gcam) ;
        m_L = gcamComputeOptimalTargetLinearTransform(gcam, m_L, .9) ;
        lin_rms = GCAMcomputeRMS(gcam, mri, lin_parms) ;
        gcamApplyLinearTransform(gcam, m_L) ;
        rms = GCAMcomputeRMS(gcam, mri, lin_parms) ;
        pct_change = 100.0*(lin_rms-rms)/lin_rms ;
        if (pct_change <= 0)
        {
          MATRIX *m_inv ;

          printf("rms increased - undoing step...\n") ;
          m_inv = MatrixInverse(m_L, NULL) ;
          gcamApplyLinearTransform(gcam, m_inv) ;
          MatrixFree(&m_inv) ;
          rms = GCAMcomputeRMS(gcam, mri, parms) ;
          MatrixIdentity(4, m_L) ;
        }
        if ((parms->write_iterations > 0) && (Gdiag & DIAG_WRITE))
        {
          write_snapshot(gcam, mri, parms, parms->start_t) ;
        }
        if (parms->log_fp)
        {
          fprintf(parms->log_fp, "%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, LIN",
                  parms->start_t, 0.0, rms, pct_change, gcam->neg, Ginvalid) ;
          if (parms->l_binary > 0)
            fprintf(parms->log_fp, ", aligned = %d (%2.3f%%)\n",
                    Galigned,
                    100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
          else
          {
            fprintf(parms->log_fp, "\n") ;
          }
          fflush(parms->log_fp) ;
        }
        printf("%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, LIN",
               parms->start_t, 0.0, rms, pct_change, gcam->neg, Ginvalid) ;
        if (parms->l_binary > 0)
          printf(", aligned = %d (%2.3f%%)\n",
                 Galigned,
                 100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          printf("\n") ;
        }
        fflush(stdout) ;
        last_rms = rms ;
        gcamCheck(gcam, mri) ;
        parms->start_t++ ;
        if (first)
        {
          /* haven't ever taken a nonlinear step, so can update metric
            properties and use this as initial conditions. */
          GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, ORIGINAL_POSITIONS) ;
          gcamComputeMetricProperties(gcam) ;
          GCAMstoreMetricProperties(gcam) ;
          GCAMcomputeOriginalProperties(gcam) ;
          if (gcam->m_affine)
          {
            double det ;
            MatrixMultiply(m_L, gcam->m_affine,gcam->m_affine) ;
            det = MatrixDeterminant(gcam->m_affine) ;
            printf("det of total affine transform updated to be %2.3f\n", det) ;
          }
        }
      }
      while (pct_change > parms->tol*0.1) ;
      //      GCAMcomputeOriginalProperties(gcam) ;
      first = 0 ;
      gcamClearGradient(gcam) ;
      last_rms = GCAMcomputeRMS(gcam, mri, parms) ;
      gcamComputeTargetGradient(gcam) ;
      gcamSmoothnessTerm(gcam, mri, parms->l_smoothness)  ;
#if 0
      gcamLikelihoodTerm(gcam, mri, mri, parms->l_likelihood, parms)  ;
      gcamLogLikelihoodTerm(gcam, mri, mri, parms->l_log_likelihood)  ;
#endif
      gcamJacobianTerm(gcam, mri, parms->l_jacobian,parms->ratio_thresh)  ;
      gcamLimitGradientMagnitude(gcam, parms, mri) ;

      gcamSmoothGradient(gcam, navgs) ;
      parms->dt = orig_dt ;
      min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
      parms->dt = min_dt ;
      gcamApplyGradient(gcam, parms) ;
      if (Gdiag & DIAG_SHOW)
      {
        gcamShowCompressed(gcam, stdout) ;
      }
      if (parms->uncompress)
      {
        gcamRemoveCompressedNodes(gcam, mri, parms, parms->ratio_thresh) ;
      }
      if ((parms->write_iterations > 0) && (Gdiag & DIAG_WRITE))
      {
        write_snapshot(gcam, mri, parms, parms->start_t) ;
      }
      rms = GCAMcomputeRMS(gcam, mri, parms) ;
      pct_change = 100.0*(last_rms-rms)/last_rms ;
      last_rms = rms ;
      if (parms->log_fp)
      {
        fprintf(parms->log_fp, "%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, PA, navgs=%d",
                parms->start_t, min_dt, rms, pct_change, gcam->neg,
                Ginvalid,navgs) ;
        if (parms->l_binary > 0)
          fprintf(parms->log_fp, ", aligned = %d (%2.3f%%)\n",
                  Galigned,
                  100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          fprintf(parms->log_fp, "\n") ;
        }
        fflush(parms->log_fp) ;
      }

      if (Gdiag & DIAG_SHOW)
      {
        printf("%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, PA, navgs=%d",
               parms->start_t, min_dt, rms, pct_change, gcam->neg, Ginvalid, navgs) ;
        if (parms->l_binary > 0)
          printf(", aligned = %d (%2.3f%%)\n",
                 Galigned,
                 100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          printf("\n") ;
        }
      }
      parms->start_t++ ;
      //      GCAMcomputeOriginalProperties(gcam) ;
    }
    while (pct_change > parms->tol) ;
    navgs /= 4 ;
    if (navgs < 1)
    {
      break ;
    }
  }
  //  GCAMcomputeOriginalProperties(gcam) ;
  parms->dt = orig_dt ;

  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;
  return(NO_ERROR) ;
}


int
gcamComputeTargetGradient(GCA_MORPH *gcam)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->status &
            GCAM_TARGET_DEFINED) // drive it towards the target position
        {
          gcamn->dx = gcamn->xs - gcamn->x ;
          gcamn->dy = gcamn->ys - gcamn->y ;
          gcamn->dz = gcamn->zs - gcamn->z ;
        }
        else  // no label specific xform - don't know where it should go
        {
          gcamn->dx = gcamn->dy = gcamn->dz = 0 ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}


MATRIX *
gcamComputeOptimalTargetLinearTransform
(GCA_MORPH *gcam, MATRIX *m_L, double reg)
{
  int            x, y, z, ncols, col ;
  GCA_MORPH_NODE *gcamn ;
  MATRIX         *m_U, *m_V, *m_Uinv, *m_I ;

  for (ncols = x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->status &
            GCAM_TARGET_DEFINED) // drive it towards the target position
        {
          ncols++ ;
        }
      }
    }
  }

  m_U = MatrixAlloc(4, ncols, MATRIX_REAL) ;
  m_V = MatrixAlloc(4, ncols, MATRIX_REAL) ;
  for (col = x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->status &
            GCAM_TARGET_DEFINED) // drive it towards the target position
        {
          col++ ;
          *MATRIX_RELT(m_U,1,col) = gcamn->x ;
          *MATRIX_RELT(m_U,2,col) = gcamn->y ;
          *MATRIX_RELT(m_U,3,col) = gcamn->z ;
          *MATRIX_RELT(m_U,4,col) = 1.0 ;

          *MATRIX_RELT(m_V,1,col) = gcamn->xs ;
          *MATRIX_RELT(m_V,2,col) = gcamn->ys ;
          *MATRIX_RELT(m_V,3,col) = gcamn->zs ;
          *MATRIX_RELT(m_V,4,col) = 1.0 ;
        }
      }
    }
  }

  // m_Uinv = MatrixPseudoInverse(m_U, NULL) ;
  m_Uinv = MatrixSVDPseudoInverse(m_U, NULL) ;
  // MatrixPrint(stdout, m_Uinv) ;
  m_L = MatrixMultiply(m_V, m_Uinv, m_L) ;
  *MATRIX_RELT(m_L, 4, 1) = 0.0 ;
  *MATRIX_RELT(m_L, 4, 2) = 0.0 ;
  *MATRIX_RELT(m_L, 4, 3) = 0.0 ;
  *MATRIX_RELT(m_L, 4, 4) = 1.0 ;

  if (DIAG_VERBOSE_ON)
  {
    MATRIX *m_tmp ;
    int    i ;
    printf("m_L:\n") ;
    MatrixPrint(stdout, m_L) ;
    m_tmp = MatrixMultiply(m_L, m_U, NULL) ;

    i = m_tmp->cols ;

    printf("U:\n") ;
    m_U->cols = 10 ;
    MatrixPrint(stdout, m_U) ;
    m_U->cols = i ;

    printf("V:\n") ;
    m_V->cols = 10 ;
    MatrixPrint(stdout, m_V) ;
    m_V->cols = i ;

    printf("Vhat:\n") ;
    m_tmp->cols = 10 ;
    MatrixPrint(stdout, m_tmp) ;
    m_tmp->cols = i ;


    MatrixFree(&m_tmp) ;
  }

  m_I = MatrixIdentity(m_L->rows, NULL) ;
  MatrixScalarMul(m_I, reg, m_I) ;
  MatrixScalarMul(m_L, 1-reg, m_L) ;
  MatrixAdd(m_I, m_L, m_L) ;
  MatrixFree(&m_I) ;
  if (DIAG_VERBOSE_ON)
  {
    MATRIX *m_tmp ;
    int    i ;
    m_tmp = MatrixMultiply(m_L, m_U, NULL) ;

    i = m_tmp->cols ;

    printf("U:\n") ;
    m_U->cols = 10 ;
    MatrixPrint(stdout, m_U) ;
    m_U->cols = i ;

    printf("V:\n") ;
    m_V->cols = 10 ;
    MatrixPrint(stdout, m_V) ;
    m_V->cols = i ;

    printf("Vhat:\n") ;
    m_tmp->cols = 10 ;
    MatrixPrint(stdout, m_tmp) ;
    m_tmp->cols = i ;


    MatrixFree(&m_tmp) ;
  }

  MatrixFree(&m_U) ;
  MatrixFree(&m_V) ;
  MatrixFree(&m_Uinv) ;
  return(m_L) ;
}


int
gcamApplyLinearTransform(GCA_MORPH *gcam, MATRIX *m_L)
{
  VECTOR         *v_src, *v_dst ;
  int             x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  v_src = VectorAlloc(4, MATRIX_REAL) ;
  v_dst = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_src,4) = 1.0 ;
  VECTOR_ELT(v_dst,4) = 1.0 ;
  // use gca information
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        V3_X(v_src) = gcamn->x ;
        V3_Y(v_src) = gcamn->y ;
        V3_Z(v_src) = gcamn->z ;
        MatrixMultiply(m_L, v_src, v_dst) ;
        gcamn->x = V3_X(v_dst) ;
        gcamn->y = V3_Y(v_dst) ;
        gcamn->z = V3_Z(v_dst) ;
      }
    }
  }
  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;
  return(NO_ERROR) ;
}


int
gcamShowCompressed(GCA_MORPH *gcam, FILE *fp)
{
  int n1, n2, n3, n4 ;
  n1 = GCAMcountCompressedNodes(gcam, 0.5) ;
  n2 = GCAMcountCompressedNodes(gcam, 0.25) ;
  n3 = GCAMcountCompressedNodes(gcam, 0.1) ;
  n4 = GCAMcountCompressedNodes(gcam, 0.025) ;
  if (n1 || n2 || n3 || n4)
    fprintf(fp,
            "%d nodes compressed > 0.5, %d > 0.25, %d > .1, %d > 0.025\n",
            n1, n2,n3,n4) ;
  return(n1+n2+n3+n4) ;
}

int
gcamSuppressNegativeGradients(GCA_MORPH *gcam, float scale)
{
  int            i, j, k /*, in, jn, kn, i1, j1, k1 */ ;
  GCA_MORPH_NODE *gcamn ;

  for (i = 0 ; i < gcam->width ; i++)
    for (j = 0 ; j < gcam->height ; j++)
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }
        if (gcamn->area1 < 0 || gcamn->area2 < 0)
        {
          gcamn->dx *= scale ;
          gcamn->dy *= scale ;
          gcamn->dz *= scale ;
        }
#if 0
        for (in = -1 ; in <= 1 ; in++)
        {
          i1 = i+in ;
          if (i1 < 0 || i1 >= gcam->width)
          {
            continue ;
          }
          for (jn = -1 ; jn <= 1 ; jn++)
          {
            j1 = j+jn ;
            if (j1 < 0 || j1 >= gcam->height)
            {
              continue ;
            }
            for (kn = -1 ; kn <= 1 ; kn++)
            {
              k1 = k+kn ;
              if (k1 < 0 || k1 >= gcam->depth)
              {
                continue ;
              }
              if (!in && !jn && !kn)
              {
                continue ;
              }
              gcamn = &gcam->nodes[i1][j1][k1] ;
              if (i1 == Gx && j1 == Gy && k1 == Gz)
              {
                DiagBreak() ;
              }
              if (gcamn->area1 < 0 || gcamn->area2 < 0)
              {
                gcamn->dx *= scale ;
                gcamn->dy *= scale ;
                gcamn->dz *= scale ;
              }
            }
          }
        }
#endif
      }
  return(NO_ERROR) ;
}


int
gcamWriteDiagnostics(GCA_MORPH *gcam)
{
  FILE  *fp ;
  static int dbg = 0 ;
  static int callno = 0 ;
  int    i, j, k ;
  GCA_MORPH_NODE *gcamn ;
  char   fname[STRLEN] ;

  if (getenv("DEBUG"))
  {
    dbg = 1 ;
  }
  if (dbg == 0)
  {
    return(NO_ERROR) ;
  }

  for (i = 0 ; i < gcam->width ; i++)
  {
    sprintf(fname, "v%d_%d.dat", callno, i) ;
    fp = fopen(fname, "w") ;
    for (j = 0 ; j < gcam->height ; j++)
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        fprintf(fp, "%f %f %f\n", gcamn->x, gcamn->y, gcamn->z) ;
      }

    fclose(fp) ;
  }

  callno++ ;
  return(NO_ERROR) ;
}


int
gcamCheck(GCA_MORPH *gcam, MRI *mri)
{
  GCA_MORPH_NODE *gcamn ;
  int            x, y, z, first = 1, num ;
  float          orig_area, dif1, dif2 ;

  if (getenv("GCAM_CHECK") == NULL)
  {
    return(NO_ERROR) ;
  }

  GCAMcomputeOriginalProperties(gcam) ;

  // compute mean orig area
  for (num = 0, orig_area = 0.0, x  = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (fabs(gcamn->x-120)<1 &&
            fabs(gcamn->y-125)<1 &&
            fabs(gcamn->z-123)<1)
        {
          DiagBreak() ;
        }
        if (gcamn->label == 0 ||
            gcamn->status & GCAM_BINARY_ZERO ||
            gcamn->invalid == GCAM_AREA_INVALID)
        {
          continue ;
        }
        num += 2;
        orig_area += (gcamn->orig_area1 + gcamn->orig_area2) ;
      }
  if (num == 0)
  {
    return(NO_ERROR) ;
  }
  orig_area /= (float)num ;
  for (x  = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == 0 || gcamn->status & GCAM_BINARY_ZERO)
        {
          continue ;
        }
        dif1 = fabs(gcamn->orig_area1-orig_area) ;
        dif2 = fabs(gcamn->orig_area2-orig_area) ;
        if ((dif1 > (0.01*orig_area)) || (dif2 > (0.01*orig_area)))
        {
          if (first)
          {
            printf("node(%d, %d, %d): "
                   "orig areas (%2.3f, %2.3f, %2.3f), should be %2.3f\n",
                   x, y, z,
                   gcamn->orig_area1,
                   gcamn->orig_area2,
                   gcamn->orig_area, orig_area) ;
            first = 0 ;
          }
          DiagBreak() ;
        }
      }
  return(NO_ERROR)  ;
}



GCA_MORPH *
GCAMupsample2(GCA_MORPH *gcam)
{
  int            xn, yn, zn, xo, yo, zo ;
  GCA_MORPH      *gcam_new ;
  GCA_MORPH_NODE *gcamn_old, *gcamn_new ;
  double           xd, yd, zd ;
  MRI            *mri_origx, *mri_origy, *mri_origz, *mri_x, *mri_y, *mri_z ;

  mri_origx = GCAMwriteMRI(gcam, NULL, GCAM_ORIGX) ;
  mri_origy = GCAMwriteMRI(gcam, NULL, GCAM_ORIGY) ;
  mri_origz = GCAMwriteMRI(gcam, NULL, GCAM_ORIGZ) ;
  mri_x = GCAMwriteMRI(gcam, NULL, GCAM_X) ;
  mri_y = GCAMwriteMRI(gcam, NULL, GCAM_Y) ;
  mri_z = GCAMwriteMRI(gcam, NULL, GCAM_Z) ;

  gcam_new = GCAMalloc(2*gcam->width, 2*gcam->height, 2*gcam->depth) ;
  *(&gcam_new->image) = *(&gcam->image) ;
  *(&gcam_new->atlas) = *(&gcam->atlas) ;
  gcam_new->ninputs = gcam->ninputs ;
  gcam_new->type = gcam->type ;
  gcam_new->status = gcam->status ;
  gcam_new->exp_k = gcam->exp_k ;
  gcam_new->spacing = gcam->spacing ;

  for (xo = 0 ; xo < gcam->width ; xo++)
  {
    for (yo = 0 ; yo < gcam->height ; yo++)
    {
      for (zo = 0 ; zo < gcam->depth ; zo++)
      {
        gcamn_old = &gcam->nodes[xo][yo][zo] ;
        for (xn = 2*xo ; xn <= 2*xo+1 ; xn++)
        {
          for (yn = 2*yo ; yn <= 2*yo+1 ; yn++)
          {
            for (zn = 2*zo ; zn <= 2*zo+1 ; zn++)
            {
              gcamn_new = &gcam_new->nodes[xn][yn][zn] ;

              // get morphed coords interpolated
              MRIsampleVolume(mri_x,
                              (float)xn/2.0,
                              (float)yn/2.0,
                              (float)zn/2.0, &xd) ;
              MRIsampleVolume(mri_y,
                              (float)xn/2.0,
                              (float)yn/2.0,
                              (float)zn/2.0, &yd) ;
              MRIsampleVolume(mri_z,
                              (float)xn/2.0,
                              (float)yn/2.0,
                              (float)zn/2.0, &zd) ;
              gcamn_new->x = xd ;
              gcamn_new->y = yd ;
              gcamn_new->z = zd ;

              // get orig coords interpolated
              MRIsampleVolume(mri_origx,
                              (float)xn/2.0,
                              (float)yn/2.0,
                              (float)zn/2.0, &xd) ;
              MRIsampleVolume(mri_origy,
                              (float)xn/2.0,
                              (float)yn/2.0,
                              (float)zn/2.0, &yd) ;
              MRIsampleVolume(mri_origz,
                              (float)xn/2.0,
                              (float)yn/2.0,
                              (float)zn/2.0, &zd) ;
              gcamn_new->origx = xd ;
              gcamn_new->origy = yd ;
              gcamn_new->origz = zd ;

              gcamn_new->xn = xd ;
              gcamn_new->yn = yd ;
              gcamn_new->zn = zd ;
              gcamn_new->label = gcamn_old->label ;
              gcamn_new->n = gcamn_old->n ;
              gcamn_new->prior = gcamn_old->prior ;
              if (gcamn_old->gc)
              {
                gcamn_new->gc = alloc_gcs(1, GCA_NO_MRF, gcam->ninputs) ;
                copy_gcs(1, gcamn_old->gc, gcamn_new->gc, gcam->ninputs) ;
              }
              gcamn_new->log_p = gcamn_old->log_p ;
              gcamn_new->orig_area = gcamn_old->orig_area / (2*2*2) ;
              gcamn_new->orig_area1 = gcamn_old->orig_area1 / (2*2*2) ;
              gcamn_new->orig_area2 = gcamn_old->orig_area2 / (2*2*2) ;
              gcamn_new->status = gcamn_old->status ;
              gcamn_new->invalid = gcamn_old->invalid ;
              gcamn_new->label_dist = gcamn_old->label_dist ;
            }
          }
        }
      }
    }
  }
  gcam_new->atlas.xsize /= 2 ;
  gcam_new->atlas.ysize /= 2 ;
  gcam_new->atlas.zsize /= 2 ;
  gcam_new->atlas.width *= 2 ;
  gcam_new->atlas.height *= 2 ;
  gcam_new->atlas.depth *= 2 ;

  GCAMfreeContents(gcam) ;
  *gcam = *gcam_new ;
  free(gcam_new) ;
  GCAMcomputeOriginalProperties(gcam) ;
  gcamComputeMetricProperties(gcam) ;
  return(NO_ERROR) ;
}

int
GCAMignoreZero(GCA_MORPH *gcam, MRI *mri_target)
{
  int             x, y, z ;
  double            val ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolume(mri_target, gcamn->x, gcamn->y, gcamn->z, &val) ;
        if (gcamn->label > 0 && !FZERO(val))
        {
          continue ;
        }
        if ((gcamn->gc == NULL || FZERO(gcamn->gc->means[0]) || FZERO(val)))
        {
          gcamn->status = gcamn->status | GCAM_NEVER_USE_LIKELIHOOD;
        }
      }
    }
  }

  return(NO_ERROR) ;
}


int
GCAMremoveIgnoreZero(GCA_MORPH *gcam, MRI *mri_target)
{
  int             x, y, z ;
  double            val ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;

        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolume(mri_target, gcamn->x, gcamn->y, gcamn->z, &val) ;
        if (gcamn->label > 0 && !FZERO(val))
        {
          continue ;
        }
        if ((gcamn->gc == NULL || FZERO(gcamn->gc->means[0]) || FZERO(val)))
          gcamn->status &= ~GCAM_NEVER_USE_LIKELIHOOD;
      }
    }
  }

  return(NO_ERROR) ;
}


MATRIX *
gcamComputeOptimalLinearTransformInRegion(GCA_MORPH *gcam, MRI *mri_mask, MATRIX *m_L, int x0, int y0, int z0, int whalf)
{
  int            ncols, col, xi, yi, zi, xk, yk, zk ;
  GCA_MORPH_NODE *gcamn ;
  MATRIX         *m_U, *m_V, *m_Uinv ;
  double            mask ;


  ncols = 0 ;    // first count the # of measurements
  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = x0+xk ;
    if (xi < 0 || xi >= gcam->width)
    {
      continue ;
    }
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = y0+yk ;
      if (yi < 0 || yi >= gcam->height)
      {
        continue ;
      }
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = z0+zk ;
        if (zi < 0 || zi >= gcam->depth)
        {
          continue ;
        }
        gcamn = &gcam->nodes[xi][yi][zi] ;
        MRIsampleVolume(mri_mask, gcamn->x, gcamn->y, gcamn->z, &mask) ;
        if (!FZERO(mask) && mask > 0)
        {
          ncols++ ;
        }
      }
    }
  }

  m_U = MatrixAlloc(4, ncols, MATRIX_REAL) ;
  m_V = MatrixAlloc(4, ncols, MATRIX_REAL) ;
  for (col = 0, xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = x0+xk ;
    if (xi < 0 || xi >= gcam->width)
    {
      continue ;
    }
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = y0+yk ;
      if (yi < 0 || yi >= gcam->height)
      {
        continue ;
      }
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = z0+zk ;
        if (zi < 0 || zi >= gcam->depth)
        {
          continue ;
        }
        gcamn = &gcam->nodes[xi][yi][zi] ;
        MRIsampleVolume(mri_mask, gcamn->x, gcamn->y, gcamn->z, &mask) ;
        if (!FZERO(mask) && mask > 0)
        {
          col++ ;
          *MATRIX_RELT(m_U,1,col) = xi ;
          *MATRIX_RELT(m_U,2,col) = yi ;
          *MATRIX_RELT(m_U,3,col) = zi ;
          *MATRIX_RELT(m_U,4,col) = 1.0 ;

          *MATRIX_RELT(m_V,1,col) = gcamn->x ;
          *MATRIX_RELT(m_V,2,col) = gcamn->y ;
          *MATRIX_RELT(m_V,3,col) = gcamn->z ;
          *MATRIX_RELT(m_V,4,col) = 1.0 ;
        }
      }
    }
  }

  // m_Uinv = MatrixPseudoInverse(m_U, NULL) ;
  m_Uinv = MatrixSVDPseudoInverse(m_U, NULL) ;
  // MatrixPrint(stdout, m_Uinv) ;
  m_L = MatrixMultiply(m_V, m_Uinv, m_L) ;
  *MATRIX_RELT(m_L, 4, 1) = 0.0 ;
  *MATRIX_RELT(m_L, 4, 2) = 0.0 ;
  *MATRIX_RELT(m_L, 4, 3) = 0.0 ;
  *MATRIX_RELT(m_L, 4, 4) = 1.0 ;

  if (DIAG_VERBOSE_ON)
  {
    MATRIX *m_tmp ;
    int    i ;
    printf("m_L:\n") ;
    MatrixPrint(stdout, m_L) ;
    m_tmp = MatrixMultiply(m_L, m_U, NULL) ;

    i = m_tmp->cols ;

    printf("U:\n") ;
    m_U->cols = 10 ;
    MatrixPrint(stdout, m_U) ;
    m_U->cols = i ;

    printf("V:\n") ;
    m_V->cols = 10 ;
    MatrixPrint(stdout, m_V) ;
    m_V->cols = i ;

    printf("Vhat:\n") ;
    m_tmp->cols = 10 ;
    MatrixPrint(stdout, m_tmp) ;
    m_tmp->cols = i ;


    MatrixFree(&m_tmp) ;
  }

#if 0
  // regularize it
  m_I = MatrixIdentity(m_L->rows, NULL) ;
  MatrixScalarMul(m_I, reg, m_I) ;
  MatrixScalarMul(m_L, 1-reg, m_L) ;
  MatrixAdd(m_I, m_L, m_L) ;
  MatrixFree(&m_I) ;
#endif

  if (DIAG_VERBOSE_ON)
  {
    MATRIX *m_tmp ;
    int    i ;
    m_tmp = MatrixMultiply(m_L, m_U, NULL) ;

    i = m_tmp->cols ;

    printf("U:\n") ;
    m_U->cols = 10 ;
    MatrixPrint(stdout, m_U) ;
    m_U->cols = i ;

    printf("V:\n") ;
    m_V->cols = 10 ;
    MatrixPrint(stdout, m_V) ;
    m_V->cols = i ;

    printf("Vhat:\n") ;
    m_tmp->cols = 10 ;
    MatrixPrint(stdout, m_tmp) ;
    m_tmp->cols = i ;


    MatrixFree(&m_tmp) ;
  }

  MatrixFree(&m_U) ;
  MatrixFree(&m_V) ;
  MatrixFree(&m_Uinv) ;
  return(m_L) ;
}

int
gcamExpansionTerm(GCA_MORPH *gcam, MRI *mri, double l_expansion)
{
  double    sse, error_0, error_n, dx, dy, dz ;
  int       x, y, z, xi, yi, zi, xk, yk, zk ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;
  float            vals[MAX_GCA_INPUTS] ;

  if (FZERO(l_expansion))
  {
    return(0.0) ;
  }

  sse = 0.0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID || gcamn->gc == NULL)
        {
          continue;
        }
        check_gcam(gcam) ;
        if (gcamn->status & (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD))
        {
          continue ;
        }
        if (IS_UNKNOWN(gcamn->label))
        {
          continue ;
        }
#if 0
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals,gcam->ninputs);
        error_0 = sqrt(GCAmahDist(gcamn->gc, vals, gcam->ninputs)) ;
        if (error_0 > 2)
        {
          continue ;  // doesn't explain intensity well
        }
#endif
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = x+xk ;
          if (xi < 0 || xi >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = y+yk ;
            if (yi < 0 || yi >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = z+zk ;
              if (xi == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak() ;
              }
              if (zi < 0 || zi >= gcam->height)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xi][yi][zi] ;
              if (gcamn_nbr->label == gcamn->label)
              {
                continue ;
              }

              load_vals(mri, gcamn_nbr->x, gcamn_nbr->y, gcamn_nbr->z,
                        vals,gcam->ninputs);
              if (gcamn_nbr->gc == NULL)
              {
                continue ;
              }
              error_0 =
                sqrt(GCAmahDist(gcamn->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn->gc, gcam->ninputs)) ;
              error_n =
                sqrt(GCAmahDist(gcamn_nbr->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn_nbr->gc, gcam->ninputs)) ;
              if (error_n <= error_0)
              {
                continue ;
              }
              dx = gcamn_nbr->x - gcamn->x ;
              dy = gcamn_nbr->y - gcamn->y ;
              dz = gcamn_nbr->z - gcamn->z ;
              dx *= l_expansion * (error_n-error_0) ;
              dy *= l_expansion * (error_n-error_0) ;
              dz *= l_expansion * (error_n-error_0) ;
              if (x == Gx && y == Gy && z == Gz)
                printf("(%d, %d, %d, %s): nbr at (%d, %d, %d, %s), I=%2.0f, D=(%2.2f, %2.2f, %2.2f)\n",
                       x, y, z, cma_label_to_name(gcamn->label),
                       xi, yi, zi, cma_label_to_name(gcamn_nbr->label),
                       vals[0], dx, dy, dz) ;
              gcamn->dx += dx ;
              gcamn->dy += dy ;
              gcamn->dz += dz ;


            }
          }
        }

        // now look for nbrs that explain central vals better than central label
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals,gcam->ninputs);
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = x+xk ;
          if (xi < 0 || xi >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = y+yk ;
            if (yi < 0 || yi >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = z+zk ;
              if (xi == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak() ;
              }
              if (zi < 0 || zi >= gcam->height)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xi][yi][zi] ;
              if (gcamn_nbr->label == gcamn->label)
              {
                continue ;
              }

              if (gcamn_nbr->gc == NULL)
              {
                continue ;
              }
              error_0 =
                sqrt(GCAmahDist(gcamn->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn->gc, gcam->ninputs)) ;
              error_n =
                sqrt(GCAmahDist(gcamn_nbr->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn_nbr->gc, gcam->ninputs)) ;
              if (error_0 <= error_n)
              {
                continue ;
              }

              // move away from nbr if it explains central intensity better than this node
              dx = gcamn->x - gcamn_nbr->x  ;
              dy = gcamn->y - gcamn_nbr->y  ;
              dz = gcamn->z - gcamn_nbr->z  ;
              dx *= l_expansion * (error_0-error_n) ;
              dy *= l_expansion * (error_0-error_n) ;
              dz *= l_expansion * (error_0-error_n) ;
              if (x == Gx && y == Gy && z == Gz)
                printf("(%d, %d, %d, %s): Inbr at (%d, %d, %d, %s), I=%2.0f, D=(%2.2f, %2.2f, %2.2f)\n",
                       x, y, z, cma_label_to_name(gcamn->label),
                       xi, yi, zi, cma_label_to_name(gcamn_nbr->label),
                       vals[0], dx, dy, dz) ;
              gcamn->dx += dx ;
              gcamn->dy += dy ;
              gcamn->dz += dz ;


            }
          }
        }

      }
  return(NO_ERROR) ;
}


double
gcamExpansionEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double    sse, error_0, error_n, sse_node ;
  int       x, y, z, xi, yi, zi, xk, yk, zk ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;
  float            vals[MAX_GCA_INPUTS] ;

  sse = 0.0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID || gcamn->gc == NULL)
        {
          continue;
        }
        check_gcam(gcam) ;
        if (gcamn->status & (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD))
        {
          continue ;
        }
        if (IS_UNKNOWN(gcamn->label))
        {
          continue ;
        }
#if 0
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals,gcam->ninputs);
        error_0 = sqrt(GCAmahDist(gcamn->gc, vals, gcam->ninputs)) ;
        if (error_0 > 2)
        {
          continue ;  // doesn't explain intensity well
        }
#endif

        sse_node = 0.0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = x+xk ;
          if (xi < 0 || xi >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = y+yk ;
            if (yi < 0 || yi >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = z+zk ;
              if (xi == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak() ;
              }
              if (zi < 0 || zi >= gcam->height)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xi][yi][zi] ;
              if (gcamn_nbr->label == gcamn->label)
              {
                continue ;
              }

              load_vals(mri, gcamn_nbr->x, gcamn_nbr->y, gcamn_nbr->z,
                        vals,gcam->ninputs);
              if (gcamn_nbr->gc == NULL)
              {
                continue ;
              }
              error_0 =
                sqrt(GCAmahDist(gcamn->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn->gc, gcam->ninputs)) ;
              error_n =
                sqrt(GCAmahDist(gcamn_nbr->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn_nbr->gc, gcam->ninputs)) ;
              if (error_n <= error_0)
              {
                continue ;
              }
              sse_node += error_n*error_n - error_0*error_0  ;
              if (!finite(sse))
              {
                DiagBreak() ;
              }
            }
          }
        }

        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals,gcam->ninputs);
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = x+xk ;
          if (xi < 0 || xi >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = y+yk ;
            if (yi < 0 || yi >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = z+zk ;
              if (xi == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak() ;
              }
              if (zi < 0 || zi >= gcam->height)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xi][yi][zi] ;
              if (gcamn_nbr->label == gcamn->label)
              {
                continue ;
              }

              if (gcamn_nbr->gc == NULL)
              {
                continue ;
              }
              error_0 =
                sqrt(GCAmahDist(gcamn->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn->gc, gcam->ninputs)) ;
              error_n =
                sqrt(GCAmahDist(gcamn_nbr->gc, vals, gcam->ninputs)) + \
                log(covariance_determinant(gcamn_nbr->gc, gcam->ninputs)) ;
              if (error_0 <= error_n)  // check to make sure nbr doesn't explain my vals better
              {
                continue ;
              }
              sse_node += error_0*error_0 - error_n*error_n  ;
              if (!finite(sse))
              {
                DiagBreak() ;
              }
            }
          }
        }
        if (sse_node < 0 || !finite(sse_node))
        {
          DiagBreak() ;
        }
        sse += sse_node ;
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("E_exp(%d, %d, %d): %2.1f\n", x, y, z,sse_node) ;
        }
      }

  return(sse) ;
}



int
GCAMmatchVentricles(GCA_MORPH *gcam, MRI *mri)
{
  int             label, ven_nbr, x, y, z, xi, yi, zi, xk, yk, zk,
                  nchanged = 0, n ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;
  float            vals[MAX_GCA_INPUTS] ;
  GCA_PRIOR      *gcap ;
  GC1D           *gc ;

  // look for regions around atlas ventricles that look like ventricle
  // and let it expand to provide a better initialization
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID || gcamn->gc == NULL)
        {
          continue;
        }
        if (IS_UNKNOWN(gcamn->label) ||
            (gcamn->label == Left_Lateral_Ventricle) ||
            (gcamn->label == Right_Lateral_Ventricle))
        {
          continue ;  // already ventricle, don't change it
        }
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals,gcam->ninputs);
        label =
          GCAcomputeMLElabelAtLocation(gcam->gca, x,y,z,vals,&n,&gcamn->log_p);
        if (label != Left_Lateral_Ventricle &&
            label != Right_Lateral_Ventricle)
        {
          continue ;
        }
        ven_nbr = 0 ;
        for (xk = -1 ; ven_nbr == 0 && xk <= 1 ; xk++)
        {
          xi = x+xk ;
          if (xi < 0 || xi >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; ven_nbr == 0 && yk <= 1 ; yk++)
          {
            yi = y+yk ;
            if (yi < 0 || yi >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; ven_nbr == 0 && zk <= 1 ; zk++)
            {
              zi = z+zk ;
              if (xi == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak() ;
              }
              if (zi < 0 || zi >= gcam->height)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xi][yi][zi] ;

              if (gcamn_nbr->gc == NULL)
              {
                continue ;
              }
              if (gcamn_nbr->label == Left_Lateral_Ventricle ||
                  gcamn_nbr->label == Right_Lateral_Ventricle)
              {
                ven_nbr = 1 ;
              }
            }
          }
        }

        if (ven_nbr == 0)
        {
          continue ;
        }
        if (n >= 0)
        {
          if (label != gcamn->label)
          {
            nchanged++ ;
          }
          gcap = &gcam->gca->priors[x][y][z] ;
          gcamn->label = label ;
          gcamn->n = n ;
          gcamn->prior = gcap->priors[n] ;
          gcamn->gc = gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;

          if (x == Gx && y == Gy && z == Gz)
            printf("RELABEL: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f, MRI=%2.0f\n",
                   x, y, z, cma_label_to_name(label), label,
                   gcamn->gc ? gcamn->gc->means[0] : 0.0,
                   gcamn->gc ? sqrt(covariance_determinant(gcamn->gc, gcam->ninputs)) : 0.0,
                   gcamn->prior,vals[0]) ;
        }
        else  /* out of FOV probably */
        {
          gcamn->label = label ;
          gcamn->n = 0 ;
          gcamn->prior = 1.0 ;
          gcamn->gc = NULL ;
          if (x == Gx && y == Gy && z == Gz)
            printf("RELABEL: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f\n",
                   x, y, z, cma_label_to_name(label), label,
                   0.0, 0.0, gcamn->prior) ;
        }
      }

  printf("%d nodes changed to ventricle\n", nchanged) ;
  return(NO_ERROR) ;
}

int
GCAMdeformVentricles(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  double      rms, last_rms, orig_dt, min_dt, pct_change ;
  int         navgs, level ;

  parms->mri = mri ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.log", parms->base_name) ;
    if (parms->log_fp == NULL)
    {
      if (parms->start_t == 0)
      {
        parms->log_fp = fopen(fname, "w") ;
      }
      else
      {
        parms->log_fp = fopen(fname, "a") ;
      }
    }
  }
  else
  {
    parms->log_fp = NULL ;
  }
  orig_dt = parms->dt ;
  if (parms->start_t == 0)
  {
    if ((parms->write_iterations > 0) && (Gdiag & DIAG_WRITE))
    {
      write_snapshot(gcam, mri, parms, parms->start_t) ;
    }
    parms->start_t++ ;
  }
  for (navgs = parms->navgs, level = parms->levels ; level >= 0 ; level--)
  {
    do
    {
      gcamClearGradient(gcam) ;
      gcamComputeMetricProperties(gcam) ;
      last_rms = GCAMcomputeRMS(gcam, mri, parms) ;
      gcamComputePeriventricularWMDeformation(gcam, mri) ;
      gcamSmoothnessTerm(gcam, mri, parms->l_smoothness)  ;
#if 0
      gcamLikelihoodTerm(gcam, mri, mri, parms->l_likelihood, parms)  ;
      gcamLogLikelihoodTerm(gcam, mri, mri, parms->l_log_likelihood)  ;
#endif
      gcamJacobianTerm(gcam, mri, parms->l_jacobian,parms->ratio_thresh)  ;
      gcamLimitGradientMagnitude(gcam, parms, mri) ;
      gcamSmoothGradient(gcam, navgs) ;
      parms->dt = orig_dt ;
      min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
      parms->dt = min_dt ;
      gcamApplyGradient(gcam, parms) ;
      if (Gdiag & DIAG_SHOW)
      {
        gcamShowCompressed(gcam, stdout) ;
      }
      if (parms->uncompress)
      {
        gcamRemoveCompressedNodes(gcam, mri, parms, parms->ratio_thresh) ;
      }
      if ((parms->write_iterations > 0) && (Gdiag & DIAG_WRITE))
      {
        write_snapshot(gcam, mri, parms, parms->start_t) ;
      }
      rms = GCAMcomputeRMS(gcam, mri, parms) ;
      pct_change = 100.0*(last_rms-rms)/last_rms ;
      last_rms = rms ;
      if (parms->log_fp)
      {
        fprintf
        (parms->log_fp,
         "%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, "
         "invalid=%d, VD, navgs=%d",
         parms->start_t, min_dt, rms, pct_change, gcam->neg,
         Ginvalid,navgs) ;
        if (parms->l_binary > 0)
          fprintf
          (parms->log_fp, ", aligned = %d (%2.3f%%)\n",
           Galigned,
           100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          fprintf(parms->log_fp, "\n") ;
        }
        fflush(parms->log_fp) ;
      }

      if (Gdiag & DIAG_SHOW)
      {
        printf
        ("%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), "
         "neg=%d, invalid=%d, VD, navgs=%d",
         parms->start_t, min_dt, rms, pct_change,
         gcam->neg, Ginvalid, navgs) ;
        if (parms->l_binary > 0)
          printf
          (", aligned = %d (%2.3f%%)\n",
           Galigned,
           100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          printf("\n") ;
        }
      }
      gcamCheck(gcam, mri) ;
      parms->start_t++ ;
    }
    while (pct_change > parms->tol) ;
    navgs /= 4 ;
    if (navgs < 1)
    {
      break ;
    }
  }
  parms->dt = orig_dt ;

  return(NO_ERROR) ;
}


int
gcamComputePeriventricularWMDeformation(GCA_MORPH *gcam, MRI *mri)
{
  int    ven_nbr, x, y, z, xi, yi, zi, xk, yk, zk ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;
  float  vals[MAX_GCA_INPUTS] ;
  double norm, xv, yv, zv, dx, dy, dz, d, mah_dist, min_mah_dist, min_d;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID || gcamn->gc == NULL)
        {
          continue;
        }
#if 0
        if (gcamn->label != Left_Cerebral_White_Matter &&
            gcamn->label != Right_Cerebral_White_Matter)
        {
          continue ;
        }
#endif
        if (gcamn->label == Left_Lateral_Ventricle ||
            gcamn->label == Right_Lateral_Ventricle)
        {
          continue ;
        }

        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals,gcam->ninputs);
        ven_nbr = 0 ;
        dx = dy = dz = 0 ;
        for (xk = -1 ; ven_nbr == 0 && xk <= 1 ; xk++)
        {
          xi = x+xk ;
          if (xi < 0 || xi >= gcam->width)
          {
            continue ;
          }
          for (yk = -1 ; ven_nbr == 0 && yk <= 1 ; yk++)
          {
            yi = y+yk ;
            if (yi < 0 || yi >= gcam->height)
            {
              continue ;
            }
            for (zk = -1 ; ven_nbr == 0 && zk <= 1 ; zk++)
            {
              zi = z+zk ;
              if (xi == Gx && yi == Gy && zi == Gz)
              {
                DiagBreak() ;
              }
              if (zi < 0 || zi >= gcam->height)
              {
                continue ;
              }
              gcamn_nbr = &gcam->nodes[xi][yi][zi] ;

              if (gcamn_nbr->gc == NULL)
              {
                continue ;
              }
              if (gcamn_nbr->label == Left_Lateral_Ventricle ||
                  gcamn_nbr->label == Right_Lateral_Ventricle)
              {
                ven_nbr++ ;
                dx += (gcamn->x - gcamn_nbr->x) ;
                dy += (gcamn->y - gcamn_nbr->y) ;
                dz += (gcamn->z - gcamn_nbr->z) ;
              }
            }
          }
        }

        if (ven_nbr == 0)
        {
          continue ;
        }

        dx /= ven_nbr ;
        dy /= ven_nbr ;
        dz /= ven_nbr ;
        norm = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (FZERO(norm))
        {
          continue ;
        }

        dx /= norm ;
        dy /= norm ;
        dz /= norm ;
        if (Gx == x && Gy == y && Gz == z)
          printf("searching node (%d, %d, %d) along "
                 "direction (%2.1f, %2.1f, %2.1f)\n",
                 Gx, Gy, Gz, dx, dy, dz) ;

        // now move away from avg ventricle position until we get
        // to an intensity that looks like wm
#define MAX_WM_DISTANCE 10
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals,gcam->ninputs);
        min_mah_dist = GCAmahDist(gcamn->gc, vals, gcam->ninputs) ;
        min_d = 0 ;
        for (d = -1.0 ; d < MAX_WM_DISTANCE ; d += 0.5)
        {
          xv = gcamn->x + d*dx ;
          yv = gcamn->y + d*dy ;
          zv = gcamn->z + d*dz ;
          load_vals(mri, xv, yv, zv, vals,gcam->ninputs);
          mah_dist = GCAmahDist(gcamn->gc, vals, gcam->ninputs) ;
          if (mah_dist < min_mah_dist)
          {
            min_mah_dist = mah_dist ;
            min_d = d ;
          }
        }
        if (Gx == x && Gy == y && Gz == z)
          printf("MLE location at min_d=%2.1f, X = (%2.0f, %2.0f, %2.0f)\n",
                 min_d,
                 gcamn->x+min_d*dx,
                 gcamn->y+min_d*dy,
                 gcamn->z+min_d*dz) ;
        gcamn->dx += min_d*dx ;
        gcamn->dy += min_d*dy ;
        gcamn->dz += min_d*dz ;
      }
  return(NO_ERROR) ;
}

HISTOGRAM *
gcamJacobianHistogram(GCA_MORPH *gcam, HISTOGRAM *h)
{
  int            i, j, k, bin_no, nbins ;
  GCA_MORPH_NODE *gcamn ;
  double         jmin, jmax, jac, orig_area ;
  static       double bin_size = -1 ;

  orig_area = gcam->spacing*gcam->spacing*gcam->spacing ;

  // compute bounds on jacobian
  jmin = 10000000 ;
  jmax = -jmin ;
  for (i = 0 ; i < gcam->width ; i++)
    for (j = 0 ; j < gcam->height ; j++)
      for (k = 0 ; k < gcam->depth ; k++)
      {
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid || gcamn->gc == NULL ||
            IS_UNKNOWN(gcamn->label) ||
            FZERO(gcamn->area) || FZERO(gcamn->orig_area1) ||
            FZERO(gcamn->orig_area2))
        {
          continue;
        }
        jac = (gcamn->area / orig_area) ;
        if (jac > jmax)
        {
          jmax = jac ;
        }
        if (jac < jmin)
        {
          jmin = jac ;
        }
      }

  if (bin_size <= 0)
  {
    nbins = 512 ;
    bin_size = (jmax - jmin + 1) / (float)nbins ;
    printf("setting bin_size to %2.2f\n", bin_size) ;
  }
  else
  {
    nbins = ceil((jmax - jmin + 1) / bin_size) ;
    printf("using %d bins\n", nbins) ;
  }
  if (h == NULL)
  {
    h = HISTOalloc(nbins) ;
  }
  else
  {
    HISTOrealloc(h, nbins) ;
  }


  for (bin_no = 0 ; bin_no < nbins ; bin_no++)
  {
    h->bins[bin_no] = (bin_no)*bin_size+jmin ;
  }
  for (i = 0 ; i < gcam->width ; i++)
    for (j = 0 ; j < gcam->height ; j++)
      for (k = 0 ; k < gcam->depth ; k++)
      {
        if (i == Gx && j == Gy && k == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[i][j][k] ;

        if (gcamn->invalid || gcamn->gc == NULL ||
            IS_UNKNOWN(gcamn->label) ||
            FZERO(gcamn->area) || FZERO(gcamn->orig_area1) ||
            FZERO(gcamn->orig_area2))
        {
          continue;
        }
        jac = (gcamn->area / orig_area) ;
        bin_no = nint((float)(jac - jmin) / (float)bin_size) ;
        h->counts[bin_no]++ ;
      }
  return(h) ;
}
int
GCAMnormalizeIntensities(GCA_MORPH *gcam, MRI *mri_target)
{
  int             x, y, z, num, peak ;
  double            val, mean, std, low_thresh, hi_thresh ;
  GCA_MORPH_NODE  *gcamn ;
  HISTOGRAM        *h ;

  h = HISTOalloc(200) ;
  HISTOinit(h, 200, 0, 5) ;

  std = mean = 0.0 ;
  num = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->gc == NULL || FZERO(gcamn->gc->means[0]))
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        MRIsampleVolume(mri_target, gcamn->x, gcamn->y, gcamn->z, &val) ;
        if (FZERO(val))
        {
          continue ;
        }
        if (mri_target->depth > 5 && MRIlabelsInNbhd(mri_target, nint(gcamn->x), nint(gcamn->y),
                            nint(gcamn->z),1, 0) > 0)
        {
          continue ;  // stay away from skull stripped areas
        }
        val /= gcamn->gc->means[0] ;
	HISTOaddSample(h, val, 0, 0) ;
        mean += val ;
        std += val*val ;
        num++ ;
      }
    }
  }

  if (Gdiag & DIAG_WRITE)
    HISTOplot(h, "hratio.plt") ;

  if ( !num )
    printf(" GCAMnormalizeIntensities - no valid "
           "gc->means information found.\n");
  mean /= (float)num ;
  std = sqrt(std / (float)num - mean*mean) ;
  low_thresh = mean-2*std ;
  hi_thresh = mean+2*std ;
  printf("trimming scaling estimates to be in [%2.4f %2.4f]\n", low_thresh,
         hi_thresh) ;

  // now trim means to get more accurate estimate (ignoring outliers)
  mean = 0.0 ;
  num = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->gc == NULL || FZERO(gcamn->gc->means[0]))
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        MRIsampleVolume(mri_target, gcamn->x, gcamn->y, gcamn->z, &val) ;
        if (FZERO(val))
        {
          continue ;
        }
        if (mri_target->depth > 5 && MRIlabelsInNbhd(mri_target, nint(gcamn->x), nint(gcamn->y),
                            nint(gcamn->z),1, 0) > 0)
        {
          continue ;  // stay away from skull stripped areas
        }
        val /= gcamn->gc->means[0] ;
        if (val < low_thresh || val > hi_thresh)
        {
          continue ;
        }
        mean += val ;
        std += val*val ;
        num++ ;
      }
    }
  }

  mean /= (float)num ;
  printf("renormalizing gcam intensities by %2.4f\n", mean) ;

  peak = HISTOfindHighestPeakInRegion(h, 2, h->nbins-2) ;
  if (peak >= 0)
    mean = h->bins[peak] ;
  HISTOfree(&h) ;
  printf("renormalizing gcam intensities using histo estimate %2.4f\n", mean) ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->gc == NULL || FZERO(gcamn->gc->means[0]))
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn->gc->means[0] *= mean ;
      }
    }
  }
  return(NO_ERROR) ;
}

int
GCAMsmoothConditionalDensities(GCA_MORPH *gcam, float sigma)
{
  MRI *mri, *mri_kernel, *mri_smooth ;

  mri = GCAMwriteMRI(gcam, NULL, GCAM_MEANS) ;
  mri_kernel = MRIgaussian1d(sigma, -1) ;
  mri_smooth = MRIconvolveGaussian(mri, NULL, mri_kernel) ;
  gcamReadMRI(gcam, mri_smooth, GCAM_MEANS) ;
  MRIfree(&mri) ;
  MRIfree(&mri_kernel) ;
  MRIfree(&mri_smooth) ;
  return(NO_ERROR) ;
}


int
dtrans_label_to_frame(GCA_MORPH_PARMS *mp, int label)
{
  int frame ;

#if 1
  for (frame = 0 ; frame < mp->ndtrans ; frame++)
    if (mp->dtrans_labels[frame] == label)
    {
      return(frame) ;
    }
  return(-1) ;
#else
  switch (label)
  {
  case Left_Thalamus_Proper:
    frame = 0 ;
    break ;
  case Right_Thalamus_Proper:
    frame = 1 ;
    break ;
  case Left_Putamen:
    frame = 2 ;
    break ;
  case Right_Putamen:
    frame = 3 ;
    break ;
  case Left_Pallidum:
    frame = 4 ;
    break ;
  case Right_Pallidum:
    frame = 5 ;
    break ;
  case Left_Lateral_Ventricle:
    frame = 6 ;
    break ;
  case Right_Lateral_Ventricle:
    frame = 7 ;
    break ;
  case Left_Caudate:
    frame = 8 ;
    break ;
  case Right_Caudate:
    frame = 9 ;
    break ;
  case Left_Cerebral_White_Matter:
    frame = 10 ;
    break ;
  case Right_Cerebral_White_Matter:
    frame = 11 ;
    break ;
  case Left_Hippocampus:
    frame = 12 ;
    break ;
  case Right_Hippocampus:
    frame = 13 ;
    break ;
  case Left_Amygdala:
    frame = 14 ;
    break ;
  case Right_Amygdala:
    frame = 15 ;
    break ;
  case Left_VentralDC:
    frame = 16 ;
    break ;
  case Right_VentralDC:
    frame = 17 ;
    break ;
  case Brain_Stem:
    frame = 18 ;
    break ;
  case Left_Inf_Lat_Vent:
    frame = 19 ;
    break ;
  case Right_Inf_Lat_Vent:
    frame = 20 ;
    break ;
  default:
    frame = -1 ;  // not one of the labels that has a dtrans frame
  }
#endif
  return(frame) ;
}


int
gcamDistanceTransformTerm(GCA_MORPH *gcam, MRI *mri_source, MRI *mri, double l_dtrans, GCA_MORPH_PARMS *mp)
{
  int             x, y, z, frame, label ;
  GCA_MORPH_NODE  *gcamn ;
  double          dist, dx, dy, dz, error ;

  if (DZERO(l_dtrans))
  {
    return(NO_ERROR) ;
  }
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
          continue;

        check_gcam(gcam) ;
        if (gcamn->status & (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD))
          continue ;

        frame = dtrans_label_to_frame(mp, gcamn->label) ;
        if (frame >= 0)  //  one that distance transform applies to
	{
	  MRIsampleVolumeGradientFrame(mri, gcamn->x, gcamn->y, gcamn->z, &dx, &dy, &dz, frame) ;
	  MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z, frame, SAMPLE_TRILINEAR,&dist) ;
	  error = (dist - gcamn->target_dist) ;
	  if (x == Gx && y == Gy && z == Gz)
	    printf("l_dtrans: node(%d,%d,%d, %s) -> (%2.1f,%2.1f,%2.1f), dist=%2.2f, T=%2.2f, D=(%2.1f,%2.1f,%2.1f)\n",
		   x, y, z, cma_label_to_name(gcamn->label), gcamn->x, gcamn->y, gcamn->z,
		   dist, gcamn->target_dist, -error*dx, -error*dy, -error*dz) ;
	  gcamn->dx -= l_dtrans*error*dx ;
	  gcamn->dy -= l_dtrans*error*dy ;
	  gcamn->dz -= l_dtrans*error*dz ;
	}

	label = MRIgetVoxVal(mri_source, x, y, z, 0) ;
        frame = dtrans_label_to_frame(mp, label) ;
        if (frame >= 0)   // add in gradients from voxels in source that are this label
	{
	  double target_dist ;
	  MRIsampleVolumeGradientFrame(mp->mri_atlas_dist_map, x, y, z, &dx, &dy, &dz, frame) ;
	  MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z, frame, SAMPLE_TRILINEAR,&dist) ;
	  MRIsampleVolumeFrameType(mp->mri_atlas_dist_map, x, y, z, frame, SAMPLE_TRILINEAR,&target_dist) ;
	  error = (target_dist - dist) ;
	  if (x == Gx && y == Gy && z == Gz)
	    printf("l_dtrans: node(%d,%d,%d, %s) -> (%2.1f,%2.1f,%2.1f), dist=%2.2f, T=%2.2f, D=(%2.1f,%2.1f,%2.1f)\n",
		   x, y, z, cma_label_to_name(label), gcamn->x, gcamn->y, gcamn->z,
		   dist, target_dist, -error*dx, -error*dy, -error*dz) ;

	  gcamn->dx += l_dtrans*error*dx ;
	  gcamn->dy += l_dtrans*error*dy ;
	  gcamn->dz += l_dtrans*error*dz ;
	}

        check_gcam(gcam) ;
      }
  return(NO_ERROR) ;
}


double
gcamDistanceTransformEnergy(GCA_MORPH *gcam, MRI *mri_source, MRI *mri, GCA_MORPH_PARMS *mp)
{
  double          sse = 0.0 ;
  int             x, y, z, frame, label ;
  GCA_MORPH_NODE  *gcamn ;
  double          dist ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;

        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
          continue;

        check_gcam(gcam) ;
        if (gcamn->status & (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD))
          continue ;


        frame = dtrans_label_to_frame(mp, gcamn->label) ;
        if (frame >= 0) // there is a distance transform for this label
	{
	  MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z,
				   frame, SAMPLE_TRILINEAR,&dist) ;
	  if (x == Gx && y == Gy && z == Gz)
	    printf("E_dtrans: node(%d,%d,%d, %s) -> (%2.1f,%2.1f,%2.1f), dist=%2.2f, target=%2.2f\n",
		   x, y, z, cma_label_to_name(gcamn->label),
		   gcamn->x, gcamn->y, gcamn->z, dist, gcamn->target_dist) ;
	  
	  check_gcam(gcam) ;
	  dist -= (gcamn->target_dist) ;
	  sse += (dist*dist) ;
	  if (!finitep(sse))
	    DiagBreak() ;
	}

	label = MRIgetVoxVal(mri_source, x, y, z, 0) ;
        frame = dtrans_label_to_frame(mp, label) ;
        if (frame >= 0)   // add in gradients from voxels in source that are this label
	{
	  double target_dist ;
	  MRIsampleVolumeFrameType(mri, gcamn->x, gcamn->y, gcamn->z, frame, SAMPLE_TRILINEAR,&dist) ;
	  MRIsampleVolumeFrameType(mp->mri_atlas_dist_map, x, y, z, frame, SAMPLE_TRILINEAR,&target_dist) ;

	  if (x == Gx && y == Gy && z == Gz)
	    printf("E_dtrans: node(%d,%d,%d, %s) -> (%2.1f,%2.1f,%2.1f), dist=%2.2f, target=%2.2f\n",
		   x, y, z, cma_label_to_name(gcamn->label),
		   gcamn->x, gcamn->y, gcamn->z, dist, target_dist) ;
	  
	  check_gcam(gcam) ;
	  dist -= (target_dist) ;
	  sse += (dist*dist) ;
	}

      }
  return(sse) ;
}

int
GCAMsetTargetDistancesForLabel(GCA_MORPH *gcam, MRI *mri_labels, MRI *mri_dist, int label)
{
  int             x, y, z, l, xv, yv, zv ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcam->gca)
        {
          GCApriorToVoxel(gcam->gca, mri_labels, x, y, z, &xv, &yv, &zv) ;
        }
        else
        {
          xv = x ;
          yv = y;
          zv = z ;
        }
        l = MRIgetVoxVal(mri_labels, xv, yv, zv, 0) ;
        if (l != label)  // only sample for this label
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn->target_dist = MRIgetVoxVal(mri_dist, xv, yv, zv, 0) ;
      }
    }
  }
  return(NO_ERROR) ;
}
int
GCAMsetMeansForLabel(GCA_MORPH *gcam, MRI *mri_labels, MRI *mri_vals,int label)
{

  int             x, y, z, l, r, xv, yv, zv ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcam->gca)
        {
          GCApriorToVoxel(gcam->gca, mri_labels, x, y, z, &xv, &yv, &zv) ;
        }
        else
        {
          xv = x ;
          yv = y ;
          zv = z ;
        }
        l = MRIgetVoxVal(mri_labels, xv, yv, zv, 0) ;
        if (l != label)  // only sample for this label
        {
          continue ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (gcamn->gc == NULL)
        {
          continue ;
        }
        for (r = 0 ; r < gcam->ninputs ; r++)
        {
          gcamn->gc->means[r] = MRIgetVoxVal(mri_vals, xv, yv, zv, r) ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}
MRI *
GCAMbuildLabelVolume(GCA_MORPH *gcam, MRI *mri)
{
  int            x,  y, z, xn, yn, zn, width, depth, height, label ;
  GCA_MORPH_NODE *gcamn ;

  // error check
  if (!mri)
    ErrorExit
    (ERROR_BADPARM, "GCAbuildLabelVolume called with null MRI.\n");
  if (mri->width != gcam->atlas.width
      || mri->height != gcam->atlas.height
      || mri->depth != gcam->atlas.depth)
    ErrorExit(ERROR_BADPARM, \
              "GCAbuildMostLikelyVolume called with mri dimension "
              "being different from M3D.\n");

  // set direction cosines etc.
  useVolGeomToMRI(&gcam->atlas, mri);

  width = mri->width ;
  depth = mri->depth ;
  height = mri->height ;

  if (gcam->gca)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          if (!GCAvoxelToPrior(gcam->gca, mri, x, y, z, &xn, &yn, &zn))
          {
            if (xn == Gx && yn == Gy && zn == Gz)
            {
              DiagBreak() ;
            }
            gcamn = &gcam->nodes[xn][yn][zn] ;

            if (gcamn->invalid == GCAM_POSITION_INVALID)
            {
              continue;
            }

            label = GCAgetMaxPriorLabel(gcam->gca, xn, yn, zn, NULL) ;
            MRIsetVoxVal(mri, x, y, z, 0, label) ;
          }// !GCA
        }
      }
    }
  }
  else   // no gca
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }

          gcamn = &gcam->nodes[x][y][z] ;

          if (gcamn->invalid == GCAM_POSITION_INVALID)
          {
            continue;
          }

          MRIsetVoxVal(mri, x, y, z, 0, gcamn->label) ;
        }
      }
    }
  }


  return(mri) ;
}

double
gcamMultiscaleEnergy(GCA_MORPH *gcam, MRI *mri)
{
  GCAM_MS         *gcam_ms = (GCAM_MS *)gcam->vgcam_ms ;
  double          vox_pval, total_log_pval, pval ;
  int             x, y, z, s, xn, yn, zn, n ;
  GCA             *gca ;
  GCA_NODE        *gcan ;
  GCA_PRIOR       *gcap;
  GCA_MORPH_NODE  *gcamn ;
  float           vals[MAX_GCA_INPUTS] ;

  for (total_log_pval = 0.0, x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->ninputs) ;
        GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn);
        for (vox_pval = 0.0, s = 0 ; s < gcam_ms->nsigmas ; s++)
        {
          gca = gcam_ms->gcas[s] ;
          gcan = &gca->nodes[xn][yn][zn] ;
          gcap = &gca->priors[x][y][z] ;
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            pval = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs, xn,yn,zn,gca) ;
            if (!finite(pval))
            {
              DiagBreak() ;
            }
            vox_pval += pval ;
          }
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("E_ms:   node(%d,%d,%d) vox sse %2.3f\n", x, y, z, vox_pval) ;
        }
        if (DZERO(vox_pval))
        {
          DiagBreak() ;
        }
        else
        {
          total_log_pval += log(vox_pval) ;
        }
      }

  return(-total_log_pval) ;
}


int
gcamMultiscaleTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_multiscale)
{
  GCAM_MS         *gcam_ms = (GCAM_MS *)gcam->vgcam_ms ;
  double          pval, dx, dy, dz, Ix, Iy, Iz, Ierror, norm ;
  int             x, y, z, s, xn, yn, zn, n ;
  GCA             *gca ;
  GCA_NODE        *gcan ;
  GC1D            *gc ;
  GCA_PRIOR       *gcap;
  GCA_MORPH_NODE  *gcamn ;
  float           vals[MAX_GCA_INPUTS] ;


  if (FZERO(l_multiscale))
  {
    return(NO_ERROR) ;
  }
  {
    static int callno = 0 ;
    char  fname[STRLEN] ;
    MRI   *mri_sigma, *mri_aseg, *mri_pvals, *mri_s_index ;
    TRANSFORM transform;

    transform.type = MORPH_3D_TYPE;
    transform.xform = (void *)gcam ;
    if (gcam->mri_xind == NULL)
    {
      TransformInvert(&transform, mri) ;
    }
    mri_aseg = MRIcloneDifferentType(mri, MRI_SHORT) ;
    mri_s_index = MRIcloneDifferentType(mri, MRI_SHORT) ;
    mri_pvals = MRIcloneDifferentType(mri, MRI_FLOAT) ;
    mri_sigma = GCAMMScomputeOptimalScale(gcam_ms, &transform, mri, NULL, mri_aseg, mri_pvals, mri_s_index) ;
    sprintf(fname, "sigma%d.mgz", callno) ;
    printf("writing optimal sigma to %s...\n", fname) ;
    MRIwrite(mri_sigma, fname) ;
    MRIfree(&mri_sigma) ;

    sprintf(fname, "aseg%d.mgz", callno) ;
    printf("writing optimal sigma aseg to %s...\n", fname) ;
    MRIwrite(mri_aseg, fname) ;
    MRIfree(&mri_aseg) ;

    sprintf(fname, "pvals%d.mgz", callno) ;
    printf("writing optimal sigma aseg pvals to %s...\n", fname) ;
    MRIwrite(mri_pvals, fname) ;
    MRIfree(&mri_pvals) ;
    callno++ ;
  }

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;

        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }
        load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->ninputs) ;
        GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn);
        MRIsampleVolumeGradientFrame(mri_smooth,
                                     gcamn->x, gcamn->y, gcamn->z,
                                     &Ix, &Iy, &Iz, 0) ;  // should be vectorized
        norm = sqrt(Ix*Ix+Iy*Iy+Iz*Iz) ;
        if (FZERO(norm))  // no gradient - don't know which way to go
        {
          continue ;
        }
        // don't worry about magnitude of gradient
        Ix /= norm ;
        Iy /= norm ;
        Iz /= norm ;

        for (dx = dy = dz = 0.0, s = 0 ; s < gcam_ms->nsigmas ; s++)
        {
          if (Gx == x && y == Gy && z == Gz)
            printf("s = %2.3f, I=%d, delI = (%2.2f, %2.2f, %2.2f)\n",
                   gcam_ms->sigmas[s], (int)vals[0], Ix, Iy, Iz) ;
          gca = gcam_ms->gcas[s] ;
          gcap = &gca->priors[x][y][z] ;
          gcan = &gca->nodes[xn][yn][zn] ;
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]) ;
            //            gc = &gcan->gcs[n] ;
            if (gc == NULL)
            {
              continue ;
            }
            Ierror = gc->means[0] - vals[0] ;  // divided by covars[0]?
            pval = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs,xn,yn,zn,gca) ;
            if (!finite(pval))
            {
              DiagBreak() ;
            }
            if (Gx == x && y == Gy && z == Gz && !FZERO(pval))
            {
              printf("               MS(%d, %s) = %2.5f * %2.0f\n",
                     n, cma_label_to_name(gcap->labels[n]), pval, Ierror) ;
            }
            dx += pval * Ix * Ierror ;
            dy += pval * Iy * Ierror ;
            dz += pval * Iz * Ierror ;
          }
        }
        if (Gx == x && y == Gy && z == Gz)
        {
          printf("            multiscale grad=(%2.3f, %2.3f, %2.3f)\n", dx, dy, dz) ;
        }
        gcamn->dx += l_multiscale * dx ;
        gcamn->dy += l_multiscale * dy ;
        gcamn->dz += l_multiscale * dz ;
      }

  return(NO_ERROR) ;
}

MRI *
GCAMMScomputeOptimalScale(GCAM_MS *gcam_ms, TRANSFORM *transform, MRI *mri_inputs, MRI *mri_sigma,
                          MRI *mri_aseg, MRI *mri_pvals, MRI *mri_s_index)

{
  int        x, y, z, s, best_s  ;
  GCA_NODE  *gcan ;
  GCA_PRIOR *gcap ;
  int       xp, yp, zp, xn, yn, zn, width, height, depth, n, max_label, best_label ;
  double    p, best_p, total_p, max_p ;
  float     vals[MAX_GCA_INPUTS] ;
  double    sigma = 2, stest, sdist, best_dist ;
  GCA       *gca ;
  MRI       *mri_kernel = MRIgaussian1d(sigma, 100) ;

  TransformInvert(transform, mri_inputs) ;
  width = mri_inputs->width ;
  height = mri_inputs->height;
  depth = mri_inputs->depth ;

  if (mri_sigma == NULL)
  {
    mri_sigma = MRIcloneDifferentType(mri_inputs, MRI_FLOAT) ;
  }

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (x == Gx2 && y == Gy2 && z == Gz2)
        {
          DiagBreak() ;
        }

        if (!GCAsourceVoxelToNode(gcam_ms->gcas[0], mri_inputs,
                                  transform, x, y, z, &xn, &yn, &zn))
        {
          load_vals(mri_inputs, x, y, z, vals, gcam_ms->gcas[0]->ninputs);

          GCAsourceVoxelToPrior(gcam_ms->gcas[0], mri_inputs, transform, x, y, z, &xp, &yp, &zp);
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          best_p = -1 ;
          best_label = best_s = -1 ;
          if (xp == Gx && yp == Gy && zp == Gz)
          {
            DiagBreak() ;
          }
          for (s = 0 ; s < gcam_ms->nsigmas ; s++)
          {
            sigma = gcam_ms->sigmas[s] ;
            gca = gcam_ms->gcas[s] ;
            gcap = &gca->priors[xp][yp][zp] ;
            if (GCAsourceVoxelToNode(gca, mri_inputs,
                                     transform, x, y, z, &xn, &yn, &zn) != NO_ERROR)
            {
              continue ;
            }
            gcan = &gca->nodes[xn][yn][zn] ;

            // marginalize over labels
            total_p = 0 ;
            max_p = -1.0 ;
            max_label = -1 ;
            for (n = 0 ; n < gcap->nlabels ; n++)
            {
              p = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs,
                                             xn, yn, zn, gca) ;
              if (Gdiag_no == gcap->labels[n] && x == Gx2 && y == Gy2 && z == Gz2)
                printf("\ts %2.1f: l %s, p = %2.5f, xp=(%d, %d, %d), prior=%2.5f\n",
                       sigma, cma_label_to_name(Gdiag_no), p, xp, yp, zp, gcap->priors[n]) ;
              if (p > max_p)
              {
                max_p = p ;
                max_label = gcap->labels[n] ;
              }
              total_p += p ;
            }

            if (max_p > best_p)
            {
              best_p = max_p ;
              best_label = max_label ;
              best_s = s ;
            }
          }
          if (mri_aseg)
          {
            MRIsetVoxVal(mri_aseg, x, y, z, 0, best_label) ;
          }
          if (mri_s_index)
          {
            MRIsetVoxVal(mri_s_index, x, y, z, 0, best_s) ;
          }

          sigma = gcam_ms->sigmas[best_s] ;
          if ((x == Ggca_x && y == Ggca_y && z == Ggca_z) || (x == Gx2 && y == Gy2 && z == Gz2))
            printf("setting optimal scale at (%d, %d, %d) to %2.3f, l=%s, p=%2.5f\n",
                   x, y, z, sigma, cma_label_to_name(best_label), best_p) ;

          MRIsetVoxVal(mri_sigma, x, y, z, 0, sigma) ;
        }
      }
    }
  }

  // spatially smooth optimal scales, then use these to recompute best label
  MRIconvolveGaussian(mri_sigma, mri_sigma, mri_kernel) ;

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (x == Gx2 && y == Gy2 && z == Gz2)
        {
          DiagBreak() ;
        }

        if (!GCAsourceVoxelToNode(gcam_ms->gcas[0], mri_inputs,
                                  transform, x, y, z, &xn, &yn, &zn))
        {
          load_vals(mri_inputs, x, y, z, vals, gcam_ms->gcas[0]->ninputs);

          GCAsourceVoxelToPrior(gcam_ms->gcas[0], mri_inputs, transform, x, y, z, &xp, &yp, &zp);
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          if (x == Gx && y == Gy && z == Gz)
          {
            DiagBreak() ;
          }
          if (xp == Gx && yp == Gy && zp == Gz)
          {
            DiagBreak() ;
          }
          best_p = -1 ;
          best_label = -1 ;
          best_s = -1 ;
          best_dist = 1000000 ;
          sigma = MRIgetVoxVal(mri_sigma, x, y, z, 0) ;
          for (s = 0 ; s < gcam_ms->nsigmas ; s++)
          {
            stest = gcam_ms->sigmas[s] ;
            sdist = fabs(stest-sigma) ;
            if (sdist < best_dist)
            {
              best_dist = sdist ;
              best_s = s ;
            }
          }

          gca = gcam_ms->gcas[best_s] ;
          GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn) ;
          gcan = &gca->nodes[xn][yn][zn] ;
          gcap = &gca->priors[xp][yp][zp] ;

          // marginalize over labels
          total_p = 0 ;
          max_p = -1.0 ;
          max_label = -1 ;
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            p = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs,
                                           xn, yn, zn, gca) ;
            if (p > best_p)
            {
              best_p = p ;
              best_label = gcap->labels[n] ;
            }
            total_p += p ;
          }

          //          best_p /= total_p ;
          if (mri_aseg)
          {
            MRIsetVoxVal(mri_aseg, x, y, z, 0, best_label) ;
          }
          if (mri_pvals) // normalize p-value
          {
            int    s, second_label = 0, second_s = 0 ;
            double second_p = 0 ;

            MRIsetVoxVal(mri_pvals, x, y, z, 0, best_p) ;
            if (x == Gx2 && y == Gy2 && z == Gz2)
              printf("setting optimal scale at (%d, %d, %d) to %2.3f, l=%s, p=%2.5f\n",
                     x, y, z, sigma, cma_label_to_name(best_label), best_p) ;
            if (mri_aseg->nframes > 1)
            {
              MRIsetVoxVal(mri_aseg, x, y, z, 1, best_p) ;

              // find 2nd optimal label and scale
              for (s = 0 ; s < gcam_ms->nsigmas ; s++)
              {
                gca = gcam_ms->gcas[s] ;
                if (GCAsourceVoxelToNode(gca, mri_inputs,
                                         transform, x, y, z, &xn, &yn, &zn) != NO_ERROR)
                {
                  continue ;
                }
                gcan = &gca->nodes[xn][yn][zn] ;
                gcap = &gca->priors[xp][yp][zp] ;

                // marginalize over labels
                for (n = 0 ; n < gcap->nlabels ; n++)
                {
                  p = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs,
                                                 xn, yn, zn, gca) ;
                  if (( gcap->labels[n] != best_label) && (p > second_p))
                  {
                    second_p = p ;
                    second_label = gcap->labels[n] ;
                    second_s = s ;
                  }
                }
              }

              gca = gcam_ms->gcas[second_s] ;
              GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn) ;
              gcan = &gca->nodes[xn][yn][zn] ;
              gcap = &gca->priors[xp][yp][zp] ;
              total_p = 0 ;
              for (n = 0 ; n < gcap->nlabels ; n++)
              {
                p = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs,
                                               xn, yn, zn, gca) ;
                total_p += p ;
              }
              //              second_p /= total_p ;
              MRIsetVoxVal(mri_aseg, x, y, z, 2, second_label) ;
              MRIsetVoxVal(mri_aseg, x, y, z, 3, second_p) ;
              if (x == Gx2 && y == Gy2 && z == Gz2)
                printf("setting 2nd optimal scale at (%d, %d, %d) to %2.3f, l=%s, p=%2.5f\n",
                       x, y, z, gcam_ms->sigmas[second_s], cma_label_to_name(second_label), second_p) ;
            }
          }
          sigma = gcam_ms->sigmas[best_s] ;
          if ((x == Ggca_x && y == Ggca_y && z == Ggca_z) ||
              (x == Gx2 && y == Gy2 && z == Gz2))
          {
            printf("setting optimal scale at (%d, %d, %d) to %2.3f\n", x, y, z, sigma) ;
          }

          MRIsetVoxVal(mri_sigma, x, y, z, 0, sigma) ;
        }
      }
    }
  }

  MRIfree(&mri_kernel) ;
  return(mri_sigma) ;
}
#if 1
double
GCAMlogPosterior(GCA_MORPH *gcam, MRI *mri_inputs)
{
  double     total_log_posterior, logp ;
  int        x, y, z, real_label  ;
  GCA_NODE  *gcan ;
  GCA_PRIOR *gcap ;
  int       xn, yn, zn, n, best_label ;
  double    p, best_p, total_p, max_p ;
  float     vals[MAX_GCA_INPUTS] ;
  GCA       *gca ;
  GCA_MORPH_NODE *gcamn ;

  gca = gcam->gca ;

  GCAMinvert(gcam, mri_inputs) ;
  for (total_log_posterior = 0.0, x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        gcamn =  &gcam->nodes[x][y][z] ;
        if (gcamn->invalid == GCAM_POSITION_INVALID)
        {
          continue;
        }
        load_vals(mri_inputs, gcamn->x, gcamn->y, gcamn->z, vals, gca->ninputs);

        max_p = best_p = -1 ;
        best_label = -1 ;
        GCApriorToNode(gca, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        gcap = &gca->priors[x][y][z] ;

        // marginalize over labels
        total_p = 0 ;
        for (real_label = n = 0 ; n < gcap->nlabels ; n++)
        {
          if (IS_UNKNOWN(gcap->labels[n]))
          {
            continue ;
          }
          real_label = 1 ;
          p = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs,xn,yn,zn,gca) ;
          if (p > max_p)
          {
            max_p = p ;
            best_label = gcap->labels[n] ;
          }
          if (!finite(p))
          {
            p = 0 ;
          }
          if (Gx == x && Gy == y && Gz == z)
            printf("node (%d, %d, %d) --> MRI(%2.0f, %2.0f, %2.0f) = %d, label = %s, pval = %2.5f\n",
                   x, y, z, gcamn->x, gcamn->y, gcamn->z, (int)(vals[0]),
                   cma_label_to_name(gcap->labels[n]),p) ;
          total_p += p ;
        }
        if (real_label == 0)
        {
          continue ;
        }
        logp = log(total_p) ;
        total_log_posterior += logp ;
      }
    }
  }

  return(total_log_posterior) ;
}
#else
double
GCAMlogPosteriorDensity(GCA_MORPH *gcam, MRI *mri_inputs)
{
  double     total_log_posterior ;
  int        x, y, z  ;
  GCA_NODE  *gcan ;
  GCA_PRIOR *gcap ;
  int       xp, yp, zp, xn, yn, zn, width, height, depth, n, best_label ;
  double    p, best_p, total_p, max_p ;
  float     vals[MAX_GCA_INPUTS] ;
  GCA       *gca ;
  TRANSFORM _transform, *transform = &_transform ;

  gca = gcam->gca ;
  width = mri_inputs->width ;
  height = mri_inputs->height;
  depth = mri_inputs->depth ;
  transform->type = MORPH_3D_TYPE;
  transform->xform = (void *)gcam ;
  if (gcam->mri_xind == NULL)
  {
    TransformInvert(transform, mri_inputs) ;
  }

  for (total_log_posterior = 0.0, x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }

        if (!GCAsourceVoxelToNode(gca, mri_inputs,
                                  transform, x, y, z, &xn, &yn, &zn))
        {
          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

          GCAsourceVoxelToPrior(gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp);
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          max_p = best_p = -1 ;
          best_label = -1 ;
          gcan = &gca->nodes[xn][yn][zn] ;
          gcap = &gca->priors[xp][yp][zp] ;

          // marginalize over labels
          total_p = 0 ;
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            p = GCAcomputePosteriorDensity(gcap, gcan, -1, n, vals, gca->ninputs,gca,xn,yn,zn) ;
            if (p > max_p)
            {
              max_p = p ;
              best_label = gcap->labels[n] ;
            }
            total_p += p ;
          }
        }
        total_log_posterior += log(total_p) ;
      }
    }
  }

  return(total_log_posterior) ;
}
#endif
MRI *
GCAMcreateDistanceTransforms(MRI *mri_source,
			     MRI *mri_target,
			     MRI *mri_all_dtrans,
			     float max_dist,
			     GCA_MORPH *gcam,
			     MRI **pmri_atlas_dist_map)
{
  MRI   *mri_dtrans, *mri_atlas_dtrans, *mri_atlas_dist_map ;
  int   frame, label ;
  char  fname[STRLEN] ;

  mri_all_dtrans = MRIallocSequence(mri_source->width,
                                    mri_source->height,
                                    mri_source->depth,
                                    MRI_FLOAT,
                                    NDTRANS_LABELS) ;
  MRIcopyHeader(mri_target, mri_all_dtrans) ; // should this be mri_source????

  mri_atlas_dist_map = MRIallocSequence(mri_target->width, mri_target->height, mri_target->depth, MRI_FLOAT, NDTRANS_LABELS) ;

  MRIcopyHeader(mri_target, mri_all_dtrans) ;

#ifdef HAVE_OPENMP
  label = 0;
  mri_atlas_dtrans = mri_dtrans = NULL ;
#pragma omp parallel for firstprivate(label, fname, mri_atlas_dtrans, mri_source, mri_target, mri_all_dtrans, mri_dtrans, max_dist, Gdiag_no, gcam)    schedule(static,1)
#endif
  for (frame = 0 ; frame < NDTRANS_LABELS ; frame++)
  {
    label = dtrans_labels[frame] ;
    printf("creating distance transform for %s, frame %d...\n", cma_label_to_name(label), frame) ;

    mri_dtrans = MRIdistanceTransform(mri_source, NULL, label, max_dist, DTRANS_MODE_SIGNED, NULL) ;
    MRIcopyFrame(mri_dtrans, mri_all_dtrans, 0, frame) ;
    mri_atlas_dtrans = MRIdistanceTransform(mri_target, NULL, label, max_dist, DTRANS_MODE_SIGNED, NULL) ;
    MRInormalizeInteriorDistanceTransform(mri_atlas_dtrans,
                                          mri_dtrans,
                                          mri_atlas_dtrans) ;
    MRIcopyFrame(mri_atlas_dtrans, mri_atlas_dist_map, 0, frame) ;
    GCAMsetTargetDistancesForLabel(gcam,
                                   mri_target,
                                   mri_atlas_dtrans,
                                   dtrans_labels[frame]);

    if (Gdiag & DIAG_WRITE && (DIAG_VERBOSE_ON || Gdiag_no == label))
    {
      sprintf(fname, "source.%s.mgz", cma_label_to_name(label)) ;
      printf("writing distance transform to %s\n", fname) ;
      MRIwrite(mri_dtrans, fname) ;
      sprintf(fname, "atlas.%s.mgz", cma_label_to_name(label)) ;
      printf("writing distance transform to %s\n", fname) ;
      MRIwrite(mri_atlas_dtrans, fname) ;
    }
    MRIfree(&mri_dtrans) ;
    MRIfree(&mri_atlas_dtrans) ;
  }

  mri_all_dtrans->outside_val = max_dist ;
  mri_atlas_dist_map->outside_val = max_dist ;
  if (pmri_atlas_dist_map)
    *pmri_atlas_dist_map = mri_atlas_dist_map ;
  else
    MRIfree(&mri_atlas_dist_map) ;
  return(mri_all_dtrans) ;
}

#if 0
MRI *
GCAMcreateDistanceTransforms(GCA_MORPH *gcam, MRI *mri, MRI *mri_all_dtrans,
                             MRI **pmri_atlas_dtrans, float max_dist)
{
  MRI   *mri_dtrans, *mri_labels ;
  int   frame ;
  char  fname[STRLEN] ;

  if (mri_all_dtrans == NULL)
    mri_all_dtrans = MRIallocSequence(mri->width, mri->height, mri->depth,
                                      MRI_FLOAT, NDTRANS_LABELS) ;
  if (*pmri_atlas_dtrans == NULL)
    *pmri_atlas_dtrans = MRIallocSequence(mri->width, mri->height, mri->depth,
                                          MRI_FLOAT, NDTRANS_LABELS) ;

  MRIcopyHeader(mri, mri_all_dtrans) ;

  mri_labels = MRIcloneDifferentType(mri, MRI_FLOAT) ;
  GCAMbuildLabelVolume(gcam, mri_labels) ;
  MRIwrite(mri_labels, "labels.mgz") ;
#ifdef HAVE_OPENMP
#pragma omp parallel for firstprivate(fname, mri_labels, mri_all_dtrans, dtrans_labels,mri_dtrans, mri, max_dist, Gdiag_no, gcam)  schedule(static,1)
#endif
  for (frame = 0 ; frame < NDTRANS_LABELS ; frame++)
  {
    printf("creating distance transform for %s, frame %d...\n", cma_label_to_name(dtrans_labels[frame]), frame) ;
    if (dtrans_labels[frame] == Gdiag_no)
    {
      DiagBreak() ;
    }
    mri_dtrans = MRIdistanceTransform(mri, NULL, dtrans_labels[frame], max_dist, DTRANS_MODE_SIGNED, NULL) ;
    if (Gdiag & DIAG_WRITE && (DIAG_VERBOSE_ON || dtrans_labels[frame] == Gdiag_no))
    {
      sprintf(fname, "%s.mgz", cma_label_to_name(dtrans_labels[frame])) ;
      MRIwrite(mri_dtrans, fname) ;
    }
    MRIcopyFrame(mri_dtrans, mri_all_dtrans, 0, frame) ;
    MRIfree(&mri_dtrans) ;

    /*
      now set distance of each gcam node to the boundary of the label (will be negative
      for interior points)
    */
    mri_dtrans =
      MRIdistanceTransform(mri_labels, NULL, dtrans_labels[frame], max_dist, DTRANS_MODE_SIGNED, NULL) ;
    GCAMsetTargetDistancesForLabel(gcam, mri_labels, mri_dtrans, dtrans_labels[frame]);
    MRIcopyFrame(mri_dtrans, *pmri_atlas_dtrans, 0, frame) ;
    MRIfree(&mri_dtrans) ;
  }

  mri_all_dtrans->outside_val = max_dist ;
  (*pmri_atlas_dtrans)->outside_val = max_dist ;
  MRIfree(&mri_labels) ;
  return(mri_all_dtrans) ;
}
#endif

int
GCAMwriteDistanceTransforms(GCA_MORPH *gcam,  MRI *mri_source_dist_map, MRI *mri_atlas_dist_map,
                            const char *dist_name)
{
  char   fname[STRLEN] ;
  int    i, label ;

  for (i = 0 ; i < NDTRANS_LABELS ; i++)
  {
    label = dtrans_labels[i] ;
    sprintf(fname, "%s_source_frame%d_%s.mgz", dist_name, i, cma_label_to_name(label)) ;
    MRIwriteFrame(mri_source_dist_map, fname, i) ;
    sprintf(fname, "%s_atlas_frame%d_%s.mgz", dist_name, i, cma_label_to_name(label)) ;
    MRIwriteFrame(mri_atlas_dist_map, fname, i) ;
  }
  return(NO_ERROR) ;
}
MRI *
GCAMreadDistanceTransforms(GCA_MORPH *gcam, MRI *mri, MRI *mri_all_dtrans,
                           MRI **pmri_atlas_dtrans, const char *dist_name, double max_dist)
{
  MRI   *mri_tmp ;
  int    i, label ;
  char   fname[STRLEN] ;

  if (mri_all_dtrans == NULL)
    mri_all_dtrans = MRIallocSequence(mri->width, mri->height, mri->depth,
                                      MRI_FLOAT, NDTRANS_LABELS) ;
  if (*pmri_atlas_dtrans == NULL)
    *pmri_atlas_dtrans = MRIallocSequence(mri->width, mri->height, mri->depth,
                                          MRI_FLOAT, NDTRANS_LABELS) ;

  MRIcopyHeader(mri, mri_all_dtrans) ;

  for (i = 0 ; i < NDTRANS_LABELS ; i++)
  {
    label = dtrans_labels[i] ;
    sprintf(fname, "%s_source_frame%d_%s.mgz", dist_name, i, cma_label_to_name(label)) ;
    mri_tmp = MRIread(fname) ;
    if (mri_tmp == NULL)
    {
      ErrorReturn(NULL, (ERROR_NOFILE, "GCAMreadDistanceTransforms: could not read volume %s", fname)) ;
    }
    MRIcopyFrame(mri_tmp, mri_all_dtrans, 0, i) ;

    sprintf(fname, "%s_atlas_frame%d_%s.mgz", dist_name, i, cma_label_to_name(label)) ;
    mri_tmp = MRIread(fname) ;
    if (mri_tmp == NULL)
    {
      ErrorReturn(NULL, (ERROR_NOFILE, "GCAMreadDistanceTransforms: could not read volume %s", fname)) ;
    }
    MRIcopyFrame(mri_tmp, *pmri_atlas_dtrans, 0, i) ;
  }

  mri_all_dtrans->outside_val = max_dist ;
  (*pmri_atlas_dtrans)->outside_val = max_dist ;
  return(mri_all_dtrans) ;
}


/*
  compute a gauss-newton step matching the source and atlas images
  and put the resulting vector field into mri_v
*/
int
MRIcomputeGaussNewtonStep(MRI *mri_source, MRI *mri_atlas, MRI *mri_v, MRI *mri_labels, GCA_MORPH_PARMS *parms)
{
  MATRIX *m_J, *m_J_t, *m_inv, *m_I, *m_tmp = NULL, *v_r, *v_v, *m_tmp2 = NULL, *m_D ;
  int    x, y, z, frame, label ;
  double Ix, Iy, Iz, dx, dy, dz, tval, sval, error, lambda = .1 ;

  m_D = MatrixIdentity(3, NULL) ;
  MatrixScalarMul(m_D, lambda, m_D) ;
  m_J = MatrixAlloc(mri_atlas->nframes, 3, MATRIX_REAL) ;
  m_J_t = MatrixTranspose(m_J, NULL) ;
  m_inv = MatrixAlloc(3, 3, MATRIX_REAL) ;
  m_I = MatrixIdentity(3, NULL) ;
  v_r = VectorAlloc(mri_atlas->nframes, MATRIX_REAL) ;
  v_v = VectorAlloc(3, MATRIX_REAL) ;

  for (x = 0 ; x < mri_atlas->width ; x++)
    for (y = 0 ; y < mri_atlas->height ; y++)
      for (z = 0 ; z < mri_atlas->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_labels, x, y, z, 0) ;
        frame = dtrans_label_to_frame(parms, label) ;
#if 0
        if (frame < 0)
        {
          continue ;
        }
#endif
        dx = MRIgetVoxVal(mri_v, x, y, z, 0) ;
        dy = MRIgetVoxVal(mri_v, x, y, z, 1) ;
        dz = MRIgetVoxVal(mri_v, x, y, z, 2) ;

#if 0
        tval = MRIgetVoxVal(mri_atlas, x, y, z, frame) ;
        MRIsampleVolumeFrameType(mri_source, x+dx, y+dy, z+dz, frame, SAMPLE_TRILINEAR,
                                 &sval) ;
        MRIsampleVolumeGradientFrame(mri_source,
                                     x+dx, y+dy, z+dz,
                                     &Ix, &Iy, &Iz, frame) ;
        error = (tval-sval) ;
        *MATRIX_RELT(m_J, 1, 1) = error*Ix ;
        *MATRIX_RELT(m_J, 1, 2) = error*Iy ;
        *MATRIX_RELT(m_J, 1, 3) = error*Iz ;

        VECTOR_ELT(v_r, 1) = error ;

#else
        for (frame = 0 ; frame < mri_atlas->nframes ; frame++)
        {
          if (dtrans_labels[frame] == label)
          {
            DiagBreak() ;
          }
          tval = MRIgetVoxVal(mri_atlas, x, y, z, frame) ;
          MRIsampleVolumeFrameType(mri_source, x+dx, y+dy, z+dz, frame, SAMPLE_TRILINEAR,
                                   &sval) ;
          MRIsampleVolumeGradientFrame(mri_source,
                                       x+dx, y+dy, z+dz,
                                       &Ix, &Iy, &Iz, frame) ;
#if 0
          if (sval > 2.0)
          {
            error = 0 ;
          }
          else
#endif
            error = (tval-sval) ;
          *MATRIX_RELT(m_J, frame+1, 1) = error*Ix ;
          *MATRIX_RELT(m_J, frame+1, 2) = error*Iy ;
          *MATRIX_RELT(m_J, frame+1, 3) = error*Iz ;

          VECTOR_ELT(v_r, frame+1) = error ;
        }
#endif
        MatrixTranspose(m_J, m_J_t) ;
        m_tmp = MatrixMultiply(m_J_t, m_J, m_tmp) ;
        MatrixAdd(m_tmp, m_D, m_tmp) ;  // Marquardt-Levemberg
        if (MatrixInverse(m_tmp, m_inv) == NULL)
        {
          continue ;
        }
        m_tmp2 = MatrixMultiply(m_inv, m_J_t, m_tmp2) ;
        MatrixMultiply(m_tmp2, v_r, v_v) ;
        if (VectorLen(v_v) > 1e5)
        {
          DiagBreak() ;
        }
        if (VectorLen(v_v) > 1)
        {
          DiagBreak() ;
        }
        MRIsetVoxVal(mri_v, x, y, z, 0, V3_X(v_v)) ;
        MRIsetVoxVal(mri_v, x, y, z, 1, V3_Y(v_v)) ;
        MRIsetVoxVal(mri_v, x, y, z, 2, V3_Z(v_v)) ;
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("Gauss-Newton(%d, %d, %d): v = (%2.2f, %2.2f, %2.2f)\n",
                 Gx, Gy, Gz, V3_X(v_v), V3_Y(v_v), V3_Z(v_v)) ;
        }
      }

  MatrixFree(&m_J) ;
  MatrixFree(&m_J_t) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_inv);
  MatrixFree(&v_r);
  MatrixFree(&m_D) ;
  MatrixFree(&v_v) ;
  return(NO_ERROR) ;
}

#if 1
int
gcamScaleAndSquare(MRI *mri_s_old, MRI *mri_s_new)
{
  int   x, y, z ;
  double dx, dy, dz, dx_old, dy_old, dz_old ;

  for (x = 0 ; x < mri_s_old->width ; x++)
    for (y = 0 ; y < mri_s_old->height ; y++)
      for (z = 0 ; z < mri_s_old->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
#if 0
        dx_old = MRIgetVoxVal(mri_s_old, x, y, z, 0) ;
        dy_old = MRIgetVoxVal(mri_s_old, x, y, z, 1) ;
        dz_old = MRIgetVoxVal(mri_s_old, x, y, z, 2) ;
#else
        dx_old = MRIFseq_vox(mri_s_old, x, y, z, 0) ;
        dy_old = MRIFseq_vox(mri_s_old, x, y, z, 1) ;
        dz_old = MRIFseq_vox(mri_s_old, x, y, z, 2) ;
#endif
        MRIsampleVolumeFrameType(mri_s_old, x+dx_old, y+dy_old, z+dz_old,
                                 0, SAMPLE_TRILINEAR, &dx) ;
        MRIsampleVolumeFrameType(mri_s_old, x+dx_old, y+dy_old, z+dz_old,
                                 1, SAMPLE_TRILINEAR, &dy) ;
        MRIsampleVolumeFrameType(mri_s_old, x+dx_old, y+dy_old, z+dz_old,
                                 2, SAMPLE_TRILINEAR, &dz) ;
#if 0
        MRIsetVoxVal(mri_s_new, x, y, z, 0, dx_old+dx) ;
        MRIsetVoxVal(mri_s_new, x, y, z, 1, dy_old+dy) ;
        MRIsetVoxVal(mri_s_new, x, y, z, 2, dz_old+dz) ;
#else
        MRIFseq_vox(mri_s_new, x, y, z, 0)  = dx_old+dx ;
        MRIFseq_vox(mri_s_new, x, y, z, 1)  = dy_old+dy ;
        MRIFseq_vox(mri_s_new, x, y, z, 2)  = dz_old+dz ;
#endif
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("gcamScaleAndSquare(%d, %d, %d): Dx = (%2.2f, %2.2f, %2.2f)\n",
                 Gx, Gy, Gz,dx_old+dx, dy_old+dy, dz_old+dz) ;
        }
      }

  return(NO_ERROR) ;
}
#endif
MRI *
MRIcomposeWarps(MRI *mri_warp1, MRI *mri_warp2, MRI *mri_dst)
{
  int   x, y, z ;
  double dx2, dy2, dz2, dx1, dy1, dz1 ;

  if (mri_dst == NULL)
  {
    mri_dst = MRIclone(mri_warp1, NULL) ;
  }

  for (x = 0 ; x < mri_warp1->width ; x++)
    for (y = 0 ; y < mri_warp1->height ; y++)
      for (z = 0 ; z < mri_warp1->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        // should use new update to sample old warp, then do vector addition
#if 0  // used to be 0 before Thomas broke everything
        dx2 = MRIgetVoxVal(mri_warp2, x, y, z, 0) ;
        dy2 = MRIgetVoxVal(mri_warp2, x, y, z, 1) ;
        dz2 = MRIgetVoxVal(mri_warp2, x, y, z, 2) ;
        MRIsampleVolumeFrameType(mri_warp1, x+dx2, y+dy2, z+dz2,
                                 0, SAMPLE_TRILINEAR, &dx1) ;
        MRIsampleVolumeFrameType(mri_warp1, x+dx2, y+dy2, z+dz2,
                                 1, SAMPLE_TRILINEAR, &dy1) ;
        MRIsampleVolumeFrameType(mri_warp1, x+dx2, y+dy2, z+dz2,
                                 2, SAMPLE_TRILINEAR, &dz1) ;


#else
        dx1 = MRIgetVoxVal(mri_warp1, x, y, z, 0) ;
        dy1 = MRIgetVoxVal(mri_warp1, x, y, z, 1) ;
        dz1 = MRIgetVoxVal(mri_warp1, x, y, z, 2) ;
#if 0
        MRIsampleVolumeFrameType(mri_warp2, x+dx1, y+dy1, z+dz1,
                                 0, SAMPLE_TRILINEAR, &dx2) ;
        MRIsampleVolumeFrameType(mri_warp2, x+dx1, y+dy1, z+dz1,
                                 1, SAMPLE_TRILINEAR, &dy2) ;
        MRIsampleVolumeFrameType(mri_warp2, x+dx1, y+dy1, z+dz1,
                                 2, SAMPLE_TRILINEAR, &dz2) ;
#else
        dx2 = MRIgetVoxVal(mri_warp2, x, y, z, 0) ;
        dy2 = MRIgetVoxVal(mri_warp2, x, y, z, 1) ;
        dz2 = MRIgetVoxVal(mri_warp2, x, y, z, 2) ;
#endif
#endif

        MRIsetVoxVal(mri_dst, x, y, z, 0, dx1+dx2) ;
        MRIsetVoxVal(mri_dst, x, y, z, 1, dy1+dy2) ;
        MRIsetVoxVal(mri_dst, x, y, z, 2, dz1+dz2) ;
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("MRIcomposeWarps(%d, %d, %d): dx1 = (%2.2f, %2.2f, %2.2f), "
                 "dx2 = (%2.2f, %2.2f, %2.2f), dx_total = (%2.2f, %2.2f, %2.2f)\n",
                 Gx, Gy, Gz, dx1, dy1, dz1, dx2, dy2, dz2,
                 dx1+dx2, dy1+dy2, dz1+dz2) ;
        }
      }

  if (Gx2 >= 0)  // find the atlas node that most closes maps to this voxel
  {
    float  min_dist, dist, x2, y2, z2, xf, yf, zf ;

    min_dist = 20*(mri_warp1->width+mri_warp1->height+mri_warp1->depth) ;

    xf = yf = zf = -1 ;
    for (x = 0 ; x < mri_warp1->width ; x++)
      for (y = 0 ; y < mri_warp1->height ; y++)
        for (z = 0 ; z < mri_warp1->depth ; z++)
        {
          dx1 = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
          dy1 = MRIgetVoxVal(mri_dst, x, y, z, 1) ;
          dz1 = MRIgetVoxVal(mri_dst, x, y, z, 2) ;
          x2 = x+dx1 ;
          y2 = y+dy1 ;
          z2 = z+dz1 ;
          dist = sqrt(SQR(x2-Gx2) + SQR(y2-Gy2) + SQR(z2-Gz2)) ;
          if (dist < min_dist)
          {
            min_dist = dist ;
            Gx = x ;
            Gy = y ;
            Gz = z ;
            xf = x2 ;
            yf = y2;
            zf = z2 ;
            DiagBreak() ;
          }
        }
    DiagBreak() ;
    printf("I (%d, %d, %d) <=- A (%2.1f, %2.1f, %2.1f), dist=%2.1f\n", Gx2, Gy2, Gz2, xf, yf, zf, min_dist) ;
  }

  return(mri_dst) ;
}
double MRImorphSSE(MRI *mri_source, MRI *mri_atlas, MRI *mri_morph) ;

MRI *MRIcomputeDistanceTransformStep(MRI *mri_source, MRI *mri_atlas, MRI *mri_delta,
                                     MRI *mri_labels, GCA_MORPH_PARMS *parms) ;


/*
  mri_source and mri_atlas are distance transforms of labeled volume.
*/

#if 0
int
check_inverse(MRI *mri_fwd, MRI *mri_inv)
{
  int    x, y, z ;
  float  dx, dy, dz ;
  double invx, invy, invz ;

  for (x = 0 ; x < mri_fwd->width ; x++)
    for (y = 0 ; y < mri_fwd->height ; y++)
      for (z = 0 ; z < mri_fwd->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        dx = MRIgetVoxVal(mri_fwd, x, y,z, 0) ;
        dy = MRIgetVoxVal(mri_fwd, x, y,z, 1) ;
        dz = MRIgetVoxVal(mri_fwd, x, y,z, 2) ;
        MRIsampleVolumeFrame(mri_inv, x+dx, y+dy, z+dz, 0, &invx) ;
        MRIsampleVolumeFrame(mri_inv, x+dx, y+dy, z+dz, 1, &invy) ;
        MRIsampleVolumeFrame(mri_inv, x+dx, y+dy, z+dz, 2, &invz) ;
        invx += (x+dx) ;
        invy += (x+dy) ;
        invz += (x+dz) ;
      }
  return(NO_ERROR) ;
}
#endif

int
GCAMdemonsRegister(GCA_MORPH *gcam, MRI *mri_source_labels,
                   MRI *mri_atlas_labels, GCA_MORPH_PARMS *parms,
                   double max_dist, MRI *mri_pvals, MRI *mri_atlas_dtrans_orig)
{
  double vmax_allowed, dt, N, old_sse, new_sse, vmax, pct_change ;
  MRI    *mri_s_new, *mri_s_old, *mri_kernel, *mri_morphed,
         *mri_current_dtrans = NULL,
          *mri_warp, *mri_current_labels, *mri_dtrans,
          *mri_atlas_dtrans ;
  int    done, step ;
  char   fname[STRLEN] ;

  mri_s_old = MRIallocSequence(mri_atlas_labels->width,
                               mri_atlas_labels->height, mri_atlas_labels->depth,
                               MRI_FLOAT, 3) ;
  mri_s_new = MRIclone(mri_s_old, NULL) ;
  mri_warp = MRIclone(mri_s_old, NULL) ;  // start with identity warp
  GCAMwriteWarpToMRI(gcam, mri_warp) ;

  parms->mri = mri_source_labels ;

  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s.log", parms->base_name) ;
    if (parms->log_fp == NULL)
    {
      if (parms->start_t == 0)
      {
        parms->log_fp = fopen(fname, "w") ;
      }
      else
      {
        parms->log_fp = fopen(fname, "a") ;
      }
    }
  }
  else
  {
    parms->log_fp = NULL ;
  }

  if (Gx > 0)
  {
    int x, y, z ;
    x = mri_atlas_labels->xi[nint(Gx + MRIgetVoxVal(mri_warp, Gx, Gy, Gz, 0))] ;
    y = mri_atlas_labels->yi[nint(Gy + MRIgetVoxVal(mri_warp, Gx, Gy, Gz, 1))] ;
    z = mri_atlas_labels->zi[nint(Gz + MRIgetVoxVal(mri_warp, Gx, Gy, Gz, 2))] ;
    printf("GCAMdemonsRegister: (%d, %d, %d) = %s in source, and (%d, %d, %d) = %s in atlas\n",
           Gx, Gy, Gz,
           cma_label_to_name(MRIgetVoxVal(mri_source_labels,Gx,Gy,Gz,0)),
           x, y, z,
           cma_label_to_name(MRIgetVoxVal(mri_atlas_labels,x,y,z,0))) ;
  }

  if (Gx2 >= 0)
  {
    float xf, yf, zf ;

    xf = yf = zf = 0 ;
    GCAMsampleInverseMorph(gcam, Gx2, Gy2, Gz2, &xf, &yf, &zf) ;
    Gx = nint(xf) ;
    Gy = nint(yf) ;
    Gz = nint(zf) ;
    printf("GCAMdemonsRegister: (%d, %d, %d) = %s in source, and (%d, %d, %d) = %s in atlas\n",
           Gx2, Gy2, Gz2,
           cma_label_to_name(MRIgetVoxVal(mri_source_labels,Gx2,Gy2,Gz2,0)),
           Gx, Gy, Gz,
           cma_label_to_name(MRIgetVoxVal(mri_atlas_labels,Gx,Gy,Gz,0))) ;
  }

  mri_kernel = MRIgaussian1d(parms->sigma, 100) ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s_target", parms->base_name) ;
    MRIwriteImageViews(mri_atlas_labels, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_target.mgh", parms->base_name) ;
    MRIwrite(mri_atlas_labels, fname) ;
  }

  /* user specifies atlas distance transform so it
     doesn't need to be recomputed every time */
#if 0
  if (mri_atlas_dtrans_orig)
  {
    mri_atlas_dtrans = MRIcopy(mri_atlas_dtrans_orig, NULL) ;
  }
  else
    mri_atlas_dtrans = MRIcreateDistanceTransforms(mri_atlas_labels, NULL,
                       max_dist,
                       parms->dtrans_labels,
                       parms->ndtrans) ;
#else
  mri_atlas_dtrans = mri_atlas_dtrans_orig ;
#endif
  step = parms->start_t ;
  mri_current_labels = MRIapplyMorph(mri_source_labels, mri_warp, NULL, SAMPLE_NEAREST) ;
  /*
    mri_current_labels is the label map using the current warp.
    mri_source_labels are the original (unwarped labels)
    mri_dtrans is the distance transform of the set of labels
    that were current when this iteration of demons started
    mri_current_dtrans is the distance transform after the current warp.
  */

  if (Gdiag & DIAG_WRITE && step == parms->start_t)
  {
    char fname[STRLEN] ;
    MRI  *mri_morphed ;

    sprintf(fname, "%s_source%3.3d.mgh", parms->base_name, parms->start_t) ;
    MRIwrite(mri_current_labels, fname) ;
    sprintf(fname, "%s_source%3.3d", parms->base_name, parms->start_t) ;
    MRIwriteImageViews(mri_current_labels, fname, IMAGE_SIZE) ;

    if (parms->mri_diag)
    {
      MRI *mri ;
      mri = MRIapplyMorph(parms->mri_diag, mri_warp, NULL, SAMPLE_TRILINEAR) ;
      sprintf(fname, "%s_intensity%3.3d.mgh", parms->base_name, step) ;
      MRIwrite(mri, fname) ;
      sprintf(fname, "%s_intensity%3.3d", parms->base_name, step) ;
      MRIwriteImageViews(mri, fname, IMAGE_SIZE) ;
      MRIfree(&mri) ;
    }
    mri_morphed = GCAMmorphFromAtlas(mri_atlas_labels, gcam, NULL, SAMPLE_NEAREST);
    sprintf(fname, "%s_atlas%3.3d.mgh", parms->base_name, parms->start_t) ;
    MRIwrite(mri_morphed, fname) ;
    MRIfree(&mri_morphed) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRI *mri_morphed, *mri_warp2 ;

      GCAMreadWarpFromMRI(gcam, mri_warp) ;
      mri_warp2 = MRIclone(mri_warp, NULL) ;
      GCAMwriteWarpToMRI(gcam, mri_warp2) ;

      mri_morphed = GCAMmorphToAtlas(mri_source_labels, gcam, NULL, 0, SAMPLE_NEAREST) ;
      sprintf(fname, "%s_source_gcam.mgh", parms->base_name) ;
      MRIwrite(mri_morphed, fname) ;
      MRIfree(&mri_morphed) ;
      MRIfree(&mri_warp2) ;
    }
  }
  /*
    normalize the interior (<0) distances in the atlas to have the same min val (max neg)
    as the src images so that we can compute true voxel distances.
  */
  mri_dtrans = MRIcreateDistanceTransforms(mri_current_labels, NULL,
               max_dist,
               parms->dtrans_labels,
               parms->ndtrans) ;
  MRInormalizeInteriorDistanceTransform(mri_atlas_dtrans, mri_dtrans, mri_atlas_dtrans) ;
  old_sse = MRIlabelMorphSSE(mri_source_labels, mri_atlas_labels, mri_warp) ;
  mri_current_dtrans = MRIcopy(mri_dtrans, mri_current_dtrans) ;
  do  // run demons on the current and atlas distance transforms
  {
    MRIcomputeDistanceTransformStep(mri_current_dtrans, mri_atlas_dtrans,
                                    mri_s_old, mri_atlas_labels,parms) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "%s_dtstep%3.3d.mgz", parms->base_name, step) ;
      printf("writing step direction to %s\n", fname) ;
      MRIwrite(mri_s_old, fname) ;
    }

    vmax = MRImaxNorm(mri_s_old) ;
    vmax_allowed = 0.125 ;  // actually 0.5 of pixel spacing would ok
    vmax_allowed = 0.25 ;  // actually 0.5 of pixel spacing would ok
    N = ceil(log2(vmax/vmax_allowed)) ;
    N = MAX(6.0, N) ;
    dt = 1.0/pow(2,N);
    MRIscalarMul(mri_s_old, mri_s_old, dt) ;
    if (mri_pvals)
    {
      MRImultiply(mri_s_old, mri_pvals, mri_s_old) ;
    }

    for ( ; dt <= 1.0 ; dt *= 2.0)
    {
#if 1
      gcamScaleAndSquare(mri_s_old, mri_s_new) ;
#else
      MRIcomposeWarps(mri_s_old, mri_s_old, mri_s_new) ;
#endif
      MRIcopy(mri_s_new, mri_s_old) ;
    }

    step++ ;
    MRIconvolveGaussian(mri_s_new, mri_s_new, mri_kernel) ;
    if (Gx > 0)
      printf("after smoothing, step (%2.2f, %2.2f, %2.2f)\n",
             MRIgetVoxVal(mri_s_new, Gx, Gy, Gz, 0),
             MRIgetVoxVal(mri_s_new, Gx, Gy, Gz, 1),
             MRIgetVoxVal(mri_s_new, Gx, Gy, Gz, 2)) ;
    MRIcomposeWarps(mri_warp, mri_s_new, mri_warp) ;
    GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri_warp) ;

    MRIapplyMorph(mri_current_dtrans, mri_s_new, mri_dtrans, SAMPLE_TRILINEAR) ;
    MRIcopy(mri_dtrans, mri_current_dtrans) ;
    MRIapplyMorph(mri_source_labels, mri_warp, mri_current_labels, SAMPLE_NEAREST) ;
    new_sse = MRIlabelMorphSSE(mri_source_labels, mri_atlas_labels, mri_warp) ;
    pct_change = 100*(old_sse-new_sse) / (old_sse) ;
    printf("step %d: old sse %2.3f, new sse %2.3f (%2.3f%%)\n",
           step, old_sse, new_sse, pct_change) ;
    if (parms->log_fp)
    {
      fprintf(parms->log_fp, "step %d: old sse %2.3f, new sse %2.3f (%2.3f%%)\n",
              step, old_sse, new_sse, pct_change) ;
      fflush(parms->log_fp) ;
    }

    done = pct_change < parms->tol ; // (new_sse > old_sse) ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN];

      sprintf(fname, "%s_source%3.3d.mgh", parms->base_name, step) ;
      MRIwrite(mri_current_labels, fname) ;
      sprintf(fname, "%s_source%3.3d", parms->base_name, step) ;
      MRIwriteImageViews(mri_current_labels, fname, IMAGE_SIZE) ;
      if (parms->mri_diag)
      {
        MRI *mri ;
        mri = MRIapplyMorph(parms->mri_diag, mri_warp, NULL, SAMPLE_TRILINEAR) ;
        sprintf(fname, "%s_intensity%3.3d.mgh", parms->base_name, step) ;
        MRIwrite(mri, fname) ;
        sprintf(fname, "%s_intensity%3.3d", parms->base_name, step) ;
        MRIwriteImageViews(mri, fname, IMAGE_SIZE) ;
        MRIfree(&mri) ;
      }
      if (parms->write_fname)
      {
        MRI *mri ;
        GCAMreadWarpFromMRI(gcam, mri_warp) ;
        printf("writing results of level to %s...\n", parms->write_fname) ;
        //            GCAMvoxToRas(gcam) ;
        GCAMwrite(gcam, parms->write_fname) ;
        mri =
          GCAMmorphToAtlas(mri_source_labels,gcam, NULL, -1, SAMPLE_NEAREST) ;

        sprintf(fname, "%s_gcam%3.3d.mgh", parms->base_name, step) ;
        MRIwrite(mri, fname) ;
        MRIfree(&mri) ;
      }
      if (DIAG_VERBOSE_ON)
      {
        if (gcam->mri_xind)
        {
          MRIfree(&gcam->mri_xind) ;
          MRIfree(&gcam->mri_yind) ;
          MRIfree(&gcam->mri_zind) ;
        }
        GCAMreadWarpFromMRI(gcam, mri_warp) ;
        GCAMinvert(gcam, mri_source_labels);
        mri_morphed = GCAMmorphFromAtlas(mri_atlas_labels, gcam, NULL, SAMPLE_NEAREST);
        sprintf(fname, "%s_atlas%3.3d.mgh", parms->base_name, step) ;
        MRIwrite(mri_morphed, fname) ;
        MRIfree(&mri_morphed) ;
      }
    }
    old_sse = new_sse ;
  }
  while (!done);
  MRIfree(&mri_current_labels) ;
  MRIfree(&mri_current_dtrans) ;
  MRIfree(&mri_dtrans) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    char fname[STRLEN] ;
    step++ ;
    mri_morphed = MRIapplyMorph(mri_source_labels, mri_warp, NULL, SAMPLE_NEAREST) ;
    sprintf(fname, "%s_source%3.3d.mgz", parms->base_name, step) ;
    MRIwrite(mri_morphed, fname) ;
    sprintf(fname, "%s_source%3.3d", parms->base_name, step) ;
    MRIwriteImageViews(mri_morphed, fname, IMAGE_SIZE) ;
    MRIfree(&mri_morphed) ;
  }
  parms->start_t = (step);
  parms->last_sse = new_sse ;
  //  GCAMreadWarpFromMRI(gcam, mri_warp) ;
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri_warp) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRI *mri_m1, *mri_m2, *mri_warp2 ;

    mri_warp2 = MRIclone(mri_warp, NULL) ;
    GCAMwriteWarpToMRI(gcam, mri_warp2) ;

    mri_morphed = GCAMmorphToAtlas(mri_source_labels, gcam, NULL, 0, SAMPLE_NEAREST) ;
    sprintf(fname, "%s_source_gcam.mgz", parms->base_name) ;

    mri_m1 = MRIapplyMorph(mri_source_labels, mri_warp, NULL, SAMPLE_NEAREST) ;
    mri_m2 = MRIapplyMorph(mri_source_labels, mri_warp2, NULL, SAMPLE_NEAREST) ;
    MRIwrite(mri_m1, "m1.mgz") ;
    MRIwrite(mri_m2, "m2.mgz") ;

    MRIwrite(mri_morphed, fname) ;
    MRIfree(&mri_morphed) ;
    MRIfree(&mri_warp2) ;
    MRIfree(&mri_m1) ;
    MRIfree(&mri_m2) ;
  }
  MRIfree(&mri_s_new) ;
  MRIfree(&mri_s_old) ;
  MRIfree(&mri_kernel) ;
  MRIfree(&mri_warp) ;
#if 0
  MRIfree(&mri_atlas_dtrans) ;
#endif

  if (gcam->mri_xind)
  {
    MRIfree(&gcam->mri_xind) ;
    MRIfree(&gcam->mri_yind) ;
    MRIfree(&gcam->mri_zind) ;
  }
  GCAMinvert(gcam, mri_source_labels);
  gcamComputeMetricProperties(gcam) ;

  return(NO_ERROR) ;
}
int
GCAMcompose(GCA_MORPH *gcam, MRI *mri_warp)
{
  int             xp, yp, zp, xv, yv, zv ;
  double          dx, dy, dz ;
  GCA_MORPH_NODE  *gcamn ;

  for (xp = 0 ; xp < gcam->width; xp++)
  {
    for (yp = 0 ; yp < gcam->height ; yp++)
    {
      for (zp = 0 ; zp < gcam->depth ; zp++)
      {
        GCApriorToVoxel(gcam->gca, mri_warp, xp, yp, zp, &xv, &yv, &zv) ;
        gcamn = &gcam->nodes[xp][yp][zp] ;
        MRIsampleVolumeFrameType(mri_warp, gcamn->x, gcamn->y, gcamn->z, 0, SAMPLE_TRILINEAR, &dx);
        MRIsampleVolumeFrameType(mri_warp, gcamn->x, gcamn->y, gcamn->z, 1, SAMPLE_TRILINEAR, &dy);
        MRIsampleVolumeFrameType(mri_warp, gcamn->x, gcamn->y, gcamn->z, 2, SAMPLE_TRILINEAR, &dz);
        gcamn->x += dx ;
        gcamn->y += dy ;
        gcamn->z += dz ;
      }
    }
  }
  return(NO_ERROR) ;
}
int
GCAMreadWarpFromMRI( GCA_MORPH *gcam, const MRI *mri_warp )
{
  int             xp, yp, zp, xv, yv, zv ;
  double          dx, dy, dz ;
  GCA_MORPH_NODE  *gcamn ;

  if (gcam->width  == mri_warp->width &&
      gcam->height == mri_warp->height &&
      gcam->depth  == mri_warp->depth)
  {
    for (xp = 0 ; xp < gcam->width ; xp++)
      for (yp = 0 ; yp < gcam->height ; yp++)
        for (zp = 0 ; zp < gcam->depth ; zp++)
        {
          if (xp == Gx && yp == Gy && zp == Gz)
          {
            DiagBreak() ;
          }
          gcamn = &gcam->nodes[xp][yp][zp] ;
          dx = MRIgetVoxVal(mri_warp, xp, yp, zp, 0) ;
          dy = MRIgetVoxVal(mri_warp, xp, yp, zp, 1) ;
          dz = MRIgetVoxVal(mri_warp, xp, yp, zp, 2) ;
          gcamn->x = gcamn->origx+dx ;
          gcamn->y = gcamn->origy+dy ;
          gcamn->z = gcamn->origz+dz ;
        }
    return(NO_ERROR) ;
  }

  for (xp = 0 ; xp < gcam->width; xp++)
  {
    xv = (xp * gcam->spacing) ;
    for (yp = 0 ; yp < gcam->height ; yp++)
    {
      yv = (yp * gcam->spacing) ;
      for (zp = 0 ; zp < gcam->depth ; zp++)
      {
        if (xp == Gx && yp == Gy && zp == Gz)
        {
          DiagBreak() ;
        }
        zv = (zp * gcam->spacing) ;
        //        GCApriorToVoxel(gcam->gca, mri_warp, xp, yp, zp, &xv, &yv, &zv) ;
        gcamn = &gcam->nodes[xp][yp][zp] ;
        MRIsampleVolumeFrameType(mri_warp, xv, yv, zv, 0, SAMPLE_TRILINEAR, &dx);
        MRIsampleVolumeFrameType(mri_warp, xv, yv, zv, 1, SAMPLE_TRILINEAR, &dy);
        MRIsampleVolumeFrameType(mri_warp, xv, yv, zv, 2, SAMPLE_TRILINEAR, &dz);
        gcamn->x = xv+dx ;
        gcamn->y = yv+dy ;
        gcamn->z = zv+dz ;
      }
    }
  }
  return(NO_ERROR) ;
}

int
GCAMreadInverseWarpFromMRI(GCA_MORPH *gcam, MRI *mri_warp)
{
  int             xv, yv, zv ;
  double          dx, dy, dz ;

  if (gcam->mri_xind)
  {
    MRIfree(&gcam->mri_xind) ;
    MRIfree(&gcam->mri_yind) ;
    MRIfree(&gcam->mri_zind) ;
  }

  gcam->mri_xind = MRIalloc(mri_warp->width, mri_warp->height, mri_warp->depth, MRI_FLOAT) ;
  gcam->mri_yind = MRIalloc(mri_warp->width, mri_warp->height, mri_warp->depth, MRI_FLOAT) ;
  gcam->mri_zind = MRIalloc(mri_warp->width, mri_warp->height, mri_warp->depth, MRI_FLOAT) ;
  for (xv = 0 ; xv < mri_warp->width; xv++)
  {
    for (yv = 0 ; yv < mri_warp->height ; yv++)
    {
      for (zv = 0 ; zv < mri_warp->depth ; zv++)
      {
        dx = MRIgetVoxVal(mri_warp, xv, yv, zv, 0) ;
        dy = MRIgetVoxVal(mri_warp, xv, yv, zv, 1) ;
        dz = MRIgetVoxVal(mri_warp, xv, yv, zv, 2) ;
        MRIsetVoxVal(gcam->mri_xind, xv, yv, zv, 0, (xv+dx)/gcam->spacing) ;
        MRIsetVoxVal(gcam->mri_yind, xv, yv, zv, 0, (yv+dy)/gcam->spacing) ;
        MRIsetVoxVal(gcam->mri_zind, xv, yv, zv, 0, (zv+dz)/gcam->spacing) ;
      }
    }
  }
  return(NO_ERROR) ;
}

MRI *
MRIapplyMorph(MRI *mri_source, MRI *mri_warp, MRI *mri_dst, int sample_type)
{
  int    x, y, z, f ;
  double dx, dy, dz, xd, yd, zd, val ;

  if (mri_dst == NULL)
  {
    mri_dst = MRIclone(mri_source, NULL) ;
  }

  MRI_BSPLINE * bspline = NULL;
  if (sample_type == SAMPLE_CUBIC_BSPLINE)
  {
    bspline = MRItoBSpline(mri_source,NULL,3);
  }

  for (x = 0 ; x < mri_dst->width ; x++)
    for (y = 0 ; y < mri_dst->height ; y++)
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        dx = MRIgetVoxVal(mri_warp, x, y, z, 0) ;
        dy = MRIgetVoxVal(mri_warp, x, y, z, 1) ;
        dz = MRIgetVoxVal(mri_warp, x, y, z, 2) ;
        xd = x+dx ;
        yd = y+dy ;
        zd = z+dz ;
        for (f = 0 ; f < mri_source->nframes ; f++)
        {
          if (sample_type == SAMPLE_CUBIC_BSPLINE)
            // recommended to externally call this and keep mri_coeff
            // if image is resampled often (e.g. in registration algo)
          {
            MRIsampleBSpline(bspline, xd, yd, zd, f, &val);
          }
          else
          {
            MRIsampleVolumeFrameType(mri_source, xd, yd, zd, f, sample_type, &val);
          }
          MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
        }
      }

  mri_dst->c_r = mri_warp->c_r ;
  mri_dst->c_a = mri_warp->c_a ;
  mri_dst->c_s = mri_warp->c_s ;
  if (bspline)
  {
    MRIfreeBSpline(&bspline);
  }
  return(mri_dst) ;
}

double
MRImorphSSE(MRI *mri_source, MRI *mri_atlas, MRI *mri_warp)
{
  int    x, y, z, f ;
  double dx, dy, dz, xd, yd, zd, sval, tval, error, sse ;

  for (sse = 0.0, x = 0 ; x < mri_atlas->width ; x++)
    for (y = 0 ; y < mri_atlas->height ; y++)
      for (z = 0 ; z < mri_atlas->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        dx = MRIgetVoxVal(mri_warp, x, y, z, 0) ;
        dy = MRIgetVoxVal(mri_warp, x, y, z, 1) ;
        dz = MRIgetVoxVal(mri_warp, x, y, z, 2) ;
        xd = x+dx ;
        yd = y+dy ;
        zd = z+dz ;
        for (f = 0 ; f < mri_source->nframes ; f++)
        {
          tval = MRIgetVoxVal(mri_atlas, x, y, z, f) ;
          MRIsampleVolumeFrameType(mri_source, xd, yd, zd, f, SAMPLE_TRILINEAR, &sval);
          error = tval-sval ;
          sse += (error*error) ;
        }
      }


  sse = sse / (mri_atlas->width*mri_atlas->height*mri_atlas->depth*mri_atlas->nframes);
  return(sse) ;
}
MRI *
MRIcreateDistanceTransforms(MRI *mri, MRI *mri_all_dtrans, float max_dist,
                            int *labels, int nlabels)
{
  MRI   *mri_dtrans ;
  int   frame ;
  char  fname[STRLEN] ;

  if (mri_all_dtrans == NULL)
    mri_all_dtrans = MRIallocSequence(mri->width, mri->height, mri->depth,
                                      MRI_FLOAT, nlabels) ;

  MRIcopyHeader(mri, mri_all_dtrans) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri, "labels.mgz") ;
  }
  for (frame = 0 ; frame < nlabels ; frame++)
  {
    printf("creating distance transform for %s, frame %d...\n",
           cma_label_to_name(labels[frame]), frame) ;
    if (dtrans_labels[frame] == Gdiag_no)
    {
      DiagBreak() ;
    }
    mri_dtrans = MRIdistanceTransform(mri, NULL, labels[frame], max_dist,
                                      DTRANS_MODE_SIGNED, NULL) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "%s.mgz", cma_label_to_name(dtrans_labels[frame])) ;
      MRIwrite(mri_dtrans, fname) ;
    }
    MRIcopyFrame(mri_dtrans, mri_all_dtrans, 0, frame) ;
    MRIfree(&mri_dtrans) ;
  }

  mri_all_dtrans->outside_val = max_dist ;
  return(mri_all_dtrans) ;
}
MRI *
MRInormalizeInteriorDistanceTransform(MRI *mri_src_dist, MRI *mri_ref_dist, MRI *mri_dst_dist)
{
  int    frame ;
  float  min_ref_dist, max_val, min_src_dist ;

  if (mri_dst_dist == NULL)
  {
    mri_dst_dist = MRIcopy(mri_src_dist, NULL) ;
  }
  for (frame = 0 ; frame < mri_src_dist->nframes ; frame++)
  {
    MRIvalRangeFrame(mri_src_dist, &min_src_dist, &max_val, frame) ;
    MRIvalRangeFrame(mri_ref_dist, &min_ref_dist, &max_val, frame) ;
    if (FZERO(min_src_dist))
    {
      continue ;
    }
    MRIscalarMulFrame(mri_src_dist, mri_dst_dist, min_ref_dist / min_src_dist, frame) ;
  }
  return(mri_dst_dist) ;
}
MRI *
MRIcomputeDistanceTransformStep(MRI *mri_source, MRI *mri_atlas, MRI *mri_delta,
                                MRI *mri_atlas_labels, GCA_MORPH_PARMS *parms)
{
  int    x, y, z, frame, label ;
  double Ix, Iy, Iz, dx, dy, dz, tval, sval, error, norm ;

  if (mri_delta == NULL)
  {
    mri_delta = MRIallocSequence(mri_atlas->width, mri_atlas->height, mri_atlas->depth,MRI_FLOAT,3);
  }
  else
  {
    MRIclear(mri_delta) ;
  }

  for (x = 0 ; x < mri_atlas->width ; x++)
    for (y = 0 ; y < mri_atlas->height ; y++)
      for (z = 0 ; z < mri_atlas->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_atlas_labels, x, y, z, 0) ;
        frame = dtrans_label_to_frame(parms, label) ;
        if (frame < 0)
        {
          continue ;
        }
        dx = MRIgetVoxVal(mri_delta, x, y, z, 0) ;
        dy = MRIgetVoxVal(mri_delta, x, y, z, 1) ;
        dz = MRIgetVoxVal(mri_delta, x, y, z, 2) ;

        tval = MRIgetVoxVal(mri_atlas, x, y, z, frame) ;
#if 0
        // used to be if 1 before Thomas broke everything
        MRIsampleVolumeFrameType(mri_source, x+dx, y+dy, z+dz, frame, SAMPLE_TRILINEAR,
                                 &sval) ;
        MRIsampleVolumeGradientFrame(mri_source,
                                     x+dx, y+dy, z+dz,
                                     &Ix, &Iy, &Iz, frame) ;
#else
        sval = MRIgetVoxVal(mri_source, x, y, z, 0) ;
        MRIsampleVolumeGradientFrame(mri_source, x, y, z,
                                     &Ix, &Iy, &Iz, frame) ;
#endif
        norm = sqrt(Ix*Ix + Iy*Iy + Iz*Iz) ;
        if (!FZERO(norm))
        {
          Ix /= norm ;
          Iy /= norm ;
          Iz /= norm ;
        }
        error = (tval-sval) ;
        MRIsetVoxVal(mri_delta, x, y, z, 0, error*Ix) ;
        MRIsetVoxVal(mri_delta, x, y, z, 1, error*Iy) ;
        MRIsetVoxVal(mri_delta, x, y, z, 2, error*Iz) ;

        if (x == Gx && y == Gy && z == Gz)
        {
          printf("MRIcomputeDistanceTransformStep: v(%d, %d, %d) %s = (%2.2f, %2.2f, %2.2f), A=%2.1f, I=%2.1f\n",
                 Gx, Gy, Gz, cma_label_to_name(label),error*Ix, error*Iy, error*Iz, tval, sval);
        }
      }

  return(mri_delta) ;
}
double
MRIlabelMorphSSE(MRI *mri_source, MRI *mri_atlas, MRI *mri_warp)
{
  int    x, y, z, src_label, atlas_label, nvox ;
  double dx, dy, dz, xd, yd, zd, sval, sse ;

  for (nvox = 0, sse = 0.0, x = 0 ; x < mri_atlas->width ; x++)
    for (y = 0 ; y < mri_atlas->height ; y++)
      for (z = 0 ; z < mri_atlas->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        atlas_label = (int)MRIgetVoxVal(mri_atlas, x, y, z, 0) ;
        if (!IS_UNKNOWN(atlas_label))
        {
          nvox++ ;
        }
        dx = MRIgetVoxVal(mri_warp, x, y, z, 0) ;
        dy = MRIgetVoxVal(mri_warp, x, y, z, 1) ;
        dz = MRIgetVoxVal(mri_warp, x, y, z, 2) ;
        xd = x+dx ;
        yd = y+dy ;
        zd = z+dz ;
        src_label = (int)MRIgetVoxVal(mri_source, x, y, z, 0) ;
        MRIsampleVolumeFrameType(mri_source, xd, yd, zd, 0, SAMPLE_NEAREST, &sval);
        src_label = nint(sval) ;
        if (src_label != atlas_label)
        {
          sse += 1 ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("atlas: %s (%d, %d, %d) + (%2.1f, %2.1f, %2.1f) = %s (%2.1f, %2.1f, %2.1f)\n",
                 cma_label_to_name(atlas_label), x, y, z, dx, dy, dz,
                 cma_label_to_name(src_label), xd, yd, zd) ;
      }


  //  sse = sse / (mri_atlas->width*mri_atlas->height*mri_atlas->depth*mri_atlas->nframes);
  sse = sse / nvox ;
  return(sse) ;
}
MRI *
GCAMwriteWarpToMRI( const GCA_MORPH *gcam, MRI *mri_warp)
{
  int    x, y, z;
  double xp, yp, zp ;
  float  xw, yw, zw ;

  if (mri_warp == NULL)  // warp is same dimensions as gcam by construction
  {
    GCA_MORPH_NODE *gcamn ;

    mri_warp = MRIallocSequence(gcam->width, gcam->height, gcam->depth, MRI_FLOAT, 3) ;
    MRIcopyVolGeomToMRI(mri_warp, &gcam->atlas) ;
    MRIsetResolution(mri_warp, gcam->atlas.xsize*gcam->spacing, 
		     gcam->atlas.ysize*gcam->spacing, 
		     gcam->atlas.zsize*gcam->spacing) ;
    for (x = 0 ; x < mri_warp->width ; x++)
      for (y = 0 ; y < mri_warp->height ; y++)
        for (z = 0 ; z < mri_warp->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          MRIsetVoxVal(mri_warp, x, y, z, 0, gcamn->x - gcamn->origx) ;
          MRIsetVoxVal(mri_warp, x, y, z, 1, gcamn->y - gcamn->origy) ;
          MRIsetVoxVal(mri_warp, x, y, z, 2, gcamn->z - gcamn->origz) ;
          if (x == Gx && y == Gy && z == Gz)
            printf("voxel (%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
                   gcamn->origx, gcamn->origy, gcamn->origz, gcamn->x, gcamn->y, gcamn->z) ;
        }
    return(mri_warp) ;
  }

  MRIcopyVolGeomToMRI(mri_warp, &gcam->atlas) ;
  for (x = 0 ; x < mri_warp->width ; x++)
    for (y = 0 ; y < mri_warp->height ; y++)
      for (z = 0 ; z < mri_warp->depth ; z++)
      {
        if (gcam->gca)
        {
          GCAvoxelToPriorReal(gcam->gca, mri_warp, x, y, z, &xp, &yp, &zp) ;
        }
        else
        {
          xp = x ;
          yp = y ;
          zp = z; // same resolution
        }
        GCAMsampleMorph(gcam, x, y, z, &xw, &yw, &zw) ;
        MRIsetVoxVal(mri_warp, x, y, z, 0, xw-x) ;
        MRIsetVoxVal(mri_warp, x, y, z, 1, yw-y) ;
        MRIsetVoxVal(mri_warp, x, y, z, 2, zw-z) ;
        if (x == Gx && y == Gy && z == Gz)
          printf("voxel (%d, %d, %d) --> (%2.1f, %2.1f, %2.1f)\n",
                 x, y, z, xw, yw, zw) ;
      }
  return(mri_warp) ;
}

MRI *
GCAMwriteInverseWarpToMRI(GCA_MORPH *gcam, MRI *mri_warp)
{
  int     x, y, z;
  double  xp, yp, zp ;
  double  xw, yw, zw ;

  if (mri_warp == NULL)
  {
    mri_warp = MRIallocSequence(gcam->width, gcam->height, gcam->depth, MRI_FLOAT, 3) ;
  }

  MRIcopyVolGeomToMRI(mri_warp, &gcam->image) ;
  for (x = 0 ; x < mri_warp->width ; x++)
    for (y = 0 ; y < mri_warp->height ; y++)
      for (z = 0 ; z < mri_warp->depth ; z++)
      {
        GCAvoxelToPriorReal(gcam->gca, mri_warp, x, y, z, &xp, &yp, &zp) ;
        MRIsampleVolume(gcam->mri_xind, x, y, z, &xw) ;
        MRIsampleVolume(gcam->mri_yind, x, y, z, &yw) ;
        MRIsampleVolume(gcam->mri_zind, x, y, z, &zw) ;
        xw *= gcam->spacing ;
        yw *= gcam->spacing ;
        zw *= gcam->spacing ;
        MRIsetVoxVal(mri_warp, x, y, z, 0, xw-x) ;
        MRIsetVoxVal(mri_warp, x, y, z, 1, yw-y) ;
        MRIsetVoxVal(mri_warp, x, y, z, 2, zw-z) ;
        if (x == Gx && y == Gy && z == Gz)
          printf("voxel (%d, %d, %d) --> (%2.1f, %2.1f, %2.1f)\n",
                 x, y, z, xw, yw, zw) ;
      }
  return(mri_warp) ;
}

extern int gca_write_iterations ;
extern char *gca_write_fname ;

#if 1
MRI  *
GCAMSreclassifyUsingGibbsPriors(MRI *mri_inputs,
                                GCAM_MS *gcam_ms,
                                MRI *mri_dst,
                                TRANSFORM *transform,
                                int max_iter, int restart,
                                void (*update_func)(MRI *), MRI *mri_s_index)
{
  int      x, y, z, width, height, depth, label, val, iter,
           n, nchanged, min_changed, index, nindices, old_label, s, max_s ;
  short    *x_indices, *y_indices, *z_indices ;
  GCA_PRIOR *gcap ;
  double   ll, lcma = 0.0, old_posterior, new_posterior, max_posterior ;
  MRI      *mri_changed, *mri_probs /*, *mri_zero */ ;
  GCA      *gca ;


  nindices = mri_dst->width * mri_dst->height * mri_dst->depth ;
  x_indices = (short *)calloc(nindices, sizeof(short)) ;
  y_indices = (short *)calloc(nindices, sizeof(short)) ;
  z_indices = (short *)calloc(nindices, sizeof(short)) ;
  if (!x_indices || !y_indices || !z_indices)
    ErrorExit(ERROR_NOMEMORY, "GCAMSreclassifyUsingGibbsPriors: "
              "could not allocate index set") ;


  width = mri_inputs->width ;
  height = mri_inputs->height;
  depth = mri_inputs->depth ;
  iter = 0 ;
  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
    if (!mri_dst)
    {
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    }
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  mri_changed = MRIclone(mri_dst, NULL) ;

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  // mark changed location
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
      {
        MRIvox(mri_changed,x,y,z) = 1 ;
      }

  if (!restart)
  {
    // calculate the statistics
    old_posterior = GCAMMSgibbsImageLogPosterior(gcam_ms, mri_dst,
                    mri_inputs, transform, mri_s_index) ;
    // get the per voxel value
    old_posterior /= (double)(width*depth*height) ;
  }
  else
  {
    old_posterior = 0 ;
  }


  do
  {
    if (restart)
    {
      for (index = x = 0 ; x < width ; x++)
        for (y = 0 ; y < height ; y++)
          for (z = 0 ; z < depth ; z++)
          {
            // not fixed voxel, but changed
            if (MRIvox(mri_changed,x,y,z) > 0)
            {
              x_indices[index] = x ;
              y_indices[index] = y ;
              z_indices[index] = z ;
              index++ ;
            }
          }
      nindices = index ;
    }
    else if (iter == 0)
    {
      /*      char  fname[STRLEN], *cp ;*/

      if (gca_write_iterations)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s%03d.mgz", gca_write_fname, iter) ;
        printf("writing snapshot to %s...\n", fname) ;
        MRIwrite(mri_dst, fname) ;
      }
      // probs has 0 to 255 values
      mri_probs = GCAlabelProbabilities(mri_inputs,
                                        gcam_ms->gcas[0], NULL, transform) ;
      // sorted according to ascending order of probs
      MRIorderIndices(mri_probs, x_indices, y_indices, z_indices) ;
      MRIfree(&mri_probs) ;
    }
    else
      // randomize the indices value ((0 -> width*height*depth)
      MRIcomputeVoxelPermutation(mri_inputs, x_indices, y_indices,
                                 z_indices) ;

    ///////////////  relabel with neighborhood likelihood  //////////////
    nchanged = 0 ;
    if (G_write_probs && !mri_probs)
    {
      mri_probs = MRIalloc(mri_inputs->width, mri_inputs->height,
                           mri_inputs->depth, MRI_FLOAT) ;
      MRIcopyHeader(mri_inputs, mri_probs) ;
    }
    for (index = 0 ; index < nindices ; index++)
    {
      x = x_indices[index] ;
      y = y_indices[index] ;
      z = z_indices[index] ;
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
      {
        DiagBreak() ;
      }

      // if not marked, don't do anything
      if (MRIvox(mri_changed, x, y, z) == 0)
      {
        continue ;
      }

      // get the grey value
      val = MRIgetVoxVal(mri_inputs, x, y, z, 0) ;

      // save the current label
      // calculate neighborhood likelihood
      max_s = MRIgetVoxVal(mri_s_index, x, y, z, 0) ;
      gca = gcam_ms->gcas[max_s] ;
      max_posterior = GCAnbhdGibbsLogPosterior(gca, mri_dst,
                      mri_inputs, x, y,z,transform, 1);
      label = old_label = nint(MRIgetVoxVal(mri_dst, x, y, z,0)) ;
      /* find the node associated with this coordinate and classify */
      //      for (s = 0 ; s < gcam_ms->nsigmas ; s++)
      for (s = max_s ; s < max_s+1 ; s++)
      {
        gca = gcam_ms->gcas[s] ;
        gcap = getGCAP(gca, mri_inputs, transform, x, y, z) ;
        // it is not in the right place
        if (gcap==NULL)
        {
          continue;
        }
#if 1
        // only one label associated, don't do anything
        if (gcap->nlabels == 1)
        {
          continue ;
        }
#endif

        // go through all labels at this point
        for (n = 0 ; n < gcap->nlabels ; n++)
        {
#if 1
          // skip the current label
          if (gcap->labels[n] == old_label)
          {
            continue ;
          }
#endif
          // assign the new label
          MRIsetVoxVal(mri_dst, x, y, z, 0,gcap->labels[n]) ;
          // calculate neighborhood likelihood
          new_posterior =
            GCAnbhdGibbsLogPosterior(gca, mri_dst,
                                     mri_inputs, x, y,z,transform, 1);
          // if it is bigger than the old one, then replace the label
          // and change max_posterior
          if (new_posterior > max_posterior)
          {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z &&
                (label == Ggca_label || old_label ==
                 Ggca_label || Ggca_label < 0))
              fprintf(stdout,
                      "NbhdGibbsLogLikelihood at (%d, %d, %d):"
                      " old = %d (ll=%.2f) new = %d (ll=%.2f)\n",
                      x, y, z, old_label, max_posterior,
                      gcap->labels[n], new_posterior);

            max_posterior = new_posterior ;
            label = gcap->labels[n] ;
            max_s = s ;
          }
        }

        /*#ifndef __OPTIMIZE__*/
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z &&
            (label == Ggca_label || old_label ==
             Ggca_label || Ggca_label < 0))
        {
          int       xn, yn, zn ;
          GCA_NODE *gcan ;

          if (!GCAsourceVoxelToNode(gca, mri_inputs, transform,
                                    x, y, z, &xn, &yn, &zn))
          {
            gcan = &gca->nodes[xn][yn][zn] ;
            printf("sigma=%2.0f, (%d, %d, %d): old label %s (%d), "
                   "new label %s (%d) (log(p)=%2.3f)\n",
                   gcam_ms->sigmas[max_s], x, y, z, cma_label_to_name(old_label), old_label,
                   cma_label_to_name(label), label, max_posterior) ;
          }
        }
      }

      MRIsetVoxVal(mri_s_index, x, y, z, 0, max_s) ;

      // if label changed
      if (label != old_label)
      {
        nchanged++ ;
        // mark it as changed
        MRIvox(mri_changed, x, y, z) = 1 ;
      }
      else
      {
        MRIvox(mri_changed, x, y, z) = 0 ;
      }
      // assign new label
      MRIsetVoxVal(mri_dst, x, y, z, 0, label) ;
      if (mri_probs)
      {
        MRIsetVoxVal(mri_probs, x, y, z, 0, -max_posterior) ;
      }
    }
    if (mri_probs)
    {
      char fname[STRLEN] ;

      sprintf(fname, "%s%03d.mgz", G_write_probs, iter) ;
      printf("writing probabilities to %s...\n", fname) ;
      MRIwrite(mri_probs, fname) ;
      MRIfree(&mri_probs) ;
    }

    /////////////////////////////////////////////////////////////////////////
    // print info
    if (nchanged > 10000 && iter < 2 && !restart)
    {
      ll = GCAMMSgibbsImageLogPosterior(gcam_ms, mri_dst,
                                        mri_inputs, transform, mri_s_index) ;
      // get the average value
      ll /= (double)(width*depth*height) ;
      if (!FZERO(lcma))
        printf("pass %d: %d changed. image ll: %2.3f "
               "(CMA=%2.3f)\n",
               iter+1, nchanged, ll, lcma) ;
      else // currently this is executed
        printf("pass %d: %d changed. image ll: %2.3f\n",
               iter+1, nchanged, ll) ;
    }
    else
    {
      printf("pass %d: %d changed.\n", iter+1, nchanged) ;
    }

    // get the largest 6 neighbor values,
    // that is, originally 0 could become 1
    MRIdilate(mri_changed, mri_changed) ;
    // if unpdate_func is present, use it
    if (update_func)
    {
      (*update_func)(mri_dst) ;
    }

#define MIN_CHANGED 5000
    min_changed = restart ? 0 : MIN_CHANGED ;
    if (nchanged <= min_changed ||
        (restart && iter >= max_iter))
    {
      if (restart)
      {
        iter = 0 ;
      }

      if (!restart)  /* examine whole volume next time */
      {
        for (x = 0 ; x < width ; x++)
          for (y = 0 ; y < height ; y++)
            for (z = 0 ; z < depth ; z++)
            {
              MRIvox(mri_changed,x,y,z) = 1 ;
            }
      }
    }
    if (gca_write_iterations > 0)
    {
      char fname[STRLEN] ;

      sprintf(fname, "%s_iter%d.mgz", gca_write_fname, iter+1) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRIwrite(mri_dst, fname) ;
    }
    if ((gca_write_iterations > 0) && !(iter % gca_write_iterations))
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s%03d.mgz", gca_write_fname, iter+1) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRIwrite(mri_dst, fname) ;
    }
  }
  while ((nchanged > MIN_CHANGED) && (iter++ < max_iter)) ;



  if (mri_probs)
  {
    MRIfree(&mri_probs) ;
  }

  free(x_indices) ;
  free(y_indices) ;
  free(z_indices) ;
  MRIfree(&mri_changed) ;

  return(mri_dst) ;
}
#endif
double
GCAMMSgibbsImageLogPosterior(GCAM_MS *gcam_ms, MRI *mri_labels, MRI *mri_inputs,
                             TRANSFORM *transform, MRI *mri_s_index)
{
  int    x, y, z, width, depth, height, s ;
  double total_log_posterior, log_posterior ;
  GCA    *gca ;

  width = mri_labels->width ;
  height = mri_labels->height ;
  depth = mri_labels->depth ;

  for (total_log_posterior = 0.0, x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        s = MRIgetVoxVal(mri_s_index, x, y, z, 0) ;
        gca = gcam_ms->gcas[s] ;
        log_posterior =
          GCAvoxelGibbsLogPosterior(gca, mri_labels,
                                    mri_inputs, x, y, z,transform,1.0);
        if (check_finite("gcaGibbsImageLogposterior",
                         log_posterior) == 0)
        {
          DiagBreak() ;
        }
        total_log_posterior += log_posterior ;
        if (total_log_posterior > 1e10)
        {
          DiagBreak() ;
        }
      }
    }
  }
  return(total_log_posterior) ;
}
int
GCAMremoveSingularitiesAndReadWarpFromMRI(GCA_MORPH *gcam, MRI *mri_warp)
{
  int             xp, yp, zp, xv, yv, zv, last_neg, wsize, iter, max_iter,
                  xk, yk, zk, nbhd = 1, max_noprogress, noprogress, min_neg ;
  double          dx, dy, dz, max_nbhd ;
  GCA_MORPH_NODE  *gcamn ;
  MRI             *mri_warp_tmp = NULL, *mri_neg, *mri_neg_orig, *mri_neg_atlas ;

  wsize = 2*gcam->spacing+1 ;
  wsize = 3 ;
  GCAMreadWarpFromMRI(gcam, mri_warp) ;
  gcamComputeMetricProperties(gcam) ;
  if (gcam->neg == 0)
  {
    return(NO_ERROR) ;
  }

  mri_neg = MRIalloc(gcam->image.width, gcam->image.height, gcam->image.depth, MRI_UCHAR) ;
  MRIcopyVolGeomToMRI(mri_warp, &gcam->image);

  mri_neg_orig = MRIalloc(gcam->image.width, gcam->image.height, gcam->image.depth, MRI_UCHAR) ;
  MRIcopyVolGeomToMRI(mri_neg_orig, &gcam->image);

  mri_neg_atlas = MRIalloc(mri_warp->width, mri_warp->height, mri_warp->depth, MRI_UCHAR) ;
  MRIcopyVolGeomToMRI(mri_neg_atlas, &gcam->atlas);

  iter = 0 ;
  max_iter = 200 ;
  max_noprogress = 4 ;
  noprogress = 0 ;
  min_neg = gcam->neg ;
  max_nbhd = 4 ;

  nbhd = MAX(0,gcam->spacing-2) ;
  printf("iter %d, gcam->neg = %d\n", iter, gcam->neg) ;
  do
  {
    MRIclear(mri_neg) ;
    mri_warp_tmp = MRIcopy(mri_warp, mri_warp_tmp) ;
    last_neg = gcam->neg ;
    for (xp = 0 ; xp < gcam->width; xp++)
    {
      for (yp = 0 ; yp < gcam->height ; yp++)
      {
        for (zp = 0 ; zp < gcam->depth ; zp++)
        {
          if (xp == Gx && yp == Gy && zp == Gz)
          {
            DiagBreak() ;
          }
          xv = (xp * gcam->spacing) ;
          yv = (yp * gcam->spacing) ;
          zv = (zp * gcam->spacing) ;
          //        GCApriorToVoxel(gcam->gca, mri_warp, xp, yp, zp, &xv, &yv, &zv) ;
          gcamn = &gcam->nodes[xp][yp][zp] ;
          if (gcamn->area1 < 0 || gcamn->area2 < 0)
          {
            int xv, yv, zv ;
            xv = nint(gcamn->x) ;
            yv = nint(gcamn->y) ;
            zv = nint(gcamn->z) ;
            if (xv >= 0 && xv < mri_neg->width &&
                yv >= 0 && yv < mri_neg->height &&
                zv >= 0 && zv < mri_neg->depth)
            {
              MRIsetVoxVal(mri_neg, nint(gcamn->x), nint(gcamn->y), nint(gcamn->z), 0, 128) ;
            }

            xv = nint(gcamn->origx) ;
            yv = nint(gcamn->origy) ;
            zv = nint(gcamn->origz) ;
            if (xv >= 0 && xv < mri_neg->width &&
                yv >= 0 && yv < mri_neg->height &&
                zv >= 0 && zv < mri_neg->depth)
            {
              MRIsetVoxVal(mri_neg_orig, xv, yv, zv, 0, 128) ;
            }
            {
              MRIsetVoxVal(mri_neg_atlas, xp, yp, zp, 0, 128) ;
            }
          }

          if ((gcamn->area1 < 0 || gcamn->area2 < 0) && gcamn->invalid == GCAM_VALID)
          {
            for (xk = -nbhd ; xk <= nbhd ; xk++)
            {
              xv = xp + xk ;
              if (xv < 0 || xv >= mri_warp->width)
              {
                continue ;
              }
              for (yk = -nbhd ; yk <= nbhd ; yk++)
              {
                yv = yp + yk ;
                if (yv < 0 || yv >= mri_warp->height)
                {
                  continue ;
                }
                for (zk = -nbhd ; zk <= nbhd ; zk++)
                {
                  zv = zp + zk ;
                  if (zv < 0 || zv >= mri_warp->depth)
                  {
                    continue ;
                  }

                  dx = MRIvoxelMean(mri_warp_tmp, xv, yv, zv, wsize, 0) ;
                  dy = MRIvoxelMean(mri_warp_tmp, xv, yv, zv, wsize, 1);
                  dz = MRIvoxelMean(mri_warp_tmp, xv, yv, zv, wsize, 2) ;
                  MRIsetVoxVal(mri_warp, xv, yv, zv, 0, dx) ;
                  MRIsetVoxVal(mri_warp, xv, yv, zv, 1, dy) ;
                  MRIsetVoxVal(mri_warp, xv, yv, zv, 2, dz) ;
                }
              }
            }
          }
        }
      }
    }

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_neg, "neg.mgz") ;
      MRIwrite(mri_neg_atlas, "neg.atlas.mgz") ;
      MRIwrite(mri_neg_orig, "neg.orig.mgz") ;
    }
    GCAMreadWarpFromMRI(gcam, mri_warp) ;
    gcamComputeMetricProperties(gcam) ;
    if (Gdiag & DIAG_SHOW)
    {
      printf("iter %d, gcam->neg = %d, nbhd=%d\n", iter+1, gcam->neg, nbhd) ;
    }
    if (gcam->neg >= min_neg)
    {
      if (noprogress++ >= max_noprogress)
      {
        min_neg = gcam->neg+1 ; // reset min neg found
        nbhd++ ;
        if (nbhd > max_nbhd)  // try it starting with just nearest nbrs again
        {
          nbhd = 1 ;
#if 0
          if (max_nbhd < 7) // bigger than 7-nbrs takes forever
          {
            max_nbhd++ ;
          }
          else
#endif
            DiagBreak() ;
        }
        noprogress = 0 ;
        if (Gdiag & DIAG_SHOW)
        {
          printf("neg = %d, increasing nbhd size to %d\n", gcam->neg,nbhd) ;
        }
      }
    }
    else
    {
      noprogress = 0 ;
      min_neg = gcam->neg ;
    }
  }
  while (gcam->neg > 0 && (++iter < max_iter || gcam->neg < last_neg)) ;
  printf("after %d iterations, nbhd size=%d, neg = %d\n", iter, nbhd, gcam->neg) ;
  if (gcam->neg > 0)
  {
    DiagBreak() ;
  }
  MRIfree(&mri_warp_tmp) ;
  MRIfree(&mri_neg) ;
  MRIfree(&mri_neg_atlas) ;
  MRIfree(&mri_neg_orig) ;
  return(NO_ERROR) ;
}
MRI *
replace_labels(MRI *mri_src_labels, MRI *mri_dst_labels,
               int combine_labels[8][2], int ncombine_labels,
               GCA_MORPH_PARMS *parms)
{
  int n, *dtrans_labels, ndtrans, n2 ;

  mri_dst_labels = MRIcopy(mri_src_labels, mri_dst_labels) ;
  if (parms)
  {
    ndtrans = parms->ndtrans ;
  }
  else
  {
    ndtrans = 0 ;
  }
  for (n = 0 ; n < ncombine_labels ; n++)
  {
    printf("replacing label %s with %s\n",
           cma_label_to_name(combine_labels[n][0]),
           cma_label_to_name(combine_labels[n][1])) ;
    MRIreplaceValuesOnly(mri_src_labels, mri_dst_labels, combine_labels[n][0], combine_labels[n][1]) ;
    if (parms)
    {
      for (n2 = 0 ; n2 < parms->ndtrans ; n2++)
        if (parms->dtrans_labels[n2] == combine_labels[n][0])
        {
          parms->dtrans_labels[n2] = -1 ;  // remove it
          ndtrans-- ;
          break ;
        }
    }
  }

  if (parms)
  {
    dtrans_labels = (int *)calloc(ndtrans, sizeof(int)) ;
    for (n = n2 = 0 ; n < parms->ndtrans ; n++)
      if (parms->dtrans_labels[n] >= 0)
      {
        dtrans_labels[n2++] = parms->dtrans_labels[n] ;
      }
    parms->ndtrans = ndtrans ;
    parms->dtrans_labels = dtrans_labels ;
  }
  return(mri_dst_labels) ;
}

int
fix_borders(GCA_MORPH *gcam)
{
  int              x, y, z ;
  GCA_MORPH_NODE   *gcamn0, *gcamn1 ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
    {
      gcamn0 = &gcam->nodes[x][y][0] ;
      gcamn1 = &gcam->nodes[x][y][1] ;
      gcamn0->dx = gcamn1->dx ;
      gcamn0->dy = gcamn1->dy ;
      gcamn0->dz = gcamn1->dz ;

      gcamn0 = &gcam->nodes[x][y][gcam->depth-1] ;
      gcamn1 = &gcam->nodes[x][y][gcam->depth-2] ;
      gcamn0->dx = gcamn1->dx ;
      gcamn0->dy = gcamn1->dy ;
      gcamn0->dz = gcamn1->dz ;
    }


  for (x = 0 ; x < gcam->width ; x++)
    for (z = 0 ; z < gcam->depth ; z++)
    {
      gcamn0 = &gcam->nodes[x][0][z] ;
      gcamn1 = &gcam->nodes[x][1][z] ;
      gcamn0->dx = gcamn1->dx ;
      gcamn0->dy = gcamn1->dy ;
      gcamn0->dz = gcamn1->dz ;

      gcamn0 = &gcam->nodes[x][gcam->height-1][z] ;
      gcamn1 = &gcam->nodes[x][gcam->height-2][z] ;
      gcamn0->dx = gcamn1->dx ;
      gcamn0->dy = gcamn1->dy ;
      gcamn0->dz = gcamn1->dz ;
    }


  for (y = 0 ; y < gcam->height ; y++)
    for (z = 0 ; z < gcam->depth ; z++)
    {
      gcamn0 = &gcam->nodes[0][y][z] ;
      gcamn1 = &gcam->nodes[1][y][z] ;
      gcamn0->dx = gcamn1->dx ;
      gcamn0->dy = gcamn1->dy ;
      gcamn0->dz = gcamn1->dz ;

      gcamn0 = &gcam->nodes[gcam->width-1][y][z] ;
      gcamn1 = &gcam->nodes[gcam->width-2][y][z] ;
      gcamn0->dx = gcamn1->dx ;
      gcamn0->dy = gcamn1->dy ;
      gcamn0->dz = gcamn1->dz ;
    }

  return(NO_ERROR) ;
}
MRI *
GCAMestimateLameConstants(GCA_MORPH *gcam)
{
  MRI    *mri_lame,*mri_warp, *mri_sobel[3] ;
  int    x, y, z, dim1, dim2 ;
  double energy, val1, val2, lambda, mu, rigid_energy, volume_change ;

  lambda = 0.57692 ;
  mu = 0.38462 ; // matches CVS
  mri_warp = GCAMwriteWarpToMRI(gcam, NULL) ;
  mri_lame = MRIallocSequence(mri_warp->width, mri_warp->height, mri_warp->depth, MRI_FLOAT, 5) ;
  mri_sobel[0] = MRIsobelFrame(mri_warp, NULL, NULL, 0) ;
  mri_sobel[1] = MRIsobelFrame(mri_warp, NULL, NULL, 1) ;
  mri_sobel[2] = MRIsobelFrame(mri_warp, NULL, NULL, 2) ;

  for (x = 0 ; x < mri_warp->width ; x++)
    for (y = 0 ; y < mri_warp->height ; y++)
      for (z = 0 ; z < mri_warp->depth ; z++)
      {
        rigid_energy = volume_change = energy = 0 ;
        for (dim1 = 0 ; dim1 < 3 ; dim1++)
          for (dim2 = 0 ; dim2 < 3 ; dim2++)
          {
            val1 = MRIgetVoxVal(mri_sobel[dim1], x, y, z, dim1) ;
            val2 = MRIgetVoxVal(mri_sobel[dim2], x, y, z, dim2) ;
            rigid_energy += val1*val2 ;

            val1 = MRIgetVoxVal(mri_sobel[dim1], x, y, z, dim2) ;
            val2 = MRIgetVoxVal(mri_sobel[dim2], x, y, z, dim1) ;
            volume_change += SQR(val1 + val2);
          }
        energy = rigid_energy * mu/4 + volume_change * lambda/2 ;
        MRIsetVoxVal(mri_lame, x, y, z, 0, energy) ;
        MRIsetVoxVal(mri_lame, x, y, z, 1, rigid_energy) ;
        MRIsetVoxVal(mri_lame, x, y, z, 2, volume_change) ;
        MRIsetVoxVal(mri_lame, x, y, z, 3, lambda) ;
        MRIsetVoxVal(mri_lame, x, y, z, 4, mu) ;
      }

  for (dim1 = 0 ; dim1 < 3 ; dim1++)
  {
    MRIfree(&mri_sobel[dim1]) ;
  }
  MRIfree(&mri_warp) ;
  MRIcopyVolGeomToMRI(mri_lame, &gcam->atlas) ;
  return(mri_lame) ;
}

int
GCAMconcatenate(GCA_MORPH *gcam1, GCA_MORPH *gcam2, GCA_MORPH *gcam_composed)
{

  // printf("Entering GCAMconcatenate \n") ;

  // x, y, z are the col, row, and slice (and xyz) in the gcam/target volume
  int x, y, z, width, height, depth;
  int out_of_gcam1, out_of_gcam2;
  float xd, yd, zd, xdd, ydd, zdd ;
  GCA_MORPH_NODE  *gcamn, *gcam1n;
  if (gcam1->width != gcam2->width || gcam1->height != gcam2->height || gcam1->depth != gcam2->depth)
  {
    ErrorExit(ERROR_BADPARM, "GCAMconcatenate: the morphs to be concatenated need to be of the same size!\n") ;
  }

  if(! gcam_composed)
  {
    printf("Allocating gcam...(%d, %d, %d)\n", gcam1->width, gcam1->height, gcam1->depth);
    gcam_composed = GCAMalloc(gcam1->width, gcam1->height, gcam1->depth) ;
  }

  // get geometry information
  gcam_composed->image = gcam1->image;

  width  = gcam_composed->width ;
  height = gcam_composed->height ;
  depth  = gcam_composed->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        // Should not divide by src thick
        //out_of_gcam = GCAMsampleMorph(gcam, (float)x*mri_src->thick,
        //   (float)y*mri_src->thick,
        //   (float)z*mri_src->thick,
        //   &xd, &yd, &zd);

        // Convert target-crs to input-crs
        out_of_gcam2 = GCAMsampleMorph(gcam2,
                                       (float)x, (float)y, (float)z,
                                       &xd, &yd, &zd);

        if (!out_of_gcam2)
        {
          out_of_gcam1 = GCAMsampleMorph(gcam1,
                                         (float)xd, (float)yd, (float)zd,
                                         &xdd, &ydd, &zdd);
          if (!out_of_gcam1)
          {

            gcamn = &gcam_composed->nodes[x][y][z] ;
            gcam1n = &gcam2->nodes[x][y][z] ;
            if (nint(xdd) == Gx && nint(ydd) == Gy && nint(zdd) == Gz)
            {
              DiagBreak() ;
            }

            if (xdd > -1 && ydd > -1 && zdd > 0 &&
                xdd < gcam1->image.width &&
                ydd < gcam1->image.height &&
                zdd < gcam1->image.depth)
            {

              gcamn->label = gcam1n->label;

              gcamn->origx =  xdd;
              gcamn->origy =  ydd;
              gcamn->origz =  zdd;

              gcamn->x =  xdd;
              gcamn->y =  ydd;
              gcamn->z =  zdd;

              gcamn->xn =  xdd;
              gcamn->yn =  ydd;
              gcamn->zn =  zdd;
            }

          }
        }
      }
    }
  }

  // printf("Exiting GCAMconcatenate \n") ;
  return(NO_ERROR) ;
}

// >>>>>>>>>>>>>>>>>>

GCA_MORPH *
GCAMfillInverse(GCA_MORPH* gcam)
{
  MRI  *mri;
  char tmpstr[2000];
  int width, height, depth;
  int x, y, z;
  GCA_MORPH_NODE  *gcamn, *invgcamn;
  GCA_MORPH* inv_gcam;


  printf("Allocating inv_gcam...(%d, %d, %d)\n", gcam->width, gcam->height, gcam->depth);
  inv_gcam = GCAMalloc(gcam->width, gcam->height, gcam->depth) ; // NOTE: forces same moving and target coordinate spaces!!
  inv_gcam->image = gcam->atlas;
  inv_gcam->atlas = gcam->image;

  sprintf(tmpstr, "%s", (gcam->image).fname);
  mri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
  if (mri==NULL)
  {
    printf("ERROR: reading %s\n",tmpstr);
    return(NULL);
  }
  if (( gcam->mri_xind == NULL) ||  ( gcam->mri_yind == NULL)  || ( gcam->mri_zind == NULL) )
  {
    // Must invert explicitly
    printf("GCAMfillInverse: Must invert gcam explicitely! \n");
    GCAMinvert(gcam, mri);
  }

  ////////
  width  = inv_gcam->width ;
  height = inv_gcam->height ;
  depth  = inv_gcam->depth ;


  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }

        gcamn = &gcam->nodes[x][y][z] ;
        invgcamn = &inv_gcam->nodes[x][y][z] ;
        // missing: all the checks

        invgcamn->origx =  MRIgetVoxVal(gcam->mri_xind, x, y, z, 0) ;
        invgcamn->origy =  MRIgetVoxVal(gcam->mri_yind, x, y, z, 0) ;
        invgcamn->origz =  MRIgetVoxVal(gcam->mri_zind, x, y, z, 0) ;

        invgcamn->x =  MRIgetVoxVal(gcam->mri_xind, x, y, z, 0) ;
        invgcamn->y =  MRIgetVoxVal(gcam->mri_yind, x, y, z, 0) ;
        invgcamn->z =  MRIgetVoxVal(gcam->mri_zind, x, y, z, 0) ;

        invgcamn->xn =  MRIgetVoxVal(gcam->mri_xind, x, y, z, 0) ;
        invgcamn->yn =  MRIgetVoxVal(gcam->mri_yind, x, y, z, 0) ;
        invgcamn->zn =  MRIgetVoxVal(gcam->mri_zind, x, y, z, 0) ;

        invgcamn->label = gcamn->label;
      }
    }
  }
  return(inv_gcam) ;
}
GCA_MORPH *
GCAMdownsample2(GCA_MORPH *gcam)
{
  int            xs, ys, zs, xd, yd, zd, labels[MAX_CMA_LABEL+1], l, max_l, max_count ;
  GCA_MORPH_NODE *gcamn_src, *gcamn_dst ;
  GCA_MORPH      *gcam_dst ;

  gcam_dst = GCAMalloc(gcam->width/2, gcam->height/2, gcam->depth/2) ;
  *(&gcam_dst->image) = *(&gcam->image) ;
  *(&gcam_dst->atlas) = *(&gcam->atlas) ;
  gcam_dst->spacing = 2*gcam->spacing ;
  gcam_dst->ninputs = gcam->ninputs ;
  gcam_dst->gca = gcam->gca ;
  gcam_dst->status = gcam->status ;
  gcam_dst->type = gcam->type ;
  gcam_dst->m_affine = gcam->m_affine ;
  gcam_dst->det = gcam->det ;

  for (xd = 0 ; xd < gcam_dst->width ; xd++)
  {
    for (yd = 0 ; yd < gcam_dst->height ; yd++)
    {
      for (zd = 0 ; zd < gcam_dst->depth ; zd++)
      {
        gcamn_dst = &gcam->nodes[xd][yd][zd] ;
        memset(labels, 0, sizeof(labels)) ;

        for (xs = xd*2 ; xs <= xd*2+1 ; xs++)
          for (ys = yd*2 ; ys <= yd*2+1 ; ys++)
            for (zs = zd*2 ; zs <= zd*2+1 ; zs++)
            {
              gcamn_src = &gcam->nodes[xs][ys][zs] ;
              labels[gcamn_src->label]++ ;

              gcamn_dst->x += gcamn_src->x ;
              gcamn_dst->y += gcamn_src->y ;
              gcamn_dst->z += gcamn_src->z ;

              gcamn_dst->origx += gcamn_src->origx ;
              gcamn_dst->origy += gcamn_src->origy ;
              gcamn_dst->origz += gcamn_src->origz ;

              gcamn_dst->xs2 += gcamn_src->xs2 ;
              gcamn_dst->ys2 += gcamn_src->ys2 ;
              gcamn_dst->zs2 += gcamn_src->zs2 ;

              gcamn_dst->xs += gcamn_src->xs ;
              gcamn_dst->ys += gcamn_src->ys ;
              gcamn_dst->zs += gcamn_src->zs ;

              gcamn_dst->xn += gcamn_src->xn ;
              gcamn_dst->yn += gcamn_src->yn ;
              gcamn_dst->zn += gcamn_src->zn ;

              gcamn_dst->saved_origx += gcamn_src->saved_origx ;
              gcamn_dst->saved_origy += gcamn_src->saved_origy ;
              gcamn_dst->saved_origz += gcamn_src->saved_origz ;

              gcamn_dst->prior += gcamn_src->prior ;
              gcamn_dst->area += gcamn_src->area ;
              gcamn_dst->area1 += gcamn_src->area1 ;
              gcamn_dst->area2 += gcamn_src->area2 ;
              gcamn_dst->orig_area1 += gcamn_src->orig_area1 ;
              gcamn_dst->orig_area2 += gcamn_src->orig_area2 ;
              if (gcamn_src->invalid)
              {
                gcamn_dst->invalid = 1 ;
              }
              if (gcamn_src->status > 0)
              {
                gcamn_dst->status = gcamn_src->status ;
              }
            }

        gcamn_dst->x /= 8 ;
        gcamn_dst->y /= 8 ;
        gcamn_dst->z /= 8 ;

        gcamn_dst->origx /= 8 ;
        gcamn_dst->origy /= 8 ;
        gcamn_dst->origz /= 8 ;

        gcamn_dst->xs2 /= 8 ;
        gcamn_dst->ys2 /= 8 ;
        gcamn_dst->zs2 /= 8 ;

        gcamn_dst->xs /= 8 ;
        gcamn_dst->ys /= 8 ;
        gcamn_dst->zs /= 8 ;

        gcamn_dst->xn /= 8 ;
        gcamn_dst->yn /= 8 ;
        gcamn_dst->zn /= 8 ;


        max_count = labels[0] ;
        max_l = 0 ;
        for (l = 1 ; l <= MAX_CMA_LABELS ; l++)
          if (labels[l] > max_count)
          {
            max_count = labels[l] ;
            max_l = l ;
          }
        gcamn_dst->label = max_l ;
      }
    }
  }
  return(gcam_dst) ;
}


int
GCAMdilateUseLikelihood(GCA_MORPH *gcam, int ndilations)
{
  MRI             *mri ;
  int             x, y, z, i ;
  GCA_MORPH_NODE *gcamn ;

  mri = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_UCHAR) ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->status & GCAM_IGNORE_LIKELIHOOD ||
            gcamn->status & GCAM_NEVER_USE_LIKELIHOOD)
        {
          continue ;
        }
        MRIsetVoxVal(mri, x, y, z, 0, 1) ;
      }
  for (i = 0 ; i < ndilations ; i++)
  {
    MRIdilate(mri, mri) ;
  }
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->status & GCAM_IGNORE_LIKELIHOOD ||
            gcamn->status & GCAM_NEVER_USE_LIKELIHOOD)
        {
          if (MRIgetVoxVal(mri, x, y, z, 0))
          {
            gcamn->status = GCAM_USE_LIKELIHOOD ;
          }
        }
      }

  MRIfree(&mri) ;
  return(NO_ERROR) ;
}



static MRI *
gcamLabelToMRI(GCA_MORPH *gcam, MRI *mri, int label)
{
  int             x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  if (mri == NULL)
  {
    mri = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_SHORT) ;
//    useVolGeomToMRI(&gcam->atlas, mri) ;
  }

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (gcamn->label == label)
        {
          MRIsetVoxVal(mri, x, y, z, 0, label) ;
        }
      }

  return(mri) ;
}

int
GCAMcomputeVentricleExpansionGradient(GCA_MORPH *gcam)
{
  int             x, y, z ;
  GCA_MORPH_NODE *gcamn ;
  MRI            *mri_vent_dist, *mri_vent, *mri_vent_dist_grad ;

  mri_vent = gcamLabelToMRI(gcam, NULL, Left_Lateral_Ventricle) ;
  gcamLabelToMRI(gcam, mri_vent, Right_Lateral_Ventricle) ;
  gcamLabelToMRI(gcam, mri_vent, Right_Inf_Lat_Vent) ;
  gcamLabelToMRI(gcam, mri_vent, Left_Inf_Lat_Vent) ;
  gcamLabelToMRI(gcam, mri_vent, Third_Ventricle) ;
  gcamLabelToMRI(gcam, mri_vent, Fourth_Ventricle) ;
  MRIbinarize(mri_vent, mri_vent, 1, 0, 1) ;
  mri_vent_dist =  MRIdistanceTransform(mri_vent, NULL, 1, -1, DTRANS_MODE_SIGNED, NULL);
  mri_vent_dist_grad = MRIsobel(mri_vent_dist, NULL, NULL) ;

  MRIfree(&mri_vent) ;
  MRIfree(&mri_vent_dist) ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        gcamn = &gcam->nodes[x][y][z] ;
        if (IS_VENTRICLE(gcamn->label))
        {
          if (gcamn->label == Left_Lateral_Ventricle)
          {
//      Gx = x ; Gy= y ; Gz=z;
            DiagBreak() ;
          }
          gcamn->dx = MRIgetVoxVal(mri_vent_dist_grad, x, y, z, 0) ;
          gcamn->dy = MRIgetVoxVal(mri_vent_dist_grad, x, y, z, 1) ;
          gcamn->dz = MRIgetVoxVal(mri_vent_dist_grad, x, y, z, 2) ;
        }
      }

  MRIfree(&mri_vent_dist_grad) ;
  return(NO_ERROR) ;
}

int
GCAMregisterVentricles(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  char   fname[STRLEN] ;
  int    navgs, which, fixed_steps = 0, n ;
  double base_sigma, pct_change, rms, last_rms = 0.0, orig_dt, l_orig_smooth, min_dt, l_orig_elastic, l_elastic, tol ;

  if (FZERO(parms->min_sigma))
  {
    parms->min_sigma = 0.4 ;
  }

  navgs = parms->navgs ;
  orig_dt = parms->dt ;
  l_orig_smooth = parms->l_smoothness;
  l_orig_elastic = l_elastic = parms->l_elastic;
  if (DZERO(parms->exp_k))
  {
    parms->exp_k = EXP_K ;
  }
  if (parms->levels < 0)
  {
    parms->levels = DEFAULT_PYRAMID_LEVELS ;
  }
  else if (parms->levels >= MAX_PYRAMID_LEVELS)
  {
    parms->levels = MAX_PYRAMID_LEVELS ;
  }

  gcam->exp_k = parms->exp_k ;
  parms->mri = mri ;
  //////////////////////////////////////////////////////////
  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s.log", parms->base_name) ;
    if (parms->log_fp == NULL)
    {
      if (parms->start_t == 0)
      {
        parms->log_fp = fopen(fname, "w") ;
      }
      else
      {
        parms->log_fp = fopen(fname, "a") ;
      }
    }
    log_integration_parms(parms->log_fp, parms) ;
  }
  else
  {
    parms->log_fp = NULL ;
  }

  if (Gdiag & DIAG_SHOW)
  {
    log_integration_parms(stdout, parms) ;
  }
  base_sigma = parms->sigma ;

  // GCAMinit() did this at the end
  // make node to have the max_prior label values
  if (parms->relabel_avgs >= parms->navgs && parms->relabel)
  {
    GCAMcomputeLabels(mri, gcam) ;
  }
  else
  {
    GCAMcomputeMaxPriorLabels(gcam) ;
  }

  tol = parms->tol ;
  if (parms->write_iterations && (Gdiag & DIAG_WRITE) && parms->start_t == 0)
  {
    write_snapshot(gcam, mri, parms, 0) ;
  }
  n = 0 ;
  for (navgs = parms->navgs ; parms->navgs >= 0 ; parms->navgs /= 4)
  {
    parms->l_smoothness = l_orig_smooth / (sqrt(parms->navgs)+1) ;
    if (parms->integration_type == GCAM_INTEGRATE_OPTIMAL || parms->integration_type == GCAM_INTEGRATE_BOTH)
    {
      which = GCAM_INTEGRATE_OPTIMAL ;
    }
    else
    {
      which = GCAM_INTEGRATE_FIXED ;
    }

    last_rms = GCAMcomputeRMS(gcam, mri, parms) ;
    if (FZERO(parms->start_t) && n==0)
    {
      if (parms->log_fp)
      {
        fprintf(parms->log_fp, "%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, PA, navgs=%d",
                parms->start_t, 0.0, last_rms, 0.0, gcam->neg, Ginvalid,navgs) ;
        if (parms->l_binary > 0)
          fprintf(parms->log_fp, ", aligned = %d (%2.3f%%)\n",
                  Galigned, 100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          fprintf(parms->log_fp, "\n") ;
        }
        fflush(parms->log_fp) ;
      }

      if (Gdiag & DIAG_SHOW)
      {
        printf("%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, PA, navgs=%d",
               parms->start_t, 0.0, last_rms, 0.0, gcam->neg, Ginvalid, navgs) ;
        if (parms->l_binary > 0)
          printf(", aligned = %d (%2.3f%%)\n",
                 Galigned, 100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          printf("\n") ;
        }
      }
    }
    do
    {
      GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, SAVED2_POSITIONS) ;
      gcamClearGradient(gcam) ;
      GCAMcomputeVentricleExpansionGradient(gcam) ;
      gcamSmoothnessTerm(gcam, mri, parms->l_smoothness)  ;
      gcamLSmoothnessTerm(gcam, mri, parms->l_lsmoothness)  ;
      gcamSpringTerm(gcam, parms->l_spring, parms->ratio_thresh)  ;
      gcamJacobianTerm(gcam, mri, parms->l_jacobian,parms->ratio_thresh)  ;
      // The following appears to be a null operation, based on current #ifdefs
      gcamLimitGradientMagnitude(gcam, parms, mri) ;
      gcamSmoothGradient(gcam, parms->navgs) ;


      if (which == GCAM_INTEGRATE_OPTIMAL)
      {
        parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ; /* will search around
                 this value */
        min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
        if (FZERO(min_dt)) // try a momentum step
        {
          fixed_steps++ ;
          which = GCAM_INTEGRATE_FIXED ;
          min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
        }
        else
        {
          fixed_steps = 0 ;
        }
        parms->dt = min_dt ;
      }
      else
      {
        fixed_steps++ ;
        min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
      }

      gcamApplyGradient(gcam, parms) ;
      gcamComputeMetricProperties(gcam) ;
      gcamCheck(gcam, mri) ;
      if (parms->constrain_jacobian)
      {
        gcamConstrainJacobian(gcam, mri, parms) ;
      }
      if (Gdiag & DIAG_SHOW)
      {
        gcamShowCompressed(gcam, stdout) ;
      }

      gcamRemoveNegativeNodes(gcam, mri, parms) ;
      if (gcam->neg > 0)   // couldn't unfold everything
      {
        printf("---------- unfolding failed - restoring original position --------------------\n");
        GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
        gcamComputeMetricProperties(gcam) ;
      }
      if (parms->uncompress)
      {
        gcamRemoveCompressedNodes(gcam, mri, parms, parms->ratio_thresh) ;
      }
      rms = GCAMcomputeRMS(gcam, mri, parms) ;
      pct_change = 100.0*(last_rms-rms)/last_rms ;
      if (pct_change < 0)
      {
        printf("rms increased (dt=%2.4f) - undoing step...\n", min_dt) ;
        GCAMcopyNodePositions(gcam, SAVED2_POSITIONS, CURRENT_POSITIONS) ;
        gcamComputeMetricProperties(gcam) ;
        gcamClearMomentum(gcam) ;
        continue ;
      }
      if (pct_change < 2*tol && parms->integration_type == GCAM_INTEGRATE_BOTH && fixed_steps > 5)
      {
        printf("switching integration type back to optimal\n") ;
        which = GCAM_INTEGRATE_OPTIMAL ;
      }

      if (!FZERO(parms->start_t))
      {
        DiagBreak() ;
      }
      if (parms->write_iterations > 0 &&
          !((n+1) % parms->write_iterations) && (Gdiag & DIAG_WRITE))
      {
        write_snapshot(gcam, mri, parms, parms->start_t+n+1) ;
      }


      last_rms = rms ;
      if (parms->log_fp)
      {
        fprintf(parms->log_fp, "%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, PA, navgs=%d",
                parms->start_t+n+1, min_dt, rms, pct_change, gcam->neg,
                Ginvalid,parms->navgs) ;
        if (parms->l_binary > 0)
          fprintf(parms->log_fp, ", aligned = %d (%2.3f%%)\n",
                  Galigned,
                  100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          fprintf(parms->log_fp, "\n") ;
        }
        fflush(parms->log_fp) ;
      }

      if (Gdiag & DIAG_SHOW)
      {
        printf("%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d, PA, navgs=%d",
               parms->start_t+n+1, min_dt, rms, pct_change, gcam->neg, Ginvalid, parms->navgs) ;
        if (parms->l_binary > 0)
          printf(", aligned = %d (%2.3f%%)\n",
                 Galigned,
                 100.0*Galigned/((gcam->width*gcam->height*gcam->depth)-Ginvalid));
        else
        {
          printf("\n") ;
        }
      }
      n++ ;
    }
    while (pct_change > tol) ;
    if (parms->navgs <= parms->min_avgs)
    {
      break ;
    }
  }

  parms->navgs = navgs ;
  parms->sigma = base_sigma ;
  parms->l_smoothness = l_orig_smooth ;
  parms->start_t += n ;

  if (parms->log_fp)
  {
    fclose(parms->log_fp) ;
    parms->log_fp = NULL ;
  }

  return(NO_ERROR) ;
}

#if 0
static int
GCAMsetNegativeNodeStatus(GCA_MORPH *gcam, int status)
{
  int             x, y, z ;
  GCA_MORPH_NODE *gcamn ;
  double         ratio ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
	gcamn = &gcam->nodes[x][y][z] ;
	if (gcamn->area1 <= 0 || gcamn->area2 <= 0)
	  gcamn->status |= status ;
	if (FZERO(gcamn->orig_area))
	  continue ;
        ratio = gcamn->area / gcamn->orig_area ;
	if (ratio <  .1)
	  gcamn->status |= status ;
      }
  return(NO_ERROR) ;
}
#endif
