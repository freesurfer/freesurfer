/**
 * @brief include file for ROI utilies.
 *
 * structures, macros, constants and prototypes for the 
 * manipulation, creation and 
 * I/O for labels (lists of vertices and/or voxels)
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

#ifndef LABEL_H
#define LABEL_H

#include "minc.h"

#include "matrix.h"
#include "const.h"
#include "mri.h"

#include "mrisurf_aaa.h"    // incomplete MRIS, MHT
#include "transform.h"

struct LABEL_VERTEX
{
  int           vno ;
  float         x ;
  float         y ;
  float         z ;
  int           xv ;  // voxel coords into label->mri struct (may not be initialzed if label->mri == NULL)
  int           yv ;
  int           zv ;
  unsigned char deleted ;
  float         stat ;     /* statistic (might not be used) */
};

struct LABEL
{
  int    max_points ;         /* # of points allocated */
  int    coords ;             // one of the LABEL_COORDS* constants below
  int    n_points ;           /* # of points in area */
  char   name[STRLEN] ;       /* name of label file */
  char   subject_name[STRLEN] ;  /* name of subject */
  LV     *lv ;                /* labeled vertices */
  General_transform transform ;   /* the next two are from this struct */
  Transform         *linear_transform ;
  Transform         *inverse_linear_transform ;
  char   space[100];          /* space description of the coords */
  double avg_stat ;
  int    *vertex_label_ind ; // mris->nvertices long - < 0 means it isn't in the label
  MRI    *mri_template ;
  MHT    *mht ;
  MRIS   *mris ; 
};

#define LABEL_COORDS_NONE         FS_COORDS_UNKNOWN
#define LABEL_COORDS_TKREG_RAS    FS_COORDS_TKREG_RAS
#define LABEL_COORDS_SCANNER_RAS  FS_COORDS_SCANNER_RAS
#define LABEL_COORDS_VOXEL        FS_COORDS_VOXEL
//#define LABEL_COORDS_SURFACE_RAS  4 // same as TKREG_RAS


LABEL *LabelToScannerRAS(LABEL *lsrc, MRI *mri, LABEL *ldst) ;
LABEL *LabelToVoxel(LABEL *lsrc, MRI *mri, LABEL *ldst) ;
LABEL *LabelVoxelToSurfaceRAS(LABEL *lsrc, MRI *mri, LABEL *ldst) ;
LABEL *LabelToSurfaceRAS(LABEL *lsrc, MRI *mri, LABEL *ldst) ;
#define LabelFromScannerRAS LabelToSurfaceRAS   

int     LabelIsCompletelyUnassigned(LABEL *area, int *unassigned);
int     LabelFillUnassignedVertices(MRI_SURFACE *mris,
                                    LABEL *area,
                                    int coords);
double  LabelMeanIntensity(LABEL *area, MRI *mri) ;
int     LabelFree(LABEL **parea) ;
int     LabelDump(FILE *fp, LABEL *area) ;
LABEL   *LabelRead(const char *subject_name,const char *label_name) ;
LABEL   *LabelReadFrom(const char *subject_name, FILE *fp) ;
int     LabelWriteInto(LABEL *area, FILE *fp) ;
int     LabelWrite(LABEL *area,const char *fname) ;
int     LabelToCurrent(LABEL *area, MRI_SURFACE *mris) ;
int     LabelToCanonical(LABEL *area, MRI_SURFACE *mris) ;
int     LabelThreshold(LABEL *area, float thresh) ;
int     LabelMarkWithThreshold(LABEL *area, MRI_SURFACE *mris, float thresh);
int     LabelMarkSurface(LABEL *area, MRI_SURFACE *mris) ;
int     LabelAddToSurfaceMark(LABEL *area, MRI_SURFACE *mris, int mark_to_add)  ;
int     LabelToOriginal(LABEL *area, MRI_SURFACE *mris) ;
int     LabelToWhite(LABEL *area, MRI_SURFACE *mris) ;
int     LabelFromCanonical(LABEL *area, MRI_SURFACE *mris) ;
int     LabelFromTalairach(LABEL *area, MRI_SURFACE *mris) ;
int     LabelToFlat(LABEL *area, MRI_SURFACE *mris) ;
int     LabelRipRestOfSurface(LABEL *area, MRI_SURFACE *mris) ;
int     LabelRipRestOfSurfaceWithThreshold(LABEL *area, MRI_SURFACE *mris, float thresh) ;
int     LabelRip(MRI_SURFACE *mris, const LABEL *area, const int Outside);
int     LabelRemoveOverlap(LABEL *area1, LABEL *area2) ;
int     LabelIntersect(LABEL *area1, LABEL *area2) ;
LABEL  *LabelRemoveAlmostDuplicates(LABEL *area, double dist, LABEL *ldst);
LABEL   *LabelCompact(LABEL *lsrc, LABEL *ldst) ;
int     LabelRemoveDuplicates(LABEL *area) ;
int     LabelHasVertex(int vtxno, LABEL *lb);
LABEL   *LabelAlloc(int max_points, const char *subject_name, const char *label_name) ;
LABEL   *LabelRealloc(LABEL *lb, int max_points);
int     LabelCurvFill(LABEL *area, int *vertex_list, int nvertices,
                      int max_vertices, MRI_SURFACE *mris) ;
int     LabelFillMarked(LABEL *area, MRI_SURFACE *mris) ;
int     LabelFillAnnotated(LABEL *area, MRI_SURFACE *mris) ;
int     LabelFillAll(LABEL *area, int *vertex_list, int nvertices,
                     int max_vertices, MRI_SURFACE *mris) ;
int     LabelTalairachTransform(LABEL *area, MRI_SURFACE *mris) ;
int     LabelSphericalTransform(LABEL *area, MRI_SURFACE *mris) ;
MATRIX  *LabelCovarianceMatrix(LABEL *area, MATRIX *mat) ;
LABEL   *LabelCombine(LABEL *area, LABEL *area_dst) ;

LABEL   *LabelTranslate(LABEL *area,
                        LABEL *area_offset,
                        float dx, float dy, float dz) ;
LABEL   *LabelCopy(LABEL *asrc, LABEL *adst) ;
LABEL   *LabelCombine(LABEL *area, LABEL *adst) ;
double  LabelArea(LABEL *area, MRI_SURFACE *mris) ;
double  LabelVariance(LABEL *area, double ux, double uy, double uz) ;
int     LabelMean(LABEL *area, double *px, double *py, double *pz) ;
int     LabelMark(LABEL *area, MRI_SURFACE *mris) ;
int     LabelMark2(LABEL *area, MRI_SURFACE *mris);
int     LabelMarkUndeleted(LABEL *area, MRI_SURFACE *mris) ;
float   LabelMaxStat(LABEL *area) ;
int     LabelMarkStats(LABEL *area, MRI_SURFACE *mris) ;
int     LabelUnmark(LABEL *area, MRI_SURFACE *mris) ;
LABEL   *LabelFromMarkedSurface(MRI_SURFACE *mris) ;
LABEL   *LabelFromMarkValue(MRI_SURFACE *mris, int mark);
int     LabelNormalizeStats(LABEL *area, float norm) ;
LABEL   *MaskSurfLabel(LABEL *lbl,
                       MRI *SurfMask,
                       float thresh, const char *masksign, int frame);

int     LabelErode(LABEL *area, MRI_SURFACE *mris, int num_times);
int     LabelDilate(LABEL *area, MRI_SURFACE *mris, int num_times, int coords);

int   LabelSetStat(LABEL *area, float stat) ;
int   LabelCopyStatsToSurface(LABEL *area, MRI_SURFACE *mris, int which) ;
LABEL *LabelFillHoles(LABEL *area_src, MRI_SURFACE *mris, int coords) ;
LABEL *LabelFillHolesWithOrig(LABEL *area_src, MRI_SURFACE *mris) ;
LABEL *LabelfromASeg(MRI *aseg, int segcode);
int   LabelFillVolume(MRI *mri, LABEL *label, int fillval) ;

MATRIX *LabelFitXYZ(LABEL *label, int order);
LABEL *LabelApplyMatrix(LABEL *lsrc, MATRIX *m, LABEL *ldst) ;
LABEL *LabelBoundary(LABEL *label, MRIS *surf);
int VertexIsInLabel(int vtxno, LABEL *label);
LABEL *LabelInFOV(MRI_SURFACE *mris, MRI *mri, float pad) ;
int   LabelUnassign(LABEL *area) ;
LABEL *MRISlabelInvert(MRIS *surf, LABEL *label);
int LabelMaskSurface(LABEL *label, MRI_SURFACE *mris) ;
int LabelMaskSurfaceValues(LABEL *label, MRI_SURFACE *mris) ;
int LabelMaskSurfaceCurvature(LABEL *label, MRI_SURFACE *mris) ;
int LabelMaskSurfaceVolume(LABEL *label, MRI *mri, float nonmask_val) ;

LABEL   *LabelSphericalCombine(MRI_SURFACE *mris, LABEL *area,
                               MRIS_HASH_TABLE *mht,
                               MRI_SURFACE *mris_dst, LABEL *area_dst);

LABEL *LabelClone(LABEL *a)  ;
int LabelCropPosterior(LABEL *area, float anterior_dist) ;
int LabelCropAnterior(LABEL *area, float anterior_dist) ;
int LabelCentroid(LABEL *area, MRI_SURFACE *mris, double *px, double *py, double *pz, int *pvno) ;
int LabelSetVals(MRI_SURFACE *mris, LABEL *area, float fillval) ;
int LabelAddToMark(LABEL *area, MRI_SURFACE *mris, int val_to_add) ;
LABEL *LabelTransform(LABEL *area_in, TRANSFORM *xform, MRI *mri, LABEL *area_out) ;
LABEL *LabelFromScannerRAS(LABEL *lsrc, MRI *mri, LABEL *ldst) ;
LABEL *LabelBaryFill(MRIS *mris, LABEL *srclabel, double delta);

LABEL *LabelSampleToSurface(MRI_SURFACE *mris, LABEL *area, MRI *mri_template, int coords) ;
int   LabelInit(LABEL *lsrc, MRI *mri_template, MRI_SURFACE *mris, int coords) ;
int   LabelAddVoxel(LABEL *area, int xv, int yv, int zv, int coords, int *vertices, int *pnvertices) ;
int   LabelDeleteVoxel(LABEL *area, int xv, int yv, int zv, int *vertices, int *pnvertices) ;
LABEL *LabelAddPoint(LABEL *label, LV *lv);
int   LabelAddVertex(LABEL *area, int vno, int coords) ;
int   LabelDeleteVertex(LABEL *area, int vno, int coords) ;
double LabelAverageVal(LABEL *area, MRI_SURFACE *mris) ;
LABEL  *LabelFromSurface(MRI_SURFACE *mris, int which, double thresh) ;
LABEL *LabelRemoveIslandsSurf(MRIS *surf, LABEL *lb);
LABEL *LabelRemoveHolesSurf(MRIS *surf, LABEL *lb);
LABEL *LabelRemoveHolesAndIslandsSurf(MRIS *surf, LABEL *lb);

int LabelVertexTableSortByVtxno(const void *p1, const void *p2);

#endif
