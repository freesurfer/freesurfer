/**
 * @brief utilities for boundary deformations
 *
 * Reference:
 * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
 * Structures in the Human Brain", Fischl et al.
 * (2002) Neuron, 33:341-355.
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


#ifndef GCABOUNDARY_H
#define GCABOUNDARY_H

#include "gca.h"
#include "icosahedron.h"
#include "const.h"
#include "histo.h"

typedef struct
{
  double       ***pdfs ;        // one per vertex, and each is 2d (I_in, I_out )
  HISTOGRAM    **h_grad_pdfs ;  // one per vertex, gradient histogram
  HISTOGRAM    **h_Iin_pdfs ;   // one per vertex, interior intensity histogram
  HISTOGRAM    **h_Iout_pdfs ;  // one per vertex, exterior intensity histogram
  float        *ntraining ;     // one per vertex
} GCA_BOUNDARY_STATS, GCABS ;

typedef struct
{
  int                  target_label ;    // what label is this the model of the boundary of
  int                  spacing ;
  int                  nvertices ;       // # of vertices in the icosahedron
  int                  ico_order ;       // redundant with previous (but makes things simpler)
  int                  nintensity_bins ;
  int                  ngrad_bins ;
  int                  width ;
  int                  height ;
  int                  depth ;
  float                max_intensity ;
  float                min_intensity ;
  float                min_grad ;
  float                max_grad ;
  float                iscale, gscale ;  // for histograming
  int                  ntraining ;
  MRI_REGION           bounding_box ;

  GCA                  *gca ;
  GCA_BOUNDARY_STATS   ***bs ;
  ICOSAHEDRON         *ico ;  // for normal directions
  char                fname[STRLEN] ;
  char                 gca_fname[STRLEN] ;
} GCA_BOUNDARY, GCAB ;


double GCABgetProbability(GCAB *gcab, MRI *mri_int, MRI *mri_dist, TRANSFORM *transform,
                          float x0, float y0, float z0,
                          float nx, float ny, float nz) ;
double GCABgetPout(GCAB *gcab, MRI *mri_int, MRI *mri_dist, TRANSFORM *transform,
                   float x0, float y0, float z0,
                   float nx, float ny, float nz) ;
double GCABgetPin(GCAB *gcab, MRI *mri_int, MRI *mri_dist, TRANSFORM *transform,
                  float x0, float y0, float z0,
                  float nx, float ny, float nz) ;
double GCABgetPgrad(GCAB *gcab, MRI *mri_int, MRI *mri_dist, TRANSFORM *transform,
                    float x0, float y0, float z0,
                    float nx, float ny, float nz) ;
GCAB *GCABalloc(GCA *gca, int spacing, int ico_order,
                int intensity_bins, int grad_bins, int target_label) ;

int  GCABfree(GCAB **pgcab) ;
int  GCABwrite(GCAB *gcab, char *fname) ;
GCAB *GCABread(char *fname, GCA *gca) ;
int  GCABtrain(GCAB *gcab, MRI *mri_int, MRI *mri_seg, TRANSFORM *transform,
               int label);
int  GCABcompleteTraining(GCAB *gcab) ;
int  GCABsourceVoxelToNode(GCAB *gcab, MRI *mri, TRANSFORM *transform,
                           float xv, float yv, float zv,
                           float *pxn, float *pyn, float *pzn);
int  GCABvoxelToNode(GCAB *gcab, float xv, float yv, float zv,
                     float *pxn, float *pyn, float *pzn) ;
int GCABsmoothPDFs(GCAB *gcab, float sigma) ;


#define GCAB_VERSION 1.0

#endif
