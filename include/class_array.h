/**
 * @brief prototypes and structs for using an array of classifiers to reclassify aseg voxels
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

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#ifndef CLASS_ARRAY_H
#define CLASS_ARRAY_H

#include "mri.h"
#include "voxlist.h"
#include "icosahedron.h"
#include "svm.h"

typedef struct
{
  int    type ;
  SVM    *svm ;
  double *c1_means ;
  double *c1_vars ;
  double *c2_means ;
  double *c2_vars ;
  int    c1_ntraining ;
  int    c2_ntraining ;
} CLASSIFIER ;


// c1 > 0, c2 < 0 in classification
typedef struct
{
  int        width ;
  int        height ;
  int        depth ;
  int        ninputs ;
  int        nscales ;
  int        wsize ;
  int        which_inputs ;
  MATRIX     *m_vox2index ;
  CLASSIFIER ****classifiers ;  // x, y, z, orientation
  char       c1_name[STRLEN] ;
  char       c2_name[STRLEN] ;
  float      *sigmas ;  // nscales long array of blurring sigmas
  double     svm_C ;
  double     tol ;
  int        max_iter ;
  int        type ;    // what kind of classifiers
  int        c1_label ;
  int        c2_label ;
} CLASSIFIER_ATLAS, CA ;

#define CA_SVM      0
#define CA_GAUSSIAN 1

CA *CAalloc(int width, int height, int depth, MATRIX *m_vox2index, int type, int which_inputs, 
            int wsize, int nscales, char *c1_name, char *c2_name, float *sigmas) ;

float **CAbuildInputs(VOXEL_LIST *vl_total, MRI *mri_intensity, MRI *mri_labels, 
                      int target_label, int which_inputs, int wsize, int nscales, float *sigmas) ;
float **CAbuildTrainingData(VOXEL_LIST *vl_total, int target_label, 
                            float *classes, MRI **mri_smooth, 
                            MRI **mri_grad, MRI **mri_laplacian, MRI *mri_dtrans, 
                            MRI **mri_2nd_deriv_s, int wsize, int nscales, 
                            int which_inputs);
float *CAbuildInputsAtVoxel(VOXEL_LIST *vl, int i, 
                            MRI **mri_smooth, MRI **mri_grad, MRI **mri_laplacian, MRI *mri_dtrans, 
                            MRI *mri_dtrans_grad, MRI **mri_2nd_deriv_s, 
                            int wsize, int nscales, float *ca_inputs, int which_inputs) ;
int ca_ninputs(int which_inputs) ;
int CAtrain(CA *ca, VOXEL_LIST *vl_total, MRI *mri_norm, MRI *mri_aseg, 
            int source_label, int target_label) ;
int CAvoxelToIndex(CA *ca, double x, double y, double z, 
                   double *pxd, double *pyd, double *pzd) ;
int  CAsetSVMparms(CA *ca, double svm_C, double svm_tol, int svm_max_iter) ;
int CAcompleteTraining(CA *ca) ;
int CAwrite(CA *ca, char *fname) ;
CA  *CAread(char *fname) ;
float CAclassifyVoxel(CA *ca, MRI *mri_normals, int x, int y, int z, float *inputs) ;
MRI *CAclassifyBorder(CA *ca, MRI *mri, MRI *mri_aseg, MRI *mri_aseg_edited,
                      int border, int label) ;


#define NINPUTS(which_inputs, wsize, nscales)    ((wsize*wsize*wsize)*nscales*ca_ninputs(which_inputs))

#define CA_INPUT_INTENSITY      0x000001
#define CA_INPUT_LABEL          0x000002
#define CA_INPUT_LAPLACIAN      0x000004
#define CA_INPUT_D2I_R          0x000008
#define CA_INPUT_D2I_A          0x000010
#define CA_INPUT_D2I_S          0x000020
#define CA_INPUT_GRADIENT       0x000040
#define CA_INPUT_DTRANS_GRAD    0x000080
#define CA_INPUT_DTRANS         0x000100


#endif
