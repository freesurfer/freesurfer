/**
 * @brief utility functions contributed by Graham Wideman
 *
 */
/*
 * Original Author: Graham Wideman
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


#ifndef GW_UTILS_H
#define GW UTILS_H

#include "mrisurf.h"
// ..which includes mrishash.h
#include "mri.h"


//=============================================================================
// 
typedef struct {
  float x;
  float y;
  float z;
} GWUTILS_VERTEX;
typedef struct
{
  int vno[3] ; // zero based
}
GWUTILS_FACE ;

MRI_SURFACE *GWU_make_surface_from_lists(GWUTILS_VERTEX * vertices,
                                         int vertexcount,
                                         GWUTILS_FACE   * faces ,
                                         int facecount ) ;
//---------------------------------------------
// Write MHT data to an MRI volume for visualization
// HMT volume will be 16-bit "short" voxel data
//---------------------------------------------
typedef enum {
    MFMM_None       = 1,  // Set mht'ed voxels to 1
    MFMM_Num           ,  // Set mht'ed voxels to first face or
                          // vertex "no" (number) encountered
                          // (one-based)
    MFMM_NumDiv16      ,  // Same as MFMM_Num, but (no >> 4) + 1;
    MFMM_Count            // Count of faces/vertices listed in this mht voxel
} MFMM_Option_t;

MRI *MRIFromMHTandMRIS(MHT * mht, MRIS * mris, MFMM_Option_t mfmm_option);

// Some simple log functions for test programs.
extern int gw_log_init(char * AProgname, 
                       char * AProgversion, 
                       char * ALogfilepath, 
                       int newfile);
extern void gw_log_message(const char * msg);
extern void gw_log_timestamp(const char * label);
extern void gw_log_begin(void);
extern void gw_log_end(void);
//=============================================================================

#endif
