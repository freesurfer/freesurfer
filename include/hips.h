/**
 * @file  hips.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.9 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


/*
 *       FILE NAME:   hips.h
 *
 *       DESCRIPTION: hips2 prototypes
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
 *
*/

#ifndef HIPS_H
#define HIPS_H

/*-----------------------------------------------------
                    TYPEDEFS & STRUCTURES
-------------------------------------------------------*/
#define MaxDisks 8
typedef struct
{
  /*
    indices (i,j) refer to which logmap point the current image point
    belongs (MAX number of logmap points is MaxDisks)
  */
  int    xi_i[MaxDisks];
  int    eta_j[MaxDisks];
  float  weight[MaxDisks]; /* differs with the type of nonuniform filetring */
}
LMtable;




/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <hipl_format.h>
#include "hipsh.h"
#include "hipsu.h"

#if 0
int h_mulscale(struct header *hdi, struct header *hdo, Pixelval *b) ;
int h_divscale(struct header *hdi, struct header *hdo, Pixelval *b) ;
int h_add(struct header *hdi1, struct header *hdi2, struct header *hdo) ;
int h_linscale(struct header *hdi, struct header * hdo, float b, float c) ;
int h_div(struct header *hdi1, struct header * hdi2, struct header *hdo) ;
int h_tof(struct header *hdi, struct header *hdo) ;
int h_toc(struct header *hdi, struct header *hdo) ;
int h_reduce(struct header *hdi, struct header *hdo, int xf, int yf) ;
int h_flipquad(struct header *hdi, struct header *hdo) ;
int h_rot180(struct header *hdi, struct header *hdo) ;
int h_minmax(struct header *hd, Pixelval *minval, Pixelval *maxval,int nzflag);
int h_discedge(struct header *hdi, struct header *hdo,int size, float varcrit);
int h_canny(struct header *Isrc, struct header *Idst, double sigma,
            int mask_size,double lfrac,double hfrac,int dothin);
int h_morphdil(struct header *hdi, struct header *hde, struct header *hdo,
               int centerr, int centerc,int gray) ;
int h_logmap(struct header *phd,struct header *phd_LM,float a,
             int R,int rows,int cols,int Fx,int Fy,float disk_size) ;
int h_logmap_f(struct header *phd,struct header *phd_LM,float a,
               int R,int rows,int cols,int Fx,int Fy,float disk_size) ;
int h_LMtable_uniform(LMtable *ptable,struct header *phd_w,int Irows,int Icols,
                      float a,int R,int Fx,int Fy,float disk_size) ;
#endif

#endif
