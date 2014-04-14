/**
 * @file  gtm.h
 * @brief Routines to create and analyze the Geometric Transfer Matrix (GTM)
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/04/14 22:13:46 $
 *    $Revision: 1.2 $
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


#ifndef GTM_INC
#define GTM_INC

#include "matrix.h"
#include "mri.h"
#include "transform.h"

#ifdef X
#undef X
#endif

typedef struct
{
  char *subject;
  int USF;
  char *apasfile;
  char *ctxannotfile;
  int ctxlhbase,ctxrhbase;
  int KeepHypo;
  int KeepCC;
  int SubSegWM;
  char *wmannotfile;
  int wmlhbase,wmrhbase;
  float dmax;
  int nlist,srclist[300],targlist[300];
  MRI *seg;
  LTA *anat2seg;
} GTMSEG;

int MRIgtmSeg(GTMSEG *gtmseg);
int GTMSEGprint(GTMSEG *gtmseg, FILE *fp);
int GTMdefaultSegReplacmentList(GTMSEG *gtmseg);

#endif
