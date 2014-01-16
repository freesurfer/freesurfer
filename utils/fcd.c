/**
 * @file  fcd.c
 * @brief I/O and analysis algorithms for FCDs (focal cortical dysplasias)
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2014/01/16 15:45:38 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include "fcd.h"
#include "error.h"
#include "diag.h"


FCD_DATA   *
FCDloadData(char *sdir, char *subject)
{
  FCD_DATA *fcd ;

  fcd = (FCD_DATA *)calloc(1, sizeof(FCD_DATA)) ;

  return(fcd) ;
}

int
FCDcomputeThicknessLabels(FCD_DATA *fcd, double thickness_thresh, double sigma) 
{
  fcd->nlabels = 1 ;
  fcd->labels[0] = LabelAlloc(10, "no one", "FCD label") ;
  sprintf(fcd->label_names[0], "FCD label") ;
  return(1) ;
}


