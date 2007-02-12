/**
 * @file  field_code.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/02/12 18:55:20 $
 *    $Revision: 1.5 $
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


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "error.h"

#include "field_code.h"

char *ReturnFieldName(int which_field)
{
  switch (which_field)
  {
  case OVERLAY_FRAME:
    return(OVERLAY_NAME) ;
    break ;
  case INFLATED_CURV_CORR_FRAME:
    return INFLATED_CURVATURE_NAME;
    break;
  case SULC_CORR_FRAME :
    return SULC_NAME;
    break;
  case CURVATURE_CORR_FRAME:
    return CURVATURE_NAME  ;
    break;
  case  GRAYMID_CORR_FRAME:
    return GRAYMID_NAME;
    break;
  case  T1MID_CORR_FRAME:
    return T1MID_NAME;
    break;
  case  T2MID_CORR_FRAME :
    return T2MID_NAME;
    break;
  case  PDMID_CORR_FRAME :
    return  PDMID_NAME;
    break;
  case  AMYGDALA_CORR_FRAME  :
    return  AMYGDALA_DIST_NAME;
    break;
  case  HIPPOCAMPUS_CORR_FRAME  :
    return HIPPOCAMPUS_DIST_NAME;
    break;
  case  PALLIDUM_CORR_FRAME :
    return PALLIDUM_DIST_NAME;
    break;
  case  PUTAMEN_CORR_FRAME :
    return PUTAMEN_DIST_NAME;
    break;
  case  CAUDATE_CORR_FRAME  :
    return CAUDATE_DIST_NAME;
    break;
  case  LAT_VENTRICLE_CORR_FRAME :
    return LAT_VENTRICLE_DIST_NAME;
    break;
  case  INF_LAT_VENTRICLE_CORR_FRAME :
    return INF_LAT_VENTRICLE_DIST_NAME;
    break;
  default:
    fprintf(stderr,"\nERROR: unknown field\n");
    exit(-1);
    break;
  }
  return NULL;
}

int IsDistanceField(int which_field)
{
  switch (which_field)
  {
  case OVERLAY_FRAME:
    return(0) ;
    break ;
  case INFLATED_CURV_CORR_FRAME:
    return 0;
    break;
  case SULC_CORR_FRAME :
    return 0;
    break;
  case CURVATURE_CORR_FRAME:
    return 0;
    break;
  case  GRAYMID_CORR_FRAME:
    return 0;
    break;
  case  T1MID_CORR_FRAME:
    return 0;
    break;
  case  T2MID_CORR_FRAME :
    return 0;
    break;
  case  PDMID_CORR_FRAME :
    return 0;
    break;
  case  AMYGDALA_CORR_FRAME  :
    return 1;
    break;
  case  HIPPOCAMPUS_CORR_FRAME  :
    return 1;
    break;
  case  PALLIDUM_CORR_FRAME :
    return 1;
    break;
  case  PUTAMEN_CORR_FRAME :
    return 1;
    break;
  case  CAUDATE_CORR_FRAME  :
    return 1;
    break;
  case  LAT_VENTRICLE_CORR_FRAME :
    return 1;
    break;
  case  INF_LAT_VENTRICLE_CORR_FRAME :
    return 1;
    break;
  default:
    fprintf(stderr,"\nERROR: unknown field\n");
    exit(-1);
    break;
  }
  return 0;
}

int InitFieldLabel(FIELD_LABEL *label)
{
  memset(label,0,sizeof(FIELD_LABEL));
  return NO_ERROR;
}

int SetFieldLabel(FIELD_LABEL *label,int field,int frame,float l_corr,float l_pcorr, int navgs, int which_norm)
{
  label->field=field;
  label->frame=frame;
  label->type=IsDistanceField(field);
  label->l_corr=MAX_OF_TWO(0.0,l_corr);
  label->l_pcorr=MAX_OF_TWO(0.0,l_pcorr);
  label->sse=0.0f;
  label->name = NULL ;
  label->navgs = navgs ;
  label->which_norm = which_norm ;
  return NO_ERROR;
}

int
SetFieldName(FIELD_LABEL *label, char *name)
{
  label->name = name ;
  return(NO_ERROR) ;
}
