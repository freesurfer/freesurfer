#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "error.h"

#include "field_code.h"

char *ReturnFieldName(int which_field){
	switch(which_field){
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
		return 	PDMID_NAME;
		break; 
	case  AMYGDALA_CORR_FRAME  :
		return 	AMYGDALA_DIST_NAME;
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

int IsDistanceField(int which_field){
	switch(which_field){
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
