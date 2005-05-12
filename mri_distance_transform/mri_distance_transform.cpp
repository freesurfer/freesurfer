#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

extern "C" {
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "timer.h"
}
#include "fastmarching.h"

char *Progname ;

int main(int argc, char *argv[])
{
	MRI *mri,*mri_distance;
	int label,mode;
	float max_distance;

	max_distance=10;
	mode=1;

	Progname=argv[0];

	fprintf(stderr,"mri_distance_transform input_distance label max_distance mode[=1] output\n");
	fprintf(stderr,"mode : 1 = outside , mode : 2 = inside , mode : 3 = both\n");

	mri=MRIread(argv[1]);
	label=atoi(argv[2]);
	max_distance=atof(argv[3]);
	mode=atoi(argv[4]);

	fprintf(stderr,"label=%d distance=%f mode=%d\n",label,max_distance,mode);

	mri_distance=MRIalloc(mri->width,mri->height,mri->depth,MRI_FLOAT);

	mri_distance=MRIextractDistanceMap(mri,mri_distance,label, max_distance, mode);

	MRIwrite(mri_distance,argv[5]);
	
	MRIfree(&mri);
	MRIfree(&mri_distance);

	return 0;
}
