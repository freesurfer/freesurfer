/**
 * @file   dummy.c
 * @author Yasunari Tosa
 * @date   Wed Oct 13 11:47:18 2004
 * 
 * @brief  sample dummy program
 * 
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "fio.h"
#include "const.h"
#include "diag.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "MRIio_old.h"
#include "mri.h"
#include "mrisurf.h"
#include "gca.h"
#include "MC.h"

char *Progname;

int main(int argc, char *argv[])
{
	MRIS *mris_in,*mris_out;

	Progname=argv[0];

	if(argc < 3) {
		fprintf(stderr,"\n\nUSAGE: mris_extract_main_component input_surface output_surface\n\n");
		exit(-1);
	} 

	mris_in=MRISread(argv[1]);

	mris_out=MRISextractMainComponent(mris_in,0);
	
	MRISwrite(mris_out,argv[2]);

	MRISfree(&mris_out);
	MRISfree(&mris_in);
	fprintf(stderr,"\ndone\n\n");
  return 0;
}


