///////////////////////////////////////////
// mris_spherical_wavelets.c
// 
// written by Peng Yu
// date: 11/10/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: pengyu $
// Revision Date  : $Date: 2004/11/11 01:05:23 $
// Revision       : $Revision: 1.1 $
////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "icosahedron.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "matrix.h"


//static char vcid[] = "$Id: mris_spherical_wavelets.c,v 1.1 2004/11/11 01:05:23 pengyu Exp $";


int             main(int argc, char *argv[]) ; 
static int      get_option(int argc, char *argv[]) ; 
char            *Progname ;             

static int      ANALYSIS=0;
static int      SYNTHESIS=0;

int 
main(int argc, char *argv[]) 
{ 
	char        fname[STRLEN];
	int         nargs, msec; 
	struct timeb  then ;
	MRIS        *mris_in; *mris_out; 

	Progname = argv[0] ; 
	DiagInit(NULL, NULL, NULL) ; 
	ErrorInit(NULL, NULL, NULL) ; 
	
	for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) 
	{ 
		nargs = get_option(argc, argv) ; 
		argc -= nargs ; 
		argv += nargs ; 
	} 
	
	if (argc < 2) 
		ErrorExit(ERROR_BADPARM, 
							"usage: %s <input surface> <output surface>", Progname); 

	TimerStart(&then) ; 

	mris_in = MRISread(argv[1]) ;
	if (!mris_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, argv[1]) ;    
	mris_out=ReadIcoByOrder(4, 100);
	MRISwrite(mris_out, argv[2]) ; 

	MRISfree(&mris_in) ;
	MRISfree(&mris_out);
	msec = TimerStop(&then) ; 
	fprintf(stdout, "spherical wavelet took %2.1f minutes\n", (float)msec/(1000.0f*60.0f)); 
	exit(0) ; 
	return(0) ; 
} 

/*----------------------------------------------------------------------
	
	Parameters: 
	
	Description: 
	----------------------------------------------------------------------*/


static int 
get_option(int argc, char *argv[]) 
{ 
	int  nargs = 0 ; 
	char *option ; 
	
	option = argv[1] + 1 ;            /* past '-' */ 
	
	if (!stricmp(option, "A"))
	{
		ANALYSIS = 1 ;
		fprintf(stdout,"spherical wavelet analysis\n");
	}
	else if (!stricmp(option, "S"))
	{
		SYNTHESIS = 1 ;
		fprintf(stdout,"Spherical wavelt synthesis\n");
	}
	else	switch (toupper(*option)) 
	{ 
	case '?': 
	case 'U': 
		fprintf(stdout, 
						"usage: %s <input volumes> <output volume>\n", 
						Progname) ; 
		exit(1) ; 
		break ;
	default: 
		fprintf(stdout, "unknown option %s\n", argv[1]) ; 
		exit(1) ; 
		break ; 
	}
	
	return(nargs) ; 
}






