///////////////////////////////////////////
// mris_spharm.c
// 
// written by Peng Yu
// date: 12/09/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: pengyu $
// Revision Date  : $Date: 2004/12/10 01:01:20 $
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


//static char vcid[] = "$Id: mris_spharm.c,v 1.1 2004/12/10 01:01:20 pengyu Exp $";


int             main(int argc, char *argv[]) ; 
static int      get_option(int argc, char *argv[]) ; 
char            *Progname ;             
static CPRT *   SphericalHarmonics(int L, int M, int theta, int phi);
static float    legendre(int l, int, m, float x) ;

#define         LEVEL = 5;

int 
main(int argc, char *argv[]) 
{ 
	char          fname[STRLEN];
	int           nargs, msec, order, i, j, k, nvertices, dimension, count; 
	float         phi, d, theta;
	struct timeb  then ;
	MRIS          *mris_in, *mris_out; 
	MRI_SP        *mrisp ;
	MATRIX        *m_Z, *m_V, *m_c;
	VERTEX        *v;
	CPTR          *z;

	Progname = argv[0] ; 
	DiagInit(NULL, NULL, NULL) ; 
	ErrorInit(NULL, NULL, NULL) ; 
	
	for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) 
	{ 
		nargs = get_option(argc, argv) ; 
		argc -= nargs ; 
		argv += nargs ; 
	} 
	
	if (argc < 4) 
		ErrorExit(ERROR_BADPARM, 
							"usage: %s <input surface> <orig surface> <finest order> <output surface>", Progname); 

	TimerStart(&then) ; 

	mris_in = MRISread(argv[1]) ;
	if (!mris_in)
		ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
							Progname, argv[1]) ; 
	fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);
	MRISreadOriginalProperties(mris_in, argv[2]) ;
	fprintf(stdout, "Reading original surface from %s\n", argv[2]);

	order = atoi (argv[3]);
	fprintf(stdout, "SPHARM decomposition up tp level %s\n", argv[3]);

	/*Sample the original surface*/

	mris_out = ReadIcoByOrder(LEVEL, 100);
	for (m = 0; m<mris_out->nvertices; m++)
		mris_out->vertices[m].nsize=1;
	mrisp = MRISPalloc(1, 3); 
	MRIScoordsToParameterization(mris_in, mrisp, 1) ;
	MRISPblur(mrisp, mrisp, 2, 0);
	MRISPblur(mrisp, mrisp, 2, 1);
	MRISPblur(mrisp, mrisp, 2, 2);
	MRIScoordsFromParameterization(mrisp, mris_out) ;
	
#if 1 /*just to test if the parameterization is correct */
	MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
	MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
	MRISwrite(mris_out, "/autofs/space/dijon_004/ksong/BIRN_processed_data/prospective/buckner/CORTICAL_THINNING/RECONS/001009_vc5398/surf/lh.hippocampus.recovered") ; 
	MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
#endif

	/*Initialize Matrix*/
	dimension = (order+1)*(order+1);
	nvertices = mris_out->nvertices;
  v_v = VectorAlloc(3, MATRIX_REAL) ;
	m_c = MatrixAlloc(dimension, 3, MATRIX_COMPLEX) ; 
	m_Z = MatrixAlloc(nvertices, dimension, MATRIX_COMPLEX) ; 
	m_V = MatrixAlloc(nvertices, 3, MATRIX_REAL) ; 

	for (i=0; i<nvertices, i++)
	{
    v = &mris->vertices[i] ;
		phi = atan2(v->y, v->x) ;
    if (phi < 0.0)  phi = 2 * M_PI + phi ; 
    d = 100*100 - v->z*v->z ;
    if (d < 0.0)   d = 0 ;
    theta = atan2(sqrt(d), v->z) ;    
    count = 0;
		*MATRIX_RELT(m_V,i+1,1) = v->orig_x;
		*MATRIX_RELT(m_V,i+1,2) = v->orig_y;
		*MATRIX_RELT(m_V,i+1,3) = v->orig_z;

    for (j=0; j<=order; j++)
			for ( k= -1*j; k<=j; k++)
			{
				count = count+1 ;
				z = SphericalHarmonics(j, k, theta, phi);
				*MATRIX_CELT_REAL(m_Z,i,count) = z->real;
				*MATRIX_CELT_IMAG(m_Z,i,count) = z->imag;
			}
	}

	MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
	MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
	fprintf(stdout, "Writing wavelets coefficient of original surface to %s\n", argv[4]);
	MRISwrite(mris_out,argv[4] ) ; 
	MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
	
	MatrixFree(&m_Z) ;
	MatrixFree(&m_V) ;
	MatrixFree(&m_c) ;
	MRISPfree(&mrisp) ; 
	MRISfree(&mris_in) ;
	MRISfree(&mris_out);
	msec = TimerStop(&then) ; 
	fprintf(stdout, "spherical wavelet took %2.1f minutes\n", (float)msec/(1000.0f*60.0f)); 
	exit(0) ; 
	return(0) ; 
} 


static CPRT *
SphericalHarmonics(int L, int M, int theta, int phi)
{
	CPTR    *value;
	float   P, real, imag, factor;	

	P = legendre(L, abs(M), cos(theta));
	factor = sqrt((L+0.5)*factorial(L-abs(M))/factorial(L+abs(M)));
	real = factor * P * cos(abs(M)*phi) /sqrt(2 * M_PI);
	imag = factor * P * sin(abs(M)*phi) /sqrt(2 * M_PI);
	if (M < 0)
	{
		real = power(-1,abs(M)) * real;
		imag = power(-1,abs(M)) *(-1) * imag;
	}
	value->real = real;
	value->imag = imag;
	return(value);
}

static float
legendre(int l, int, m, float x)
{
	float  fact, pll, pmm, pmmp1, somx2;
	int    i, ll;

	if ( m<0 || m>1 || fabs(x)>1.0 )
		fprintf(stdout, "Bad arguments in routine legendre");
	pmm = 1.0;
	if ( m > 0)
	{
		somx2 = sqrt((1.0-x)*(1.0+x));
		fact = 1.0;
		for ( i=1 ; i<=m ; i++)
		{
			pmm *= -fact * somx2;
			fact += 2.0;
		}
	}

	if ( l==m )
		return pmm;
	else
	{
		pmmp1 = x*(2*m+1)*pmm;
		if ( l==(m+1) )
			return pmmp1;
		else
		{
			for ( ll=m+2 ; ll<=l ; ll++ )
			{
				pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm = pmmp1;
				pmmp1 = pll;
			}
			return pll;
		}
	}
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
	
	if (!stricmp(option, "L"))
	{
		LEVEL = atoi(argv[2]) ;
		fprintf(stdout,"Use %d order icosahedron \n", LEVEL);
		nargs=1;
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






