///////////////////////////////////////////
// mris_spharm.c
// 
// written by Peng Yu
// date: 12/09/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: pengyu $
// Revision Date  : $Date: 2005/01/26 22:52:12 $
// Revision       : $Revision: 1.3 $
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


//static char vcid[] = "$Id: mris_spharm.c,v 1.3 2005/01/26 22:52:12 pengyu Exp $";


int             main(int argc, char *argv[]) ; 
static int      get_option(int argc, char *argv[]) ; 
char            *Progname ;             
static int      SphericalHarmonics(int L, int M, float theta, float phi) ;
static double   legendre(int l, int m, float x) ;
static double   factorial(int L, int M) ;
static MATRIX   *ComplexMatrixTranspose(MATRIX *mIn, MATRIX *mOut);
static MATRIX   *ComplexMatrixPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv);
static MRI_SURFACE *center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
static int      LEVEL = 7;
static double   REAL, IMAG;
static float    threshold = 1;
static int      COMPARE = 0;
static char     *cfname;

int 
main(int argc, char *argv[]) 
{ 
	char          fname[STRLEN];
	int           nargs, msec, order, i, j, k, nvertices, dimension, count, fno; 
	float         phi, d, theta, area;
	struct timeb  then ;
	MRIS          *mris_in, *mris_out; 
	MRI_SP        *mrisp ;
	MATRIX        *m_Z, *m_V, *m_c, *m_Z_inv, *m_V_new;
	VERTEX        *v;

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
	for (i = 0; i<mris_out->nvertices; i++)
		mris_out->vertices[i].nsize=1;
	mrisp = MRISPalloc(1, 3); 
#if 1
	MRIScoordsToParameterization(mris_in, mrisp, 1) ;
	//MRISPblur(mrisp, mrisp, 0.5, 0);
	//MRISPblur(mrisp, mrisp, 0.5, 1);
	//MRISPblur(mrisp, mrisp, 0.5, 2);
	MRIScoordsFromParameterization(mrisp, mris_out) ;
#else
	MRISreadOriginalProperties(mris_out, argv[2]) ;
#endif
#if 1 /*just to test if the parameterization is correct */
	MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
	MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
	MRISupdateSurface(mris_out);
	fprintf(stderr, "original area becomes %f\n", mris_out->total_area);
	center_brain(mris_out, mris_out);
	MRISscaleBrain(mris_out, mris_out, sqrt(100000.0f/mris_out->total_area)) ;
	MRISupdateSurface(mris_out);
	for (fno=0; fno<mris_out->nfaces; fno++)
		area += mris_out->faces[fno].area;
	fprintf(stderr, "original area becomes %f\n", area);
	//MRISwrite(mris_out, "/space/xrt/1/users/btquinn/buckner_paper/011121_vc8048/surf/lh.sampled") ; 
	MRISsaveVertexPositions(mris_out, ORIGINAL_VERTICES) ;
	MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
#endif	

	/*Initialize Matrix*/
	dimension = (order+1)*(order+1);
	nvertices = mris_out->nvertices;
	m_c = MatrixAlloc(dimension, 3, MATRIX_COMPLEX) ; 
	m_Z = MatrixAlloc(nvertices, dimension, MATRIX_COMPLEX) ; 
	m_V = MatrixAlloc(nvertices, 3, MATRIX_REAL) ; 
	m_Z_inv = MatrixAlloc(dimension, nvertices, MATRIX_COMPLEX) ; 
	m_V_new = MatrixAlloc(nvertices, 3, MATRIX_COMPLEX) ; 

	for (i=0; i<nvertices; i++)
	{
    v = &mris_out->vertices[i] ;
		phi = atan2(v->y, v->x) ;
    if (phi < 0.0)  phi = 2 * M_PI + phi ; 
    d = 100*100 - v->z*v->z ;
    if (d < 0.0)   d = 0 ;
    theta = atan2(sqrt(d), v->z) ;    
    count = 0;
		*MATRIX_RELT(m_V,i+1,1) = v->origx;
		*MATRIX_RELT(m_V,i+1,2) = v->origy;
		*MATRIX_RELT(m_V,i+1,3) = v->origz;

    for (j=0; j<=order; j++)
			for ( k= -1*j; k<=j; k++)
			{
				count++ ;
				SphericalHarmonics(j, k, theta, phi);
				MATRIX_CELT_REAL(m_Z,i+1,count) = REAL;
				MATRIX_CELT_IMAG(m_Z,i+1,count) = IMAG;
			}
	}

	m_Z_inv=ComplexMatrixPseudoInverse(m_Z, NULL) ;
	MatrixMultiply(m_Z_inv, m_V, m_c) ;
	MatrixPrint(stdout,m_c);

	if (COMPARE)
	{
		MATRIX *m_cc; 
		double r1, r2, m1, m2, diff;
		fprintf(stdout, "%s\n", cfname);
		m_cc = MatrixAlloc(dimension, 3, MATRIX_COMPLEX) ; 
		mris_out = ReadIcoByOrder(LEVEL, 100);
		for (i = 0; i<mris_out->nvertices; i++)
			mris_out->vertices[i].nsize=1;
		mrisp = MRISPalloc(1, 3); 
#if 1
		MRIScoordsToParameterization(mris_in, mrisp, 1) ;
		//MRISPblur(mrisp, mrisp, 0.5, 0);
		//MRISPblur(mrisp, mrisp, 0.5, 1);
		//MRISPblur(mrisp, mrisp, 0.5, 2);
		MRIScoordsFromParameterization(mrisp, mris_out) ;
#else
		MRISreadOriginalProperties(mris_out, cfname) ;
#endif
#if 1 /*just to test if the parameterization is correct */
		MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
		MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
		MRISupdateSurface(mris_out);
		fprintf(stderr, "original area becomes %f\n", mris_out->total_area);
		center_brain(mris_out, mris_out);
		MRISscaleBrain(mris_out, mris_out, sqrt(100000.0f/mris_out->total_area)) ;
		MRISupdateSurface(mris_out);
		for (fno=0; fno<mris_out->nfaces; fno++)
			area += mris_out->faces[fno].area;
		fprintf(stderr, "original area becomes %f\n", area);
		//MRISwrite(mris_out, "/space/xrt/1/users/btquinn/buckner_paper/011121_vc8048/surf/lh.sampled") ; 
		MRISsaveVertexPositions(mris_out, ORIGINAL_VERTICES) ;
		MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
#endif	

		for (i=0; i<nvertices; i++)
		{
			v = &mris_out->vertices[i] ;
			phi = atan2(v->y, v->x) ;
			if (phi < 0.0)  phi = 2 * M_PI + phi ; 
			d = 100*100 - v->z*v->z ;
			if (d < 0.0)   d = 0 ;
			theta = atan2(sqrt(d), v->z) ;    
			count = 0;
			*MATRIX_RELT(m_V,i+1,1) = v->origx;
			*MATRIX_RELT(m_V,i+1,2) = v->origy;
			*MATRIX_RELT(m_V,i+1,3) = v->origz;
			
			for (j=0; j<=order; j++)
				for ( k= -1*j; k<=j; k++)
				{
					count++ ;
					SphericalHarmonics(j, k, theta, phi);
					MATRIX_CELT_REAL(m_Z,i+1,count) = REAL;
					MATRIX_CELT_IMAG(m_Z,i+1,count) = IMAG;
				}
		}
		
		m_Z_inv=ComplexMatrixPseudoInverse(m_Z, NULL) ;
		MatrixMultiply(m_Z_inv, m_V, m_cc) ;
		MatrixPrint(stdout,m_cc);

		/*check for difference*/
		for (count=0, j=0; j<=order; j++)
			for ( k= -1*j; k<=j; k++)
			{
				count++;
				for (i=1;i<=3;i++)
				{
					r1 = MATRIX_CELT_REAL(m_c,count,i);
					r2 = MATRIX_CELT_REAL(m_cc,count,i);
					m1 = MATRIX_CELT_IMAG(m_c,count,i);
					m2 = MATRIX_CELT_IMAG(m_cc,count,i);
					if (sqrt(r1*r1+m1*m1)==0)
						diff = fabs(sqrt(r1*r1+m1*m1)-sqrt(r2*r2+m2*m2));
					else
						diff = fabs(sqrt(r1*r1+m1*m1)-sqrt(r2*r2+m2*m2))/sqrt(r1*r1+m1*m1);
					if (diff>threshold)
					{	MATRIX_CELT_REAL(m_c,count,i) = r2;
					MATRIX_CELT_IMAG(m_c,count,i) = m2;
					fprintf(stdout, "%d %d %f\n",count, i, diff);}
				}
			}
		MatrixFree(&m_cc) ;
	} 

	/*Reconstruct Surface using only part of the coefficients*/
#if 0
	for (i=5; i<=m_c->rows; i++)
		for (j=1; j<=m_c->cols;j++)
		{
			MATRIX_CELT_REAL(m_c,i,j) = 0;
			MATRIX_CELT_IMAG(m_c,i,j) = 0;
		}
#endif
	MatrixMultiply(m_Z, m_c, m_V_new) ;

	for (i=0; i<nvertices; i++)
	{
		v = &mris_out->vertices[i] ;
		v->origx = MATRIX_CELT_REAL(m_V_new,i+1,1) ;
		v->origy = MATRIX_CELT_REAL(m_V_new,i+1,2) ;
		v->origz = MATRIX_CELT_REAL(m_V_new,i+1,3) ;
	}

	MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
	MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
	fprintf(stdout, "Writing reconstructed surface to %s\n", argv[4]);
	MRISwrite(mris_out,argv[4] ) ; 
	MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;


	MatrixFree(&m_V_new) ;
	MatrixFree(&m_Z_inv) ;
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


static int
SphericalHarmonics(int L, int M, float theta, float phi)
{
	double   P, factor;	

	P = legendre(L, abs(M), cos(theta));
	factor = sqrt((L+0.5)*factorial(L,abs(M)));
	REAL = factor * P * cos(abs(M)*phi) /sqrt(2 * M_PI);
	IMAG = factor * P * sin(abs(M)*phi) /sqrt(2 * M_PI);
	if (M < 0)
	{
		REAL = pow(-1,abs(M)) * REAL;
		IMAG = pow(-1,abs(M)) *(-1) * IMAG;
	}
	return(NO_ERROR);
}

static double
legendre(int l, int m, float x)
{
	double  fact, pll, pmm, pmmp1, somx2;
	int    i, ll;

	if ( m<0 || m>l || fabs(x)>1.0 )
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

static double 
factorial(int L, int M)
{
	double ratio=1;
	int i;
	
	for (i=L+M; i>L-M; i--)
		ratio *= (double) 1/i;
	return(ratio);
}

MATRIX *
ComplexMatrixPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv)
{
  MATRIX  *mT, *mTm, *mTm_inv ;

  /* build (mT m)-1 mT */
  mT = ComplexMatrixTranspose(m, NULL) ;
  mTm = MatrixMultiply(mT, m, NULL) ;
  mTm_inv = MatrixInverse(mTm, NULL) ;
  if (!mTm_inv)
  {
    MatrixFree(&mT) ; MatrixFree(&mTm) ;
    return(NULL) ;
  }
  m_pseudo_inv = MatrixMultiply(mTm_inv, mT, m_pseudo_inv) ;

  MatrixFree(&mT) ; MatrixFree(&mTm) ; MatrixFree(&mTm_inv) ; 
  return(m_pseudo_inv) ;
}

MATRIX  *
ComplexMatrixTranspose(MATRIX *mIn, MATRIX *mOut)
{
  int  row, col, rows, cols ;

  if (!mOut)
  {
    mOut = MatrixAlloc(mIn->cols, mIn->rows, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
		{
			MATRIX_CELT_REAL(mOut,col,row) = MATRIX_CELT_REAL(mIn,row,col)  ;
			MATRIX_CELT_IMAG(mOut,col,row) = -1*MATRIX_CELT_IMAG(mIn,row,col)  ;
		}
  }

  return(mOut) ;
}


MRI_SURFACE *
center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int         fno, vno ;
  FACE        *face ;
	VERTEX      *vdst;
  float       x, y, z, x0, y0, z0 ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x0 = y0 = z0 = 0 ;   /* silly compiler warning */

  for (fno = 0 ; fno < mris_src->nfaces ; fno++)
  {
    face = &mris_dst->faces[fno] ;
    if (face->ripflag)
      continue ;
    x = mris_dst->vertices[face->v[0]].x;
    y = mris_dst->vertices[face->v[0]].y;
    z = mris_dst->vertices[face->v[0]].z;
    x += mris_dst->vertices[face->v[1]].x;
    y += mris_dst->vertices[face->v[1]].y;
    z += mris_dst->vertices[face->v[1]].z;
    x += mris_dst->vertices[face->v[2]].x;
    y += mris_dst->vertices[face->v[2]].y;
    z += mris_dst->vertices[face->v[2]].z;
		x /= face->area; y/= face->area; z/= face->area;
		x0 += x; y0 += y; z0 += z;
  }
  x0 /= mris_dst->total_area ;  y0 /= mris_dst->total_area  ; z0 /= mris_dst->total_area ;
 
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag)
      continue ;
    vdst->x -= x0 ; vdst->y -= y0 ; vdst->z -= z0 ;
  }

  mris_dst->xctr = mris_dst->yctr = mris_dst->zctr = 0 ;
  return(mris_dst) ;
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
	else if (!stricmp(option, "C"))
	{
		cfname = argv[2] ;
		COMPARE = 1;
		fprintf(stdout,"Compare with %s \n", cfname);
		nargs=1;
	}
	else if (!stricmp(option, "T"))
	{
		threshold = atof(argv[2]) ;
		fprintf(stdout,"Set threshold as %f \n", threshold);
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






