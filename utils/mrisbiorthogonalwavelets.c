/**
 * @file  mrisbiorthogonalwavelets.c
 * @brief Routines for biothogonal wavelets transformation.
 *
 * References:
 * P. Schroder and W. Sweldens. Spherical wavelets: Texture processing. 
 * In Rendering Techniques '95. Springer Verlag, 1995.
 * P. Schroder and W. Sweldens. Spherical wavelets: Efficiently representing 
 * functions on the sphere. Computer Graphics Proceedings (SIGGRAPH 95), 
 * pages 161-172, 1995.
 */
/*
 * Original Author: Peng Yu
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:46 $
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

MRI_SURFACE  *wavelet_analysis_curv(MRI_SURFACE *mris_out, int order) ;
MRI_SURFACE  *wavelet_analysis_vec(MRI_SURFACE *mris_out, int order);
MRI_SURFACE  *wavelet_synthesis_curv(MRI_SURFACE *mris_out, int order) ;
MRI_SURFACE  *wavelet_synthesis_vec(MRI_SURFACE *mris_out, int order);


MRI_SURFACE *
wavelet_analysis_curv(MRI_SURFACE *mris_out, int order)
{
	int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0; 
	MRIS          *mris_high; 
	VERTEX        *vm_out, *vm_high, *v;
	double         s_jkm;

	/* Initialize Ij,k*/
	for (vno = 0 ; vno<mris_out->nvertices; vno++)
		{
			vm_out = &mris_out->vertices[vno];
			vm_out->val = 1;
		}		
	
	/*Iteratively compute Ij,k*/
	for (i=order;i>0;i--)
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;
			number = IcoNVtxsFromOrder(i-1); //the start of m vertices
			for (m = number; m<mris_high->nvertices; m++) 
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{					
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						v->val += 0.5*vm_out->val ;
					}
				for(; nnum<vm_high->v2num; nnum++)
					if ( vm_high->v[nnum]<number ) //B(j,m)
					{						
						k = vm_high->v[nnum];
						if (flag==0) b1=k;
						else b2=k;
						flag++;
						v = &mris_out->vertices[k];	
						v->val += 0.125*vm_out->val ;
					}
				for(; nnum<vm_high->v3num; nnum++)
					if ( vm_high->v[nnum]<number ) //C(j,m)
					{					
						k = vm_high->v[nnum];
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];	
							v->val -= 0.0625*vm_out->val ;
						}
					}
			}
		}
	
		/*Analysis Stage I:*/
		for (i=order;i>0;i--) 
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;

			number = IcoNVtxsFromOrder(i-1); //the start of m vertices
			/* compute Yj,m for each m vertices */
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m]; 
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
					{
						k = vm_high->v[nnum] ; 
						v = &mris_out->vertices[k];
						vm_out->curv -= 0.5*v->curv;
						//if(m==67770) fprintf(stdout, "%f, %d, %f\n", v->curv, k, vm_out->curv);
					}
				for(; nnum<vm_high->v2num; nnum++) //second order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
					{
						k = vm_high->v[nnum] ; 
						if (flag==0) b1=k;
						else b2=k;
						flag++;
						v = &mris_out->vertices[k];
						vm_out->curv -= 0.125*v->curv;		
						//if(m==67770) fprintf(stdout, "%f, %d, %f\n", v->curv, k, vm_out->curv);
					}
				for(; nnum<vm_high->v3num; nnum++)
					if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
					{
						k = vm_high->v[nnum] ; 
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];
							vm_out->curv += 0.0625*v->curv;	
							//if(m==67770) fprintf(stdout, "%f, %d, %f\n", v->curv, k, vm_out->curv);
						}
					}
			}


			/*Analysis Stage II: */
			/*Compute Lamda(j,k) using the Yita(j,m)*/
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						s_jkm = vm_out->val/2/v->val; 
						//if(k==6642) fprintf(stdout, "%f, %f, %f, %f, %d\n", vm_out->val, v->val, s_jkm, vm_out->curv, m);
						v->curv += s_jkm*vm_out->curv;
					}

			}
		} 
	MRISfree(&mris_high) ;		
	return(mris_out);
}

MRI_SURFACE *
wavelet_analysis_vec(MRI_SURFACE *mris_out, int order)
{
	int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0; 
	MRIS          *mris_high; 
	VERTEX        *vm_out, *vm_high, *v;
	double         s_jkm;

		/* Initialize Ij,k*/
		for (vno = 0 ; vno<mris_out->nvertices; vno++)
		{
			vm_out = &mris_out->vertices[vno];
			vm_out->val = 1;
		}		

		/*Iteratively compute Ij,k*/
		for (i=order;i>0;i--)
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;
			number = IcoNVtxsFromOrder(i-1); //the start of m vertices
			for (m = number; m<mris_high->nvertices; m++) 
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{					
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						v->val += 0.5*vm_out->val ;
					}
				for(; nnum<vm_high->v2num; nnum++)
					if ( vm_high->v[nnum]<number ) //B(j,m)
					{						
						k = vm_high->v[nnum];
						if (flag==0) b1=k;
						else b2=k;
						flag++;
						v = &mris_out->vertices[k];	
						v->val += 0.125*vm_out->val ;
					}
				for(; nnum<vm_high->v3num; nnum++)
					if ( vm_high->v[nnum]<number ) //C(j,m)
					{					
						k = vm_high->v[nnum];
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];	
							v->val -= 0.0625*vm_out->val ;
						}
					}
			}
		}


		/*Analysis Stage I:*/
		for (i=order;i>0;i--) 
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;

			number = IcoNVtxsFromOrder(i-1); //the start of m vertices
			/* compute Yj,m for each m vertices */
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m]; 
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
					{
						k = vm_high->v[nnum] ; 
						v = &mris_out->vertices[k];
						vm_out->origx -= 0.5*v->origx;
						vm_out->origy -= 0.5*v->origy;
						vm_out->origz -= 0.5*v->origz;	
					}
				for(; nnum<vm_high->v2num; nnum++) //second order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
					{
						k = vm_high->v[nnum] ; 
						if (flag==0) b1=k;
						else b2=k;
						flag++;
						v = &mris_out->vertices[k];
						vm_out->origx -= 0.125*v->origx;
						vm_out->origy -= 0.125*v->origy;
						vm_out->origz -= 0.125*v->origz;			
					}
				for(; nnum<vm_high->v3num; nnum++)
					if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
					{
						k = vm_high->v[nnum] ; 
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];
							vm_out->origx += 0.0625*v->origx;
							vm_out->origy += 0.0625*v->origy;
							vm_out->origz += 0.0625*v->origz;		
						}
					}
			}


			/*Analysis Stage II: */
			/*Compute Lamda(j,k) using the Yita(j,m)*/
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						s_jkm = vm_out->val/2/v->val; 
						v->origx += s_jkm*vm_out->origx;
						v->origy += s_jkm*vm_out->origy;
						v->origz += s_jkm*vm_out->origz;
					}

			}

#if 0
		
			{	//optional output
			  char fname[200];
			  MRISfree(&mris_high);
			  mris_high = ReadIcoByOrder(i-1, 100); 
			  for (m = 0; m<mris_high->nvertices; m++)
				{
				  mris_high->vertices[m].x = mris_out->vertices[m].origx;
				  mris_high->vertices[m].y = mris_out->vertices[m].origy;
				  mris_high->vertices[m].z = mris_out->vertices[m].origz;
				}
			  sprintf(fname,"/space/birn/42/users/pengyu/simulation/010223_61223/surf/lh.wavelet%d",i-1);
			  MRISwrite(mris_high, fname) ;	   
			}
#endif
		} 
	MRISfree(&mris_high) ;
	return(mris_out);
}

MRI_SURFACE *
wavelet_synthesis_curv(MRI_SURFACE *mris_out, int order)
{
	int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0; 
	MRIS          *mris_high; 
	VERTEX        *vm_out, *vm_high, *v;
	double         s_jkm;

		/*Initialize Ij,k*/
		for (vno = 0; vno<mris_out->nvertices; vno++)
		{
			vm_out = &mris_out->vertices[vno];
			vm_out->val = 1;
		}			

		/*Iteratively compute Ij,k*/
		for (i=order;i>0;i--)
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;
			number = IcoNVtxsFromOrder(i-1); //the start of m vertices
			for (m = number; m<mris_high->nvertices; m++) 
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{					
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						v->val += 0.5*vm_out->val ;
					}
				for(; nnum<vm_high->v2num; nnum++)
					if ( vm_high->v[nnum]<number ) //B(j,m)
					{						
						k = vm_high->v[nnum];
						if (flag==0) b1=k;
						else b2=k;
						flag++;
						v = &mris_out->vertices[k];	
						v->val += 0.125*vm_out->val ;
					}
				for(; nnum<vm_high->v3num; nnum++)
					if ( vm_high->v[nnum]<number ) //C(j,m)
					{					
						k = vm_high->v[nnum];
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];	
							v->val -= 0.0625*vm_out->val ;
						}
					}
			}
		}


		for (i=1;i<=order;i++)
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;
			number = IcoNVtxsFromOrder(i-1); //the start of m vertices

			/* Synthesis Stage I */
			/* Compute Lamda(j+1,k) using the Yita(j,m) */
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						s_jkm = vm_out->val/2/v->val; 
						v->curv -= s_jkm*vm_out->curv;
					}
			}

			/* compute Lamda(j+1,m) for each m vertices */
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m]; 
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
					{
						k = vm_high->v[nnum] ; 
						v = &mris_out->vertices[k];	
						vm_out->curv += 0.5*v->curv;
					}
				for(; nnum<vm_high->v2num; nnum++) //second order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
					{
						k = vm_high->v[nnum] ;
						if (flag==0) b1=k;
						else b2=k;
						flag++;
 						v = &mris_out->vertices[k];	
						vm_out->curv += 0.125*v->curv;
					}
				for(; nnum<vm_high->v3num; nnum++) //third order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
					{
						k = vm_high->v[nnum] ; 
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];	
							vm_out->curv -= 0.0625*v->curv;
						}
					}
			}
		} 
	MRISfree(&mris_high) ;
	return(mris_out);
}

MRI_SURFACE *
wavelet_synthesis_vec(MRI_SURFACE *mris_out, int order)
{
	int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0; 
	MRIS          *mris_high; 
	VERTEX        *vm_out, *vm_high, *v;
	double         s_jkm;

				
		/*Initialize Ij,k*/
		for (vno = 0; vno<mris_out->nvertices; vno++)
		{
			vm_out = &mris_out->vertices[vno];
			vm_out->val = 1;
		}			

		/*Iteratively compute Ij,k*/
		for (i=order;i>0;i--)
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;
			number = IcoNVtxsFromOrder(i-1); //the start of m vertices
			for (m = number; m<mris_high->nvertices; m++) 
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{					
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						v->val += 0.5*vm_out->val ;
					}
				for(; nnum<vm_high->v2num; nnum++)
					if ( vm_high->v[nnum]<number ) //B(j,m)
					{						
						k = vm_high->v[nnum];
						if (flag==0) b1=k;
						else b2=k;
						flag++;
						v = &mris_out->vertices[k];	
						v->val += 0.125*vm_out->val ;
					}
				for(; nnum<vm_high->v3num; nnum++)
					if ( vm_high->v[nnum]<number ) //C(j,m)
					{					
						k = vm_high->v[nnum];
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];	
							v->val -= 0.0625*vm_out->val ;
						}
					}
			}
		}


		for (i=1;i<=order;i++)
		{
			mris_high = ReadIcoByOrder(i, 100); //higher order surface
			for (m = 0; m<mris_high->nvertices; m++)
				mris_high->vertices[m].nsize=1;
			MRISsetNeighborhoodSize(mris_high, 3) ;
			number = IcoNVtxsFromOrder(i-1); //the start of m vertices

			/* Synthesis Stage I */
			/* Compute Lamda(j+1,k) using the Yita(j,m) */
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m];
				vm_high = &mris_high->vertices[m];
				for (nnum=0; nnum<vm_high->vnum; nnum++)
					if ( vm_high->v[nnum]<number ) //A(j,m)
					{
						k = vm_high->v[nnum];
						v = &mris_out->vertices[k];	
						s_jkm = vm_out->val/2/v->val; 
						v->origx -= s_jkm*vm_out->origx;
						v->origy -= s_jkm*vm_out->origy;
						v->origz -= s_jkm*vm_out->origz;
					}
			}

			/* compute Lamda(j+1,m) for each m vertices */
			for (m = number; m<mris_high->nvertices; m++)
			{
				vm_out = &mris_out->vertices[m]; 
				vm_high = &mris_high->vertices[m];
				flag=0;
				for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
					{
						k = vm_high->v[nnum] ; 
						v = &mris_out->vertices[k];	
						vm_out->origx += 0.5*v->origx;
						vm_out->origy += 0.5*v->origy;
						vm_out->origz += 0.5*v->origz;	
					}
				for(; nnum<vm_high->v2num; nnum++) //second order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
					{
						k = vm_high->v[nnum] ;
						if (flag==0) b1=k;
						else b2=k;
						flag++;
 						v = &mris_out->vertices[k];	
						vm_out->origx += 0.125*v->origx;
						vm_out->origy += 0.125*v->origy;
						vm_out->origz += 0.125*v->origz;
					}
				for(; nnum<vm_high->v3num; nnum++) //third order neighborhood
					if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
					{
						k = vm_high->v[nnum] ; 
						flag=0; //C has to be a second-order neighbor of B
						for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
							if(mris_high->vertices[b1].v[cno]==k) flag=1;
						for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
							if(mris_high->vertices[b2].v[cno]==k) flag=1;
						if(flag)
						{
							v = &mris_out->vertices[k];	
							vm_out->origx -= 0.0625*v->origx;
							vm_out->origy -= 0.0625*v->origy;
							vm_out->origz -= 0.0625*v->origz;		
						}
					}
			}
		} 
	MRISfree(&mris_high) ;
	return(mris_out);
}
