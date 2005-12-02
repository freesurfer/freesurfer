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
#include "tags.h"
#include "version.h"
#include "MC.h"

#define USE_WM -1
#define SQR(x) ((x)*(x))

char *Progname;


static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;
static int corners = 1;

static int mriRemoveEdgeConfiguration(MRI *mri_seg, MRI *mri_orig, int label){
	static int niter=0;
	int i,j,k;
	int ntotal=0,nmodified,nfound,npass;
	
	niter++;
	fprintf(stderr,"\nIteration Number : %d",niter);

	/* dealing with xy-plane */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 0 ; k < mri_seg->depth ; k++)
			for(j = 0 ; j < mri_seg->height-1 ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)!=label) continue;
					if(MRIvox(mri_seg,i+1,j+1,k)!=label) continue;
					if((MRIvox(mri_seg,i,j+1,k)==label) || (MRIvox(mri_seg,i+1,j,k)==label)) continue;
					/* select the brigther voxel */
					if(MRIgetVoxVal(mri_orig,i,j+1,k,0)>MRIgetVoxVal(mri_orig,i+1,j,k,0))
						MRIvox(mri_seg,i,j+1,k)=label;
					else
						MRIvox(mri_seg,i+1,j,k)=label;
					nfound++;
				} 
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (xy+): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}
	/* dealing with xy-plane */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 0 ; k < mri_seg->depth ; k++)
			for(j = 0 ; j < mri_seg->height-1 ; j++)
				for(i = 1 ; i < mri_seg->width ; i++){
					if(MRIvox(mri_seg,i,j,k)!=label) continue;
					if(MRIvox(mri_seg,i-1,j+1,k)!=label) continue;
					if((MRIvox(mri_seg,i,j+1,k)==label) || (MRIvox(mri_seg,i-1,j,k)==label)) continue;
					/* select the brigther voxel */
					if(MRIgetVoxVal(mri_orig,i,j+1,k,0)>MRIgetVoxVal(mri_orig,i-1,j,k,0))
						MRIvox(mri_seg,i,j+1,k)=label;
					else
						MRIvox(mri_seg,i-1,j,k)=label;
					nfound++;
				} 
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (xy-): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}

  /* dealing with yz-plane */
  nfound=1;
  nmodified=0;
	npass=0;
  while(nfound){
    nfound=0;
		npass++;
    for(k = 0 ; k < mri_seg->depth-1 ; k++)
      for(j = 0 ; j < mri_seg->height-1 ; j++)
        for(i = 0 ; i < mri_seg->width ; i++){
          if(MRIvox(mri_seg,i,j,k)!=label) continue;
          if(MRIvox(mri_seg,i,j+1,k+1)!=label) continue;
          if((MRIvox(mri_seg,i,j+1,k)==label) || (MRIvox(mri_seg,i,j,k+1)==label)) continue;
          /* select the brigther voxel */
          if(MRIgetVoxVal(mri_orig,i,j+1,k,0)>MRIgetVoxVal(mri_orig,i,j,k+1,0))
            MRIvox(mri_seg,i,j+1,k)=label;
          else
            MRIvox(mri_seg,i,j,k+1)=label;
          nfound++;
        }
    nmodified+=nfound;
    ntotal += nfound;
    fprintf(stderr,"\npass %3d (yz+): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
  }
	/* dealing with yz-plane */
  nfound=1;
  nmodified=0;
	npass=0;
  while(nfound){
    nfound=0;
		npass++;
    for(k = 1 ; k < mri_seg->depth ; k++)
      for(j = 0 ; j < mri_seg->height-1 ; j++)
        for(i = 0 ; i < mri_seg->width ; i++){
					if(MRIvox(mri_seg,i,j,k)!=label) continue;
          if(MRIvox(mri_seg,i,j+1,k-1)!=label) continue;
          if((MRIvox(mri_seg,i,j+1,k)==label) || (MRIvox(mri_seg,i,j,k-1)==label)) continue;
          /* select the brigther voxel */
          if(MRIgetVoxVal(mri_orig,i,j+1,k,0)>MRIgetVoxVal(mri_orig,i,j,k-1,0))
            MRIvox(mri_seg,i,j+1,k)=label;
          else
            MRIvox(mri_seg,i,j,k-1)=label;
          nfound++;
        }
    nmodified+=nfound;
    ntotal += nfound;
    fprintf(stderr,"\npass %3d (yz-): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
  }

	/* dealing with xz-plane */
  nfound=1;
  nmodified=0;
	npass=0;
  while(nfound){
    nfound=0;
    npass++;
    for(k = 0 ; k < mri_seg->depth-1 ; k++)
      for(j = 0 ; j < mri_seg->height ; j++)
        for(i = 0 ; i < mri_seg->width-1 ; i++){
          if(MRIvox(mri_seg,i,j,k)!=label) continue;
          if(MRIvox(mri_seg,i+1,j,k+1)!=label) continue;
          if((MRIvox(mri_seg,i+1,j,k)==label) || (MRIvox(mri_seg,i,j,k+1)==label)) continue;
          /* select the brigther voxel */
          if(MRIgetVoxVal(mri_orig,i+1,j,k,0)>MRIgetVoxVal(mri_orig,i,j,k+1,0))
            MRIvox(mri_seg,i+1,j,k)=label;
          else
            MRIvox(mri_seg,i,j,k+1)=label;
          nfound++;
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,"\npass %3d (xz+): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
  }
	/* dealing with xz-plane */
  nfound=1;
  nmodified=0;
	npass=0;
  while(nfound){
    nfound=0;
    npass++;
    for(k = 0 ; k < mri_seg->depth-1 ; k++)
      for(j = 0 ; j < mri_seg->height ; j++)
        for(i = 1 ; i < mri_seg->width ; i++){
          if(MRIvox(mri_seg,i,j,k)!=label) continue;
          if(MRIvox(mri_seg,i-1,j,k+1)!=label) continue;
          if((MRIvox(mri_seg,i-1,j,k)==label) || (MRIvox(mri_seg,i,j,k+1)==label)) continue;
          /* select the brigther voxel */
          if(MRIgetVoxVal(mri_orig,i-1,j,k,0)>MRIgetVoxVal(mri_orig,i,j,k+1,0))
            MRIvox(mri_seg,i-1,j,k)=label;
          else
            MRIvox(mri_seg,i,j,k+1)=label;
          nfound++;
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,"\npass %3d (xz-): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
  }

	return ntotal;
}

static int mriRemoveCornerConfiguration(MRI *mri_seg, MRI *mri_orig, int label){
	static int niter=0;
	int i,j,k,p,refp,ind1_i[6],ind1_j[6],ind1_k[6],ind2_i[6],ind2_j[6],ind2_k[6];
	int ntotal=0,nmodified,nfound,npass;
	float dist[6],maxdist;
	
	
	niter++;
	fprintf(stderr,"\nIteration Number : %d",niter);

	/* dealing with i+1,j+1,k+1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 0 ; k < mri_seg->depth-1; k++)
			for(j = 0 ; j < mri_seg->height-1 ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)!=label) continue;
					if(MRIvox(mri_seg,i+1,j+1,k+1)!=label) continue;

					/* find problematic configuration */
					if(((MRIvox(mri_seg,i+1,j,k)==label) 
							&& ((MRIvox(mri_seg,i+1,j+1,k)==label)
									||(MRIvox(mri_seg,i+1,j,k+1)==label))) || 
						 ((MRIvox(mri_seg,i,j+1,k)==label) 
							&& ((MRIvox(mri_seg,i,j+1,k+1)==label)
									||(MRIvox(mri_seg,i+1,j+1,k)==label))) || 
						 ((MRIvox(mri_seg,i,j,k+1)==label) 
							&& ((MRIvox(mri_seg,i+1,j,k+1)==label)
									||(MRIvox(mri_seg,i,j+1,k+1)==label)))) continue;
					
					/* select the brigther path */
					dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
					ind1_i[0]=i+1;ind1_j[0]=j;ind1_k[0]=k;
					ind2_i[0]=i+1;ind2_j[0]=j+1;ind2_k[0]=k;
					
					dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
					ind1_i[1]=i+1;ind1_j[1]=j;ind1_k[1]=k;
					ind2_i[1]=i+1;ind2_j[1]=j;ind2_k[1]=k+1;
					
					dist[2]=MRIgetVoxVal(mri_orig,i,j+1,k,0)+ MRIgetVoxVal(mri_orig,i,j+1,k+1,0);
					ind1_i[2]=i;ind1_j[2]=j+1;ind1_k[2]=k;
					ind2_i[2]=i;ind2_j[2]=j+1;ind2_k[2]=k+1;
					
					dist[3]=MRIgetVoxVal(mri_orig,i,j+1,k,0)+ MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
					ind1_i[3]=i;ind1_j[3]=j+1;ind1_k[3]=k;
					ind2_i[3]=i+1;ind2_j[3]=j+1;ind2_k[3]=k;
					
					dist[4]=MRIgetVoxVal(mri_orig,i,j,k+1,0)+ MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
					ind1_i[4]=i;ind1_j[4]=j;ind1_k[4]=k+1;
					ind2_i[4]=i+1;ind2_j[4]=j;ind2_k[4]=k+1;

					dist[5]=MRIgetVoxVal(mri_orig,i,j,k+1,0)+ MRIgetVoxVal(mri_orig,i,j+1,k+1,0);
					ind1_i[5]=i;ind1_j[5]=j;ind1_k[5]=k+1;
					ind2_i[5]=i;ind2_j[5]=j+1;ind2_k[5]=k+1;

					/* find max path */
					refp=0;maxdist=dist[0];
					for(p = 1 ; p < 6 ; p++){
						if(maxdist<dist[p]){
							maxdist=dist[p];
							refp=p;
						}
					}
					/* assign value */
					
					if(MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label){
						MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
						//						fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
						nfound++;
					}
					if(MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label){
						MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
						//						fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
						nfound++;
					}
					//					if(nfound) exit(-1);
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}
	/* dealing with i+1,j+1,k-1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 1 ; k < mri_seg->depth; k++)
			for(j = 0 ; j < mri_seg->height-1 ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)!=label) continue;
					if(MRIvox(mri_seg,i+1,j+1,k-1)!=label) continue;

					/* find problematic configuration */
					if(((MRIvox(mri_seg,i+1,j,k)==label) 
							&& ((MRIvox(mri_seg,i+1,j+1,k)==label)
									||(MRIvox(mri_seg,i+1,j,k-1)==label))) || 
						 ((MRIvox(mri_seg,i,j+1,k)==label) 
							&& ((MRIvox(mri_seg,i,j+1,k-1)==label)
									||(MRIvox(mri_seg,i+1,j+1,k)==label))) || 
						 ((MRIvox(mri_seg,i,j,k-1)==label) 
							&& ((MRIvox(mri_seg,i+1,j,k-1)==label)
									||(MRIvox(mri_seg,i,j+1,k-1)==label)))) continue;
					
					/* select the brigther path */
					dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
					ind1_i[0]=i+1;ind1_j[0]=j;ind1_k[0]=k;
					ind2_i[0]=i+1;ind2_j[0]=j+1;ind2_k[0]=k;
					
					dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
					ind1_i[1]=i+1;ind1_j[1]=j;ind1_k[1]=k;
					ind2_i[1]=i+1;ind2_j[1]=j;ind2_k[1]=k-1;
					
					dist[2]=MRIgetVoxVal(mri_orig,i,j+1,k,0)+ MRIgetVoxVal(mri_orig,i,j+1,k-1,0);
					ind1_i[2]=i;ind1_j[2]=j+1;ind1_k[2]=k;
					ind2_i[2]=i;ind2_j[2]=j+1;ind2_k[2]=k-1;
					
					dist[3]=MRIgetVoxVal(mri_orig,i,j+1,k,0)+ MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
					ind1_i[3]=i;ind1_j[3]=j+1;ind1_k[3]=k;
					ind2_i[3]=i+1;ind2_j[3]=j+1;ind2_k[3]=k;
					
					dist[4]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+ MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
					ind1_i[4]=i;ind1_j[4]=j;ind1_k[4]=k-1;
					ind2_i[4]=i+1;ind2_j[4]=j;ind2_k[4]=k-1;

					dist[5]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+ MRIgetVoxVal(mri_orig,i,j+1,k-1,0);
					ind1_i[5]=i;ind1_j[5]=j;ind1_k[5]=k-1;
					ind2_i[5]=i;ind2_j[5]=j+1;ind2_k[5]=k-1;

					/* find max path */
					refp=0;maxdist=dist[0];
					for(p = 1 ; p < 6 ; p++){
						if(maxdist<dist[p]){
							maxdist=dist[p];
							refp=p;
						}
					}
					/* assign value */
					
					if(MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label){
						MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
						//						fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
						nfound++;
					}
					if(MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label){
						MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
						//						fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
						nfound++;
					}
					//					if(nfound) exit(-1);
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}
	

	/* dealing with i+1,j-1,k-1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 1 ; k < mri_seg->depth; k++)
			for(j = 1 ; j < mri_seg->height ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)!=label) continue;
					if(MRIvox(mri_seg,i+1,j-1,k-1)!=label) continue;

					/* find problematic configuration */
					if(((MRIvox(mri_seg,i+1,j,k)==label) 
							&& ((MRIvox(mri_seg,i+1,j-1,k)==label)
									||(MRIvox(mri_seg,i+1,j,k-1)==label))) || 
						 ((MRIvox(mri_seg,i,j-1,k)==label) 
							&& ((MRIvox(mri_seg,i,j-1,k-1)==label)
									||(MRIvox(mri_seg,i+1,j-1,k)==label))) || 
						 ((MRIvox(mri_seg,i,j,k-1)==label) 
							&& ((MRIvox(mri_seg,i+1,j,k-1)==label)
									||(MRIvox(mri_seg,i,j-1,k-1)==label)))) continue;
					
					/* select the brigther path */
					dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
					ind1_i[0]=i+1;ind1_j[0]=j;ind1_k[0]=k;
					ind2_i[0]=i+1;ind2_j[0]=j-1;ind2_k[0]=k;
					
					dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
					ind1_i[1]=i+1;ind1_j[1]=j;ind1_k[1]=k;
					ind2_i[1]=i+1;ind2_j[1]=j;ind2_k[1]=k-1;
					
					dist[2]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+ MRIgetVoxVal(mri_orig,i,j-1,k-1,0);
					ind1_i[2]=i;ind1_j[2]=j-1;ind1_k[2]=k;
					ind2_i[2]=i;ind2_j[2]=j-1;ind2_k[2]=k-1;
					
					dist[3]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+ MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
					ind1_i[3]=i;ind1_j[3]=j-1;ind1_k[3]=k;
					ind2_i[3]=i+1;ind2_j[3]=j-1;ind2_k[3]=k;
					
					dist[4]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+ MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
					ind1_i[4]=i;ind1_j[4]=j;ind1_k[4]=k-1;
					ind2_i[4]=i+1;ind2_j[4]=j;ind2_k[4]=k-1;

					dist[5]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+ MRIgetVoxVal(mri_orig,i,j-1,k-1,0);
					ind1_i[5]=i;ind1_j[5]=j;ind1_k[5]=k-1;
					ind2_i[5]=i;ind2_j[5]=j-1;ind2_k[5]=k-1;

					/* find max path */
					refp=0;maxdist=dist[0];
					for(p = 1 ; p < 6 ; p++){
						if(maxdist<dist[p]){
							maxdist=dist[p];
							refp=p;
						}
					}
					/* assign value */
					
					if(MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label){
						MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
						//						fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
						nfound++;
					}
					if(MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label){
						MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
						//						fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
						nfound++;
					}
					//					if(nfound) exit(-1);
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}

	/* dealing with i+1,j-1,k+1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 0 ; k < mri_seg->depth-1; k++)
			for(j = 1 ; j < mri_seg->height ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)!=label) continue;
					if(MRIvox(mri_seg,i+1,j-1,k+1)!=label) continue;

					/* find problematic configuration */
					if(((MRIvox(mri_seg,i+1,j,k)==label) 
							&& ((MRIvox(mri_seg,i+1,j-1,k)==label)
									||(MRIvox(mri_seg,i+1,j,k+1)==label))) || 
						 ((MRIvox(mri_seg,i,j-1,k)==label) 
							&& ((MRIvox(mri_seg,i,j-1,k+1)==label)
									||(MRIvox(mri_seg,i+1,j-1,k)==label))) || 
						 ((MRIvox(mri_seg,i,j,k+1)==label) 
							&& ((MRIvox(mri_seg,i+1,j,k+1)==label)
									||(MRIvox(mri_seg,i,j-1,k+1)==label)))) continue;
					
					/* select the brigther path */
					dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
					ind1_i[0]=i+1;ind1_j[0]=j;ind1_k[0]=k;
					ind2_i[0]=i+1;ind2_j[0]=j-1;ind2_k[0]=k;
					
					dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+ MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
					ind1_i[1]=i+1;ind1_j[1]=j;ind1_k[1]=k;
					ind2_i[1]=i+1;ind2_j[1]=j;ind2_k[1]=k+1;
					
					dist[2]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+ MRIgetVoxVal(mri_orig,i,j-1,k+1,0);
					ind1_i[2]=i;ind1_j[2]=j-1;ind1_k[2]=k;
					ind2_i[2]=i;ind2_j[2]=j-1;ind2_k[2]=k+1;
					
					dist[3]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+ MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
					ind1_i[3]=i;ind1_j[3]=j-1;ind1_k[3]=k;
					ind2_i[3]=i+1;ind2_j[3]=j-1;ind2_k[3]=k;
					
					dist[4]=MRIgetVoxVal(mri_orig,i,j,k+1,0)+ MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
					ind1_i[4]=i;ind1_j[4]=j;ind1_k[4]=k+1;
					ind2_i[4]=i+1;ind2_j[4]=j;ind2_k[4]=k+1;

					dist[5]=MRIgetVoxVal(mri_orig,i,j,k+1,0)+ MRIgetVoxVal(mri_orig,i,j-1,k+1,0);
					ind1_i[5]=i;ind1_j[5]=j;ind1_k[5]=k+1;
					ind2_i[5]=i;ind2_j[5]=j-1;ind2_k[5]=k+1;

					/* find max path */
					refp=0;maxdist=dist[0];
					for(p = 1 ; p < 6 ; p++){
						if(maxdist<dist[p]){
							maxdist=dist[p];
							refp=p;
						}
					}
					/* assign value */
					
					if(MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label){
						MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
						//						fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
						nfound++;
					}
					if(MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label){
						MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
						//						fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
						nfound++;
					}
					//					if(nfound) exit(-1);
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}

	return ntotal;
}

static int mriRemoveBackgroundCornerConfiguration(MRI *mri_seg, MRI *mri_orig, int label){
	static int niter=0;
	int i,j,k;
	int ntotal=0,nmodified,nfound,npass;
	
	niter++;
	fprintf(stderr,"\nIteration Number : %d",niter);

	/* dealing with i+1,j+1,k+1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 0 ; k < mri_seg->depth-1; k++)
			for(j = 0 ; j < mri_seg->height-1 ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)==label) continue;
					if(MRIvox(mri_seg,i+1,j+1,k+1)==label) continue;
						 /* find problematic configuration */
						 if((MRIvox(mri_seg,i+1,j,k)==label) && (MRIvox(mri_seg,i+1,j+1,k)==label) && 
								(MRIvox(mri_seg,i,j+1,k)==label) && (MRIvox(mri_seg,i,j+1,k+1)==label) &&
								(MRIvox(mri_seg,i+1,j,k)==label) && 
								(MRIvox(mri_seg,i+1,j,k+1)==label)){
							 if(MRIgetVoxVal(mri_orig,i,j,k,0)>MRIgetVoxVal(mri_orig,i+1,j+1,k+1,0))
								 MRIvox(mri_seg,i,j,k)=label;
							 else
								 MRIvox(mri_seg,i+1,j+1,k+1)=label;
							 nfound++;
						 }
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (++): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}
		/* dealing with i+1,j+1,k-1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 1 ; k < mri_seg->depth; k++)
			for(j = 0 ; j < mri_seg->height-1 ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)==label) continue;
					if(MRIvox(mri_seg,i+1,j+1,k-1)==label) continue;
						 /* find problematic configuration */
						 if((MRIvox(mri_seg,i+1,j,k)==label) && (MRIvox(mri_seg,i+1,j+1,k)==label) && 
								(MRIvox(mri_seg,i,j+1,k)==label) && (MRIvox(mri_seg,i,j+1,k-1)==label) &&
								(MRIvox(mri_seg,i+1,j,k)==label) && 
								(MRIvox(mri_seg,i+1,j,k-1)==label)){
							 if(MRIgetVoxVal(mri_orig,i,j,k,0)>MRIgetVoxVal(mri_orig,i+1,j+1,k-1,0))
								 MRIvox(mri_seg,i,j,k)=label;
							 else
								 MRIvox(mri_seg,i+1,j+1,k-1)=label;
							 nfound++;
						 }
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (+-): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}
	/* dealing with i+1,j-1,k-1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 1 ; k < mri_seg->depth; k++)
			for(j = 1 ; j < mri_seg->height ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)==label) continue;
					if(MRIvox(mri_seg,i+1,j-1,k-1)==label) continue;
						 /* find problematic configuration */
						 if((MRIvox(mri_seg,i+1,j,k)==label) && (MRIvox(mri_seg,i+1,j-1,k)==label) && 
								(MRIvox(mri_seg,i,j-1,k)==label) && (MRIvox(mri_seg,i,j-1,k-1)==label) &&
								(MRIvox(mri_seg,i+1,j,k)==label) && 
								(MRIvox(mri_seg,i+1,j,k-1)==label)){
							 if(MRIgetVoxVal(mri_orig,i,j,k,0)>MRIgetVoxVal(mri_orig,i+1,j-1,k-1,0))
								 MRIvox(mri_seg,i,j,k)=label;
							 else
								 MRIvox(mri_seg,i+1,j-1,k-1)=label;
							 nfound++;
						 }
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (--): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}
		/* dealing with i+1,j-1,k+1 */
	nfound=1;
	nmodified=0;
	npass=0;
	while(nfound){
		nfound=0;
		npass++;
		for(k = 0 ; k < mri_seg->depth-1; k++)
			for(j = 1 ; j < mri_seg->height ; j++)
				for(i = 0 ; i < mri_seg->width-1 ; i++){
					if(MRIvox(mri_seg,i,j,k)==label) continue;
					if(MRIvox(mri_seg,i+1,j-1,k+1)==label) continue;
						 /* find problematic configuration */
						 if((MRIvox(mri_seg,i+1,j,k)==label) && (MRIvox(mri_seg,i+1,j-1,k)==label) && 
								(MRIvox(mri_seg,i,j-1,k)==label) && (MRIvox(mri_seg,i,j-1,k+1)==label) &&
								(MRIvox(mri_seg,i+1,j,k)==label) && 
								(MRIvox(mri_seg,i+1,j,k+1)==label)){
							 if(MRIgetVoxVal(mri_orig,i,j,k,0)>MRIgetVoxVal(mri_orig,i+1,j-1,k+1,0))
								 MRIvox(mri_seg,i,j,k)=label;
							 else
								 MRIvox(mri_seg,i+1,j-1,k+1)=label;
							 nfound++;
						 }
				}
		nmodified += nfound;
		ntotal += nfound;
		fprintf(stderr,"\npass %3d (-+): %3d found - %3d modified     |    TOTAL: %3d",npass,nfound,nmodified,ntotal);
	}

	return ntotal;
}

int main(int argc, char *argv[])
{
	MRI *mri_seg,*mri_orig, *mri_seg_orig ;
	int niter=10,ntotal=0,nmodified,i,j,k,nvoxels, ac ;
	int label, nargs;
	char cmdline[CMD_LINE_LEN], **av ;
	int x, y, z ,convert=0;
	
  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_pretess.c,v 1.9 2005/12/02 20:44:45 segonne Exp $", "$Name:  $");

  make_cmd_version_string (argc, argv, "$Id: mri_pretess.c,v 1.9 2005/12/02 20:44:45 segonne Exp $", "$Name:  $", cmdline);

	Progname=argv[0];

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

	if(argc < 5) {
		fprintf(stderr,"\n\nUSAGE: mri_pretess filled label normalized newfilled\n\n");
		exit(-1);
	} 
	
	mri_seg=MRIread(argv[1]);
	if (!stricmp(argv[2], "wm"))
		label = USE_WM ;
	else
		label=atoi(argv[2]);
	if (mri_seg->type != MRI_UCHAR)
	{
		MRI *mri_tmp ;
		mri_tmp = MRIchangeType(mri_seg, MRI_UCHAR, 0, 255, 1) ;
		MRIfree(&mri_seg) ; mri_seg = mri_tmp ;
	}
	mri_orig=MRIread(argv[3]);
	if (mri_orig == NULL)
		ErrorExit(ERROR_NOFILE, "%s: could not open %s", Progname, argv[3]) ;

	if (label == USE_WM)
	{
		printf("binarizing input wm segmentation...\n") ;
		label = 128 ;
		mri_seg_orig = mri_seg ;
		mri_seg=MRIalloc(mri_seg_orig->width,mri_seg_orig->height,mri_seg_orig->depth,MRI_UCHAR);
		MRIbinarize(mri_seg_orig, mri_seg, WM_MIN_VAL, 0, label) ;
	}else if(mri_seg->type!=MRI_UCHAR){
		printf("converting input segmentation to MRI_UCHAR...\n") ;
		mri_seg_orig = mri_seg ;
		mri_seg=MRIalloc(mri_seg_orig->width,mri_seg_orig->height,mri_seg_orig->depth,MRI_UCHAR);
		for (x = 0 ; x < mri_seg->width ; x++)
			for (y = 0 ; y < mri_seg->height ; y++)
				for (z = 0 ; z < mri_seg->depth ; z++)
					if(((int)MRIgetVoxVal(mri_seg_orig,x,y,z,0))==label)
						MRIvox(mri_seg,x,y,z)=label;
					else
						MRIvox(mri_seg,x,y,z)=0;
		convert=1;
	}else
		mri_seg_orig = NULL ;
	

	fprintf(stderr,"\nAmbiguous edge configurations...");	
	while(niter--){
		nmodified=0;

		if(1) nmodified += mriRemoveEdgeConfiguration(mri_seg,mri_orig,label);
		if(corners) nmodified += mriRemoveCornerConfiguration(mri_seg,mri_orig,label);
		if(1) nmodified += mriRemoveBackgroundCornerConfiguration(mri_seg,mri_orig,label);
		if(nmodified==0) break;
		ntotal += nmodified;
	}

	nvoxels=0;
	for(k=0;k<mri_seg->depth;k++)
		for(j=0;j<mri_seg->height;j++)
			for(i=0;i<mri_seg->width;i++)
				if(MRIvox(mri_seg,i,j,k)==label) nvoxels++;

	fprintf(stderr,"\n\nTotal Number of Modified Voxels = %d (out of %d: %f)\n",ntotal,nvoxels,100.0*ntotal/nvoxels);

	fprintf(stderr,"\nWriting out volume...");

	if (mri_seg_orig)
	{
		if(!convert){
			for (x = 0 ; x < mri_seg->width ; x++)
				for (y = 0 ; y < mri_seg->height ; y++)
					for (z = 0 ; z < mri_seg->depth ; z++)
						if (MRIvox(mri_seg,x,y,z) > WM_MIN_VAL && MRIgetVoxVal(mri_seg_orig,x,y,z,0) < WM_MIN_VAL) 
							MRIsetVoxVal(mri_seg_orig, x, y, z,0,PRETESS_FILL);
			MRIfree(&mri_seg) ;
			mri_seg = mri_seg_orig ;
		}else{
			for (x = 0 ; x < mri_seg->width ; x++)
				for (y = 0 ; y < mri_seg->height ; y++)
					for (z = 0 ; z < mri_seg->depth ; z++)
						if (MRIvox(mri_seg,x,y,z)==label)
							MRIsetVoxVal(mri_seg_orig, x, y, z,0,label);
			MRIfree(&mri_seg) ;
			mri_seg = mri_seg_orig ;
		}
	}


	MRIaddCommandLine(mri_seg, cmdline) ;
	MRIwrite(mri_seg,argv[4]);

	fprintf(stderr,"\n\n");
  return 0;
}


static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!stricmp(option, "debug_voxel"))
	{
		Gx = atoi(argv[2]) ;
		Gy = atoi(argv[3]) ;
		Gz = atoi(argv[4]) ;
		nargs = 3 ;
		printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
	}
  else if (!stricmp(option, "nocorners"))
  {
		corners = 0 ;
		printf("no removal of corner configurations in addition to edge ones\n") ;
  }
	else switch (*option)
	{
	case 'W':
		Gdiag |= DIAG_WRITE ;
		break ;
  case '?':
  case 'U':
    usage_exit(1);
    break ;
	default:
    printf("unknown option %s\n", argv[1]) ;
		usage_exit(1) ;
    break ;
	}
	return(nargs) ;
}

static void 
usage_exit(int ecode)
{
	printf(
				 "usage: %s <wm vol> <label> <T1 vol> <output wm>\n",
				 Progname) ;
	printf("       where <wm vol> is usually wm.mgz\n") ;
	printf("             <label> is usually wm\n") ;
	printf("             <T1 vol> is brain.mgz\n") ;
	exit(ecode) ;
}
