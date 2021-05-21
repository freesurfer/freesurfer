/**
 * @brief make sure a filled volume can be tessellated
 *
 * Changes white matter (WM) segmentation so that the neighbors of all
 * voxels labeled as WM have a face in common - no edges or corners allowed.
 */
/*
 * Original Author: Florent Segonne
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

const char *Progname;

int MRIaddEdgeVoxel(MRI *mri, int SetVal);

static void  print_help(void);
static int get_option(int argc, char *argv[]) ;
static int corners = 1, keep_edits = 0,test_edge = 0, test_edge_added;

/*-------------------------------------------------------------------*/
static int mriRemoveEdgeConfiguration(MRI *mri_seg, MRI *mri_orig, int label)
{
  static int niter=0;
  int i,j,k;
  int ntotal=0,nmodified,nfound,npass;

  niter++;
  fprintf(stderr,"\nIteration Number : %d",niter);

  /* dealing with xy-plane */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth ; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j+1,k)!=label)
          {
            continue;
          }
          if ((MRIvox(mri_seg,i,j+1,k)==label) ||
              (MRIvox(mri_seg,i+1,j,k)==label))
          {
            continue;
          }
          /* select the brigther voxel */
          if (MRIgetVoxVal(mri_orig,i,j+1,k,0) >
              MRIgetVoxVal(mri_orig,i+1,j,k,0))
          {
            MRIvox(mri_seg,i,j+1,k)=label;
          }
          else
          {
            MRIvox(mri_seg,i+1,j,k)=label;
          }
          nfound++;
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf
    (stderr,
     "\npass %3d (xy+): %3d found - %3d modified     |    TOTAL: %3d",
     npass,nfound,nmodified,ntotal);
  }
  /* dealing with xy-plane */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth ; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 1 ; i < mri_seg->width ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i-1,j+1,k)!=label)
          {
            continue;
          }
          if ((MRIvox(mri_seg,i,j+1,k)==label) ||
              (MRIvox(mri_seg,i-1,j,k)==label))
          {
            continue;
          }
          /* select the brigther voxel */
          if (MRIgetVoxVal(mri_orig,i,j+1,k,0) >
              MRIgetVoxVal(mri_orig,i-1,j,k,0))
          {
            MRIvox(mri_seg,i,j+1,k)=label;
          }
          else
          {
            MRIvox(mri_seg,i-1,j,k)=label;
          }
          nfound++;
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf
    (stderr,
     "\npass %3d (xy-): %3d found - %3d modified     |    TOTAL: %3d",
     npass,nfound,nmodified,ntotal);
  }

  /* dealing with yz-plane */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth-1 ; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 0 ; i < mri_seg->width ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i,j+1,k+1)!=label)
          {
            continue;
          }
          if ((MRIvox(mri_seg,i,j+1,k)==label) ||
              (MRIvox(mri_seg,i,j,k+1)==label))
          {
            continue;
          }
          /* select the brigther voxel */
          if (MRIgetVoxVal(mri_orig,i,j+1,k,0) >
              MRIgetVoxVal(mri_orig,i,j,k+1,0))
          {
            MRIvox(mri_seg,i,j+1,k)=label;
          }
          else
          {
            MRIvox(mri_seg,i,j,k+1)=label;
          }
          nfound++;
        }
    nmodified+=nfound;
    ntotal += nfound;
    fprintf
    (stderr,
     "\npass %3d (yz+): %3d found - %3d modified     |    TOTAL: %3d",
     npass,nfound,nmodified,ntotal);
  }
  /* dealing with yz-plane */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 1 ; k < mri_seg->depth ; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 0 ; i < mri_seg->width ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i,j+1,k-1)!=label)
          {
            continue;
          }
          if ((MRIvox(mri_seg,i,j+1,k)==label) ||
              (MRIvox(mri_seg,i,j,k-1)==label))
          {
            continue;
          }
          /* select the brigther voxel */
          if (MRIgetVoxVal(mri_orig,i,j+1,k,0) >
              MRIgetVoxVal(mri_orig,i,j,k-1,0))
          {
            MRIvox(mri_seg,i,j+1,k)=label;
          }
          else
          {
            MRIvox(mri_seg,i,j,k-1)=label;
          }
          nfound++;
        }
    nmodified+=nfound;
    ntotal += nfound;
    fprintf
    (stderr,
     "\npass %3d (yz-): %3d found - %3d modified     |    TOTAL: %3d",
     npass,nfound,nmodified,ntotal);
  }

  /* dealing with xz-plane */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth-1 ; k++)
      for (j = 0 ; j < mri_seg->height ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j,k+1)!=label)
          {
            continue;
          }
          if ((MRIvox(mri_seg,i+1,j,k)==label) ||
              (MRIvox(mri_seg,i,j,k+1)==label))
          {
            continue;
          }
          /* select the brigther voxel */
          if (MRIgetVoxVal(mri_orig,i+1,j,k,0) >
              MRIgetVoxVal(mri_orig,i,j,k+1,0))
          {
            MRIvox(mri_seg,i+1,j,k)=label;
          }
          else
          {
            MRIvox(mri_seg,i,j,k+1)=label;
          }
          nfound++;
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (xz+): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }
  /* dealing with xz-plane */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth-1 ; k++)
      for (j = 0 ; j < mri_seg->height ; j++)
        for (i = 1 ; i < mri_seg->width ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i-1,j,k+1)!=label)
          {
            continue;
          }
          if ((MRIvox(mri_seg,i-1,j,k)==label) ||
              (MRIvox(mri_seg,i,j,k+1)==label))
          {
            continue;
          }
          /* select the brigther voxel */
          if (MRIgetVoxVal(mri_orig,i-1,j,k,0) >
              MRIgetVoxVal(mri_orig,i,j,k+1,0))
          {
            MRIvox(mri_seg,i-1,j,k)=label;
          }
          else
          {
            MRIvox(mri_seg,i,j,k+1)=label;
          }
          nfound++;
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (xz-): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }
  return ntotal;
}

static int mriRemoveCornerConfiguration(MRI *mri_seg,
                                        MRI *mri_orig,
                                        int label)
{
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
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth-1; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j+1,k+1)!=label)
          {
            continue;
          }

          /* find problematic configuration */
          if (((MRIvox(mri_seg,i+1,j,k)==label)
               && ((MRIvox(mri_seg,i+1,j+1,k)==label)
                   ||(MRIvox(mri_seg,i+1,j,k+1)==label))) ||
              ((MRIvox(mri_seg,i,j+1,k)==label)
               && ((MRIvox(mri_seg,i,j+1,k+1)==label)
                   ||(MRIvox(mri_seg,i+1,j+1,k)==label))) ||
              ((MRIvox(mri_seg,i,j,k+1)==label)
               && ((MRIvox(mri_seg,i+1,j,k+1)==label)
                   ||(MRIvox(mri_seg,i,j+1,k+1)==label))))
          {
            continue;
          }

          /* select the brigther path */
          dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0) +
                  MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
          ind1_i[0]=i+1;
          ind1_j[0]=j;
          ind1_k[0]=k;
          ind2_i[0]=i+1;
          ind2_j[0]=j+1;
          ind2_k[0]=k;

          dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0) +
                  MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
          ind1_i[1]=i+1;
          ind1_j[1]=j;
          ind1_k[1]=k;
          ind2_i[1]=i+1;
          ind2_j[1]=j;
          ind2_k[1]=k+1;

          dist[2]=MRIgetVoxVal(mri_orig,i,j+1,k,0) +
                  MRIgetVoxVal(mri_orig,i,j+1,k+1,0);
          ind1_i[2]=i;
          ind1_j[2]=j+1;
          ind1_k[2]=k;
          ind2_i[2]=i;
          ind2_j[2]=j+1;
          ind2_k[2]=k+1;

          dist[3]=MRIgetVoxVal(mri_orig,i,j+1,k,0) +
                  MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
          ind1_i[3]=i;
          ind1_j[3]=j+1;
          ind1_k[3]=k;
          ind2_i[3]=i+1;
          ind2_j[3]=j+1;
          ind2_k[3]=k;

          dist[4]=MRIgetVoxVal(mri_orig,i,j,k+1,0) +
                  MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
          ind1_i[4]=i;
          ind1_j[4]=j;
          ind1_k[4]=k+1;
          ind2_i[4]=i+1;
          ind2_j[4]=j;
          ind2_k[4]=k+1;

          dist[5]=MRIgetVoxVal(mri_orig,i,j,k+1,0) +
                  MRIgetVoxVal(mri_orig,i,j+1,k+1,0);
          ind1_i[5]=i;
          ind1_j[5]=j;
          ind1_k[5]=k+1;
          ind2_i[5]=i;
          ind2_j[5]=j+1;
          ind2_k[5]=k+1;

          /* find max path */
          refp=0;
          maxdist=dist[0];
          for (p = 1 ; p < 6 ; p++)
          {
            if (maxdist<dist[p])
            {
              maxdist=dist[p];
              refp=p;
            }
          }
          /* assign value */

          if (MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label)
          {
            MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
            //      fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",
            // i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
            nfound++;
          }
          if (MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label)
          {
            MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
            //      fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",
            // i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
            nfound++;
          }
          //     if(nfound) exit(-1);
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }
  /* dealing with i+1,j+1,k-1 */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 1 ; k < mri_seg->depth; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j+1,k-1)!=label)
          {
            continue;
          }

          /* find problematic configuration */
          if (((MRIvox(mri_seg,i+1,j,k)==label)
               && ((MRIvox(mri_seg,i+1,j+1,k)==label)
                   ||(MRIvox(mri_seg,i+1,j,k-1)==label))) ||
              ((MRIvox(mri_seg,i,j+1,k)==label)
               && ((MRIvox(mri_seg,i,j+1,k-1)==label)
                   ||(MRIvox(mri_seg,i+1,j+1,k)==label))) ||
              ((MRIvox(mri_seg,i,j,k-1)==label)
               && ((MRIvox(mri_seg,i+1,j,k-1)==label)
                   ||(MRIvox(mri_seg,i,j+1,k-1)==label))))
          {
            continue;
          }

          /* select the brigther path */
          dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
          ind1_i[0]=i+1;
          ind1_j[0]=j;
          ind1_k[0]=k;
          ind2_i[0]=i+1;
          ind2_j[0]=j+1;
          ind2_k[0]=k;

          dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
          ind1_i[1]=i+1;
          ind1_j[1]=j;
          ind1_k[1]=k;
          ind2_i[1]=i+1;
          ind2_j[1]=j;
          ind2_k[1]=k-1;

          dist[2]=MRIgetVoxVal(mri_orig,i,j+1,k,0)+
                  MRIgetVoxVal(mri_orig,i,j+1,k-1,0);
          ind1_i[2]=i;
          ind1_j[2]=j+1;
          ind1_k[2]=k;
          ind2_i[2]=i;
          ind2_j[2]=j+1;
          ind2_k[2]=k-1;

          dist[3]=MRIgetVoxVal(mri_orig,i,j+1,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j+1,k,0);
          ind1_i[3]=i;
          ind1_j[3]=j+1;
          ind1_k[3]=k;
          ind2_i[3]=i+1;
          ind2_j[3]=j+1;
          ind2_k[3]=k;

          dist[4]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+
                  MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
          ind1_i[4]=i;
          ind1_j[4]=j;
          ind1_k[4]=k-1;
          ind2_i[4]=i+1;
          ind2_j[4]=j;
          ind2_k[4]=k-1;

          dist[5]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+
                  MRIgetVoxVal(mri_orig,i,j+1,k-1,0);
          ind1_i[5]=i;
          ind1_j[5]=j;
          ind1_k[5]=k-1;
          ind2_i[5]=i;
          ind2_j[5]=j+1;
          ind2_k[5]=k-1;

          /* find max path */
          refp=0;
          maxdist=dist[0];
          for (p = 1 ; p < 6 ; p++)
          {
            if (maxdist<dist[p])
            {
              maxdist=dist[p];
              refp=p;
            }
          }
          /* assign value */

          if (MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label)
          {
            MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
            //      fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",
            // i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
            nfound++;
          }
          if (MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label)
          {
            MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
            //      fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",
            // i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
            nfound++;
          }
          //     if(nfound) exit(-1);
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }

  /* dealing with i+1,j-1,k-1 */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 1 ; k < mri_seg->depth; k++)
      for (j = 1 ; j < mri_seg->height ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j-1,k-1)!=label)
          {
            continue;
          }

          /* find problematic configuration */
          if (((MRIvox(mri_seg,i+1,j,k)==label)
               && ((MRIvox(mri_seg,i+1,j-1,k)==label)
                   ||(MRIvox(mri_seg,i+1,j,k-1)==label))) ||
              ((MRIvox(mri_seg,i,j-1,k)==label)
               && ((MRIvox(mri_seg,i,j-1,k-1)==label)
                   ||(MRIvox(mri_seg,i+1,j-1,k)==label))) ||
              ((MRIvox(mri_seg,i,j,k-1)==label)
               && ((MRIvox(mri_seg,i+1,j,k-1)==label)
                   ||(MRIvox(mri_seg,i,j-1,k-1)==label))))
          {
            continue;
          }

          /* select the brigther path */
          dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
          ind1_i[0]=i+1;
          ind1_j[0]=j;
          ind1_k[0]=k;
          ind2_i[0]=i+1;
          ind2_j[0]=j-1;
          ind2_k[0]=k;

          dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
          ind1_i[1]=i+1;
          ind1_j[1]=j;
          ind1_k[1]=k;
          ind2_i[1]=i+1;
          ind2_j[1]=j;
          ind2_k[1]=k-1;

          dist[2]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+
                  MRIgetVoxVal(mri_orig,i,j-1,k-1,0);
          ind1_i[2]=i;
          ind1_j[2]=j-1;
          ind1_k[2]=k;
          ind2_i[2]=i;
          ind2_j[2]=j-1;
          ind2_k[2]=k-1;

          dist[3]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
          ind1_i[3]=i;
          ind1_j[3]=j-1;
          ind1_k[3]=k;
          ind2_i[3]=i+1;
          ind2_j[3]=j-1;
          ind2_k[3]=k;

          dist[4]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+
                  MRIgetVoxVal(mri_orig,i+1,j,k-1,0);
          ind1_i[4]=i;
          ind1_j[4]=j;
          ind1_k[4]=k-1;
          ind2_i[4]=i+1;
          ind2_j[4]=j;
          ind2_k[4]=k-1;

          dist[5]=MRIgetVoxVal(mri_orig,i,j,k-1,0)+
                  MRIgetVoxVal(mri_orig,i,j-1,k-1,0);
          ind1_i[5]=i;
          ind1_j[5]=j;
          ind1_k[5]=k-1;
          ind2_i[5]=i;
          ind2_j[5]=j-1;
          ind2_k[5]=k-1;

          /* find max path */
          refp=0;
          maxdist=dist[0];
          for (p = 1 ; p < 6 ; p++)
          {
            if (maxdist<dist[p])
            {
              maxdist=dist[p];
              refp=p;
            }
          }
          /* assign value */

          if (MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label)
          {
            MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
            //      fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",
            // i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
            nfound++;
          }
          if (MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label)
          {
            MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
            //      fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",
            // i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
            nfound++;
          }
          //     if(nfound) exit(-1);
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }

  /* dealing with i+1,j-1,k+1 */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth-1; k++)
      for (j = 1 ; j < mri_seg->height ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)!=label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j-1,k+1)!=label)
          {
            continue;
          }

          /* find problematic configuration */
          if (((MRIvox(mri_seg,i+1,j,k)==label)
               && ((MRIvox(mri_seg,i+1,j-1,k)==label)
                   ||(MRIvox(mri_seg,i+1,j,k+1)==label))) ||
              ((MRIvox(mri_seg,i,j-1,k)==label)
               && ((MRIvox(mri_seg,i,j-1,k+1)==label)
                   ||(MRIvox(mri_seg,i+1,j-1,k)==label))) ||
              ((MRIvox(mri_seg,i,j,k+1)==label)
               && ((MRIvox(mri_seg,i+1,j,k+1)==label)
                   ||(MRIvox(mri_seg,i,j-1,k+1)==label))))
          {
            continue;
          }

          /* select the brigther path */
          dist[0]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
          ind1_i[0]=i+1;
          ind1_j[0]=j;
          ind1_k[0]=k;
          ind2_i[0]=i+1;
          ind2_j[0]=j-1;
          ind2_k[0]=k;

          dist[1]=MRIgetVoxVal(mri_orig,i+1,j,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
          ind1_i[1]=i+1;
          ind1_j[1]=j;
          ind1_k[1]=k;
          ind2_i[1]=i+1;
          ind2_j[1]=j;
          ind2_k[1]=k+1;

          dist[2]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+
                  MRIgetVoxVal(mri_orig,i,j-1,k+1,0);
          ind1_i[2]=i;
          ind1_j[2]=j-1;
          ind1_k[2]=k;
          ind2_i[2]=i;
          ind2_j[2]=j-1;
          ind2_k[2]=k+1;

          dist[3]=MRIgetVoxVal(mri_orig,i,j-1,k,0)+
                  MRIgetVoxVal(mri_orig,i+1,j-1,k,0);
          ind1_i[3]=i;
          ind1_j[3]=j-1;
          ind1_k[3]=k;
          ind2_i[3]=i+1;
          ind2_j[3]=j-1;
          ind2_k[3]=k;

          dist[4]=MRIgetVoxVal(mri_orig,i,j,k+1,0)+
                  MRIgetVoxVal(mri_orig,i+1,j,k+1,0);
          ind1_i[4]=i;
          ind1_j[4]=j;
          ind1_k[4]=k+1;
          ind2_i[4]=i+1;
          ind2_j[4]=j;
          ind2_k[4]=k+1;

          dist[5]=MRIgetVoxVal(mri_orig,i,j,k+1,0)+
                  MRIgetVoxVal(mri_orig,i,j-1,k+1,0);
          ind1_i[5]=i;
          ind1_j[5]=j;
          ind1_k[5]=k+1;
          ind2_i[5]=i;
          ind2_j[5]=j-1;
          ind2_k[5]=k+1;

          /* find max path */
          refp=0;
          maxdist=dist[0];
          for (p = 1 ; p < 6 ; p++)
          {
            if (maxdist<dist[p])
            {
              maxdist=dist[p];
              refp=p;
            }
          }
          /* assign value */

          if (MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])!=label)
          {
            MRIvox(mri_seg,ind1_i[refp],ind1_j[refp],ind1_k[refp])=label;
            //      fprintf(stderr,"(%d,%d,%d)&-(%d,%d,%d)",
            // i,j,k,ind1_i[refp],ind1_j[refp],ind1_k[refp]);
            nfound++;
          }
          if (MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])!=label)
          {
            MRIvox(mri_seg,ind2_i[refp],ind2_j[refp],ind2_k[refp])=label;
            //      fprintf(stderr,"+(%d,%d,%d)&-(%d,%d,%d)-",
            // i,j,k,ind2_i[refp],ind2_j[refp],ind2_k[refp]);
            nfound++;
          }
          //     if(nfound) exit(-1);
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (+++): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }

  return ntotal;
}

static int mriRemoveBackgroundCornerConfiguration(MRI *mri_seg,
    MRI *mri_orig,
    int label)
{
  static int niter=0;
  int i,j,k;
  int ntotal=0,nmodified,nfound,npass;

  niter++;
  fprintf(stderr,"\nIteration Number : %d",niter);

  /* dealing with i+1,j+1,k+1 */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth-1; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)==label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j+1,k+1)==label)
          {
            continue;
          }
          /* find problematic configuration */
          if ((MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j+1,k)==label) &&
              (MRIvox(mri_seg,i,j+1,k)==label) &&
              (MRIvox(mri_seg,i,j+1,k+1)==label) &&
              (MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j,k+1)==label))
          {
            if (MRIgetVoxVal(mri_orig,i,j,k,0) >
                MRIgetVoxVal(mri_orig,i+1,j+1,k+1,0))
            {
              MRIvox(mri_seg,i,j,k)=label;
            }
            else
            {
              MRIvox(mri_seg,i+1,j+1,k+1)=label;
            }
            nfound++;
          }
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (++): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }
  /* dealing with i+1,j+1,k-1 */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 1 ; k < mri_seg->depth; k++)
      for (j = 0 ; j < mri_seg->height-1 ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)==label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j+1,k-1)==label)
          {
            continue;
          }
          /* find problematic configuration */
          if ((MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j+1,k)==label) &&
              (MRIvox(mri_seg,i,j+1,k)==label) &&
              (MRIvox(mri_seg,i,j+1,k-1)==label) &&
              (MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j,k-1)==label))
          {
            if (MRIgetVoxVal(mri_orig,i,j,k,0) >
                MRIgetVoxVal(mri_orig,i+1,j+1,k-1,0))
            {
              MRIvox(mri_seg,i,j,k)=label;
            }
            else
            {
              MRIvox(mri_seg,i+1,j+1,k-1)=label;
            }
            nfound++;
          }
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (+-): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }
  /* dealing with i+1,j-1,k-1 */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 1 ; k < mri_seg->depth; k++)
      for (j = 1 ; j < mri_seg->height ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)==label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j-1,k-1)==label)
          {
            continue;
          }
          /* find problematic configuration */
          if ((MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j-1,k)==label) &&
              (MRIvox(mri_seg,i,j-1,k)==label) &&
              (MRIvox(mri_seg,i,j-1,k-1)==label) &&
              (MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j,k-1)==label))
          {
            if (MRIgetVoxVal(mri_orig,i,j,k,0) >
                MRIgetVoxVal(mri_orig,i+1,j-1,k-1,0))
            {
              MRIvox(mri_seg,i,j,k)=label;
            }
            else
            {
              MRIvox(mri_seg,i+1,j-1,k-1)=label;
            }
            nfound++;
          }
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (--): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }
  /* dealing with i+1,j-1,k+1 */
  nfound=1;
  nmodified=0;
  npass=0;
  while (nfound)
  {
    nfound=0;
    npass++;
    for (k = 0 ; k < mri_seg->depth-1; k++)
      for (j = 1 ; j < mri_seg->height ; j++)
        for (i = 0 ; i < mri_seg->width-1 ; i++)
        {
          if (MRIvox(mri_seg,i,j,k)==label)
          {
            continue;
          }
          if (MRIvox(mri_seg,i+1,j-1,k+1)==label)
          {
            continue;
          }
          /* find problematic configuration */
          if ((MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j-1,k)==label) &&
              (MRIvox(mri_seg,i,j-1,k)==label) &&
              (MRIvox(mri_seg,i,j-1,k+1)==label) &&
              (MRIvox(mri_seg,i+1,j,k)==label) &&
              (MRIvox(mri_seg,i+1,j,k+1)==label))
          {
            if (MRIgetVoxVal(mri_orig,i,j,k,0) >
                MRIgetVoxVal(mri_orig,i+1,j-1,k+1,0))
            {
              MRIvox(mri_seg,i,j,k)=label;
            }
            else
            {
              MRIvox(mri_seg,i+1,j-1,k+1)=label;
            }
            nfound++;
          }
        }
    nmodified += nfound;
    ntotal += nfound;
    fprintf(stderr,
            "\npass %3d (-+): %3d found - %3d modified     |    TOTAL: %3d",
            npass,nfound,nmodified,ntotal);
  }

  return ntotal;
}

/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  MRI *mri_seg,*mri_orig, *mri_seg_orig, *mri_old = NULL;
  int niter=100,ntotal=0,nmodified,i,j,k,nvoxels, ac ;
  int label, nargs,nvox;
  char **av ;
  int x, y, z ,convert=0;

  nargs = handleVersionOption(argc, argv, "mri_pretess");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }

  std::string cmdline = getAllInfo(argc, argv, "mri_pretess");

  Progname=argv[0];

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
  {
    print_help();
    exit(-1);
  }

  mri_seg=MRIread(argv[1]);
  if (mri_seg == NULL)
  {
    ErrorExit(ERROR_NOFILE,"%s: could not open %s", Progname, argv[1]) ;
  }
  if (test_edge)
  {
    printf("Adding WM edge to input segmentation.\n");
    test_edge_added = MRIaddEdgeVoxel(mri_seg, WM_EDITED_ON_VAL);
    //test_edge_added = MRIaddEdgeVoxel(mri_seg, 110);
    if (!test_edge_added)
    {
      printf("ERROR: could not add test edge\n");
      exit(1);
    }
  }
  if (!stricmp(argv[2], "wm"))
  {
    label = USE_WM ;
  }
  else
  {
    label=atoi(argv[2]);
  }
  if (mri_seg->type != MRI_UCHAR)
  {
    MRI *mri_tmp ;
    mri_tmp = MRIchangeType(mri_seg, MRI_UCHAR, 0, 255, 1) ;
    MRIfree(&mri_seg) ;
    mri_seg = mri_tmp ;
  }
  mri_orig=MRIread(argv[3]);
  if (mri_orig == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open %s", Progname, argv[3]) ;
  }

  if (mri_orig->width != mri_seg->width ||
      mri_orig->height != mri_seg->height ||
      mri_orig->depth != mri_seg->depth)
    ErrorExit(ERROR_BADPARM, 
	      "intensity (%d, %d, %d) and segmentation volume"
	      " (%d x %d x %d) don't match",
	      mri_orig->width, mri_orig->height, mri_orig->depth,
	      mri_seg->width, mri_seg->height, mri_seg->depth) ;
  if (label == USE_WM)
  {
    printf("binarizing input wm segmentation...\n") ;
    label = 128 ;
    mri_seg_orig = mri_seg ;
    mri_seg=MRIalloc(mri_seg_orig->width,
                     mri_seg_orig->height,
                     mri_seg_orig->depth,
                     MRI_UCHAR);
    MRIbinarize(mri_seg_orig, mri_seg, WM_MIN_VAL, 0, label) ;
  }
  else if (mri_seg->type!=MRI_UCHAR)
  {
    printf("converting input segmentation to MRI_UCHAR...\n") ;
    mri_seg_orig = mri_seg ;
    mri_seg=MRIalloc(mri_seg_orig->width,
                     mri_seg_orig->height,
                     mri_seg_orig->depth,
                     MRI_UCHAR);
    for (x = 0 ; x < mri_seg->width ; x++)
      for (y = 0 ; y < mri_seg->height ; y++)
        for (z = 0 ; z < mri_seg->depth ; z++)
          if (((int)MRIgetVoxVal(mri_seg_orig,x,y,z,0))==label)
          {
            MRIvox(mri_seg,x,y,z)=label;
          }
          else
          {
            MRIvox(mri_seg,x,y,z)=0;
          }
    convert=1;
  }
  else
  {
    mri_seg_orig = NULL ;
  }

  printf("Ambiguous edge configurations... \n");
  while (niter--)
  {
    nmodified=0;

    if (1)
    {
      nmodified += mriRemoveEdgeConfiguration(mri_seg,mri_orig,label);
    }
    if (corners) nmodified +=
        mriRemoveCornerConfiguration(mri_seg,mri_orig,label);
    if (1) nmodified +=
        mriRemoveBackgroundCornerConfiguration(mri_seg,mri_orig,label);
    if (nmodified==0)
    {
      break;
    }
    ntotal += nmodified;
  }

  nvoxels=0;
  for (k=0; k<mri_seg->depth; k++)
    for (j=0; j<mri_seg->height; j++)
      for (i=0; i<mri_seg->width; i++)
        if (MRIvox(mri_seg,i,j,k)==label)
        {
          nvoxels++;
        }

  fprintf(stderr,"\n\nTotal Number of Modified Voxels = %d (out of %d: %f)\n",
          ntotal,nvoxels,100.0*ntotal/nvoxels);
  printf("\n");

  //fprintf(stderr,"\nWriting out volume...");
  //if(argc==6) MRIwrite(mri_seg,argv[5]);
  if (mri_seg_orig)
  {
    if (!convert)
    {
      for (x = 0 ; x < mri_seg->width ; x++)
        for (y = 0 ; y < mri_seg->height ; y++)
          for (z = 0 ; z < mri_seg->depth ; z++)
            if (MRIvox(mri_seg,x,y,z) > WM_MIN_VAL &&
                MRIgetVoxVal(mri_seg_orig,x,y,z,0) < WM_MIN_VAL)
            {
              MRIsetVoxVal(mri_seg_orig, x, y, z,0,PRETESS_FILL);
            }
      MRIfree(&mri_seg) ;
      mri_seg = mri_seg_orig ;
    }
    else
    {
      for (x = 0 ; x < mri_seg->width ; x++)
        for (y = 0 ; y < mri_seg->height ; y++)
          for (z = 0 ; z < mri_seg->depth ; z++)
            if (MRIvox(mri_seg,x,y,z)==label)
            {
              MRIsetVoxVal(mri_seg_orig, x, y, z,0,label);
            }
      MRIfree(&mri_seg) ;
      mri_seg = mri_seg_orig ;
    }
  }

  /*-------------- Keep Edits ---------------------------*/
  if (keep_edits)
  {
    printf("Searching for edits to keep ...\n");

    mri_old = MRIread(argv[4]);
    if (!mri_old)
    {
      ErrorPrintf(ERROR_NOFILE, "%s: could not read file %s to preserve edits",
                  Progname, argv[4]);
      exit(1);
    }

    nvox = MRIcopyLabel(mri_old, mri_seg, WM_EDITED_ON_VAL) ;
    printf("  kept %d WM ON voxels\n",nvox);
    nvox = MRIcopyLabel(mri_old, mri_seg, WM_EDITED_OFF_VAL) ;
    printf("  kept %d WM OFF voxels\n",nvox);
    MRIfree(&mri_old) ;
    printf("\n");
  }

  MRIaddCommandLine(mri_seg, cmdline) ;

  if (! test_edge)
  {
    MRIwrite(mri_seg,argv[4]);
  }
  else
  {
    printf("\n");
    printf("NOT saving output because called with -test.\n");
    printf("\n");
  }

  printf("mri_pretess done\n");
  printf("\n");
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
  else if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help();
    exit(1);
  }
  else if (!stricmp(option, "nocorners"))
  {
    corners = 0 ;
    printf("no removal of corner configurations in addition to edge ones\n") ;
  }
  else if (!stricmp(option, "keep"))
  {
    keep_edits = 1;
    printf("keeping edits\n") ;
  }
  else if (!stricmp(option, "test"))
  {
    test_edge = 1;
    printf("testing, output will not be saved\n") ;
  }
  else switch (*option)
    {
    case 'W':
      Gdiag |= DIAG_WRITE ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_help();
      exit(1);
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }
  return(nargs) ;
}

#include "mri_pretess.help.xml.h"
static void print_help(void)
{
  outputHelpXml(mri_pretess_help_xml,
                mri_pretess_help_xml_len);
}

/*------------------------------------------------------------------------
  MRIaddEdgeVoxel() - creates a segmentation situation which should be
  removed by pretess. Finds a voxel that is not WM and has no
  face-neighbors that are WM but does have edge neighbors that are WM.
  It then sets this voxel to SetVal. This voxel should be removed by
  mri_pretess (unless maybe if keep_edits is set). Returns 1 if it
  found a voxel to add, 0 otherwise.
  ------------------------------------------------------------------------*/
int MRIaddEdgeVoxel(MRI *mri, int SetVal)
{
  int set,c,r,s;
  double v;

  for (c=1; c < mri->width-1; c++)
  {
    for (r=1; r < mri->height-1; r++)
    {
      for (s=1; s < mri->depth-1; s++)
      {
        v = MRIgetVoxVal(mri,c,r,s,0);

        // Skip if center or face-neighbor is labeled as WM
        if (v >= WM_MIN_VAL)
        {
          continue;
        }
        if (MRIgetVoxVal(mri,c-1,r,s,0) >= WM_MIN_VAL)
        {
          continue;
        }
        if (MRIgetVoxVal(mri,c+1,r,s,0) >= WM_MIN_VAL)
        {
          continue;
        }
        if (MRIgetVoxVal(mri,c,r-1,s,0) >= WM_MIN_VAL)
        {
          continue;
        }
        if (MRIgetVoxVal(mri,c,r+1,s,0) >= WM_MIN_VAL)
        {
          continue;
        }
        if (MRIgetVoxVal(mri,c,r,s-1,0) >= WM_MIN_VAL)
        {
          continue;
        }
        if (MRIgetVoxVal(mri,c,r,s+1,0) >= WM_MIN_VAL)
        {
          continue;
        }

        // Set center as WM if any edge-neighbor is labeled as WM
        set=0;
        if (MRIgetVoxVal(mri,c-1,r-1,s,0) >= WM_MIN_VAL)
        {
          set=1;
        }
        if (MRIgetVoxVal(mri,c+1,r-1,s,0) >= WM_MIN_VAL)
        {
          set=1;
        }
        if (MRIgetVoxVal(mri,c+1,r+1,s,0) >= WM_MIN_VAL)
        {
          set=1;
        }
        if (MRIgetVoxVal(mri,c-1,r+1,s,0) >= WM_MIN_VAL)
        {
          set=1;
        }
        // Should go thru all the other edges, but I'm lazy

        if (set)
        {
          printf("  Adding test edge at %d %d %d, val=%d \n",c,r,s,SetVal);
          MRIsetVoxVal(mri,c,r,s,0,SetVal);
          return(1);
        }
      }
    }
  }
  printf("Could not find a test voxel to set\n");
  return(0);
}
