#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"

#define IMGSIZE 256

/* MRIribbon determines the space between the inner and outer MRI surfaces provided, */
/* and creates a volume in mri_dst corresponding to the input format mri_src */
MRI *MRISribbon(MRI_SURFACE *inner_mris,MRI_SURFACE *outer_mris,MRI *mri_src,MRI *mri_dst)
{
  MRI *mri_inter;

  /* Allocate new MRI structures as needed */
  printf("MRIribbon: Creating new (_inter)MRI of %d, %d, %d...\n",mri_src->width,mri_src->height,mri_src->depth);
  mri_inter=MRIalloc(mri_src->width,mri_src->height,mri_src->depth,mri_src->type);
  MRIcopyHeader(mri_src, mri_inter);
  if (!mri_dst) {
    printf("MRIribbon: Creating new (_dst)MRI...\n");
    mri_dst=MRIalloc(mri_src->width,mri_src->height,mri_src->depth,mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  printf("Creating volume inside outer shell...\n");
  /* Create volume inside outer shell */
  /* Create shell corresponding to surface in MRI volume (includes outer shell in surface) */
  MRISshell(mri_src,outer_mris,mri_inter);
  MRISfloodoutside(mri_inter,mri_inter);
  MRISaccentuate(mri_inter,mri_inter,1,254);
  MRIcomplement(mri_inter,mri_inter);

  printf("Creating volume outside inner shell...\n");
  /* Create volume outside inner shell */
  MRISshell(mri_src,inner_mris,mri_dst);
  MRISfloodoutside(mri_dst,mri_dst);

  printf("Finding intersection of volumes...\n");
  /* Find intersection of volumes to create ribbon */
  MRIintersect(mri_inter,mri_dst,mri_dst);
  MRISaccentuate(mri_dst,mri_dst,1,255);

  printf("Done with intersection...\n");
  MRIfree(&mri_inter); /* Why segmentation fault? */

  return mri_dst;
}

/* MRISshell needs mri_src as an example of the format for the output, to match size and type.  The */
/* surface is recreated in the MRI space (mri_dst) from the tesselated surface (mris) as voxels of 255 */
MRI *MRISshell(MRI *mri_src,MRI_SURFACE *mris,MRI *mri_dst)
{
  int width,height,depth,i,j,imnr,fno,numu,numv,u,v;
  int imnr0;
  float ps,st,xx0,xx1,yy0,yy1,zz0,zz1,x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax;
  float px0,py0,pz0,px1,py1,pz1,px,py,pz;
  VERTEX *v_0,*v_1,*v_2;
  FACE *f;

  imnr0=mri_src->imnr0;
  st=mri_src->thick; /* slice thickness */
  ps=mri_src->ps;
  xx0=mri_src->xstart;
  xx1=mri_src->xend;
  yy0=mri_src->ystart;
  yy1=mri_src->yend;
  zz0=mri_src->zstart;
  zz1=mri_src->zend;

  /* Create new blank MRI or clear existing destination MRI */
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    printf("MRISshell: Creating new (_dst)MRI...\n");
    mri_dst = MRIalloc(width, height, depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }
  MRIclear(mri_dst);

  /* Fill each face in MRI volume */
  for (fno=0; fno<mris->nfaces; fno++) {
    /* Calculate (x,y,z) for each vertex for face */
    f = &mris->faces[fno];
    v_0 = &mris->vertices[f->v[0]];
    v_1 = &mris->vertices[f->v[1]];
    v_2 = &mris->vertices[f->v[2]];
    x0 = v_0->x;
    y0 = v_0->y;
    z0 = v_0->z;
    x1 = v_1->x;
    y1 = v_1->y;
    z1 = v_1->z;
    x2 = v_2->x;
    y2 = v_2->y;
    z2 = v_2->z;

    /* Calculate triangle side lengths */
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    /* Divide space between sides into numv parallel lines */
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
    numu = ceil(2*d0);
    numv = ceil(2*dmax);
    /* Fill each line in MRI volume */
    for (v=0; v<=numv; v++) {
      px0 = x0 + (x2-x0)*(float)v/(float)numv;
      py0 = y0 + (y2-y0)*(float)v/(float)numv;
      pz0 = z0 + (z2-z0)*(float)v/(float)numv;
      px1 = x1 + (x2-x1)*(float)v/(float)numv;
      py1 = y1 + (y2-y1)*(float)v/(float)numv;
      pz1 = z1 + (z2-z1)*(float)v/(float)numv;
      /* Fill each voxel on line in MRI volume */
      for (u=0; u<=numu; u++) {
        px = px0 + (px1-px0)*(float)u/(float)numu;
        py = py0 + (py1-py0)*(float)u/(float)numu;
        pz = pz0 + (pz1-pz0)*(float)u/(float)numu;
        /* Note mapping (x,y,z)<->(i,j,k) */
        imnr = (int)((py-yy0)/st+1.5-imnr0);
        i = (int)((xx1-px)/ps+0.5);
        j = (int)((zz1-pz)/ps+0.5);
        if (i>=0 && i<IMGSIZE && j>=0 && j<IMGSIZE && imnr>=0 && imnr<depth)
          MRIvox(mri_dst,i,j,imnr)=255;
      }
    }
  }
/* quick voxel test */
/*  MRIvox(mri_dst,32,32,32)=255;*/

  return mri_dst;
}

/* Floods MRI volume from outermost corners inward */
/* Fill with 1, boundary is anything but 0 and 1 */
MRI *MRISfloodoutside(MRI *mri_src,MRI *mri_dst)
{
  int newfilled,width,height,depth,i,j,k;

  /* Set MRI size */
  width=mri_src->width;
  height=mri_src->height;
  depth=mri_src->depth;

  /* Set seed voxel in corner of box */
  MRIvox(mri_dst,1,1,1)=1;

  newfilled=1;
  while (newfilled>0) {
    newfilled=0;

    for (i=1;i<width-1;i++)
      for (j=1;j<height-1;j++)
        for (k=1;k<depth-1;k++)
          if (MRIvox(mri_dst,i,j,k)==0)
            if (MRIvox(mri_dst,i,j,k-1)==1||
                MRIvox(mri_dst,i-1,j,k)==1||
                MRIvox(mri_dst,i,j-1,k)==1) {
                  MRIvox(mri_dst,i,j,k)=1;
                  newfilled++;
            }
    for (i=width-2;i>=1;i--)
      for (j=height-2;j>=1;j--)
        for (k=depth-2;k>=1;k--)
          if (MRIvox(mri_dst,i,j,k)==0)
            if (MRIvox(mri_dst,i,j,k+1)==1||
                MRIvox(mri_dst,i+1,j,k)==1||
                MRIvox(mri_dst,i,j+1,k)==1) {
                  MRIvox(mri_dst,i,j,k)=1;
                  newfilled++;
            }
  }

  return mri_dst;
}

MRI *MRISaccentuate(MRI *mri_src,MRI *mri_dst,int lo_thresh,int hi_thresh)
{
  int width,height,depth,i,j,k;

  printf("In MRISaccentuate.\n");
  width=mri_src->width;
  height=mri_src->height;
  depth=mri_src->depth;
  if (!mri_dst) {
    printf("MRISaccentuate: Creating new (_dst)MRI...\n");
    mri_dst = MRIalloc(width,height,depth,mri_src->type);
    MRIcopyHeader(mri_src,mri_dst);
  }

  printf("Loop.\n");
  for (k=0; k<depth; k++)
    for (j=0; j<height; j++)
      for (i=0; i<width; i++) {
        int vox=MRIvox(mri_src,i,j,k);
        MRIvox(mri_dst,i,j,k)=((vox>=lo_thresh)&&(vox<=hi_thresh))?255:0;
      }
  printf("Return.\n");
  return mri_dst;
}
