#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "cma.h"

#define IMGSIZE 256

#define subvoxmask(i,j,k) (1<<((i)+((j)<<1)+((k)<<2)))

unsigned char SubFill(unsigned char vox,int i,int j,int k);
MRI *MRIbitwiseand(MRI *mri1,MRI *mri2,MRI *mri_dst);
MRI *MRIbitwiseor(MRI *mri1,MRI *mri2,MRI *mri_dst);
MRI *MRIbitwisenot(MRI *mri_src,MRI *mri_dst);
MRI *MRImajority(MRI *mri_src,MRI *mri_dst);
MRI *MRImergecortexwhitecma(MRI *mri_cortex,MRI *mri_white,MRI *mri_cma,MRI *mri_dst,int hemisphere);

/* MRIribbon determines the space between the inner and outer MRI surfaces provided, */
/* and creates a volume in mri_dst corresponding to the input format mri_src */
MRI *MRISribbon(MRI_SURFACE *inner_mris,MRI_SURFACE *outer_mris,MRI *mri_src,MRI *mri_dst)
{
  MRI *mri_inter;

  /* Allocate new MRI structures as needed */
  printf("MRISribbon: Creating new (_inter)MRI of %d, %d,%d...\n",mri_src->width,mri_src->height,mri_src->depth);  
  mri_inter=MRIalloc(mri_src->width,mri_src->height,mri_src->depth,mri_src->type);
  MRIcopyHeader(mri_src, mri_inter);
  if (!mri_dst) {
    printf("MRISribbon: Creating new (_dst)MRI...\n");
    mri_dst=MRIalloc(mri_src->width,mri_src->height,mri_src->depth,mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  printf("Creating volume inside outer shell...\n");
  /* Create volume inside outer shell */
  /* Create shell corresponding to surface in MRI volume (includes outer shell in surface) */
  MRISpartialshell(mri_src,outer_mris,mri_inter);
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
  MRIfree(&mri_inter);

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
/*    printf("MRISshell: Creating new (_dst)MRI...\n");*/
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

  width=mri_src->width;
  height=mri_src->height;
  depth=mri_src->depth;
  if (!mri_dst) {
/*    printf("MRISaccentuate: Creating new (_dst)MRI...\n");*/
    mri_dst = MRIalloc(width,height,depth,mri_src->type);
    MRIcopyHeader(mri_src,mri_dst);
  }

  for (k=0; k<depth; k++)
    for (j=0; j<height; j++)
      for (i=0; i<width; i++) {
        int vox=MRIvox(mri_src,i,j,k);
        MRIvox(mri_dst,i,j,k)=((vox>=lo_thresh)&&(vox<=hi_thresh))?255:0;
      }
  return mri_dst;
}

/* Set mri_dst voxel to 255 for every mri_src voxel for which the majority
   of subvoxels is set. */
MRI *MRImajority(MRI *mri_src,MRI *mri_dst)
{
  int width,height,depth,i,j,k,vox;

  width=mri_src->width;
  height=mri_src->height;
  depth=mri_src->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width,height,depth,mri_src->type);
    MRIcopyHeader(mri_src,mri_dst);
  }

  for (k=0; k<depth; k++)
    for (j=0; j<height; j++)
      for (i=0; i<width; i++) {
        vox=MRIvox(mri_src,i,j,k);
        MRIvox(mri_dst,i,j,k)=((((vox>>7)&1)+((vox>>6)&1)+((vox>>5)&1)+
          ((vox>>4)&1)+((vox>>3)&1)+((vox>>2)&1)+((vox>>1)&1)+(vox&1))>4)?255:0;
      }

  return mri_dst;
}

MRI *MRIbitwisenot(MRI *mri_src,MRI *mri_dst)
{
  int width,height,depth,i,j,k,vox;

  width=mri_src->width;
  height=mri_src->height;
  depth=mri_src->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width,height,depth,mri_src->type);
    MRIcopyHeader(mri_src,mri_dst);
  }

  for (k=0; k<depth; k++)
    for (j=0; j<height; j++)
      for (i=0; i<width; i++) {
        vox=MRIvox(mri_src,i,j,k);
        MRIvox(mri_dst,i,j,k)=~vox;
      }
  return mri_dst;
}

/* The following is the same as above, but adapted for 2 x 2 x 2 resolution */
/* per voxel. */

/* Partial shell fills the voxel space with a superresolution surface - */
/* each voxel is divided into 8 subvoxels represented by the 8 bits of the */
/* unsigned char in the MRI structure. */
MRI *MRISpartialshell(MRI *mri_src,MRI_SURFACE *mris,MRI *mri_dst)
{
  int width,height,depth,i,j,imnr,isub,jsub,isubmnr,fno,numu,numv,u,v;
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
    for (v=0; v<=numv*2; v++) {
      px0 = x0 + (x2-x0)*(float)v/(float)numv/2;
      py0 = y0 + (y2-y0)*(float)v/(float)numv/2;
      pz0 = z0 + (z2-z0)*(float)v/(float)numv/2;
      px1 = x1 + (x2-x1)*(float)v/(float)numv/2;
      py1 = y1 + (y2-y1)*(float)v/(float)numv/2;
      pz1 = z1 + (z2-z1)*(float)v/(float)numv/2;
      /* Fill each voxel on line in MRI volume */
      for (u=0; u<=numu*2; u++) {
        px = px0 + (px1-px0)*(float)u/(float)numu/2;
        py = py0 + (py1-py0)*(float)u/(float)numu/2;
        pz = pz0 + (pz1-pz0)*(float)u/(float)numu/2;
        /* Note mapping (x,y,z)<->(i,j,k) */
        imnr = (int)((py-yy0)/st+1.5-imnr0);
        i = (int)((xx1-px)/ps+0.5);
        j = (int)((zz1-pz)/ps+0.5);
        if (i>=0 && i<IMGSIZE && j>=0 && j<IMGSIZE && imnr>=0 && imnr<depth) {
          /* Each voxel has 8 subvoxels, represented by the 8 bits of the unsigned char. */
          isubmnr = ((int)(((py-yy0)/st+1.5-imnr0)*2))-((int)((py-yy0)/st+1.5-imnr0))*2;
          isub = ((int)(((xx1-px)/ps+0.5)*2))-((int)((xx1-px)/ps+0.5))*2;
          jsub = ((int)(((zz1-pz)/ps+0.5)*2))-((int)((zz1-pz)/ps+0.5))*2;
    /* (isubmnr, isub, jsub) should be in the range (0..1, 0..1, 0..1) */
          /* Assume that the initial value for all voxels is zero */
    MRIvox(mri_dst,i,j,imnr)=MRIvox(mri_dst,i,j,imnr)|subvoxmask(isub,jsub,isubmnr);
        }
      }
    }
  }

  return mri_dst;
}

MRI *MRIbitwiseor(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int width,height,depth,x,y,z;
  BUFTYPE *p1,*p2,*pdst,v1,v2;

  width = mri1->width;
  height = mri1->height;
  depth = mri1->depth;

  if (!mri_dst)
    mri_dst = MRIclone(mri1,NULL);
 
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);
      p1 = &MRIvox(mri1, 0, y, z);
      p2 = &MRIvox(mri2, 0, y, z);
      for (x = 0 ; x < width ; x++) {
        v1 = *p1++;
        v2 = *p2++;
        *pdst++ = v1 | v2;
      }
    }
  }
  return(mri_dst);
}

MRI *MRIbitwiseand(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int width,height,depth,x,y,z;
  BUFTYPE *p1,*p2,*pdst,v1,v2;

  width = mri1->width;
  height = mri1->height;
  depth = mri1->depth;

  if (!mri_dst)
    mri_dst = MRIclone(mri1,NULL);
 
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);
      p1 = &MRIvox(mri1, 0, y, z);
      p2 = &MRIvox(mri2, 0, y, z);
      for (x = 0 ; x < width ; x++) {
        v1 = *p1++;
        v2 = *p2++;
        *pdst++ = v1 & v2;
      }
    }
  }
  return(mri_dst);
}

/* Floods MRI volume from outermost corners inward. */
/* Upon return, mri_dst contains the filled volume NOT including the shell. */
/* mri_src and mri_dst cannot be the same volume! */
MRI *MRISpartialfloodoutside(MRI *mri_src,MRI *mri_dst)
{
  int newfilled,width,height,depth,i,j,k,is,js,ks,isub,jsub,ksub;

  /* Set MRI size */
  width=mri_src->width;
  height=mri_src->height;
  depth=mri_src->depth;

  /* Set seed voxel in corner of voxel in corner of box */
  MRIvox(mri_dst,1,1,1)=MRIvox(mri_dst,1,1,1)|subvoxmask(0,0,0);

  newfilled=1;
  while (newfilled>0) {
    newfilled=0;

printf("    (left to right)\n");
    for (is=2;is<2*width-1;is++) {
      for (js=2;js<2*height-1;js++)
        for (ks=2;ks<2*depth-1;ks++) {
          i=is/2; j=js/2; k=ks/2;
          isub=is%2; jsub=js%2; ksub=ks%2;
/*printf("i,j,k,isub,jsub,ksub = %d,%d,%d,%d,%d,%d\n",i,j,k,isub,jsub,ksub);
printf("curr. vox = %d,%d\n",(MRIvox(mri_src,i,j,k)&subvoxmask(isub,jsub,ksub)),(MRIvox(mri_dst,i,j,k)&subvoxmask(isub,jsub,ksub)));
printf("adj. voxels = %d,%d,%d\n",(MRIvox(mri_dst,i,j,(ks-1)/2)&subvoxmask(isub,jsub,(ks-1)%2)),
  (MRIvox(mri_dst,(is-1)/2,j,k)&subvoxmask((is-1)%2,jsub,ksub)),
  (MRIvox(mri_dst,(is-1)/2,j,k)&subvoxmask((is-1)%2,jsub,ksub)));*/
          if (((MRIvox(mri_src,i,j,k)&subvoxmask(isub,jsub,ksub))==0)
             &&((MRIvox(mri_dst,i,j,k)&subvoxmask(isub,jsub,ksub))==0))
            if (((MRIvox(mri_dst,i,j,(ks-1)/2)&subvoxmask(isub,jsub,(ks-1)%2))>0) ||
                ((MRIvox(mri_dst,(is-1)/2,j,k)&subvoxmask((is-1)%2,jsub,ksub))>0) ||
                ((MRIvox(mri_dst,i,(js-1)/2,k)&subvoxmask(isub,(js-1)%2,ksub))>0)) {
                  MRIvox(mri_dst,i,j,k)=MRIvox(mri_dst,i,j,k)|subvoxmask(isub,jsub,ksub);
                  newfilled++;
          }
        }
      }
printf("    (right to left)\n");
    for (is=2*width-2;is>=1;is--) {
      for (js=2*height-2;js>=1;js--)
        for (ks=2*depth-2;ks>=1;ks--) {
          i=is/2; j=js/2; k=ks/2;
          isub=is%2; jsub=js%2; ksub=ks%2;
          if (((MRIvox(mri_src,i,j,k)&subvoxmask(isub,jsub,ksub))==0)
             &&((MRIvox(mri_dst,i,j,k)&subvoxmask(isub,jsub,ksub))==0))
            if (((MRIvox(mri_dst,i,j,(ks+1)/2)&subvoxmask(isub,jsub,(ks+1)%2))>0) ||
                ((MRIvox(mri_dst,(is+1)/2,j,k)&subvoxmask((is+1)%2,jsub,ksub))>0) ||
                ((MRIvox(mri_dst,i,(js+1)/2,k)&subvoxmask(isub,(js+1)%2,ksub))>0)) {
                  MRIvox(mri_dst,i,j,k)=MRIvox(mri_dst,i,j,k)|subvoxmask(isub,jsub,ksub);
                  newfilled++;
            }
        }
    }
printf("    (filled %d voxels)\n",newfilled);
  }

  return mri_dst;
}

MRI *MRISpartialribbon(MRI_SURFACE *inner_mris,MRI_SURFACE *outer_mris,MRI *mri_src,MRI *mri_dst,MRI *mri_mask,int hemisphere)
{
  MRI *mri_inter1,*mri_inter2;

  /* Allocate new MRI structures as needed */
  mri_inter1=MRIalloc(mri_src->width,mri_src->height,mri_src->depth,mri_src->type);
  MRIcopyHeader(mri_src, mri_inter1);
  mri_inter2=MRIalloc(mri_src->width,mri_src->height,mri_src->depth,mri_src->type);
  MRIcopyHeader(mri_src, mri_inter2);
  if (!mri_dst) {
    mri_dst=MRIalloc(mri_src->width,mri_src->height,mri_src->depth,mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  printf("Creating partial volume inside outer shell...\n");
  /* Create volume inside outer shell */
  /* Create shell corresponding to surface in MRI volume (includes outer shell in surface) */
  printf("  - creating partial shell...\n");
  MRISpartialshell(mri_src,outer_mris,mri_inter1); /* partial shell in mri_inter1 */
  printf("  - flooding outside shell...\n");
  MRISpartialfloodoutside(mri_inter1,mri_inter2); /* flooded outside shell in mri_inter2 */
  printf("  - inverting flooded region...\n");
  MRIbitwisenot(mri_inter2,mri_inter2); /* flooded inside shell and shell in mri_inter2 */

  printf("Creating partial volume outside inner shell...\n");
  /* Create volume outside inner shell */
  printf("  - creating partial shell...\n");
  MRISpartialshell(mri_src,inner_mris,mri_inter1);
  printf("  - flooding outside shell...\n");
  MRISpartialfloodoutside(mri_inter1,mri_dst);
  printf("  - finding union of shell and filled volume...\n");
  MRIbitwiseor(mri_inter1,mri_dst,mri_dst);

  printf("Finding intersection of volumes...\n");
  /* Find bitwise and of volumes to create ribbon */
  MRIbitwiseand(mri_inter2,mri_dst,mri_dst);

  /* Clear up the partial edges. */
  printf("Converting volume to original resolution...\n");
  MRImajority(mri_dst,mri_dst);

  /* If masked, change to CMA labels, add white matter label and apply to mask. */
  if (mri_mask) {
    /* Create white matter volume in mri_inter1, including shell */
    printf("Creating full volume outside inner shell...\n");
    MRISshell(mri_src,inner_mris,mri_inter1);
    MRISfloodoutside(mri_inter1,mri_inter2);

    printf("Inverting volume...\n");
    MRISaccentuate(mri_inter2,mri_inter2,1,254);
    /* White matter volume includes inner surface (note that complement of 0 is 1, therefore accentuate). */
    MRIcomplement(mri_inter2,mri_inter2);
    MRISaccentuate(mri_inter2,mri_inter2,1,254);

    /* mri_dst contains cerebral cortex, mri_inter2 contains white matter and some of the gray inner shell */
    printf("Merging labelled volumes...\n");
    MRIcopy(mri_dst,mri_inter1);
    MRImergecortexwhitecma(mri_inter1,mri_inter2,mri_mask,mri_dst,hemisphere);
  }

  MRIfree(&mri_inter1);
  MRIfree(&mri_inter2);

  return mri_dst;
}

MRI *MRImergecortexwhitecma(MRI *mri_cortex,MRI *mri_white,MRI *mri_cma,MRI *mri_dst,int hemisphere)
{
  /* mri_cortex = cerebral cortex is labelled as 255, all else is 0.
     mri_white  = white matter and some of the inner gray matter shell is labelled as 255, all else is 0.
     mri_cma    = CMA labelled volume.

     This function takes the mri_cma volume and replaces all instances of:
       Left_Cerebral_Cortex;
       Left_Cerebral_White_Matter;
       Right_Cerebral_Cortex;
       Right_Cerebral_White_Matter;
       Uknown,
     with:
       Left/Right_Cerebral_Cortex if labelled in mri_cortex;
       Left/Right_Cerebral_White_Matter if labelled in mri_white (and not also labelled in mri_cortex),
       where left/right hemisphere is given by hemisphere;
       Unknown otherwise.
  */

  int width,height,depth,i,j,k,vox;

  width=mri_cma->width;
  height=mri_cma->height;
  depth=mri_cma->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width,height,depth,mri_cma->type);
    MRIcopyHeader(mri_cma,mri_dst);
  }

  for (k=0; k<depth; k++)
    for (j=0; j<height; j++)
      for (i=0; i<width; i++) {
        vox=MRIvox(mri_cma,i,j,k);
        MRIvox(mri_dst,i,j,k)=vox;
        if (hemisphere==MRI_LEFT_HEMISPHERE) {
          if ((vox==Left_Cerebral_Cortex)||(vox==Left_Cerebral_White_Matter)||(vox==Unknown)) {
            if (MRIvox(mri_cortex,i,j,k)==255)
              MRIvox(mri_dst,i,j,k)=Left_Cerebral_Cortex;
            else {
              if (MRIvox(mri_white,i,j,k)==255)
                MRIvox(mri_dst,i,j,k)=Left_Cerebral_White_Matter;
              else
                MRIvox(mri_dst,i,j,k)=Unknown;
            }
          }
        }
        if (hemisphere==MRI_RIGHT_HEMISPHERE) {
          if ((vox==Right_Cerebral_Cortex)||(vox==Right_Cerebral_White_Matter)||(vox==Unknown)) {
            if (MRIvox(mri_cortex,i,j,k)==255)
              MRIvox(mri_dst,i,j,k)=Right_Cerebral_Cortex;
            else {
              if (MRIvox(mri_white,i,j,k)==255)
                MRIvox(mri_dst,i,j,k)=Right_Cerebral_White_Matter;
              else
                MRIvox(mri_dst,i,j,k)=Unknown;
            }
          }
        }
      }

  return mri_dst;
}

