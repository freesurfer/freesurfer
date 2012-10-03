/**
 * @file  mrisutils.c
 * @brief more surface processing utils
 *
 */
/*
 * Original Authors: Segonne and Greve 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/10/03 21:28:35 $
 *    $Revision: 1.42 $
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
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

#include "fio.h"
#include "const.h"
#include "diag.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "mri.h"
#include "mrisurf.h"
#include "matrix.h"
#include "stats.h"
#include "timer.h"
#include "const.h"
#include "mrishash.h"
#include "icosahedron.h"
#include "tritri.h"
#include "timer.h"
#include "chklc.h"
#include "mrisutils.h"
#include "cma.h"
#include "gca.h"
#include "sig.h"
#include "annotation.h"
#include "mrisegment.h"

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
////                    USEFUL ROUTINES               ////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

//smooth a surface 'niter' times with a step (should be around 0.5)
void MRISsmoothSurface(MRI_SURFACE *mris,int niter,float step)
{
  VERTEX *v;
  int iter,k,m,n;
  float x,y,z;

  if (step>1)
    step=1.0f;

  for (iter=0;iter<niter;iter++)
  {
    MRIScomputeMetricProperties(mris) ;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
    }

    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      n=0;
      x = y = z = 0;
      for (m=0;m<v->vnum;m++)
      {
        x += mris->vertices[v->v[m]].tx;
        y += mris->vertices[v->v[m]].ty;
        z += mris->vertices[v->v[m]].tz;
        n++;
      }
      x/=n;
      y/=n;
      z/=n;

      v->x=v->x+step*(x-v->x);
      v->y=v->y+step*(y-v->y);
      v->z=v->z+step*(z-v->z);
    }
  }
}



void MRIScenterCOG(MRIS* mris)
{
  MRIScenterCOG2(mris,NULL,NULL,NULL);
}

// translate the COG of a surface to (0,0,0)
void MRIScenterCOG2(MRI_SURFACE *mris,double *xCOG,double *yCOG,double *zCOG)
{
  int k;
  double x,y,z;
  x=0;
  y=0;
  z=0;
  for (k=0;k<mris->nvertices;k++)
  {
    x+=mris->vertices[k].x;
    y+=mris->vertices[k].y;
    z+=mris->vertices[k].z;
  }
  x/=mris->nvertices;
  y/=mris->nvertices;
  z/=mris->nvertices;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].x-=x;
    mris->vertices[k].y-=y;
    mris->vertices[k].z-=z;
  }
  if (xCOG && yCOG && zCOG)
  {
    (*xCOG)=x;
    (*yCOG)=y;
    (*zCOG)=z;
  }
  /*       fprintf(stderr,"\nCOG Centered at x=%f y=%f z=%f",
     (float)x,(float)y,(float)z);*/
}

//peel a volume from a surface (World coordinates)
// *** type == 0:  keep the inside
//     type == 1:  peel the outside surface and set the inside value to 'val'
//     type == 2: keep the outside
//     type == 3: peel the inside surface and set the outside value to 'val'
// *** NbVoxels corresponds to:
//               - the number of kept voxels if (type>=0)
//               - the number of removed voxels if (type<0)
MRI* MRISpeelVolume(MRIS *mris,
                    MRI *mri_src,
                    MRI *mri_dst,
                    int type,
                    unsigned char val,
                    unsigned long *NbVoxels)
{
  int i,j,k,imnr;
  float x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax,u,v;
  float px,py,pz,px0,py0,pz0,px1,py1,pz1;
  int numu,numv,totalfilled,newfilled;
  double tx,ty,tz;
  unsigned char tmpval;
  unsigned long size;
  int width, height,depth;
  MRI *mri_buff;

  width=mri_src->width;
  height=mri_src->height;
  depth=mri_src->depth;

  if (mri_dst==NULL)
    mri_dst=MRIalloc(width,height,depth,mri_src->type);


  mri_buff= MRIalloc(width, height, depth, MRI_UCHAR) ;


  for (k=0;k<mris->nfaces;k++)
  {
    x0 =mris->vertices[mris->faces[k].v[0]].x;
    y0 =mris->vertices[mris->faces[k].v[0]].y;
    z0 =mris->vertices[mris->faces[k].v[0]].z;
    x1 =mris->vertices[mris->faces[k].v[1]].x;
    y1 =mris->vertices[mris->faces[k].v[1]].y;
    z1 =mris->vertices[mris->faces[k].v[1]].z;
    x2 =mris->vertices[mris->faces[k].v[2]].x;
    y2 =mris->vertices[mris->faces[k].v[2]].y;
    z2 =mris->vertices[mris->faces[k].v[2]].z;
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
    numu = ceil(2*d0);
    numv = ceil(2*dmax);


    for (v=0;v<=numv;v++)
    {
      px0 = x0 + (x2-x0)*v/numv;
      py0 = y0 + (y2-y0)*v/numv;
      pz0 = z0 + (z2-z0)*v/numv;
      px1 = x1 + (x2-x1)*v/numv;
      py1 = y1 + (y2-y1)*v/numv;
      pz1 = z1 + (z2-z1)*v/numv;
      for (u=0;u<=numu;u++)
      {
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;

        // MRIworldToVoxel(mri_src,px,py,pz,&tx,&ty,&tz);
        MRIsurfaceRASToVoxel(mri_src,px,py,pz,&tx,&ty,&tz);
        imnr=(int)(tz+0.5);
        j=(int)(ty+0.5);
        i=(int)(tx+0.5);


        if (i>=0 && i<width && j>=0 && j<height && imnr>=0 && imnr<depth)
          MRIvox(mri_buff,i,j,imnr) = 255;

      }
    }
  }

  MRIvox(mri_buff,1,1,1)= 64;
  totalfilled = newfilled = 1;
  while (newfilled>0)
  {
    newfilled = 0;
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++)
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,k-1)==64||MRIvox(mri_buff,i,j-1,k)==64||
                MRIvox(mri_buff,i-1,j,k)==64)
            {
              MRIvox(mri_buff,i,j,k)= 64;
              newfilled++;
            }
    for (k=depth-2;k>=1;k--)
      for (j=height-2;j>=1;j--)
        for (i=width-2;i>=1;i--)
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,k+1)==64||MRIvox(mri_buff,i,j+1,k)==64||
                MRIvox(mri_buff,i+1,j,k)==64)
            {
              MRIvox(mri_buff,i,j,k) = 64;
              newfilled++;
            }
    totalfilled += newfilled;
  }

  size=0;
  switch (type)
  {
  case 0:
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++)
        {
          if (MRIvox(mri_buff,i,j,k)==64)
            MRIvox(mri_dst,i,j,k) = 0 ;
          else
          {
            tmpval=MRIvox(mri_src,i,j,k);
            MRIvox(mri_dst,i,j,k) = tmpval;
            size++;
          }
        }
    break;
  case 1:
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++)
        {
          if (MRIvox(mri_buff,i,j,k)==64)
            MRIvox(mri_dst,i,j,k) = 0 ;
          else
          {
            MRIvox(mri_dst,i,j,k) = val;
            size++;
          }
        }
    break;
  case 2:
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++)
        {
          if (MRIvox(mri_buff,i,j,k)==64)
          {
            tmpval=MRIvox(mri_src,i,j,k);
            MRIvox(mri_dst,i,j,k) =  tmpval ;
          }
          else
          {
            MRIvox(mri_dst,i,j,k) = 0;
            size++;
          }
        }
    break;
  case 3:
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++)
        {
          if (MRIvox(mri_buff,i,j,k)==64)
            MRIvox(mri_dst,i,j,k) =  val;
          else
          {
            MRIvox(mri_dst,i,j,k) = 0;
            size++;
          }
        }
    break;
  }
  if (NbVoxels)
    (*NbVoxels)=size;

  MRIfree(&mri_buff);
  return mri_dst;
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
////   ROUTINES FOR MATCHING A SURFACE TO A VOLUME LABEL //////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


static int mrisClearMomentum(MRI_SURFACE *mris);
static int mrisClearGradient(MRI_SURFACE *mris);
//static int mrisClearExtraGradient(MRI_SURFACE *mris);
static MRI* mriIsolateLabel(MRI* mri_seg,int label,MRI_REGION* bbox);
static int mrisAverageSignedGradients(MRI_SURFACE *mris, int num_avgs);
static double mrisRmsValError(MRI_SURFACE *mris, MRI *mri);
static void mrisSetVal(MRIS *mris,float val);
static double mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, 
                                          float momentum,
                                          float delta_t, 
                                          MHT *mht, 
                                          float max_mag);
static int mrisLimitGradientDistance(MRI_SURFACE *mris, MHT *mht, int vno);
static int mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno);
static int mrisRemoveNormalGradientComponent(MRI_SURFACE *mris, int vno);
static int mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv);
static int mrisComputeTangentPlanes(MRI_SURFACE *mris);
static int mrisComputeIntensityTerm(MRI_SURFACE *mris, 
                                    double l_intensity, 
                                    MRI *mri,
                                    double sigma);
static int mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring);
static int mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring);

static int
mrisClearMomentum(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->odx = 0 ;
    v->ody = 0 ;
    v->odz = 0 ;
  }
  return(NO_ERROR) ;
}

static int
mrisClearGradient(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->dx = 0 ;
    v->dy = 0 ;
    v->dz = 0 ;
  }
  return(NO_ERROR) ;
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
mrisClearExtraGradient(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (mris->dx2)
      mris->dx2[vno] = mris->dy2[vno] = mris->dz2[vno] = 0 ;
  }
  return(NO_ERROR) ;
}
#endif


static int
mrisAverageSignedGradients(MRI_SURFACE *mris, int num_avgs)
{
  int    i, vno, vnb, *pnb, vnum ;
  float  dx, dy, dz, num, sigma, dot ;
  VERTEX *v, *vn ;
  MRI_SP *mrisp, *mrisp_blur ;

  if (num_avgs <= 0)
    return(NO_ERROR) ;

  if (Gdiag_no >= 0)
  {
    v = &mris->vertices[Gdiag_no] ;
    fprintf(stdout, "before averaging dot = %2.2f ",
            v->dx*v->nx+v->dy*v->ny+v->dz*v->nz) ;
  }
  if (0 && mris->status == MRIS_PARAMETERIZED_SPHERE)  /* use convolution */
  {
    sigma = sqrt((float)num_avgs) / 4.0 ;
    mrisp = MRISgradientToParameterization(mris, NULL, 1.0) ;
    mrisp_blur = MRISPblur(mrisp, NULL, sigma, -1) ;
    MRISgradientFromParameterization(mrisp_blur, mris) ;
    MRISPfree(&mrisp) ;
    MRISPfree(&mrisp_blur) ;
  }
  else for (i = 0 ; i < num_avgs ; i++)
    {
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        dx = v->dx ;
        dy = v->dy ;
        dz = v->dz ;
        pnb = v->v ;
        /*      vnum = v->v2num ? v->v2num : v->vnum ;*/
        vnum = v->vnum ;
        for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
        {
          vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
          if (vn->ripflag)
            continue ;
          dot = vn->dx * v->dx + vn->dy * v->dy + vn->dz*v->dz ;
          if (dot < 0)
            continue ;  /* pointing in opposite directions */

          num++ ;
          dx += vn->dx ;
          dy += vn->dy ;
          dz += vn->dz ;
#if 0
          if (vno == Gdiag_no)
          {
            float dot ;
            dot = vn->dx*v->dx + vn->dy*v->dy + vn->dz*v->dz ;
            if (dot < 0)
              fprintf(stdout, 
                      "vn %d: dot = %2.3f, dx = (%2.3f, %2.3f, %2.3f)\n",
                      v->v[vnb], dot, vn->dx, vn->dy, vn->dz) ;
          }
#endif
        }
        num++ ;
        v->tdx = dx / num ;
        v->tdy = dy / num ;
        v->tdz = dz / num ;
      }
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        v->dx = v->tdx ;
        v->dy = v->tdy ;
        v->dz = v->tdz ;
      }
    }
  if (Gdiag_no >= 0)
  {
    float dot ;
    v = &mris->vertices[Gdiag_no] ;
    dot = v->nx*v->dx + v->ny*v->dy + v->nz*v->dz ;
    fprintf(stdout, " after dot = %2.2f\n",dot) ;
    if (fabs(dot) > 50)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}


MRI_REGION* MRIlocateRegion(MRI *mri,int label)
{
  int i,j,k;
  int xmin,xmax,ymin,ymax,zmin,zmax;
  MRI_REGION *mri_region=(MRI_REGION*)malloc(sizeof(MRI_REGION));

  zmax=ymax=xmax=0;
  zmin=ymin=xmin=10000;

  for (k=0;k<mri->depth;k++)
    for (j=0;j<mri->height;j++)
      for (i=0;i<mri->width;i++)
        if (MRIvox(mri,i,j,k)==label)
        {
          if (k<zmin)
            zmin=k;
          if (j<ymin)
            ymin=j;
          if (i<xmin)
            xmin=i;
          if (k>zmax)
            zmax=k;
          if (j>ymax)
            ymax=j;
          if (i>xmax)
            xmax=i;
        }

  mri_region->x=xmin;
  mri_region->y=ymin;
  mri_region->z=zmin;
  mri_region->dx=xmax-xmin;
  mri_region->dy=ymax-ymin;
  mri_region->dz=zmax-zmin;

  return mri_region;
}

static MRI* mriIsolateLabel(MRI* mri_seg,int label,MRI_REGION* bbox)
{
  int i,j,k;
  int xplusdx,yplusdy,zplusdz;
  MRI *mri=MRIalloc(mri_seg->width,
                    mri_seg->height,
                    mri_seg->depth,
                    mri_seg->type);

  xplusdx=bbox->x+bbox->dx+1;
  yplusdy=bbox->y+bbox->dy+1;
  zplusdz=bbox->z+bbox->dz+1;

  for (k=bbox->z;k<zplusdz;k++)
    for (j=bbox->y;j<yplusdy;j++)
      for (i=bbox->x;i<xplusdx;i++)
        if (MRIvox(mri_seg,i,j,k)==label)
          MRIvox(mri,i,j,k)=1;

  return mri;
}

static double
mrisRmsValError(MRI_SURFACE *mris, MRI *mri)
{
  int     vno, n, xv, yv, zv ;
  double    val, total, delta, x, y, z ;
  VERTEX  *v ;

  for (total = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    n++ ;
    MRISvertexToVoxel(mris,v, mri, &x, &y, &z) ;
    xv = nint(x) ;
    yv = nint(y) ;
    zv = nint(z) ;
    MRIsampleVolume(mri, x, y, z, &val) ;
    delta = (val - v->val) ;
    total += delta*delta ;
  }
  return(sqrt(total / (double)n)) ;
}

static void mrisSetVal(MRIS *mris,float val)
{
  int n;
  for (n=0;n<mris->nvertices;n++)
    mris->vertices[n].val=val;
}

#define MIN_NBR_DIST  (0.25)


static int
mrisRemoveNormalGradientComponent(MRI_SURFACE *mris, int vno)
{
  VERTEX   *v ;
  float    dot ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
    return(NO_ERROR) ;

  dot = v->nx*v->odx + v->ny*v->ody + v->nz*v->odz ;
  v->odx -= dot*v->nx ;
  v->ody -= dot*v->ny ;
  v->odz -= dot*v->nz ;

  return(NO_ERROR) ;
}

static int
mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno)
{
  VERTEX   *v, *vn ;
  int      n ;
  float    dx, dy, dz, dot, x, y, z, dist ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
    return(NO_ERROR) ;

  x = v->x ;
  y = v->y ;
  z = v->z ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    dx = vn->x - x ;
    dy = vn->y - y ;
    dz = vn->z - z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;

    /* too close - take out gradient component in this dir. */
    if (dist <= MIN_NBR_DIST)
    {
      dx /= dist ;
      dy /= dist ;
      dz /= dist ;
      dot = dx*v->odx + dy*v->ody + dz*v->odz ;
      if (dot > 0.0)
      {
        v->odx -= dot*dx ;
        v->ody -= dot*dy ;
        v->odz -= dot*dz ;
      }
    }
  }

  return(NO_ERROR) ;
}


static int
mrisLimitGradientDistance(MRI_SURFACE *mris, MHT *mht, int vno)
{
  VERTEX   *v ;

  v = &mris->vertices[vno] ;

  mrisRemoveNeighborGradientComponent(mris, vno) ;
  if (MHTisVectorFilled(mht, mris, vno, v->odx, v->ody, v->odz))
  {
    mrisRemoveNormalGradientComponent(mris, vno) ;
    if (MHTisVectorFilled(mht, mris, vno, v->odx, v->ody, v->odz))
    {
      v->odx = v->ody = v->odz = 0.0 ;
      return(NO_ERROR) ;
    }
  }

  return(NO_ERROR) ;
}


static double
mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, float momentum,
                            float delta_t, MHT *mht, float max_mag)
{
  static int direction = 1 ;
  double  mag ;
  int     vno, i ;
  VERTEX  *v ;

  for (i = 0 ; i < mris->nvertices ; i++)
  {
    if (direction < 0)
      vno = mris->nvertices - i - 1 ;
    else
      vno = i ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->odx = delta_t * v->dx + momentum*v->odx ;
    v->ody = delta_t * v->dy + momentum*v->ody ;
    v->odz = delta_t * v->dz + momentum*v->odz ;
    mag = sqrt(v->odx*v->odx + v->ody*v->ody + v->odz*v->odz) ;
    if (mag > max_mag) /* don't let step get too big */
    {
      mag = max_mag / mag ;
      v->odx *= mag ;
      v->ody *= mag ;
      v->odz *= mag ;
    }

    /* erase the faces this vertex is part of */

    if (mht)
      MHTremoveAllFaces(mht, mris, v) ;

    if (mht)
      mrisLimitGradientDistance(mris, mht, vno) ;

    v->x += v->odx ;
    v->y += v->ody ;
    v->z += v->odz ;

    if ((fabs(v->x) > 128.0f) ||
        (fabs(v->y) > 128.0f) ||
        (fabs(v->z) > 128.0f))
      DiagBreak() ;

    if (mht)
      MHTaddAllFaces(mht, mris, v) ;
  }

  direction *= -1 ;
  return(delta_t) ;
}

#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y,\
                                               v1->z-v0->z)

static int
mrisComputeTangentPlanes(MRI_SURFACE *mris)
{
  VECTOR  *v_n, *v_e1, *v_e2, *v ;
  int     vno ;
  VERTEX  *vertex ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    /* now find some other non-parallel vector */
#if 0
    if (!FZERO(vertex->nx) || !FZERO(vertex->ny))
    {
      VECTOR_LOAD(v, 0.0, 0.0, 1.0) ;
    }
    else
    {
      VECTOR_LOAD(v, 0.0, 1.0, 0.0) ;
    }
#else
    VECTOR_LOAD(v, vertex->ny, vertex->nz, vertex->nx) ;
#endif
    V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    if ((V3_LEN_IS_ZERO(v_e1)))  /* happened to pick a parallel vector */
    {
      VECTOR_LOAD(v, vertex->ny, -vertex->nz, vertex->nx) ;
      V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    }

    if ((V3_LEN_IS_ZERO(v_e1)) && 
        DIAG_VERBOSE_ON)  /* happened to pick a parallel vector */
      fprintf(stderr, "vertex %d: degenerate tangent plane\n", vno) ;
    V3_CROSS_PRODUCT(v_n, v_e1, v_e2) ;
    V3_NORMALIZE(v_e1, v_e1) ;
    V3_NORMALIZE(v_e2, v_e2) ;
    vertex->e1x = V3_X(v_e1) ;
    vertex->e2x = V3_X(v_e2) ;
    vertex->e1y = V3_Y(v_e1) ;
    vertex->e2y = V3_Y(v_e2) ;
    vertex->e1z = V3_Z(v_e1) ;
    vertex->e2z = V3_Z(v_e2) ;
  }

  VectorFree(&v) ;
  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Fit a 1-d quadratic to the surface locally and move the
          vertex in the normal direction to improve the fit.
------------------------------------------------------*/
static int
mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv)
{
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n ;
  VERTEX   *v, *vn ;
  float    ui, vi, rsq, a, b ;

  if (FZERO(l_curv))
    return(NO_ERROR) ;

  mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(2, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 2, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < v->vtotal ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ;
      vi = V3_DOT(v_e2, v_nbr) ;
      rsq = ui*ui + vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsq ;
      *MATRIX_RELT(m_R, n+1, 2) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    if (!m_R_inv)
    {
      MatrixFree(&m_R) ;
      VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
    b *= l_curv ;
    v->dx += b * v->nx ;
    v->dy += b * v->ny ;
    v->dz += b * v->nz ;

    if (vno == Gdiag_no)
      fprintf(stdout, "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "a=%2.2f, b=%2.1f\n",
              vno, b*v->nx, b*v->ny, b*v->nz, a, b) ;
    MatrixFree(&m_R) ;
    VectorFree(&v_Y) ;
    MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_nbr) ;
  VectorFree(&v_A) ;
  return(NO_ERROR) ;
}


static int
mrisComputeIntensityTerm(MRI_SURFACE *mris, 
                         double l_intensity, 
                         MRI *mri, 
                         double sigma)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z, nx, ny, nz, dx, dy, dz ;
  double    val0, xw, yw, zw, del, val_outside, val_inside, delI, delV;
  //int k,ktotal ;

  if (FZERO(l_intensity))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;
    y = v->y ;
    z = v->z ;

    // MRIworldToVoxel(mri, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolume(mri, xw, yw, zw, &val0) ;
    //    sigma = v->val2 ;

    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;

#if 1
    {
      double val;
      xw = x + nx ;
      yw = y + ny ;
      zw = z + nz ;
      // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
      MRIsampleVolume(mri, xw, yw, zw, &val) ;
      val_outside = val ;



      xw = x - nx ;
      yw = y - ny ;
      zw = z - nz ;
      // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
      MRIsampleVolume(mri, xw, yw, zw, &val) ;
      val_inside = val ;
    }
#else
    /* compute intensity gradient using smoothed volume */
    {
      double dist, val, step_size ;
      int  n ;

      step_size = MIN(sigma/2, 0.5) ;
      ktotal = 0.0 ;
      for (n = 0, val_outside = val_inside = 0.0, dist = step_size ;
           dist <= 2*sigma;
           dist += step_size, n++)
      {
        k = exp(-dist*dist/(2*sigma*sigma)) ;
        ktotal += k ;
        xw = x + dist*nx ;
        yw = y + dist*ny ;
        zw = z + dist*nz ;
        // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
        MRIsampleVolume(mri, xw, yw, zw, &val) ;
        val_outside += k*val ;

        xw = x - dist*nx ;
        yw = y - dist*ny ;
        zw = z - dist*nz ;
        // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
        MRIsampleVolume(mri, xw, yw, zw, &val) ;
        val_inside += k*val ;
      }
      val_inside /= (double)ktotal ;
      val_outside /= (double)ktotal ;
    }
#endif


    delV = v->val - val0 ;
    delI = (val_outside - val_inside) / 2.0 ;

    if (!FZERO(delI))
      delI /= fabs(delI) ;
    else if (delV<0) /*we are inside*/
      delI = -1 ;
    else
      delI = 1 ;   /* intensities tend to increase inwards */

    del = l_intensity * delV * delI ;

    dx = nx * del ;
    dy = ny * del ;
    dz = nz * del ;

    v->dx += dx ;
    v->dy += dy ;
    v->dz += dz ;

  }

  return(NO_ERROR) ;
}

static int
mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *vertex, *vn ;
  float   sx, sy, sz, nx, ny, nz, nc, x, y, z ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = vertex->nx ;
    ny = vertex->ny ;
    nz = vertex->nz ;
    x = vertex->x ;
    y = vertex->y ;
    z = vertex->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < vertex->vnum ; m++)
    {
      vn = &mris->vertices[vertex->v[m]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
    nc = sx*nx+sy*ny+sz*nz;   /* projection onto normal */
    sx = l_spring*nc*nx ;              /* move in normal direction */
    sy = l_spring*nc*ny ;
    sz = l_spring*nc*nz ;

    vertex->dx += sx ;
    vertex->dy += sy ;
    vertex->dz += sz ;
    //if (vno == Gdiag_no)
    //fprintf(stdout, "v %d spring normal term:  (%2.3f, %2.3f, %2.3f)\n",
    //        vno, sx, sy, sz) ;
  }


  return(NO_ERROR) ;
}

static int
mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int     vno, n, m ;
  VERTEX  *v, *vn ;
  float   sx, sy, sz, x, y, z, nc ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    if (v->border && !v->neg)
      continue ;

    x = v->x ;
    y = v->y ;
    z = v->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
#if 0
    n = 4 ;  /* avg # of nearest neighbors */
#endif
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }

    nc = sx*v->nx+sy*v->ny+sz*v->nz;   /* projection onto normal */
    sx -= l_spring*nc*v->nx ;                   /* remove  normal component */
    sy -= l_spring*nc*v->ny ;
    sz -= l_spring*nc*v->nz;

    v->dx += sx ;
    v->dy += sy ;
    v->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring tangent term: (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }


  return(NO_ERROR) ;
}


#define MOMENTUM 0.8
#define MAX_ASYNCH_MM       0.3
#define MAX_REDUCTIONS     0
#define REDUCTION_PCT      0.5
#define AVERAGES_NBR 1
#define BASE_DT_SCALE    1.0

MRIS *MRISmatchSurfaceToLabel(MRIS *mris,
                              MRI *mri_seg,
                              int label,
                              MRI_REGION *mri_region,
                              INTEGRATION_PARMS *integration_parms,
                              int connectivity)
{
  int bbox_indicator=0,done,niterations,n,nreductions=0;
  MRI_REGION *bbox;
  MRI *mri;
  INTEGRATION_PARMS *parms;
  int parms_indicator=0;
  double sse,last_sse,rms,last_rms,base_dt,dt,delta_t=0.0;
  double tol;
  MHT *mht = NULL;
  int avgs=AVERAGES_NBR;
  struct timeb then;
  int msec;


  TimerStart(&then);

  if (integration_parms==NULL)
  {
    parms=(INTEGRATION_PARMS*)calloc(1,sizeof(INTEGRATION_PARMS));
    parms_indicator=1;
    parms->projection = NO_PROJECTION ;
    parms->niterations=5;
    parms->dt=0.5f;
    parms->base_dt = BASE_DT_SCALE*parms->dt ;
    parms->tol=1e-3;
    parms->l_spring = 0.0f ;
    parms->l_curv = 1.0 ;
    parms->l_intensity = 1.0f ;
    parms->l_tspring = 1.0f ;
    parms->l_nspring = 0.2f ;
    parms->momentum = MOMENTUM;
    parms->dt_increase = 1.0 /* DT_INCREASE */;
    parms->dt_decrease = 0.50 /* DT_DECREASE*/ ;
    parms->error_ratio = 50.0 /*ERROR_RATIO */;
    /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
    parms->l_surf_repulse = 1.0 ;
    parms->l_repulse = 0.0;// ;
    parms->sigma=2.0f;
    parms->mri_brain=mri_seg; /*necessary for using MRIcomputeSSE*/
  }
  else
    parms=integration_parms;

  niterations=parms->niterations;
  base_dt=parms->base_dt;
  tol=parms->tol;

  if (mri_region)
    bbox=mri_region;
  else
  {
    bbox=MRIlocateRegion(mri_seg,label);
    bbox_indicator=1;
  }

  mri=mriIsolateLabel(mri_seg,label,bbox);

  mrisClearMomentum(mris);
  MRIScomputeMetricProperties(mris);
  MRISstoreMetricProperties(mris);
  MRIScomputeNormals(mris);

  switch (connectivity)
  {
  case 1:
    mrisSetVal(mris,0.65);
    break;
  case 2:
    mrisSetVal(mris,0.5);
    break;
  case 3:
    mrisSetVal(mris,0.75);
    break;
  case 4:
    mrisSetVal(mris,0.3);
    break;
  default:
    mrisSetVal(mris,0.5);
    break;
  }

  last_rms=rms=mrisRmsValError(mris,mri);

  mris->noscale=1;
  last_sse=sse=MRIScomputeSSE(mris,parms)/(double)mris->nvertices;


  if (1)//Gdiag & DIAG_SHOW)
    fprintf(stdout,"%3.3d: dt: %2.4f, sse=%2.4f, rms=%2.4f\n",
            0,0.0f,(float)sse,(float)rms);

  for (n=0;n<niterations;n++)
  {
    dt=base_dt;
    nreductions=0;

    mht=MHTfillTable(mris,mht);

    mrisClearGradient(mris);
    MRISstoreMetricProperties(mris);
    MRIScomputeNormals(mris);

    /*intensity based terms*/
    mrisComputeIntensityTerm(mris,parms->l_intensity,mri, parms->sigma);

    mrisAverageSignedGradients(mris,avgs);

    /*smoothness terms*/
    mrisComputeQuadraticCurvatureTerm(mris,parms->l_curv);
    mrisComputeNormalSpringTerm(mris,parms->l_nspring);
    mrisComputeTangentialSpringTerm(mris,parms->l_tspring);

    //      sprintf(fname,"./zurf%d",n);
    //MRISwrite(mris,fname);

    do
    {
      MRISsaveVertexPositions(mris,TMP_VERTICES);
      delta_t=mrisAsynchronousTimeStepNew(mris,
                                          parms->momentum,
                                          dt,
                                          mht,
                                          MAX_ASYNCH_MM);

      MRIScomputeMetricProperties(mris);
      rms=mrisRmsValError(mris,mri);
      sse=MRIScomputeSSE(mris,parms)/(double)mris->nvertices;

      done=1;
      if (rms>last_rms-tol) /*error increased - reduce step size*/
      {
        nreductions++;
        dt*=REDUCTION_PCT;
        if (0)//Gdiag & DIAG_SHOW)
          fprintf(stdout,"      sse=%2.1f, last_sse=%2.1f,\n"
                  "      ->  time setp reduction %d of %d to %2.3f...\n",
                  sse,last_sse,nreductions,MAX_REDUCTIONS+1,dt);

        mrisClearMomentum(mris);
        if (rms>last_rms) /*error increased - reject step*/
        {
          MRISrestoreVertexPositions(mris,TMP_VERTICES);
          MRIScomputeMetricProperties(mris);
          done=(nreductions>MAX_REDUCTIONS);
        }
      }
    }
    while (!done);
    last_sse=sse;
    last_rms=rms;

    rms=mrisRmsValError(mris,mri);
    sse=MRIScomputeSSE(mris,parms)/(double)mris->nvertices;
    fprintf(stdout,"%3d, sse=%2.1f, last sse=%2.1f,\n"
            "     rms=%2.4f, last_rms=%2.4f\n",
            n,sse,last_sse,rms,last_rms);
  }

  msec=TimerStop(&then);
  if (1)//Gdiag & DIAG_SHOW)
    fprintf(stdout,"positioning took %2.2f minutes\n",
            (float)msec/(60*1000.0f));

  if (bbox_indicator)
    free(bbox);
  if (parms_indicator)
    free(parms);
  MRIfree(&mri);
  MHTfree(&mht);

  return mris;
}



//smooth a surface 'niter' times with a step (should be around 0.5)
void MRISsmoothSurface2(MRI_SURFACE *mris,int niter,float step,int avrg)
{
  VERTEX *v;
  int iter,k,m,n;
  float x,y,z;

  if (step>1)
    step=1.0f;

  for (iter=0;iter<niter;iter++)
  {
    MRIScomputeMetricProperties(mris) ;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
    }

    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      n=0;
      x = y = z = 0;
      for (m=0;m<v->vnum;m++)
      {
        x += mris->vertices[v->v[m]].tx;
        y += mris->vertices[v->v[m]].ty;
        z += mris->vertices[v->v[m]].tz;
        n++;
      }
      x/=n;
      y/=n;
      z/=n;

      v->dx=step*(x-v->x);
      v->dy=step*(y-v->y);
      v->dz=step*(z-v->z);

    }
    mrisAverageSignedGradients(mris,avrg);
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->x+=v->dx;
      v->y+=v->dy;
      v->z+=v->dz;
    }
  }
}


/*--------------------------------------------------------------------*/
MRIS *MRISloadSurfSubject(char *subj, char *hemi, char *surfid,
                          char *SUBJECTS_DIR)
{
  MRIS *Surf;
  char fname[2000];

  if (SUBJECTS_DIR == NULL)
  {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL)
    {
      printf("ERROR: SUBJECTS_DIR not defined in environment.\n");
      return(NULL);
    }
  }

  sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subj,hemi,surfid);
  printf("  INFO: loading surface  %s\n",fname);
  fflush(stdout);
  Surf = MRISread(fname) ;
  if (Surf == NULL)
  {
    printf("ERROR: could not load registration surface\n");
    exit(1);
  }
  printf("nvertices = %d\n",Surf->nvertices);
  fflush(stdout);

  return(Surf);
}


/*-----------------------------------------------------------------
  MRISfdr2vwth() - computes the voxel-wise (or vertex-wise) threshold
    needed to realize the given False Discovery Rate (FDR) based on the
    values in the val field. The val field is copied to the val2 field,
    and then the val2 field is thresholded.

  fdr - false dicovery rate, between 0 and 1, eg: .05
  signid -
      0 = use all values regardless of sign
     +1 = use only positive values
     -1 = use only negative values
     If a vertex does not meet the sign criteria, its val2 is 0
  log10flag - interpret val field as -log10(p)
  maskflag - use the undefval field as a mask. If the undefval of
     a vertex is 1, then its val will be used to compute the threshold
     (if the val also meets the sign criteria). If the undefval is
     0, then val2 will be set to 0.
  vwth - voxel-wise threshold between 0 and 1. If log10flag is set,
     then vwth = -log10(vwth). Vertices with p values
     GREATER than vwth have val2=0. Note that this is the same
     as requiring -log10(p) > -log10(vwth).

  So, for the val2 to be set to something non-zero, the val must
  meet the sign, mask, and threshold criteria. If val meets all
  the criteria, then val2=val (ie, no log10 transforms). The val
  field itself is not altered.

  Return Values:
    0 - evertying is OK
    1 - no vertices met the mask and sign criteria

  Ref: http://www.sph.umich.edu/~nichols/FDR/FDR.m
  Thresholding of Statistical Maps in Functional Neuroimaging Using
  the False Discovery Rate.  Christopher R. Genovese, Nicole A. Lazar,
  Thomas E. Nichols (2002).  NeuroImage 15:870-878.

  See also: fdr2vwth() in sig.c
  ---------------------------------------------------------------*/
int MRISfdr2vwth(MRIS *surf, double fdr, int signid,
                 int log10flag, int maskflag, double *vwth)
{
  double *p=NULL, val, val2null;
  int vtxno, np;

  if (log10flag) val2null = 0;
  else          val2null = 1;

  p = (double *) calloc(surf->nvertices,sizeof(double));
  np = 0;
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++)
  {
    if ((maskflag && !surf->vertices[vtxno].undefval) ||
        surf->vertices[vtxno].ripflag) continue;
    val = surf->vertices[vtxno].val;
    if (signid == -1 && val > 0) continue;
    if (signid == +1 && val < 0) continue;
    val = fabs(val);
    if (log10flag) val = pow(10,-val);
    p[np] = val;
    np++;
  }

  // Check that something met the match criteria,
  // otherwise return an error
  if (np==0)
  {
    printf("WARNING: MRISfdr2vwth(): no vertices met threshold\n");
    free(p);
    return(1);
  }

  *vwth = fdr2vwth(p,np,fdr);

  // Perform the thresholding
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++)
  {
    if ((maskflag && !surf->vertices[vtxno].undefval) ||
        surf->vertices[vtxno].ripflag)
    {
      // Set to null if masking and not in the mask
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }
    val = surf->vertices[vtxno].val;
    if (signid == -1 && val > 0)
    {
      // Set to null if wrong sign
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }
    if (signid == +1 && val < 0)
    {
      // Set to null if wrong sign
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }

    val = fabs(val);
    if (log10flag) val = pow(10,-val);

    if (val > *vwth)
    {
      // Set to null if greather than thresh
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }

    // Otherwise, it meets all criteria, so
    // pass the original value through
    surf->vertices[vtxno].val2 = surf->vertices[vtxno].val;
  }

  // Change the vwth to log10 if needed
  if (log10flag) *vwth = -log10(*vwth);
  free(p);

  printf("MRISfdr2vwth(): np = %d, nv = %d, fdr = %g, vwth=%g\n",
	 np,surf->nvertices,fdr,*vwth);

  return(0);
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*---------------------------------------------------------------*/
int MRISfwhm2niters(double fwhm, MRIS *surf)
{
  double avgvtxarea, gstd;
  int niters;

  MRIScomputeMetricProperties(surf);
  avgvtxarea = surf->total_area/surf->nvertices;

  if (surf->group_avg_surface_area > 0)  {
    // This should be ok even if metric properties have been scaled ??
    // Always do this now (4/9/10)
    avgvtxarea *= (surf->group_avg_surface_area/surf->total_area);
    //if(getenv("FIX_VERTEX_AREA") != NULL)    {
    //printf("INFO: fwhm2niters: Fixing group surface area\n");
    // avgvtxarea *= (surf->group_avg_surface_area/surf->total_area);
    //}
    //else printf("INFO: fwhm2niters: NOT fixing group surface area\n");
  }

  gstd = fwhm/sqrt(log(256.0));
  //1.14 is a fudge factor based on empirical fit of nearest neighbor
  niters = floor(1.14*(4*PI*(gstd*gstd))/(7*avgvtxarea) + 0.5);
  return(niters);
}


/*----------------------------------------------------------------------*/
double MRISniters2fwhm(int niters, MRIS *surf)
{
  double avgvtxarea, gstd, fwhm;

  MRIScomputeMetricProperties(surf);
  avgvtxarea = surf->total_area/surf->nvertices;

  if (surf->group_avg_surface_area > 0)  {
    // This should be ok even if metric properties have been scaled ??
    // Always do this now (4/9/10)
    avgvtxarea *= (surf->group_avg_surface_area/surf->total_area);
    //if (getenv("FIX_VERTEX_AREA") != NULL)    {
    //printf("INFO: niters2fwhm: Fixing group surface area\n");
    // avgvtxarea *= (surf->group_avg_surface_area/surf->total_area);
    //}
    //else printf("INFO: fwhm2niters: NOT fixing group surface area\n");
  }

  gstd = sqrt(7*avgvtxarea*niters/(1.14*4*PI));
  fwhm = gstd*sqrt(log(256.0));
  return(fwhm);
}


/*---------------------------------------------------------------*/
int MRISfwhm2nitersSubj(double fwhm, char *subject, char *hemi, char *surfname)
{
  char *SUBJECTS_DIR, surfpath[2000];
  MRIS *surf;
  int niters;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL)  {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    return(-1);
  }
  sprintf(surfpath,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);

  surf = MRISread(surfpath);
  if (surf == NULL)  {
    printf("ERROR: could not read %s\n",surfpath);
    return(-1);
  }
  MRIScomputeMetricProperties(surf);

  niters = MRISfwhm2niters(fwhm,surf);

  return(niters);
}


/*----------------------------------------------------------------------
  MRISfwhm() - estimates fwhm from global ar1 mean
  ----------------------------------------------------------------------*/
double MRISfwhmFromAR1(MRIS *surf, double ar1)
{
  double fwhm, gstd,InterVertexDistAvg;

  MRIScomputeMetricProperties(surf);
  InterVertexDistAvg = surf->avg_vertex_dist;
  if (surf->group_avg_surface_area > 0)  {
    // This should be ok even if metric properties have been scaled ??
    // Always do this now (4/9/10)
    InterVertexDistAvg *= 
      sqrt(surf->group_avg_surface_area/surf->total_area);
    //if (getenv("FIX_VERTEX_AREA") != NULL)    {
    //printf("INFO: fwhmFromAR1: Fixing group surface area\n");
    //InterVertexDistAvg *= 
    //  sqrt(surf->group_avg_surface_area/surf->total_area);
    //}
    //else printf("INFO: fwhmFromAR1: NOT fixing group surface area\n");
  }

  if (ar1 > 0.0)  {
    gstd = InterVertexDistAvg/sqrt(-4*log(ar1));
    fwhm = gstd*sqrt(log(256.0));
  }
  else
  {
    printf("WARNING: ar1 = %g <= 0. Setting fwhm to 0.\n",ar1);
    fwhm = 0.0;
  }

  return(fwhm);
}


/*----------------------------------------------------------------------
  MRIscale() - scales vertex XYZ by scale.
  ----------------------------------------------------------------------*/
int MRISscale(MRIS *mris, double scale)
{
  int    vno ;
  VERTEX *v ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->x *= scale;
    v->y *= scale;
    v->z *= scale;
  }
  return(0);
}


/*-------------------------------------------------------------------------
  MRISsmoothingArea() - computes the area coverted by the give number
  of iterations at the given vertex. Essentially it's the area of the
  neighborhood.
  -------------------------------------------------------------------------*/
double MRISsmoothingArea(MRIS *mris, int vtxno, int niters)
{
  MRI *mri;
  int n,nhits;
  double val,area;

  // alloc a surface mri, zeros for everybody
  mri = MRIalloc(mris->nvertices,1,1,MRI_FLOAT);

  // create a delta function at the target vertex
  MRIsetVoxVal(mri,vtxno,0,0,0,1);

  // smooth it by the number of iterations
  MRISsmoothMRI(mris, mri, niters, NULL, mri);

  // find the non-zero vertices. these are the vertices in the neighborhood
  // add up the number of vertices and area
  nhits = 0;
  area = 0.0;
  for (n=0; n < mris->nvertices; n++)
  {
    val = MRIgetVoxVal(mri,n,0,0,0);
    if (val > 0.0)
    {
      nhits ++;
      if (mris->group_avg_vtxarea_loaded)
        area += mris->vertices[n].group_avg_area;
      else
        area += mris->vertices[n].area;
    }
  }
  MRIfree(&mri);
  printf("%6d  %3d   %4d %7.1f\n",vtxno,niters,nhits,area);
  fflush(stdout);

  return(area);
}


/*-------------------------------------------------------------------
  MRISseg2annot() - converts a segmentation in a surface-encoded
  MRI struct to an annotation in an MRIS struct. If ctab is null,
  it tries to use mris->ct. See also MRISannotIndex2Seg().
  -------------------------------------------------------------------*/
int MRISseg2annot(MRIS *mris, MRI *surfseg, COLOR_TABLE *ctab)
{
  int vtxno, segid, ano;

  if (ctab == NULL)
  {
    if (mris->ct == NULL)
    {
      printf("ERROR: MRISseg2annot: both ctab and mris->ct are NULL\n");
      return(1);
    }
  }
  else mris->ct = ctab;
  set_atable_from_ctable(mris->ct);

  for (vtxno=0; vtxno < mris->nvertices; vtxno++)
  {
    segid = MRIgetVoxVal(surfseg,vtxno,0,0,0);
    ano = index_to_annotation(segid);
    mris->vertices[vtxno].annotation = ano;
    //printf("%5d %2d %2d %s\n",vtxno,segid,ano,index_to_name(segid));
  }

  return(0);
}


/*----------------------------------------------------------------
  MRISannotIndex2Seg() - creates an MRI struct where each voxel is the
  annotation index. The struct has nvertices columns, 1 row, 1 slice,
  and 1 frame. It should not be misconstrued as a volume. See also
  MRISseg2annot().
  ----------------------------------------------------------------*/
MRI *MRISannotIndex2Seg(MRIS *mris)
{
  MRI *seg;
  int vno, annot, annotid;

  annotid = 0;
  seg = MRIalloc(mris->nvertices,1,1,MRI_FLOAT);
  for (vno = 0; vno < mris->nvertices; vno ++)
  {
    annot = mris->vertices[vno].annotation;
    if (mris->ct)      CTABfindAnnotation(mris->ct, annot, &annotid);
    else        annotid = annotation_to_index(annot);
    MRIsetVoxVal(seg,vno,0,0,0,annotid);
  }
  return(seg);
}



/*!/
  \fn double MRISvolumeInSurf(MRIS *mris)
  \brief Computes the volume inside a surface. (Xiao)
*/
double MRISvolumeInSurf(MRIS *mris)
{
  int fno;
  double total_volume, face_area;
  VECTOR *v_a, *v_b, *v_n, *v_cen;
  VERTEX  *v0, *v1, *v2;
  FACE *face;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_cen = VectorAlloc(3, MATRIX_REAL) ;     /* centroid vector */

  total_volume = 0;
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;

    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;

    VERTEX_EDGE(v_a, v0, v1) ;
    VERTEX_EDGE(v_b, v0, v2) ;

    /* face normal vector */
    V3_CROSS_PRODUCT(v_a, v_b, v_n) ;
    face_area = V3_LEN(v_n) * 0.5f ;

    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */

    /* compute face centroid */
    V3_X(v_cen) = (v0->x + v1->x + v2->x)/3.0;
    V3_Y(v_cen) = (v0->y + v1->y + v2->y)/3.0;
    V3_Z(v_cen) = (v0->z + v1->z + v2->z)/3.0;

    total_volume += V3_DOT(v_cen, v_n)*face_area;
  }

  MatrixFree(&v_cen);
  MatrixFree(&v_a);
  MatrixFree(&v_b);
  MatrixFree(&v_n);

  total_volume /= 3.0;
  return(total_volume);
}

/*
  note that if the v->marked2 fields are set in vertices (e.g. from
  a .annot file), then these vertices will not be considered in
  the search for non-cortical vertices (that is, they will be labeled
  cortex).
*/
LABEL *MRIScortexLabel(MRI_SURFACE *mris, MRI *mri_aseg, int min_vertices) 
{
  LABEL      *lcortex ;
  int        vno, label, nvox, total_vox, adjacent, x, y, z, target_label,l ;
  VERTEX     *v ;
  double     xv, yv, zv, val, xs, ys, zs, d ;
  MRI_REGION box ;

  printf("generating cortex label...\n") ;

  mri_aseg = MRIcopy(mri_aseg, NULL) ; // so we can mess with it

  // remove the posterior few mm of the ventricles to prevent
  // them poking into the calcarine
#define ERASE_MM 3
  for (l = 0 ; l < 2 ; l++)
  {
    if (l == 0)
      target_label = Left_Lateral_Ventricle ;
    else
      target_label = Right_Lateral_Ventricle ; 
    MRIlabelBoundingBox(mri_aseg, target_label, &box) ;
    for (z = box.z+box.dz-(ERASE_MM+1) ; z < box.z+box.dz ; z++)
      for (y = 0 ; y < mri_aseg->height ; y++)
        for (x = 0 ; x < mri_aseg->width ; x++)
        {
          label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
          if (label == target_label)
            MRIsetVoxVal(mri_aseg, x, y, z, 0, 0) ; // erase it
        }
  }
  MRISsetMarks(mris, 1) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->marked2 > 0)  // already must be cortex
      continue ;
    if (vno == Gdiag_no )
      DiagBreak() ;

    // don't sample inside here due to thin parahippocampal wm.
    // The other interior labels
    // will already have been ripped (and hence not marked)
    for (d = 0 ; d <= 2 ; d += 0.5) {
      xs = v->x + d*v->nx ;
      ys = v->y + d*v->ny ;
      zs = v->z + d*v->nz ;
      if (mris->useRealRAS)
        MRIworldToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      else
        MRIsurfaceRASToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      label = nint(val) ;
      if (label == Left_Lateral_Ventricle ||
          label == Right_Lateral_Ventricle ||
          label == Third_Ventricle ||
          label == Left_Accumbens_area ||
          label == Right_Accumbens_area ||
          label == Left_Caudate ||
          label == Right_Caudate ||
          IS_CC(label) || 
          label == Left_Pallidum ||
          label == Right_Pallidum ||
          IS_HIPPO(label) ||
          IS_AMYGDALA(label) ||
          IS_LAT_VENT(label) ||
          label == Third_Ventricle ||
          label == Right_Thalamus_Proper ||
          label == Left_Thalamus_Proper ||
          label == Brain_Stem ||
          label == Left_VentralDC ||
          label == Right_VentralDC)
      {
        if (label == Left_Putamen || label == Right_Putamen)
          DiagBreak() ;
        if (vno == Gdiag_no)
          DiagBreak() ;
        v->marked = 0 ;
      }
    }
    // putamen can be adjacent to insula in aseg, but shouldn't be inferior
    /* now check for putamen superior to this point. If there's a lot
       of it there, then we are in basal forebrain and not cortex. */
    for (adjacent = total_vox = nvox = 0, d = 0 ; 
         d <= 10 ; d += 0.5, total_vox++) 
    {
      xs = v->x ; ys = v->y ;
      zs = v->z + d ;  // sample superiorly
      if (mris->useRealRAS)
        MRIworldToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      else
        MRIsurfaceRASToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      label = nint(val) ;
      if (label == Left_Putamen || label == Right_Putamen)
      {
        nvox++ ;
        if (d < 1.5)
          adjacent = 1 ;
      }
    }
    if (adjacent &&
        (double)nvox/(double)total_vox > 0.5) // more than 50% putamen
      v->marked = 0 ;
  }

  // remove small holes that shouldn't be non-cortex
  {
    LABEL **label_array ;
    int   nlabels, n, i ;

    MRISinvertMarks(mris) ; // marked->not cortex now
    MRISsegmentMarked(mris, &label_array, &nlabels, 0) ;
    printf("%d non-cortical segments detected\n", nlabels) ;
    if (min_vertices < 0)  // means only keep max segment
    {
      for (n = 0 ; n < nlabels ; n++)
        if (label_array[n]->n_points > min_vertices)
          min_vertices = label_array[n]->n_points ;
      printf("only using segment with %d vertices\n", min_vertices) ;
    }

    for (n = 0 ; n < nlabels ; n++)
    {
      if (label_array[n]->n_points < min_vertices)
      {
        printf("erasing segment %d (vno[0] = %d)\n", 
               n, label_array[n]->lv[0].vno) ;
        for (i = 0 ; i < label_array[n]->n_points ; i++)
          mris->vertices[label_array[n]->lv[i].vno].marked = 0 ; // mark it as cortex
      }
      LabelFree(&label_array[n]) ;
    }
    free(label_array) ;
  }
  
  MRISinvertMarks(mris) ;  // marked --> is cortex again
  lcortex = LabelFromMarkedSurface(mris) ;

  MRIfree(&mri_aseg) ;  // a locally edited copy, not the original
  return(lcortex) ;
}


/*-------------------------------------------------------*/
/*!
  \fn int MRISsphericalCoords(MRIS *surf)
  \brief Replaces x,y,z with theta,phi,radius. Assumes
    that the surface xyz are already on the sphere.
    Note: this is not realated to the surface-based
    spherical coords.
 */
int MRISsphericalCoords(MRIS *surf)
{
  int k;
  double x,y,z,d2,d,r,theta,phi;

  for(k=0; k < surf->nvertices; k++){
    x = surf->vertices[k].x;
    y = surf->vertices[k].y;
    z = surf->vertices[k].z;
    d2 = x*x + y*y;
    d = sqrt(d2);
    r = sqrt(d2 + z*z);
    theta = atan2(y,x);
    phi = atan2(z,d);
    surf->vertices[k].x = theta;
    surf->vertices[k].y = phi;
    surf->vertices[k].z = r;
  }
  return(0);
}


/*-------------------------------------------------------*/
/*!
  \fn int MRISripZeros(MRIS *surf, MRI *mri)
  \brief Sets ripflag=1 for vertices where the mri value is 0 
    (actually less than 1e-5). If mri is null, then uses the
    val field. No change to a vertex if ripflag already = 1.
*/
int MRISripZeros(MRIS *surf, MRI *mri)
{
  int k;
  double v;

  if(mri){
    if(mri->width != surf->nvertices){
      printf("ERROR: MRISripZeros(): dimension mismatch\n");
      return(1);
    }
  }

  for(k=0; k < surf->nvertices; k++){
    if(surf->vertices[k].ripflag) continue;
    if(mri) v = MRIgetVoxVal(mri,k,0,0,0);
    else    v = surf->vertices[k].val;
    if(fabs(v) < 1e-5) surf->vertices[k].ripflag = 1;
  }
  return(0);
}

/*!
  \fn int MRISfindPath ( int *vert_vno, int num_vno, int max_path_length,
                         int *path, int *path_length, MRIS *mris ).
  \brief Finds a path that connects all the vertices listed in vert_vno
  \params vert_vno[] - array of vertex numbers to connect
  \params num_vno - number in array of vertex numbers to connect
  \params max_path_length - max number in path
  \params path - array of connected vertex numbers
  \params path_length - pointer to length of array of connected vertex numbers
  \params mris - surface
  This was copied from tksurfer.c find_path() and modified slightly.
*/
int MRISfindPath ( int *vert_vno, int num_vno, int max_path_length,
		   int *path, int *path_length, MRIS *mris ) {
  int cur_vert_vno;
  int src_vno;
  int dest_vno;
  int vno;
  char* check;
  float* dist;
  int* pred;
  char done;
  VERTEX* v;
  VERTEX* u;
  float closest_dist;
  int closest_vno;
  int neighbor;
  int neighbor_vno;
  float dist_uv;
  int path_vno;
  int num_path = 0;
  int num_checked;
  float vu_x, vu_y, vu_z;
  int flag2d = 0; // for flattend surface?

  dist = (float*) calloc (mris->nvertices, sizeof(float));
  pred = (int*) calloc (mris->nvertices, sizeof(int));
  check = (char*) calloc (mris->nvertices, sizeof(char));
  num_path = 0;
  num_checked = 0;
  (*path_length) = 0;

  for (cur_vert_vno = 0; cur_vert_vno < num_vno-1; cur_vert_vno++) {
    /* clear everything */
    for (vno = 0; vno < mris->nvertices; vno++) {
      dist[vno] = 999999;
      pred[vno] = -1;
      check[vno] = FALSE;
    }

    /* Set src and dest */
    src_vno = vert_vno[cur_vert_vno+1];
    dest_vno = vert_vno[cur_vert_vno];

    /* make sure both are in range. */
    if (src_vno < 0 || src_vno >= mris->nvertices ||
        dest_vno < 0 || dest_vno >= mris->nvertices)
      continue;

    if (src_vno == dest_vno)
      continue;

    /* pull the src vertex in. */
    dist[src_vno] = 0;
    pred[src_vno] = vno;
    check[src_vno] = TRUE;

    done = FALSE;
    while (!done) {

      /* find the vertex with the shortest edge. */
      closest_dist = 999999;
      closest_vno = -1;
      for (vno = 0; vno < mris->nvertices; vno++)
        if (check[vno])
          if (dist[vno] < closest_dist) {
            closest_dist = dist[vno];
            closest_vno = vno;
          }
      v = &(mris->vertices[closest_vno]);
      check[closest_vno] = FALSE;

      /* if this is the dest node, we're done. */
      if (closest_vno == dest_vno) {
        done = TRUE;
      } else {
        /* relax its neighbors. */
        for (neighbor = 0; neighbor < v->vnum; neighbor++) {
          neighbor_vno = v->v[neighbor];
          u = &(mris->vertices[neighbor_vno]);

          /* calc the vector from u to v. */
          vu_x = u->x - v->x;
          vu_y = u->y - v->y;
	  if (flag2d)	    vu_z = 0;
	  else     	    vu_z = u->z - v->z;

          /* recalc the weight. */
	  if (flag2d)
	    dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
			   ((v->y - u->y) * (v->y - u->y)));
	  else
	    dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
			   ((v->y - u->y) * (v->y - u->y)) +
			   ((v->z - u->z) * (v->z - u->z)));

          /* if this is a new shortest path, update the predecessor,
             weight, and add it to the list of ones to check next. */
          if (dist_uv + dist[closest_vno] < dist[neighbor_vno]) {
            pred[neighbor_vno] = closest_vno;
            dist[neighbor_vno] = dist_uv + dist[closest_vno];
            check[neighbor_vno] = TRUE;
          }
        }
      }
      num_checked++;
      if ((num_checked % 100) == 0) {
        printf (".");
        fflush (stdout);
      }
    }

    /* add the predecessors from the dest to the src to the path. */
    path_vno = dest_vno;
    path[(*path_length)++] = dest_vno;
    while (pred[path_vno] != src_vno &&
           (*path_length) < max_path_length ) {
      path[(*path_length)++] = pred[path_vno];
      path_vno = pred[path_vno];
    }
  }
  printf (" done\n");
  fflush (stdout);

  free (dist);
  free (pred);
  free (check);

  return (ERROR_NONE);
}

MRI *MRIScomputeFlattenedVolume(MRI_SURFACE *mris,
                                MRI *mri,
                                double res,
                                int nsamples,
                                int normalize,
                                MRI **pmri_vertices,
                                int smooth_iters,
                                double wm_dist,
                                double outside_dist)
{
  MRI    *mri_flat, *mri_mask, *mri_counts, *mri_vno ;
  int    vno, width, height, u, v, w, fno, num ;
  int    uk, vk, ui, vi, whalf = 3, nv, wm_samples, outside_samples ;
  double xmin, xmax, ymin, ymax, fdist, x, y, z, dx, dy, dz, norm, xf, yf, val, xv, yv,zv,
         oval, max_change ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *face ;
  MHT    *mht ;

  wm_samples = nint(wm_dist / res) ;
  outside_samples = nint(outside_dist / res) ;
  mht = MHTfillTableAtResolution(mris, NULL, FLATTENED_VERTICES, 2.0) ;
  ymax = xmax = -1e10 ;
  ymin = xmin = 1e10 ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v0 = &mris->vertices[vno] ;
    if (v0->ripflag)
      continue ;
    if (v0->fx < xmin)
      xmin = v0->fx ;
    if (v0->fy < ymin)
      ymin = v0->fy ;

    if (v0->fx > xmax)
      xmax = v0->fx ;
    if (v0->fy > ymax)
      ymax = v0->fy ;
  }

  width = ceil((xmax - xmin)/res) ;
  height = ceil((ymax - ymin)/res) ;
  mri_flat = MRIalloc(width, height, nsamples, MRI_FLOAT) ;
  mri_mask = MRIalloc(width, height, nsamples, MRI_UCHAR) ;
  mri_counts = MRIalloc(width, height, nsamples, MRI_INT) ;

/*
  the first frame of mri_vno contains the # of vertices mapped to that (i,j) position, then
  the subsequent frames contain the vertex numbers
*/
  mri_vno = MRIalloc(width, height, nsamples, MRI_INT) ;
  MRIsetResolution(mri_flat, res, res, 3.0/(float)nsamples) ;
  MRIsetResolution(mri_mask, res, res, 3.0/(float)nsamples) ;
  MRIsetResolution(mri_vno, res, res, 3.0/(float)nsamples) ;
  num = 0 ;
  whalf = ceil(2.0 / res) ;
  printf("using mask window size = %d\n", 2*whalf+1) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v0 = &mris->vertices[vno] ;
    u = ceil(v0->fx - xmin) / res ;
    v = ceil(v0->fy - ymin) / res ;
    for (uk = -whalf ; uk <= whalf ; uk++)
    {
      ui = mri_mask->xi[u+uk] ;
      for (vk = -whalf ; vk <= whalf ; vk++)
      {
        vi = mri_mask->yi[v+vk] ;
        if (MRIgetVoxVal(mri_mask, ui, vi, 0, 0) == 0)
          num++ ;
        MRIsetVoxVal(mri_mask, ui, vi, 0, 0, 1) ;
      }
    }
  }
  printf("%d voxels set in mask\n", num);

  for (num = u = 0 ; u < width ; u++)
  {
    if (!(u % 100))
    {
      printf("u = %d of %d\n", u, width) ;
      fflush(stdout) ;
    }
    for (v = 0 ; v < height ; v++)
    {
      if (u == Gx && v == Gy)
        DiagBreak() ;
#if 1
      if (MRIgetVoxVal(mri_mask, u, v, 0, 0) == 0)
        continue ;
#endif
      num++ ;
      xf = u*res + xmin ;
      yf = v*res + ymin ;
      MHTfindClosestFaceGeneric(mht, mris, xf,  yf, 0, 1000, -1, -1, &face, &fno, &fdist) ;
      v0 = &mris->vertices[face->v[0]] ;
      v1 = &mris->vertices[face->v[1]] ;
      v2 = &mris->vertices[face->v[2]] ;
      if (v0->ripflag || v1->ripflag || v2->ripflag /* || fdist > 2*/)
        continue ;
      if (v0-mris->vertices == Gdiag_no ||
          v1-mris->vertices == Gdiag_no ||
          v2-mris->vertices == Gdiag_no)
        DiagBreak() ;

      dx = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->nx, v1->nx, v2->nx) ;
      dy = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->ny, v1->ny, v2->ny) ;
      dz = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->nz, v1->nz, v2->nz) ;
      norm = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
      if (FZERO(norm))
        continue ;
      dx /= norm ; dy /= norm ; dz /= norm ;  // make it unit length
      if (normalize)   // divide the ribbon at this point into equally spaced samples
      {
        norm =
          (sqrt(SQR(v0->pialx - v0->whitex) +SQR(v0->pialy - v0->whitey) +SQR(v0->pialz - v0->whitez)) +
           sqrt(SQR(v1->pialx - v1->whitex) +SQR(v1->pialy - v1->whitey) +SQR(v1->pialz - v1->whitez)) +
           sqrt(SQR(v2->pialx - v2->whitex) +SQR(v2->pialy - v2->whitey) +SQR(v2->pialz-v2->whitez)))/3;
	norm += wm_dist + outside_dist ;
        norm /= nsamples ;  // divide average thickness into this many samples
      }
      else   // use uniform spacing
        norm = 1.0 / nsamples ;

      dx *= norm; dy *= norm ; dz *= norm ;

      x = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf,0,v0->whitex, v1->whitex, v2->whitex) ;
      y = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf,0,v0->whitey, v1->whitey, v2->whitey) ;
      z = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf,0,v0->whitez, v1->whitez, v2->whitez) ;
      nv = MRIgetVoxVal(mri_vno, u, v, 0, 0) ;   // # of vertices mapped to this location
      vno = v0-mris->vertices ;
      for (w = 0 ; w < nv ; w++)
      {
	int vno2 ;
	vno2 = (int)MRIgetVoxVal(mri_vno, u, v, w+1, 0) ;
	if (vno2 == vno)  // already in the list
	  break ;
      }
      MRIsetVoxVal(mri_vno, u, v, nv+1, 0, vno) ;
      if (w == nv)
	MRIsetVoxVal(mri_vno, u, v, 0, 0, nv+1) ;  // increment # of vertices mapping here

      x -= dx * wm_dist ; y -= dy * wm_dist ; z -= dz * wm_dist ;
      for (w = 0 ; w < nsamples ; w++)
      {
        if (u == Gx && y == Gy && w == Gz)
          DiagBreak() ;
        MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
          DiagBreak() ;
        MRIsetVoxVal(mri_flat, u, v, w, 0, val) ;
        MRIsetVoxVal(mri_counts, u, v, w, 0, 1) ;
        if (v0-mris->vertices == Gdiag_no ||
            v1-mris->vertices == Gdiag_no ||
            v2-mris->vertices == Gdiag_no)
          printf("(%2.1f %2.1f %2.1f) --> (%d, %d, %d) : %2.1f\n",
                 xv, yv, zv, u, v, w, val) ;

        x += dx ; y += dy ; z += dz ;
      }
      if (v0-mris->vertices == Gdiag_no ||
          v1-mris->vertices == Gdiag_no ||
          v2-mris->vertices == Gdiag_no)
        DiagBreak() ;
    }
  }
  printf("%d voxel visited - %2.1f %% of total %d\n", num, 100.0*num/(width*height), width*height) ;
  MRIremoveNaNs(mri_vno, mri_vno) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    char fname[STRLEN] ;
    sprintf(fname, "vno.mgz") ;
    printf("writing vertex numbers to %s\n", fname) ;
    MRIwrite(mri_vno, fname) ;
  }
  // fill in holes in the flatmap where no vertices mapped by averaging with neighboring filled locations
#define MAX_ITERS 50
  for (num = 0 ; num < MAX_ITERS ; num++)
  {
    max_change = 0.0 ;
    for (u = 0 ; u < width ; u++)
    {
      for (v = 0 ; v < height ; v++)
      {
        for (w = 0 ; w < nsamples ; w++)
        {
          if (MRIgetVoxVal(mri_flat, u, v, w, 0) > 30000)
          {
            fprintf(stderr, "voxel (%d, %d, %d) = %f\n", u, v, w, MRIgetVoxVal(mri_flat,u,v,w,0)) ;
            DiagBreak() ;
          }
          if (MRIgetVoxVal(mri_counts, u, v, w, 0) > 0)  // not a hole
            continue ;
          for (val = 0.0, uk = -1 ; uk <= 1 ; uk++)
          {
            ui = mri_flat->xi[u+uk] ;
            for (vk = -1 ; vk <= 1 ; vk++)
            {
              vi = mri_flat->yi[v+vk] ;
              val += MRIgetVoxVal(mri_flat, ui, vi, w, 0) ;
            }
          }
          oval = MRIgetVoxVal(mri_flat, u, v, w, 0) ;
          val /= 9.0 ;   // # of voxels visited
          if (fabs(oval-val) > max_change)
            max_change = fabs(oval-val) ;
          MRIsetVoxVal(mri_flat, u, v, w, 0, val) ;
#if 1
          if (fabs(oval-val) < 1 && val > 50)
            MRIsetVoxVal(mri_counts, u, v, w, 0, 1) ;
#endif
        }
      }
    }
    if (max_change < 1)
      break ;
    printf("%d of %d: max change %2.1f\n", num+1, MAX_ITERS, max_change) ;
    fflush(stdout) ;
  }

  // now apply a bit of tangential smoothing
  for (num = 0 ; num < smooth_iters ; num++)
  {
    for (u = 0 ; u < width ; u++)
    {
      for (v = 0 ; v < height ; v++)
      {
        for (w = 0 ; w < nsamples ; w++)
        {
          for (val = 0.0, uk = -1 ; uk <= 1 ; uk++)
          {
            ui = mri_flat->xi[u+uk] ;
            for (vk = -1 ; vk <= 1 ; vk++)
            {
              vi = mri_flat->yi[v+vk] ;
              val += MRIgetVoxVal(mri_flat, ui, vi, w, 0) ;
            }
          }
          val /= 9.0 ;   // # of voxels visited
          MRIsetVoxVal(mri_flat, u, v, w, 0, val) ;
        }
      }
    }
  }

  MHTfree(&mht) ; MRIfree(&mri_mask) ; MRIfree(&mri_counts) ; 
  if (pmri_vertices)
    *pmri_vertices = mri_vno ;
  else
    MRIfree(&mri_vno) ;
  return(mri_flat) ;
}

/*!
  \fn MRI_SP *MRISmakeTemplate(int nsubjects, char **subjlist, 
      int nhemis, char **hemilist, char *surfregname)
  \brief Creates a surface registration template. Overlaps with 
    mris_make_template. Produces the same result as with -norot -aparc.
    Can be used on both lh and rh.
  \params number of subjects in subjlist
  \params subjlist - list of subjects
  \params nhemis - number of hemispheres in hemilist
  \params hemilist - list of hemispheres to use
  \params surfregname - name of surface registration
  Example 1:
    subjlist[0] = "s02.ghent";subjlist[1] = "s05.ghent";
    hemilist[0] = "lh";
    mrisp_template=MRISmakeTemplate(2, subjlist, 1, hemilist, "sphere.reg");
  Example 2:
    subjlist[0] = "s02.ghent";subjlist[1] = "s05.ghent";
    hemilist[0] = "rh";  hemilist[1] = "lh";
    mrisp_template=MRISmakeTemplate(2, subjlist, 2, hemilist, "sphere.left_right");
*/
MRI_SP *MRISmakeTemplate(int nsubjects, char **subjlist, int nhemis, char **hemilist, char *surfregname)
{
  static char *surface_names[] = {"inflated","smoothwm","smoothwm"} ;
  static char *curvature_names[] =  {"inflated.H","sulc",NULL} ;
  char tmpstr[2000];
  int images_per_surface = 3;
  int nsurfaces = sizeof(curvature_names) / sizeof(curvature_names[0]);
  int nparam_images = images_per_surface*nsurfaces;
  float scale = 1 ;
  char *annot_name = "aparc", *SUBJECTS_DIR, *hemi, *subject ;
  INTEGRATION_PARMS parms ;
  int which_norm = NORM_MEAN;
  int navgs = 0, nthhemi, sno, nthsubject ;
  int nbrs = 3, err ;
  MRI_SP  *mrisp, /* *mrisp_aligned,*/ *mrisp_template ;
  MRIS *mris;

  /* default template fields*/
  memset(&parms, 0, sizeof(parms)) ;
  parms.nfields=3;
  SetFieldLabel(&parms.fields[0], INFLATED_CURV_CORR_FRAME,0,0.0,0.0,0,which_norm);
  /* only use sulc for rigid registration */
  SetFieldLabel(&parms.fields[1],SULC_CORR_FRAME,1,1.0,0.0,0,which_norm);
  SetFieldLabel(&parms.fields[2],CURVATURE_CORR_FRAME,2,0.0,0.0,0,which_norm);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  mrisp_template = MRISPalloc(scale, nparam_images);
  for(nthsubject = 0; nthsubject < nsubjects; nthsubject++){
    subject = subjlist[nthsubject];

    for(nthhemi = 0; nthhemi < nhemis; nthhemi++){
      hemi = hemilist[nthhemi];
      printf("subject %s hemi %s\n",subject,hemi);
      sprintf(tmpstr, "%s/%s/surf/%s.%s",SUBJECTS_DIR, subject, hemi, surfregname);
      printf("   reading surface %s...\n",tmpstr);
      mris = MRISread(tmpstr);
      if(mris == NULL){
	printf("ERROR: could not load %s\n",tmpstr);
	return(NULL);
      }
      
      err = MRISreadAnnotation(mris, annot_name);
      if(err){
	printf("ERROR: could not load %s\n",annot_name);
	return(NULL);
      }

      MRISripMedialWall(mris) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      MRISstoreMetricProperties(mris) ;
      
      for (sno = 0; sno < nsurfaces ; sno++){
	if(curvature_names[sno]){
	  /* read in precomputed curvature file */
	  sprintf(tmpstr, "%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,curvature_names[sno]) ;
	  err = MRISreadCurvatureFile(mris, tmpstr);
	  if(err){
	    printf("ERROR: could not load %s\n",tmpstr);
	    return(NULL);
	  }
	  MRISaverageCurvatures(mris, navgs) ;
	  MRISnormalizeCurvature(mris, which_norm) ;
	}
	else {
	  sprintf(tmpstr, "%s/%s/surf/%s.%s",SUBJECTS_DIR,subject, hemi,surface_names[sno]) ;
	  err = MRISreadVertexPositions(mris, tmpstr);
	  if(err){
	    printf("ERROR: could not load %s\n",tmpstr);
	    return(NULL);
	  }
	  if (nbrs > 1) MRISsetNeighborhoodSize(mris, nbrs) ;
	  MRIScomputeMetricProperties(mris) ;
	  MRIScomputeSecondFundamentalForm(mris) ;
	  MRISuseMeanCurvature(mris) ;
	  MRISaverageCurvatures(mris, navgs) ;
	  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
	  MRISnormalizeCurvature(mris, which_norm) ;
	}
	printf("  computing parameterization for surface %s...\n",tmpstr);
	mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
	MRISPcombine(mrisp, mrisp_template, sno*3) ;
	MRISPfree(&mrisp) ;
      }
      MRISfree(&mris) ;
    }
  }
  return(mrisp_template);
}
