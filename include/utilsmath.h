/**
 * @brief utilities for distance field algorithm and the OBB Tree algorithm
 */
/*
 * Original Author: Krish Subramanium
 *
 * Geometric Tools, LLC
 * Copyright (c) 1998-2010
 * Distributed under the Boost Software License, Version 1.0.
 * http://www.boost.org/LICENSE_1_0.txt
 * http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef utilsmath_h
#define utilsmath_h

#include "mri.h"
#include "mrisurf.h"
#include "diag.h"
#include "proto.h"

#include <iostream>
#include <limits.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace Math
{
  //! A simple templated 3D Point class 
  template <typename T>
    class Point
  {
  public:      // a very simple Point class which has 3 elements 
    // where v[j] = x, y or z coord depending on j = 0, 1, 2
    T v[3]; // this is easier because VERTEX has x, y and z 
            // which doesn't lend itself to computational ease
  };
    
  //! A simple axis aligned bounding box class for a MRIS surface
  class AABB
  {
  public:   // Axis Aligned Bounding Box to trivially
    // reject points in inclusion tests
    double mincorner[3];
    double maxcorner[3];

    //! A constructor which constructs AABB using the surface verts
    //! Results in a bounding box which just bounds the surface
    //! a leeway of delta units is added to both the corners of AABB
    //! which results in a looser bounding box
    //! Essentially a tighter bounding box just touches the surface,
    //! so rejection test might 
    //! exclude voxels which are in the surface.. so use a looser AABB
    AABB(MRIS *mris, double delta)
    {
      VERTEX* v;
      mincorner[0] = mincorner[1] = mincorner[2] =
        std::numeric_limits<float>::max();
      maxcorner[0] = maxcorner[1] = maxcorner[2] =
        -std::numeric_limits<float>::max();

      for ( int i=0; i< mris->nvertices; i++)
      {
        v = &mris->vertices[i];
        if ( v->x < mincorner[0] )
          mincorner[0] = v->x ;
        if ( v->y < mincorner[1] )
          mincorner[1] = v->y ;
        if ( v->z < mincorner[2] )
          mincorner[2] = v->z ;
        
        if ( v->x > maxcorner[0] )
          maxcorner[0] = v->x ;
        if ( v->y > maxcorner[1] )
          maxcorner[1] = v->y ;
        if ( v->z > maxcorner[2] )
          maxcorner[2] = v->z ;
      }
      mincorner[0] -= delta;
      mincorner[1] -= delta;
      mincorner[2] -= delta;
      maxcorner[0] += delta;
      maxcorner[1] += delta;
      maxcorner[2] += delta;
    }

    ~AABB(){}

    //! Check to see if the point given is outside.. 
    /*! 
      \param x,y,z - the point
      \returns 1 if the point is inside and 0 otherwise
    */
    bool CheckOutside(double x, double y, double z)
    {
      if ( x < mincorner[0] || y < mincorner[1] || z < mincorner[2] || \
           x > maxcorner[0] || y > maxcorner[1] || z > maxcorner[2] )
        return 1;
      return 0;
    }
  };

  //! A simple templated bounding box class
  template <typename T>
    class BoundingBox
  {
  public:
    T minc[3]; //left bottom ( minimum corner)
    T maxc[3]; //top right ( max corner )

    //! Adjust bounds
    /*!
      \param maxdist - the maximum distance by which the 
      BB's bounds have to be extended
    */
    void adjustbounds(double maxdist)
    {
      minc[0] = minc[0] - maxdist;
      minc[1] = minc[1] - maxdist;
      minc[2] = minc[2] - maxdist;
      maxc[0] = maxc[0] + maxdist;
      maxc[1] = maxc[1] + maxdist;
      maxc[2] = maxc[2] + maxdist;
    }
  };

  /** The surface vertices are in RAS. They need to be converted to Vox space
   * because the test voxels are in Vox space.
   * For this we need a mri_template image which is in the same
   * space as the MRIS surface
   * Note that for every voxel we check, we can convert them
   * to RAS and eliminate this method
   * But that is a order of magnitude slower because 
   * the number of voxels to check is usually 
   * much greater than the number of vertex points 
   \param mris - the surface whose vertices have to be converted to vox space
   \param mri_template - the MRI template of the same subject which is needed for the ras2vox call
  */
  void ConvertSurfaceRASToVoxel(MRIS *mris, MRI *mri_template)
  {
    MRISfreeDistsButNotOrig(mris);
      // MRISsetXYZ will invalidate all of these,
      // so make sure they are recomputed before being used again!

    for (int i=0; i< mris->nvertices; i++)
    {
      VERTEX const * const vertex = &mris->vertices[i];
      double cx = vertex->x;
      double cy = vertex->y;
      double cz = vertex->z;

      // for every surface vertex do this call
      double vx, vy, vz;
      MRISsurfaceRASToVoxelCached(mris,
                                  mri_template,
                                  cx, cy, cz,
                                  &vx, &vy, &vz);
      // and reassign the vertices
      MRISsetXYZ(mris,i, (float)vx,(float)vy,(float)vz);
    }
  }

  //! Return the dot product of v1 and v2
  inline double Dot(const double v1[3], const double v2[3])
  {
    return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
  }

  //! Cross product of two 3-vectors. Result vector in z[3].
  inline void Cross(const double x[3], const double y[3], double z[3])
  {
    z[0] = x[1]*y[2] - x[2]*y[1]; 
    z[1] = x[2]*y[0] - x[0]*y[2];
    z[2] = x[0]*y[1] - x[1]*y[0];
  }

  //! Returns the length of the vector x
  inline double Norm(const double x[3])
  {
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  }

  //! Returns the squared-length of the vector
  inline double SqrNorm(const double x[3])
  {
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  }
  
  //! Normalizes the vector so that it's length is 1
  inline double Normalize(double x[3])
  {
    double den; 
    if ( (den = Norm(x)) != 0.0 )
    {
      for (int i=0; i < 3; i++)
      {
        x[i] /= den;
      }
    }
    return den;
  }

  //! Return the normal of the triangle formed by v1,v2 and v3 in 'n'
  inline void ComputeNormal(const double v1[3],
                            const double v2[3],
                            const double v3[3], 
                            double n[3])
  {
    double ax, ay, az, bx, by, bz;
    double length;
    ax = v3[0] - v2[0]; ay = v3[1] - v2[1]; az = v3[2] - v2[2];
    bx = v1[0] - v2[0]; by = v1[1] - v2[1]; bz = v1[2] - v2[2];

    n[0] = (ay * bz - az * by);
    n[1] = (az * bx - ax * bz);
    n[2] = (ax * by - ay * bx);
    if ( (length = Norm(n)) != 0.0 )
    {
      n[0] /= length;
      n[1] /= length;
      n[2] /= length;
    }
  }

  //! Sort the eigen vectors of a 3x3 matrix in descending order based on its eigen values
  //! and return the eigen vectors in max, mid and min arrays
  /*!
    \param eigenSystem - the input vnl_symmetric_eigensystem
    \param evalues - <deprecated>
    \returns max - the maximum eigenvector
    \returns mid - the mid eigenvector
    \returns min - the min eigenvector
  */
  void GetSortedEigenVectors ( vnl_symmetric_eigensystem<double> &eigenSystem,
                               double *evalues,
                               double *max,
                               double *mid,
                               double *min)
  {
    vnl_vector<double> _ev(3);
    int _max=0, _min=0, _mid; 
    double maxval=std::numeric_limits<float>::min(),
      minval=std::numeric_limits<float>::max();

    // sort the eigenvector indices according to the eigenvalues
    for ( int i=0; i<3; i++)
    {
      double curval = eigenSystem.get_eigenvalue(i);
      if ( curval > maxval ) 
      {
        maxval = curval;
        _max = i;
      }
      if ( curval < minval )
      {
        minval = curval;
        _min = i;
      }
    }
    _mid = 3 - _max - _min; // since there are only 3 indices

    _ev = eigenSystem.get_eigenvector(_max);
    max[0] = _ev[0]; max[1] = _ev[1]; max[2] = _ev[2];
    _ev = eigenSystem.get_eigenvector(_mid);
    mid[0] = _ev[0]; mid[1] = _ev[1]; mid[2] = _ev[2];
    _ev = eigenSystem.get_eigenvector(_min);
    min[0] = _ev[0]; min[1] = _ev[1]; min[2] = _ev[2];
  }

  //! Compute distance from x to line segment p1p2.
  // Returns parametric coordinate t
  inline double DistanceToLine( double x[3], double p1[3], double p2[3] )
  {
    double p21[3], denom, num;
    // vectors 
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    //get parametric location
    num   = p21[0]*(x[0] - p1[0]) + p21[1]*(x[1]-p1[1]) + p21[2]*(x[2]-p1[2]);
    denom = p21[0] * p21[0] + p21[1]*p21[1] + p21[2]*p21[2]; //dot product
    
    return (num/denom); 
  }

  // Geometric Tools, LLC
  // Copyright (c) 1998-2010
  // Distributed under the Boost Software License, Version 1.0.
  // http://www.boost.org/LICENSE_1_0.txt
  // http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
  //
  // File Version: 4.10.0 (2009/11/18)
  //! Given the face with vertices v0, v1 and v2.. Find the closest distance
  //! from the face to a point pt
  //! Algorithm - http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
  /*! 
    \param v0, v1, v2 - the vertices of the face
    \param pt - the point 
    \returns the distance in double
  */
  double  DistancePointToFace(VERTEX *v0,
                              VERTEX *v1,
                              VERTEX *v2,
                              double pt[3]) 
  {
    double kDiff[3];
    kDiff[0]= v0->x - pt[0];
    kDiff[1]= v0->y - pt[1];
    kDiff[2]= v0->z - pt[2];
    double kEdge0[3];
    kEdge0[0] = v1->x - v0->x;
    kEdge0[1] = v1->y - v0->y;
    kEdge0[2] = v1->z - v0->z;
    double kEdge1[3];
    kEdge1[0] = v2->x - v0->x;
    kEdge1[1] = v2->y - v0->y;
    kEdge1[2] = v2->z - v0->z;

    double fA00 = SqrNorm(kEdge0);
    double fA01 = Dot(kEdge0, kEdge1);
    double fA11 = SqrNorm(kEdge1);
    double fB0  = Dot(kDiff, kEdge0);
    double fB1  = Dot(kDiff, kEdge1);
    double fC   = SqrNorm(kDiff);
    double fDet = fabs(fA00*fA11-fA01*fA01);
    double fS   = fA01*fB1-fA11*fB0;
    double fT   = fA01*fB0-fA00*fB1;
    double fSqrDistance;

    if (fS + fT <= fDet)
    {
      if (fS < 0.0)
      {
        if (fT < 0.0)  // region 4
        {
          if (fB0 < 0.0)
          {
            fT = 0.0;
            if (-fB0 >= fA00)
            {
              fS = 1.0;
              fSqrDistance = fA00+(2.0)*fB0+fC;
            }
            else
            {
              fS = -fB0/fA00;
              fSqrDistance = fB0*fS+fC;
            }
          }
          else
          {
            fS = 0.0;
            if (fB1 >= 0.0)
            {
              fT = 0.0;
              fSqrDistance = fC;
            }
            else if (-fB1 >= fA11)
            {
              fT = 1.0;
              fSqrDistance = fA11+(2.0)*fB1+fC;
            }
            else
            {
              fT = -fB1/fA11;
              fSqrDistance = fB1*fT+fC;
            }
          }
        }
        else  // region 3
        {
          fS = 0.0;
          if (fB1 >= 0.0)
          {
            fT = 0.0;
            fSqrDistance = fC;
          }
          else if (-fB1 >= fA11)
          {
            fT = 1.0;
            fSqrDistance = fA11+(2.0)*fB1+fC;
          }
          else
          {
            fT = -fB1/fA11;
            fSqrDistance = fB1*fT+fC;
          }
        }
      }
      else if (fT < 0.0)  // region 5
      {
        fT = 0.0;
        if (fB0 >= 0.0)
        {
          fS = 0.0;
          fSqrDistance = fC;
        }
        else if (-fB0 >= fA00)
        {
          fS = 1.0;
          fSqrDistance = fA00+(2.0)*fB0+fC;
        }
        else
        {
          fS = -fB0/fA00;
          fSqrDistance = fB0*fS+fC;
        }
      }
      else  // region 0
      {
        // minimum at interior point
        double fInvDet = (1.0)/fDet;
        fS *= fInvDet;
        fT *= fInvDet;
        fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
          fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
      }
    }
    else
    {
      double fTmp0, fTmp1, fNumer, fDenom;

      if (fS < 0.0)  // region 2
      {
        fTmp0 = fA01 + fB0;
        fTmp1 = fA11 + fB1;
        if (fTmp1 > fTmp0)
        {
          fNumer = fTmp1 - fTmp0;
          fDenom = fA00-2.0f*fA01+fA11;
          if (fNumer >= fDenom)
          {
            fS = 1.0;
            fT = 0.0;
            fSqrDistance = fA00+(2.0)*fB0+fC;
          }
          else
          {
            fS = fNumer/fDenom;
            fT = 1.0 - fS;
            fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
              fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
          }
        }
        else
        {
          fS = 0.0;
          if (fTmp1 <= 0.0)
          {
            fT = 1.0;
            fSqrDistance = fA11+(2.0)*fB1+fC;
          }
          else if (fB1 >= 0.0)
          {
            fT = 0.0;
            fSqrDistance = fC;
          }
          else
          {
            fT = -fB1/fA11;
            fSqrDistance = fB1*fT+fC;
          }
        }
      }
      else if (fT < 0.0)  // region 6
      {
        fTmp0 = fA01 + fB1;
        fTmp1 = fA00 + fB0;
        if (fTmp1 > fTmp0)
        {
          fNumer = fTmp1 - fTmp0;
          fDenom = fA00-(2.0)*fA01+fA11;
          if (fNumer >= fDenom)
          {
            fT = 1.0;
            fS = 0.0;
            fSqrDistance = fA11+(2.0)*fB1+fC;
          }
          else
          {
            fT = fNumer/fDenom;
            fS = 1.0 - fT;
            fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
              fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
          }
        }
        else
        {
          fT = 0.0;
          if (fTmp1 <= 0.0)
          {
            fS = 1.0;
            fSqrDistance = fA00+(2.0)*fB0+fC;
          }
          else if (fB0 >= 0.0)
          {
            fS = 0.0;
            fSqrDistance = fC;
          }
          else
          {
            fS = -fB0/fA00;
            fSqrDistance = fB0*fS+fC;
          }
        }
      }
      else  // region 1
      {
        fNumer = fA11 + fB1 - fA01 - fB0;
        if (fNumer <= 0.0)
        {
          fS = 0.0;
          fT = 1.0;
          fSqrDistance = fA11+(2.0)*fB1+fC;
        }
        else
        {
          fDenom = fA00-2.0f*fA01+fA11;
          if (fNumer >= fDenom)
          {
            fS = 1.0;
            fT = 0.0;
            fSqrDistance = fA00+(2.0)*fB0+fC;
          }
          else
          {
            fS = fNumer/fDenom;
            fT = 1.0 - fS;
            fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
              fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
          }
        }
      }
    }

    // account for numerical round-off error
    if (fSqrDistance < 0.0)
    {
      fSqrDistance = 0.0;
      //fSqrDistance = -fSqrDistance;
    }
    return sqrt(fSqrDistance);
  }

  /**
   * This version of finding a shell does the following
   * Take the max length of edges at a vertex =  r_i . 
   * For a given vertex mark all the voxels less than r_i distance apart
   * as voxels in the shell
   */
  MRI* MRISmaxedgeshell(MRI *mrisrc, MRIS *mris, MRI* mridst, int clearflag)
  {
    /* For each vertex, find the edge with the 
       maximum length and store it in an array */
    for (int vno = 0; vno < mris->nvertices; vno++)
    {
      VERTEX_TOPOLOGY const * const vert_topology = &mris->vertices_topology[vno];
      VERTEX          const * const vert          = &mris->vertices         [vno];
      /* neighbor m */
      double distance = 0.0; // An edge length obviously has to be atleast 0.0
      for (int nc = 0; nc < vert_topology->vnum; nc++)
      {
        VERTEX const * const nvert = &mris->vertices[ vert_topology->v[nc] ];
        double const xd    = vert->x - nvert->x;
        double const yd    = vert->y - nvert->y;
        double const zd    = vert->z - nvert->z;
        double const _dist = sqrt( xd*xd + yd*yd + zd*zd);
        // length of the edge with max length among its neighbors
        if ( _dist > distance )
          distance = _dist;
      }
      for (double xk = -distance; xk <=distance ; xk+=distance/5.0)
        for (double yk = -distance; yk <=distance ; yk+=distance/5.0)
          for (double zk = -distance; zk <=distance ; zk+=distance/5.0)
          {
            if ( xk*xk + yk*yk + zk*zk > distance*distance ) continue;

    	    double fi, fj, fimnr;
            MRISsurfaceRASToVoxelCached(mris, mrisrc,
                                        vert->x + xk,vert->y+yk,vert->z+zk,
                                        &fi,&fj,&fimnr);
            int i=nint(fi);
            int j=nint(fj);
            int imnr=nint(fimnr);
            if (i>=0 && i<mridst->width && 
                j>=0 && j<mridst->height && 
                imnr>=0 && imnr<mridst->depth)
              MRIsetVoxVal(mridst,i,j,imnr,0,255);
          }
    }
    return(mridst);
  }

  /* MRISshell needs mri_src as an example of 
   * the format for the output, to match size and type.
   * The surface is recreated in the MRI space (mri_dst) 
   * from the tesselated surface (mris) as voxels of 255 */
  MRI *MRISbshell(MRI *mri_src,MRI_SURFACE *mris,MRI *mri_dst,int clearflag)
  {
    int width,height,depth,j,imnr, fno;
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    double px,py,pz;
    double a, b, c;
    double fi,fj,fimnr;
    VERTEX *v_0,*v_1,*v_2;
    FACE *f;


    /* Create new blank MRI or clear existing destination MRI */
    width = mri_src->width;
    height = mri_src->height;
    depth = mri_src->depth;
    if (!mri_dst)
    {
      /*    printf("MRISshell: Creating new (_dst)MRI...\n");*/
      mri_dst = MRIalloc(width, height, depth, mri_src->type);
      MRIcopyHeader(mri_src, mri_dst);
    }

    if (clearflag)
      MRIclear(mri_dst);

    /* Fill each face in MRI volume */
    for (fno=0; fno<mris->nfaces; fno++)
    {
      /* Calculate (x,y,z) for each vertex for face */
      f = &mris->faces[fno];
      v_0 = &mris->vertices[f->v[0]];
      v_1 = &mris->vertices[f->v[1]];
      v_2 = &mris->vertices[f->v[2]];
      if (f->v[0] == Gdiag_no ||
          f->v[1] == Gdiag_no ||
          f->v[2] == Gdiag_no)
        DiagBreak();
      x0 = v_0->x;
      y0 = v_0->y;
      z0 = v_0->z;
      x1 = v_1->x;
      y1 = v_1->y;
      z1 = v_1->z;
      x2 = v_2->x;
      y2 = v_2->y;
      z2 = v_2->z;

      //Points generation using Barycentric coordinates
      for (int i=0; i < 50; i++ )
      {
        a = randomNumber(0,1);
        b = randomNumber(0,1);
        if ( a+b > 1)
        {
          a = 1 - a;
          b = 1 - b;
        }
        c = 1 - a - b;
        
        px = a * x0 + b * x1 + c * x2;
        py = a * y0 + b * y1 + c * y2;
        pz = a * z0 + b * z1 + c * z2;
      
        MRISsurfaceRASToVoxelCached(mris, mri_src,px,py,pz,&fi,&fj,&fimnr);
        i=nint(fi);
        j=nint(fj);
        imnr=nint(fimnr);
        if (i>=0 && i<mri_dst->width && 
            j>=0 && j<mri_dst->height && 
            imnr>=0 && imnr<mri_dst->depth)
          MRIsetVoxVal(mri_dst,i,j,imnr,0,255);
      }
    }

    return mri_dst;
  }

}

#endif //utilsmath_h
