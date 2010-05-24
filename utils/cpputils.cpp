/**
 * @file  cpputils.cpp
 * @brief various utilities which make use of C++ classes
 *
 */
/*
 * Original Author: Krish Subramaniam
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2010/05/24 15:36:11 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2010
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

char *CPPUTILS_VERSION = "$Revision: 1.1 $";

extern "C"
{
#include "diag.h"
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "cma.h"
}
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <queue>

#include "MRISOBBTree.h"
#include "MRISdistancefield.h"


typedef Math::Point<int> Pointd;
typedef Math::Point<float> Point_f;
using namespace std;

/* This version of fill interior uses the surface distance transform and an inside outside check using OBBTree datastructure to accurately fill. hence this is generally a more accurate algorithm than the old MRISfillInterior but the downside is it takes more time ( around 5 minutes ) . 
 *  returned MRI structure (mri_out) is an MRI_INT MRI volume
 */ 
extern "C" MRI *
MRISfillInterior2(MRI_SURFACE *mris, MRI *mri_interior)
{
  int res;
  MRI *mri_visited, *_mridist, *mri_out; 

  _mridist = MRIcloneDifferentType(mri_interior, MRI_FLOAT);
  mri_visited  = MRIcloneDifferentType(_mridist, MRI_INT);
  mri_out  = MRIclone(mri_interior, NULL);

  // Save the MRIS surface verts
  Point_f *surfverts = new Point_f[mris->nvertices];
  for (int vc=0; vc<mris->nvertices; vc++)
  {
    VERTEX *v = &mris->vertices[vc];
    surfverts[vc].v[0] = v->x;
    surfverts[vc].v[1] = v->y;
    surfverts[vc].v[2] = v->z;
  }

  // Convert surface vertices to vox space
  Math::ConvertSurfaceRASToVoxel(mris, _mridist);
  
  // Find the distance field 
  MRISDistanceField *distfield = new MRISDistanceField(mris, _mridist);
  distfield->SetMaxDistance(3.0);
  distfield->Generate(); //mri_dist now has the distancefield 

  // Construct the OBB Tree
  MRISOBBTree* OBBTree = new MRISOBBTree(mris);
  OBBTree->ConstructTree();
  
  std::queue<Pointd* > ptsqueue;
  // iterate through all the volume points 
  // and apply sign 
  for(int i=0; i< _mridist->width; i++)
  {
    //std::cerr << i <<" ";
    for(int j=0; j< _mridist->height; j++)
    {
      for(int k=0; k< _mridist->depth; k++)
      {
        if ( MRIIvox(mri_visited, i, j, k )) continue; 
        res = OBBTree->PointInclusionTest(i, j, k);
        Pointd *pt = new Pointd;
        pt->v[0] = i; pt->v[1] = j; pt->v[2] = k;
        ptsqueue.push( pt );

        // First serve all the points in the queue before going to the next voxel
        while ( !ptsqueue.empty() )
        {
          // serve the front and pop it
          Pointd *p = ptsqueue.front();
          const int x = p->v[0];
          const int y = p->v[1];
          const int z = p->v[2];
          delete p;
          ptsqueue.pop();

          if ( MRIIvox(mri_visited, x, y, z) ) continue; 
          MRIIvox(mri_visited, x, y, z) =  res;
          const float dist = MRIFvox(_mridist, x, y, z);
          MRIFvox(_mridist, x, y, z) =  dist*res;

          // mark its 6 neighbors if distance > 1 ( triangle inequality )
          if ( dist > 1 )
          {
            // left neighbor in x
            if ( x>0 && !MRIIvox(mri_visited, x-1, y, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x - 1;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // bottom neighbor in y
            if ( y>0 && !MRIIvox(mri_visited, x, y-1, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y - 1;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // front neighbor in z
            if ( z>0 && !MRIIvox(mri_visited, x, y, z-1))
            { 
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z - 1;
              ptsqueue.push(ptemp);
            }
            // right neighbor in x
            if ( x<mri_visited->width-1 && !MRIIvox(mri_visited, x+1, y, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x + 1;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // top neighbor in y
            if ( y<mri_visited->height-1 && !MRIIvox(mri_visited, x, y+1, z))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y + 1;
              ptemp->v[2]   = z;
              ptsqueue.push(ptemp);
            }
            // back neighbor in z
            if ( z<mri_visited->depth-1 && !MRIIvox(mri_visited, x, y, z+1))
            {
              Pointd *ptemp = new Pointd;
              ptemp->v[0]   = x;
              ptemp->v[1]   = y;
              ptemp->v[2]   = z + 1;
              ptsqueue.push(ptemp);
            }
          }
        }
      }
    }
  }
  for(int i=0; i< _mridist->width; i++)
  {
    for(int j=0; j< _mridist->height; j++)
    {
      for(int k=0; k< _mridist->depth; k++)
      {
        if ( MRIFvox(_mridist, i, j, k) > 0)
          MRIsetVoxVal(mri_out, i, j, k, 0, 1);
        else
          MRIsetVoxVal(mri_out, i, j, k, 0, 0);
      }
    }
  }

  // Restore the saved MRIS surface verts
  for (int vc=0; vc<mris->nvertices; vc++)
  {
    VERTEX *v = &mris->vertices[vc];
    v->x = surfverts[vc].v[0];
    v->y = surfverts[vc].v[1];
    v->z = surfverts[vc].v[2];
  }

  delete [] surfverts;  //no longer need
  MRIfree(&mri_visited);
  MRIfree(&_mridist);
  delete OBBTree;
  delete distfield;
  return(mri_out);

}

