/**
 * @brief Given a MRIS surface, generate a distance field where the value at a voxel is the 
 * distance between the voxel and the closest point on the surface
 *
 * Original Author: Krish Subramaniam 
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




// include guard
#ifndef mrisdistancefield_h
#define mrisdistancefield_h

#include "mri.h"
#include "mrisurf.h"
#include "diag.h"

MRI *MRISdistancefield(MRIS *mris, MRI *mri_tmp, double max_distance, int signedfield);

#include <iostream>
#include "utilsmath.h"

/** This algorithm is based on "Generating Signed Distance Fields from Triangle Meshes"
 * by J.Baerentzen and Henrik Aanaes 
 * For a given mesh with N triangles and X, Y, Z number of voxels in x, y and z direction the complexity
 * of the bruce force approach is O(N*X*Y*Z).
 * This algorithm differs in the sense it calculates distances for voxels only within a bounding box of each face. 
 * Hence its complexity is around O(N*k) where k is the size of the bounding box. ( distance factor )
 * Thus, the distance factor is critical to how quickly the algorithm executes. 
 * If you are only interested in the transition region make sure the distance factor you give is the minimum needed.
 */
class MRISDistanceField
{
  private:
    MRIS *mris;           // The MRIS surface to which the distance field is generated
    MRI *mri_distfield;   // This has the distance field output once this is run
    double max_distance;  // calculate distance field only for voxels less than max_distance far
    bool signedfield;     // TODO - whether we need a signed field or not

  public:
    //! Constructor needs the surface and the mri
    MRISDistanceField(MRIS *_mris, MRI *_mri_distfield) : mris(_mris), mri_distfield(_mri_distfield)
    {
      this->max_distance = 3.0; //default max distance
      this->signedfield = false;
    }

    //! Setter for the max distance
    void SetMaxDistance(double dist)
    {
      this->max_distance = dist;
    }

    //! Set if a signed field is desired ( NOT IMPLEMENTED )
    void SetSignedField()
    {
      this->signedfield = true;
    }

    //! Initially set all voxels in the volume as max_distance
    void InitVolume()
    {
      // If unsigned field, init the volume so that all voxels are max_distance
      for (int i=0; i < mri_distfield->width; i++)
      {
        for (int j=0; j < mri_distfield->height; j++)
        {
          for (int k=0; k < mri_distfield->depth; k++)
          {
            MRIsetVoxVal(mri_distfield, i, j, k, 0, (float)max_distance);
          }
        }
      }
    }


    //! Generate the distance field
    void Generate()
    {
      FACE *face;
      VERTEX *v0, *v1, *v2;
      Math::BoundingBox<int> bbox;
      double calcdist;

      // Init the volume so that all voxels are max_distance
      InitVolume();

      // For every face, 
      for ( int fno=0; fno < mris->nfaces; fno++)
      {
        face = &mris->faces[fno];

        // get its integer bounding box adjusted by the distance
        computeBB(face, bbox);

        // bounds check -- the way surfaces are located in the brain
        // these checks will not ever be true
        if (bbox.minc[0] < 0) bbox.minc[0] = 0;
        if (bbox.minc[1] < 0) bbox.minc[0] = 0;
        if (bbox.minc[2] < 0) bbox.minc[0] = 0;
        if (bbox.maxc[0] >= mri_distfield->width)  bbox.maxc[0] = mri_distfield->width-1;
        if (bbox.maxc[1] >= mri_distfield->height) bbox.maxc[1] = mri_distfield->height-1;
        if (bbox.maxc[2] >= mri_distfield->depth)  bbox.maxc[2] = mri_distfield->depth-1;

        // for all voxels in the bounding box
        for ( int i=bbox.minc[0]; i<=bbox.maxc[0]; i++)
        {
          for ( int j=bbox.minc[1]; j<=bbox.maxc[1]; j++)
          {
            for ( int k=bbox.minc[2]; k<=bbox.maxc[2]; k++)
            {
              // calculate distance to the face
              double pt[3];
              pt[0] = (double)i; pt[1] = (double)j; pt[2] = (double)k;
              v0 = &mris->vertices[face->v[0]];
              v1 = &mris->vertices[face->v[1]];
              v2 = &mris->vertices[face->v[2]];
              calcdist = Math::DistancePointToFace(v0, v1, v2, pt); 
              
              // if distance to face is smaller than previous voxel-value, update the voxel value
              if ( calcdist < MRIgetVoxVal(mri_distfield, i, j, k, 0) )
                MRIsetVoxVal(mri_distfield, i, j, k, 0, (float)calcdist);
            }
          }
        }
      }
    }

    //! Compute the bounding box of the face scaled by the max_distance
    //! return the integer bounding box ( so that it's easy to iterate )
    void computeBB(FACE *face, Math::BoundingBox<int>& bbox)
    {
      VERTEX *v1, *v2, *v3;
      Math::BoundingBox<double> _bb;

      v1 = &mris->vertices[face->v[0]];
      v2 = &mris->vertices[face->v[1]];
      v3 = &mris->vertices[face->v[2]];

      _bb.maxc[0] = std::max(std::max(v1->x, v2->x), v3->x);
      _bb.maxc[1] = std::max(std::max(v1->y, v2->y), v3->y);
      _bb.maxc[2] = std::max(std::max(v1->z, v2->z), v3->z);
      _bb.minc[0] = std::min(std::min(v1->x, v2->x), v3->x);
      _bb.minc[1] = std::min(std::min(v1->y, v2->y), v3->y);
      _bb.minc[2] = std::min(std::min(v1->z, v2->z), v3->z);

      _bb.adjustbounds(this->max_distance); //scale the bounding box

      bbox.minc[0] = nint(_bb.minc[0]);
      bbox.minc[1] = nint(_bb.minc[1]);
      bbox.minc[2] = nint(_bb.minc[2]);
      bbox.maxc[0] = nint(_bb.maxc[0]);
      bbox.maxc[1] = nint(_bb.maxc[1]);
      bbox.maxc[2] = nint(_bb.maxc[2]);
    }

};

#endif

