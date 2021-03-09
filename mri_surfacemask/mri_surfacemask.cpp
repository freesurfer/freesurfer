/**
 * @brief taking mri volume and surface and create volume
 *
 * where all pixels outside the surface has been zeroed
 *
 */
/*
 * Original Author: Yasunari Tosa
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

#include <iostream>
#include <iomanip>

 
#include "mri.h"
#include "mrisurf.h"
#include "version.h"
  const char *Progname = "mri_surfacemask";


using namespace std;

void markSurfaceVoxels(MRI *mri_filled, MRIS *mris) {
  double x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax,u,v;
  double px,py,pz,px0,py0,pz0,px1,py1,pz1;
  int numu,numv;
  double tx,ty,tz;
  int xv, yv, zv;

  int width = mri_filled->width;
  int height = mri_filled->height;
  int depth = mri_filled->depth;

  for (int k=0; k < mris->nfaces; k++) {
    // get three vertices
    x0 =mris->vertices[mris->faces[k].v[0]].x;
    y0 =mris->vertices[mris->faces[k].v[0]].y;
    z0 =mris->vertices[mris->faces[k].v[0]].z;

    x1 =mris->vertices[mris->faces[k].v[1]].x;
    y1 =mris->vertices[mris->faces[k].v[1]].y;
    z1 =mris->vertices[mris->faces[k].v[1]].z;

    x2 =mris->vertices[mris->faces[k].v[2]].x;
    y2 =mris->vertices[mris->faces[k].v[2]].y;
    z2 =mris->vertices[mris->faces[k].v[2]].z;
    // calculate sides
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
    numu = int(ceil(2*d0));
    numv = int(ceil(2*dmax));

    /*         x0
    //        /  \
    //       /    \ px0
    //      /     /\
    //     /   px/  \
    //  x1/_____/____\x2
    //         px1
    */
    for (v=0;v<=numv;v++) {
      // px0 spans x0 to x2
      px0 = x0 + (x2-x0)*v/numv;
      py0 = y0 + (y2-y0)*v/numv;
      pz0 = z0 + (z2-z0)*v/numv;
      // px1 spans x1 to x2
      px1 = x1 + (x2-x1)*v/numv;
      py1 = y1 + (y2-y1)*v/numv;
      pz1 = z1 + (z2-z1)*v/numv;
      for (u=0;u<=numu;u++) {
        // px spans px0 to px1
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;

        // calculate voxel value of a point in the triangle
        // C function don't know about const
        MRIsurfaceRASToVoxel(const_cast<MRI *> (mri_filled),px,py,pz,&tx,&ty,&tz);
        // voxel position is floating point
        // change it to int.  we use the convention of -0.5 to 0.5 as voxel 0
        // thus floor(tx + 0.5) i.e. -0.5 -> -0.5+0.5 = 0
        xv = int(floor(tx + 0.5));
        yv = int(floor(ty + 0.5));
        zv = int(floor(tz + 0.5));
        if (xv >= 0 && xv < width && 
            yv >=0 && yv < height && 
            zv >= 0 && zv < depth) {
          MRIvox(mri_filled,xv,yv,zv)  = 255;
        }
      }  // numu
    } // numv
  } // nsurface
}

// mark outside voxels to be 64.  Only modify 0'ed voxels
void markOutsideVoxels(MRI *mri_filled) {
  // now we mark outside voxel as 64
  int totalfilled, newfilled;
  int k,j,i;
  int width = mri_filled->width;
  int height = mri_filled->height;
  int depth = mri_filled->depth;

  MRIvox(mri_filled,1,1,1)= 64;
  totalfilled = newfilled = 1;
  while (newfilled>0) {
    newfilled = 0;
    // going from left to right
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++)
          if (MRIvox(mri_filled,i,j,k)==0)
            if (MRIvox(mri_filled,i,j,mri_filled->zi[k-1])==64||
                MRIvox(mri_filled,i,mri_filled->yi[j-1],k)==64||
                MRIvox(mri_filled,mri_filled->xi[i-1],j,k)==64) {
              MRIvox(mri_filled,i,j,k)= 64;
              newfilled++;
            }
    // going from right to left
    for (k=depth-1;k>=0;k--)
      for (j=height-1;j>=0;j--)
        for (i=width-1;i>=0;i--)
          if (MRIvox(mri_filled,i,j,k)==0)
            if (MRIvox(mri_filled,i,j,mri_filled->zi[k+1])==64||
                MRIvox(mri_filled,i,mri_filled->yi[j+1],k)==64||
                MRIvox(mri_filled,mri_filled->xi[i+1],j,k)==64) {
              MRIvox(mri_filled,i,j,k) = 64;
              newfilled++;
            }
    totalfilled += newfilled;
  }
  // fill all volume surface boundary voxels to be 64 (there are 6 faces)
  for (k=0; k < depth;k++)
    for (j=0; j < height; j++) {
      MRIvox(mri_filled,       0, j, k ) = 64;
      MRIvox(mri_filled, width-1, j, k ) = 64;
    }

  for (k=0; k < depth;k++)
    for (i=0; i < width ; i++) {
      MRIvox(mri_filled, i,        0, k ) = 64;
      MRIvox(mri_filled, i, height-1, k ) = 64;
    }

  for (i=0; i < width ;i++)
    for (j=0; j < height; j++) {
      MRIvox(mri_filled, i, j,      0 ) = 64;
      MRIvox(mri_filled, i, j, depth-1) = 64;
    }
}

int main(int argc, char *argv[]) {
  int nargs = handleVersionOption(argc, argv, "mri_surfacemask");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;
  fprintf(stderr,"\n");

  if (argc < 4) {
    cout << "  Usage: " << Progname << " <volume> <surface> <surface-masked volume>" << endl;
    cout << "Purpose: Produce a new volume where all pixels outside the surface are set to zero."
    << endl;
    exit(-1);
  }

  // readin volume
  MRI *mri_input = MRIread(argv[1]);
  if (!mri_input) {
    cout << "Could not read volume " << argv[1] << endl;
    exit(-1);
  }
  // readin surface
  MRIS *mris = MRISread(argv[2]);
  if (!mris) {
    cout << "Could not read surface " << argv[2] << endl;
    MRIfree(&mri_input);
    exit(-1);
  }
  int width = mri_input->width;
  int height = mri_input->height;
  int depth = mri_input->depth;

  // copy the input volume as the output. initialized to be zero
  MRI *mri_filled = MRIalloc(width, height, depth, MRI_UCHAR);
  if (!mri_filled) {
    MRIfree(&mri_input);
    MRISfree(&mris);
    cout << "could not allocate memory for working buffer" << endl;
    exit(-1);
  }
  MRIcopyHeader(mri_input, mri_filled);

  // we marked all voxels which touches the surface to be 255
  markSurfaceVoxels(mri_filled, mris);

  // free surface memory
  MRISfree(&mris);

  // mark outside voxels to be 64
  markOutsideVoxels(mri_filled);

  // prepare the out put
  MRI *mri_out = MRIcopy(mri_input, NULL);
  if (!mri_out) {
    MRIfree(&mri_input);
    MRIfree(&mri_filled);
    cout << "could not allocate memory for output buffer" << endl;
    exit(-1);
  }
  // mask it (make voxel with value of 64 to be made 0)
  MRImask(mri_out, mri_filled, mri_out, 64, 0);

  // write it
  MRIwrite(mri_out, argv[3]);
  //
  MRIfree(&mri_input);
  MRIfree(&mri_filled);
  MRIfree(&mri_out);
}
