/*
 * Original Author: Ruopeng Wang
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
#ifndef _Track_h
#define _Track_h
#include <QList>

struct Track
{
  int   nNum;
  float*  fPts;
  QList<float*> fScalars;
  float*  fProperty;

  float fLength;

  short* nVoxels;
  int    nNumberOfVoxels;
  unsigned char charColor[3];

  Track()
  {
    Reset();
  }

  // NOTE: destructor does NOT delete allocated data
  //       Have to call Delete() explicitly to free the allocated data
  void Delete();

  void Reset();

  bool Update(const short* dim, const float* voxel_size);

  bool AddScalars(const short* dim, const float* voxel_size, const float* value_ptr);

  bool RemoveScalars(int nIndex);

  static float GetTrackLength(const float* track_pts, int n);

  static void GetVoxels(const float* track_pts_in, int n_pts_in, const short* dim, const float* voxel_size,
                        short** indices_out, int* num_out);

  void GetVoxels(const short* dim, const float* voxel_size, short** indices_out, int* num_out);

  float GetMeanScalar(int nIndex);

  inline float GetProperty(int nIndex)
  {
    return fProperty[nIndex];
  }


  static double DistanceBetween2Points( float* track_pts, int n1, int n2 );

  /*
  static bool IntersectWithPlane(float* track_pts, int nNum, double* normal, double* pt,
              double* x, int* index = 0);
  static bool IntersectWithDisk(float* track_pts, int nNum, double* normal, double* pt, double r,
                 double* x, int* indices, double angle_end, double angle_start);

  bool IntersectWithPlane(double* normal, double* pt, double* x, int* index);

  bool IntersectWithDisk(double* normal, double* pt, double r,
                 double* x, double angle_end, double angle_start);
  */
};


#endif
