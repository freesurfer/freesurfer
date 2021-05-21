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
#include "Track.h"
#include "MyUtils.h"

// NOTE: destructor doesn't delete allocated data
//       Have to call Delete() explicitly to free the allocated data
void Track::Delete()
{
  delete[] fPts;
  delete[] nVoxels;
  delete[] fProperty;
  for (int i = 0; i < fScalars.size(); i++)
  {
    delete[] fScalars[i];
  }

  Reset();
}

void Track::Reset()
{
  nNum = 0;
  fPts = 0;
  nVoxels = 0;
  nNumberOfVoxels = 0;
  fLength = 0;
  fProperty = 0;
  fScalars.clear();
}

bool Track::Update(const short* dim, const float* voxel_size)
{
  fLength = GetTrackLength(fPts, nNum);

  if (nVoxels)
  {
    delete[] nVoxels;
    nVoxels = NULL;
    nNumberOfVoxels = 0;
  }

  GetVoxels(dim, voxel_size, &nVoxels, &nNumberOfVoxels);

  return true;
}

bool Track::AddScalars(const short* dim, const float* voxel_size, const float* value_ptr)
{
  float* p = fPts;
  int x, y, z;
  float* v_ptr = new float[nNum];

  for (int i = 0; i < nNum; i++)
  {
    x = qMax(0, qMin(dim[0]-1, (int)(p[0]/voxel_size[0])));
    y = qMax(0, qMin(dim[1]-1, (int)(p[1]/voxel_size[1])));
    z = qMax(0, qMin(dim[2]-1, (int)(p[2]/voxel_size[2])));

    v_ptr[i] = value_ptr[z*dim[0]*dim[1]+y*dim[0]+x];
    //  qDebug() << v_ptr[i];
    p += 3;
  }

  fScalars.push_back(v_ptr);

  return true;
}

bool Track::RemoveScalars(int nIndex)
{
  if (nIndex >= fScalars.size())
  {
    return false;
  }
  delete[] fScalars[nIndex];
  fScalars.erase(fScalars.begin()+nIndex);
  return true;
}

float Track::GetTrackLength(const float* track_pts, int n)
{
  double len = 0;

  for (int i = 1; i < n; i++)
  {
    len += MyUtils::GetDistance<float>(track_pts+i*3, track_pts+i*3-3);
  }
  return (float)len;
}

void GetVoxels( short* n1, short* n2, QList<short>& indices_out )
{
  size_t k;
  size_t n_out = indices_out.size()/3;
  for (k = 0; k < n_out; k++)
  {
    if (n1[0] == indices_out[k*3] && n1[1] == indices_out[k*3+1] && n1[2] == indices_out[k*3+2])
    {
      break;
    }
  }
  if ( k >= n_out )
  {
    for ( int i = 0; i < 3; i++ )
    {
      indices_out.push_back( n1[i] );
    }
  }

  if ( qAbs( n1[0]-n2[0]) > 1 || qAbs( n1[1]-n2[1]) > 1 || qAbs( n1[2]-n2[2]) > 1 )
  {
    short m[3];
    for ( int i = 0; i < 3; i++ )
    {
      m[i] = (n1[i] + n2[i] ) /2;
    }

    GetVoxels( n1, m, indices_out );
    GetVoxels( m, n2, indices_out );
  }
  else
  {
    if ( n1[0] != n2[0] || n1[1] != n2[1] || n1[2] != n2[2] )
    {
      for ( int i = 0; i < 3; i++ )
      {
        indices_out.push_back( n2[i] );
      }
    }
  }
}

void Track::GetVoxels(const float* track_pts_in, int n_pts_in, const short* dim, const float* voxel_size,
                      short** indices_out, int* n_out)
{
  short n[3], m[3];
  QList<short> indices;
  for (int i = 0; i < n_pts_in-1; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      n[j] = qMax(0, qMin(dim[j]-1, (int)(track_pts_in[i*3+j]/voxel_size[j])));
      m[j] = qMax(0, qMin(dim[j]-1, (int)(track_pts_in[i*3+3+j]/voxel_size[j])));
    }
    ::GetVoxels( n, m, indices );
  }
  *n_out = indices.size() / 3;
  *indices_out = new short[indices.size()];
  for ( int i = 0; i < indices.size(); i++ )
  {
    (*indices_out)[i] = indices[i];
  }
}

void Track::GetVoxels(const short* dim, const float* voxel_size, short** indices_out, int* n_out)
{
  GetVoxels(fPts, nNum, dim, voxel_size, indices_out, n_out);
}

/*
bool Track::IntersectWithPlane(float* track_pts, int nNum, double* normal, double* pt, double* x, int* index)
  {
    for (int i = 1; i < nNum; i++)
    {
      double t;
      if (CMyMath::PlaneIntersectWithLine(track_pts+i*3, track_pts+(i-1)*3, normal, pt, t, x))
      {
        if (index)
          *index = i-1;
        return true;
      }
    }

    return false;
  }

bool Track::IntersectWithDisk(float* track_pts, int nNum, double* normal, double* pt, double r,
                   double* x, int* indices, double angle_end, double angle_start)
{
  for ( int i = 0; i < nNum-1; i++ )
  {
    double t;
    if (CMyMath::PlaneIntersectWithLine(track_pts+i*3, track_pts+(i+1)*3, normal, pt, t, x))
    {
      if ( CMyMath::Distance2BetweenPoints(pt, x) <= r*r)
      {
        if ( angle_end < 90 || angle_start > 0)
        {
          float v[3];
          CMyMath::PtsToVector(track_pts+i*3, track_pts+i*3+3, v);
          float a = CMyMath::AngleBetweenVectors(v, normal)*180/PI;
          if (a > 90)
            a = 180-a;
          indices[0] = i;
          indices[1] = i+1;
          return (a >= angle_start && a <= angle_end);
        }
        else
        {
          indices[0] = i;
          indices[1] = i+1;
          return true;
        }
      }
    }
  }
  return false;
}

bool Track::IntersectWithPlane(double* normal, double* pt, double* x, int* index)
{
  return IntersectWithPlane(fPts, nNum, normal, pt, x, index);
}

bool Track::IntersectWithDisk(double* normal, double* pt, double r,
                   double* x, double angle_end, double angle_start)
{
  int n[2];
  return IntersectWithDisk(fPts, nNum, normal, pt, r, x, n, angle_end, angle_start);
}

*/

float Track::GetMeanScalar(int nIndex)
{
  float* p = fScalars[nIndex];
  double len = 0;
  double dMean = 0;
  for (int i = 0; i < nNum-1; i++)
  {
    //  qDebug() << p[i] << p[i+1];
    double dLen = MyUtils::GetDistance<float>(fPts + i*3, fPts + i*3+3);
    //  qDebug() << dLen;
    dMean += dLen * (p[i] + p[i+1]) / 2;
    len += dLen;
  }
  return (float)(dMean/len);
}

double Track::DistanceBetween2Points( float* track_pts, int n1, int n2 )
{
  double d = 0;
  for ( int i = n1; i < n2-1; i++ )
  {
    d += MyUtils::GetDistance<float>(track_pts + i*3, track_pts + i*3+3);
  }
  return d;
}
