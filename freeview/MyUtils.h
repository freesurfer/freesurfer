/**
 * @brief Misc utility class.
 *
 */
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
 *
 */

#ifndef MyUtils_h
#define MyUtils_h

#include <math.h>
#include <vector>
#include <QStringList>

class MyUtils
{
public:
  MyUtils()
  {}

  static void FloodFill( char** data, int x, int y, int min_x, int min_y, int max_x, int max_y, int fill_value, int border_value );
  static void FloodFill( char* data, int x, int y, int dim_x, int dim_y, int fill_value, int border_value );

  template <class T> static double GetDistance( const T* pt1, const T* pt2 );
  template <class T> static void GetVector( T* pt1, T* pt2, double* v_out, bool bNormalize = true);
  template <class T> static bool Equal( T* pt1, T* p2, int n = 3 );
  template <class T> static T** AllocateMatrix(int ny, int nx);
  template <class T> static void FreeMatrix(T** p, int ny);
  template <class T> static double Dot( T* v1, T* v2 );


  static double RoundToGrid( double dvalue );

  static bool CalculateOptimalVolume( int* vox, int nsize1, int* vox2, int nsize2,
                                      std::vector<void*> input_volumes, float* output_volume, int vol_size );
  static bool CalculateOptimalVolume( int* vox, int nsize1, int* vox2, int nsize2,
                                      std::vector<void*> input_volumes, double* output_volume, int vol_size );
  static bool CalculateOptimalVolume( int* vox, int nsize1, int* vox2, int nsize2,
                                      std::vector<void*> input_volumes, int* output_volume, int vol_size );
  static bool CalculateOptimalVolume( int* vox, int nsize1, int* vox2, int nsize2,
                                      std::vector<void*> input_volumes, short* output_volume, int vol_size );
  static bool CalculateOptimalVolume( int* vox, int nsize1, int* vox2, int nsize2,
                                      std::vector<void*> input_volumes, unsigned char* output_volume, int vol_size );
  static bool CalculateOptimalVolume( int* vox, int nsize1, int* vox2, int nsize2,
                                      std::vector<void*> input_volumes, long* output_volume, int vol_size );

  static double CalculateCorrelationCoefficient(float* x, float* y, int n);

  static bool IsIdentity( double m[4][4], double scale = 1.0 );

  static bool IsOblique( double m[4][4] );

  static QStringList SplitString( const QString& strg, const QString& divider, int nIgnoreStart = 0, int nIgnoreLength = 0 );

  static QString CygwinPathToWin32Path(const QString& path_in);

  static QString Win32PathToCygwinPath(const QString& path_in);

  static QString NormalizeCygwinPath(const QString& path_in);

  static QString CygwinPathProof(const QString& path_in);

  static QString Win32PathProof(const QString& path_in);

  static bool FindIntersection(const double& x0, const double& y0,
                               const double& x1, const double& y1,
                               const double& a0, const double& b0,
                               const double& a1, const double& b1,
                               double* x, double* y);

  static bool FindIntersection(std::vector < std::vector < double > >& line0,
                               std::vector < std::vector < double > >& line1,
                               double* x, double* y, int* n0 = NULL, int* n1 = NULL);

  static QString RealToNumber(qreal val, int nPrecision);
};

template <class T>
double MyUtils::GetDistance( const T* pt1, const T* pt2 )
{
  double dRes = 0;
  for ( int i = 0; i < 3; i++ )
  {
    dRes += ( pt2[i] - pt1[i] ) * ( pt2[i] - pt1[i] );
  }

  return sqrt( dRes );
}

template <class T>
double MyUtils::Dot( T* v1, T* v2 )
{
  return ( double )( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] );
}

template <class T>
void MyUtils::GetVector( T* pt1_in, T* pt2_in, double* v_out, bool bNormalize )
{
  T pt1[3], pt2[3];
  for ( int i = 0; i < 3; i++ )
  {
    pt1[i] = pt1_in[i];
    pt2[i] = pt2_in[i];
  }

  double dist = GetDistance( pt1, pt2 );
  if ( dist == 0 )
  {
    return;
  }

  for ( int i = 0; i < 3; i++ )
  {
    v_out[i] = ( (double)pt2[i] - (double)pt1[i] );
    if (bNormalize)
      v_out[i] /= dist;
  }
}

template <class T>
bool MyUtils::Equal( T* pt1, T* pt2, int nLength )
{
  for ( int i = 0; i < nLength; i++ )
  {
    if ( pt1[i] != pt2[i] )
    {
      return false;
    }
  }

  return true;
}

template <class T>
T** MyUtils::AllocateMatrix(int ny, int nx)
{
  T** p = new T*[ny];
  if (!p)
  {
    return 0;
  }

  for (int i = 0; i < ny; i++)
  {
    p[i] = new T[nx];
    if (!p[i])
    {
      for (int j = 0; j < i; j++)
      {
        delete[] p[j];
      }
      delete[] p;
      return 0;
    }
  }

  return p;
}

template <class T>
void MyUtils::FreeMatrix(T** p, int ny)
{
  if (p == 0)
  {
    return;
  }
  for (int i = 0; i < ny; i++)
  {
    delete[] p[i];
  }
  delete[] p;
  p = 0;
}

inline bool MyUtils::IsIdentity( double m[4][4], double scale )
{
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      if (i == j && qAbs(m[i][j]-1) > 1e-6)
        return false;
      else if (i != j && qAbs(m[i][j]) > scale*1e-4)
        return false;
    }
  }
  return true;
}

inline bool MyUtils::IsOblique( double m[4][4] )
{
  double near_zero = 1e-5;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (fabs(m[i][j]) > near_zero && fabs(fabs(m[i][j])-1) > near_zero)
      {
        return true;
      }
    }
  }
  return false;
}

#endif
