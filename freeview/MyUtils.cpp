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

#include "MyUtils.h"
#include <math.h>
#include <stddef.h>
#include <QFileInfo>
#include <QDir>
#include <QDebug>
#include <iostream>


#include "matrix.h"


using namespace std;

typedef struct
{
  int xl, xr, y, dy;
}
LINESEGMENT;

#define MAXDEPTH 10000

#define PUSH(XL, XR, Y, DY) \
  if( sp < stack+MAXDEPTH && Y+(DY) >= min_y && Y+(DY) <= max_y ) \
{ sp->xl = XL; sp->xr = XR; sp->y = Y; sp->dy = DY; sp++; }

#define POP(XL, XR, Y, DY) \
{ sp--; XL = sp->xl; XR = sp->xr; Y = sp->y + (DY = sp->dy); }

void MyUtils::FloodFill(char** data, int x, int y,
                        int min_x, int min_y,
                        int max_x, int max_y,
                        int fill_value, int border_value)
{
  int left, x1, x2, dy;
  LINESEGMENT stack[MAXDEPTH], *sp = stack;

  if (data[y][x] == border_value || data[y][x] == fill_value)
  {
    return;
  }

  if (x < min_x || x > max_x || y < min_y || y > max_y)
  {
    return;
  }

  PUSH(x, x, y, 1);        /* needed in some cases */
  PUSH(x, x, y+1, -1);    /* seed segment (popped 1st) */

  while (sp > stack )
  {
    POP(x1, x2, y, dy);

    for (x = x1; x >= min_x &&
         data[y][x] != border_value && data[y][x] != fill_value; x--)
    {
      data[y][x] = fill_value;
    }

    if ( x >= x1 )
    {
      goto SKIP;
    }

    left = x+1;
    if ( left < x1 )
    {
      PUSH(left, x1-1, y, -dy);  /* leak on left? */
    }

    x = x1+1;

    do
    {
      for (; x<=max_x &&
           data[y][x] != border_value && data[y][x] != fill_value; x++)
      {
        data[y][x] = fill_value;
      }

      PUSH(left, x-1, y, dy);

      if (x > x2+1)
      {
        PUSH(x2+1, x-1, y, -dy);  /* leak on right? */
      }

SKIP:
      for (x++; x <= x2 &&
           (data[y][x] == border_value || data[y][x] == fill_value); x++)
      {
        ;
      }

      left = x;
    }
    while (x <= x2);
  }
}

void MyUtils::FloodFill(char* data, int x, int y,
                        int dim_x, int dim_y,
                        int fill_value, int border_value)
{
  int left, x1, x2, dy;
  LINESEGMENT stack[MAXDEPTH], *sp = stack;

  if (data[y*dim_x+x] == border_value || data[y*dim_x+x] == fill_value)
  {
    return;
  }

  int min_x = 0, max_x = dim_x-1, min_y = 0, max_y = dim_y-1;
  if (x < min_x || x > max_x || y < min_y || y > max_y)
  {
    return;
  }

  PUSH(x, x, y, 1);        /* needed in some cases */
  PUSH(x, x, y+1, -1);    /* seed segment (popped 1st) */

  while (sp > stack )
  {
    POP(x1, x2, y, dy);

    for (x = x1; x >= min_x &&
         data[y*dim_x+x] != border_value && data[y*dim_x+x] != fill_value; x--)
    {
      data[y*dim_x+x] = fill_value;
    }

    if ( x >= x1 )
    {
      goto SKIP;
    }

    left = x+1;
    if ( left < x1 )
    {
      PUSH(left, x1-1, y, -dy);  /* leak on left? */
    }

    x = x1+1;

    do
    {
      for (; x<=max_x &&
           data[y*dim_x+x] != border_value && data[y*dim_x+x] != fill_value; x++)
      {
        data[y*dim_x+x] = fill_value;
      }

      PUSH(left, x-1, y, dy);

      if (x > x2+1)
      {
        PUSH(x2+1, x-1, y, -dy);  /* leak on right? */
      }

SKIP:
      for (x++; x <= x2 &&
           (data[y*dim_x+x] == border_value || data[y*dim_x+x] == fill_value); x++)
      {
        ;
      }

      left = x;
    }
    while (x <= x2);
  }
}


template <class T> bool CalculateOptimalVolume_t( int* vox1,
                                                  int nsize1,
                                                  int* vox2,
                                                  int nsize2,
                                                  std::vector<void*> input_volumes,
                                                  T* output_volume,
                                                  int vol_size )
{
  int nvars = input_volumes.size();
  MATRIX* m1 = MatrixAlloc( nsize1, nvars, MATRIX_REAL );
  MATRIX* m2 = MatrixAlloc( nsize2, nvars, MATRIX_REAL );
  for ( int i = 0; i < nvars; i++ )
  {
    T* input_vol = (T*)input_volumes[i];
    for ( int j = 0; j < nsize1; j++ )
    {
      *MATRIX_RELT( m1, j+1, i+1 ) = input_vol[vox1[j]];
    }
    for ( int j = 0; j < nsize2; j++ )
    {
      *MATRIX_RELT( m2, j+1, i+1 ) = input_vol[vox2[j]];
    }
  }
  VECTOR* mean1 = VectorAlloc( nvars, m1->type );
  VECTOR* mean2 = VectorAlloc( nvars, m2->type );
  MATRIX* cov1 = MatrixCovariance( m1, NULL, mean1 );
  if ( cov1 == NULL )
  {
    return false;
  }
  MATRIX* cov2 = MatrixCovariance( m2, NULL, mean2 );
  if ( cov2 == NULL )
  {
    return false;
  }
  MATRIX* scov1 =
      MatrixScalarMul( cov1, (float)nsize1 / (nsize1 + nsize2 ), NULL );
  MATRIX* scov2 =
      MatrixScalarMul( cov2, (float)nsize2 / (nsize1 + nsize2 ), NULL );
  MATRIX* cov = MatrixAdd( scov1, scov2, NULL );
  MATRIX* cov_inv = MatrixInverse( cov, NULL );
  if ( cov_inv == NULL )
  {
    return false;
  }
  MATRIX* mean_sub = MatrixSubtract( mean1, mean2, NULL );
  MATRIX* weight = MatrixMultiply( cov_inv, mean_sub, NULL );
  cout << "condition number: " << MatrixConditionNumber( cov ) << "\n";
  // MATRIX* weight = MatrixCopy( mean_sub, NULL );

  double* w = new double[nvars];
  double sum = 0;
  for ( int i = 0; i < nvars; i++ )
  {
    w[i] = *MATRIX_RELT( weight, i+1, 1 );
    sum += fabs( w[i] );
  }
  cout << "Weight: ";
  for ( int i = 0; i < nvars; i++ )
  {
    w[i] /= sum;
    cout << w[i] << "\n";
  }

  double tmp = 0;
  for ( int i = 0; i < vol_size; i++ )
  {
    tmp = 0;
    for ( int j = 0; j < nvars; j++ )
    {
      tmp +=  ((T*)input_volumes[j])[i] * w[j];
    }

    output_volume[i] = (T)tmp;
  }

  MatrixFree( &m1 );
  MatrixFree( &m2 );
  MatrixFree( &mean1 );
  MatrixFree( &mean2 );
  MatrixFree( &cov1 );
  MatrixFree( &cov2 );
  MatrixFree( &scov1 );
  MatrixFree( &scov2 );
  MatrixFree( &cov );
  MatrixFree( &cov_inv );
  MatrixFree( &mean_sub );
  MatrixFree( &weight );
  delete[] w;

  return true;
}


bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1,
                                      int* vox2, int nsize2,
                                      std::vector<void*> input_volumes,
                                      float* output_volume, int vol_size )
{
  return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2,
                                   input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1,
                                      int* vox2, int nsize2,
                                      std::vector<void*> input_volumes,
                                      double* output_volume, int vol_size )
{
  return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2,
                                   input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1,
                                      int* vox2, int nsize2,
                                      std::vector<void*> input_volumes,
                                      int* output_volume, int vol_size )
{
  return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2,
                                   input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1,
                                      int* vox2, int nsize2,
                                      std::vector<void*> input_volumes,
                                      short* output_volume, int vol_size )
{
  return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2,
                                   input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1,
                                      int* vox2, int nsize2,
                                      std::vector<void*> input_volumes,
                                      unsigned char* output_volume,
                                      int vol_size )
{
  return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2,
                                   input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1,
                                      int* vox2, int nsize2,
                                      std::vector<void*> input_volumes,
                                      long* output_volume, int vol_size )
{
  return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2,
                                   input_volumes, output_volume, vol_size );
}


double MyUtils::RoundToGrid(double x)
{
  int n = 0;
  if (x == 0)
    return 0;
  int sign = 1;
  if (x < 0)
  {
    x = -x;
    sign = -1;
  }
  if (x < 1)
  {
    do
    {
      x *= 10;
      n++;
    }
    while (x < 1 || x > 10);
  }
  else if (x > 10)
  {
    do
    {
      x /= 10;
      n--;
    }
    while (x < 1 || x > 10);
  }

  if (x > 5)
  {
    x = 10;
  }
  else if (x > 2)
  {
    x = 5;
  }
  else if (x > 1)
  {
    x = 2;
  }

  return x/pow(10, n)*sign;
}

template <class T> double calculate_correlation_coefficient(float* x, T* y, int n)
{
  double sxy = 0, sx = 0, sy=0, sx2 = 0, sy2 = 0;
  for (int i = 0; i < n; i++)
  {
    sxy += x[i]*y[i];
    sx += x[i];
    sy += y[i];
    sx2 += x[i]*x[i];
    sy2 += y[i]*y[i];
  }
  double sq2 = (n*sx2-sx*sx)*(n*sy2-sy*sy);
  if (sq2 > 0)
    return (n*sxy-sx*sy)/sqrt(sq2);
  else
    return 0;
}

double MyUtils::CalculateCorrelationCoefficient(float *x, float *y, int n)
{
  return calculate_correlation_coefficient(x, y, n);
}

bool IsBetween(const double& x0, const double& x, const double& x1)
{
  return (x >= x0) && (x <= x1);
}

bool MyUtils::FindIntersection(const double& x0, const double& y0,
                      const double& x1, const double& y1,
                      const double& a0, const double& b0,
                      const double& a1, const double& b1,
                      double* x, double* y)
{
  // four endpoints are x0, y0 & x1,y1 & a0,b0 & a1,b1

  double xy, ab;
  bool partial = false;
  double denom = (b0 - b1) * (x0 - x1) - (y0 - y1) * (a0 - a1);
  if (denom == 0)
  {
    xy = -1;
    ab = -1;
  }
  else
  {
    xy = (a0 * (y1 - b1) + a1 * (b0 - y1) + x1 * (b1 - b0)) / denom;
    partial = IsBetween(0, xy, 1);
    if (partial)
    {
      // no point calculating this unless xy is between 0 & 1
      ab = (y1 * (x0 - a1) + b1 * (x1 - x0) + y0 * (a1 - x1)) / denom;
    }
  }
  if ( partial && IsBetween(0, ab, 1))
  {
    ab = 1-ab;
    xy = 1-xy;
    *x = x0 + (x1-x0)*xy;
    *y = y0 + (y1-y0)*xy;
    return true;
  }
  else
    return false;
}


bool MyUtils::FindIntersection(std::vector<std::vector<double> > &line0,
                               std::vector<std::vector<double> > &line1, double *x, double *y,
                               int* n0, int* n1)
{
  for (size_t i = 0; i < line0.size()-1; i++)
  {
    for (size_t j = 0; j < line1.size()-1; j++)
    {
      if (FindIntersection(line0[i][0], line0[i][1], line0[i+1][0], line0[i+1][1],
                           line1[j][0], line1[j][1], line1[j+1][0], line1[j+1][1], x, y))
      {
        if (n0)
          *n0 = i;
        if (n1)
          *n1 = j;
        return true;
      }
    }
  }
  return false;
}

QStringList MyUtils::SplitString( const QString& strg_to_split,
                                  const QString& divider,
                                  int nIgnoreStart,
                                  int nIgnoreLength )
{
  QStringList sa;
  QString strg = strg_to_split.trimmed();
  int n = strg.indexOf( divider );
  int nMark = n + divider.length();
  while ( n != -1 )
  {
    if ( nMark < nIgnoreStart || nMark >= nIgnoreStart + nIgnoreLength )
    {
      QString substr = strg.left( n ).trimmed();
      if ( !substr.isEmpty() )
      {
        sa << substr;
      }
      strg = strg.mid( n + divider.length() );
      n = strg.indexOf( divider );
      if ( n != -1 )
      {
        nMark += n + divider.length();
      }
    }
    else
    {
      nMark -= ( n + divider.length() );
      int nStart = 0;
      n = strg.indexOf( divider, nStart );
      while ( n != -1 &&
              (nMark + n + (int)divider.length()) >= nIgnoreStart &&
              (nMark + n + (int)divider.length()) < (nIgnoreStart + nIgnoreLength) )
      {
        nStart = n + divider.length();
        n = strg.indexOf( divider, nStart );
      }

      if ( n != -1 )
      {
        nMark += n + divider.length();
      }
    }
  }
  if ( strg.length() > 0 )
  {
    strg = strg.trimmed();
    if ( !strg.isEmpty() )
    {
      sa << strg;
    }
  }

  return sa;
}

QString MyUtils::CygwinPathToWin32Path(const QString &path_in)
{
  QString strg = path_in;
  if (strg.left(10).toLower() == "/cygdrive/")
  {
    strg = strg.mid(10).insert(1, ':');
  }
  strg.replace( '/', '\\');
  if ( strg[0] == '\\')
  {
    strg = QString("c:") + strg;
  }
  return strg;
}

QString MyUtils::Win32PathToCygwinPath(const QString &path_in)
{
  QString strg = path_in;
  if (strg[1] == ':')
  {
    strg = QString("/cygdrive/") + strg.remove(1,1);
  }
  strg.replace('\\', '/');
  return strg;
}

QString MyUtils::NormalizeCygwinPath(const QString &path_in)
{
  return Win32PathToCygwinPath(CygwinPathToWin32Path(path_in));
}

QString MyUtils::CygwinPathProof(const QString &path_in)
{
#ifdef Q_CYGWIN_WIN
  return Win32PathToCygwinPath(path_in);
#else
  return path_in;
#endif
}

QString MyUtils::Win32PathProof(const QString &path_in)
{
#ifdef Q_CYGWIN_WIN
  return CygwinPathToWin32Path(path_in);
#else
  return path_in;
#endif
}

QString MyUtils::RealToNumber(qreal val, int nPrecision)
{
  if (qAbs(val) >= pow(10, nPrecision))
    return QString("%1").arg(val, 0, 'f', 0);
  else
    return QString("%1").arg(val, 0, 'g', nPrecision);
}
