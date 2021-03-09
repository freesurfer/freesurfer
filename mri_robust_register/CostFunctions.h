/**
 * @brief A class that makes available many different cost functions for images
 *   and to combine multiple volumes by mean or median
 *   MRIiterator iterates through MRI (readonly)
 */

/*
 * Original Author: Martin Reuter
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

// written by Martin Reuter
// March 20th ,2009
//
#ifndef CostFunctions_H
#define CostFunctions_H

#include "mri.h"
#include "mriBSpline.h"
#include "matrix.h"

//#include <utility>
//#include <string>
#include <vector>
#include <iostream>

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/vnl_matrix_fixed.h>
#undef export

#include "JointHisto.h"

/** \class CostFunctions
 * \brief Static class implementing several cost functions and other image statistics
 */
class CostFunctions
{
public:
  //! Get min and max intensity value
  static std::pair<float, float> minmax(MRI *i);
  //! Get max intensity value
  static float max(MRI *i);
  //! Get mean intensity value
  static double mean(MRI *i, int frame);
  static double mean(MRI *i, const vnl_matrix_fixed<double, 4, 4>& Mi,
                     int frame,
                     int d1, int d2, int d3);
  static double mean(MRI * mri,
    const vnl_matrix_fixed<double, 4, 4>& Mi,
    int frame,
    int d1, int d2, int d3,
    int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
  //! Get variance of intensity values
  static double var(MRI *i, int frame);
  //! Get standard deviation of intensity values
  static double sdev(MRI *i, int frame)
  {
    return sqrt(var(i,frame));
  }

  //! Get median of intensity values
  static float median(MRI *i);
  //! Get d times median absolute deviation of intensity values (robust std)
  static float mad(MRI *i, float d = 1.4826);
  //! Get moment of image
  static double moment(MRI * i, int x, int y, int z);
  //! Get weighted centroid of intensity image
  static std::vector<double> centroid(MRI * i);
  //! Get principal orientation
  static vnl_matrix_fixed<double, 3, 3> orientation(MRI * i);

  //! Sum of Squared Difference (SSD)
  static double leastSquares(MRI * si, MRI * ti,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3);
  static double leastSquares(MRI * si, MRI * ti,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3,
    const double & s1, const double & s2);
  static double leastSquares(MRI_BSPLINE * sbi, MRI_BSPLINE * tbi,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3);

  //! Sum of Absolute Difference (SAD)  
  static double absDiff(MRI * si, MRI * ti,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3,
    const double & s1, const double & s2);
    
  //! never really tested
  static double tukeyBiweight(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3, double sat = 4.6851);
  //! never really tested (old)
  static double tukeyBiweight(MRI *i1, MRI * i2 = NULL, double sat = 4.6851);
  //! never really tested
  static double normalizedCorrelation(MRI * i1, MRI * i2);

  //! Normalized Correlation Coefficient
  static double NCC(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3 );


  //! Local Normalized Correlation Coefficient
  static double localNCC(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3 );


  //! Mutual Information (joint histograms)
  static double mutualInformation(MRI * i1, MRI * i2, double fwhm = 7)
  {
    JointHisto H(i1, i2);
    H.smooth(fwhm);
    return -H.computeMI();
  }

  //! Normalized Mutual Information (joint histograms)
  static double normalizedMutualInformation(MRI * i1, MRI * i2, double fwhm = 7)
  {
    JointHisto H(i1, i2);
    H.smooth(fwhm);
    return -H.computeNMI();
  }

  //! Entropy Correlation (joint histograms) (does it work?)
  static double entropyCorrelationCoefficient(MRI * i1, MRI * i2, double fwhm =
      7)
  {
    JointHisto H(i1, i2);
    H.smooth(fwhm);
    return -H.computeECC();
  }

  //! Normalized Cross Correlation (joint histograms) (does it work?)
  static double normalizedCrossCorrelation(MRI * i1, MRI * i2, double fwhm = 7)
  {
    JointHisto H(i1, i2);
    H.smooth(fwhm);
    return -H.computeNCC();
  }


  //! not implemented and not sure where they are from? Flirt?
  static float woods(MRI * i1, MRI * i2 = NULL);
  //! not implemented and not sure where they are from? Flirt?
  static float correlationRatio(MRI * i1, MRI * i2);

protected:
  inline static double rhoTukeyBiweight(double d, double sat = 4.685);
  inline static double rhoSquare(double d)
  {
    return d * d;
  }

  static MRI * mapMRI(MRI* mri, const vnl_matrix_fixed<double, 4, 4>& Mi,
    int d1, int d2, int d3);

};

inline double CostFunctions::rhoTukeyBiweight(double d, double sat)
{
  if (d > sat || d < -sat)
    return sat * sat / 2.0;
  else
  {
    double a = d / sat;
    double b = 1.0 - a * a;
    return (sat * sat / 2.0) * (1.0 - b * b * b);
  }
}

/** \class MRIiterator
 * \brief Class for iterating through an MRI
 */
class MRIiterator
{
public:
  //! Set current position to begin and initialize end pointer
  MRIiterator(MRI * i);

  //! Set position to first element
  void begin();
  //! Check wether end position is reached
  bool isEnd();
  //! Increase to next position
  MRIiterator& operator++(int);
  //! Return value at current position as float
  float operator*();

protected:

  float fromUCHAR(void);
  float fromSHORT(void);
  float fromINT(void);
  float fromLONG(void);
  float fromFLOAT(void);

  MRIiterator& opincchunk(int);
  MRIiterator& opincnochunk(int);

  MRI * img;
  unsigned char * pos;
  unsigned char * end;
  float (MRIiterator::*getVal)(void);
  MRIiterator& (MRIiterator::*opinc)(int);
  int x, y, z;
  int bytes_per_voxel;
};

inline MRIiterator::MRIiterator(MRI * i) :
    img(i)
{
  // set current pos to begin
  // and initialize end pointer
  begin();

  switch (img->type)
  {
  case MRI_UCHAR:
    getVal = &MRIiterator::fromUCHAR;
    bytes_per_voxel = sizeof(unsigned char);
    break;
  case MRI_SHORT:
    getVal = &MRIiterator::fromSHORT;
    bytes_per_voxel = sizeof(short);
    break;
  case MRI_INT:
    getVal = &MRIiterator::fromINT;
    bytes_per_voxel = sizeof(int);
    break;
  case MRI_LONG:
    getVal = &MRIiterator::fromLONG;
    bytes_per_voxel = sizeof(long);
    break;
  case MRI_FLOAT:
    getVal = &MRIiterator::fromFLOAT;
    bytes_per_voxel = sizeof(float);
    break;
  default:
    std::cerr << "MRIiterator: Type not supported: " << img->type << std::endl;
    exit(1);
  }

}

inline void MRIiterator::begin()
// set pos to first element
{
  if (img->ischunked)
  {
    pos = (unsigned char*) img->chunk;
    end = (unsigned char*) img->chunk + img->bytes_total;
    opinc = &MRIiterator::opincchunk;
  }
  else
  {
    x = 0;
    y = 0;
    z = 0;
    pos = (unsigned char*) img->slices[0][0];
    end = NULL;
    opinc = &MRIiterator::opincnochunk;
  }
}

inline bool MRIiterator::isEnd()
{
  //    if(pos > end && end != 0)
// {
//    std::cerr << "MRIiterator::isEnd outside data???" << std::endl;
//    exit(1);
// }
  return (pos == end);
}

inline MRIiterator& MRIiterator::operator++(int i)
{
  return (this->*opinc)(i);
}

inline MRIiterator& MRIiterator::opincchunk(int)
{
//   if (pos < end)
  pos += img->bytes_per_vox;
  return *this;
}

inline MRIiterator& MRIiterator::opincnochunk(int)
{
  x++;
  if (x == img->width)
  {
    x = 0;
    y++;
    if (y == img->height)
    {
      y = 0;
      z++;
      if (z == img->depth)
      {
        z = 0;
        pos = NULL;
        return *this;
      }
    }
    pos = (unsigned char*) img->slices[z][y];
  }
  else
    pos += bytes_per_voxel;
  return *this;
}

inline float MRIiterator::fromUCHAR()
{
  return ((float) *(unsigned char *) pos);
}

inline float MRIiterator::fromSHORT()
{
  return ((float) *(short *) pos);
}

inline float MRIiterator::fromINT()
{
  return ((float) *(int *) pos);
}

inline float MRIiterator::fromLONG()
{
  return ((float) *(long *) pos);
}

inline float MRIiterator::fromFLOAT()
{
  return ((float) *(float *) pos);
}

inline float MRIiterator::operator*()
{
//   if (pos < end && pos >= img->chunk)
  return (this->*getVal)();
}

// example:
// MRIiterator it(mri);
// for (it.begin(); !it.isEnd(); it++)
// {
//    std::cout << *it << std::endl;
////    *it = 0;
// }

#endif
