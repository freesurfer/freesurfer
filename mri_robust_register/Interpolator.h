/**
 * @brief Classes for different interpolation methods
 *
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
//
// written by Martin Reuter
// Sep. 4th ,2013
//
#ifndef Interpolator_H
#define Interpolator_H

#include "mri.h"
#include "mriBSpline.h"

/** \class Interpolator
 * \brief Base class for different interpolation methods
 */
class Interpolator
{
public:
  virtual ~Interpolator()
  {
  }
  ;
  //! Sample image at specific location and frame
  virtual double sample(double x, double y, double z, int frame) const =0;
  //! Resample full image (in mri_dst space), after applying voxel transform
  virtual MRI * transformVOX(MRI *mri_dst, MATRIX *mA) const=0;
  //! Resample full image (in mri_dst space), after applying voxel transform
  virtual MRI * transformVOX(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const=0;
  //! Resample full image (in mri_dst space), after applying RAS transform
  virtual MRI * transformRAS(MRI *mri_dst, MATRIX *mA) const =0;
  //! Resample full image (in mri_dst space), after applying RAS transform
  virtual MRI * transformRAS(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const =0;
  //! Resample full image (in mri_dst space), after applying LTA
  virtual MRI * transformLTA(MRI *mri_dst, LTA *lta) const =0;
  //! Get the degree (0=nearest, 1=linear, higher = BSpline)
  inline virtual int getDegree() const
  {
    return degree;
  }
  ;
protected:
  //! The degree (0=nearest, 1=linear, higher = BSpline)
  int degree;
};

/** \class InterpolatorNearest
 * \brief Class for nearest neighbor interpolation
 */
class InterpolatorNearest: public Interpolator
{
public:
  //! Constructor creates a copy of the image
  InterpolatorNearest(MRI* tmri) :
      degree(0)
  {
    mri = MRIcopy(tmri, NULL);
  }
  ;
  //! Destructor frees copy of the image
  ~InterpolatorNearest()
  {
    MRIfree(&mri);
  }
  ;

  inline virtual double sample(double x, double y, double z, int frame) const
  {
    double val;
    return MRIsampleVolumeFrameType(mri, x, y, z, frame, SAMPLE_NEAREST, &val);
    return val;
  }
  ;
  inline virtual MRI * transformVOX(MRI *mri_dst, MATRIX *mA) const
  {
    return MRIlinearTransformInterp(mri, mri_dst, mA, SAMPLE_NEAREST);
  }
  ;
  inline virtual MRI * transformVOX(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const
  {
    MATRIX * mA = MyMatrix::convertVNL2MATRIX(m, NULL);
    MRI* mriT = ::MRIlinearTransformInterp(mri, mri_dst, mA, SAMPLE_NEAREST);
    MatrixFree(&mA);
    return mriT;
  }
  ;
  inline virtual MRI * transformRAS(MRI *mri_dst, MATRIX *mA) const
  {
    return MRIapplyRASlinearTransformInterp(mri, mri_dst, mA, SAMPLE_NEAREST);
  }
  ;
  inline virtual MRI * transformRAS(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const
  {
    MATRIX * mA = MyMatrix::convertVNL2MATRIX(m, NULL);
    MRI* mriT = ::MRIapplyRASlinearTransformInterp(mri, mri_dst, mA,
        SAMPLE_NEAREST);
    MatrixFree(&mA);
    return mriT;
  }
  ;
  inline virtual MRI * transformLTA(MRI *mri_dst, LTA *lta) const
  {
    return LTAtransformInterp(mri, mri_dst, lta, SAMPLE_NEAREST);
  }
  ;

private:
  MRI * mri;
};

/** \class InterpolatorLinear
 * \brief Class for tri-linear interpolation
 */
class InterpolatorLinear: public Interpolator
{
public:
  //! Constructor creates a copy of the image
  InterpolatorLinear(MRI* tmri) :
      degree(1)
  {
    mri = MRIcopy(tmri, NULL);
  }
  ;
  //! Destructor frees copy of the image
  ~InterpolatorLinear()
  {
    MRIfree(&mri);
  }
  ;

  inline virtual double sample(double x, double y, double z, int frame) const
  {
    double val;
    return MRIsampleVolumeFrameType(mri, x, y, z, frame, SAMPLE_TRILINEAR, &val);
    return val;
  }
  ;
  inline virtual MRI * transformVOX(MRI *mri_dst, MATRIX *mA) const
  {
    return MRIlinearTransformInterp(mri, mri_dst, mA, SAMPLE_TRILINEAR);
  }
  ;
  inline virtual MRI * transformVOX(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const
  {
    MATRIX * mA = MyMatrix::convertVNL2MATRIX(m, NULL);
    MRI* mriT = ::MRIlinearTransformInterp(mri, mri_dst, mA, SAMPLE_TRILINEAR);
    MatrixFree(&mA);
    return mriT;
  }
  ;
  inline virtual MRI * transformRAS(MRI *mri_dst, MATRIX *mA) const
  {
    return MRIapplyRASlinearTransformInterp(mri, mri_dst, mA, SAMPLE_TRILINEAR);
  }
  ;
  inline virtual MRI * transformRAS(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const
  {
    MATRIX * mA = MyMatrix::convertVNL2MATRIX(m, NULL);
    MRI* mriT = ::MRIapplyRASlinearTransformInterp(mri, mri_dst, mA,
        SAMPLE_TRILINEAR);
    MatrixFree(&mA);
    return mriT;
  }
  ;
  inline virtual MRI * transformLTA(MRI *mri_dst, LTA *lta) const
  {
    return LTAtransformInterp(mri, mri_dst, lta, SAMPLE_TRILINEAR);
  }
  ;

private:
  MRI * mri;
};

/** \class InterpolatorBSpline
 * \brief Class for BSpline interpolation
 */
class InterpolatorBSpline: public Interpolator
{
public:
  //! Constructor creates bspline coefficient image
  InterpolatorBSpline(MRI* tmri, int tdegree) :
      degree(tdegree)
  {
    bspline = MRItoBSpline(tmri, NULL, tdegree);
  }
  ;
  //! Destructor frees coefficient image
  virtual ~InterpolatorBSpline()
  {
    MRIfreeBSpline(&bspline);
  }
  ;

  inline virtual double sample(double x, double y, double z, int frame) const
  {
    double val;
    return MRIsampleBSpline(bspline, x, y, z, frame, &val);
    return val;
  }
  ;
  inline virtual MRI * transformVOX(MRI *mri_dst, MATRIX *mA) const
  {
    return MRIlinearTransformBSpline(bspline, mri_dst, mA);
  }
  ;
  inline virtual MRI * transformVOX(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const
  {
    MATRIX * mA = MyMatrix::convertVNL2MATRIX(m, NULL);
    MRI* mriT = ::MRIlinearTransformBSpline(bspline, mri_dst, mA);
    MatrixFree(&mA);
    return mriT;
  }
  ;
  inline virtual MRI * transformRAS(MRI *mri_dst, MATRIX *mA) const
  {
    return MRIapplyRASlinearTransformBSpline(bspline, mri_dst, mA);
  }
  ;
  inline virtual MRI * transformRAS(MRI *mri_dst,
      const vnl_matrix_fixed<double, 4, 4>& m) const
  {
    MATRIX * mA = MyMatrix::convertVNL2MATRIX(m, NULL);
    MRI* mriT = ::MRIapplyRASlinearTransformBSpline(bspline, mri_dst, mA);
    MatrixFree(&mA);
    return mriT;
  }
  ;
  inline virtual MRI * transformLTA(MRI *mri_dst, LTA *lta) const
  {
    return LTAtransformBSpline(bspline, mri_dst, lta);
  }
  ;
//  inline virtual MRI * downsample(MRI *mri_dst) const {return MRIdownsample2BSpline(bspline,mri_dst);};
//  inline virtual MRI * upsample(MRI *mri_dst) const {return MRIupsample2BSpline(bspline,mri_dst);};

private:
  MRI_BSPLINE * bspline;
};

/** \class InterpolatorCubic
 * \brief Class for cubic BSpline interpolation
 */
class InterpolatorCubic: public InterpolatorBSpline
{
public:
  //! Constructor creates bspline coefficient image
  InterpolatorCubic(MRI* tmri) :
      InterpolatorBSpline(tmri, 3)
  {
  }
  ;

};

#endif
