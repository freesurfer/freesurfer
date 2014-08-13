//TensorCubicSmoothing.h
//
//Created 7/17/2014
//By: Benjamin Lewin
//

#ifndef TENSOR_CUBIC_SMOOTHING_INCLUDED
#define TENSOR_CUBIC_SMOOTHING_INCLUDED

#include "matrix3d.h"
#include <vnl/vnl_vector.h>
//namespace ctl 
//{

  class TensorCubicSmoothing
  {
    public:

      TensorCubicSmoothing();
      //TensorCubicSmoothing(const Matrix &Bx, const Matrix &By,
      //  const Matrix &Bz, const Matrix &W, const Matrix &indexMap);
      ~TensorCubicSmoothing();

      //Copy Constructor and Assignment overload
      TensorCubicSmoothing(const TensorCubicSmoothing& other);
      TensorCubicSmoothing& operator=(const TensorCubicSmoothing& other); 
      
      /// takes in separable basis function in three dimensions and weights
      /// and constructs the matrix AtWA by exploiting the seperability
      int constructAtWA(const vnl_matrix<float> &B2x, const vnl_matrix<float> &B2y,
                        const vnl_matrix<float> &B2z, const Matrix3d &W,
                        const vnl_vector<int> &indexMap);

      /// takes in separable basis function in three dimensions and weights
      /// and constructs the matrix AtWr by exploiting the seperability
      int constructAtWr(const vnl_matrix<float> &Bx, const vnl_matrix<float> &By,
                        const vnl_matrix<float> &Bz, const Matrix3d &W,
                        const Matrix3d &r);
      
      /// solves the non-linear system by means of least squares given regularization parameters
      int solve(const Matrix3d &P, float lambda);
      /// solves the non-linear system by means of least squares without regularization
      int solve();

      /// returns the solved coefficients
      void getCoefficients(vnl_vector<float> &c) const;
      /// returns the smoothed data
      void expandCoefficients(Matrix3d &d, const vnl_matrix<float> &Bx,
                                           const vnl_matrix<float> &By,
                                           const vnl_matrix<float> &Bz) const;

      // helper function for  constructAtWA and construct AtWr
      // TODO NEEDS A NEW NAME
      // will be private once testing is done
      Matrix3d doSneakySeperabilityMultiplications(const Matrix3d &data,
                                                   const vnl_matrix<float> &Bx,
                                                   const vnl_matrix<float> &By,
                                                   const vnl_matrix<float> &Bz);
    private:

      vnl_matrix<float> AtWA;
      /// 1d  vector
      vnl_vector<float> AtWr;

      /// coefficients  (vector)
      vnl_vector<float> coefficients;

  };

//};

#endif


