/**
 * @brief A base class for different transformation models used for registration
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
// Aug. 28th ,2013
//
#ifndef Transformation_H
#define Transformation_H

#include <vector>
#include <vnl/vnl_matrix_fixed.h>

#include "Quaternion.h"
#include "MyMatrix.h"

/** \class Transformation
 * \brief Base class for all 2d and 3d transformation models
 */
class Transformation
{
public:
  virtual ~Transformation()
  {}
  
  //! Get the transformation as a 4x4 affine matrix
  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const =0;
  //! Get the degrees of freedom of the transform
  virtual unsigned int getDOF() const
  {
    return parameters.size();
  }

  //! Get the parameters
  const vnl_vector<double> & getParameters() const
  {
    return parameters;
  }

  //! Get steps for Powell
  virtual vnl_vector<double> getSteps() const =0;
  //! Get the gradient ( grad Image * grad Transform )
  virtual vnl_vector<double> getGradient(const unsigned int& x, const float& fx,
      const unsigned int& y, const float& fy, const unsigned int& z,
      const float& fz) const =0;

  //! Set the parameters from double std vector
  void setParameters(const std::vector<double> &p)
  {
    assert(p.size() >= getDOF());
    assert(parameters.size() == getDOF());
    for (unsigned int i = 0; i < getDOF(); i++)
      parameters[i] = p[i];
  }

  //! Set the parameters from float vnl vector
  void setParameters(const vnl_vector<float> &v)
  {
    assert(v.size() >= getDOF());
    assert(parameters.size() == getDOF());
    for (unsigned int i = 0; i < getDOF(); i++)
      parameters[i] = v[i];
  }

  //! Set the parameters from double vnl vector
  void setParameters(const vnl_vector<double> &v)
  {
    assert(v.size() >= getDOF());
    assert(parameters.size() == getDOF());
    parameters = v;
    parameters.set_size(getDOF());
  }

  //! Set parameters to identity transform
  virtual void setIdentity()
  {
    parameters.fill(0.0);
  }

protected:
  //! vector to store the parameters
  vnl_vector<double> parameters;
  //std::vector < double > parameters;
};

/****************************************************************************************************

 2D transforms
 
 affine, isoscale (rigid+global scaling), rigid (translation+rotation), translation only
 
 Note the transforms ending with 2 are slightly different and are mainly used for the powell method.
 
 *****************************************************************************************************/

/** \class Transform2dAffine
 * \brief Describes affine 2d transformation (6 DOF) as 2x3 matrix added to the 4x4 identity matrix
 */
class Transform2dAffine: public Transformation
{
public:
  Transform2dAffine()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform2dAffine(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform2dAffine(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform2dAffine(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform2dAffine()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 6;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[2] = 0.02;
    v[5] = 0.02;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    ret.set_identity();
    ret[0][0] += parameters[0];
    ret[0][1] = parameters[1];
    ret[0][3] = parameters[2];
    ret[1][0] = parameters[3];
    ret[1][1] += parameters[4];
    ret[1][3] = parameters[5];
    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret(6);
    ret[0] = (double) x * fx;
    ret[1] = (double) y * fx;
    ret[2] = fx;
    ret[3] = (double) x * fy;
    ret[4] = (double) y * fy;
    ret[5] = fy;
    return ret;
  }

};

/** \class Transform2dAffine2
 * \brief Describes affine 2d transformation (6 DOF) as t_x,t_y,alpha_z,scale_x,scale_y,shear
 * 
 * where alpha_z is negative euler angle around z axis (clockwise, same as spm)
 */
class Transform2dAffine2: public Transformation
{ // M = T*shear*Scale*Rot
public:
  Transform2dAffine2()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform2dAffine2(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform2dAffine2(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform2dAffine2(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform2dAffine2()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 6;
  }

  virtual void setIdentity()
  {
    std::fill(parameters.begin(), parameters.end(), 0.0);
    parameters[3] = 1.0;
    parameters[4] = 1.0;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    v[3] = 0.01;
    v[4] = 0.01;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    // M = T*shear*Scale*Rot
    vnl_matrix_fixed<double, 4, 4> ret;
    // Translation
    vnl_vector_fixed<double, 3> t(parameters[0], parameters[1], 0);
    //Rotation
    Quaternion q;
    q.importZYXAngles(-parameters[2], 0, 0);
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);
    //Scale
    vnl_matrix<double> smat(3, 3, 0.0);
    smat[0][0] = parameters[3];
    smat[1][1] = parameters[4];
    smat[2][2] = 1;
    //Shear
    vnl_matrix<double> zmat(3, 3);
    zmat.set_identity();
    zmat[0][1] = parameters[5];
    // product 3x3
    vnl_matrix<double> M3 = zmat * smat * rmat;
    // consturct 4x4 with translation also:
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
    {
      for (cc = 0; cc < 3; cc++) // copy M3
        ret[rr][cc] = M3[rr][cc];

      // copy translation into 4th column
      ret[rr][3] = t[rr];
      // set 4th row to zero
      ret[3][rr] = 0.0;
    }
    //except 4,4
    ret[3][3] = 1.0;

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    std::cerr << "Gradient for affine 2d (second type) not implemented!"
        << std::endl;
    exit(1);
    vnl_vector<double> ret(6);
    ret[0] = fx;
    ret[1] = fx;
    ret[2] = (-y * fx) + (x * fy); // or negative rotation?
    ret[3] = x * fx;
    ret[4] = y * fy;
    ret[5] = y * fx;
    return ret;
  }

};

/** \class Transform2dIsoscale
 * \brief Describes rigid 2d transformation with isotropic scaling (4 DOF) [ p -q ; q p ] + T
 * 
 */
class Transform2dIsoscale: public Transformation
{
public:
  Transform2dIsoscale()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform2dIsoscale(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform2dIsoscale(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform2dIsoscale(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform2dIsoscale()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 4;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    ret.set_identity();
    ret[0][0] += parameters[2];
    ret[0][1] = -parameters[3];
    ret[0][3] = parameters[0];
    ret[1][0] = parameters[3];
    ret[1][1] += parameters[2];
    ret[1][3] = parameters[1];
    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret(4);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = fx * x + fy * y;
    ret[3] = -fx * y + fy * x;
    return ret;
  }

};

/** \class Transform2dIsoscale2
 * \brief Describes rigid 2d transformation with isotropic scaling (4 DOF) tx,ty,alpha,s
 * 
 */
class Transform2dIsoscale2: public Transformation
{ // rigid and isotropic scaling 
public:
  Transform2dIsoscale2()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform2dIsoscale2(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform2dIsoscale2(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform2dIsoscale2(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform2dIsoscale2()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 4;
  }

  inline virtual void setIdentity()
  {
    parameters.fill(0.0);
    parameters[3] = 1.0;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    v[3] = 0.01;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    // M = T*(Scale*Rot)
    //Rot
    Quaternion q;
    q.importZYXAngles(-parameters[2], 0, 0);
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);

    // scale
    ret.set_identity();
    ret[0][0] = parameters[3] * rmat[0][0];
    ret[0][1] = parameters[3] * rmat[0][1];
    ret[1][0] = parameters[3] * rmat[1][0];
    ret[1][1] = parameters[3] * rmat[1][1];

    // translation
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    std::cerr << "Gradient for isoscale 2d (second type) not implemented!"
        << std::endl;
    exit(1);
    vnl_vector<double> ret(4);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = (-y * fx) + (x * fy); // or negative rotation?
    ret[3] = fx * x + fy * y;
    return ret;
  }

};

/** \class Transform2dRigid
 * \brief Describes rigid 2d transformation (3 DOF) tx,ty,alpha (pos around z)
 * 
 */
class Transform2dRigid: public Transformation
{ // rigid  tx, ty, r
public:
  Transform2dRigid()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform2dRigid(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform2dRigid(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform2dRigid(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform2dRigid()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 3;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    Quaternion q;
    // r1,r2,r3 are as in robust paper (axis, and length is angle)
    q.importRotVec(0, 0, parameters[2]);
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);
    vnl_matrix_fixed<double, 4, 4> ret;
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
    {
      for (cc = 0; cc < 3; cc++) // copy rot-matrix
        ret[rr][cc] = rmat[rr][cc];

      // set 4th row and col to zero
      ret[3][rr] = 0.0;
      ret[rr][3] = 0.0;
    }
    //except 4,4
    ret[3][3] = 1.0;
    // translation
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];
    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret(3);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = (fy * x - fx * y);
    return ret;
  }

};

/** \class Transform2dRigid2
 * \brief Describes rigid 2d transformation (3 DOF) tx,ty,alpha (clockwise around z)
 */
class Transform2dRigid2: public Transformation
{ // rigid  tx, ty, r  where r is rot around z axis (clockwise?, same as SPM)
public:
  Transform2dRigid2()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform2dRigid2(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform2dRigid2(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform2dRigid2(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform2dRigid2()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 3;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    // rigid: first 2 translation, next rotation (as a vector)
    // split translation and rotation:
    vnl_vector_fixed<double, 3> t(parameters[0], parameters[1], 0);
    double r = parameters[2];

    // converts rot vector (3x1) and translation vector (3x1)
    // into an affine matrix (homogeneous coord) 4x4
    // if global rtype ==2 then r1,r2,r3 are angles around x,y,z axis (order 1zrot,2yrot,3xrot)
    vnl_matrix<double> rmat;
    Quaternion q;
    // first convert rotation to quaternion (clockwise)
    //q.importZYXAngles(-r[2], -r[1], -r[0]);
    q.importZYXAngles(-r, 0, 0); // same as spm now
    // then to rotation matrix
    rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);

    vnl_matrix_fixed<double, 4, 4> ret;
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
    {
      for (cc = 0; cc < 3; cc++) // copy rot-matrix
        ret[rr][cc] = rmat[rr][cc];

      // copy translation into 4th column
      ret[rr][3] = t[rr];
      // set 4th row to zero
      ret[3][rr] = 0.0;
    }
    //except 4,4
    ret[3][3] = 1.0;

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    std::cerr << "Gradient for rigid 2d (second type) not implemented!"
        << std::endl;
    exit(1);
    vnl_vector<double> ret(3);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = (y * fx) - (x * fy); // negative rotation compared to other 2d rigid
    return ret;
  }

};

/** \class Transform2dTranslate
 * \brief Describes translation in 2d (2 DOF) tx,ty
 */
class Transform2dTranslate: public Transformation
{
public:
  Transform2dTranslate()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform2dTranslate(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform2dTranslate(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform2dTranslate(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform2dTranslate()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 2;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.02);
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    ret.set_identity();
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];
    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret(2);
    ret[0] = fx;
    ret[1] = fy;
    return ret;
  }

};

/**************************************************************

 3D transforms
 
 ***************************************************************/

/** \class Transform3dAffine
 * \brief Describes affine transformation in 3d (12 DOF) matrix add-ons
 */
class Transform3dAffine: public Transformation
{
public:
  Transform3dAffine()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform3dAffine(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform3dAffine(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform3dAffine(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform3dAffine()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 12;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[3] = 0.02;
    v[7] = 0.02;
    v[11] = 0.02;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    ret.set_identity();
    int count = 0;
    for (int rr = 0; rr < 3; rr++)
      for (int cc = 0; cc < 4; cc++)
      {
        ret[rr][cc] += parameters[count];
        count++;
      }
    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret(12);
    ret[0] = fx * x;
    ret[1] = fx * y;
    ret[2] = fx * z;
    ret[3] = fx;
    ret[4] = fy * x;
    ret[5] = fy * y;
    ret[6] = fy * z;
    ret[7] = fy;
    ret[8] = fz * x;
    ret[9] = fz * y;
    ret[10] = fz * z;
    ret[11] = fz;
    return ret;
  }

};

/** \class Transform3dAffine2
 * \brief Describes affine transformation in 3d (12 DOF) M = T*shear*Scale*Rot
 */
class Transform3dAffine2: public Transformation
{
public:
  Transform3dAffine2()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform3dAffine2(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform3dAffine2(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform3dAffine2(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform3dAffine2()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 12;
  }

  inline virtual void setIdentity()
  {
    parameters.fill(0.0);
    parameters[6] = 1.0;
    parameters[7] = 1.0;
    parameters[8] = 1.0;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    v[2] = 0.02;
    v[6] = 0.01;
    v[7] = 0.01;
    v[8] = 0.01;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    // M = T*shear*Scale*Rot

    //Rot
    Quaternion q;
    q.importZYXAngles(-parameters[5], parameters[4], -parameters[3]); // same as spm now
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);
    //Scale
    vnl_matrix<double> smat(3, 3, 0.0);
    smat[0][0] = parameters[6];
    smat[1][1] = parameters[7];
    smat[2][2] = parameters[8];
    //Shear
    vnl_matrix<double> zmat(3, 3);
    zmat.set_identity();
    zmat[0][1] = parameters[9];
    zmat[0][2] = parameters[10];
    zmat[1][2] = parameters[11];
    // product 3x3
    vnl_matrix<double> M3 = zmat * smat * rmat;
    // consturct 4x4 with translation also:
    vnl_matrix_fixed<double, 4, 4> ret;
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
    {
      for (cc = 0; cc < 3; cc++) // copy M3
        ret[rr][cc] = M3[rr][cc];

      // copy translation into 4th column
      ret[rr][3] = parameters[rr];
      // set 4th row to zero
      ret[3][rr] = 0.0;
    }
    //except 4,4
    ret[3][3] = 1.0;

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    std::cerr << " Affine in 3D (type 2): gradient not implemented!"
        << std::endl;
    exit(1);
    vnl_vector<double> ret(12);
    ret[0] = fx * x;
    ret[1] = fx * y;
    ret[2] = fx * z;
    ret[3] = fx;
    ret[4] = fy * x;
    ret[5] = fy * y;
    ret[6] = fy * z;
    ret[7] = fy;
    ret[8] = fz * x;
    ret[9] = fz * y;
    ret[10] = fz * z;
    ret[11] = fz;
    return ret;
  }

};

/** \class Transform3dIsoscale
 * \brief Describes rigid and isoscale transformation in 3d (7 DOF)  M = T*(Scale*Rot)
 */
class Transform3dIsoscale: public Transformation
{
public:
  Transform3dIsoscale()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform3dIsoscale(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform3dIsoscale(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform3dIsoscale(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform3dIsoscale()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 7;
  }

  inline virtual void setIdentity()
  {
    parameters.fill(0.0);
    parameters[6] = 1.0;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    v[2] = 0.02;
    v[6] = 0.01;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    std::cerr << " Isoscale in 3D not implemented yet, use ridig or affine"
        << std::endl;
    exit(1);
    vnl_matrix_fixed<double, 4, 4> ret;
    // M = T*(Scale*Rot)

    //Rotation
    Quaternion q;
    q.importRotVec(parameters[3], parameters[4], parameters[5]);
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);

    // scale
    rmat = ((double) parameters[6]) * rmat;

    // copy 
    ret.set_identity();
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
      for (cc = 0; cc < 3; cc++) // copy 
        ret[rr][cc] = rmat[rr][cc];

    // translation
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];
    ret[2][3] = parameters[2];

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    std::cerr << " Isoscale in 3D not implemented yet, use ridig or affine"
        << std::endl;
    exit(1);
    vnl_vector<double> ret(7);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = fz;
    ret[3] = (fz * y - fy * z);
    ret[4] = (fx * z - fz * x);
    ret[5] = (fy * x - fx * y);
    ret[6] = (fx * x + fy * y);
    return ret;
  }

};

/** \class Transform3dIsoscale2
 * \brief Describes rigid and isoscale transformation in 3d (7 DOF) M = T*(Scale*Rot)
 */
class Transform3dIsoscale2: public Transformation
{
public:
  Transform3dIsoscale2()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform3dIsoscale2(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform3dIsoscale2(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform3dIsoscale2(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform3dIsoscale2()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 7;
  }

  inline virtual void setIdentity()
  {
    parameters.fill(0.0);
    parameters[6] = 1.0;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    v[2] = 0.02;
    v[6] = 0.01;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    // M = T*(Scale*Rot)

    //Rotation
    Quaternion q;
    q.importZYXAngles(-parameters[5], parameters[4], -parameters[3]); // same as spm now
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);

    // scale
    rmat = ((double) parameters[6]) * rmat;

    // copy 
    ret.set_identity();
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
      for (cc = 0; cc < 3; cc++) // copy 
        ret[rr][cc] = rmat[rr][cc];

    // translation
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];
    ret[2][3] = parameters[2];

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    std::cerr << " Isoscale in 3D (type 2): gradient not implemented!"
        << std::endl;
    exit(1);
    vnl_vector<double> ret(7);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = fz;
    ret[3] = (fz * y - fy * z);
    ret[4] = (fx * z - fz * x);
    ret[5] = (fy * x - fx * y);
    ret[6] = (fx * x + fy * y);
    return ret;
  }

};

/** \class Transform3dRigid
 * \brief Describes rigid transformation in 3d (6 DOF) 
 */
class Transform3dRigid: public Transformation
{
public:
  Transform3dRigid()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform3dRigid(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform3dRigid(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform3dRigid(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform3dRigid()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 6;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    v[2] = 0.02;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    // M = T*(Scale*Rot)

    //Rotation
    Quaternion q;
    q.importRotVec(parameters[3], parameters[4], parameters[5]);
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);

    // copy 
    ret.set_identity();
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
      for (cc = 0; cc < 3; cc++) // copy 
        ret[rr][cc] = rmat[rr][cc];

    // translation
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];
    ret[2][3] = parameters[2];

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret(6);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = fz;
    ret[3] = (fz * y - fy * z);
    ret[4] = (fx * z - fz * x);
    ret[5] = (fy * x - fx * y);
    return ret;
  }

};

/** \class Transform3dRigid2
 * \brief Describes rigid transformation in 3d (6 DOF) 
 */
class Transform3dRigid2: public Transformation
{
public:
  Transform3dRigid2()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform3dRigid2(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform3dRigid2(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform3dRigid2(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform3dRigid2()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 6;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.001);
    v[0] = 0.02;
    v[1] = 0.02;
    v[2] = 0.02;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    // M = T*(Scale*Rot)

    //Rotation
    Quaternion q;
    q.importZYXAngles(-parameters[5], parameters[4], -parameters[3]); // same as spm now
    vnl_matrix<double> rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(), 3);

    // copy 
    ret.set_identity();
    int rr, cc;
    for (rr = 0; rr < 3; rr++)
      for (cc = 0; cc < 3; cc++) // copy 
        ret[rr][cc] = rmat[rr][cc];

    // translation
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];
    ret[2][3] = parameters[2];

    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    std::cerr << "ERROR rigid in 3D (type 2): gradient not implemented !"
        << std::endl;
    exit(1);
    vnl_vector<double> ret(6);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = fz;
    ret[3] = (fz * y - fy * z);
    ret[4] = (fx * z - fz * x);
    ret[5] = (fy * x - fx * y);
    return ret;
  }

};

/** \class Transform3dTranslate
 * \brief Describes translation in 3d (3 DOF) 
 */
class Transform3dTranslate: public Transformation
{
public:
  Transform3dTranslate()
  {
    parameters.set_size(getDOF());
    setIdentity();
  }

  Transform3dTranslate(const std::vector<double> &p)
  {
    parameters.set_size(getDOF());
    setParameters(p);
  }

  Transform3dTranslate(const vnl_vector<float> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  Transform3dTranslate(const vnl_vector<double> &v)
  {
    parameters.set_size(getDOF());
    setParameters(v);
  }

  virtual ~Transform3dTranslate()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 3;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v(getDOF(), 0.02);
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    ret.set_identity();
    // translation
    ret[0][3] = parameters[0];
    ret[1][3] = parameters[1];
    ret[2][3] = parameters[2];
    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret(3);
    ret[0] = fx;
    ret[1] = fy;
    ret[2] = fz;
    return ret;
  }

};

/** \class TransformIdentity
 * \brief Describes Identity (no geometric transform) (0 DOF) 
 */
class Transform3dIdentity: public Transformation
{
public:
  Transform3dIdentity()
  {
  }

  Transform3dIdentity(const std::vector<double> &p)
  {
  }

  Transform3dIdentity(const vnl_vector<float> &v)
  {
  }

  Transform3dIdentity(const vnl_vector<double> &v)
  {
  }

  virtual ~Transform3dIdentity()
  {}
  
  inline virtual unsigned int getDOF() const
  {
    return 0;
  }

  inline virtual vnl_vector<double> getSteps() const
  {
    vnl_vector<double> v;
    return v;
  }

  virtual vnl_matrix_fixed<double, 4, 4> getMatrix() const
  {
    vnl_matrix_fixed<double, 4, 4> ret;
    ret.set_identity();
    return ret;
  }

  inline virtual vnl_vector<double> getGradient(const unsigned int& x,
      const float& fx, const unsigned int& y, const float& fy,
      const unsigned int& z, const float& fz) const
  {
    vnl_vector<double> ret;
    return ret;
  }

};

#endif
