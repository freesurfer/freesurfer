/**
 * @file  generators.hpp
 * @brief Holds generators for GPU routines
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef GENERATORS_CUDA_H
#define GENERATORS_CUDA_H


//! Class to generate linear progressions
class LinearGenerator
{
public:
  //! Starting value for this generator
  float minVal;
  //! Maximum value for this generator
  float maxVal;
  //! Number of steps for this generator
  unsigned int n;

  //! Error value
  static const float errVal = 1e10;

  //! Operator to access the ith value
  __host__ __device__ float operator() ( const unsigned int i ) const
  {
    float res = errVal;

    if( i < this->n )
    {
      res = static_cast<float>(i) / ( this->n - 1 );
      res *= ( this->maxVal - this->minVal );
      res += this->minVal;
    }

    return( res );
  }

  //! Constructor with inputs
  __host__ __device__ LinearGenerator( const float myMin,
                   const float myMax,
                   const unsigned int myn ) : minVal(myMin),
    maxVal(myMax),
    n(myn) {};

};



class TranslationGenerator
{
public:

  LinearGenerator xGen, yGen, zGen;

  //! Operator to return the ith translation vector
  __host__ __device__ float3 operator() ( unsigned int targetIndex ) const
  {
    float3 res;

    // Sanity check
    if( targetIndex >= ( this->xGen.n * this->yGen.n * this->zGen.n ) )
    {
      res.x = res.y = res.z = LinearGenerator::errVal;
      return( res );
    }

    res.z = this->zGen( targetIndex % this->zGen.n );
    targetIndex /= this->zGen.n;
    res.y = this->yGen( targetIndex % this->yGen.n );
    targetIndex /= this->yGen.n;
    res.x = this->xGen( targetIndex );

    return( res );
  }

  //! Constructor with identical limits for all
  __host__ __device__ TranslationGenerator( const float myMin,
                        const float myMax,
                        const unsigned int myn ) : xGen(myMin,myMax,myn),
    yGen(myMin,myMax,myn),
    zGen(myMin,myMax,myn) {}
};



class TransformGenerator
{
public:

  LinearGenerator xTrans, yTrans, zTrans;
  LinearGenerator xScale, yScale, zScale;
  LinearGenerator xRot, yRot, zRot;


  //! Method to return the ith transform parameters
  __host__ __device__ void GetTransform( unsigned int targetIndex,
                                         float3& trans, float3& scale, float3& rot ) const
  {

    // Sanity check
    long tot = ( this->xTrans.n * this->yTrans.n * this->zTrans.n );
    tot *= ( this->xScale.n * this->yScale.n * this->zScale.n );
    tot *= ( this->xRot.n * this->yRot.n * this->zRot.n );
    if( targetIndex >= tot )
    {
      trans.x = trans.y = trans.z = LinearGenerator::errVal;
      scale.x = scale.y = scale.z = LinearGenerator::errVal;
      rot.x = rot.y = rot.z = LinearGenerator::errVal;
      return;
    }

    // Now go through, setting things
    trans.z = this->zTrans( targetIndex % this->zTrans.n );
    targetIndex /= this->zTrans.n;
    trans.y = this->yTrans( targetIndex % this->yTrans.n );
    targetIndex /= this->yTrans.n;
    trans.x = this->xTrans( targetIndex % this->xTrans.n );
    targetIndex /= this->xTrans.n;

    rot.z = this->zRot( targetIndex % this->zRot.n );
    targetIndex /= this->zRot.n;
    rot.y = this->yRot( targetIndex % this->yRot.n );
    targetIndex /= this->yRot.n;
    rot.x = this->xRot( targetIndex % this->xRot.n );
    targetIndex /= this->xRot.n;

    scale.z = this->zScale( targetIndex % this->zScale.n );
    targetIndex /= this->zScale.n;
    scale.y = this->yScale( targetIndex % this->yScale.n );
    targetIndex /= this->yScale.n;
    scale.x = this->xScale( targetIndex % this->xScale.n );
    targetIndex /= this->xScale.n;


  }


  //! Constructor with generators for translations, scalings and rotations
  __host__ __device__ TransformGenerator( const LinearGenerator genTrans,
                      const LinearGenerator genScale,
                      const LinearGenerator genRot ) : xTrans(genTrans),
    yTrans(genTrans),
    zTrans(genTrans),
    xScale(genScale),
    yScale(genScale),
    zScale(genScale),
    xRot(genRot),
    yRot(genRot),
    zRot(genRot) {}
};


#endif
