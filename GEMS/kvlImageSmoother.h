/**
 * @file  kvlImageSmoother.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#ifndef __kvlImageSmoother_h
#define __kvlImageSmoother_h

#include "itkImage.h"


namespace kvl
{


class ImageSmoother: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef ImageSmoother  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ImageSmoother, itk::Object );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;
  typedef itk::Image< bool, 3 >  MaskImageType;

  //
  void  SetMaskImage( const MaskImageType* maskImage )
  {
    m_MaskImage = maskImage;
    m_BasisFunctions = vnl_matrix< double >();
  }

  const MaskImageType*  GetMaskImage() const
  {
    return m_MaskImage;
  }

  //
  void  SetPolynomialOrder( int polynomialOrder )
  {
    m_PolynomialOrder = polynomialOrder;
    m_PolynomialOrderUsed = polynomialOrder;
    m_BasisFunctions = vnl_matrix< double >();
  }

  int  GetPolynomialOrder() const
  {
    return m_PolynomialOrder;
  }

  //
  void  SetPolynomialOrderUsed( int polynomialOrderUsed )
  {
    m_PolynomialOrderUsed = polynomialOrderUsed;
  }

  int  GetPolynomialOrderUsed() const
  {
    return m_PolynomialOrderUsed;
  }

  //
  void  SetImage( const ImageType* image )
  {
    m_Image = image;
  }

  const ImageType*  GetImage() const
  {
    return m_Image;
  }

  //
  void  SetWeightImage( const ImageType* weightImage )
  {
    m_WeightImage = weightImage;
  }

  const ImageType*  GetWeightImage() const
  {
    return m_WeightImage;
  }

  //
  ImageType::Pointer  GetSmoothedImage( const MaskImageType* maskImage = 0, int voxelSizeFactor = 1 ) const;

protected:
  ImageSmoother();
  virtual ~ImageSmoother();

  void Initialize() const;

  void EstimateParameters() const;

  static vnl_matrix< double >  GetNonOrthogonalizedBasisFunctions( const MaskImageType*  maskImage, int polynomialOrder,
      const int* indicesWhereUnityIsReached = 0 );

  ImageType::Pointer  ExpandPolynomialToImage( const vnl_matrix< double >& basisFunctions, const MaskImageType* maskImage ) const;


private:
  ImageSmoother(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MaskImageType::ConstPointer  m_MaskImage;
  ImageType::ConstPointer  m_Image;
  ImageType::ConstPointer  m_WeightImage;
  int  m_PolynomialOrder;
  int  m_PolynomialOrderUsed;

  // Cached stuff
  mutable vnl_matrix< double >  m_BasisFunctions;
  mutable vnl_matrix< double >  m_Orthogonalizer;
  mutable vnl_vector< double >  m_Parameters;

};


} // end namespace kvl

#endif

