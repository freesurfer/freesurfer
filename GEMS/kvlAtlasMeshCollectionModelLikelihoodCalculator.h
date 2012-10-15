/**
 * @file  kvlAtlasMeshCollectionModelLikelihoodCalculator.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#ifndef __kvlAtlasMeshCollectionModelLikelihoodCalculator_h
#define __kvlAtlasMeshCollectionModelLikelihoodCalculator_h

#include "kvlAtlasMeshCollection.h"
#include "itkImage.h"


namespace kvl
{


/**
 *
 */
class AtlasMeshCollectionModelLikelihoodCalculator: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshCollectionModelLikelihoodCalculator  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollectionModelLikelihoodCalculator, itk::Object );

  /** Some typedefs */
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
  {
    return m_LabelImages;
  }

  // Get label image
  const LabelImageType*  GetLabelImage( unsigned int labelImageNumber ) const;

  // Initialize
  void SetMeshCollection( const AtlasMeshCollection* meshCollection )
  {
    m_MeshCollection = meshCollection;
  }

  //
  const AtlasMeshCollection*  GetMeshCollection() const
  {
    return m_MeshCollection;
  }


  /** */
  void GetDataCostAndAlphasCost( float& dataCost, float& alphasCost ) const;

protected:
  AtlasMeshCollectionModelLikelihoodCalculator();
  virtual ~AtlasMeshCollectionModelLikelihoodCalculator();

private:
  AtlasMeshCollectionModelLikelihoodCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //
  AtlasMeshCollection::ConstPointer  m_MeshCollection;

  std::vector< LabelImageType::ConstPointer >  m_LabelImages;

  int  m_NumberOfLabelImages;


};




} // end namespace kvl

#endif

