/**
 * @file  kvlAtlasMeshCollectionPositionCostCalculator2.h
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
#ifndef __kvlAtlasMeshCollectionPositionCostCalculator_h
#define __kvlAtlasMeshCollectionPositionCostCalculator_h

#include "itkImage.h"
#include "kvlAtlasMeshCollection.h"
#include <vector>


namespace kvl
{


/**
 *
 */
class AtlasMeshCollectionPositionCostCalculator: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshCollectionPositionCostCalculator  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollectionPositionCostCalculator, itk::Object );

  // Some typedefs
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  //
  void SetMeshCollection( AtlasMeshCollection* meshCollection )
  {
    m_MeshCollection = meshCollection;
  }

  AtlasMeshCollection* GetMeshCollection()
  {
    return m_MeshCollection;
  }

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
  {
    m_LabelImages = labelImages;
  }

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
  {
    return m_LabelImages;
  }

  //
  float GetPositionCost()
  {
    int  numberOfProblemsDummy;
    return this->GetPositionCost( numberOfProblemsDummy );
  }

  //
  float GetPositionCost( int& numberOfProblems );

  //
  static void SetReturnZero( bool returnZero )
  {
    m_ReturnZero = returnZero;
  }

  //
  static bool  GetReturnZero()
  {
    return m_ReturnZero;
  }

protected:
  AtlasMeshCollectionPositionCostCalculator();
  virtual ~AtlasMeshCollectionPositionCostCalculator();

private:
  AtlasMeshCollectionPositionCostCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  AtlasMeshCollection::Pointer  m_MeshCollection;
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;

  static bool  m_ReturnZero;

};




} // end namespace kvl

#endif

