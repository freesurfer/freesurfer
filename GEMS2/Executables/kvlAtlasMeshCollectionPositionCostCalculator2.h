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
    { return m_MeshCollection; }

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
    { m_LabelImages = labelImages; }

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
    { return m_LabelImages; }

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
    { m_ReturnZero = returnZero; }

  //
  static bool  GetReturnZero()
    { return m_ReturnZero; }

  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    { m_mapCompToComp = mapCompToComp; }
  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }


protected:
  AtlasMeshCollectionPositionCostCalculator();
  virtual ~AtlasMeshCollectionPositionCostCalculator();

private:
  AtlasMeshCollectionPositionCostCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  AtlasMeshCollection::Pointer  m_MeshCollection;
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;

  static bool  m_ReturnZero;

  std::vector<unsigned char > *m_mapCompToComp;

};




} // end namespace kvl

#endif

