#ifndef __kvlMultiResolutionAtlasMesher_h
#define __kvlMultiResolutionAtlasMesher_h

#include "kvlAtlasParameterEstimator.h"


namespace kvl
{
 

class MultiResolutionAtlasMesher: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef MultiResolutionAtlasMesher  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiResolutionAtlasMesher, itk::Object );

  // Some typedefs
  typedef CompressionLookupTable::ImageType  LabelImageType;

  // Set up
  void SetUp( const std::vector< LabelImageType::ConstPointer >& labelImages,
              const CompressionLookupTable*  compressionLookupTable,
              const itk::Size< 3 >&  initialSize, 
              const std::vector< double >&  initialStiffnesse );

  //
  const AtlasMeshCollection*  GetCurrentMeshCollection() const
    { return m_Current; }

  //
  const AtlasParameterEstimator* GetEstimator() const
    { return m_Estimator; }

  //
  void  Go();

protected :
  // Constructor
  MultiResolutionAtlasMesher();

  // Destructor
  virtual ~MultiResolutionAtlasMesher();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;  

  
#ifdef USE_TETGEN
  //
  typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
  static void AddHexahedron( MeshSourceType* meshSource,
                             AtlasMesh::CellsContainer* hexahedra,
                             const AtlasMesh::PointType&  p0,
                             const AtlasMesh::PointType&  p1,
                             const AtlasMesh::PointType&  p2,
                             const AtlasMesh::PointType&  p3,
                             const AtlasMesh::PointType&  p4,
                             const AtlasMesh::PointType&  p5,
                             const AtlasMesh::PointType&  p6,
                             const AtlasMesh::PointType&  p7 )
    {
    AtlasMesh::PointIdentifier  p0IdDummy;
    AtlasMesh::PointIdentifier  p1IdDummy;
    AtlasMesh::PointIdentifier  p2IdDummy;
    AtlasMesh::PointIdentifier  p3IdDummy;
    AtlasMesh::PointIdentifier  p4IdDummy;
    AtlasMesh::PointIdentifier  p5IdDummy;
    AtlasMesh::PointIdentifier  p6IdDummy;
    AtlasMesh::PointIdentifier  p7IdDummy;

    AddHexahedron( meshSource, hexahedra,
                   p0, p1, p2, p3, p4, p5, p6, p7,
                   p0IdDummy, p1IdDummy, p2IdDummy, p3IdDummy, p4IdDummy,
                   p5IdDummy, p6IdDummy, p7IdDummy );

    }

  // 
  static void AddHexahedron( MeshSourceType* meshSource,
                             AtlasMesh::CellsContainer* hexahedra,
                             const AtlasMesh::PointType&  p0,
                             const AtlasMesh::PointType&  p1,
                             const AtlasMesh::PointType&  p2,
                             const AtlasMesh::PointType&  p3,
                             const AtlasMesh::PointType&  p4,
                             const AtlasMesh::PointType&  p5,
                             const AtlasMesh::PointType&  p6,
                             const AtlasMesh::PointType&  p7,
                             AtlasMesh::PointIdentifier&  p0Id,
                             AtlasMesh::PointIdentifier&  p1Id,
                             AtlasMesh::PointIdentifier&  p2Id,
                             AtlasMesh::PointIdentifier&  p3Id,
                             AtlasMesh::PointIdentifier&  p4Id,
                             AtlasMesh::PointIdentifier&  p5Id,
                             AtlasMesh::PointIdentifier&  p6Id,
                             AtlasMesh::PointIdentifier&  p7Id );


 

  template<class TCoordRep>
  static void GetUpsampledHexahedronPoints( const AtlasMesh::PointType&  p0,
                                            const AtlasMesh::PointType&  p1,
                                            const AtlasMesh::PointType&  p2,
                                            const AtlasMesh::PointType&  p3,
                                            const AtlasMesh::PointType&  p4,
                                            const AtlasMesh::PointType&  p5,
                                            const AtlasMesh::PointType&  p6,
                                            const AtlasMesh::PointType&  p7,
                                            AtlasMesh::PointType&  p01,
                                            AtlasMesh::PointType&  p02,
                                            AtlasMesh::PointType&  p13,
                                            AtlasMesh::PointType&  p23,
                                            AtlasMesh::PointType&  p15,
                                            AtlasMesh::PointType&  p37,
                                            AtlasMesh::PointType&  p26,
                                            AtlasMesh::PointType&  p04,
                                            AtlasMesh::PointType&  p57,
                                            AtlasMesh::PointType&  p67,
                                            AtlasMesh::PointType&  p46,
                                            AtlasMesh::PointType&  p45,
                                            AtlasMesh::PointType&  p0123,
                                            AtlasMesh::PointType&  p1357,
                                            AtlasMesh::PointType&  p2367,
                                            AtlasMesh::PointType&  p0246,
                                            AtlasMesh::PointType&  p0145,
                                            AtlasMesh::PointType&  p4567,
                                            AtlasMesh::PointType&  pMiddle );
#endif  


private :
  //
  MultiResolutionAtlasMesher(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

#ifdef USE_TETGEN  
  //
  AtlasMesh::CellsContainer::Pointer  GetCells( const AtlasMesh::PointsContainer* position ) const;

  //
  AtlasMeshCollection::Pointer  GetMeshCollection( AtlasMesh::PointsContainer* referencePosition,
                                                   std::vector< AtlasMesh::PointsContainer::Pointer >& positions,
                                                   double  stiffness ) const;
#endif                                                   

  //
  void Upsample();

  //
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;
  CompressionLookupTable::ConstPointer  m_CompressionLookupTable;
  itk::Size< 3 >  m_InitialSize;
  std::vector< double >  m_InitialStiffnesses;
  
  AtlasParameterEstimator::Pointer  m_Estimator;
  int  m_NumberOfClasses;
  int  m_NumberOfMeshes;
  itk::Size< 3 >  m_DomainSize;

  AtlasMeshCollection::Pointer  m_Current;

#ifdef USE_TETGEN  
  AtlasMesh::CellsContainer::Pointer  m_Hexahedra;
#endif

};



} // end namespace kvl


#endif
