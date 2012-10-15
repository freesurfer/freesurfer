/**
 * @file  kvlMapDeformationField.cxx
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
#include "itkImage.h"
#include "kvlAtlasMeshCollection.h"



#include "kvlAtlasMeshRasterizor.h"


namespace kvl
{


namespace FragmentProcessor
{

/**
 *
 */
class MapPosition
{
public:

  typedef itk::Image< AtlasMesh::PointType, 3 >  PositionImageType;

  MapPosition()
  {
    m_PositionImage = 0;
    m_Mesh = 0;
    m_MeshToMap = 0;
  }

  ~MapPosition() {};

  void  SetMeshToMap( const AtlasMesh* meshToMap )
  {
    m_MeshToMap = meshToMap;
  }

  void AllocatePositionImage( PositionImageType::SizeType  size )
  {
    //std::cout << "Allocating image of size " << size << std::endl;
    m_PositionImage = PositionImageType::New();
    m_PositionImage->SetRegions( size );
    m_PositionImage->Allocate();
    kvl::AtlasMesh::PointType  uncovered;
    for ( int i = 0; i < 3; i++ )
    {
      uncovered[ i ] = itk::NumericTraits< float >::max();
    }
    m_PositionImage->FillBuffer( uncovered );
  }

  const PositionImageType* GetPositionImage() const
  {
    return m_PositionImage;
  }

  inline void operator()( const float& pi0, const float& pi1, const float& pi2, const float& pi3 )
  {
    kvl::AtlasMesh::PointType  position;
    for ( int i = 0; i < 3; i++ )
    {
      position[ i ] = pi0 * m_PositionInVertex0[ i ] +
                      pi1 * m_PositionInVertex1[ i ] +
                      pi2 * m_PositionInVertex2[ i ] +
                      pi3 * m_PositionInVertex3[ i ];
    }

    // std::cout << "               Trying to fill in pixel with index " << m_Index << " with position " << position << std::endl;
    // std::cout << "                    pi0: " << pi0 << std::endl;
    // std::cout << "                    pi1: " << pi1 << std::endl;
    // std::cout << "                    pi2: " << pi2 << std::endl;
    // std::cout << "                    pi3: " << pi3 << std::endl;
    //std::cout << "The image itself is size " << m_PositionImage << std::endl;
    m_PositionImage->SetPixel( m_Index, position );

#if 0
    if ( ( m_Index[ 0 ] == 0 ) && ( m_Index[ 1 ] == 1 ) && ( m_Index[ 2 ] == 0 ) )
    {
      std::cout << "Filled index " << m_Index << " with position " << position << std::endl;
    }
#endif

    m_Index[ 0 ]++;
  }

  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
  {
    m_Index[ 0 ] = x;
    m_Index[ 1 ] = y;
    m_Index[ 2 ] = z;
  }

  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
  {
    //std::cout << "Starting to rasterize tetrahedron with id " << cellId << std::endl;

    // Cache the alpha of the specified class in each of the vertices of this tetrahedron
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );

    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    m_PositionInVertex0 = m_MeshToMap->GetPoints()->ElementAt( *pit );
    ++pit;
    m_PositionInVertex1 = m_MeshToMap->GetPoints()->ElementAt( *pit );
    ++pit;
    m_PositionInVertex2 = m_MeshToMap->GetPoints()->ElementAt( *pit );
    ++pit;
    m_PositionInVertex3 = m_MeshToMap->GetPoints()->ElementAt( *pit );

    //std::cout << "    m_PositionInVertex0: " << m_PositionInVertex0 << std::endl;
    //std::cout << "    m_PositionInVertex1: " << m_PositionInVertex1 << std::endl;
    //std::cout << "    m_PositionInVertex2: " << m_PositionInVertex2 << std::endl;
    //std::cout << "    m_PositionInVertex3: " << m_PositionInVertex3 << std::endl;
#if 0
    {
      AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
      std::cout << "Rasterizing tetrahedron: " << std::endl;
      for ( int i = 0; i < 4; i++, ++pit )
      {
        std::cout << "      " << m_Mesh->GetPoints()->ElementAt( *pit ) << " is mapped onto"
                  << m_MeshToMap->GetPoints()->ElementAt( *pit )  << std::endl;
      }

      std::cout << "    m_PositionInVertex0: " << m_PositionInVertex0 << std::endl;
      std::cout << "    m_PositionInVertex1: " << m_PositionInVertex1 << std::endl;
      std::cout << "    m_PositionInVertex2: " << m_PositionInVertex2 << std::endl;
      std::cout << "    m_PositionInVertex3: " << m_PositionInVertex3 << std::endl;

    }
#endif

    return true;
  }

  inline void SetMesh( const AtlasMesh* mesh )
  {
    //std::cout << "Getting mesh " << mesh << std::endl;
    m_Mesh = mesh;
  }

private:

  PositionImageType::Pointer  m_PositionImage;
  PositionImageType::IndexType  m_Index;

  AtlasMesh::PointType  m_PositionInVertex0;
  AtlasMesh::PointType  m_PositionInVertex1;
  AtlasMesh::PointType  m_PositionInVertex2;
  AtlasMesh::PointType  m_PositionInVertex3;

  AtlasMesh::ConstPointer  m_Mesh;
  AtlasMesh::ConstPointer  m_MeshToMap;


};


} // End namespace FragmentProcessor


/**
 *
 */
class AtlasMeshPositionMapper: public AtlasMeshRasterizor< FragmentProcessor::MapPosition >
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshPositionMapper  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::MapPosition >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshPositionMapper, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::PositionImageType  PositionImageType;

  /** */
  virtual void SetLabelImage( const LabelImageType*  labelImage )
  {
    // Use the label image as a template for the alpha image
    this->GetFragmentProcessor().AllocatePositionImage( labelImage->GetLargestPossibleRegion().GetSize() );

    // Invoke superclass' implementation
    Superclass::SetLabelImage( labelImage );
  }

  /** */
  void  SetMeshToMap( const AtlasMesh* meshToMap )
  {
    this->GetFragmentProcessor().SetMeshToMap( meshToMap );
  }

  /** */
  const PositionImageType*  GetPositionImage() const
  {
    return this->GetFragmentProcessor().GetPositionImage();
  }

protected:
  AtlasMeshPositionMapper() {};
  virtual ~AtlasMeshPositionMapper() {};

private:
  AtlasMeshPositionMapper(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


};


} // end namespace kvl






int main( int argc, char* argv[] )
{

  // Check input
  if ( argc != 6 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " sourceMeshCollection targetMeshCollection "
              << " numberOfMeshNodesX numberOfMeshNodesY numberOfMeshNodesZ" << std::endl;
    return( -1 );
  }

  // Read collections from file
  kvl::AtlasMeshCollection::Pointer  sourceCollection = kvl::AtlasMeshCollection::New();
  sourceCollection->Read( argv[ 1 ] );

  kvl::AtlasMeshCollection::Pointer  targetCollection = kvl::AtlasMeshCollection::New();
  targetCollection->Read( argv[ 2 ] );

#if 0
  {
    unsigned int  meshSize[] = { 2, 2, 2 };
    unsigned int  domainSize[] = { 2, 2, 2 };
    sourceCollection->Construct( meshSize, domainSize, 1000, 2, 1 );
    targetCollection->Construct( meshSize, domainSize, 1000, 2, 1 );
  }
#endif

  // Retrieve the original mesh size
  kvl::AtlasMeshPositionMapper::LabelImageType::SizeType  originalMeshSize;
  originalMeshSize[ 0 ] = 0;
  originalMeshSize[ 1 ] = 0;
  originalMeshSize[ 2 ] = 0;
  for ( kvl::AtlasMesh::PointsContainer::ConstIterator  it = targetCollection->GetReferencePosition()->Begin();
        it != targetCollection->GetReferencePosition()->End(); ++it )
  {
    //std::cout << "Found vertex with reference position " << it.Value() << std::endl;
    for ( int i = 0; i < 3; i++ )
    {
      if ( it.Value()[ i ] > originalMeshSize[ i ] )
      {
        originalMeshSize[ i ] = static_cast< kvl::AtlasMeshPositionMapper::LabelImageType::SizeType::SizeValueType >( it.Value()[ i ] );
      }
    }
  }
  originalMeshSize[ 0 ]++;
  originalMeshSize[ 1 ]++;
  originalMeshSize[ 2 ]++;
  for ( int i = 0; i < 3; i++ )
  {
    std::cout << "Got originalMeshSize[" << i << "]: " << originalMeshSize[ i ] << std::endl;
  }

  // Retrieve the target mesh size
  kvl::AtlasMeshPositionMapper::LabelImageType::SizeType  size;
  for ( int i = 0; i < 3; i++ )
  {
    std::istringstream  sizeStream( argv[ 3 + i ] );
    sizeStream >> size[ i ];
  }
  for ( int i = 0; i < 3; i++ )
  {
    std::cout << "Got size[" << i << "]: " << size[ i ] << std::endl;
  }

  // Create an empty template image for the rasterizor
  kvl::AtlasMeshPositionMapper::LabelImageType::Pointer  templateLabelImage = kvl::AtlasMeshPositionMapper::LabelImageType::New();
  templateLabelImage->SetRegions( size );
  templateLabelImage->Allocate();
  //templateLabelImage->FillBuffer( 0 );


  // Adjust the position of the mesh to map to the correct size
  for ( kvl::AtlasMesh::PointsContainer::Iterator  it = sourceCollection->GetReferencePosition()->Begin();
        it != sourceCollection->GetReferencePosition()->End(); ++it )
  {
    for ( int i = 0; i < 3; i++ )
    {
      it.Value()[ i ] /= ( static_cast< float >( originalMeshSize[ i ] - 1 ) / ( size[ i ] - 1 ) );
      if ( fabs( it.Value()[ i ] - ( size[ i ] - 1 ) ) < 1e-3 )
      {
        //std::cout << "Heeeeeeeeeee having voxel right on border that won't be rendered " << it.Value()[ i ] << std::endl;
        it.Value()[ i ] += 1e-3;
        //std::cout << "   Changed a little to " << it.Value()[ i ] << std::endl;
      }
    }
  }

  // For each mesh in the collection, rasterize the mesh node source positions onto an image,
  // and then loop over all nodes in the target and assign them the position most close
  bool  badRoundingErrors = false;
  for ( unsigned int meshNumber = 0; meshNumber < sourceCollection->GetNumberOfMeshes(); meshNumber++ )
  {

    //std::cout << "Rasterizing mesh number " << meshNumber << std::endl;
    kvl::AtlasMeshPositionMapper::Pointer  mapper = kvl::AtlasMeshPositionMapper::New();
    mapper->SetLabelImage( templateLabelImage );
    mapper->SetMeshToMap( sourceCollection->GetMesh( meshNumber ) );
    mapper->Rasterize( sourceCollection->GetReferenceMesh() );
    //std::cout << "Rasterized mesh number " << meshNumber << std::endl;

    kvl::AtlasMesh::PointsContainer::ConstIterator  refIt = targetCollection->GetReferencePosition()->Begin();
    kvl::AtlasMesh::PointsContainer::Iterator  it = targetCollection->GetPositions()[ meshNumber ]->Begin();
    for ( ; refIt != targetCollection->GetReferencePosition()->End(); ++refIt, ++it )
    {
      //std::cout << "Lookup up position for vertex with reference position " << refIt.Value() << std::endl;
      kvl::AtlasMeshPositionMapper::PositionImageType::IndexType  index;
      for ( int i = 0; i < 3; i++ )
      {
        const float  normalizedReferencePosition = ( refIt.Value()[ i ] ) / static_cast< float >( originalMeshSize[ i ] - 1 ) *
            ( size[ i ] - 1 );
        index[ i ] = static_cast< int >( normalizedReferencePosition + 0.5 );
        if ( fabs( index[ i ] - normalizedReferencePosition ) > 1e-5 )
        {
          badRoundingErrors = true;
        }

      }

      //
      //std::cout << "    found that it was at index " << index << std::endl;
      kvl::AtlasMesh::PointType  mappedPoint = mapper->GetPositionImage()->GetPixel( index );
      if ( mappedPoint[ 0 ] != itk::NumericTraits< float >::max() )
      {
        //std::cout << "    and that it maps to " << mappedPoint << std::endl;

        // Un-normalize
        // for ( int i = 0; i < 3; i++ )
        //   {
        //   if ( fabs( mappedPoint[ i ] - size[ i ] ) < 1e-1 )
        //     {
        //     std::cout << "Heeeeeeeeeee having voxel right on border that wouldn't have been rendered" << mappedPoint[ i ] << std::endl;
        //     mappedPoint[ i ] -= 1e-1;
        //     std::cout << "   Changed a little back to " << mappedPoint[ i ] << std::endl;
        //     }
        //   mappedPoint[ i ] /= ( static_cast< float >( size[ i ] - 1 ) / ( originalMeshSize[ i ] - 1 ) );
        //   }
        it.Value() = mappedPoint;
        //std::cout << "    and that it was finally " << mappedPoint << std::endl;
      }
      else
      {
        std::cout << "Oooooooooooooooops that voxel is not covered!" << std::endl;
      }


    } // End loop over all mesh nodes

  } // End loop over all meshes


  if ( badRoundingErrors )
  {
    std::cout << "Warning: bad rounding errors occurred!" << std::endl;
  }


  // Write out results
  targetCollection->Write( "mapped.txt" );

  return 0;
};

