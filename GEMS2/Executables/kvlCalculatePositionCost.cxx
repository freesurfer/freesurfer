#include "itkImage.h"
#include "kvlAtlasMeshCollection.h"
#include "itkImageFileReader.h"
#include "kvlAtlasMeshCollectionPositionCostCalculator2.h"




int main( int argc, char* argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << argv[0] << " meshFileName imageFileName1 imageFileName2 ..." << std::endl;
    return -1;
    }
  
  // Read the collection from file
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  collection->Read( argv[ 1 ] );

  // Check that as many images were provided as there are meshes in the mesh collection
  if ( argc != static_cast< int >( 2 + collection->GetNumberOfMeshes() ) )
    {
    std::cerr << "Number of meshes must match number of images" << std::endl;
    return -1;
    }

  // Read the images
  typedef itk::Image< unsigned char, 3 >  LabelImageType;
  std::vector< LabelImageType::ConstPointer >  labelImages;
  for ( unsigned int meshNumber = 0; meshNumber < collection->GetNumberOfMeshes(); meshNumber++ )
    {
    typedef itk::ImageFileReader< LabelImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( argv[ 2 + meshNumber ] );
    reader->Update();
    labelImages.push_back( reader->GetOutput() );
    }


  // Now the real stuff
  kvl::AtlasMeshCollectionPositionCostCalculator::Pointer  calculator =
                                           kvl::AtlasMeshCollectionPositionCostCalculator::New();
  calculator->SetMeshCollection( collection );
  calculator->SetLabelImages( labelImages );
  calculator->DebugOn();
  for ( int i = 0; i<3; i++ )
    {
    const float positionCost = calculator->GetPositionCost();
    std::cout << "positionCost: " << positionCost << std::endl;
    }


  return 0;
};

