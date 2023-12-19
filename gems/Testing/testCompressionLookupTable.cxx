#include "itkImageFileReader.h"
#include "itkMGHImageIOFactory.h"
#include "kvlCompressionLookupTable.h"

int main(int argc, char **argv)
{
  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );
  
  // Read the input images
  typedef kvl::CompressionLookupTable::ImageType  LabelImageType;
  std::vector< LabelImageType::ConstPointer >  labelImages;
  for ( int argumentNumber = 1; argumentNumber < argc; argumentNumber++ )
  {
    std::cout << "Reading input image: " << argv[ argumentNumber ] << std::endl;
    // Read the input image
    typedef itk::ImageFileReader< LabelImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( argv[ argumentNumber ] );
    reader->Update();
    LabelImageType::ConstPointer  labelImage = reader->GetOutput();

    // Over-ride the spacing and origin since at this point we can't deal with that
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    const_cast< LabelImageType* >( labelImage.GetPointer() )->SetSpacing( spacing );
    const_cast< LabelImageType* >( labelImage.GetPointer() )->SetOrigin( origin );

    // Remember this image
    labelImages.push_back( labelImage );
  }  

  // Build a lookup table that maps the original intensities onto class numbers starting
  // at 0 and densely packed
  kvl::CompressionLookupTable::Pointer  lookupTable = kvl::CompressionLookupTable::New();
  lookupTable->Construct( labelImages );
  lookupTable->Write( "compressionLookupTable.txt" );

  /*********** START OF SIMULATION ***********/
  printf("\nClasses Contributed To Label Cost Calculations:\n");
  // loop through each label,
  // report which classes will contribute to cost/gradient/likelihood/color calculation
  std::vector<LabelImageType::PixelType> labels = lookupTable->GetLabels();
  std::vector<LabelImageType::PixelType>::const_iterator labelIt;
  for (labelIt = labels.begin(); labelIt != labels.end(); labelIt++)
  {
    printf("label %5d <= ", *labelIt);
    const std::vector< int >& classNumbers = lookupTable->GetClassNumbers(*labelIt);
    for ( std::vector< int >::const_iterator classIt = classNumbers.begin(); 
          classIt != classNumbers.end();
          ++classIt )
      printf(" %3d, ", *classIt);

    printf("\n");
  }

  printf("\nTotal # of Labels: %d\n", labels.size());
  printf("Total # of Classes: %d\n", lookupTable->GetNumberOfClasses());

  exit(0);
}
