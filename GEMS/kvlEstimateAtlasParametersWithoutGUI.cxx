/**
 * @file  kvlEstimateAtlasParametersWithoutGUI.cxx
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
#include "kvlMultiResolutionAtlasMesher.h"
#include "itkCommand.h"
#include "kvlCompressionLookupTable.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "kvlAtlasMeshSummaryDrawer.h"
#include "itkImageFileWriter.h"


namespace kvl
{

class  MesherCommand : public itk::Command
{
public:
  /** Standard class typedefs. */
  typedef MesherCommand         Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self ) ;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MesherCommand, itk::Command );

  /** Abstract method that defines the action to be taken by the command. */
  virtual void Execute( itk::Object *caller, const itk::EventObject & event )
  {
    this->Execute(  const_cast< const Object* >( caller ), event );
  }

  /** Abstract method that defines the action to be taken by the command.
   * This variant is expected to be used when requests comes from a
   * const Object */
  virtual void Execute(const itk::Object *caller, const itk::EventObject & event )
  {

    //if ( typeid( event ) == typeid( UpsampledCollectionConstructedEvent ) )
    if ( typeid( event ) == typeid( UpsamplingEvent ) )
    {
      // Retrieve the mesher
      kvl::MultiResolutionAtlasMesher::ConstPointer  constMesher = static_cast< const kvl::MultiResolutionAtlasMesher* >(  caller );
      kvl::MultiResolutionAtlasMesher::Pointer  mesher = const_cast< kvl::MultiResolutionAtlasMesher* >(  constMesher.GetPointer() );

      if ( mesher->GetUpsamplingStepNumber() == ( mesher->GetNumberOfUpsamplingSteps() - 1 ) )
      {
        // OK, we just estimated a mesh and now on the verge of upsampling one last time. Quickly prune
        // the mesh collection before the upsampling occurs

        // Determine which points will be certainly kept: these are points that map to non-zero values in
        // *any* mask images
        std::set< AtlasMesh::PointIdentifier >  idsOfPointsToBeCertainlyKept;
        for ( int meshNumber = 0; meshNumber < mesher->GetCurrentMeshCollection()->GetNumberOfMeshes(); meshNumber++ )
        {
          for ( AtlasMesh::PointsContainer::ConstIterator  it = mesher->GetCurrentMeshCollection()->GetPositions()[ meshNumber ]->Begin();
                it != mesher->GetCurrentMeshCollection()->GetPositions()[ meshNumber ]->End(); ++it )
          {
            // Get the nearest neighbor image grid index
            MaskImageType::IndexType  index;
            for ( int i = 0; i < 3; i++ )
            {
              index[ i ] = static_cast< int >( it.Value()[ i ] + 0.5 );
            }

            if ( m_MaskImages[ meshNumber ]->GetPixel( index ) )
            {
              // std::cout << "In mesh number " << meshNumber << " position of point with id " << it.Index()
              //           << " is " << it.Value() << ", which maps to index " << index
              //           << " in the mask image, where the intensity is " << m_MaskImages[ meshNumber ]->GetPixel( index ) << std::endl;

              idsOfPointsToBeCertainlyKept.insert( it.Index() );
            }

          }
        }


        std::cout << "Starting to prune the mesher's mesh collection" << std::endl;
        std::cout << "   number of points in original mesh collection: "
                  << mesher->GetCurrentMeshCollection()->GetReferencePosition()->Size() << std::endl;
        std::cout << "   number of points to be certainly kept: "
                  << idsOfPointsToBeCertainlyKept.size() << std::endl;
        mesher->Prune( idsOfPointsToBeCertainlyKept );
        std::cout << "Finished pruning the mesher's mesh collection" << std::endl;
        mesher->GetCurrentMeshCollection()->Write( "estimatedPruned.txt" );
        std::cout << "Just wrote pruned mesh collection to estimatedPruned.txt" << std::endl;
      }
    }


  }


  //
  typedef kvl::CompressionLookupTable::ImageType  MaskImageType;
  void  SetMaskImages( const std::vector< MaskImageType::ConstPointer >&  maskImages )
  {
    m_MaskImages = maskImages;
  }

protected:
  MesherCommand() {}
  virtual ~MesherCommand() {}

private:
  MesherCommand(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::vector< MaskImageType::ConstPointer >  m_MaskImages;

};




class  MesherEstimatorCommand : public itk::Command
{
public:
  /** Standard class typedefs. */
  typedef MesherEstimatorCommand         Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self ) ;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MesherEstimatorCommand, itk::Command );

  /** Abstract method that defines the action to be taken by the command. */
  virtual void Execute( itk::Object *caller, const itk::EventObject & event )
  {
    this->Execute(  const_cast< const Object* >( caller ), event );
  }

  /** Abstract method that defines the action to be taken by the command.
   * This variant is expected to be used when requests comes from a
   * const Object */
  virtual void Execute(const itk::Object *caller, const itk::EventObject & event )
  {

    if ( ( typeid( event ) == typeid( itk::IterationEvent ) ) ||
         ( typeid( event ) == typeid( itk::EndEvent ) ) )
    {
      kvl::AtlasParameterEstimator::ConstPointer  estimator = static_cast< const kvl::AtlasParameterEstimator* >(  caller );

      // Write out the current mesh collection
      std::ostringstream  fileNameStream;
      fileNameStream << "estimatedMeshForUpsamplingNumber" << m_Mesher->GetUpsamplingStepNumber();
      if ( typeid( event ) == typeid( itk::IterationEvent ) )
      {
        fileNameStream  << "_estimatorIteration" << estimator->GetIterationNumber();
      }
      else
      {
        // The estimator updates the alphas once more after the final registration has been updated
        fileNameStream  << "_estimatorFinalEstimation";
      }
      fileNameStream << ".txt";
      std::string  fileName = fileNameStream.str();
      if ( !estimator->GetCurrentMeshCollection()->Write( fileName.c_str() ) )
      {
        std::cerr <<  "Couldn't write mesh collection to file " << fileName << std::endl;
        exit( -1 );
      }
      std::cout << "Just wrote to " << fileName << std::endl;

      // Just for illustration's sake, also write out a summary image of the atlas in reference position
      AtlasMeshSummaryDrawer::Pointer  summaryDrawer = AtlasMeshSummaryDrawer::New();
      summaryDrawer->SetLabelImage( m_TemplateImage );
      summaryDrawer->Rasterize( estimator->GetCurrentMeshCollection()->GetReferenceMesh() );
      typedef itk::ImageFileWriter< AtlasMeshSummaryDrawer::SummaryImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      fileNameStream << ".mhd";
      writer->SetFileName( fileNameStream.str().c_str() );
      writer->SetInput( summaryDrawer->GetSummaryImage() );
      writer->Update();
    }

  }

  void  SetMesher( const kvl::MultiResolutionAtlasMesher*  mesher )
  {
    m_Mesher = mesher;
  }

  //
  typedef AtlasMeshSummaryDrawer::LabelImageType  TemplateImageType;
  void  SetTemplateImage( const TemplateImageType* templateImage )
  {
    m_TemplateImage = const_cast< TemplateImageType* >( templateImage );
  }


protected:
  MesherEstimatorCommand()
  {
    m_Mesher = 0;
    m_TemplateImage = 0;
  }
  virtual ~MesherEstimatorCommand() {}

private:
  MesherEstimatorCommand(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  kvl::MultiResolutionAtlasMesher::ConstPointer  m_Mesher;
  TemplateImageType::ConstPointer  m_TemplateImage;

};



};




int main( int argc, char** argv )
{

  // Sanity check on input
  if ( argc < 9 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " numberOfUpsamplingSteps meshSizeX meshSizeY meshSizeZ stiffness useGaussians useMask imageFileName1 [ imageFileName2 ... maskImageFileName1 maskImageFile2 ... ]" << std::endl;

    return -1;
  }



  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 8; argumentNumber++ )
  {
    inputParserStream << argv[ argumentNumber ] << " ";
  }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  unsigned int  numberOfUpsamplingSteps;
  unsigned int  meshSizeX;
  unsigned int  meshSizeY;
  unsigned int  meshSizeZ;
  float  stiffness;
  bool  useGaussians;
  bool  useMask;
  inputStream >> numberOfUpsamplingSteps >>  meshSizeX >> meshSizeY >> meshSizeZ >> stiffness >> useGaussians >> useMask;

  try
  {
    // Determine which is the end of the image file names and the start of the masks
    int  originalImageFileNamesLastArgumentNumber = argc-1;
    if ( useMask )
    {
      const int  totalNumberOfImages = argc - 8;
      if ( totalNumberOfImages % 2 )
      {
        std::cerr << "You should the same number of masks as there are images" << std::endl;
        exit( -1 );
      }

      originalImageFileNamesLastArgumentNumber = 7 + totalNumberOfImages / 2;
    }

    // Read the original images
    typedef kvl::CompressionLookupTable::ImageType  InputImageType;
    std::vector< InputImageType::ConstPointer >  originalImages;
    for ( int argumentNumber = 8; argumentNumber <= originalImageFileNamesLastArgumentNumber; argumentNumber++ )
    {
      // Read the input image
      typedef itk::ImageFileReader< InputImageType >  ReaderType;
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( argv[ argumentNumber ] );
      reader->Update();
      InputImageType::ConstPointer  originalImage = reader->GetOutput();

      // Over-ride the spacing and origin since at this point we can't deal with that
      const double spacing[] = { 1, 1, 1 };
      const double origin[] = { 0, 0, 0 };
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetSpacing( spacing );
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetOrigin( origin );

      // Remember this image
      originalImages.push_back( originalImage );
    }

    // Read the maks images
    std::vector< InputImageType::ConstPointer >  maskImages;
    for ( int argumentNumber = ( originalImageFileNamesLastArgumentNumber + 1 );
          argumentNumber < argc; argumentNumber++ )
    {
      // Read the input image
      typedef itk::ImageFileReader< InputImageType >  ReaderType;
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( argv[ argumentNumber ] );
      reader->Update();
      InputImageType::ConstPointer  maskImage = reader->GetOutput();

      // Over-ride the spacing and origin since at this point we can't deal with that
      const double spacing[] = { 1, 1, 1 };
      const double origin[] = { 0, 0, 0 };
      const_cast< InputImageType* >( maskImage.GetPointer() )->SetSpacing( spacing );
      const_cast< InputImageType* >( maskImage.GetPointer() )->SetOrigin( origin );

      // Remember this image
      maskImages.push_back( maskImage );
    }


    // Build up the internal label images (in uchar - whereas the original images are in ushort)
    typedef kvl::CompressionLookupTable::CompressedImageType  OutputImageType;
    kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
    std::vector< OutputImageType::ConstPointer >  labelImages;
    if ( !useGaussians )
    {
      // Build a lookup table that maps the original intensities onto uchar starting
      // at 0 and densely packed
      compressor->Construct( originalImages );
      compressor->Write( "compressionLookupTable.txt" );

      // Collect the label images resulting from pushing the original images through the
      // lookup table
      for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
            it != originalImages.end(); ++it )
      {
        labelImages.push_back( ( compressor->CompressImage( ( *it ).GetPointer() ) ).GetPointer() );
      }
    }
    else
    {
      // Just convert the original images into uchar
      for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
            it != originalImages.end(); ++it )
      {
        typedef itk::CastImageFilter< InputImageType, OutputImageType >   CasterType;
        CasterType::Pointer  caster = CasterType::New();
        caster->SetInput( ( *it ).GetPointer() );
        caster->Update();
        labelImages.push_back( caster->GetOutput() );
      }

    }


    // Set up the multi resolution mesher
    kvl::MultiResolutionAtlasMesher::Pointer  mesher = kvl::MultiResolutionAtlasMesher::New();
    mesher->SetLabelImages( labelImages );
    mesher->SetNumberOfUpsamplingSteps( numberOfUpsamplingSteps );
    mesher->SetTryToBeSparse( false );
    unsigned int  size[ 3 ] = { meshSizeX,  meshSizeY, meshSizeZ };
    float  stiffnesses[ 5 ] = { stiffness, stiffness, stiffness, stiffness, stiffness };
    mesher->SetUp( size, stiffnesses );
    const int  numberOfClasses = mesher->GetNumberOfClasses();
    std::cout << "Set up with numberOfClasses: " << numberOfClasses << std::endl;
    if ( useMask )
    {
      kvl::MesherCommand::Pointer  mesherCommand = kvl::MesherCommand::New();
      mesherCommand->SetMaskImages( maskImages );
      //mesher->AddObserver( kvl::UpsampledCollectionConstructedEvent(), mesherCommand );
      mesher->AddObserver( kvl::UpsamplingEvent(), mesherCommand );
    }
    kvl::MesherEstimatorCommand::Pointer  mesherEstimatorCommand = kvl::MesherEstimatorCommand::New();
    mesherEstimatorCommand->SetMesher( mesher );
    mesherEstimatorCommand->SetTemplateImage( labelImages[ 0 ] );
    mesher->AddEstimatorObserver( itk::IterationEvent(), mesherEstimatorCommand );
    mesher->AddEstimatorObserver( itk::EndEvent(), mesherEstimatorCommand );
    mesher->SetUseGaussians( useGaussians );


    // Let the beast go
    mesher->Go();
  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};


