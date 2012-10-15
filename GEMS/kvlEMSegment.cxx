/**
 * @file  kvlEMSegment.cxx
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
#include "kvlEMSegmenter.h"
#include "kvlCroppedImageReader.h"
#include "itkImageFileWriter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "kvlAtlasMeshCollection.h"
#include "itkMGHImageIOFactory.h"
#include "itkShrinkImageFilter.h"
#if 0
#include "vnl/vnl_det.h"
#endif
#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"


typedef kvl::EMSegmenter::ImageType  ImageType;
typedef kvl::EMSegmenter::ClassificationImageType  ClassificationImageType;
typedef itk::Image< unsigned char, 3 >  UCharImageType;
typedef itk::Image< short, 3 >  ShortImageType;


// Somce forward declarations
ImageType::Pointer  LogTransform( const ImageType* image );

std::vector< UCharImageType::ConstPointer >  GetCrispClassifications( const std::vector< ClassificationImageType::ConstPointer >& images );

template < typename SourceImageType, typename TargetImageType >
typename TargetImageType::Pointer Convert( const SourceImageType* image,
    typename SourceImageType::PixelType sourceMinVal
    = itk::NumericTraits< typename SourceImageType::PixelType >::min(),
    typename SourceImageType::PixelType sourceMaxVal
    = itk::NumericTraits< typename SourceImageType::PixelType >::max(),
    typename TargetImageType::PixelType targetMinVal
    = itk::NumericTraits< typename TargetImageType::PixelType >::min(),
    typename TargetImageType::PixelType targetMaxVal
    = itk::NumericTraits< typename TargetImageType::PixelType >::max() );


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " imageFileName [ numberOfGaussians=2 biasFieldOrder=4 downSamplingFactor=1 meshCollectionFileName boundingFileName ]" << std::endl;
    return -1;
  }

  // Pare input
  const std::string imageFileName( argv[ 1 ] );
  int  numberOfGaussians = 2;
  int  biasFieldOrder = 4;
  int  downSamplingFactor = 1;
  std::string  meshCollectionFileName;
  std::string  boundingFileName;


  if ( argc > 2 )
  {
    std::istringstream  numberOfGaussiansStream( argv[ 2 ] );
    numberOfGaussiansStream >> numberOfGaussians;
  }

  if ( argc > 3 )
  {
    std::istringstream  biasFieldOrderStream( argv[ 3 ] );
    biasFieldOrderStream >> biasFieldOrder;
  }

  if ( argc > 4 )
  {
    std::istringstream  downSamplingFactorStream( argv[ 4 ] );
    downSamplingFactorStream >> downSamplingFactor;
  }

  if ( argc > 5 )
  {
    meshCollectionFileName = argv[ 5 ];
  }

  if ( argc > 6 )
  {
    boundingFileName = argv[ 6 ];
  }

  std::cout << "imageFileName: " << imageFileName << std::endl;
  std::cout << "numberOfGaussians: " << numberOfGaussians << std::endl;
  std::cout << "biasFieldOrder: " << biasFieldOrder << std::endl;
  std::cout << "downSamplingFactor: " << downSamplingFactor << std::endl;
  if ( meshCollectionFileName.size() )
  {
    std::cout << "meshCollectionFileName: " << meshCollectionFileName << std::endl;
  }
  if ( boundingFileName.size() )
  {
    std::cout << "boundingFileName: " << boundingFileName << std::endl;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  try
  {
    // Read image
    kvl::CroppedImageReader::Pointer  reader = kvl::CroppedImageReader::New();
    reader->SetExtraFraction( 0.1 );
    if ( boundingFileName.size() == 0 )
    {
      reader->Read( imageFileName.c_str() );
    }
    else
    {
      reader->Read( imageFileName.c_str(), boundingFileName.c_str() );
    }

    // Read mesh collection. If none if given, then make one
    kvl::AtlasMeshCollection::Pointer  meshCollection = kvl::AtlasMeshCollection::New();
    if ( meshCollectionFileName.size() == 0 )
    {
      // Create a simple mesh atlas, as required by the segmenter
      const unsigned int  meshSize[] = { 2, 2, 2 };
      unsigned int  domainSize[] = { reader->GetImage()->GetBufferedRegion().GetSize( 0 ),
                                     reader->GetImage()->GetBufferedRegion().GetSize( 1 ),
                                     reader->GetImage()->GetBufferedRegion().GetSize( 2 )
                                   };
      const float K = 1000.0f;
      const unsigned int numberOfClasses = 1;
      const unsigned int numberOfMeshes = 1;
      const bool  forceBorderVerticesToBackground = false;
      meshCollection->Construct( meshSize, domainSize, K, numberOfClasses, numberOfMeshes, forceBorderVerticesToBackground );
      //meshCollection->Write( "debugger.txt" );
    }
    else
    {
      if ( !meshCollection->Read( meshCollectionFileName.c_str() ) )
      {
        std::cerr << "Couldn't read mesh collection from file " << meshCollectionFileName << std::endl;
        exit( -1 );
      }
    }


    // Apply the correct transform to the mesh
    std::cout << "Applying transform: " << std::endl;
    reader->GetTransform()->Print( std::cout );
    meshCollection->Transform( -1, reader->GetTransform() );
    for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
    {
      meshCollection->Transform( i, reader->GetTransform() );
    }

#if 0
    const float  determinant = vnl_det( reader->GetTransform()->GetMatrix().GetVnlMatrix() );
    if ( determinant < 0 )
    {
      std::cout << "Careful here: the applied transformation will turn positive tetrahedra into negative ones." << std::endl;
      std::cout << reader->GetTransform()->GetMatrix().GetVnlMatrix() << std::endl;
      std::cout << " determinant: " << determinant << std::endl;
      std::cout << "Starting to swap the point assignments of each tetrahedron..." << std::endl;

      for ( kvl::AtlasMesh::CellsContainer::Iterator  cellIt = meshCollection->GetCells()->Begin();
            cellIt != meshCollection->GetCells()->End(); ++cellIt )
      {
        kvl::AtlasMesh::CellType*  cell = cellIt.Value();

        if( cell->GetType() != kvl::AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
          continue;
        }

        // Swap points assigned to first two vertices. This will readily turn negative tetrahedra
        //into positives ones.
        kvl::AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
        const kvl::AtlasMesh::PointIdentifier  p0Id = *pit;
        ++pit;
        const kvl::AtlasMesh::PointIdentifier  p1Id = *pit;

        pit = cell->PointIdsBegin();
        *pit = p1Id;
        ++pit;
        *pit = p0Id;
      } // End loop over all tetrahedra


      std::cout << "...done!" << std::endl;
    }
#endif

    // Write out the image that's actually being segmented
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( "imageBeingSegmented.mhd" );
    writer->SetInput( reader->GetImage() );
    writer->Update();

    // Log-transform the image
    std::cout << "Log-transforming data" << std::endl;
    ImageType::Pointer  image = LogTransform( reader->GetImage() ).GetPointer();

    // Threshold out background voxels, based on atlas prior
    if ( meshCollectionFileName.size() )
    {
      // Create an empty image that serves as a template for the alpha rasterizor
      typedef kvl::AtlasMeshAlphaDrawer::LabelImageType  LabelImageType;
      LabelImageType::Pointer  templateImage = LabelImageType::New();
      templateImage->SetRegions( image->GetLargestPossibleRegion().GetSize() );
      templateImage->Allocate();
      templateImage->FillBuffer( 0 );

      // Now rasterize the prior for the background class
      kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
      alphaDrawer->SetLabelImage( templateImage );
      alphaDrawer->SetLabelNumber( 0 );
      ( const_cast< kvl::AtlasMeshAlphaDrawer::AlphaImageType* >( alphaDrawer->GetAlphaImage() ) )->FillBuffer( 1.0 );
      alphaDrawer->Rasterize( meshCollection->GetReferenceMesh() );

      // Smooth it
      const float  backgroundPriorThresholdSmoothingSigma = 0.0f;
      typedef itk::DiscreteGaussianImageFilter< kvl::AtlasMeshAlphaDrawer::AlphaImageType,
              kvl::AtlasMeshAlphaDrawer::AlphaImageType >  SmootherType;
      SmootherType::Pointer smoother = SmootherType::New();
      smoother->SetInput( alphaDrawer->GetAlphaImage() );
      smoother->SetMaximumError( 0.1 );
      smoother->SetUseImageSpacingOff();
      smoother->SetVariance( pow( backgroundPriorThresholdSmoothingSigma, 2 ) );
      smoother->Update();
      kvl::AtlasMeshAlphaDrawer::AlphaImageType::ConstPointer  backgroundPrior = smoother->GetOutput();


      // Set all image intensities in areas < probability 1 for background to zero
      itk::ImageRegionIterator< kvl::CroppedImageReader::ImageType >  imageIt( image, image->GetLargestPossibleRegion() );
      itk::ImageRegionConstIterator< kvl::AtlasMeshAlphaDrawer::AlphaImageType >
      priorIt( backgroundPrior,
               backgroundPrior->GetLargestPossibleRegion() );
      for ( ; !imageIt.IsAtEnd(); ++imageIt, ++priorIt )
      {
        if ( priorIt.Value() > 0.99f )
        {
          imageIt.Value() = 0;
        }

      } // End loop over all voxels

    } // End test if we need to threshold out background voxels


    // Downsample the image
    ImageType::ConstPointer  nonDownsampledImage = image.GetPointer();
    kvl::AtlasMeshCollection::ConstPointer  nonDownsampledCollection = meshCollection.GetPointer();
    if ( downSamplingFactor > 1 )
    {
      // Downsample image
      typedef itk::ShrinkImageFilter< ImageType, ImageType >  ShrinkerType;
      ShrinkerType::Pointer  shrinker = ShrinkerType::New();
      shrinker->SetInput( image );
      shrinker->SetShrinkFactors( downSamplingFactor );
      shrinker->Update();
      image = shrinker->GetOutput();

      // Apply downsampling also to the mesh collection
      meshCollection = kvl::AtlasMeshCollection::New();
      meshCollection->GenerateFromSingleMesh(
        const_cast< kvl::AtlasMesh* >( nonDownsampledCollection->GetReferenceMesh().GetPointer() ), 1, 1000.0f );
      typedef kvl::AtlasMeshCollection::TransformType  TransformType;
      TransformType::Pointer  transform = TransformType::New();
      transform->Scale( 1 / static_cast< float >( downSamplingFactor ) );
      meshCollection->Transform( -1, transform );
      for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
      {
        meshCollection->Transform( i, transform );
      }
    }


    // Set up the EM segmenter. We allow more than one Gaussian to be assigned to a class in the atlas: if
    // numberOfGaussians is greater than the number of classes, the first class absorbs the extra ones.
    const unsigned int  numberOfClasses = meshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size();
    int  numberOfExtraGaussians = numberOfGaussians - numberOfClasses;
    if ( numberOfExtraGaussians < 0 )
    {
      std::cerr << "numberOfGaussians (" << numberOfGaussians
                << ") can't be smaller than the number of classes in the atlas mesh (" << numberOfClasses << ")" << std::endl;
      exit( -1 );
    }
    std::vector< unsigned int > lookupTable;
    lookupTable.push_back( 0 );
    for ( int i = 0; i < numberOfExtraGaussians; i++ )
    {
      lookupTable.push_back( 0 );
    }
    for ( unsigned int i = 1; i < numberOfClasses; i++ )
    {
      lookupTable.push_back( i );
    }
    kvl::EMSegmenter::Pointer  segmenter = kvl::EMSegmenter::New();
    segmenter->SetImage( image );
    segmenter->SetAtlasMesh( meshCollection->GetMesh( 0 ), lookupTable );
    segmenter->SetBiasFieldOrder( biasFieldOrder );
    //segmenter->SetStopCriterion( 1e-10 );
    segmenter->SetMaximumNumberOfIterations( 300 );
    segmenter->SetReestimatedRelativeWeightOfSharedClasses( false );


    // Let the beast go
    segmenter->Segment();


    ImageType::Pointer  correctedImage = 0;
    if ( downSamplingFactor > 1 )
    {
      //
      std::cout << "Parameters estimated on downsampled image; classicification at full resolution forced" << std::endl;

      correctedImage = kvl::EMSegmenter::ConvertToNativeImageType( segmenter->BiasCorrect( nonDownsampledImage, downSamplingFactor ) );
      segmenter->SetAtlasMesh( nonDownsampledCollection->GetMesh( 0 ), lookupTable );

      segmenter->SetImageAndSegmentItWithCurrentParametersAndWithoutBiasFieldCorrection( correctedImage );
    }
    else
    {
      correctedImage = kvl::EMSegmenter::ConvertToNativeImageType( segmenter->GetBiasCorrectedImage() );
    }

    // Make sure spacing, origin, and directions are copied
    correctedImage->SetSpacing( nonDownsampledImage->GetSpacing() );
    correctedImage->SetOrigin( nonDownsampledImage->GetOrigin() );
    correctedImage->SetDirection( nonDownsampledImage->GetDirection() );

    // Write out the results
    const std::string  biasCorrectedFileName = "biasCorrected.mgz";
    typedef itk::ImageFileWriter< ShortImageType >  ShortWriterType;
    ShortWriterType::Pointer  shortWriter = ShortWriterType::New();
    shortWriter->SetFileName( biasCorrectedFileName.c_str() );
#if 0
    shortWriter->SetInput( Convert< ImageType, ShortImageType >( correctedImage ) );
#else
    shortWriter->SetInput( Convert< ImageType, ShortImageType >( correctedImage,
                           itk::NumericTraits< ImageType::PixelType >::min(),
                           itk::NumericTraits< ImageType::PixelType >::max(),
                           0,
                           static_cast< ShortImageType::PixelType >( pow( 2, 12 ) - 1 ) ) );
#endif
    shortWriter->Update();
    std::cout << "Wrote " << biasCorrectedFileName << std::endl;


    std::vector< ClassificationImageType::ConstPointer >  posteriors;
    for ( unsigned int classNumber = 0; classNumber < segmenter->GetNumberOfClasses(); classNumber++ )
    {
      //
      std::ostringstream  posteriorFileNameStream;
      posteriorFileNameStream << "posterior_class" << classNumber << ".mgz";
      const std::string  posteriorFileName = posteriorFileNameStream.str();

      // Get the posterior
      ClassificationImageType::ConstPointer  posterior = segmenter->GetPosterior( classNumber );

      // Make sure spacing, origin, and directions are copied
      const_cast< ClassificationImageType* >( posterior.GetPointer() )->SetSpacing( nonDownsampledImage->GetSpacing() );
      const_cast< ClassificationImageType* >( posterior.GetPointer() )->SetOrigin( nonDownsampledImage->GetOrigin() );
      const_cast< ClassificationImageType* >( posterior.GetPointer() )->SetDirection( nonDownsampledImage->GetDirection() );

      //
      typedef itk::ImageFileWriter< UCharImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      writer->SetFileName( posteriorFileName.c_str() );
      writer->SetInput( Convert< ClassificationImageType, UCharImageType >( posterior, 0.0f, 1.0f ) );
      writer->Update();

      std::cout << "Wrote " << posteriorFileName << std::endl;

      posteriors.push_back( posterior );
    }


    // Write out crisp results
    std::vector< UCharImageType::ConstPointer >  crispClassifications = GetCrispClassifications( posteriors );
    for ( unsigned int classNumber = 0; classNumber < crispClassifications.size(); classNumber++ )
    {
      //
      std::ostringstream  crispPosteriorFileNameStream;
      crispPosteriorFileNameStream << "crispPosterior_class" << classNumber << ".mgz";
      const std::string  crispPosteriorFileName = crispPosteriorFileNameStream.str();

      //
      typedef itk::ImageFileWriter< UCharImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      writer->SetFileName( crispPosteriorFileName.c_str() );
      writer->SetInput( crispClassifications[ classNumber ] );
      writer->Update();

      std::cout << "Wrote " << crispPosteriorFileName << std::endl;
    }



  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }


  return 0;
};




//
//
//
ImageType::Pointer  LogTransform( const ImageType* image )
{
  // Allocate image
  typedef itk::Image< float, 3 >  LogImageType;
  LogImageType::Pointer  logImage = LogImageType::New();
  logImage->SetRegions( image->GetBufferedRegion() );
  logImage->Allocate();

  // Make sure spacing, origin, and directions are copied
  logImage->SetSpacing( image->GetSpacing() );
  logImage->SetOrigin( image->GetOrigin() );
  logImage->SetDirection( image->GetDirection() );

  // Loop over all voxels, log-transforming intensities and at the same
  // time calculate min/max values
  itk::ImageRegionConstIterator< ImageType >  imageIt( image,
      image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< LogImageType >  logIt( logImage,
      logImage->GetLargestPossibleRegion() );
  float  minimum = itk::NumericTraits< float >::max();
  float  maximum = itk::NumericTraits< float >::min();
  for ( ; !imageIt.IsAtEnd(); ++imageIt, ++logIt )
  {
    if ( imageIt.Value() > 0 )
    {
      logIt.Value() = log( imageIt.Value() );

      if ( logIt.Value() > maximum )
      {
        maximum = logIt.Value();
      }
      if ( logIt.Value() < minimum )
      {
        minimum = logIt.Value();
      }

    }
    else
    {
      logIt.Value() = 0;
    }
  }

  std::cout << "minimum: " << minimum << std::endl;
  std::cout << "maximum: " << maximum << std::endl;


  // Scale and clip intensities to fill range of original pixel type. Leave some room
  // for bias field correction
  const ImageType::PixelType  intensityRange = itk::NumericTraits< ImageType::PixelType >::max() - itk::NumericTraits< ImageType::PixelType >::min();
  const ImageType::PixelType  middleOfRange = ( itk::NumericTraits< ImageType::PixelType >::max() + itk::NumericTraits< ImageType::PixelType >::min() ) / 2;
  const ImageType::PixelType  lowerBound = middleOfRange - intensityRange / 2;
  const ImageType::PixelType  upperBound = middleOfRange + intensityRange / 2;
  typedef itk::IntensityWindowingImageFilter< LogImageType, ImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( logImage );
  windower->SetWindowMinimum( minimum );
  windower->SetWindowMaximum( maximum );
  windower->SetOutputMinimum( lowerBound );
  windower->SetOutputMaximum( upperBound );
  windower->Update();

  return  windower->GetOutput();

}




//
//
//
std::vector< UCharImageType::ConstPointer >  GetCrispClassifications( const std::vector< ClassificationImageType::ConstPointer >& images )
{

  // Determine, for each voxel, which class is maximal
  ClassificationImageType::Pointer  maximumProbabilityImage = ClassificationImageType::New();
  maximumProbabilityImage->SetRegions( images[ 0 ]->GetBufferedRegion() );
  maximumProbabilityImage->Allocate();
  maximumProbabilityImage->FillBuffer( 0.0f );

  UCharImageType::Pointer  bestClassImage = UCharImageType::New();
  bestClassImage->SetRegions( images[ 0 ]->GetBufferedRegion() );
  bestClassImage->Allocate();
  bestClassImage->FillBuffer( images.size() );

  for ( unsigned int classNumber = 0; classNumber < images.size(); classNumber++ )
  {
    // Loop over all voxels
    itk::ImageRegionConstIterator< ClassificationImageType >  sourceIt( images[ classNumber ],
        images[ classNumber ]->GetBufferedRegion() );
    itk::ImageRegionIterator< ClassificationImageType >  maxProbIt( maximumProbabilityImage,
        maximumProbabilityImage->GetBufferedRegion() );
    itk::ImageRegionIterator< UCharImageType >  bestClassIt( bestClassImage,
        bestClassImage->GetBufferedRegion() );
    for ( ; !sourceIt.IsAtEnd(); ++sourceIt, ++maxProbIt, ++bestClassIt )
    {
      if ( sourceIt.Value() > maxProbIt.Value() )
      {
        maxProbIt.Value() = sourceIt.Value();
        bestClassIt.Value() = classNumber;
      }
    }

  }




  // Loop over each class, each time creating an image with 0 or 255 filled in
  std::vector< UCharImageType::ConstPointer >   crispClassifications;
  for ( unsigned int classNumber = 0; classNumber < images.size(); classNumber++ )
  {
    // Create empty image
    UCharImageType::Pointer  crispClassification = UCharImageType::New();
    crispClassification->SetRegions( images[ 0 ]->GetBufferedRegion() );
    crispClassification->Allocate();
    crispClassification->FillBuffer( 0 );

    // Make sure spacing, origin, and directions are copied
    crispClassification->SetSpacing( images[ 0 ]->GetSpacing() );
    crispClassification->SetOrigin( images[ 0 ]->GetOrigin() );
    crispClassification->SetDirection( images[ 0 ]->GetDirection() );


    // Loop over all voxels
    itk::ImageRegionConstIterator< UCharImageType >  bestClassIt( bestClassImage,
        bestClassImage->GetBufferedRegion() );
    itk::ImageRegionIterator< UCharImageType >  crispIt( crispClassification,
        crispClassification->GetBufferedRegion() );
    for ( ; !bestClassIt.IsAtEnd(); ++bestClassIt, ++crispIt )
    {
      if ( bestClassIt.Value() == classNumber )
      {
        crispIt.Value() = itk::NumericTraits< UCharImageType::PixelType >::max();
      }
    }

    //
    crispClassifications.push_back( crispClassification.GetPointer() );
  }


  // Return result
  return crispClassifications;

}



//
//
//
template < typename SourceImageType, typename TargetImageType >
typename TargetImageType::Pointer Convert( const SourceImageType* image,
    typename SourceImageType::PixelType sourceMinVal
    = itk::NumericTraits< typename SourceImageType::PixelType >::min(),
    typename SourceImageType::PixelType sourceMaxVal
    = itk::NumericTraits< typename SourceImageType::PixelType >::max(),
    typename TargetImageType::PixelType targetMinVal
    = itk::NumericTraits< typename TargetImageType::PixelType >::min(),
    typename TargetImageType::PixelType targetMaxVal
    = itk::NumericTraits< typename TargetImageType::PixelType >::max() )
{
  typedef itk::IntensityWindowingImageFilter< SourceImageType, TargetImageType >   WindowerType;
  typename WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( image );
  windower->SetWindowMinimum( sourceMinVal );
  windower->SetWindowMaximum( sourceMaxVal );
  windower->SetOutputMinimum( targetMinVal );
  windower->SetOutputMaximum( targetMaxVal );
  windower->Update();

  return  windower->GetOutput();
};
