/**
 * @file  kvlMapLabels.cxx
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
#include "itkCommand.h"
#include "itkImageFileWriter.h"




/**
 *
 * Read a deformed mesh, replace its alphas with those of the mesh of label
 * probabilities, and rasterize into an image to will contain the label with
 * maximal probability
 *
 */



int main( int argc, char** argv )
{

  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " deformedMesh labelFileName1 [ labelFileName2 ... ]" << std::endl;

    return -1;
  }


  try
  {
    typedef kvl::CompressionLookupTable::ImageType  InputImageType;
    std::vector< InputImageType::ConstPointer >  originalImages;
    for ( int argumentNumber = 2; argumentNumber < argc; argumentNumber++ )
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


    // Build up the internal label images (in uchar - whereas the original images are in ushort)
    typedef kvl::CompressionLookupTable::CompressedImageType  OutputImageType;
    kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
    compressor->Construct( originalImages );
    compressor->Write( "compressionLookupTable.txt" );

    // Collect the label images resulting from pushing the original images through the
    // lookup table
    std::vector< OutputImageType::ConstPointer >  labelImages;
    for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
          it != originalImages.end(); ++it )
    {
      labelImages.push_back( ( compressor->CompressImage( ( *it ).GetPointer() ) ).GetPointer() );
    }


    // Read the mesh collection
    kvl::AtlasMeshCollection::Pointer  collection = kvl::AtlasMeshCollection::New();
    if ( !collection->Read( argv[ 1 ] ) )
    {
      std::cerr << "Couldn't read mes collection from file " << argv[ 1 ] << std::endl;
      exit( -1 );
    }


    // Replace the original alphas with flat alphas of the correct length
    const int  numberOfClasses = compressor->GetCompressionLookupTable().size();
    std::cout << "Found that we have " << numberOfClasses << " classes" << std::endl;
    kvl::AtlasAlphasType  flatAlphas( numberOfClasses );
    flatAlphas.Fill( 1.0f / static_cast< float >( numberOfClasses ) );
    for ( kvl::AtlasMesh::PointDataContainer::Iterator it = collection->GetPointParameters()->Begin();
          it != collection->GetPointParameters()->End(); ++it )
    {
      it.Value().m_Alphas = flatAlphas;
    }


    // Now estimate the alphas using an AtlasParameterEstimator
    kvl::AtlasParameterEstimator::Pointer  estimator = kvl::AtlasParameterEstimator::New();
    estimator->SetLabelImages( labelImages );
    estimator->SetInitialMeshCollection( collection );
    estimator->SetMaximumNumberOfIterations( 0 );
    kvl::EstimatorCommand::Pointer  command = kvl::EstimatorCommand::New();
    estimator->AddObserver( kvl::AlphasEstimationStartEvent(), command );
    estimator->AddObserver( kvl::AlphasEstimationIterationEvent(), command );
    estimator->AddObserver( kvl::AlphasEstimationEndEvent(), command );
    estimator->Estimate();


    // Write out
    const std::string  fileName = "calculatedLabelProbabilities.txt";
    if ( !collection->Write( fileName.c_str() ) )
    {
      std::cerr << "Couldn't write results to " << fileName << std::endl;
      exit( -1 );
    }
    std::cout << "Just wrote results to " << fileName << std::endl;

  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};


