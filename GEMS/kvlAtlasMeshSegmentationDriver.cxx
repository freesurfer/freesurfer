/**
 * @file  kvlAtlasMeshSegmentationDriver.cxx
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
#include "kvlAtlasMeshSegmentationDriver.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCommand.h"
#include "kvlAtlasMeshCollection.h"
#include "itkShiftScaleImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkAffine3DTransform.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "kvlAtlasMeshSmoother.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkMGHImageIOFactory.h"



namespace kvl
{

//
//
//
AtlasMeshSegmentationDriver
::AtlasMeshSegmentationDriver()
{

  m_Segmenter = AtlasMeshSegmenter::New();
  m_NonDownsampledImage = 0;
  m_Compressor = CompressionLookupTable::New();

  m_CoregisterToPosteriorProbabilities = false;
  for ( int i = 0; i < 3; i++ )
  {
    m_PartialVolumeUpsamplingFactors[ i ] = 1;
  }
  m_StartingMeshNumber = -1;
  m_DownSamplingFactor = 1;
  m_BackgroundPriorThresholdSmoothingSigma = itk::NumericTraits< float >::max();
  m_LogTransformData = false;
  m_BiasFieldOrder = 0;
  m_K = -1.0f;

  m_MeshSmoothingSigmas.push_back( 0.0f );
  m_ImageSmoothingSigmas.push_back( 0.0f );

  m_MultiResolutionLevel = 0;
  m_NumberOfMultiResolutionLevels = 0;


  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  m_WorldToImageTransform = TransformType::New();

  m_TrackingCollection = AtlasMeshCollection::New();

  // Keep an eye on the segmenter, so that intermediate progression can be
  // written out the the log directory
  typedef itk::MemberCommand< Self >   MemberCommandType;
  MemberCommandType::Pointer  command = MemberCommandType::New();
  command->SetCallbackFunction( this, &Self::HandleSegmenterEvent );
  m_Segmenter->AddObserver( PositionUpdatingStartEvent(), command );
  m_Segmenter->AddObserver( PositionUpdatingIterationEvent(), command );
  m_Segmenter->AddObserver( PositionUpdatingEndEvent(), command );

}



//
//
//
AtlasMeshSegmentationDriver
::~AtlasMeshSegmentationDriver()
{

}


//
//
//
void
AtlasMeshSegmentationDriver
::Initialize()
{


  // Make sure the log directory exists and is empty
  if ( itksys::SystemTools::FileIsDirectory( m_LogDirectory.c_str() ) )
  {
    // Logging directory exists already. Remove
    if ( !( itksys::SystemTools::RemoveADirectory( m_LogDirectory.c_str() ) ) )
    {
      itkExceptionMacro( << "Couldn't remove existing log directory " << m_LogDirectory );
    }
  }

  if ( !( itksys::SystemTools::MakeDirectory( m_LogDirectory.c_str() ) ) )
  {
    itkExceptionMacro( "Couldn't create log directory" << m_LogDirectory );
  }
  m_LogDirectory += "/";
  std::cout << "Having m_LogDirectory" << m_LogDirectory << std::endl;

  std::cout << "Having downSamplingFactor of " << m_DownSamplingFactor << std::endl;
  CroppedImageReader::Pointer  reader = CroppedImageReader::New();
  reader->SetExtraFraction( 0.1 );
  if ( m_BoundingFileName.size() == 0 )
  {
    reader->Read( m_ImageFileName.c_str() );
  }
  else
  {
    reader->Read( m_ImageFileName.c_str(), m_BoundingFileName.c_str() );
  }

  // Remember what parts of the original image we're actually analyzing, and the size
  // of the original image.
  m_OriginalImageCropRegion = reader->GetOriginalImageRegion();
  m_OriginalImageOriginalRegion = reader->GetOriginalImageOriginalRegion();

  // Remember the world-to-image transform of the cropped image, where the first voxel
  // is considered having index (0,0,0)^T
  m_WorldToImageTransform->Compose( reader->GetWorldToImageTransform() );

  // Let the segmenter know about our partial volume desires
  m_Segmenter->SetPartialVolumeUpsamplingFactors( m_PartialVolumeUpsamplingFactors );


  // Read the mesh collection
  AtlasMeshCollection::Pointer  meshCollection = AtlasMeshCollection::New();

  if ( m_MeshCollectionFileName.empty() )
  {
    // Generate a mesh from scratch
    const unsigned int  meshSize[] = { 2, 2, 2 };
    const unsigned int  domainSize[] = { m_OriginalImageOriginalRegion.GetSize( 0 ),
                                         m_OriginalImageOriginalRegion.GetSize( 1 ),
                                         m_OriginalImageOriginalRegion.GetSize( 2 )
                                       };
    const float K = 1000.0f;
    const unsigned int numberOfClasses = 1;
    const unsigned int numberOfMeshes = 1;
    const bool  forceBorderVerticesToBackground = false;
    meshCollection->Construct( meshSize, domainSize, K, numberOfClasses, numberOfMeshes, forceBorderVerticesToBackground );

    std::cout << "Created a mesh with " << numberOfClasses << " classes" << std::endl;
    meshCollection->Write( "debugger.txt" );
  }
  else if ( !meshCollection->Read( m_MeshCollectionFileName.c_str() ) )
  {
    itkExceptionMacro( "Couldn't read mesh collection from file " << m_MeshCollectionFileName );
  }

  // Change K if user has specified a value
  if ( m_K > 0 )
  {
    std::cout << "Setting K of mesh collection to: " << m_K << std::endl;
    meshCollection->SetK( m_K );
  }

  // Don't allow deformation of K >= 1000
  if ( meshCollection->GetK() >= 1000 )
  {
    std::cout << "Forbidding deformation" << std::endl;
    m_Segmenter->SetDontDeform( true );
  }

  // Apply the correct transform
  std::cout << "Applying transform: " << std::endl;
  reader->GetTransform()->Print( std::cout );
  meshCollection->Transform( -1, reader->GetTransform() );
  for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
  {
    meshCollection->Transform( i, reader->GetTransform() );
  }

  // Remember this transform in the segmenter, so that it can implement the correct
  // sliding boundary conditions of the mesh deformation
  m_Segmenter->SetMeshToImageTranform( reader->GetTransform() );

  const float  determinant = vnl_det( reader->GetTransform()->GetMatrix().GetVnlMatrix() );
  if ( determinant < 0 )
  {
    std::cout << "Careful here: the applied transformation will turn positive tetrahedra into negative ones." << std::endl;
    std::cout << reader->GetTransform()->GetMatrix().GetVnlMatrix() << std::endl;
    std::cout << " determinant: " << determinant << std::endl;
    std::cout << "Starting to swap the point assignments of each tetrahedron..." << std::endl;

    for ( AtlasMesh::CellsContainer::Iterator  cellIt = meshCollection->GetCells()->Begin();
          cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
      AtlasMesh::CellType*  cell = cellIt.Value();

      if( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
        continue;
      }

      // Swap points assigned to first two vertices. This will readily turn negative tetrahedra
      //into positives ones.
      AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
      const AtlasMesh::PointIdentifier  p0Id = *pit;
      ++pit;
      const AtlasMesh::PointIdentifier  p1Id = *pit;

      pit = cell->PointIdsBegin();
      *pit = p1Id;
      ++pit;
      *pit = p0Id;
    } // End loop over all tetrahedra


    std::cout << "...done!" << std::endl;
  }



  {
    // Write out the image at this point
    typedef AtlasMeshSegmenter::ImageType  ImageType;
    typedef itk::ImageFileWriter< ImageType >  ImageWriterType;
    ImageWriterType::Pointer  imageWriter = ImageWriterType::New();
    imageWriter->SetFileName( ( m_LogDirectory + "originalImageBeingSegmented.mhd" ).c_str() );
    imageWriter->SetInput( reader->GetImage() );
    imageWriter->Update();

    // Write out the mesh at this point
    //AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
    //collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( meshCollection->GetReferenceMesh().GetPointer() ), 1, 1000.0f );
    meshCollection->Write( ( m_LogDirectory + "originalMeshCollectionBeingUsed.txt" ).c_str() );
  }


  // Pass the image to the segmenter, after log transforming, downsampling, and masking it out if necessary
  typedef CroppedImageReader::ImageType  ImageType;
  ImageType::Pointer  image = const_cast< ImageType* >( reader->GetImage() );

  if ( m_LogTransformData )
  {
    // Pump through log transform
    std::cout << "Log-transforming data" << std::endl;
    typedef itk::Image< float, 3 >  LogImageType;
    LogImageType::Pointer  logImage = LogImageType::New();
    logImage->SetRegions( image->GetBufferedRegion() );
    logImage->Allocate();

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
    image = windower->GetOutput();
  }


  if ( m_BackgroundPriorThresholdSmoothingSigma != itk::NumericTraits< float >::max() )
  {
    // Create an empty image that serves as a template for the alpha rasterizor
    typedef AtlasMeshAlphaDrawer::LabelImageType  LabelImageType;
    LabelImageType::Pointer  templateImage = LabelImageType::New();
    templateImage->SetRegions( image->GetLargestPossibleRegion().GetSize() );
    templateImage->Allocate();
    templateImage->FillBuffer( 0 );

    // Now rasterize the prior for the background class
    AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
    alphaDrawer->SetLabelImage( templateImage );
    alphaDrawer->SetLabelNumber( 0 );
    ( const_cast< AtlasMeshAlphaDrawer::AlphaImageType* >( alphaDrawer->GetAlphaImage() ) )->FillBuffer( 1.0 );
    alphaDrawer->Rasterize( meshCollection->GetReferenceMesh() );

    // Smooth it
    typedef itk::DiscreteGaussianImageFilter< AtlasMeshAlphaDrawer::AlphaImageType,
            AtlasMeshAlphaDrawer::AlphaImageType >  SmootherType;
    SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetInput( alphaDrawer->GetAlphaImage() );
    smoother->SetMaximumError( 0.1 );
    smoother->SetUseImageSpacingOff();
    smoother->SetVariance( pow( m_BackgroundPriorThresholdSmoothingSigma, 2 ) );
    smoother->Update();
    AtlasMeshAlphaDrawer::AlphaImageType::ConstPointer  backgroundPrior = smoother->GetOutput();


    // Set all image intensities in areas < probability 1 for background to zero
    itk::ImageRegionIterator< CroppedImageReader::ImageType >  imageIt( image, image->GetLargestPossibleRegion() );
    itk::ImageRegionConstIterator< AtlasMeshAlphaDrawer::AlphaImageType >
    priorIt( backgroundPrior,
             backgroundPrior->GetLargestPossibleRegion() );
    for ( ; !imageIt.IsAtEnd(); ++imageIt, ++priorIt )
    {
      if ( priorIt.Value() > 0.99f )
      {
        imageIt.Value() = 0;
      }

    } // End loop over all voxels

  } // End test if we need to mask out voxels based on background prior


  m_NonDownsampledImage = image;
  if ( m_DownSamplingFactor > 1 )
  {
    // Downsample image
    typedef itk::ShrinkImageFilter< ImageType, ImageType >  ShrinkerType;
    ShrinkerType::Pointer  shrinker = ShrinkerType::New();
    shrinker->SetInput( image );
    shrinker->SetShrinkFactors( m_DownSamplingFactor );
    shrinker->Update();
    image = shrinker->GetOutput();

    // Make sure to unset spacing and origin as we don't look at it,
    // but VTK sure does!
    const float spacing[] = { 1.0f, 1.0f, 1.0f };
    const float origin[] = { 0.0f, 0.0f, 0.0f };
    image->SetSpacing( spacing );
    image->SetOrigin( origin );

    // Apply downsampling also to the mesh collection
    typedef AtlasMeshCollection::TransformType  TransformType;
    TransformType::Pointer  transform = TransformType::New();
    transform->Scale( 1 / static_cast< float >( m_DownSamplingFactor ) );
    meshCollection->Transform( -1, transform );
    for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
    {
      meshCollection->Transform( i, transform );
    }

  }


  m_Segmenter->SetImage( image );

  if ( m_StartingMeshNumber == -1 )
  {
    m_Segmenter->SetMesh( const_cast< AtlasMesh* >( meshCollection->GetReferenceMesh().GetPointer() ) );
  }
  else
  {
    m_Segmenter->SetMesh( const_cast< AtlasMesh* >( meshCollection->GetMesh( m_StartingMeshNumber ) ) );
  }

  // Remember components of the mesh collection
  m_TrackingCollection->SetCells( meshCollection->GetCells() );
  m_TrackingCollection->SetReferencePosition( meshCollection->GetReferencePosition() );
  m_TrackingCollection->SetK( meshCollection->GetK() );


  //
  std::cout << "Setting coregisterToPosteriorProbabilities to " << m_CoregisterToPosteriorProbabilities << std::endl;
  m_Segmenter->SetCoregisterToPosteriorProbabilities( m_CoregisterToPosteriorProbabilities );

  std::cout << "Setting biasFieldOrder to " << m_BiasFieldOrder << std::endl;
  m_Segmenter->SetBiasFieldOrder( m_BiasFieldOrder );

  // Read the lookup table
  if ( !m_Compressor->Read( m_CompressionLookupTableFileName.c_str() ) )
  {
    itkExceptionMacro( << "Coulnd't read the required compression lookup table file " << m_CompressionLookupTableFileName );
  }

  // Retrieve label names corresponding to compressed labels of mesh
  std::cout << "Building labels..." << std::endl;
  m_Labels.clear();
  for ( unsigned int i = 0; i < m_Segmenter->GetNumberOfClasses(); i++ )
  {
    const unsigned int  originalLabel = m_Segmenter->GetEMSegmenter()->GetLookupTable()[ i ];
    CompressionLookupTable::LabelStringLookupTableType::const_iterator it
    = m_Compressor->GetLabelStringLookupTable().find( originalLabel );
    if ( it == m_Compressor->GetLabelStringLookupTable().end() )
    {
      // Not found, which should really not be possible...
      m_Labels.push_back( "No label found" );
    }
    else
    {
      m_Labels.push_back( it->second.c_str() );
    }
  }
  std::cout << "...done! (m_Labels.size(): " << m_Labels.size() << ")" << std::endl;


  // Write out the image you're working with
  typedef AtlasMeshSegmenter::ImageType  ImageType;
  typedef itk::ImageFileWriter< ImageType >  ImageWriterType;
  ImageWriterType::Pointer  imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( ( m_LogDirectory + "imageBeingSegmented.mhd" ).c_str() );
  imageWriter->SetInput( m_Segmenter->GetImage() );
  imageWriter->Update();

  // Write out the mesh
  AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
  collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( m_Segmenter->GetMesh() ), 1, 1000.0f );
  collection->Write( ( m_LogDirectory + "meshBeingUsed.txt" ).c_str() );




  // Create a new, internal labeling where classes that are split into more than one Gaussian are
  // treated as seperate "classes" (i.e., each internal label corresponds to one Gaussian distribution)
  m_InternalToExternalLookupTable.clear();
  CompressionLookupTable::CompressionLookupTableType  compressionLookupTable = m_Compressor->GetCompressionLookupTable();
  // Loop over all classes, each time incrementing the internal label by one unless user specified more than one
  // for the class at hand.
  for ( CompressionLookupTable::CompressionLookupTableType::const_iterator  it = compressionLookupTable.begin();
        it != compressionLookupTable.end(); ++it )
  {
    int  numberOfGaussians = 1;
    SplitClassesContainerType::const_iterator  it2 = m_SplitClasses.find( it->first );
    if ( it2 != m_SplitClasses.end() )
    {
      numberOfGaussians = it2->second;
    }

    for ( int i = 0; i < numberOfGaussians; i++ )
    {
      m_InternalToExternalLookupTable[ m_InternalToExternalLookupTable.size() ] = it->second;
    }
  }

  std::cout << "m_InternalToExternalLookupTable: " << std::endl;
  for ( MoreToLessLabelsLookupTableType::const_iterator it = m_InternalToExternalLookupTable.begin();
        it != m_InternalToExternalLookupTable.end(); ++it )
  {
    const std::string  freeSurferLabelString = m_Compressor->GetLabelStringLookupTable().find( it->second )->second;
    std::cout << "    " << static_cast< int >( it->first )
              << "  -> " << static_cast< int >( it->second )
              << " (" << freeSurferLabelString << ")" << std::endl;
  }



  // Create a new, reduced labeling where internal classes that share the same Gaussian distribution parameters,
  // as specified by the user, are added up to one single class. We do this because it will significantly save
  // the time needed by the model parameter estimation process. Once these parameters have been estimated, we
  // will then pass them on to the internal labeling segmentation procedure.
  m_InternalToReducedLookupTable.clear();
  unsigned char  reducedLabel = 0;
  // First add the ones that are grouped into same Gaussians
  for ( unsigned int i = 0; i < m_SameGaussianParameters.size(); i++, reducedLabel++ )
  {
    for ( unsigned int j = 0; j < m_SameGaussianParameters[ i ].size(); j++ )
    {
      const unsigned int  freeSurferLabel = ( m_SameGaussianParameters[ i ] )[ j ];
      CompressionLookupTable::CompressionLookupTableType::const_iterator  it = compressionLookupTable.find( freeSurferLabel );
      if ( it == compressionLookupTable.end() )
      {
        itkExceptionMacro( "You specified a freeSurferLabel " << static_cast< int >( freeSurferLabel )
                           << " that doesn't exist in the provide lookupTable" );
      }
      const unsigned char  externalLabel = it->second;
      unsigned char  internalLabel = 0;
      bool  foundInternalLabel = false;
      for ( MoreToLessLabelsLookupTableType::const_iterator it2 = m_InternalToExternalLookupTable.begin();
            it2 != m_InternalToExternalLookupTable.end(); ++it2 )
      {
        if ( it2->second == externalLabel )
        {

          if ( foundInternalLabel )
          {
            // Ooops... you found more than one internal label that maps to the same external label that
            // is supposed to be split into two or more classes. Don't know how to do that; bailing out
            itkExceptionMacro( << "You specified freeSurferLabel " <<  static_cast< int >( freeSurferLabel )
                               << " as belonging BOTH to a group with the same Gaussian parameters AND as"
                               " being represented by more than one Gaussian!!!???" );
          }
          internalLabel = it2->first;
          foundInternalLabel= true;

        }
      }

      if ( m_InternalToReducedLookupTable.find( internalLabel ) != m_InternalToReducedLookupTable.end() )
      {
        // This guy already existed! That's not how it's supposed to be...
        itkExceptionMacro( << "You specified freeSurferLabel " <<  static_cast< int >( freeSurferLabel )
                           << " as belonging to several groups with the same Gaussian parameters" );
      }
      m_InternalToReducedLookupTable[ internalLabel ] = reducedLabel;
    }
  }
  // Now add all the rest. Do this by looping over all internal labels, and adding
  // entries only for who's not already there
  for ( MoreToLessLabelsLookupTableType::const_iterator  it = m_InternalToExternalLookupTable.begin();
        it != m_InternalToExternalLookupTable.end(); ++it )
  {
    const unsigned char  internalLabel = it->first;
    if ( m_InternalToReducedLookupTable.find( internalLabel ) == m_InternalToReducedLookupTable.end() )
    {
      // Not found, so insert it
      m_InternalToReducedLookupTable[ internalLabel ] = reducedLabel;
      ++reducedLabel;
    }
  }

  std::cout << "m_InternalToReducedLookupTable: " << std::endl;
  for ( MoreToLessLabelsLookupTableType::const_iterator  it = m_InternalToReducedLookupTable.begin();
        it != m_InternalToReducedLookupTable.end(); ++it )
  {
    const unsigned char  externalLabel = m_InternalToExternalLookupTable[ it->first ];
    const std::string  freeSurferLabelString = ( m_Compressor->GetLabelStringLookupTable().find( externalLabel ) )->second;
    std::cout << "     " << static_cast< int >( it->first )
              << " (" << freeSurferLabelString << ")"
              << " -> " << static_cast< int >( it->second ) << std::endl;
  }


  // Remember multiresolution
  if ( m_MeshSmoothingSigmas.size() != m_ImageSmoothingSigmas.size() )
  {
    itkExceptionMacro( << "size of m_MeshSmoothingSigmas must be the same as size of m_ImageSmoothingSigmas" );
  }
  m_NumberOfMultiResolutionLevels = m_MeshSmoothingSigmas.size();


}


//
//
//
void
AtlasMeshSegmentationDriver
::Segment( bool useAffine )
{

  this->InvokeEvent( MultiResolutionStartEvent() );

  // Save original point parameters
  AtlasMesh::Pointer  mesh = const_cast< AtlasMesh* >( m_Segmenter->GetMesh() );
  AtlasMesh::PointDataContainer::Pointer  originalPointParameters = mesh->GetPointData();

  // Get the internal point parameters and save them for later
  AtlasMesh::ConstPointer  internalMesh = this->GetDistributedMesh( mesh.GetPointer(), m_InternalToExternalLookupTable );
  AtlasMesh::PointDataContainer::ConstPointer  internalPointParameters = internalMesh->GetPointData();

  // Set up segmenter with a reduced mesh
  mesh = const_cast< AtlasMesh* >( this->GetReducedMesh( internalMesh.GetPointer(), m_InternalToReducedLookupTable ).GetPointer() );
  AtlasMesh::PointDataContainer::Pointer  reducedPointParameters = mesh->GetPointData();
  m_Segmenter->SetMesh( mesh );
  mesh = const_cast< AtlasMesh* >( m_Segmenter->GetMesh() );

  // Run segmentation for a number of times, with decreasing smoothing of the mesh
  AtlasMeshSegmenter::ImageType::ConstPointer  unsmoothedImage = m_Segmenter->GetImage();

#if 1
  {
    AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
    collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( m_Segmenter->GetMesh() ), 1, 1000.0f );
    collection->Write( ( m_LogDirectory + "debug_BeforeFirstLevel.txt" ).c_str() );
  }
#endif

  for ( m_MultiResolutionLevel = 0;
        m_MultiResolutionLevel < m_NumberOfMultiResolutionLevels;
        m_MultiResolutionLevel++ )
  {
    const float  sigma = m_MeshSmoothingSigmas[ m_MultiResolutionLevel ];
    const float  imageSigma = m_ImageSmoothingSigmas[ m_MultiResolutionLevel ];

    // Get smoothed alphas
    AtlasMesh::PointDataContainer::Pointer  smoothedPointParameters = reducedPointParameters;
    if ( sigma != 0 )
    {
#if 0
      mesh->SetPointData( reducedPointParameters );
      AtlasMeshCollection::Pointer  originalCollection = AtlasMeshCollection::New();
      originalCollection->GenerateFromSingleMesh( mesh, 1, 1000.0f );
      AtlasMeshSmoother::Pointer  smoother = AtlasMeshSmoother::New();
      smoother->SetMeshCollection( originalCollection );
      smoother->SetSigma( sigma / static_cast< float >( m_DownSamplingFactor ) );
      smoothedPointParameters = smoother->GetSmoothedMeshCollection()->GetPointParameters();
#else
      AtlasMeshCollection::Pointer  collectionToBeSmoothed = AtlasMeshCollection::New();
      collectionToBeSmoothed->SetReferencePosition( m_TrackingCollection->GetReferencePosition() );
      collectionToBeSmoothed->SetPointParameters( reducedPointParameters );

      AtlasMeshSmoother::Pointer  smoother = AtlasMeshSmoother::New();
      smoother->SetMeshCollection( collectionToBeSmoothed );
      smoother->SetSigma( sigma / static_cast< float >( m_DownSamplingFactor ) );
      smoothedPointParameters = smoother->GetSmoothedMeshCollection()->GetPointParameters();

#endif
    }

    // Get smoothed image
    AtlasMeshSegmenter::ImageType::ConstPointer  smoothedImage = unsmoothedImage;
    if ( imageSigma != 0 )
    {
      typedef itk::Image< float, 3 >  FloatImageType;
      typedef itk::CastImageFilter< AtlasMeshSegmenter::ImageType, FloatImageType >   CasterType;
      typedef itk::DiscreteGaussianImageFilter< FloatImageType, FloatImageType >  SmootherType;
      typedef itk::CastImageFilter< FloatImageType, AtlasMeshSegmenter::ImageType >  BackCasterType;

      CasterType::Pointer caster = CasterType::New();
      caster->SetInput( unsmoothedImage );
      SmootherType::Pointer smoother = SmootherType::New();
      smoother->SetInput( caster->GetOutput() );
      smoother->SetMaximumError( 0.1 );
      smoother->SetUseImageSpacingOff();
      smoother->SetVariance( imageSigma * imageSigma );
      BackCasterType::Pointer  backCaster = BackCasterType::New();
      backCaster->SetInput( smoother->GetOutput() );
      backCaster->Update();
      smoothedImage = backCaster->GetOutput();
    }

    // Assign smoothed alphas to the mesh instead of the original alphas
    mesh->SetPointData( smoothedPointParameters );

    // Assign smoothed image to segmenter
    m_Segmenter->SetImage( smoothedImage );

    // Initialize the tracking collection for this level
    m_TrackingCollection->SetPointParameters( smoothedPointParameters );
    //m_TrackingCollection->GetPositions().clear();

    this->InvokeEvent( MultiResolutionIterationEvent() );

    std::cout << "Running for smoothing sigma " << sigma << std::endl;
    m_Segmenter->Segment( useAffine );

  } // End loop over multi-resolution

  // Signal that we're done
  this->InvokeEvent( MultiResolutionEndEvent() );





  {
    // Write out a summary as well
    typedef itk::ImageFileWriter< SummaryImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( this->GetSummaryImage() );
    writer->SetFileName( ( m_LogDirectory + "mappedSummary.mhd" ).c_str() );
    writer->Write();
    writer->SetInput( this->GetSummaryImage( false, true ) );
    writer->SetFileName( ( m_LogDirectory + "mappedSummaryCrisp.mhd" ).c_str() );
    writer->Write();
  }




  // Write out a summary as well
  typedef itk::ImageFileWriter< SummaryImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetInput( this->GetSummaryImage() );
  writer->SetFileName( ( m_LogDirectory + "summary.mhd" ).c_str() );
  writer->Write();
  writer->SetInput( this->GetSummaryImage( false, true ) );
  writer->SetFileName( ( m_LogDirectory + "summaryCrisp.mhd" ).c_str() );
  writer->Write();

  // Write out bias field corrected image
  typedef itk::ImageFileWriter< EMSegmenter::BiasCorrectedImageType >  BiasCorrectedWriterType;
  BiasCorrectedWriterType::Pointer  biasCorrectedWriter = BiasCorrectedWriterType::New();
  biasCorrectedWriter->SetInput( m_Segmenter->GetEMSegmenter()->GetBiasCorrectedImage() );
  biasCorrectedWriter->SetFileName( ( m_LogDirectory + "biasCorrectedImage.mhd" ).c_str() );
  biasCorrectedWriter->Write();

#if 0
  typedef itk::ImageFileWriter< ExponentiatedImageType >  ExponentiatedImageWriterType;
  ExponentiatedImageWriterType::Pointer  exponentiatedWriter = ExponentiatedImageWriterType::New();
  exponentiatedWriter->SetFileName( "exponentiatedBiasCorrectedImage.mhd" );
  exponentiatedWriter->SetInput( AtlasMeshSegmenterConsole::GetExponentiatedImage( m_Segmenter->GetEMSegmenter()->GetBiasCorrectedImage() ) );
  exponentiatedWriter->Write();
#endif

  ImageType::Pointer  nativeImage = 0;
  if ( m_Segmenter->GetEMSegmenter()->GetBiasField() )
  {
    // Write out estimated bias field
    typedef itk::ImageFileWriter< EMSegmenter::BiasFieldImageType >  BiasWriterType;
    BiasWriterType::Pointer  biasWriter = BiasWriterType::New();
    biasWriter->SetInput( m_Segmenter->GetEMSegmenter()->GetBiasField() );
    biasWriter->SetFileName( ( m_LogDirectory + "bias.mhd" ).c_str() );
    biasWriter->Write();


    // Write out non-downsampled bias corrected image
    EMSegmenter::BiasCorrectedImageType::ConstPointer  nonDownsampledBiasCorrectedImage
    = m_Segmenter->GetEMSegmenter()->BiasCorrect( m_NonDownsampledImage, m_DownSamplingFactor ).GetPointer();
    biasCorrectedWriter->SetInput( nonDownsampledBiasCorrectedImage );
    biasCorrectedWriter->SetFileName( ( m_LogDirectory + "nonDownsampledBiasCorrectedImage.mhd" ).c_str() );
    biasCorrectedWriter->Write();

    // Convert it to native pixel type, simply by casting. Make sure we don't exceed range of what can be represented
    nativeImage = ImageType::New();
    nativeImage->SetRegions( nonDownsampledBiasCorrectedImage->GetBufferedRegion() );
    nativeImage->Allocate();
    itk::ImageRegionConstIterator< EMSegmenter::BiasCorrectedImageType >  nonNativeIt( nonDownsampledBiasCorrectedImage,
        nonDownsampledBiasCorrectedImage->GetBufferedRegion() );
    itk::ImageRegionIterator< ImageType >  nativeIt( nativeImage,
        nativeImage->GetBufferedRegion() );
    for ( ; !nonNativeIt.IsAtEnd(); ++nonNativeIt, ++nativeIt )
    {
      if ( ( nonNativeIt.Value() + 0.5 ) > itk::NumericTraits< ImageType::PixelType >::max() )
      {
        nativeIt.Value() = itk::NumericTraits< ImageType::PixelType >::max();
      }
      else if ( ( nonNativeIt.Value() + 0.5 ) < itk::NumericTraits< ImageType::PixelType >::min() )
      {
        nativeIt.Value() = itk::NumericTraits< ImageType::PixelType >::min();
      }
      else
      {
        nativeIt.Value() = static_cast< ImageType::PixelType >( nonNativeIt.Value() + 0.5 );
      }
    }


    // Write out non-downsampled bias corrected image in short format
    typedef itk::ImageFileWriter< ImageType >  NativeWriterType;
    NativeWriterType::Pointer  nativeWriter = NativeWriterType::New();
    nativeWriter->SetInput( nativeImage );
    nativeWriter->SetFileName( ( m_LogDirectory + "nonDownsampledBiasCorrectedImageNonFloat.mhd" ).c_str() );
    nativeWriter->Write();
  }
  else
  {
    nativeImage = m_NonDownsampledImage;
  }


  // Write out the mesh, undoing downsampling
  AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
  collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( m_Segmenter->GetMesh() ), 1, 1000.0f );
#if 0
  collection->SetReferencePosition( m_OriginalMesh->GetPoints() );
#endif

  typedef AtlasMeshCollection::TransformType  TransformType;
  TransformType::Pointer  transform = TransformType::New();
  transform->Scale( m_DownSamplingFactor );
  collection->Transform( -1, transform );
  collection->Transform( 0, transform );
  std::ostringstream  fileNameStream;
  collection->Write( ( m_LogDirectory + "mesh.txt" ).c_str() );


#if 1
  if ( !m_MeshCollectionFileName.empty() )
  {
    std::cout << "Mapping warped mesh back to original space" << std::endl;

    // Read again the original mesh collection in. This is needed because we may have messed around
    // with the ordering of the vertices in each tetrahedron, in order to keep tetrahedra "positive"
    // irrespective of the affine meshToImage transformation applied
    AtlasMeshCollection::Pointer  warpedOriginalMeshCollection = AtlasMeshCollection::New();
    if ( !warpedOriginalMeshCollection->Read( m_MeshCollectionFileName.c_str() ) )
    {
      itkExceptionMacro( "Couldn't read mesh collection from file " << m_MeshCollectionFileName );
    }

    // Replace positions in the original mesh with the ones we have now. Make sure
    // these are actual *copies* of the data, and not just pointers to our valuable
    // data that should be protected
    AtlasMesh::PointsContainer::Pointer  warpedPosition = AtlasMesh::PointsContainer::New();
    for ( AtlasMesh::PointsContainer::ConstIterator  it = collection->GetReferencePosition()->Begin();
          it != collection->GetReferencePosition()->End(); ++it )
    {
      warpedPosition->InsertElement( it.Index(), it.Value() );
    }
    warpedOriginalMeshCollection->SetPositions( std::vector< AtlasMesh::PointsContainer::Pointer >( 1, warpedPosition ) );

    // Undo the mesh to image transform we applied
    TransformType::Pointer  inverseMeshToImagetransform = TransformType::New();
    m_Segmenter->GetMeshToImageTranform()->GetInverse( inverseMeshToImagetransform );
    //warpedOriginalMeshCollection->Transform( -1, inverseMeshToImagetransform );
    warpedOriginalMeshCollection->Transform( 0, inverseMeshToImagetransform );

    // Write out
    warpedOriginalMeshCollection->Write( ( m_LogDirectory + "warpedOriginalMesh.txt" ).c_str() );
  }
#endif


  // Now classify native according to estimated parameters
  std::cout << "Abusing EM segmenter now..." << std::endl;
  EMSegmenter::Pointer  myAbusedEMSegmenter = const_cast< EMSegmenter* >( m_Segmenter->GetEMSegmenter() );
  mesh->SetPointData( const_cast< AtlasMesh::PointDataContainer* >( internalPointParameters.GetPointer() ) );
  mesh->SetPoints( collection->GetReferencePosition() );
  m_Segmenter->SetMesh( mesh );

#if 0
  {
    std::ostringstream  fileNameStream;
    fileNameStream << m_LogDirectory << "internalMesh.txt";

    AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
    collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( m_Segmenter->GetMesh() ), 1, 1000.0f );
    collection->Write( fileNameStream.str().c_str() );
    //exit( -1 );
  }
#endif

  const int  numberOfInternalLabels = m_InternalToReducedLookupTable.size();
  std::vector< float >  internalMeans( numberOfInternalLabels );
  std::vector< float >  internalVariances( numberOfInternalLabels );
  for ( MoreToLessLabelsLookupTableType::const_iterator  it = m_InternalToReducedLookupTable.begin();
        it != m_InternalToReducedLookupTable.end(); ++it )
  {
    internalMeans[ it->first ] = myAbusedEMSegmenter->GetMeans()[ it->second ];
    internalVariances[ it->first ] = myAbusedEMSegmenter->GetVariances()[ it->second ];
  }
  std::vector< float >  internalPriorWeights( numberOfInternalLabels, 1.0f );

  //
  const_cast< std::vector< float >& >( myAbusedEMSegmenter->GetMeans() ) = internalMeans;
  const_cast< std::vector< float >& >( myAbusedEMSegmenter->GetVariances() ) = internalVariances;
  const_cast< std::vector< float >& >( myAbusedEMSegmenter->GetPriorWeights() ) = internalPriorWeights;

  std::cout << "myAbusedEMSegmenter->GetMeans(): [";
  for ( unsigned int i = 0; i < myAbusedEMSegmenter->GetMeans().size(); i++ )
  {
    std::cout << " " << myAbusedEMSegmenter->GetMeans()[ i ];
  }
  std::cout << "]" << std::endl;

  std::cout << "myAbusedEMSegmenter->GetVariances(): [";
  for ( unsigned int i = 0; i < myAbusedEMSegmenter->GetVariances().size(); i++ )
  {
    std::cout << " " << myAbusedEMSegmenter->GetVariances()[ i ];
  }
  std::cout << "]" << std::endl;

  const_cast< std::vector< float >& >( myAbusedEMSegmenter->GetPriorWeights() ) = internalPriorWeights;

  std::cout << "myAbusedEMSegmenter->GetPriorWeights(): [";
  for ( unsigned int i = 0; i < myAbusedEMSegmenter->GetPriorWeights().size(); i++ )
  {
    std::cout << " " << myAbusedEMSegmenter->GetPriorWeights()[ i ];
  }
  std::cout << "]" << std::endl;


  myAbusedEMSegmenter->SetImageAndSegmentItWithCurrentParametersAndWithoutBiasFieldCorrection( nativeImage );

  std::cout << "Still here..." << std::endl;
#if 0
  {
    std::ostringstream  fileNameStream;
    fileNameStream << m_LogDirectory << "internalMeshOfAbusedEMSegmenter.txt";

    AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
    collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( myAbusedEMSegmenter->GetAtlasMesh() ), 1, 1000.0f );
    collection->Write( fileNameStream.str().c_str() );
    //exit( -1 );
  }
#endif
#if 0
  {
    for ( int i = 0; i < myAbusedEMSegmenter->GetNumberOfClasses(); i++ )
    {
      typedef itk::ImageFileWriter< EMSegmenter::ClassificationImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      std::ostringstream  fileNameStream;

      fileNameStream << "myAbusedEMSegmenter_Prior" << i << ".mhd";
      writer->SetFileName( fileNameStream.str() );
      writer->SetInput( myAbusedEMSegmenter->GetPrior( i) );
      writer->Write();

      fileNameStream.clear();
      fileNameStream << "myAbusedEMSegmenter_Posterior" << i << ".mhd";
      writer->SetFileName( fileNameStream.str() );
      writer->SetInput( myAbusedEMSegmenter->GetPosterior( i) );
      writer->Write();

    }


  }
#endif





  // Write out posteriors
  TransformType::Pointer  posteriorImageToWorldTransform = TransformType::New();
  m_WorldToImageTransform->GetInverse( posteriorImageToWorldTransform );
  ImageType::SpacingType  posteriorSpacing;
  ImageType::PointType  posteriorOrigin;
  ImageType::DirectionType  posteriorDirection;
  GetSpacingOriginAndDirection( posteriorImageToWorldTransform,
                                posteriorSpacing,
                                posteriorOrigin,
                                posteriorDirection );
  for ( unsigned int i = 0; i < myAbusedEMSegmenter->GetNumberOfClasses(); i++ )
  {
    // Allocate empty image
    typedef  itk::Image< unsigned char, 3 >  MyImageType;
    MyImageType::Pointer  myImage = MyImageType::New();
    myImage->SetRegions( myAbusedEMSegmenter->GetPosterior( i )->GetLargestPossibleRegion().GetSize() );
    myImage->Allocate();

    // Set spacing, origin, and direction
    myImage->SetSpacing( posteriorSpacing );
    myImage->SetOrigin( posteriorOrigin );
    myImage->SetDirection( posteriorDirection );

    // Fill in
    itk::ImageRegionConstIterator< EMSegmenter::ClassificationImageType >
    sourceIt( myAbusedEMSegmenter->GetPosterior( i ),
              myAbusedEMSegmenter->GetPosterior( i )->GetBufferedRegion() );
    itk::ImageRegionIterator< MyImageType >  targetIt( myImage,
        myImage->GetBufferedRegion() );
    for ( ; !sourceIt.IsAtEnd(); ++sourceIt, ++targetIt )
    {
      targetIt.Value() = static_cast< MyImageType::PixelType >( 255.0f * sourceIt.Value() + 0.5  );
    }


    // Write out
    const unsigned char  externalLabel = m_InternalToExternalLookupTable[ i ];
    const std::string  freeSurferLabelString = ( m_Compressor->GetLabelStringLookupTable().find( externalLabel ) )->second;


    typedef itk::ImageFileWriter< MyImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    std::ostringstream  fileNameStream;

    fileNameStream << m_LogDirectory << "posterior_" << freeSurferLabelString << ".mgz";
    writer->SetFileName( fileNameStream.str() );
    writer->SetInput( myImage );
    writer->Write();
  }


  if ( myAbusedEMSegmenter->GetSuperResolutionPosteriors() )
  {
    std::cout << "Having super resolution posterior!!!" << std::endl;


    // Determine spacing, origin, and direction of superresolution images
    // w.r.t. original image
    TransformType::ParametersType  imageToSuperResolutionImageTransformParams( 12 );
    imageToSuperResolutionImageTransformParams.Fill( 0.0f );
    for ( int i = 0; i < 3; i++ )
    {
      // Diagonal (scaling)
      imageToSuperResolutionImageTransformParams( i + i * 3 ) = m_PartialVolumeUpsamplingFactors[ i ];

      // Translation
      imageToSuperResolutionImageTransformParams( 9 + i ) = ( m_PartialVolumeUpsamplingFactors[ i ] - 1.0f ) / 2.0f;
    }
    TransformType::Pointer  imageToSuperResolutionImageTransform = TransformType::New();
    imageToSuperResolutionImageTransform->SetParameters( imageToSuperResolutionImageTransformParams );
    std::cout << "imageToSuperResolutionImageTransform: " << std::endl;
    imageToSuperResolutionImageTransform->Print( std::cout );

    TransformType::Pointer  superResolutionWorldToImageTransform = TransformType::New();
    superResolutionWorldToImageTransform->Compose( m_WorldToImageTransform );
    superResolutionWorldToImageTransform->Compose( imageToSuperResolutionImageTransform );

    TransformType::Pointer  superResolutionImageToWorldTransform = TransformType::New();
    superResolutionWorldToImageTransform->GetInverse( superResolutionImageToWorldTransform );
    ImageType::SpacingType  superResolutionSpacing;
    ImageType::PointType  superResolutionOrigin;
    ImageType::DirectionType  superResolutionDirection;
    GetSpacingOriginAndDirection( superResolutionImageToWorldTransform,
                                  superResolutionSpacing,
                                  superResolutionOrigin,
                                  superResolutionDirection );


    // Write out super resolution image
    typedef  itk::ImageFileWriter< EMSegmenter::ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( myAbusedEMSegmenter->GetSuperResolutionImage() );
    writer->SetFileName( ( m_LogDirectory + "superResolutionImage.mhd" ).c_str() );
    writer->Update();

    // Retrieve the super resolution posterior
    EMSegmenter::SuperResolutionProbabilityImageType::ConstPointer  superResolutionPosteriors =
      myAbusedEMSegmenter->GetSuperResolutionPosteriors().GetPointer();

    // Split up into different classes and write out
    for ( unsigned int i = 0; i < myAbusedEMSegmenter->GetNumberOfClasses(); i++ )
    {
      // Allocate empty image
      typedef  itk::Image< unsigned char, 3 >  MyImageType;
      MyImageType::Pointer  myImage = MyImageType::New();
      myImage->SetRegions( superResolutionPosteriors->GetLargestPossibleRegion().GetSize() );
      myImage->Allocate();

      // Set spacing, origin, and direction
      myImage->SetSpacing( superResolutionSpacing );
      myImage->SetOrigin( superResolutionOrigin );
      myImage->SetDirection( superResolutionDirection );

      // Fill in
      itk::ImageRegionConstIterator< EMSegmenter::SuperResolutionProbabilityImageType >
      sourceIt( superResolutionPosteriors,
                superResolutionPosteriors->GetBufferedRegion() );
      itk::ImageRegionIterator< MyImageType >  targetIt( myImage,
          myImage->GetBufferedRegion() );
      for ( ; !sourceIt.IsAtEnd(); ++sourceIt, ++targetIt )
      {
        targetIt.Value() = static_cast< MyImageType::PixelType >( 255.0f * sourceIt.Value()[ i ] + 0.5  );
      }


      // Write out
      const unsigned char  externalLabel = m_InternalToExternalLookupTable[ i ];
      const std::string  freeSurferLabelString = ( m_Compressor->GetLabelStringLookupTable().find( externalLabel ) )->second;


      typedef itk::ImageFileWriter< MyImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      std::ostringstream  fileNameStream;

      fileNameStream << m_LogDirectory << "superResolutionPosterior_" << freeSurferLabelString << ".mgz";
      writer->SetFileName( fileNameStream.str() );
      writer->SetInput( myImage );
      writer->Write();
    }

  } // End test if we have super resolution posterior



  writer->SetInput( this->GetSummaryImage() );
  writer->SetFileName( ( m_LogDirectory + "nonDownsampledSummary.mhd" ).c_str() );
  writer->Write();
  writer->SetInput( this->GetSummaryImage( false, true ) );
  writer->SetFileName( ( m_LogDirectory + "nonDownsampledSummaryCrisp.mhd" ).c_str() );
  writer->Write();
  writer->SetInput( this->GetSummaryImage( false, true, false ) );
  writer->SetFileName( ( m_LogDirectory + "nonDownsampledSummaryCrispNotDemangled.mhd" ).c_str() );
  writer->Write();


#if 0
  //
  std::cout << "Writing FreeSurfer image..." << std::flush;
  FreeSurferLabelImageType::Pointer  freeSurferImage
  = this->MapCrispSummaryBackToFreeSurfer( this->GetSummaryImage( false, true ) );
  typedef itk::ImageFileWriter< FreeSurferLabelImageType >  FreeSurferWriterType;
  FreeSurferWriterType::Pointer  freeSurferWriter = FreeSurferWriterType::New();
  freeSurferWriter->SetInput( freeSurferImage );
  freeSurferWriter->SetFileName( ( m_LogDirectory + "freeSurfer.mhd" ).c_str() );
  freeSurferWriter->Update();
#endif
  std::cout << "done!" << std::endl;


}



//
//
//
AtlasMeshSegmentationDriver::SummaryImageType::Pointer
AtlasMeshSegmentationDriver
::GetSummaryImage( bool  usePrior, bool  useCrisp, bool joinSplitClasses, const std::vector< float >&  means ) const
{

  // Create an empty image to hold the reconstruction
  SummaryImageType::Pointer  summary = SummaryImageType::New();
  summary->SetRegions( m_Segmenter->GetEMSegmenter()->GetPosterior( 0 )->GetLargestPossibleRegion() );
  summary->SetSpacing( m_Segmenter->GetEMSegmenter()->GetPosterior( 0 )->GetSpacing() );
  summary->SetOrigin( m_Segmenter->GetEMSegmenter()->GetPosterior( 0 )->GetOrigin() );
  summary->Allocate();
  summary->FillBuffer( 0 );

  // Create an empty image to hold the maximum probability in each voxel
  SummaryImageType::Pointer  maximumProbability = SummaryImageType::New();
  maximumProbability->SetRegions( m_Segmenter->GetEMSegmenter()->GetPosterior( 0 )->GetLargestPossibleRegion() );
  maximumProbability->SetSpacing( m_Segmenter->GetEMSegmenter()->GetPosterior( 0 )->GetSpacing() );
  maximumProbability->SetOrigin( m_Segmenter->GetEMSegmenter()->GetPosterior( 0 )->GetOrigin() );
  maximumProbability->Allocate();
  maximumProbability->FillBuffer( 0 );

  // Loop over all the posterior images
  for ( unsigned int  classNumber = 0; classNumber < m_Segmenter->GetNumberOfClasses(); classNumber++ )
  {
    // Decide what mean to use
    float  mean;
    if ( means.size() != 0 )
    {
      mean = means[ classNumber ];
    }
    else
    {
      if ( joinSplitClasses )
      {
        //mean = m_Segmenter->GetEMSegmenter()->GetLookupTable()[ classNumber ];
        mean = m_InternalToExternalLookupTable.find( classNumber )->second;
      }
      else
      {
        mean = classNumber;
      }
    }

    // Loop over all voxels
    EMSegmenter::ClassificationImageType::ConstPointer  classification = 0;
    if ( usePrior )
    {
      classification =  m_Segmenter->GetEMSegmenter()->GetPrior( classNumber );
    }
    else
    {
      classification =  m_Segmenter->GetEMSegmenter()->GetPosterior( classNumber );
    }
    itk::ImageRegionConstIterator< EMSegmenter::ClassificationImageType >
    classificationIt( classification,
                      classification->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< SummaryImageType >  maximumIt( maximumProbability,
        maximumProbability->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< SummaryImageType >  summaryIt( summary,
        summary->GetLargestPossibleRegion() );

    for ( ; !classificationIt.IsAtEnd(); ++classificationIt, ++maximumIt, ++summaryIt )
    {
      if ( useCrisp )
      {
        // Only override voxel value if this is the MAP so far
        if ( classificationIt.Value() > maximumIt.Value() )
        {
          maximumIt.Value() = classificationIt.Value();
          summaryIt.Value() = mean;
        }

      }
      else
      {
        summaryIt.Value() += mean * classificationIt.Value();
      }


    } // End loop over all voxels

  } // End loop over all labels


  return summary;


}


//
//
//
AtlasMeshSegmentationDriver::FreeSurferLabelImageType::Pointer
AtlasMeshSegmentationDriver
::MapCrispSummaryBackToFreeSurfer( const SummaryImageType*  crispSummary ) const
{

  // Sanity check
  if ( m_OriginalImageCropRegion.GetSize() !=  crispSummary->GetLargestPossibleRegion().GetSize() )
  {
    itkExceptionMacro( << "MapCrispSummaryBackToFreeSurfer failed because of image size mismatch" );
  }


  // The compressor tells how to map from freeSurfer to compressed. Create an inverse mapping
  // for going the other way around
  std::map< CompressionLookupTable::CompressedImageType::PixelType,
      CompressionLookupTable::ImageType::PixelType >   inverseLookupTable;
  for ( CompressionLookupTable::CompressionLookupTableType::const_iterator  it = m_Compressor->GetCompressionLookupTable().begin();
        it != m_Compressor->GetCompressionLookupTable().end(); ++it )
  {
    inverseLookupTable[ it->second ] = it->first;
  }


  // Create an empty image, and fill it up
  FreeSurferLabelImageType::Pointer  freeSurferImage = FreeSurferLabelImageType::New();
  freeSurferImage->SetRegions( m_OriginalImageOriginalRegion );
  freeSurferImage->Allocate();
  freeSurferImage->FillBuffer( 0 );

  itk::ImageRegionConstIterator< SummaryImageType >  summaryIt( crispSummary, crispSummary->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< FreeSurferLabelImageType >  freeSurferIt( freeSurferImage, m_OriginalImageCropRegion );
  for ( ; !summaryIt.IsAtEnd(); ++summaryIt, ++freeSurferIt )
  {
    freeSurferIt.Value() = inverseLookupTable[ static_cast< CompressionLookupTable::CompressedImageType::PixelType >( summaryIt.Value() ) ];
  }

  return freeSurferImage;
}


//
//
//
void
AtlasMeshSegmentationDriver
::ParseSetUpFile( const std::string& setUpFileName )
{

  // Open the file for reading
  std::ifstream  in( setUpFileName.c_str() );
  if ( in.bad() )
  {
    itkExceptionMacro( << "Couldn't open " << setUpFileName << " for reading!" );
  }

  // Loop over all lines, and fill in the correct field
  const int size = 255;
  char buffer[ size ];
  while ( in.getline( buffer, size ) )
  {
    if ( ( buffer[0] == '#' ) || ( buffer[0] == 10 ) || ( buffer[0] == 13 ) ) // '\n' corresponds to 10 or 13, depending on platform
    {
      //std::cout << "skipping the line: " << buffer << std::endl;
      continue;
    }

    const std::string  line = buffer;

    // Look for the colon symbol, and split the line up into a key and a value string
    std::size_t foundPosition;
    if ( ( foundPosition = line.find( ":" ) ) == std::string::npos )
    {
      //std::cout << "skipping the line: " << line << std::endl;
      continue;
    }
    std::string  key( line, 0, foundPosition );
    std::string  value( line, foundPosition+1 );

    // Remove white spaces from key
    std::istringstream  keyStream( key );
    keyStream >> key;

    // Depending on the key, do the Right Thing
    if ( !key.compare( "logDirectory" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_LogDirectory;
      std::cout << "m_LogDirectory: " << m_LogDirectory << std::endl;
    }
    else if ( !key.compare( "imageFileName" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_ImageFileName;
      std::cout << "m_ImageFileName: " << m_ImageFileName << std::endl;
    }
    else if ( !key.compare( "meshCollectionFileName" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_MeshCollectionFileName;
      std::cout << "m_MeshCollectionFileName: " << m_MeshCollectionFileName << std::endl;
    }
    else if ( !key.compare( "startingMeshNumber" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_StartingMeshNumber;
      std::cout << "m_StartingMeshNumber: " << m_StartingMeshNumber << std::endl;
    }
    else if ( !key.compare( "boundingFileName" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_BoundingFileName;
      std::cout << "m_BoundingFileName: " << m_BoundingFileName << std::endl;
    }
    else if ( !key.compare( "compressionLookupTableFileName" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_CompressionLookupTableFileName;
      std::cout << "m_CompressionLookupTableFileName: " << m_CompressionLookupTableFileName << std::endl;
    }
    else if ( !key.compare( "coregisterToPosteriorProbabilities" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_CoregisterToPosteriorProbabilities;
      std::cout << "m_CoregisterToPosteriorProbabilities: " << m_CoregisterToPosteriorProbabilities << std::endl;
    }
    else if ( !key.compare( "partialVolumeUpsamplingFactors" ) )
    {
      std::istringstream  valueStream( value );
      std::cout << "m_PartialVolumeUpsamplingFactors: [";
      for ( int i = 0; i < 3; i++ )
      {
        valueStream >> m_PartialVolumeUpsamplingFactors[ i ];
        std::cout << " " << m_PartialVolumeUpsamplingFactors[ i ];
      }
      std::cout << "]" << std::endl;
    }
    else if ( !key.compare( "downSamplingFactor" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_DownSamplingFactor;
      std::cout << "m_DownSamplingFactor: " << m_DownSamplingFactor << std::endl;
    }
    else if ( !key.compare( "backgroundPriorThresholdSmoothingSigma" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_BackgroundPriorThresholdSmoothingSigma;
      std::cout << "m_BackgroundPriorThresholdSmoothingSigma: " << m_BackgroundPriorThresholdSmoothingSigma << std::endl;
    }
    else if ( !key.compare( "logTransformData" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_LogTransformData;
      std::cout << "m_LogTransformData: " << m_LogTransformData << std::endl;
    }
    else if ( !key.compare( "biasFieldOrder" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_BiasFieldOrder;
      std::cout << "m_BiasFieldOrder: " << m_BiasFieldOrder << std::endl;
    }
    else if ( !key.compare( "K" ) )
    {
      std::istringstream  valueStream( value );
      valueStream >> m_K;
      std::cout << "m_K: " << m_K << std::endl;
    }
    else if ( !key.compare( "splitClass" ) )
    {
      std::istringstream  valueStream( value );
      unsigned int  label;
      int  numberOfClasses;
      valueStream >> label >> numberOfClasses;
      m_SplitClasses[ label ] = numberOfClasses;
      std::cout << "m_SplitClasses[ " << label << " ]: " << m_SplitClasses[ label ] << std::endl;
    }
    else if ( !key.compare( "sameGaussianParameters" ) )
    {
      std::istringstream  valueStream( value );
      unsigned int  label;
      std::vector< unsigned int >  bag;
      while ( valueStream >> label )
      {
        bag.push_back( label );
      }
      m_SameGaussianParameters.push_back( bag );

      std::cout << "m_SameGaussianParameters[ " << m_SameGaussianParameters.size() - 1 << " ]: ";
      for ( unsigned int i = 0; i < ( m_SameGaussianParameters[ m_SameGaussianParameters.size() - 1 ] ).size(); i++ )
      {
        std::cout << " " << ( m_SameGaussianParameters[ m_SameGaussianParameters.size() - 1 ] )[ i ];
      }
      std::cout << std::endl;
    }
    else if ( !key.compare( "meshSmoothingSigmas" ) )
    {
      std::istringstream  valueStream( value );
      float  sigma;
      m_MeshSmoothingSigmas.clear();
      while ( valueStream >> sigma )
      {
        m_MeshSmoothingSigmas.push_back( sigma );
      }

      std::cout << "m_MeshSmoothingSigmas: ";
      for ( SmoothingSigmasContainterType::const_iterator  it = m_MeshSmoothingSigmas.begin();
            it != m_MeshSmoothingSigmas.end(); ++it )
      {
        std::cout << " " << *it;
      }
      std::cout << std::endl;
    }
    else if ( !key.compare( "imageSmoothingSigmas" ) )
    {
      std::istringstream  valueStream( value );
      float  sigma;
      m_ImageSmoothingSigmas.clear();
      while ( valueStream >> sigma )
      {
        m_ImageSmoothingSigmas.push_back( sigma );
      }

      std::cout << "m_ImageSmoothingSigmas: ";
      for ( SmoothingSigmasContainterType::const_iterator  it = m_ImageSmoothingSigmas.begin();
            it != m_ImageSmoothingSigmas.end(); ++it )
      {
        std::cout << " " << *it;
      }
      std::cout << std::endl;
    }

  } // End loop over all lines


}


//
//
//
AtlasMesh::ConstPointer
AtlasMeshSegmentationDriver
::GetDistributedMesh( const AtlasMesh*  originalMesh,
                      const MoreToLessLabelsLookupTableType& lookupTable ) const
{

  // Pre-calculate in how much classes each atlas prior will be divided
  std::map< unsigned char, int >  numberOfDistributedLabelsPerOriginalLabel;
  for ( MoreToLessLabelsLookupTableType::const_iterator  it = lookupTable.begin();
        it != lookupTable.end(); ++it )
  {
    std::map< unsigned char, int >::iterator  it2 = numberOfDistributedLabelsPerOriginalLabel.find( it->second );
    if ( it2 == numberOfDistributedLabelsPerOriginalLabel.end() )
    {
      numberOfDistributedLabelsPerOriginalLabel[ it->second ] = 1;
    }
    else
    {
      ( it2->second )++;
    }
  }

  for ( std::map< unsigned char, int >::iterator  it = numberOfDistributedLabelsPerOriginalLabel.begin();
        it != numberOfDistributedLabelsPerOriginalLabel.end(); ++it )
  {
    std::cout << "original label " << static_cast< int >( it->first )
              << " is distributed over " << it->second << " distributed labels" << std::endl;
  }



  // Create a new, distributed point parameters container in which the original alphas
  // are distributed properly:
  // First copy the original point parameters
  AtlasMesh::PointDataContainer::Pointer  distributedPointParameters = AtlasMesh::PointDataContainer::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  it = originalMesh->GetPointData()->Begin();
        it != originalMesh->GetPointData()->End(); ++it )
  {
    distributedPointParameters->InsertElement( it.Index(), it.Value() );
  }

  // Now adjust the alphas in each vertex
  for ( AtlasMesh::PointDataContainer::Iterator  it = distributedPointParameters->Begin();
        it != distributedPointParameters->End(); ++it )
  {
    // Retrieve the original alphas
    AtlasAlphasType  originalAlphas = it.Value().m_Alphas;

    // Loop over all classes in the final configuration
    AtlasAlphasType  finalAlphas( lookupTable.size() );
    for ( MoreToLessLabelsLookupTableType::const_iterator  it2 = lookupTable.begin();
          it2 != lookupTable.end(); ++it2 )
    {
      const int  distributionFactor = numberOfDistributedLabelsPerOriginalLabel[ it2->second ];
      finalAlphas[ it2->first ] = originalAlphas[ it2->second ] /
                                  static_cast< float >( distributionFactor );
    }

    // Set the alpha in this vertex to the final alphas
    it.Value().m_Alphas = finalAlphas;
  }


  // Construct the final atlas mesh from its components
  AtlasMesh::Pointer  distributedMesh = AtlasMesh::New();
  distributedMesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( originalMesh->GetPoints() ) );
  distributedMesh->SetCells(
    const_cast< AtlasMesh::CellsContainer* >( originalMesh->GetCells() ) );
  distributedMesh->SetPointData( distributedPointParameters );
  distributedMesh->SetCellData(
    const_cast< AtlasMesh::CellDataContainer* >( originalMesh->GetCellData() ) );

  return distributedMesh.GetPointer();

}



//
//
//
AtlasMesh::ConstPointer
AtlasMeshSegmentationDriver
::GetReducedMesh( const AtlasMesh*  originalMesh,
                  const MoreToLessLabelsLookupTableType& lookupTable ) const
{

  // Determine number of labels in the reduced mesh
  int  maximalReducedLabel = 0;
  for ( MoreToLessLabelsLookupTableType::const_iterator  it = lookupTable.begin();
        it != lookupTable.end(); ++it )
  {
    if ( it->second > maximalReducedLabel )
    {
      maximalReducedLabel = it->second;
    }
  }
  const int  numberOfReducedLabels = maximalReducedLabel + 1;
  std::cout << "numberOfReducedLabels: " << numberOfReducedLabels << std::endl;


  // Create a new, reduced point parameters container in which the original alphas
  // are summed up properly:
  // First copy the original point parameters
  AtlasMesh::PointDataContainer::Pointer  reducedPointParameters = AtlasMesh::PointDataContainer::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  it = originalMesh->GetPointData()->Begin();
        it != originalMesh->GetPointData()->End(); ++it )
  {
    reducedPointParameters->InsertElement( it.Index(), it.Value() );
  }

  // Now adjust the alphas in each vertex
  for ( AtlasMesh::PointDataContainer::Iterator  it = reducedPointParameters->Begin();
        it != reducedPointParameters->End(); ++it )
  {
    // Retrieve the original alphas
    AtlasAlphasType  originalAlphas = it.Value().m_Alphas;

    // Loop over all labels in the original configuration, and add their contribution
    // to the correct reduced label
    AtlasAlphasType  finalAlphas( numberOfReducedLabels );
    finalAlphas.Fill( 0.0f );
    for ( unsigned char  originalLabel = 0; originalLabel < originalAlphas.Size(); originalLabel++ )
    {
      finalAlphas[ lookupTable.find( originalLabel )->second ] += originalAlphas[ originalLabel ];
    }

    // Set the alpha in this vertex to the final alphas
    it.Value().m_Alphas = finalAlphas;
  }


  // Construct the final atlas mesh from its components
  AtlasMesh::Pointer  reducedMesh = AtlasMesh::New();
  reducedMesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( originalMesh->GetPoints() ) );
  reducedMesh->SetCells(
    const_cast< AtlasMesh::CellsContainer* >( originalMesh->GetCells() ) );
  reducedMesh->SetPointData( reducedPointParameters );
  reducedMesh->SetCellData(
    const_cast< AtlasMesh::CellDataContainer* >( originalMesh->GetCellData() ) );

  return reducedMesh.GetPointer();

}


//
//
//
void
AtlasMeshSegmentationDriver
::GetSpacingOriginAndDirection( const TransformType* transform,
                                ImageType::SpacingType& spacing,
                                ImageType::PointType& origin,
                                ImageType::DirectionType& direction )
{

  //
  std::cout << "Splitting affine matrix up into spacing, origin, and direction components" << std::endl;
  transform->Print( std::cout );

  TransformType::ParametersType  parameters = transform->GetParameters();

  origin[ 0 ] = parameters[ 9 ];
  origin[ 1 ] = parameters[ 10 ];
  origin[ 2 ] = parameters[ 11 ];

  for ( int i = 0; i < 3; i++ )
  {
    spacing[ i ] = sqrt( pow( parameters[ i ], 2 ) + pow( parameters[ i + 3 ], 2 ) + pow( parameters[ i + 6 ], 2 ) );

    direction( 0 , i ) =  parameters[ i ] / spacing[ i ];
    direction( 1 , i ) =  parameters[ i + 3 ] / spacing[ i ];
    direction( 2 , i ) =  parameters[ i + 6 ] / spacing[ i ];
  }

  std::cout << "spacing: " << spacing << std::endl;
  std::cout << "origin: " << origin << std::endl;
  std::cout << "direction: " << direction << std::endl;

}


//
//
//
void
AtlasMeshSegmentationDriver
::HandleSegmenterEvent( itk::Object* object, const itk::EventObject & event )
{

  // Clear the tracking collection if a new position updating round is starting
  if ( typeid( event ) == typeid( PositionUpdatingStartEvent ) )
  {
    // Clear the tracking collection
    //std::cout << "PositionUpdatingStartEvent: clearing m_TrackingCollection" << std::endl;
    m_TrackingCollection->GetPositions().clear();
  }


  // Add the current position to what we have already. We're not testing what
  // event is occurring because we're only listening to PositionUpdatingStartEvent,
  // PositionUpdatingIterationEvent, and PositionUpdatingEndEvent only (see constructor)
  AtlasMesh::PointsContainer::ConstPointer  source = m_Segmenter->GetMesh()->GetPoints();
  AtlasMesh::PointsContainer::Pointer  target = AtlasMesh::PointsContainer::New();
  for ( AtlasMesh::PointsContainer::ConstIterator  sourceIt = source->Begin();
        sourceIt != source->End(); ++sourceIt )
  {
    target->InsertElement( sourceIt.Index(), sourceIt.Value() );
  }
  m_TrackingCollection->GetPositions().push_back( target );
  //std::cout << "Writing m_TrackingCollection to m_TrackingCollection.txt" << std::endl;
  //m_TrackingCollection->Write( "m_TrackingCollection.txt" );


  // Write out the tracking collection if a position updating round is ending
  if ( typeid( event ) == typeid( PositionUpdatingEndEvent ) )
  {
    // Write the tracking collection out
    std::ostringstream  fileNameStream;
    fileNameStream << m_LogDirectory << "evolutionOfMeshNodePositions"
                   << "_MultiResolutionLevel" << m_MultiResolutionLevel
                   << "_Iteration";
    if ( m_Segmenter->GetIterationNumber() < 100 )
    {
      fileNameStream << "0";
      if ( m_Segmenter->GetIterationNumber() < 10 )
      {
        fileNameStream << "0";
      }
    }
    fileNameStream << m_Segmenter->GetIterationNumber()
                   << ".txt";
    //std::cout << "PositionUpdatingEndEvent: writing m_TrackingCollection to "
    //          << fileNameStream.str() << std::endl;
    m_TrackingCollection->Write( fileNameStream.str().c_str() );
  }



}



} // end namespace kvl
