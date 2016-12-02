#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkShrinkImageFilter.h"
#include "kvlEMSegmenter.h"


namespace kvl
{

class BiasFieldCorrectImage : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef BiasFieldCorrectImage         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( BiasFieldCorrectImage, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    std::cout << "I am " << this->GetNameOfClass() 
              << " and I'm running! " << std::endl;
              
              
    // correctedImage = kvlBiasFieldCorrectImage( images, mesh, biasFieldOrder, downSamplingFactor ) //Note that the bias correction in the C++ code corrects for each channel independently, i.e., diagonal covariances
 
    // Some typedefs
    typedef itk::Image< unsigned short, 3 >   ImageType;

              
    // Make sure input arguments are correct
    if ( ( nrhs < 3 ) || !mxIsInt64( prhs[ 0 ] ) || !mxIsInt64( prhs[ 1 ] ) || 
         !mxIsDouble( prhs[ 2 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input image(s)
    typedef itk::Image< unsigned short, 3 >  ImageType;
    //typedef AtlasMeshDeformationConjugateGradientOptimizer::ImageType  ImageType;

    const int  N = mxGetN( prhs[ 1 ] );
    const int  M = mxGetM( prhs[ 1 ] );
    int numberOfImages = 0;
    
    if(N<M)
    {
      numberOfImages = M;
    }
    else
    {
      numberOfImages = N;
    }

    mexPrintf("numberOfImages = %d\n",numberOfImages);
    std::vector<ImageType::Pointer> images;
    uint64_T *  imagesHandle =  static_cast< uint64_T * >( mxGetData( prhs[ 1 ] ) );

    for(unsigned int nima = 0; nima < numberOfImages; nima++, imagesHandle++)
       {
	 
	 const int handle = *(imagesHandle);
	 std::cout<<"Image: "<<handle<<std::endl;
	 itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( handle );
	 if ( typeid(*(object)  ) != typeid( ImageType ) )
	   {
	     mexErrMsgTxt( "image doesn't refer to the correct ITK object type" );
	   }
	 ImageType::ConstPointer constImage = static_cast< const ImageType* >( object.GetPointer() );
         ImageType::Pointer image = const_cast< ImageType* >( constImage.GetPointer() );
         images.push_back(image); 

       }
    
    // Get the mesh
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 1 ] ) ) );
    itk::Object::ConstPointer meshObject = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    if ( typeid( *meshObject ) != typeid( AtlasMesh ) )
      {
      mexErrMsgTxt( "mesh is not correct ITK object type" );
      }
    AtlasMesh::ConstPointer  mesh
            = static_cast< const AtlasMesh* >( meshObject.GetPointer() );
    
    // Get the biasFieldOrder
    const int biasFieldOrder = static_cast< int >( *( mxGetPr( prhs[ 2 ] ) ) );
    
    // Get the downSamplingFactor
    int downSamplingFactor = 1;
    if ( nrhs > 3 )
      {
      downSamplingFactor = static_cast< int >( *( mxGetPr( prhs[ 3 ] ) ) );
      }

    // Show what we have so far
    std::cout << "image: " << images[0] << std::endl;
    std::cout << "mesh: " << mesh << std::endl;
    std::cout << "biasFieldOrder: " << biasFieldOrder << std::endl;
    std::cout << "downSamplingFactor: " << downSamplingFactor << std::endl;
    
    
    // Do the actual work in ITK
    AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
    collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( mesh.GetPointer() ), 1, 1000.0f );
 
    // Downsample the image and the mesh
    std::vector<ImageType::Pointer>  nonDownsampledImages = images;
    if ( downSamplingFactor > 1 )
      {
      // Downsample image
      for(unsigned int nima = 0; nima < numberOfImages; nima++)
      {
        typedef itk::ShrinkImageFilter< ImageType, ImageType >  ShrinkerType;
        ShrinkerType::Pointer  shrinker = ShrinkerType::New();
        shrinker->SetInput( images[nima] );
        shrinker->SetShrinkFactors( downSamplingFactor );
        shrinker->Update();
        images[nima] = shrinker->GetOutput();
      }

      // Apply downsampling also to the mesh
      AtlasMeshCollection::Pointer  meshCollection = kvl::AtlasMeshCollection::New();
      meshCollection->GenerateFromSingleMesh(
            const_cast< kvl::AtlasMesh* >( mesh.GetPointer() ), 1, 1000.0f );
      typedef AtlasMeshCollection::TransformType  TransformType;
      TransformType::Pointer  transform = TransformType::New();
      transform->Scale( 1 / static_cast< float >( downSamplingFactor ) );
      meshCollection->Transform( -1, transform );
      for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
        {
        meshCollection->Transform( i, transform );
        }
      mesh = meshCollection->GetMesh( 0 );
      }


    // Do the work
    EMSegmenter::Pointer  segmenter = EMSegmenter::New();
    segmenter->SetImages( images );
    segmenter->SetAtlasMesh( mesh );
    segmenter->SetBiasFieldOrder( biasFieldOrder );
    //segmenter->SetStopCriterion( 1e-10 );
    //segmenter->SetMaximumNumberOfIterations( 300 );
    segmenter->Segment();


    std::vector<ImageType::Pointer>  correctedImages;
    if ( downSamplingFactor > 1 )
      {
      //
        for(unsigned int nima = 0; nima < numberOfImages; nima++)
          {
          std::cout << "Parameters estimated on downsampled image; applying bias field in full resolution" << std::endl;
          correctedImages.push_back(EMSegmenter::ConvertToNativeImageType( segmenter->BiasCorrect( nonDownsampledImages[nima], downSamplingFactor ) ));
          }
      }
    else
      {
       for(unsigned int nima = 0; nima < numberOfImages; nima++)
          {
            correctedImages.push_back(EMSegmenter::ConvertToNativeImageType( segmenter->GetBiasCorrectedImage(nima) ));
           }
      }

    // Make sure spacing, origin, and directions are copied
    for(unsigned int nima = 0; nima < numberOfImages; nima++)
      {
        correctedImages[nima]->SetSpacing( nonDownsampledImages[nima]->GetSpacing() );
        correctedImages[nima]->SetOrigin( nonDownsampledImages[nima]->GetOrigin() );
        correctedImages[nima]->SetDirection( nonDownsampledImages[nima]->GetDirection() );
      }
    // Store the corrected image in persistent memory
    std::vector<int> correctedImageHandles(numberOfImages,0);
    for(unsigned int nima = 0; nima < numberOfImages; nima++)
    {
      correctedImageHandles[nima] = kvl::MatlabObjectArray::GetInstance()->AddObject( correctedImages[nima] );
    }
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = numberOfImages;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    double* z;
    z = mxGetPr(plhs[0]);    
    for(unsigned int nima = 0; nima < numberOfImages; nima++, z++)
    {
      *( z ) = correctedImageHandles[nima];
    }
    
    }
  
protected:
  BiasFieldCorrectImage() {};
  virtual ~BiasFieldCorrectImage() {};


  BiasFieldCorrectImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl







