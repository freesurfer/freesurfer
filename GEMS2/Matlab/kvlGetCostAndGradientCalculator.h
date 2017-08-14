#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "kvlConditionalGaussianEntropyCostAndGradientCalculator.h"
#include "kvlMutualInformationCostAndGradientCalculator.h"


namespace kvl
{

class GetCostAndGradientCalculator : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetCostAndGradientCalculator         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetCostAndGradientCalculator, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // calculator = kvlGetCostAndGradientCalculator( typeName, image(s), boundaryCondition, transform )
  
    // Make sure input arguments are correct
    const  std::string  usageString = "Usage: calculator = kvlGetCostAndGradientCalculator( typeName, image(s), boundaryCondition, [ transform ], [ means ], [ precisions ] )\n where typeName = {'AtlasMeshToIntensityImage','ConditionalGaussianEntropy','MutualInformation'}\n and boundaryCondition = {'Sliding', 'Affine', 'Translation', 'None'}";
    if ( ( nrhs < 3 ) || 
         !mxIsChar( prhs[ 0 ] ) || 
         !mxIsInt64( prhs[ 1 ] ) || 
         !mxIsChar( prhs[ 2 ] ) ) 
      {
      mexErrMsgTxt( usageString.c_str() ); 
      }

      
    // Retrieve input image(s)
    typedef AtlasMeshToIntensityImageCostAndGradientCalculator::ImageType  ImageType;

    const int  N = mxGetN( prhs[ 1 ] );
    const int  M = mxGetM( prhs[ 1 ] );
    int numberOfContrasts = 0;
    if ( N < M )
      {
      numberOfContrasts = M;
      }
    else
      {
      numberOfContrasts = N;
      }

    mexPrintf( "numberOfContrasts = %d\n", numberOfContrasts );
    std::vector< ImageType::ConstPointer >  images;
    uint64_T *  imagesHandle =  static_cast< uint64_T * >( mxGetData( prhs[ 1 ] ) );
    for ( int imageNummber = 0; imageNummber < numberOfContrasts; imageNummber++, imagesHandle++)
       {
       const int handle = *(imagesHandle);
       std::cout << "Image: " << handle << std::endl;
       itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( handle );
       if ( typeid( *(object) ) != typeid( ImageType ) )
         {
         mexErrMsgTxt( "image doesn't refer to the correct ITK object type" );
         }
       ImageType::ConstPointer constImage = static_cast< const ImageType* >( object.GetPointer() );
       //ImageType::Pointer image = const_cast< ImageType* >( constImage.GetPointer() );
       images.push_back( constImage ); 
       }
    
    
    // Retrieve transform if one is provided
    typedef CroppedImageReader::TransformType  TransformType;
    TransformType::ConstPointer  constTransform = 0;
    if ( nrhs > 3 )
      {
      // Sanity check
      if ( !mxIsInt64( prhs[ 3 ] ) )
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      const int transformHandle = *( static_cast< int* >( mxGetData( prhs[ 3 ] ) ) );
      itk::Object::ConstPointer  object = kvl::MatlabObjectArray::GetInstance()->GetObject( transformHandle );
      if ( typeid( *object ) != typeid( TransformType ) )
        {
        mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
        }
      constTransform = static_cast< const TransformType* >( object.GetPointer() );
      }
      
    
    // Retrieve means if they are provided
    std::vector< vnl_vector< float > >  means;
    if ( nrhs > 4 )
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 4 ] ) ) 
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      //
      const int  numberOfClasses = mxGetDimensions( prhs[ 4 ] )[ 0 ];
      const int  numberOfContrasts  = mxGetDimensions( prhs[ 4 ] )[ 1 ];
      //mexPrintf("numberOfClasses = %d\n",numberOfClasses);
      //mexPrintf("numberOfContrasts = %d\n",numberOfContrasts);
      for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        vnl_vector< float >  mean( numberOfContrasts, 0.0f );
        for ( int imageNumber = 0; imageNumber < numberOfContrasts; imageNumber++ )
          {
          mean[ imageNumber ] = (mxGetPr( prhs[ 4 ] ))[ classNumber + numberOfClasses*imageNumber ];
          }
        means.push_back( mean );
        }
        
      // Print what we recovered
      for (unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        vnl_vector<float> miini = means[classNumber];
        for ( unsigned int nima = 0; nima < numberOfContrasts; nima++ )
          {
          //mexPrintf("means[%d][%d] = %f\n",nima, classNumber, miini[nima]);
          }
        }

      } // End test if means are provided 
        
        
        
    // Retrieve precisions if they are provided
    std::vector< vnl_matrix< float > >  precisions;
    if ( nrhs > 5 )
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 5 ] ) ) 
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      //
      const int  numberOfClasses = mxGetDimensions( prhs[ 4 ] )[ 0 ];
      const int  numberOfContrasts  = mxGetDimensions( prhs[ 4 ] )[ 1 ];
      //mexPrintf("numberOfClasses = %d\n",numberOfClasses);
      //mexPrintf("numberOfContrasts = %d\n",numberOfContrasts);
      
      // Does not really matter which way you read these in, because the precisions are symmetric matrices
      //transpose wont do any harm. 
      for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        vnl_matrix< float >  precision( numberOfContrasts, numberOfContrasts, 0.0f );
        for ( unsigned int row = 0; row < numberOfContrasts; row++ )
          {
          for ( unsigned int col = 0; col < numberOfContrasts; col++ )
            {
            precision[ row ][ col ] = mxGetPr( prhs[ 5 ] )[ classNumber + row * numberOfClasses + col * numberOfClasses * numberOfContrasts ];
            }
          }
        precisions.push_back( precision );
        }

      // Print what we've recovered
      for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        vnl_matrix<float> precMat = precisions[classNumber];
        for ( unsigned int row = 0; row < numberOfContrasts; row++ )
          {
          for ( unsigned int col = 0; col < numberOfContrasts; col++ )
            {
            //mexPrintf("precisions[%d][%d][%d] = %f\n",row,col,classNumber,precMat[row][col]);
            }
          }
        }
      
      } // End test if precisions are provided 
        
        
        
    // Construct the correct type of calculator
    AtlasMeshPositionCostAndGradientCalculator::Pointer  calculator = 0; 
    const std::string  typeName = mxArrayToString( prhs[ 0 ] );
    switch( typeName[ 0 ] ) 
      {
      case 'A': 
        {
        std::cout << "AtlasMeshToIntensityImage" << std::endl;
        AtlasMeshToIntensityImageCostAndGradientCalculator::Pointer  myCalculator 
                          = AtlasMeshToIntensityImageCostAndGradientCalculator::New();
        myCalculator->SetImages( images );
        myCalculator->SetParameters( means, precisions );
        calculator = myCalculator;
        break;
        } 
      case 'C': 
        {
        std::cout << "ConditionalGaussianEntropy" << std::endl;
        ConditionalGaussianEntropyCostAndGradientCalculator::Pointer  myCalculator 
                          = ConditionalGaussianEntropyCostAndGradientCalculator::New();
        myCalculator->SetImage( images[ 0 ] );
        calculator = myCalculator;
        break;
        } 
      case 'M': 
        {
        std::cout << "MutualInformation" << std::endl;
        MutualInformationCostAndGradientCalculator::Pointer  myCalculator 
                          = MutualInformationCostAndGradientCalculator::New();
        myCalculator->SetImage( images[ 0 ] );
        calculator = myCalculator;
        break;
        } 
      default:
        {
        mexErrMsgTxt( "typeName not understood" );
        break;
        }
      }

    
    // Specify the correct type of boundary condition  
    const std::string  boundaryCondition = mxArrayToString( prhs[ 2 ] );
    switch( boundaryCondition[ 0 ] ) 
      {
      case 'S': 
        {
        std::cout << "SLIDING" << std::endl;
        calculator->SetBoundaryCondition( AtlasMeshPositionCostAndGradientCalculator::SLIDING );
        if ( constTransform.GetPointer() )
          {
          calculator->SetMeshToImageTransform( constTransform );
          }
        break;
        } 
      case 'A': 
        {
        std::cout << "AFFINE" << std::endl;
        calculator->SetBoundaryCondition( AtlasMeshPositionCostAndGradientCalculator::AFFINE );
        break;
        } 
      case 'T': 
        {
        std::cout << "TRANSLATION" << std::endl;
        calculator->SetBoundaryCondition( AtlasMeshPositionCostAndGradientCalculator::TRANSLATION );
        break;
        } 
      case 'N': 
        {
        std::cout << "NONE" << std::endl;
        calculator->SetBoundaryCondition( AtlasMeshPositionCostAndGradientCalculator::NONE );
        break;
        } 
      default:
        {
        mexErrMsgTxt( "boundaryCondition not understood" );
        break;
        }
      }

                             
    // Store the calculator in persistent memory
    const int calculatorHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( calculator );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = calculatorHandle;
    
    }
  
protected:
  GetCostAndGradientCalculator() {};
  virtual ~GetCostAndGradientCalculator() {};


  GetCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


