#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "kvlConditionalGaussianEntropyCostAndGradientCalculator.h"
#include "kvlMutualInformationCostAndGradientCalculator.h"
#include "kvlAtlasMeshToPointSetCostAndGradientCalculator.h"


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
    const  std::string  usageString = "Usage: calculator = kvlGetCostAndGradientCalculator( typeName, image(s), boundaryCondition, [ transform ], [ means ], [ variances ], [ mixtureWeights ], [ numberOfGaussiansPerClass ], [ targetPoints ] )\n where typeName = {'AtlasMeshToIntensityImage','ConditionalGaussianEntropy','MutualInformation','PointSet'}\n and boundaryCondition = {'Sliding', 'Affine', 'Translation', 'None'}";
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
    std::vector< vnl_vector< double > >  means;
    if ( nrhs > 4 )
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 4 ] ) ) 
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      //
      const int  numberOfGaussians = mxGetDimensions( prhs[ 4 ] )[ 0 ];
      const int  numberOfContrasts  = mxGetDimensions( prhs[ 4 ] )[ 1 ];
      //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
      //mexPrintf("numberOfContrasts = %d\n",numberOfContrasts);
      for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
        {
        vnl_vector< double >  mean( numberOfContrasts, 0.0f );
        for ( int contrastNumber = 0; contrastNumber < numberOfContrasts; contrastNumber++ )
          {
          mean[ contrastNumber ] = (mxGetPr( prhs[ 4 ] ))[ gaussianNumber + numberOfGaussians*contrastNumber ];
          }
        means.push_back( mean );
        }
        
      if ( false )  
        {
        // Print what we recovered
        for (unsigned int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
          {
          vnl_vector<double> miini = means[gaussianNumber];
          for ( unsigned int nima = 0; nima < numberOfContrasts; nima++ )
            {
            mexPrintf("means[%d][%d] = %f\n",nima, gaussianNumber, miini[nima]);
            }
          }
        } // End test if printing  

      } // End test if means are provided 
        
        
        
    // Retrieve variances if they are provided
    std::vector< vnl_matrix< double > >  variances;
    if ( nrhs > 5 )
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 5 ] ) ) 
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      //
      const int  numberOfGaussians = mxGetDimensions( prhs[ 4 ] )[ 0 ];
      const int  numberOfContrasts  = mxGetDimensions( prhs[ 4 ] )[ 1 ];
      //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
      //mexPrintf("numberOfContrasts = %d\n",numberOfContrasts);
      
      // Does not really matter which way you read these in, because the variances are symmetric matrices
      // transpose won't do any harm. 
      for ( unsigned int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
        {
        vnl_matrix< double >  variance( numberOfContrasts, numberOfContrasts, 0.0f );
        for ( unsigned int row = 0; row < numberOfContrasts; row++ )
          {
          for ( unsigned int col = 0; col < numberOfContrasts; col++ )
            {
            variance[ row ][ col ] = mxGetPr( prhs[ 5 ] )[ gaussianNumber + row * numberOfGaussians + col * numberOfGaussians * numberOfContrasts ];
            }
          }
        variances.push_back( variance );
        }

      if ( false )
        {
        // Print what we've recovered
        for ( unsigned int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
          {
          vnl_matrix< double >  varMat = variances[gaussianNumber];
          for ( unsigned int row = 0; row < numberOfContrasts; row++ )
            {
            for ( unsigned int col = 0; col < numberOfContrasts; col++ )
              {
              mexPrintf("variances[%d][%d][%d] = %f\n",row,col,gaussianNumber,varMat[row][col]);
              }
            }
          }
        } // End test if printing out
        
      
      } // End test if variances are provided 
        
        
        
    // Retrieve mixtureWeights if they are provided
    std::vector< double >  mixtureWeights;
    if ( nrhs > 6 )
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 6 ] ) ) 
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      //
      const int  numberOfGaussians = mxGetDimensions( prhs[ 4 ] )[ 0 ];
      mixtureWeights = std::vector< double >( numberOfGaussians, 0.0f );
      for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
        {
        mixtureWeights[ gaussianNumber ] = (mxGetPr( prhs[ 6 ] ))[ gaussianNumber ];
        }
    
      if ( false )
        {
        // Print what we've recovered
        for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
          {
          mexPrintf("mixtureWeights[%d] = %f\n", gaussianNumber, mixtureWeights[ gaussianNumber ] );
          }  
        } // End test printing
        
      } // End test if mixtureWeights are provided

      
    // Retrieve numberOfGaussiansPerClass if they are provided
    std::vector< int >  numberOfGaussiansPerClass;
    if ( nrhs > 7 )
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 7 ] ) ) 
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      //
      const int  numberOfClasses = mxGetNumberOfElements( prhs[ 7 ] );
      numberOfGaussiansPerClass = std::vector< int >( numberOfClasses, 0 );
      for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        numberOfGaussiansPerClass[ classNumber ] = static_cast< int >( (mxGetPr( prhs[ 7 ] ))[ classNumber ] );
        }
    
      if ( false )
        {
        // Print what we've recovered
        for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
          {
          mexPrintf("numberOfGaussiansPerClass[%d] = %d\n", classNumber, numberOfGaussiansPerClass[ classNumber ] );
          }  
        } // End test printing
        
      } // End test if numberOfGaussiansPerClass are provided
        
        
    // Retrieve targetPoints if they are provided
    AtlasMesh::PointsContainer::Pointer  targetPoints = AtlasMesh::PointsContainer::New();
    if ( nrhs > 8 )
      {
      // Sanity check
      if ( ( !mxIsDouble( prhs[ 8 ] ) ) || ( mxGetDimensions( prhs[ 8 ] )[ 1 ] != 3 ) )
        {
        mexErrMsgTxt( usageString.c_str() ); 
        }
        
      //
      const int  numberOfPoints = mxGetDimensions( prhs[ 8 ] )[ 0 ];
      mexPrintf("numberOfPoints = %d\n",numberOfPoints );
      for ( int pointNumber = 0; pointNumber < numberOfPoints; pointNumber++ )
        {
        AtlasMesh::PointType  point;
        for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )
          {
          point[ dimensionNumber ] =  (mxGetPr( prhs[ 8 ] ))[ pointNumber + numberOfPoints*dimensionNumber ]; 
          }
          
        targetPoints->InsertElement( targetPoints->Size(), point );
        } // End loop over all points  

      if ( true )  
        {
        // Print what we recovered
        for ( AtlasMesh::PointsContainer::ConstIterator  it = targetPoints->Begin(); 
              it != targetPoints->End(); ++it )
          {
          std::cout << it.Value() << std::endl; 
          }
          
        } // End test if printing  

      } // End test if targetPoints are provided 
        
        
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
        myCalculator->SetParameters( means, variances, mixtureWeights, numberOfGaussiansPerClass );
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
      case 'P': 
        {
        std::cout << "PointSet" << std::endl;
        AtlasMeshToPointSetCostAndGradientCalculator::Pointer  myCalculator 
                          = AtlasMeshToPointSetCostAndGradientCalculator::New();
        myCalculator->SetTargetPoints( targetPoints );
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


