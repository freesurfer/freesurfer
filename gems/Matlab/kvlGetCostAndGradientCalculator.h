#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageLogDomainCostAndGradientCalculator.h"
#include "kvlAtlasMeshToWishartGaussMixtureCostAndGradientCalculator.h"
#include "kvlAtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator.h"
#include "kvlAtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator.h"
#include "kvlConditionalGaussianEntropyCostAndGradientCalculator.h"
#include "kvlMutualInformationCostAndGradientCalculator.h"
#include "kvlAtlasMeshToPointSetCostAndGradientCalculator.h"
#include "kvlCroppedImageReader.h"


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
    std::string  usageString = "Usage: calculator = kvlGetCostAndGradientCalculator( typeName, image(s),";
    usageString.append(" boundaryCondition, [ transform ], [ means ], [ variances ], [ mixtureWeights ],");
    usageString.append(" [ numberOfGaussiansPerClass ], [ targetPoints ] )\n where typeName =");
    usageString.append(" {'AtlasMeshToIntensityImage','ConditionalGaussianEntropy','MutualInformation'");
    usageString.append(",'PointSet','WishartMixtureToAtlasMesh','FrobeniusMixutureToAtlasMesh'}\n and boundaryCondition = {'Sliding', 'Affine', 'Translation', 'None'}");
    usageString.append("\n\n When the type is WishartMixtureToAtlasMesh the targetPoints argument is skipped");
    usageString.append(" and the 9th-12th parameters are :\n\ncalculator = kvlGetCostAndGradientCalculator(...,");
    usageString.append(" DTIimages , degreesOfFreedom , scaleMatrices , wmmMixtureWeights , numberOfWishartsPerClass)");


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
    int numberOfGMMContrasts = 0;
    if ( N < M )
      {
      numberOfGMMContrasts = M;
      }
    else
      {
      numberOfGMMContrasts = N;
      }

    mexPrintf( "numberOfGMMContrasts = %d\n", numberOfGMMContrasts );
    std::vector< ImageType::ConstPointer >  images;
    uint64_T *  imagesHandle =  static_cast< uint64_T * >( mxGetData( prhs[ 1 ] ) );
    for ( int imageNummber = 0; imageNummber < numberOfGMMContrasts; imageNummber++, imagesHandle++)
       {
       const int handle = *(imagesHandle);
       std::cout << "Image: " << handle << std::endl;
       itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( handle );
       // if ( typeid( *(object) ) != typeid( ImageType ) )
       if ( strcmp(typeid( *object ).name(), typeid( ImageType ).name()) )  // Eugenio: MAC compatibility
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
      // if ( typeid( *object ) != typeid( TransformType ) )
      if ( strcmp(typeid( *object ).name(), typeid( TransformType ).name()) )  // Eugenio: MAC compatibility
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
        
        
    // Retrieve targetPoints if they are provided or DTimages if the typeName is WishartMixtureToAtlasMesh
    const std::string  typeName = mxArrayToString( prhs[ 0 ] ); // check calculator type first
    std::vector< ImageType::ConstPointer >  DTIimages;
    AtlasMesh::PointsContainer::Pointer  targetPoints = AtlasMesh::PointsContainer::New();
    if ( ( nrhs > 8 ) && ( typeName[0] != 'W' ) && ( typeName[0] != 'F' ) && ( typeName[0] != 'D' ))
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
    else if (( nrhs > 8 ) && (typeName[0] == 'W'))
    {
        const int  N = mxGetN( prhs[ 8 ] );
        const int  M = mxGetM( prhs[ 8 ] );
        int numberOfWMMContrasts = 0;
        if ( N < M )
        {
            numberOfWMMContrasts = M;
        }
        else
        {
            numberOfWMMContrasts = N;
        }

        if (numberOfWMMContrasts != 7)
        {
            mexErrMsgTxt( "Seven DTI images required for WMM. One log determinant, 3 diagonal entries and 3 off-diagonal." );
        }

        mexPrintf( "numberOfWMMContrasts = %d\n", numberOfWMMContrasts );
        uint64_T *  imagesHandleWMM =  static_cast< uint64_T * >( mxGetData( prhs[ 8 ] ) );
        for ( int imageNummberWMM = 0; imageNummberWMM < numberOfWMMContrasts; imageNummberWMM++, imagesHandleWMM++)
        {
            const int handle = *(imagesHandleWMM);
            std::cout << "Image: " << handle << std::endl;
            itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( handle );
            // if ( typeid( *(object) ) != typeid( ImageType ) )
            if ( strcmp(typeid( *object ).name(), typeid( ImageType ).name()) )  // Eugenio: MAC compatibility
            {
                mexErrMsgTxt( "image doesn't refer to the correct ITK object type" );
            }
            ImageType::ConstPointer constImage = static_cast< const ImageType* >( object.GetPointer() );
            //ImageType::Pointer image = const_cast< ImageType* >( constImage.GetPointer() );
            DTIimages.push_back( constImage );
        }
    }
    else if (( nrhs > 8 ) && (typeName[0] == 'F'))
    {
        const int  N = mxGetN( prhs[ 8 ] );
        const int  M = mxGetM( prhs[ 8 ] );
        int numberOfFMMContrasts = 0;
        if ( N < M )
        {
            numberOfFMMContrasts = M;
        }
        else
        {
            numberOfFMMContrasts = N;
        }

        if (numberOfFMMContrasts != 6)
        {
            mexErrMsgTxt( "Six log-DTI images required for WMM. 3 diagonal entries and 3 off-diagonal in order 11,12,13,22,23,33." );
        }

        mexPrintf( "numberOfWMMContrasts = %d\n", numberOfFMMContrasts );
        uint64_T *  imagesHandleWMM =  static_cast< uint64_T * >( mxGetData( prhs[ 8 ] ) );
        for ( int imageNummberWMM = 0; imageNummberWMM < numberOfFMMContrasts; imageNummberWMM++, imagesHandleWMM++)
        {
            const int handle = *(imagesHandleWMM);
            std::cout << "Image: " << handle << std::endl;
            itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( handle );
            // if ( typeid( *(object) ) != typeid( ImageType ) )
            if ( strcmp(typeid( *object ).name(), typeid( ImageType ).name()) )  // Eugenio: MAC compatibility
            {
                mexErrMsgTxt( "image doesn't refer to the correct ITK object type" );
            }
            ImageType::ConstPointer constImage = static_cast< const ImageType* >( object.GetPointer() );
            //ImageType::Pointer image = const_cast< ImageType* >( constImage.GetPointer() );
            DTIimages.push_back( constImage );
        }
    }
    else if (( nrhs > 8 ) && (typeName[0] == 'D'))
    {
        const int  N = mxGetN( prhs[ 8 ] );
        const int  M = mxGetM( prhs[ 8 ] );
        int numberOfFMMContrasts = 0;
        if ( N < M )
        {
            numberOfFMMContrasts = M;
        }
        else
        {
            numberOfFMMContrasts = N;
        }

        if (numberOfFMMContrasts != 4)
        {
            mexErrMsgTxt( "Four images required for MM. 1 FA and 3 elements of the primary eigenvector." );
        }

        mexPrintf( "numberOfWMMContrasts = %d\n", numberOfFMMContrasts );
        uint64_T *  imagesHandleWMM =  static_cast< uint64_T * >( mxGetData( prhs[ 8 ] ) );
        for ( int imageNummberWMM = 0; imageNummberWMM < numberOfFMMContrasts; imageNummberWMM++, imagesHandleWMM++)
        {
            const int handle = *(imagesHandleWMM);
            std::cout << "Image: " << handle << std::endl;
            itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( handle );
            // if ( typeid( *(object) ) != typeid( ImageType ) )
            if ( strcmp(typeid( *object ).name(), typeid( ImageType ).name()) )  // Eugenio: MAC compatibility
            {
                mexErrMsgTxt( "image doesn't refer to the correct ITK object type" );
            }
            ImageType::ConstPointer constImage = static_cast< const ImageType* >( object.GetPointer() );
            //ImageType::Pointer image = const_cast< ImageType* >( constImage.GetPointer() );
            DTIimages.push_back( constImage );
        }
    }// end test if DTIimages are provided


    // Retrieve degreesOfFreedom if they are provided
    std::vector< double >  degreesOfFreedom;
    std::vector< double >  frobVariance;
    std::vector< double >  DSWbetaConcentration;
    if ( ( nrhs > 9 ) && (typeName[0] == 'W'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 9 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfWisharts = mxGetDimensions( prhs[ 9 ] )[ 0 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        for ( int wishartNumber = 0; wishartNumber < numberOfWisharts; wishartNumber++ )
        {
            double  dof = (mxGetPr( prhs[ 9 ] ))[ wishartNumber ];

            degreesOfFreedom.push_back( dof );
        }

        if ( false )
        {
            // Print what we recovered
            for (unsigned int wishartNumber = 0; wishartNumber < numberOfWisharts; wishartNumber++ )
            {
                double miini = degreesOfFreedom[wishartNumber];
                mexPrintf("degreesOfFreedom[%d] = %f\n", wishartNumber, miini);
            }
        } // End test if printing

    } // End test if degreesOfFreedom are provided
    else if ( ( nrhs > 9 ) && (typeName[0] == 'F'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 9 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfFrobenius = mxGetDimensions( prhs[ 9 ] )[ 0 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        for ( int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
        {
            double  fvar = (mxGetPr( prhs[ 9 ] ))[ wishartNumber ];

            frobVariance.push_back( fvar );
        }

        if ( false )
        {
            // Print what we recovered
            for (unsigned int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
            {
                double miini = frobVariance[wishartNumber];
                mexPrintf("frobVariance[%d] = %f\n", wishartNumber, miini);
            }
        } // End test if printing

    } // End test if degreesOfFreedom are provided
    else if ( ( nrhs > 9 ) && (typeName[0] == 'D'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 9 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfFrobenius = mxGetDimensions( prhs[ 9 ] )[ 0 ];
        mexPrintf("numberOfDSWbetas = %d\n",numberOfFrobenius);
        for ( int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
        {
            double  fvar = (mxGetPr( prhs[ 9 ] ))[ wishartNumber ];

            DSWbetaConcentration.push_back( fvar );
        }

        if ( false )
        {
            // Print what we recovered
            for (unsigned int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
            {
                double miini = DSWbetaConcentration[wishartNumber];
                mexPrintf("DSWbetaConcentration[%d] = %f\n", wishartNumber, miini);
            }
        } // End test if printing

    } // End test if DSWbetaConcentration are provided


    // Retrieve scaleMatrices if they are provided
    std::vector< vnl_matrix< double > >  scaleMatrices;
    std::vector< vnl_vector< double > >  frobMeans;
    std::vector< vnl_vector< double > >  DSWbetaMeans;
    if ( ( nrhs > 10 ) && (typeName[0] == 'W'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 10 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfWisharts = mxGetDimensions( prhs[ 10 ] )[ 0 ];
        const int  numberOfDimensions  = mxGetDimensions( prhs[ 10 ] )[ 1 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        //mexPrintf("numberOfContrasts = %d\n",numberOfContrasts);

        // Does not really matter which way you read these in, because the scales are symmetric matrices
        // transpose won't do any harm.
        for ( unsigned int wishartNumber = 0; wishartNumber < numberOfWisharts; wishartNumber++ )
        {
            vnl_matrix< double >  scale( numberOfDimensions, numberOfDimensions, 0.0f );
            for ( unsigned int row = 0; row < numberOfDimensions; row++ )
            {
                for ( unsigned int col = 0; col < numberOfDimensions; col++ )
                {
                    scale[ row ][ col ] = mxGetPr( prhs[ 10 ] )[ wishartNumber + row * numberOfWisharts + col * numberOfWisharts * numberOfDimensions ];
                }
            }
            scaleMatrices.push_back( scale );
        }

        if ( false )
        {
            // Print what we've recovered
            for ( unsigned int wishartNumber = 0; wishartNumber < numberOfWisharts; wishartNumber++ )
            {
                vnl_matrix< double >  varMat = scaleMatrices[wishartNumber];
                for ( unsigned int row = 0; row < numberOfDimensions; row++ )
                {
                    for ( unsigned int col = 0; col < numberOfDimensions; col++ )
                    {
                        mexPrintf("scaleMatrices[%d][%d][%d] = %f\n",row,col,wishartNumber,varMat[row][col]);
                    }
                }
            }
        } // End test if printing out
    } // End test if scaleMatrices are provided
    else  if ( ( nrhs > 10 ) && (typeName[0] == 'F') )
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 10 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfFrobenius = mxGetDimensions( prhs[ 10 ] )[ 0 ];
        const int  numberOfDimensions  = mxGetDimensions( prhs[ 10 ] )[ 1 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        //mexPrintf("numberOfContrasts = %d\n",numberOfContrasts);
        for ( int frobNumber = 0; frobNumber < numberOfFrobenius; frobNumber++ )
          {
          vnl_vector< double >  mean( numberOfDimensions, 0.0f );
          for ( int dimnNumber = 0; dimnNumber < numberOfDimensions; dimnNumber++ )
            {
            mean[ dimnNumber ] = (mxGetPr( prhs[ 10 ] ))[ frobNumber + numberOfFrobenius*dimnNumber ];
            }
          frobMeans.push_back( mean );
          }

        if ( false )
          {
          // Print what we recovered
          for (unsigned int frobNumber = 0; frobNumber < numberOfFrobenius; frobNumber++ )
            {
            vnl_vector<double> miini = frobMeans[frobNumber];
            for ( unsigned int nima = 0; nima < numberOfDimensions; nima++ )
              {
              mexPrintf("means[%d][%d] = %f\n",nima, frobNumber, miini[nima]);
              }
            }
          } // End test if printing
    }
    else  if ( ( nrhs > 10 ) &&  (typeName[0] == 'D'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 10 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfFrobenius = mxGetDimensions( prhs[ 10 ] )[ 0 ];
        const int  numberOfDimensions  = mxGetDimensions( prhs[ 10 ] )[ 1 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        //mexPrintf("numberOfContrasts = %d\n",numberOfContrasts);
        for ( int frobNumber = 0; frobNumber < numberOfFrobenius; frobNumber++ )
          {
          vnl_vector< double >  mean( numberOfDimensions, 0.0f );
          for ( int dimnNumber = 0; dimnNumber < numberOfDimensions; dimnNumber++ )
            {
            mean[ dimnNumber ] = (mxGetPr( prhs[ 10 ] ))[ frobNumber + numberOfFrobenius*dimnNumber ];
            }
          DSWbetaMeans.push_back( mean );
          }

        if ( false )
          {
          // Print what we recovered
          for (unsigned int frobNumber = 0; frobNumber < numberOfFrobenius; frobNumber++ )
            {
            vnl_vector<double> miini = DSWbetaMeans[frobNumber];
            for ( unsigned int nima = 0; nima < numberOfDimensions; nima++ )
              {
              mexPrintf("DSWbetaMeans[%d][%d] = %f\n",nima, frobNumber, miini[nima]);
              }
            }
          } // End test if printing
    }




    // Retrieve wmmMixtureWeights if they are provided
    std::vector< double >  wmmMixtureWeights;
    if ( ( nrhs > 11 ) && ((typeName[0] == 'W')||(typeName[0] == 'F')||(typeName[0] == 'D')))
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 11 ] ) )
        {
        mexErrMsgTxt( usageString.c_str() );
        }

      //
      const int  numberOfWisharts = mxGetDimensions( prhs[ 11 ] )[ 0 ];
      wmmMixtureWeights = std::vector< double >( numberOfWisharts, 0.0f );
      for ( int wishartNumber = 0; wishartNumber < numberOfWisharts; wishartNumber++ )
        {
        wmmMixtureWeights[ wishartNumber ] = (mxGetPr( prhs[ 11 ] ))[ wishartNumber ];
        }

      if ( false )
        {
        // Print what we've recovered
        for ( int gaussianNumber = 0; gaussianNumber < numberOfWisharts; gaussianNumber++ )
          {
          mexPrintf("mixtureWeightsWMM[%d] = %f\n", gaussianNumber, wmmMixtureWeights[ gaussianNumber ] );
          }
        } // End test printing

      } // End test if mixtureWeights are provided


    // Retrieve numberOfWishartsPerClass if they are provided
    std::vector< int >  numberOfWishartsPerClass;
    if ( nrhs > 12 )
      {
      // Sanity check
      if ( !mxIsDouble( prhs[ 12 ] ) )
        {
        mexErrMsgTxt( usageString.c_str() );
        }

      //
      const int  numberOfWMMClasses = mxGetNumberOfElements( prhs[ 12 ] );
      numberOfWishartsPerClass = std::vector< int >( numberOfWMMClasses, 0 );
      for ( int classNumber = 0; classNumber < numberOfWMMClasses; classNumber++ )
        {
        numberOfWishartsPerClass[ classNumber ] = static_cast< int >( (mxGetPr( prhs[ 12 ] ))[ classNumber ] );
        }

      if ( false )
        {
        // Print what we've recovered
        for ( int classNumber = 0; classNumber < numberOfWMMClasses; classNumber++ )
          {
          mexPrintf("numberOfWishartsPerClass[%d] = %d\n", classNumber, numberOfWishartsPerClass[ classNumber ] );
          }
        } // End test printing

      } // End test if numberOfWishartsPerClass are provided


    std::vector< double >  DSWbetaAlpha;
    double voxratio = 1.0;
    if ( ( nrhs > 13 ) && (typeName[0] == 'D'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 13 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfFrobenius = mxGetDimensions( prhs[ 13 ] )[ 0 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        for ( int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
        {
            double  fvar = (mxGetPr( prhs[ 13 ] ))[ wishartNumber ];

            DSWbetaAlpha.push_back( fvar );
        }

        if ( false )
        {
            // Print what we recovered
            for (unsigned int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
            {
                double miini = DSWbetaAlpha[wishartNumber];
                mexPrintf("DSWbetaAlpha[%d] = %f\n", wishartNumber, miini);
            }
        } // End test if printing

    } // End test if DSWbetaAlpha are provided start test for voxratio
    else  if ( ( nrhs > 13 ) && ((typeName[0] == 'F')  || (typeName[0] == 'W')) )
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 13 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        // get voxratio
        voxratio = (mxGetPr( prhs[ 13 ] ))[ 0 ];

        if ( false )
          {
            // Print what we recovered
            mexPrintf("voxratio= %f\n",voxratio);
          } // End test if printing
    } // End test for voxratio


    std::vector< double >  DSWbetaBeta;
    if ( ( nrhs > 14 ) && (typeName[0] == 'D'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 14 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfFrobenius = mxGetDimensions( prhs[ 14 ] )[ 0 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        for ( int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
        {
            double  fvar = (mxGetPr( prhs[ 14 ] ))[ wishartNumber ];

            DSWbetaBeta.push_back( fvar );
        }

        if ( false )
        {
            // Print what we recovered
            for (unsigned int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
            {
                double miini = DSWbetaBeta[wishartNumber];
                mexPrintf("DSWbetaBeta[%d] = %f\n", wishartNumber, miini);
            }
        } // End test if printing

    } // End test if DSWbetaBeta are provided


    std::vector< double >  logKummerSamples;
    if ( ( nrhs > 15 ) && (typeName[0] == 'D'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 15 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }

        //
        const int  numberOfFrobenius = mxGetDimensions( prhs[ 15 ] )[ 0 ];
        //mexPrintf("numberOfGaussians = %d\n",numberOfGaussians);
        for ( int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
        {
            double  fvar = (mxGetPr( prhs[ 15 ] ))[ wishartNumber ];

            logKummerSamples.push_back( fvar );
        }

        if ( false )
        {
            // Print what we recovered
            for (unsigned int wishartNumber = 0; wishartNumber < numberOfFrobenius; wishartNumber++ )
            {
                double miini = logKummerSamples[wishartNumber];
                mexPrintf("logKummerSamples[%d] = %f\n", wishartNumber, miini);
            }
        } // End test if printing

    } // End test if logKummerSamples are provided


    double negLogKummerIncrement;
    if ( ( nrhs > 16 ) && (typeName[0] == 'D'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 16 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }
            negLogKummerIncrement = (mxGetPr( prhs[ 16 ] ))[ 0 ];

        if ( false )
        {
                mexPrintf("DSWbetaBeta[%d] = %f\n", 0, negLogKummerIncrement);
        } // End test if printing

    } // End test if negLogKummerIncrement are provided


    if ( ( nrhs > 17 ) && (typeName[0] == 'D'))
    {
        // Sanity check
        if ( !mxIsDouble( prhs[ 17 ] ) )
        {
            mexErrMsgTxt( usageString.c_str() );
        }
            voxratio = (mxGetPr( prhs[ 17 ] ))[ 0 ];

        if ( false )
        {
            mexPrintf("voxratio= %f\n",voxratio);
        } // End test if printing

    } // End test if voxratio is provided


    // Construct the correct type of calculator
    AtlasMeshPositionCostAndGradientCalculator::Pointer  calculator = 0;
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
      case 'L': 
        {
        std::cout << "AtlasMeshToIntensityImageLogDomain" << std::endl;
        AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator::Pointer  myCalculator
                          = AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator::New();
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
      case 'W':
        {
        std::cout << "WishartMixtureToAtlasMesh" << std::endl;
        AtlasMeshToWishartGaussMixtureCostAndGradientCalculator::Pointer  myCalculator
                          = AtlasMeshToWishartGaussMixtureCostAndGradientCalculator::New();
        //myCalculator->SetGaussianImages( images );
        myCalculator->SetImages( images );
        //myCalculator->SetWishartImages( DTIimages );
        myCalculator->SetDiffusionImages( DTIimages );
        myCalculator->SetParameters( means, variances, mixtureWeights, numberOfGaussiansPerClass );
        myCalculator->SetDiffusionParameters( means.size(), 
                                              wmmMixtureWeights, numberOfWishartsPerClass, voxratio, degreesOfFreedom,
                                              scaleMatrices);
        calculator = myCalculator;
        break;
        }
      case 'F':
        {
        std::cout << "FrobeniusMixtureToAtlasMesh" << std::endl;
        AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator::Pointer  myCalculator
                          = AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator::New();
        //myCalculator->SetGaussianImages( images );
        myCalculator->SetImages( images );
        //myCalculator->SetFrobeniusImages( DTIimages );
        myCalculator->SetDiffusionImages( DTIimages );
        myCalculator->SetParameters( means, variances, mixtureWeights, numberOfGaussiansPerClass );
        myCalculator->SetDiffusionParameters( means.size(),
                                              wmmMixtureWeights, numberOfWishartsPerClass, voxratio, frobVariance,
                                              frobMeans);

        std::cout<<"frobMixtureWeights size = "<<wmmMixtureWeights.size()<<std::endl;

        calculator = myCalculator;
        break;
        }
      case 'D':
        {
        std::cout << "DSWbetaMixtureToAtlasMesh" << std::endl;
        AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator::Pointer  myCalculator
                        = AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator::New();
        //myCalculator->SetGaussianImages( images );
        myCalculator->SetImages( images );
        //myCalculator->SetDSWbetaImages( DTIimages );
        myCalculator->SetDiffusionImages( DTIimages );
        myCalculator->SetParameters( means, variances, mixtureWeights, numberOfGaussiansPerClass );
        myCalculator->SetDiffusionParameters( means.size(),
                                              wmmMixtureWeights, numberOfWishartsPerClass, voxratio, DSWbetaAlpha,
                                              DSWbetaMeans, DSWbetaBeta, DSWbetaConcentration,
                                              logKummerSamples, negLogKummerIncrement);
                                   //frobVariance, frobMeans, wmmMixtureWeights, numberOfWishartsPerClass);

        std::cout<<"DSWbMixtureWeights size = "<<wmmMixtureWeights.size()<<std::endl;

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


