#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "kvlAtlasMeshMultiAlphaDrawer.h"



namespace kvl
{

class RasterizeAtlasMesh : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef RasterizeAtlasMesh         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( RasterizeAtlasMesh, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // alphaBuffer = kvlRasterizeAtlasMesh( mesh, DIM, classNumber )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 2 ) || !mxIsInt64( prhs[ 0 ] ) || 
         ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Some typedefs
    typedef AtlasMeshAlphaDrawer::ImageType  AlphaImageType;
    typedef AlphaImageType::SizeType  SizeType;
    typedef AtlasMeshMultiAlphaDrawer::ImageType  MultiAlphasImageType;
   
    // Retrieve input arguments
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMesh ) )
    if ( strcmp(typeid( *object ).name(), typeid( kvl::AtlasMesh ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }
    kvl::AtlasMesh::ConstPointer mesh = static_cast< const kvl::AtlasMesh* >( object.GetPointer() );
    double*  tmp = mxGetPr( prhs[ 1 ] );
    SizeType  imageSize;
    for ( int i = 0; i < 3; i++, tmp++ )
      {
      imageSize[ i ] = static_cast< int >( *tmp );
      }
      
    int classNumber = -1;
    if ( nrhs > 2 )
      {
      double* tmp = mxGetPr( prhs[ 2 ] ); 
      classNumber = static_cast< int >( *tmp );
      }
    
    //std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    //std::cout << "imageSize: " << imageSize << std::endl;
    //std::cout << "classNumber: " << classNumber << std::endl;
    

    if ( classNumber >= 0 )
      {
      // Rasterize the specified prior. If the class number is 0, then pre-fill everything
      // so that parts not overlayed by the mesh are still considered to the background
      kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
      alphaDrawer->SetRegions( imageSize );
      alphaDrawer->SetClassNumber( classNumber );
      if ( classNumber == 0 )
        {
        ( const_cast< AlphaImageType* >( alphaDrawer->GetImage() ) )->FillBuffer( 1.0 );
        }
      alphaDrawer->Rasterize( mesh );

      
      // Finally, copy the buffer contents into a Matlab matrix
      mwSize  dims[ 3 ];
      for ( int i = 0; i < 3; i++ )
        {
        dims[ i ] = imageSize[ i ];
        }
      //plhs[ 0 ] = mxCreateNumericArray( 3, dims, mxSINGLE_CLASS, mxREAL );
      //float*  data = static_cast< float* >( mxGetData( plhs[ 0 ] ) ); 
      plhs[ 0 ] = mxCreateNumericArray( 3, dims, mxUINT16_CLASS, mxREAL );
      unsigned short*  data = static_cast< unsigned short* >( mxGetData( plhs[ 0 ] ) ); 

      itk::ImageRegionConstIterator< AlphaImageType >  
                        it( alphaDrawer->GetImage(),
                            alphaDrawer->GetImage()->GetBufferedRegion() );
      for ( ;!it.IsAtEnd(); ++it, ++data )
        {
        float  probability = it.Value();
        if ( probability < 0 )
          { 
          probability = 0.0f;  
          }
        else if ( probability > 1 )
          {
          probability = 1.0f;
          }  
        *data = static_cast< unsigned short >( probability * 65535 + .5 );
        }
      }
    else
      {
      // Rasterize all priors simultaneously
      const unsigned int  numberOfClasses = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
      //std::cout << "numberOfClasses: " << numberOfClasses << std::endl;

      //std::cout << "Rasterizing mesh..." << std::flush;
      kvl::AtlasMeshMultiAlphaDrawer::Pointer  drawer = kvl::AtlasMeshMultiAlphaDrawer::New();
      drawer->SetRegions( imageSize );
      //std::cout << "here: " << numberOfClasses << std::endl;
      drawer->Rasterize( mesh );
      MultiAlphasImageType::ConstPointer  alphasImage = drawer->GetImage();
      //std::cout << "done" << std::endl;
      
      
      // Convert into 4-D Matlab matrix
      //std::cout << "Converting into Matlab format..." << std::flush;

      mwSize  dims[ 4 ];
      dims[ 0 ] = imageSize[ 0 ];
      dims[ 1 ] = imageSize[ 1 ];
      dims[ 2 ] = imageSize[ 2 ];
      dims[ 3 ] = numberOfClasses;
      // The following is much faster than 
      // 
      //   plhs[ 0 ] = mxCreateNumericArray( 4, dims, mxSINGLE_CLASS, mxREAL );
      //
      // because mxCreateNumericArray initializes all the matrix values to zero, which is something
      // we don't need
      int  numberOfElements = 1;
      for ( int i = 0; i < 4; i++ )
        {
        numberOfElements *= dims[ i ];  
        }
      //plhs[ 0 ] = mxCreateNumericMatrix( 0, 0, mxSINGLE_CLASS, mxREAL );
      plhs[ 0 ] = mxCreateNumericMatrix( 0, 0, mxUINT16_CLASS, mxREAL );
      mxSetDimensions( plhs[ 0 ], dims, 4 );
      //mxSetData( plhs[ 0 ], mxMalloc( sizeof(float) * numberOfElements ) );
      mxSetData( plhs[ 0 ], mxMalloc( sizeof(unsigned short) * numberOfElements ) );

      //float*  data = static_cast< float* >( mxGetData( plhs[ 0 ] ) ); 
      unsigned short*  data = static_cast< unsigned short* >( mxGetData( plhs[ 0 ] ) ); 
      for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        // Loop over all voxels
        itk::ImageRegionConstIterator< MultiAlphasImageType >  alphasIt( alphasImage,
                                                                         alphasImage->GetBufferedRegion() );
        for ( ;!alphasIt.IsAtEnd(); ++alphasIt, ++data )
          {
          float  probability = 0.0f;
          if ( alphasIt.Value().sum() == 0 )
            {
            // Outside of mesh so not computed. Let's put that to background
            if ( classNumber == 0 )
              {
              //probability = 1.0f;  
              }
            }
          else  
            {
            //*data = alphasIt.Value()[ classNumber ];
            probability = alphasIt.Value()[ classNumber ];
            if ( probability < 0 )
              { 
              probability = 0.0f;  
              }
            else if ( probability > 1 )
              {
              probability = 1.0f;
              }
            } // End test if outside mesh area
            
          *data = static_cast< unsigned short >( probability * 65535 + .5 );
          }
          
        //std::cout << "Done" << std::endl;  
        } // End loop over all labels
      //std::cout << "done" << std::endl;
      
      }  // End test if rasterizing one prior or all priors simultaneously


    }
  
protected:
  RasterizeAtlasMesh() {};
  virtual ~RasterizeAtlasMesh() {};


  RasterizeAtlasMesh(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



