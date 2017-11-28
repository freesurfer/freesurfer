#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIO.h"


namespace kvl
{

class WriteImage : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef WriteImage         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( WriteImage, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    // std::cout << "I am " << this->GetNameOfClass() 
    //           << " and I'm running! " << std::endl;
              
              
      // kvlWriteImage( image, fileName, transform )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 2 ) || !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsChar( prhs[ 1 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input arguments
    const int imageHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    const std::string  fileName = mxArrayToString( prhs[1] );

    // Retrieve the image
    itk::Object::Pointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( imageHandle );
    typedef itk::Image< float, 3 >   FloatImageType;
    // if ( typeid( *object ) != typeid( FloatImageType ) )
    if ( strcmp(typeid( *object ).name(), typeid( FloatImageType ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "image doesn't refer to the correct ITK object type" );
      }
    FloatImageType::Pointer  image 
            = static_cast< FloatImageType* >( object.GetPointer() );
    
    
    // If transform is given, retrieve and apply it
    if ( nrhs > 2 )
      {
      // Retrieve the transform
      if ( !mxIsInt64( prhs[ 2 ] ) )
        {
        mexErrMsgTxt( "Incorrect arguments" );
        }
      typedef CroppedImageReader::TransformType  TransformType;
      const int transformHandle = *( static_cast< int* >( mxGetData( prhs[ 2 ] ) ) );
      object = kvl::MatlabObjectArray::GetInstance()->GetObject( transformHandle );
      // if ( typeid( *object ) != typeid( TransformType ) )
      if ( strcmp(typeid( *object ).name(), typeid( TransformType ).name()) )  // Eugenio: MAC compatibility
        {
        mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
        }
      TransformType::ConstPointer  transform = static_cast< const TransformType* >( object.GetPointer() );
    
#if 0      
      // If this is MGH format, fiddle with the directions (rotation around Z-axis)
      itk::ImageIOBase::Pointer  io = itk::ImageIOFactory::CreateImageIO(  fileName.c_str(), 
                                                                           itk::ImageIOFactory::WriteMode );
      if ( io )
        {
        if ( dynamic_cast< itk::MGHImageIO* >( io.GetPointer() ) )
          {
          std::cout << "==========================================" << std::endl;  
          std::cout << "Dealing with MGH format here - rotating orientation around Z-axis!" << std::endl;
          std::cout << "==========================================" << std::endl;
          TransformType::OutputVectorType  scaling;
          scaling[ 0 ] = -1.0;
          scaling[ 1 ] = -1.0;
          scaling[ 2 ] = 1.0;
          const_cast< TransformType* >( transform.GetPointer() )->Scale( scaling );
          }  
        }
#endif
        
      
      // In order not to modify the original image, we create a new one. The proper way of doing this
      // would be to only copy the header information and of course not the pixel intensities, but I'm
      // too lazy now to figure out how to do it in ITK
      typedef itk::CastImageFilter< FloatImageType, FloatImageType >  CasterType;
      CasterType::Pointer  caster = CasterType::New();
      caster->SetInput( image );
      caster->Update();
      image = caster->GetOutput();
      
      
      // Apply the transform
      FloatImageType::PointType   newOrigin;
      FloatImageType::SpacingType  newSpacing;
      FloatImageType::DirectionType  newDirection;
      for ( int i = 0; i < 3; i++ )
        {
        // Offset part
        newOrigin[ i ] = transform->GetOffset()[ i ];

        // For every column, determine norm (which will be voxel spacing), and normalize direction
        double  normOfColumn = 0.0;
        for ( int j = 0; j < 3; j++ )
          {
          normOfColumn += pow( transform->GetMatrix()[ j ][ i ], 2 );
          }
        normOfColumn = sqrt( normOfColumn );
        newSpacing[ i ] = normOfColumn;
        for ( int j = 0; j < 3; j++ )
          {
          newDirection[ j ][ i ] = transform->GetMatrix()[ j ][ i ] / normOfColumn;
          }
        }
      image->SetOrigin( newOrigin );
      image->SetSpacing( newSpacing );
      image->SetDirection( newDirection );

      } // End test if transform is given


    // Write it out
    typedef itk::ImageFileWriter< FloatImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( fileName.c_str() );
    writer->SetInput( image );
    writer->Update();
    std::cout << "Wrote image to file " << fileName << std::endl;
    
    }
  
protected:
  WriteImage() {};
  virtual ~WriteImage() {};


  WriteImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


