#include "pyKvlImage.h"
#include "itkCastImageFilter.h"
#include "itkMGHImageIOFactory.h"
#include "pyKvlNumpy.h"
#include "pyKvlTransform.h"
#include <pybind11/numpy.h>
#include <itkImageFileWriter.h>


KvlImage::KvlImage(const std::string &imageFileName) {
    // Add support for MGH file format to ITK. An alternative way to add this by default would be
    // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
    itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

    // Read the image
    kvl::CroppedImageReader::Pointer reader = kvl::CroppedImageReader::New();
    reader->Read( imageFileName.c_str() );

    // Convert the image to float
    typedef itk::CastImageFilter< kvl::CroppedImageReader::ImageType, ImageType > CasterType;
    CasterType::Pointer  caster = CasterType::New();
    caster->SetInput( reader->GetImage() );
    caster->Update();

    // Store the image and transform in persistent memory
    m_image = caster->GetOutput();
    m_transform = TransformType::New();
    reader->GetWorldToImageTransform()->GetInverse( m_transform );
    std::cout << "Read image: " << imageFileName << std::endl;
}

KvlImage::KvlImage(const py::array_t<float> &buffer) {
    typedef typename ImageType::PixelType  PixelType;

    // Determine the size of the image to be created
    typedef typename ImageType::SizeType  SizeType;
    SizeType  imageSize;
    for ( int i = 0; i < 3; i++ )
    {
        imageSize[ i ] = buffer.shape( i );
        //std::cout << "imageSize[ i ]: " << imageSize[ i ] << std::endl;
    }

    // Construct an ITK image
    m_image = ImageType::New();
    m_image->SetRegions( imageSize );
    m_image->Allocate();
    m_transform = TransformType::New();
    // Loop over all voxels and copy contents
    itk::ImageRegionIterator< ImageType >  it( m_image,
                                               m_image->GetBufferedRegion() );
    const float *bufferPointer = buffer.data(0);
    for ( ;!it.IsAtEnd(); ++it, ++bufferPointer )
    {
        it.Value() = *bufferPointer;
    }
}

std::unique_ptr<KvlTransform> KvlImage::GetTransform() {
    return std::unique_ptr<KvlTransform>(new KvlTransform(m_transform));
}

py::array_t<float> KvlImage::GetImageBuffer() {
    auto region = m_image->GetLargestPossibleRegion();
    auto shape = region.GetSize();
    auto* const buffer = new float[region.GetNumberOfPixels()];

    itk::ImageRegionConstIterator< ImageType > it( m_image, region );
    float* buffer_p = buffer;
    for ( ;!it.IsAtEnd(); ++it, ++buffer_p ){
        *buffer_p = it.Value();
    }
    return createNumpyArrayFStyle({shape[0], shape[1], shape[2]}, buffer);
}

void KvlImage::Write(std::string fileName, KvlTransform &transform) {

    // If transform is given, retrieve and apply it
    ImageType *image = m_image;
    if ( true )
    {
        typedef kvl::CroppedImageReader::TransformType  TransformType;
        // In order not to modify the original image, we create a new one. The proper way of doing this
        // would be to only copy the header information and of course not the pixel intensities, but I'm
        // too lazy now to figure out how to do it in ITK
        typedef itk::CastImageFilter< ImageType, ImageType >  CasterType;
        CasterType::Pointer  caster = CasterType::New();
        caster->SetInput( m_image );
        caster->Update();
        ImageType *image = caster->GetOutput();


        // Apply the transform
        ImageType::PointType   newOrigin;
        ImageType::SpacingType  newSpacing;
        ImageType::DirectionType  newDirection;
        for ( int i = 0; i < 3; i++ )
        {
            // Offset part
            newOrigin[ i ] = transform.m_transform->GetOffset()[ i ];

            // For every column, determine norm (which will be voxel spacing), and normalize direction
            double  normOfColumn = 0.0;
            for ( int j = 0; j < 3; j++ )
            {
                normOfColumn += pow( transform.m_transform->GetMatrix()[ j ][ i ], 2 );
            }
            normOfColumn = sqrt( normOfColumn );
            newSpacing[ i ] = normOfColumn;
            for ( int j = 0; j < 3; j++ )
            {
                newDirection[ j ][ i ] = transform.m_transform->GetMatrix()[ j ][ i ] / normOfColumn;
            }
        }
        image->SetOrigin( newOrigin );
        image->SetSpacing( newSpacing );
        image->SetDirection( newDirection );

    } // End test if transform is given


    // Write it out
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( fileName.c_str() );
    writer->SetInput( image );
    writer->Update();
    std::cout << "Wrote image to file " << fileName << std::endl;

}
