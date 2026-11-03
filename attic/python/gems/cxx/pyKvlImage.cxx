#include "pyKvlImage.h"
#include "itkCastImageFilter.h"
#include "itkMGHImageIOFactory.h"
#include "pyKvlNumpy.h"
#include "pyKvlTransform.h"
#include <pybind11/numpy.h>
#include <itkImageFileWriter.h>
#include "itkDiscreteGaussianImageFilter.h"

static bool mgh_factory_is_loaded;
void load_mgh_factory() {
    if (!mgh_factory_is_loaded) {
        mgh_factory_is_loaded = true;
        // Add support for MGH file format to ITK. An alternative way to add this by default would be
        // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
        itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );
    }
}

py::array_t<float> KvlImage::image_to_numpy(ImagePointer image) {
    auto region = image->GetLargestPossibleRegion();
    auto shape = region.GetSize();
    auto* const buffer = new float[region.GetNumberOfPixels()];

    itk::ImageRegionConstIterator< ImageType > it( image, region );
    float* buffer_p = buffer;
    for ( ;!it.IsAtEnd(); ++it, ++buffer_p ){
        *buffer_p = it.Value();
    }
    return createNumpyArrayFStyle({shape[0], shape[1], shape[2]}, buffer);
}

ImagePointer KvlImage::numpy_to_image(const py::array_t<float> &buffer) {
    typedef typename ImageType::PixelType  PixelType;

    // Determine the size of the image to be created
    typedef typename ImageType::SizeType  SizeType;
    SizeType  imageSize;
    for ( int i = 0; i < 3; i++ ) {
        imageSize[ i ] = buffer.shape( i );
    }

    // Construct an ITK image
    ImagePointer image = ImageType::New();
    image->SetRegions( imageSize );
    image->Allocate();

    // Loop over all voxels and copy contents
    itk::ImageRegionIterator< ImageType >  it( image, image->GetBufferedRegion() );
    const float *bufferPointer = buffer.data(0);
    for ( ;!it.IsAtEnd(); ++it, ++bufferPointer ) {
        it.Value() = *bufferPointer;
    }
    return image;
}

KvlImage::KvlImage(const std::string &imageFileName) {
    load_mgh_factory();

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

KvlImage::KvlImage(const std::string &imageFileName, const std::string &boundingFileName):
        m_non_cropped_image_size(3), m_crop_slices(3) {
    load_mgh_factory();

    // Read the image
    kvl::CroppedImageReader::Pointer reader = kvl::CroppedImageReader::New();
    reader->SetExtraFraction(0.0);
    reader->Read(imageFileName.c_str(), boundingFileName.c_str());

    // Convert the image to float
    typedef itk::Image<float, 3> ImageType;
    typedef itk::CastImageFilter<kvl::CroppedImageReader::ImageType, ImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(reader->GetImage());
    caster->Update();

    // Store the image and transform in persistent memory
    m_image = caster->GetOutput();
    m_transform = reader->GetTransform()->Clone();

    for ( int index = 0; index < 3; index++ ) {
        int start = reader->GetOriginalImageRegion().GetIndex(index) - reader->GetCroppedImageRegion().GetIndex(index);
        int stop = start + reader->GetCroppedImageRegion().GetSize(index);
        m_crop_slices[index] = py::slice(start, stop, 1);
        m_non_cropped_image_size[index] = reader->GetOriginalImageOriginalRegion().GetSize(index);
    }
    std::cout << "Read image: " << imageFileName << " cropped by " << boundingFileName << std::endl;
}

KvlImage::KvlImage(const py::array_t<float> &buffer) {
    m_image = numpy_to_image(buffer);
    m_transform = TransformType::New();
}

std::vector<int> KvlImage::GetNonCroppedImageSize() {
    return m_non_cropped_image_size;
}

py::tuple KvlImage::GetCropSlices() {
    return m_crop_slices;
}

std::unique_ptr<KvlTransform> KvlImage::GetTransform() {
    return std::unique_ptr<KvlTransform>(new KvlTransform(m_transform));
}

py::array_t<float> KvlImage::GetImageBuffer() {
    return image_to_numpy(m_image);
}

void KvlImage::Write(std::string fileName, KvlTransform &transform) {

    // If transform is given, retrieve and apply it
    ImagePointer image = m_image;
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
        image = caster->GetOutput();


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
    std::cout << "Wrote image and transform to file " << fileName << std::endl;

}

void KvlImage::WriteImage(std::string fileName) {

    // If transform is given, retrieve and apply it
    ImagePointer image = m_image;

    // Write it out
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( fileName.c_str() );
    writer->SetInput( image );
    writer->Update();
    std::cout << "Wrote image to file " << fileName << std::endl;

}

py::array_t<float> KvlImage::smoothImageBuffer(const py::array_t<float>& imageBuffer, std::vector<double> sigmas)
{
    ImagePointer image = numpy_to_image(imageBuffer);

    // Do the actual work in ITK
    typedef itk::Image< float, 3 >  InternalImageType;
    typedef itk::CastImageFilter< ImageType, InternalImageType >   CasterType;
    typedef itk::DiscreteGaussianImageFilter< InternalImageType, InternalImageType >  SmootherType;
    typedef itk::CastImageFilter< InternalImageType, ImageType >  BackCasterType;

    typename CasterType::Pointer caster = CasterType::New();
    caster->SetInput( image );
    SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetInput( caster->GetOutput() );
    smoother->SetMaximumError( 0.1 );
    smoother->SetUseImageSpacingOff();
    double  variances[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
        variances[ i ] = sigmas[ i ] * sigmas[ i ];
    }
    smoother->SetVariance( variances );
    typename BackCasterType::Pointer  backCaster = BackCasterType::New();
    backCaster->SetInput( smoother->GetOutput() );
    backCaster->Update();
    auto smoothedImage = backCaster->GetOutput();

    return image_to_numpy(smoothedImage);
}
