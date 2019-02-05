#ifndef _MeshToImageFilter_h_
#define _MeshToImageFilter_h_

#include "itkImageSource.h"
#include "itkGaussianSpatialFunction.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShiftScaleImageFilter.h"

#include "itkContinuousIndex.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"
#include "itkEuclideanDistanceMetric.h"

using namespace itk;
 
  template<class TInputMesh, class TOutputImage>
    class ITK_EXPORT MeshToImageFilter :
  public ImageSource<TOutputImage>
  {
    
  public:
    
    typedef MeshToImageFilter     Self;
    typedef ImageSource<TOutputImage>  Superclass;
    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;

    itkNewMacro  (Self);
    itkTypeMacro (MeshToImageFilter, ImageSource);
    
    itkStaticConstMacro(ImageDimension, unsigned int,
			TOutputImage::ImageDimension);

    typedef TInputMesh                              MeshType;
    typedef typename MeshType::Pointer              MeshTypePointer;


    /** Image typedefs */
    typedef TOutputImage                            OutputImageType;
    typedef typename OutputImageType::PixelType     OutputPixelType;
    typedef typename OutputImageType::RegionType    OutputImageRegionType;
    typedef typename OutputImageType::PointType     PointType;
    typedef typename OutputImageType::IndexType     IndexType;
    typedef typename OutputImageType::SpacingType   SpacingType;
    typedef typename OutputImageType::SizeType      SizeType;
    typedef typename OutputImageType::DirectionType DirectionType;

    /** Typedefs for base image. */
    typedef ImageBase< itkGetStaticConstMacro( ImageDimension ) > ImageBaseType;
    

    /** Set the size of the output image. */
    virtual void SetOutputSize( const SizeType & size );

    /** Get the size of the output image. */
    virtual const SizeType & GetOutputSize();

    /** Set the start index of the output largest possible region. 
     * The default is an index of all zeros. */
    virtual void SetOutputIndex( const IndexType & index );
    
    /** Get the start index of the output largest possible region. */
    virtual const IndexType & GetOutputIndex();
    
    /** Set the region of the output image. */
    itkSetMacro( OutputRegion, OutputImageRegionType );

    /** Get the region of the output image. */
    itkGetConstReferenceMacro( OutputRegion, OutputImageRegionType );
     
    /** Set the output image spacing. */
    itkSetMacro( OutputSpacing, SpacingType );
    virtual void SetOutputSpacing( const double* values );

    /** Get the output image spacing. */
    itkGetConstReferenceMacro( OutputSpacing, SpacingType );

    /** Set the output image origin. */
    itkSetMacro( OutputOrigin, PointType );
    virtual void SetOutputOrigin( const double* values);

    /** Get the output image origin. */
    itkGetConstReferenceMacro( OutputOrigin, PointType );

    /** Set the output direction cosine matrix. */
    itkSetMacro( OutputDirection, DirectionType );
    itkGetConstReferenceMacro( OutputDirection, DirectionType );

    /** Helper method to set the output parameters based on this image */
    void SetOutputParametersFromImage( const ImageBaseType * image );
  
    /** TransformToVelocityFieldSource produces a vector image. */
    virtual void GenerateOutputInformation( void );

    /** Set/Get the vector of positions */
    void SetInput(MeshTypePointer mesh)
    { m_Input = mesh; }
    MeshTypePointer GetInput()
    { return this->m_Input; }
    
   float BinaryImageOfLabels(int label, int flip);
  protected:
    MeshToImageFilter();
    ~MeshToImageFilter(){};
    
    /** Threaded implementation */
    void GenerateData();
    
    void PrintSelf(std::ostream& os, Indent indent) const
    {
      Superclass::PrintSelf(os,indent);
    }

  private:
    MeshToImageFilter (const Self&);
    void operator=(const Self&);

    OutputImageRegionType   m_OutputRegion;      // region of the output image
    SpacingType             m_OutputSpacing;     // output image spacing
    PointType               m_OutputOrigin;      // output image origin
    DirectionType           m_OutputDirection;   // output image direction cosines

    MeshTypePointer         m_Input;
	bool m_usingLabels;
  };


#include "MeshToImageFilter.txx"
#endif
