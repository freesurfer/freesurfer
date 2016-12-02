#ifndef __kvlAtlasMeshDeformationOptimizer_h
#define __kvlAtlasMeshDeformationOptimizer_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "kvlAtlasMesh.h"


namespace kvl
{

// Events generated
itkEventMacro( DeformationStartEvent, itk::UserEvent );
itkEventMacro( DeformationIterationEvent, itk::UserEvent );
itkEventMacro( DeformationEndEvent, itk::UserEvent );


/**
 *
 */
class AtlasMeshDeformationOptimizer: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationOptimizer  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  //itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationOptimizer, itk::Object );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;
  typedef itk::Image< AtlasAlphasType, 3 >  ProbabilityImageType;
  typedef itk::Image< unsigned char, 3 >  SegmentedImageType;
  typedef itk::AffineTransform< double, 3 >  TransformType;


  /** */
  /** */
  void  SetImages( std::vector<itk::Image< float, 3 >::Pointer> images )
    {
    m_Images = images;
    m_Initialized = false;
    }

  /** */
  const itk::Image< float, 3 >::Pointer GetImage(int ind) const
    { 
    return m_Images[ind]; 
    }
  /** */
  void  SetProbabilityImage( const ProbabilityImageType*  probabilityImage )
    {
    m_ProbabilityImage = probabilityImage;
    m_Initialized = false;
    }

  /** */
  const ProbabilityImageType*  GetProbabilityImage() const
    {
    return m_ProbabilityImage;
    }

  /** */
  void  SetSegmentedImage( const SegmentedImageType*  segmentedImage )
    {
    m_SegmentedImage = segmentedImage;
    m_Initialized = false;
    }

  /** */
  const SegmentedImageType*  GetSegmentedImage() const
    {
    return m_SegmentedImage;
    }


  /** */
  void  SetMesh( AtlasMesh* mesh )
    {
    m_Mesh = mesh;
    m_Initialized = false;
    }

  /** */
  const AtlasMesh*  GetMesh() const
    {
    return m_Mesh;
    }


  /** */
  void SetMeans( std::vector< vnl_vector<float> >& means )
    { 
    m_Means = means;
    m_Initialized = false; 
    }

  /** */
  void SetPrecisions( std::vector< vnl_matrix<float> >& precisions )
    {
      m_Precisions = precisions;
      m_Initialized = false; 
    }

  /** */
  void  SetMeshToImageTransform( const TransformType* meshToImageTransform )
    {
    m_MeshToImageTransform = meshToImageTransform;
    m_Initialized = false;  
    }

  const TransformType* GetMeshToImageTransform() const
    {
    return m_MeshToImageTransform;
    }

  // 
  unsigned int  GetIterationNumber() const
    { return m_IterationNumber; }
    
  //
  unsigned int  GetMaximumNumberOfIterations() const
    { return m_MaximumNumberOfIterations; }

  void  SetMaximumNumberOfIterations( unsigned int maximumNumberOfIterations )
    { m_MaximumNumberOfIterations = maximumNumberOfIterations; }

  //
  unsigned int  GetIterationEventResolution() const
    { return m_IterationEventResolution; }
    
  //
  void SetIterationEventResolution( unsigned int  iterationEventResolution )
    { m_IterationEventResolution = iterationEventResolution; 
      // std::cout << "Now setting m_IterationEventResolution: " << m_IterationEventResolution << std::endl;
    }
 
  /** */
  virtual double GetMinLogLikelihoodTimesPrior() const = 0;
  
  /** */
  virtual bool Go() = 0;

  /** */
  void SetNumberOfThreads( int numberOfThreads )
    { m_NumberOfThreads = numberOfThreads; }
    
  /** */
  int  GetNumberOfThreads() const
    { return m_NumberOfThreads; }

  //
  void  SetVerbose( bool verbose )
    { m_Verbose = verbose; }
    
  //
  bool  GetVerbose() const
    { return m_Verbose; }

  // Set/Get mapping of collapsed labels.
  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    { m_mapCompToComp = mapCompToComp; }
  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }

protected:
  AtlasMeshDeformationOptimizer();
  virtual ~AtlasMeshDeformationOptimizer();
  
  virtual void Initialize()
    {
    m_Initialized = true;  
    }

  //
  int  m_IterationNumber;
  int  m_MaximumNumberOfIterations;
  int  m_IterationEventResolution;


  std::vector<itk::Image< float, 3 >::Pointer> m_Images;
  ProbabilityImageType::ConstPointer  m_ProbabilityImage;
  SegmentedImageType::ConstPointer  m_SegmentedImage;
  AtlasMesh::Pointer  m_Mesh;
  TransformType::ConstPointer  m_MeshToImageTransform;

  std::vector< vnl_vector<float> > m_Means;
  std::vector< vnl_matrix<float> > m_Precisions;
  
  std::vector<unsigned char > *m_mapCompToComp;

  bool  m_Initialized;
  int  m_NumberOfThreads;
  bool  m_Verbose;

private:
  AtlasMeshDeformationOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
};


} // end namespace kvl

#endif

