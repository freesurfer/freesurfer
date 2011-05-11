#ifndef __kvlAtlasMeshAlphaDrawer_h
#define __kvlAtlasMeshAlphaDrawer_h

#include "kvlAtlasMeshRasterizor.h"


namespace kvl
{


namespace FragmentProcessor 
{

/**
 *
 */
class DrawAlpha
{
public:
  
  typedef itk::Image< float, 3 >  ImageType;

  DrawAlpha() 
    {
    m_Image = 0;
    m_AlphaInVertex0 = 1.0f;
    m_AlphaInVertex1 = 0.0f;
    m_AlphaInVertex2 = 0.0f;
    m_AlphaInVertex3 = 0.0f;
    m_LabelNumber = 0;
    m_Mesh = 0;
    }

  ~DrawAlpha() {};

  void AllocateImage( ImageType::SizeType  size )
    {
    m_Image = ImageType::New();
    m_Image->SetRegions( size );
    m_Image->Allocate();
    m_Image->FillBuffer( 0 );
    }
    
  void SetImage( ImageType* image )
    { m_Image = image; }

  const ImageType* GetImage() const
    { return m_Image; }

  void SetLabelNumber( unsigned char labelNumber )
    { m_LabelNumber = labelNumber; }
  
  inline void operator()( const float& pi0, const float& pi1, const float& pi2, const float& pi3 )
    {
    float  alpha = pi0 * m_AlphaInVertex0 +
                   pi1 * m_AlphaInVertex1 +
                   pi2 * m_AlphaInVertex2 +
                   pi3 * m_AlphaInVertex3;
#if 1                   
    m_Image->SetPixel( m_Index, alpha );
#else
    const unsigned long  offset = m_Image->ComputeOffset( m_Index );
    *( m_Image->GetBufferPointer()+offset ) = alpha;
#endif
    
    
    
      
    m_Index[ 0 ]++;
    }
    
  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
    {
    m_Index[ 0 ] = x;
    m_Index[ 1 ] = y;
    m_Index[ 2 ] = z;
    }
    
  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
    {
    // Cache the alpha of the specified class in each of the vertices of this tetrahedron
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );
          
    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    AtlasAlphasType  alphas0;
    AtlasAlphasType  alphas1;
    AtlasAlphasType  alphas2;
    AtlasAlphasType  alphas3;
    
    alphas0 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    alphas1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    alphas2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    alphas3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    
    m_AlphaInVertex0 = alphas0[ m_LabelNumber ];
    m_AlphaInVertex1 = alphas1[ m_LabelNumber ];
    m_AlphaInVertex2 = alphas2[ m_LabelNumber ];
    m_AlphaInVertex3 = alphas3[ m_LabelNumber ];

    return true;
    }
    
  inline void SetMesh( const AtlasMesh* mesh )
    {
    m_Mesh = mesh;
    }
    
private:

  ImageType::Pointer  m_Image;
  ImageType::IndexType  m_Index;
  
  float  m_AlphaInVertex0;
  float  m_AlphaInVertex1;
  float  m_AlphaInVertex2;
  float  m_AlphaInVertex3;
  
  unsigned char  m_LabelNumber;
  
  AtlasMesh::ConstPointer  m_Mesh;
  
    
};


} // End namespace FragmentProcessor


/**
 *
 */
class AtlasMeshAlphaDrawer: public AtlasMeshRasterizor< FragmentProcessor::DrawAlpha >
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshAlphaDrawer  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::DrawAlpha >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshAlphaDrawer, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::ImageType  AlphaImageType;

  /** */
  void SetLabelNumber( unsigned char labelNumber )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin(); 
          it != this->GetFragmentProcessors().end(); ++it )
      {  
      it->SetLabelNumber( labelNumber );
      }
    }
    
  /** */
  virtual void SetLabelImage( const LabelImageType*  labelImage )
    {
    // Allocate the result image for the first fragment processor, and point to the other
    // fragment processors to that 
    std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();  
    it->AllocateImage( labelImage->GetLargestPossibleRegion().GetSize() );
    AlphaImageType::Pointer  image = const_cast< AlphaImageType* >( it->GetImage() );
    for ( ; it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetImage( image );
      }
    // Invoke superclass' implementation
    Superclass::SetLabelImage( labelImage );
    }
  
  /** */
  const AlphaImageType*  GetAlphaImage() const
    { return this->GetFragmentProcessor().GetImage(); }
    
  /**  */
  void Rasterize( const AtlasMesh* mesh )
    {
    // Rasterize using multithreading
    Superclass::Rasterize( mesh, true );
    }
  
protected:
  AtlasMeshAlphaDrawer() {};
  virtual ~AtlasMeshAlphaDrawer() {};

private:
  AtlasMeshAlphaDrawer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
};


} // end namespace kvl

#endif

