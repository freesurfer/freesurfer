#ifndef __kvlAtlasMeshMinLogLikelihoodCalculator_h
#define __kvlAtlasMeshMinLogLikelihoodCalculator_h

#include "kvlAtlasMeshRasterizor.h"


namespace kvl
{


namespace FragmentProcessor 
{

class CalculateMinLogLikelihood
{
public:
  
  CalculateMinLogLikelihood() 
    {
    m_SourcePointer = 0;
    m_Mesh = 0;
    m_MinLogLikelihood = 0.0f;
    m_mapCompToComp = 0;
    }

  ~CalculateMinLogLikelihood() {};

  inline void operator()( const float& pi0, const float& pi1, const float& pi2, const float& pi3 )
    {
    
    float  denominator = 0;
    if(m_mapCompToComp == 0)  // no collapsed labels
      {

      // Calculate the numerator of the Labelification weights W0, W1, and W2 in this point
      float  W0 = m_AlphasInVertex0[ *m_SourcePointer ] * pi0;
      float  W1 = m_AlphasInVertex1[ *m_SourcePointer ] * pi1;
      float  W2 = m_AlphasInVertex2[ *m_SourcePointer ] * pi2;
      float  W3 = m_AlphasInVertex3[ *m_SourcePointer ] * pi3;
    
      // Calculate the sum
      denominator = W0 + W1 + W2 + W3 + 1e-15;
    
      }
    else   // version with collapsed labels
      { 
      float denominator = 1e-15;
      for(int ind=0; ind<m_mapCompToComp[*m_SourcePointer].size(); ind++)
        {
        const unsigned char k = m_mapCompToComp[*m_SourcePointer][ind];

// // std::cout << " ammllc " << ((int) *m_SourcePointer) << " " << ((int) k) << std::endl;

        float  W0 = m_AlphasInVertex0[ k ] * pi0;
        float  W1 = m_AlphasInVertex1[ k ] * pi1;
        float  W2 = m_AlphasInVertex2[ k ] * pi2;
        float  W3 = m_AlphasInVertex3[ k ] * pi3;
    
        // Calculate the sum
        denominator += (W0 + W1 + W2 + W3);
        }
      }

    m_MinLogLikelihood -= log( denominator );

    // Move on to the next pixel      
    m_SourcePointer++;
    }
    
  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
    {
    m_SourcePointer = sourcePointer;
    }
    
  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
    {
    // Cache relevant elements of the vertices of this triangle
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );
          
    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    m_AlphasInVertex0 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    m_AlphasInVertex1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    m_AlphasInVertex2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    m_AlphasInVertex3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;

    return true;
    }
    
  inline void SetMesh( const AtlasMesh* mesh )
    {
    m_Mesh = mesh;
    
    m_MinLogLikelihood = 0.0f;
    }
    
    
  float GetMinLogLikelihood() const
    {
    return m_MinLogLikelihood;
    }

  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    { m_mapCompToComp = mapCompToComp; } 
  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }
    
private:

  const unsigned char*  m_SourcePointer;
  
  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;
  AtlasAlphasType  m_AlphasInVertex3;
  
  AtlasMesh::ConstPointer  m_Mesh;
  
  float  m_MinLogLikelihood;

   std::vector<unsigned char > *m_mapCompToComp;
};


 
 

} // End namespace FragmentProcessor





class AtlasMeshMinLogLikelihoodCalculator: public AtlasMeshRasterizor< FragmentProcessor::CalculateMinLogLikelihood >
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshMinLogLikelihoodCalculator  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::CalculateMinLogLikelihood >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshMinLogLikelihoodCalculator, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  
  /** */
  float GetMinLogLikelihood() const
    { return this->GetFragmentProcessor().GetMinLogLikelihood(); }

  /** */
  void  SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetMapCompToComp( mapCompToComp );  
      }
    }

  /** */
  std::vector<unsigned char > * GetMapCompToComp()
    {
    return this->GetFragmentProcessor().GetMapCompToComp();
    }
  
protected:
  AtlasMeshMinLogLikelihoodCalculator() {};
  virtual ~AtlasMeshMinLogLikelihoodCalculator() {};

private:
  AtlasMeshMinLogLikelihoodCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
};


} // end namespace kvl

#endif

