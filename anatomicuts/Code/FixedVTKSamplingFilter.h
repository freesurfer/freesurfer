#ifndef _FixedVTKSamplingFilter_h_
#define _FixedVTKSamplingFilter_h_


#include "itkMeshToMeshFilter.h"
using namespace itk;
  template <class TInputMesh, class TOutputMesh>
    class FixedVTKSamplingFilter :
  public MeshToMeshFilter<TInputMesh, TOutputMesh>
  {
  public:
    typedef FixedVTKSamplingFilter Self;
    typedef MeshToMeshFilter<TInputMesh, TOutputMesh>       Superclass;
    typedef SmartPointer<Self>                              Pointer;
    typedef SmartPointer<const Self>                        ConstPointer;

    itkNewMacro  (Self);
    itkTypeMacro (FixedVTKSamplingFilter, MeshToMeshFilter);

    typedef TOutputMesh                         OutputMeshType;
    typedef typename OutputMeshType::MeshTraits OutputMeshTraits;
    typedef typename OutputMeshType::PointType  OutputPointType;
    typedef typename OutputMeshType::PixelType  OutputPixelType;

    
    typedef TInputMesh                         InputMeshType;
    typedef typename InputMeshType::MeshTraits InputMeshTraits;
    typedef typename InputMeshType::PointType  InputPointType;
    typedef typename InputMeshType::PixelType  InputPixelType;

    
    /** Some convenient typedefs. */
    typedef typename OutputMeshType::Pointer         OutputMeshPointer;
    typedef typename OutputMeshType::CellTraits      OutputCellTraits;
    typedef typename OutputMeshType::CellIdentifier  OutputCellIdentifier;
    typedef typename OutputMeshType::CellType        OutputCellType;
    typedef typename OutputMeshType::CellAutoPointer OutputCellAutoPointer;
    typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
    typedef typename OutputCellTraits::PointIdIterator     OutputPointIdIterator;
    
    typedef typename OutputMeshType::PointsContainerPointer
      OutputPointsContainerPointer;
    
    typedef typename OutputMeshType::PointsContainer
      OutputPointsContainer;

    typedef typename OutputMeshType::CellsContainer
      OutputCellsContainer;

    typedef typename OutputMeshType::CellsContainerPointer
      OutputCellsContainerPointer;

    typedef PolylineCell<OutputCellType>                      PolylineCellType;

    typedef typename InputMeshType::Pointer         InputMeshPointer;
    typedef typename InputMeshType::CellTraits      InputMeshCellTraits;
    typedef typename InputMeshType::CellIdentifier  InputMeshCellIdentifier;
    typedef typename InputMeshType::CellType        InputMeshCellType;
    typedef typename InputMeshType::CellAutoPointer InputMeshCellAutoPointer;
    typedef typename InputMeshType::PointIdentifier InputMeshPointIdentifier;
    typedef typename InputMeshCellTraits::PointIdIterator   InputMeshPointIdIterator;
    
    typedef typename InputMeshType::PointsContainerPointer
      InputPointsContainerPointer;
    
    typedef typename InputMeshType::PointsContainer
      InputPointsContainer;

    typedef typename InputMeshType::CellsContainer
      InputCellsContainer;

    typedef typename InputMeshType::CellsContainerPointer
      InputCellsContainerPointer;

    int GetSampling()
    {
      return this->sampling;
    }

    void SetSampling(int sampling)
    {
       this->sampling = sampling;
    }
/*betweeen 0 and 1*/
	void SetExtendPercentage(float percentage)
	{
		this->m_extendPercentage = percentage;
	}
    void SampleToFiberLenght();
  protected:
    FixedVTKSamplingFilter();
    ~FixedVTKSamplingFilter() {};

    virtual void GenerateData (void);
    
  private:
    int sampling;
    FixedVTKSamplingFilter (const Self&);
    void operator=(const Self&);
	float m_extendPercentage;
  };
  


#endif
