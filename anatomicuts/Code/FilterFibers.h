#ifndef _FilterFibers_h_
#define _FilterFibers_h_

#include "itkMeshToMeshFilter.h"
using namespace itk;

  template <class TInputMesh, class TOutputMesh>
    class FilterFibers :
  public MeshToMeshFilter<TInputMesh, TOutputMesh>
  {
  public:
    typedef FilterFibers Self;
    typedef MeshToMeshFilter<TInputMesh, TOutputMesh>       Superclass;
    typedef SmartPointer<Self>                              Pointer;
    typedef SmartPointer<const Self>                        ConstPointer;

    itkNewMacro  (Self);
    itkTypeMacro (FilterFibers, MeshToMeshFilter);

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

    float GetPercentage()
	{
      return this->m_percentage;
    }

    void SetPercentage(float percentage)
    {
       
	if(percentage <0.0 || percentage > 1.0 )
	{	
		std::cout << "Percentage should be between 0.0 and 1.0: dividing by 100. " << std::endl;
		percentage /=100;
	}
	this->m_percentage = percentage;
    }

  protected:
    FilterFibers(){};
    ~FilterFibers() {};

    virtual void GenerateData (void);

    
  private:
    FilterFibers (const Self&);
    void operator=(const Self&);
    float m_percentage;
  };
  


#endif
