#ifndef __LabelPerPointVariableLengthVector_h
#define __LabelPerPointVariableLengthVector_h


#include "itkVariableLengthVector.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"
#include <unordered_map>

template <typename TValueType , class TMesh>
class LabelPerPointVariableLengthVector : public VariableLengthVector<TValueType>
{
public:
 
  /** The element type stored at each location in the Array. */
  typedef TValueType                                    ValueType;
  typedef TValueType                                    ComponentType;
  typedef typename NumericTraits< ValueType >::RealType RealValueType;
  typedef LabelPerPointVariableLengthVector                          Self;
  typedef VariableLengthVector<TValueType>              Superclass;
  typedef TMesh						MeshType;
  typedef typename MeshType::Pointer			MeshPointerType;
  typedef typename MeshType::ConstPointer		ConstMeshPointerType;
  typedef typename MeshType::PixelType			PixelType; //??
  typedef typename MeshType::CellPixelType		CellType; //??
  typedef typename MeshType::CellAutoPointer		CellAutoPointerType; //??
  typedef std::unordered_map<int,float> LabelsMapType;
  typedef typename std::vector<LabelsMapType> LabelsDirectionType;
  LabelPerPointVariableLengthVector():Superclass(){
	;
}; 

  /** Constructor with size. Size can only be changed by assignment */
	LabelPerPointVariableLengthVector(unsigned int dimension):Superclass(dimension){;};
	LabelPerPointVariableLengthVector( ValueType* data, unsigned int sz, bool LetArrayManageMemory = false):Superclass(data, sz, LetArrayManageMemory){};

 	void SetCell(MeshPointerType mesh, int cellID);

	const std::vector<CellType>* GetLabels() const 
	{
		return &this->m_labels;
	}
	int GetCellId() const
	{
		return this->m_cellId;
	}
	int GetLength() const
	{
		return this->m_length;
	}
	const LabelsDirectionType GetLabelsPerDirection() const
	{
		return this->m_labelsPerDirection;
	}
	int GetNumberOfPoints()const{ return this->m_numberOfPoints;}
void Print() const;
private:
	std::vector<CellType> m_labels;
	LabelsDirectionType m_labelsPerDirection;
	int m_cellId;
	int m_numberOfPoints;
	float  m_length;
};
#include "LabelPerPointVariableLengthVector.txx"
#endif
