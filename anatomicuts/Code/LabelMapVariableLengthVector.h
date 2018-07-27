#ifndef __LabelMapVariableLengthVector_h
#define __LabelMapVariableLengthVector_h

#include "itkVariableLengthVector.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_transpose.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;
using namespace itk;

template <typename TValueType , class TMesh>
class LabelMapVariableLengthVector : public VariableLengthVector<TValueType>
{
public:
 
  /** The element type stored at each location in the Array. */
  typedef TValueType                                    ValueType;
  typedef TValueType                                    ComponentType;
  typedef typename NumericTraits< ValueType >::RealType RealValueType;
  typedef LabelMapVariableLengthVector                          Self;
  typedef VariableLengthVector<TValueType>              Superclass;
  typedef TMesh						MeshType;
  typedef typename MeshType::Pointer			MeshPointerType;
  typedef typename MeshType::CellPixelType		CellType; //??
  typedef typename MeshType::CellAutoPointer		CellAutoPointerType; //??

  LabelMapVariableLengthVector():Superclass(){
	;
}; 

  /** Constructor with size. Size can only be changed by assignment */
  LabelMapVariableLengthVector(unsigned int dimension):Superclass(dimension){
	;
};
  LabelMapVariableLengthVector( ValueType* data, unsigned int sz, 
                                        bool LetArrayManageMemory = false):Superclass(data, sz, LetArrayManageMemory){};
 	void SetCell(MeshPointerType mesh, int cellID);
	const CellType* GetLabels() const 
	{
		return &this->m_labels;
	}
	const CellType* GetDirections() const 
	{
		return &this->m_directions;
	}
void Print() const{};
private:
	CellType m_labels;
	CellType m_directions;
};


#endif
