#ifndef __LabelVariableLengthVector_h
#define __LabelVariableLengthVector_h

#include "itkVariableLengthVector.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_transpose.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>

using namespace itk;

template <typename TValueType , class TMesh>
class LabelVariableLengthVector : public VariableLengthVector<TValueType>
{
public:
 
  /** The element type stored at each location in the Array. */
  typedef TValueType                                    ValueType;
  typedef TValueType                                    ComponentType;
  typedef typename NumericTraits< ValueType >::RealType RealValueType;
  typedef LabelVariableLengthVector                          Self;
  typedef VariableLengthVector<TValueType>              Superclass;
  typedef TMesh						MeshType;
  typedef typename MeshType::Pointer			MeshPointerType;
  typedef typename MeshType::CellType		CellType; //??
  typedef typename MeshType::CellAutoPointer		CellAutoPointerType; //??

  LabelVariableLengthVector():Superclass(){
	;
}; 

  /** Constructor with size. Size can only be changed by assignment */
  LabelVariableLengthVector(unsigned int dimension):Superclass(dimension){
	;
};
  LabelVariableLengthVector( ValueType* data, unsigned int sz, 
                                        bool LetArrayManageMemory = false):Superclass(data, sz, LetArrayManageMemory){};
 	void SetCell(MeshPointerType mesh, int cellID);
void Print() const{};
private:
};




#endif
