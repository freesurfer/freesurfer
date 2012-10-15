#ifndef __itkAffine3DTransform_h
#define __itkAffine3DTransform_h

#include <iostream>
#include "itkTransform.h"
#include "itkExceptionObject.h"
#include "itkMatrix.h"
#include "itkVersor.h"

namespace itk
{

/** \brief Affine3DTransform of a vector space (e.g. space coordinates)
 *
 *
 * \ingroup Transforms
 */
template < class TScalarType=double >    // Data type for scalars (float or double)
class ITK_EXPORT Affine3DTransform :
  public Transform< TScalarType, 3, 3> // Dimensions of input and output spaces
{
public:
  /** Standard class typedefs. */
  typedef Affine3DTransform Self;
  typedef Transform< TScalarType, 3, 3 > Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( Affine3DTransform, Transform );

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro( Self );

  /** Dimension of the space. */
  enum  { SpaceDimension  = 3 };
  enum  { ParametersDimension  = 12 };

  /** */
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);

  /** Scalar type. */
  typedef typename Superclass::ScalarType  ScalarType;

  /** Parameters type. */
  typedef typename Superclass::ParametersType  ParametersType;

  /** Jacobian type. */
  typedef typename Superclass::JacobianType  JacobianType;

  /// Standard matrix type for this class
  typedef Matrix<ScalarType, SpaceDimension, SpaceDimension> MatrixType;

  /// Standard vector type for this class
  typedef Vector<TScalarType, SpaceDimension> OffsetType;

  /// Standard vector type for this class
  typedef Vector<TScalarType, SpaceDimension> InputVectorType;
  typedef Vector<TScalarType, SpaceDimension> OutputVectorType;

  /// Standard covariant vector type for this class
  typedef CovariantVector<TScalarType, SpaceDimension> InputCovariantVectorType;
  typedef CovariantVector<TScalarType, SpaceDimension> OutputCovariantVectorType;

  /// Standard vnl_vector type for this class
  typedef vnl_vector_fixed<TScalarType, SpaceDimension> InputVnlVectorType;
  typedef vnl_vector_fixed<TScalarType, SpaceDimension> OutputVnlVectorType;

  /// Standard coordinate point type for this class
  typedef Point<TScalarType, SpaceDimension>    InputPointType;
  typedef Point<TScalarType, SpaceDimension>    OutputPointType;

  /** Standard offset type for this class   */
  //typedef     OutputVectorType    OffsetType;


  /**
   * Get offset of an Affine3DTransform
   **/
  itkGetConstReferenceMacro( Offset, OffsetType );

  /**
   * Get matrix from an Affine3DTransform
   **/
  itkGetConstReferenceMacro( Matrix, MatrixType );


  /**
   * Set offset of a Affine3D Transform
   **/
  void SetOffset(const OffsetType &offset)
  {
    m_Offset = offset;
    this->ConstructParametersFromTransform();
    this->Modified();
    return;
  }

  /**
   * Set the matrix of a Rigid3D Transform
   **/
  void SetMatrix(const MatrixType &matrix)
  {
    m_Matrix = matrix;
    this->ConstructParametersFromTransform();
    this->Modified();
    return;
  }


  /**
   * Set the transformation from a container of parameters
  **/
  void SetParameters( const ParametersType & parameters );

  /** Get the Transformation Parameters. */
  itkGetConstReferenceMacro( Parameters, ParametersType );


  /**
   * Transform points, vectors, etc.
   **/
  OutputPointType     TransformPoint(const InputPointType  &point ) const;
  OutputVectorType    TransformVector(const InputVectorType &vector) const;
  OutputVnlVectorType    TransformVector(const InputVnlVectorType &vector) const;

  OutputCovariantVectorType TransformCovariantVector(
    const InputCovariantVectorType &vector) const;

  /**
   * Print contents of an Affine3DTransform
   **/
  void PrintSelf(std::ostream &os, Indent indent) const;

  /**
   * Find inverse of an affine transformation
   *
   * This method creates and returns a new Affine3DTransform object
   * which is the inverse of self.  If self is not invertible,
   * an exception is thrown.
   **/
  Pointer GetInverse( void ) const;

  /** Set the parameters to the IdentityTransform */
  void SetIdentity(void);

  /** Compute the Jacobian Matrix of the transformation at one point */
  const JacobianType & GetJacobian(const InputPointType  &point ) const;

protected:
  Affine3DTransform();
  ~Affine3DTransform();


private:
  Affine3DTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Function to calculate determinant of m_Matrix
  double GetDeterminantOfMatrix( void ) const;

  // Function to provide a set of parameters that generate the
  // current transform. Since such a set is not unique, the extra
  // constraint is introduced that all angles lie in the interval
  // [ -pi/2, pi/2 ]
  void ConstructParametersFromTransform( void );


  MatrixType          m_Matrix;
  OffsetType          m_Offset;

}; //class Affine3DTransform


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAffine3DTransform.txx"
#endif

#endif /* __itkAffine3DTransform_h */
