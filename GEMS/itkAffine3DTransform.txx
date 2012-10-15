#ifndef _itkAffine3DTransform_txx
#define _itkAffine3DTransform_txx

#include "itkAffine3DTransform.h"
#include "vnl/vnl_math.h"


namespace itk
{

// Constructor with default arguments
template<class TScalarType>
Affine3DTransform<TScalarType>::
Affine3DTransform():Superclass(OutputSpaceDimension,ParametersDimension)
{
  m_Offset.Fill( 0 );
  m_Matrix.SetIdentity();
  
  //m_Jacobian = JacobianType( SpaceDimension, ParametersDimension );
  this->m_Jacobian.Fill( NumericTraits< ITK_TYPENAME JacobianType::ValueType >::Zero );
  
  this->m_Parameters = ParametersType( ParametersDimension );
  this->m_Parameters.Fill( NumericTraits< ITK_TYPENAME ParametersType::ValueType >::Zero );
  for (unsigned int i=6; i<9; i++)
    {
    this->m_Parameters[i] = 1.0;
    }


}


// Destructor
template<class TScalarType>
Affine3DTransform<TScalarType>::
~Affine3DTransform()
{
}


// Print self
template<class TScalarType>
void
Affine3DTransform<TScalarType>::
PrintSelf(std::ostream &os, Indent indent) const
{

  Superclass::PrintSelf(os,indent);

  os << indent << "Offset: " << m_Offset   << std::endl;
  os << indent << "Matrix: " << m_Matrix   << std::endl;
}



  /**
   * Set the transformation from a container of parameters
  **/
template<class TScalarType>
void
Affine3DTransform<TScalarType>::
SetParameters( const ParametersType & parameters )
{
  // First test for decent input
  if ( parameters.Size() != ParametersDimension )
    {
      itkExceptionMacro( << "Size of parameters must be " << ParametersDimension <<  "!" );
    }

  this->m_Parameters = parameters;

  double sinX = sin( this->m_Parameters[ 0 ] );
  double cosX = cos( this->m_Parameters[ 0 ] );
  double sinY = sin( this->m_Parameters[ 1 ] );
  double cosY = cos( this->m_Parameters[ 1 ] );
  double sinZ = sin( this->m_Parameters[ 2 ] );
  double cosZ = cos( this->m_Parameters[ 2 ] );
  double translationX = this->m_Parameters[ 3 ];
  double translationY = this->m_Parameters[ 4 ];
  double translationZ = this->m_Parameters[ 5 ];
  double scaleX = this->m_Parameters[ 6 ];
  double scaleY = this->m_Parameters[ 7 ];
  double scaleZ = this->m_Parameters[ 8 ];
  double gantryX = this->m_Parameters[ 9 ];
  double gantryY = this->m_Parameters[ 10 ];
  double gantryZ = this->m_Parameters[ 11 ];

  m_Matrix[ 0 ][ 0 ] = scaleX * cosY * ( cosZ + gantryY * sinZ );
  m_Matrix[ 0 ][ 1 ] = scaleY * ( cosY * ( sinZ + gantryX * gantryZ * cosZ ) - gantryZ * sinY );
  m_Matrix[ 0 ][ 2 ] = scaleZ * ( gantryX * cosY * cosZ - sinY );
  m_Matrix[ 1 ][ 0 ] = scaleX * ( sinX * sinY * ( cosZ + gantryY * sinZ ) - cosX * ( sinZ - gantryY * cosZ ) );
  m_Matrix[ 1 ][ 1 ] = scaleY * ( sinX * sinY * ( sinZ + gantryX * gantryZ * cosZ ) +
                                  cosX * ( cosZ - gantryX * gantryZ * sinZ ) + gantryZ * sinX * cosY );
  m_Matrix[ 1 ][ 2 ] = scaleZ * ( sinX * cosY + gantryX * ( sinX * sinY * cosZ - cosX * sinZ ) );
  m_Matrix[ 2 ][ 0 ] = scaleX * ( cosX * sinY * ( cosZ + gantryY * sinZ ) + sinX * ( sinZ - gantryY * cosZ ) );
  m_Matrix[ 2 ][ 1 ] = scaleY * ( cosX * sinY * ( sinZ + gantryX * gantryZ * cosZ ) -
                                  sinX * ( cosZ - gantryX * gantryZ * sinZ ) + gantryZ * cosX * cosY );
  m_Matrix[ 2 ][ 2 ] = scaleZ * ( cosX * cosY + gantryX * ( cosX * sinY * cosZ + sinX * sinZ ) );

  m_Offset[ 0 ] = translationX;
  m_Offset[ 1 ] = translationY;
  m_Offset[ 2 ] = translationZ;

  this->Modified();
}



// Transform a point
template<class TScalarType>
typename Affine3DTransform<TScalarType>::OutputPointType
Affine3DTransform<TScalarType>::
TransformPoint(const InputPointType &point) const
{
  return m_Matrix * point + m_Offset;
}


// Transform a vector
template<class TScalarType>
typename Affine3DTransform<TScalarType>::OutputVectorType
Affine3DTransform<TScalarType>::
TransformVector(const InputVectorType &vect) const
{
  return  m_Matrix * vect;
}


// Transform a vnl_vector_fixed
template<class TScalarType>
typename Affine3DTransform<TScalarType>::OutputVnlVectorType
Affine3DTransform<TScalarType>::
TransformVector(const InputVnlVectorType &vect) const
{
  return  m_Matrix * vect;
}


// Transform a CovariantVector
template<class TScalarType>
typename Affine3DTransform<TScalarType>::OutputCovariantVectorType
Affine3DTransform<TScalarType>::
TransformCovariantVector(const InputCovariantVectorType &vect) const
{
  // Covariant vectors are transformed like contravariant
  // vectors under orthogonal transformations.
  return  m_Matrix * vect;
}


// Create and return an inverse transformation
template<class TScalarType>
typename Affine3DTransform<TScalarType>::Pointer
Affine3DTransform<TScalarType>::
GetInverse( void ) const
{

  // If matrix singular, throw an exception
  if ( this->GetDeterminantOfMatrix() <= NumericTraits< ScalarType >::epsilon() )
    {
      itkExceptionMacro( << "Failed to invert singular matrix!" );
    }

  MatrixType inverse;
  inverse = m_Matrix.GetInverse();
  Pointer result = New();
  result->m_Matrix   =  inverse;
  result->m_Offset   = -(inverse * m_Offset);
  result->ConstructParametersFromTransform();
  return result;
}



// Set the transform to identity
template<class TScalarType >
void
Affine3DTransform< TScalarType >::
SetIdentity( void )
{
  ParametersType parameters( ParametersDimension );
  parameters.Fill( NumericTraits< ITK_TYPENAME ParametersType::ValueType >::Zero );
  for (unsigned int i=6; i<9; i++)
    {
      parameters[i] = 1.0;
    }
  
  this->SetParameters( parameters );

}



// Compute the Jacobian in one position
template<class TScalarType >
const typename Affine3DTransform<TScalarType>::JacobianType &
Affine3DTransform< TScalarType >::
GetJacobian( const InputPointType & p ) const
{

  // Some constants to use...
  double sinX = sin( this->m_Parameters[ 0 ] );
  double cosX = cos( this->m_Parameters[ 0 ] );
  double sinY = sin( this->m_Parameters[ 1 ] );
  double cosY = cos( this->m_Parameters[ 1 ] );
  double sinZ = sin( this->m_Parameters[ 2 ] );
  double cosZ = cos( this->m_Parameters[ 2 ] );
  //double translationX = this->m_Parameters[ 3 ];
  //double translationY = this->m_Parameters[ 4 ];
  //double translationZ = this->m_Parameters[ 5 ];
  double scaleX = this->m_Parameters[ 6 ];
  double scaleY = this->m_Parameters[ 7 ];
  double scaleZ = this->m_Parameters[ 8 ];
  double gantryX = this->m_Parameters[ 9 ];
  double gantryY = this->m_Parameters[ 10 ];
  double gantryZ = this->m_Parameters[ 11 ];

  // First row
  this->m_Jacobian[0][0] = 0;
  this->m_Jacobian[0][1] = -scaleX * sinY * ( cosZ + gantryY * sinZ ) * p[0]
                             -scaleY * ( sinY * ( sinZ + gantryX * gantryZ * cosZ ) + gantryZ * cosY ) * p[1]
                             -scaleZ * ( gantryX * sinY * cosZ + cosY ) * p[2];
  this->m_Jacobian[0][2] = scaleX * cosY * ( -sinZ + gantryY * cosZ ) * p[0] +
                             scaleY * cosY * ( cosZ - gantryX * gantryZ * sinZ ) * p[1] -
                             scaleZ * gantryX * cosY * sinZ * p[2];
  this->m_Jacobian[0][3] = 1;
  this->m_Jacobian[0][4] = 0;
  this->m_Jacobian[0][5] = 0;
  this->m_Jacobian[0][6] = cosY * ( cosZ + gantryY * sinZ ) * p[0];
  this->m_Jacobian[0][7] = ( cosY * ( sinZ + gantryX * gantryZ * cosZ ) - gantryZ * sinY ) * p[1];
  this->m_Jacobian[0][8] = ( gantryX * cosY * cosZ - sinY ) * p[2];
  this->m_Jacobian[0][9] = scaleY * cosY * gantryZ * cosZ * p[1] +
                             scaleZ * cosY * cosZ * p[2];
  this->m_Jacobian[0][10] = scaleX * cosY * sinZ * p[0];
  this->m_Jacobian[0][11] = scaleY * ( cosY * gantryX * cosZ - sinY ) * p[1];

  // Second row
  this->m_Jacobian[1][0] = scaleX * ( cosX * sinY * ( cosZ + gantryY * sinZ ) + sinX * ( sinZ - gantryY * cosZ ) ) * p[0] +
                             scaleY * ( cosX * sinY * ( sinZ + gantryX * gantryZ * cosZ ) - 
                                            sinX * ( cosZ - gantryX * gantryZ * sinZ ) + gantryZ * cosX * cosY ) * p[1] +
                             scaleZ * ( cosX * cosY + gantryX * ( cosX * sinY * cosZ + sinX * sinZ ) ) * p[2];
  this->m_Jacobian[1][1] = scaleX * ( sinX * cosY * ( cosZ + gantryY * sinZ ) ) * p[0] +
                             scaleY * ( sinX * cosY * ( sinZ + gantryX * gantryZ * cosZ ) - gantryZ * sinX * sinY ) * p[1] +
                             scaleZ * ( -sinX * sinY + gantryX *sinX * cosY * cosZ ) * p[2];
  this->m_Jacobian[1][2] = scaleX * ( sinX * sinY * ( -sinZ + gantryY * cosZ ) - cosX * ( cosZ + gantryY * sinZ ) ) * p[0] +
                             scaleY * ( sinX * sinY * ( cosZ - gantryX * gantryZ * sinZ ) + 
                                            cosX * ( -sinZ - gantryX * gantryZ * cosZ ) ) * p[1] +
                             scaleZ * gantryX * ( -sinX * sinY * sinZ - cosX * cosZ ) * p[2];
  this->m_Jacobian[1][3] = 0;
  this->m_Jacobian[1][4] = 1;
  this->m_Jacobian[1][5] = 0;
  this->m_Jacobian[1][6] = ( sinX * sinY * ( cosZ + gantryY * sinZ ) - cosX * ( sinZ - gantryY * cosZ ) ) * p[0];
  this->m_Jacobian[1][7] = ( sinX * sinY * ( sinZ + gantryX * gantryZ * cosZ ) + cosX * ( cosZ - gantryX * gantryZ * sinZ )
                               + gantryZ * sinX * cosY ) * p[1];
  this->m_Jacobian[1][8] = ( sinX * cosY + gantryX * ( sinX * sinY * cosZ - cosX * sinZ ) ) * p[2];
  this->m_Jacobian[1][9] = scaleY * ( sinX * sinY * gantryZ * cosZ - cosX * gantryZ * sinZ ) * p[1] +
                             scaleZ * ( sinX * sinY * cosZ - cosX * sinZ ) * p[2];
  this->m_Jacobian[1][10] = scaleX * ( sinX * sinY * sinZ + cosX * cosZ ) * p[0];
  this->m_Jacobian[1][11] = scaleY * ( sinX * sinY * gantryX * cosZ - cosX * gantryX * sinZ + sinX * cosY ) * p[1];

  // Third row
  this->m_Jacobian[2][0] = scaleX * ( -sinX * sinY * ( cosZ + gantryY * sinZ ) + cosX * ( sinZ - gantryY * cosZ ) ) * p[0] +
                             scaleY * ( -sinX * sinY * ( sinZ + gantryX * gantryZ * cosZ )
                                            - cosX * ( cosZ - gantryX * gantryZ * sinZ )  - gantryZ * sinX * cosY ) * p[1] +
                             scaleZ * ( -sinX * cosY + gantryX * ( -sinX * sinY * cosZ + cosX * sinZ ) ) * p[2];
  this->m_Jacobian[2][1] = scaleX * cosX * cosY * ( cosZ + gantryY * sinZ ) * p[0] +
                              scaleY * ( cosX * cosY * ( sinZ + gantryX * gantryZ * cosZ ) - gantryZ * cosX * sinY ) * p[1] +
                              scaleZ * ( -cosX * sinY + gantryX * cosX * cosY * cosZ ) * p[2];
  this->m_Jacobian[2][2] = scaleX * ( cosX * sinY * ( -sinZ + gantryY * cosZ ) + sinX * ( cosZ + gantryY * sinZ ) ) * p[0] +
                             scaleY * ( cosX * sinY * ( cosZ - gantryX * gantryZ * sinZ ) - 
                                            sinX * ( -sinZ - gantryX * gantryZ * cosZ ) ) * p[1] +
                             scaleZ * gantryX * ( -cosX * sinY * sinZ + sinX * cosZ ) * p[2];
  this->m_Jacobian[2][3] = 0;
  this->m_Jacobian[2][4] = 0;
  this->m_Jacobian[2][5] = 1;
  this->m_Jacobian[2][6] = ( cosX * sinY * ( cosZ + gantryY * sinZ ) + sinX * ( sinZ - gantryY * cosZ ) ) * p[0];
  this->m_Jacobian[2][7] = ( cosX * sinY * ( sinZ + gantryX * gantryZ * cosZ ) - sinX * ( cosZ - gantryX * gantryZ * sinZ )
                               + gantryZ * cosX * cosY ) * p[1];
  this->m_Jacobian[2][8] = ( cosX * cosY + gantryX * ( cosX * sinY * cosZ + sinX * sinZ ) ) * p[2];
  this->m_Jacobian[2][9] = scaleY * ( cosX * sinY * gantryZ * cosZ + sinX * gantryZ * sinZ ) * p[1] +
                             scaleZ * ( cosX * sinY * cosZ + sinX * sinZ ) * p[2];
  this->m_Jacobian[2][10] = scaleX * ( cosX * sinY * sinZ - sinX * cosZ ) * p[0];
  this->m_Jacobian[2][11] = scaleY * ( cosX * sinY * gantryX * cosZ + sinX * gantryX * sinZ + cosX * cosY ) * p[1];


  return this->m_Jacobian;

}



// Construct a set of parameters that would result in the current transformation
template<class TScalarType >
void 
Affine3DTransform< TScalarType >::
ConstructParametersFromTransform( void )
{
  MatrixType transpose;
  transpose = m_Matrix.GetTranspose();
  MatrixType p = transpose * m_Matrix;
  double determinant = GetDeterminantOfMatrix();

  double scaleYSquare = p[1][1] - vnl_math_sqr( p[1][2] ) / p[2][2];
  double scaleXSquare = p[0][0] - vnl_math_sqr( p[0][1] * p[2][2] - p[0][2] * p[1][2]) / vnl_math_sqr( p[2][2] ) / scaleYSquare;
  double scaleZSquare = p[2][2] - vnl_math_sqr( p[0][2] ) / scaleXSquare;

  double scaleXAbs = vcl_sqrt(  scaleXSquare );
  double scaleYAbs = vcl_sqrt(  scaleYSquare );
  double scaleZAbs = vcl_sqrt(  scaleZSquare );

  double gantryXAbs = p[0][2] * scaleYAbs / determinant;
  double gantryYAbs = ( p[0][1] * p[2][2] - p[0][2] * p[1][2] ) / p[2][2] * scaleZAbs / determinant;
  double gantryZAbs = p[1][2] / p[2][2] * scaleZSquare * scaleXAbs / determinant;

  double s1 = m_Matrix[0][0] / scaleXAbs -
                   gantryYAbs * m_Matrix[0][1] / scaleYAbs +
                   gantryYAbs * gantryZAbs * m_Matrix[0][2] / scaleZAbs;
  double s2 = - gantryXAbs * m_Matrix[2][0] / scaleXAbs
                    + gantryXAbs * gantryYAbs * m_Matrix[2][1] / scaleYAbs
                    + ( 1 - gantryXAbs * gantryYAbs * gantryZAbs ) * m_Matrix[2][2] / scaleZAbs;

  double signScaleX = vnl_math_sgn0( s1 );
  double signScaleY = vnl_math_sgn0( s1 ) * vnl_math_sgn0( s2 ) * vnl_math_sgn0( determinant );
  double signScaleZ = vnl_math_sgn0( s2 );

  double scaleX = signScaleX * scaleXAbs;
  double scaleY = signScaleY * scaleYAbs;
  double scaleZ = signScaleZ * scaleZAbs;

  double gantryX = signScaleY * gantryXAbs;
  double gantryY = signScaleZ * gantryYAbs;
  double gantryZ = signScaleX * gantryZAbs;

  MatrixType scale;
  scale.Fill( 0.0 );
  scale[ 0 ][ 0 ] = scaleX;
  scale[ 1 ][ 1 ] = scaleY;
  scale[ 2 ][ 2 ] = scaleZ;
//   std::cout << "scale: " << scale << std::endl;

  MatrixType gantry;
  gantry.Fill( 0.0 );
  gantry.SetIdentity();
  gantry[ 0 ][ 1 ] = gantryX * gantryZ;
  gantry[ 0 ][ 2 ] = gantryX;
  gantry[ 1 ][ 0 ] = gantryY;
  gantry[ 2 ][ 1 ] = gantryZ;
//   std::cout << "gantry: " << gantry << std::endl;

  MatrixType rotation;
  rotation = m_Matrix * scale.GetInverse() * gantry.GetInverse();

  double rotationY = - asin( rotation[0][2] );
  double rotationX = asin( rotation[1][2] / cos( rotationY ) );
  double rotationZ = asin( rotation[0][1] / cos( rotationY ) );

  // Fill in in this->m_Parameters
  this->m_Parameters[0] = rotationX;
  this->m_Parameters[1] = rotationY;
  this->m_Parameters[2] = rotationZ;
  this->m_Parameters[3] = m_Offset[0];
  this->m_Parameters[4] = m_Offset[1];
  this->m_Parameters[5] = m_Offset[2];
  this->m_Parameters[6] = scaleX;
  this->m_Parameters[7] = scaleY;
  this->m_Parameters[8] = scaleZ;
  this->m_Parameters[9] = gantryX;
  this->m_Parameters[10] = gantryY;
  this->m_Parameters[11] = gantryZ;

}




// Compute the determinant of the matrix
template<class TScalarType >
double
Affine3DTransform< TScalarType >::
GetDeterminantOfMatrix( void ) const
{
  double determinant = m_Matrix[ 0 ][ 0 ] * ( m_Matrix[ 1 ][ 1 ] * m_Matrix[ 2 ][ 2 ] - m_Matrix[ 2 ][ 1 ] * m_Matrix[ 1 ][ 2 ] ) -
                                m_Matrix[ 0 ][ 1 ] * ( m_Matrix[ 1 ][ 0 ] * m_Matrix[ 2 ][ 2 ] - m_Matrix[ 2 ][ 0 ] * m_Matrix[ 1 ][ 2 ] ) +
                                m_Matrix[ 0 ][ 2 ] * ( m_Matrix[ 1 ][ 0 ] * m_Matrix[ 2 ][ 1 ] - m_Matrix[ 2 ][ 0 ] * m_Matrix[ 1 ][ 1 ] );

  return determinant;
}




} // namespace

#endif
