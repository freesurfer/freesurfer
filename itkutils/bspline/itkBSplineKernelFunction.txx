/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkBSplineKernelFunction_txx
#define _itkBSplineKernelFunction_txx

#include "itkBSplineKernelFunction.h"

namespace itk
{

template <unsigned int VSplineOrder>
BSplineKernelFunction<VSplineOrder>
::BSplineKernelFunction() 
{ 
  this->m_SplineOrder = VSplineOrder;
  this->GenerateBSplineShapeFunctions( this->m_SplineOrder+1 );
}

template <unsigned int VSplineOrder>
BSplineKernelFunction<VSplineOrder>
::~BSplineKernelFunction() 
{ 
}

template <unsigned int VSplineOrder>
void
BSplineKernelFunction<VSplineOrder>
::SetSplineOrder( unsigned int order )
{
  if ( order != this->m_SplineOrder )
    {
    this->m_SplineOrder = order;
    this->GenerateBSplineShapeFunctions( this->m_SplineOrder+1 );
    this->Modified();
    }
}

template <unsigned int VSplineOrder>
void
BSplineKernelFunction<VSplineOrder>
::GenerateBSplineShapeFunctions( unsigned int order )
{
  unsigned int NumberOfPieces = static_cast<unsigned int>( 0.5*( order+1 ) );
  this->m_BSplineShapeFunctions.set_size( NumberOfPieces, order );

  VectorType knots( order+1 );
  for ( unsigned int i = 0; i < knots.size(); i++)
    {
    knots[i] = -0.5*static_cast<RealType>( order ) + static_cast<RealType>( i );
    }

  for ( unsigned int i = 0; i < NumberOfPieces; i++ )
    {
    PolynomialType poly = this->CoxDeBoor(order, knots, 
             0, static_cast<unsigned int>( 0.5*( order ) ) + i );
    this->m_BSplineShapeFunctions.set_row( i, poly.coefficients() );
    }   
}

template <unsigned int VSplineOrder>
typename BSplineKernelFunction<VSplineOrder>::PolynomialType 
BSplineKernelFunction<VSplineOrder>
::CoxDeBoor( unsigned short order, VectorType knots, 
             unsigned int whichBasisFunction, unsigned int whichPiece )
{
  VectorType tmp(2);
  PolynomialType poly1(0.0), poly2(0.0);
  RealType den;
  unsigned short p = order-1;
  unsigned short i = whichBasisFunction;   

  if ( p == 0 && whichBasisFunction == whichPiece )
    {
    PolynomialType poly(1.0);
    return poly;
    }          

  // Term 1
  den = knots(i+p)-knots(i);
  if ( den == NumericTraits<RealType>::Zero )  
    {
    PolynomialType poly(0.0);
    poly1 = poly;
    }
  else
    {
    tmp(0) = 1.0;
    tmp(1) = -knots(i);
    tmp /= den;
    poly1 = PolynomialType(tmp) * this->CoxDeBoor( order-1, knots, i, whichPiece );   
    }

  // Term 2
  den = knots(i+p+1)-knots(i+1);
  if ( den == NumericTraits<RealType>::Zero )  
    {
    PolynomialType poly(0.0);
    poly2 = poly;
    }
  else
    {
    tmp(0) = -1.0;
    tmp(1) = knots(i+p+1);
    tmp /= den;
    poly2 = PolynomialType(tmp) * this->CoxDeBoor( order-1, knots, i+1, whichPiece );
    }    
  return ( poly1 + poly2 );
}

template <unsigned int VSplineOrder>
typename BSplineKernelFunction<VSplineOrder>::MatrixType 
BSplineKernelFunction<VSplineOrder>
::GetShapeFunctionsInZeroToOneInterval()
{
  int order = this->m_SplineOrder+1;
  unsigned int NumberOfPieces = static_cast<unsigned int>( order );
  MatrixType ShapeFunctions( NumberOfPieces, order );

  VectorType knots( 2*order );
  for ( unsigned int i = 0; i < knots.size(); i++ )
    {
    knots[i] = -static_cast<RealType>( this->m_SplineOrder ) 
               + static_cast<RealType>( i );
    }

  for ( unsigned int i = 0; i < NumberOfPieces; i++ )
    {
    PolynomialType poly = this->CoxDeBoor( order, knots, i, order-1 );
    ShapeFunctions.set_row( i, poly.coefficients() ); 
    }   
  return ShapeFunctions;
}

template <unsigned int VSplineOrder>
void
BSplineKernelFunction<VSplineOrder>
::PrintSelf( std::ostream& os, Indent indent ) const
{ 
  Superclass::PrintSelf( os, indent ); 
  os << indent  << "Spline Order: " << this->m_SplineOrder << std::endl;
  os << indent  << "Piecewise Polynomial Pieces: " << std::endl;  
  RealType a;
  RealType b = 0;
  for ( unsigned int i = 0; i < this->m_BSplineShapeFunctions.rows(); i++ )
    {
    os << indent << indent;
    PolynomialType( this->m_BSplineShapeFunctions.get_row( i ) ).print( os );
    if ( i == 0 )
      {
      a = 0.0;
      if ( this->m_SplineOrder % 2 == 0 )
        {
        b = 0.5;
        }
      else
        {
        b = 1.0;
        }
      }
    else
      {
      a = b;
      b += 1.0;
      }  
    os << ",  X \\in [" << a << ", " << b << "]" << std::endl;
    }  
}  


} // namespace itk

#endif
