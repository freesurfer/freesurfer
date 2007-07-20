#ifndef __itkBSplineDataPointSetToImageFilter_txx
#define __itkBSplineDataPointSetToImageFilter_txx

#include "itkBSplineDataPointSetToImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageDuplicator.h"
#include "itkCastImageFilter.h"
#include "itkNumericTraits.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_vector.h"

namespace itk
{

template <class TInputPointSet, class TOutputImage>
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::BSplineDataPointSetToImageFilter()
{
  this->m_SplineOrder.Fill( 3 );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    { 
    this->m_NumberOfControlPoints[i] = ( this->m_SplineOrder[i]+1 );
    this->m_Kernel[i] = KernelType::New();
    this->m_Kernel[i]->SetSplineOrder( this->m_SplineOrder[i] );
    }

  this->m_CloseDimension.Fill( 0 );
  this->m_DoMultilevel = false;
  this->m_GenerateOutputImage = true;
  this->m_NumberOfLevels.Fill( 1 );
  this->m_MaximumNumberOfLevels = 1;
   
  this->m_PhiLattice = PointDataImageType::New();
  this->m_PsiLattice = PointDataImageType::New();  
  this->m_InputPointData = PointDataContainerType::New();
  this->m_OutputPointData = PointDataContainerType::New();
  
  this->m_PointWeights = WeightsContainerType::New();
  this->m_UsePointWeights = false;
}

template <class TInputPointSet, class TOutputImage>
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::~BSplineDataPointSetToImageFilter()
{  
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetSplineOrder( unsigned int order )
{ 
  this->m_SplineOrder.Fill( order );
  this->SetSplineOrder( this->m_SplineOrder );
} 

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetSplineOrder( ArrayType order )
{ 
  itkDebugMacro(<< "Setting m_SplineOrder to " << order);

  this->m_SplineOrder = order;
  for ( int i = 0; i < ImageDimension; i++ )
    { 
    if ( this->m_SplineOrder[i] == 0 )
      {
      itkExceptionMacro( << "The spline order in each dimension must be greater than 0" );
      }

    this->m_Kernel[i] = KernelType::New();
    this->m_Kernel[i]->SetSplineOrder( this->m_SplineOrder[i] );
  
    if ( this->m_DoMultilevel )
      {
      typename KernelType::MatrixType R;
      typename KernelType::MatrixType C;
      C = this->m_Kernel[i]->GetShapeFunctionsInZeroToOneInterval();
      R = C;
      for ( unsigned int j = 0; j < C.cols(); j++ )
        {
        RealType c = pow( 2.0, static_cast<RealType>( C.cols()-j-1 ) );
        for ( unsigned int k = 0; k < C.rows(); k++)
          {
          R(k, j) *= c;
          }
        }  
      R = R.transpose();  
      R.flipud();
      C = C.transpose();  
      C.flipud();      
      this->m_RefinedLatticeCoefficients[i] 
          = ( vnl_svd<RealType>( R ).solve( C ) ).extract( 2, C.cols() );
      }  
    } 
  this->Modified();
}  

template <class TInputPointSet, class TOutputImage>
void 
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetNumberOfLevels( unsigned int levels )
{
  this->m_NumberOfLevels.Fill( levels );
  this->SetNumberOfLevels( this->m_NumberOfLevels );
}     

template <class TInputPointSet, class TOutputImage>
void 
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetNumberOfLevels( ArrayType levels )
{
  this->m_NumberOfLevels = levels;
  this->m_MaximumNumberOfLevels = 1;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( this->m_NumberOfLevels[i] == 0 )
      {
      itkExceptionMacro( << "The number of levels in each dimension must be greater than 0" );
      }
    if ( this->m_NumberOfLevels[i] > this->m_MaximumNumberOfLevels )
      {
      this->m_MaximumNumberOfLevels = this->m_NumberOfLevels[i];
      } 
    }

  itkDebugMacro( << "Setting m_NumberOfLevels to " 
                 << this->m_NumberOfLevels );
  itkDebugMacro( << "Setting m_MaximumNumberOfLevels to " 
                 << this->m_MaximumNumberOfLevels );

  if ( this->m_MaximumNumberOfLevels > 1 ) 
    {
    this->m_DoMultilevel = true;
    }
  else
    {  
    this->m_DoMultilevel = false;
    }
  this->SetSplineOrder( this->m_SplineOrder );
  this->Modified();
}     

template <class TInputPointSet, class TOutputImage>
void 
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetPointWeights( typename WeightsContainerType::Pointer weights )
{
  this->m_UsePointWeights = true; 
  this->m_PointWeights = weights;
  this->Modified();
}
  
template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::GenerateData()
{

  Superclass::GenerateData();

  this->m_InputPointData->Initialize();
  this->m_OutputPointData->Initialize();

  if ( this->m_UsePointWeights && 
       ( this->m_PointWeights->Size() != this->GetInput()->GetNumberOfPoints() ) )
    {
    itkExceptionMacro( << "The number of weight points and input points must be equal." );   
    }

  typename PointSetType::PointDataContainer::ConstIterator It;
  It = this->GetInput()->GetPointData()->Begin();
  while ( It != this->GetInput()->GetPointData()->End() )
    {
    if ( !this->m_UsePointWeights )
      {
      this->m_PointWeights->InsertElement( It.Index(), 1.0 );
      }
    this->m_InputPointData->InsertElement( It.Index(), It.Value() );
    this->m_OutputPointData->InsertElement( It.Index(), It.Value() );
    ++It;
    }

  this->m_CurrentLevel = 0; 
  this->m_CurrentNumberOfControlPoints = this->m_NumberOfControlPoints;    
  this->GenerateControlLattice();
  this->UpdatePointSet();

  itkDebugMacro( << "Current Level = " << this->m_CurrentLevel );
  itkDebugMacro( << "  Number of control points = " 
                 << this->m_CurrentNumberOfControlPoints );

  if ( this->m_DoMultilevel ) 
    {   
    this->m_PsiLattice->SetRegions( this->m_PhiLattice->GetLargestPossibleRegion() );
    this->m_PsiLattice->Allocate();
    PointDataType P( 0.0 );
    this->m_PsiLattice->FillBuffer( P );
    }
  
  for ( this->m_CurrentLevel = 1; 
        this->m_CurrentLevel < this->m_MaximumNumberOfLevels; this->m_CurrentLevel++ )
    {
  
    ImageRegionIterator<PointDataImageType> ItPsi( this->m_PsiLattice, 
            this->m_PsiLattice->GetLargestPossibleRegion() );
    ImageRegionIterator<PointDataImageType> ItPhi( this->m_PhiLattice, 
            this->m_PhiLattice->GetLargestPossibleRegion() );
    for ( ItPsi.GoToBegin(), ItPhi.GoToBegin(); 
          !ItPsi.IsAtEnd(); ++ItPsi, ++ItPhi )
      {
      ItPsi.Set(ItPhi.Get()+ItPsi.Get());
      }  
    this->RefineControlLattice();

    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if ( this->m_CurrentLevel < this->m_NumberOfLevels[i] )
        {
        this->m_CurrentNumberOfControlPoints[i] = 
      2*this->m_CurrentNumberOfControlPoints[i]-this->m_SplineOrder[i];
        }
      }

    itkDebugMacro( << "Current Level = " << this->m_CurrentLevel );
    itkDebugMacro( << "  No. of control points = " 
                   << this->m_CurrentNumberOfControlPoints );

    RealType avg_p = NumericTraits<RealType>::Zero;

    typename PointDataContainerType::Iterator  ItIn, ItOut;
    ItIn = this->m_InputPointData->Begin();
    ItOut = this->m_OutputPointData->Begin();
    while ( ItIn != this->m_InputPointData->End() )
      {
      this->m_InputPointData->InsertElement( ItIn.Index(), ItIn.Value()-ItOut.Value() );

      if ( this->GetDebug() )
       { 
       avg_p += ( ItIn.Value()-ItOut.Value() ).GetNorm(); 
        }  
      ++ItIn;
      ++ItOut;
      }
    itkDebugMacro( << "The average difference norm of the point set is " 
                   << avg_p / static_cast<RealType>( this->GetInput()->GetNumberOfPoints() ) );

    this->GenerateControlLattice();
    this->UpdatePointSet();
    }

  if (this->m_DoMultilevel) 
    {   
    ImageRegionIterator<PointDataImageType> 
      ItPsi( this->m_PsiLattice, this->m_PsiLattice->GetLargestPossibleRegion() );
    ImageRegionIterator<PointDataImageType> 
      ItPhi( this->m_PhiLattice, this->m_PhiLattice->GetLargestPossibleRegion() );
    for ( ItPsi.GoToBegin(), ItPhi.GoToBegin(); 
          !ItPsi.IsAtEnd(); ++ItPsi, ++ItPhi )
      {
      ItPsi.Set( ItPhi.Get()+ItPsi.Get() );
      } 

    typedef ImageDuplicator<PointDataImageType> ImageDuplicatorType;
    typename ImageDuplicatorType::Pointer Duplicator = ImageDuplicatorType::New();  
    Duplicator->SetInputImage( this->m_PsiLattice );
    Duplicator->Update();
    this->m_PhiLattice = Duplicator->GetOutput();      
    this->UpdatePointSet();
    }

  if ( this->m_GenerateOutputImage )
    {  
    this->GenerateOutputImage();
    }
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::RefineControlLattice() 
{
  ArrayType NumberOfNewControlPoints = this->m_CurrentNumberOfControlPoints;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( this->m_CurrentLevel < this->m_NumberOfLevels[i] )
      {  
      NumberOfNewControlPoints[i] = 2*NumberOfNewControlPoints[i]-this->m_SplineOrder[i];
      }  
    }
  typename RealImageType::RegionType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( this->m_CloseDimension[i] )
      {
      size[i] = NumberOfNewControlPoints[i] - this->m_SplineOrder[i];
      }
    else
      {
      size[i] = NumberOfNewControlPoints[i];
      }
    }
  
  typename PointDataImageType::Pointer RefinedLattice = PointDataImageType::New();
  RefinedLattice->SetRegions( size );
  RefinedLattice->Allocate();
  PointDataType data;
  data.Fill( 0.0 );
  RefinedLattice->FillBuffer( data );  

  typename PointDataImageType::IndexType idx;
  typename PointDataImageType::IndexType idx_Psi;
  typename PointDataImageType::IndexType tmp;
  typename PointDataImageType::IndexType tmp_Psi;
  typename PointDataImageType::IndexType off;
  typename PointDataImageType::IndexType off_Psi;
  typename PointDataImageType::RegionType::SizeType size_Psi;

  size.Fill( 2 );
  int N = 1;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    N *= ( this->m_SplineOrder[i]+1 );
    size_Psi[i] = this->m_SplineOrder[i] + 1;  
    }

  ImageRegionIteratorWithIndex<PointDataImageType> 
    It( RefinedLattice, RefinedLattice->GetLargestPossibleRegion() );

  It.GoToBegin();
  while ( !It.IsAtEnd() )
    {
    idx = It.GetIndex();
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if ( this->m_CurrentLevel < this->m_NumberOfLevels[i] )
        {
        idx_Psi[i] = static_cast<unsigned int>( 0.5*idx[i] );
       }
      else
        {
        idx_Psi[i] = static_cast<unsigned int>( idx[i] );    
        }
      }
    for ( unsigned int i = 0; 
          i < pow( 2.0, static_cast<RealType>( ImageDimension ) ); i++ )
      {
      PointDataType sum( 0.0 );
      PointDataType val;
      off = this->IndexToSubscript( i, size );

      bool OutOfBoundary = false;
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        tmp[j] = idx[j] + off[j];
        if ( static_cast< unsigned int >( tmp[j] )
           >= NumberOfNewControlPoints[j] && !this->m_CloseDimension[j] )
          {
             OutOfBoundary = true;
             break;
          }
        if ( this->m_CloseDimension[j] )
          {
          tmp[j] %=  RefinedLattice->GetLargestPossibleRegion().GetSize()[j];
          } 
        }
      if ( OutOfBoundary )
        {
        continue;
        }      
 
      for ( int j = 0; j < N; j++ )
        {
        off_Psi = this->IndexToSubscript( j, size_Psi );

        bool OutOfBoundary = false;
        for ( unsigned int k = 0; k < ImageDimension; k++ )
          {
          tmp_Psi[k] = idx_Psi[k] + off_Psi[k];
          if ( static_cast< unsigned int>( tmp_Psi[k] )
             >= this->m_CurrentNumberOfControlPoints[k] 
                  && !this->m_CloseDimension[k] )
            {
            OutOfBoundary = true;
            break;
            }
          if ( this->m_CloseDimension[k] )
            {
           tmp_Psi[k] %= this->m_PsiLattice->GetLargestPossibleRegion().GetSize()[k];
            } 
          }
          if ( OutOfBoundary )
             {
               continue;
             }      
          RealType coeff = 1.0;
          for ( unsigned int k = 0; k < ImageDimension; k++ )
           {
            coeff *= this->m_RefinedLatticeCoefficients[k]( off[k], off_Psi[k] );
           }  
          val = this->m_PsiLattice->GetPixel( tmp_Psi );  
         val *= coeff;
          sum += val;
          }
        RefinedLattice->SetPixel( tmp, sum );
        }  

    bool IsEvenIndex = false;
    while ( !IsEvenIndex && !It.IsAtEnd() )
      {      
      ++It;  
      idx = It.GetIndex();
      IsEvenIndex = true;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        if ( idx[i] % 2 )
          {
          IsEvenIndex = false;
          }
        }
      }    
    }

  typedef ImageDuplicator<PointDataImageType> ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer Duplicator = ImageDuplicatorType::New();  
  Duplicator->SetInputImage( RefinedLattice );
  Duplicator->Update();
  this->m_PsiLattice = Duplicator->GetOutput();  
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::GenerateControlLattice() 
{
  typename RealImageType::RegionType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( this->m_CloseDimension[i] )
      {  
      size[i] = this->m_CurrentNumberOfControlPoints[i]-this->m_SplineOrder[i];
      }
    else
      {
      size[i] = this->m_CurrentNumberOfControlPoints[i];
      }
    }

  this->m_PhiLattice->SetRegions( size );
  this->m_PhiLattice->Allocate();
  this->m_PhiLattice->FillBuffer( 0.0 );

  typename RealImageType::Pointer omega;
  omega = RealImageType::New();
  omega->SetRegions( size );
  omega->Allocate();
  omega->FillBuffer( 0.0 );

  typename PointDataImageType::Pointer delta;
  delta = PointDataImageType::New();  
  delta->SetRegions( size );
  delta->Allocate();
  delta->FillBuffer( 0.0 );
  
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    size[i] = this->m_SplineOrder[i]+1;  
    }

  typename RealImageType::Pointer w;
  w = RealImageType::New();  
  w->SetRegions( size );
  w->Allocate();

  typename PointDataImageType::Pointer phi;
  phi = PointDataImageType::New();  
  phi->SetRegions( size );
  phi->Allocate();

  ImageRegionIteratorWithIndex<RealImageType> 
     Itw( w, w->GetBufferedRegion() );
  ImageRegionIteratorWithIndex<PointDataImageType> 
     Itp( phi, phi->GetBufferedRegion() );

  vnl_vector<RealType> p( ImageDimension );
  vnl_vector<RealType> r( ImageDimension );  
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    r[i] = static_cast<RealType>( this->m_CurrentNumberOfControlPoints[i]-this->m_SplineOrder[i] )
           /static_cast<RealType>( this->m_Size[i]*this->m_Spacing[i]+1e-10 );
    }  
    
  typename PointDataContainerType::ConstIterator It;
  for ( It = this->m_InputPointData->Begin(); 
        It != this->m_InputPointData->End(); ++It )
    {
    
    PointType point;
    this->GetInput()->GetPoint( It.Index(), &point );

    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      p[i] = ( point[i] - this->m_Origin[i] ) * r[i];
      }

    RealType w2_sum = 0.0;
    for ( Itw.GoToBegin(); !Itw.IsAtEnd(); ++Itw )
      {
      RealType B = 1.0;
      typename RealImageType::IndexType idx = Itw.GetIndex();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        { 
    RealType u = static_cast<RealType>( p[i] - 
             static_cast<unsigned>( p[i]) - idx[i] ) 
             + 0.5*static_cast<RealType>( this->m_SplineOrder[i]-1 );
        B *= this->m_Kernel[i]->Evaluate( u );
        }  
      Itw.Set( B );
      w2_sum += B*B;
      }

    for ( Itp.GoToBegin(), Itw.GoToBegin(); !Itp.IsAtEnd(); ++Itp, ++Itw )
      {
      typename RealImageType::IndexType idx = Itw.GetIndex();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        idx[i] += static_cast<unsigned>( p[i] );
        if ( this->m_CloseDimension[i] )
          {
          idx[i] %= delta->GetLargestPossibleRegion().GetSize()[i];
          }  
        }
        RealType wc = this->m_PointWeights->GetElement( It.Index() );
        RealType t = Itw.Get();
        PointDataType data = It.Value();
        data *= ( t/w2_sum ); 
        Itp.Set( data );
        data *= ( t*t*wc );       
        delta->SetPixel( idx, delta->GetPixel( idx ) + data );
        omega->SetPixel( idx, omega->GetPixel( idx ) + wc*t*t );
      }
    }  

  ImageRegionIterator<PointDataImageType> 
      Itl( this->m_PhiLattice, this->m_PhiLattice->GetBufferedRegion() );
  ImageRegionIterator<PointDataImageType> 
      Itd( delta, delta->GetBufferedRegion() );  
  ImageRegionIterator<RealImageType> 
      Ito( omega, omega->GetBufferedRegion() );
  
  for ( Itl.GoToBegin(), Ito.GoToBegin(), Itd.GoToBegin(); 
           !Itl.IsAtEnd(); ++Itl, ++Ito, ++Itd )
    {   
    PointDataType P;   
    P.Fill( 0 );
    if ( Ito.Get() != 0 )
      {
      P = Itd.Get()/Ito.Get();
      for ( unsigned int i = 0; i < P.Size(); i++ )
        {
        if ( vnl_math_isnan( P[i] ) || vnl_math_isinf( P[i] ) )
          {
          P[i] = 0;
          } 
        }         
      Itl.Set( P );
      }
    }
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::UpdatePointSet() 
{
  typename PointDataContainerType::ConstIterator ItIn;

  ItIn = this->m_InputPointData->Begin();
  while ( ItIn != this->m_InputPointData->End() )
    {
    PointType point;
    PointDataType dataOut;

    this->GetInput()->GetPoint( ItIn.Index(), &point );
    this->EvaluateAtPoint( point, dataOut );
    this->m_OutputPointData->InsertElement( ItIn.Index(), dataOut );
    ++ItIn;
    }
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::GenerateOutputImage() 
{
  ImageRegionIteratorWithIndex<ImageType> 
     It( this->GetOutput(), this->GetOutput()->GetBufferedRegion() );

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::IndexType idx = It.GetIndex();
    PointType point;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      point[i] = static_cast<RealType>( idx[i] ) 
               /( static_cast<RealType>( this->m_Size[i]-1 ) );
      }
    PointDataType data;
    this->Evaluate( point, data );
    It.Set( data );
    }
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::EvaluateAtPoint( PointType point, PointDataType &data ) 
{
  for ( unsigned int i = 0; i < ImageDimension; i++)
    {
    point[i] -= this->m_Origin[i];
    point[i] /= ( static_cast<RealType>( this->m_Size[i]*this->m_Spacing[i] ) );
    }  
  this->Evaluate( point, data );
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::EvaluateAtIndex( IndexType idx, PointDataType &data ) 
{
  PointType point;
  this->GetOutput()->TransformIndexToPhysicalPoint( idx, point );
  this->EvaluateAtPoint( point, data );
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::EvaluateAtContinuousIndex( ContinuousIndexType idx, PointDataType &data ) 
{
  PointType point;
  this->GetOutput()->TransformContinuousIndexToPhysicalPoint( idx, point );
  this->EvaluateAtPoint( point, data );
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::Evaluate( PointType params, PointDataType &data ) 
{
  vnl_vector<RealType> p( ImageDimension );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( params[i] >= 1.0 )
      {
      params[i] = 0.99999;
      }
    if ( params[i] < 0.0 )
      {
      params[i] = 0.0;
      }
    p[i] = static_cast<RealType>( params[i] ) 
         * static_cast<RealType>( this->m_CurrentNumberOfControlPoints[i]
                          - this->m_SplineOrder[i] );
    }

  typename RealImageType::RegionType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    size[i] = this->m_SplineOrder[i]+1;
    } 
  typename RealImageType::Pointer w;
  w = RealImageType::New();  
  w->SetRegions( size );
  w->Allocate(); 

  ImageRegionIteratorWithIndex<RealImageType> 
     Itw( w, w->GetLargestPossibleRegion() );

  PointDataType val;
  data.Fill( 0.0 );
  
  for ( Itw.GoToBegin(); !Itw.IsAtEnd(); ++Itw )
    {
    RealType B = 1.0;
    typename RealImageType::IndexType idx = Itw.GetIndex();
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {  
      RealType u = p[i] - static_cast<RealType>( static_cast<unsigned>(p[i]) + idx[i] ) 
                        + 0.5*static_cast<RealType>( this->m_SplineOrder[i] - 1 );
      B *= this->m_Kernel[i]->Evaluate( u );
      }  
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      idx[i] += static_cast<unsigned int>( p[i] );
      if ( this->m_CloseDimension[i] )
        {
        idx[i] %= this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[i];
        }  
      }      
    if ( this->m_PhiLattice->GetLargestPossibleRegion().IsInside( idx ) )
      {
      val = this->m_PhiLattice->GetPixel( idx );  
     val *= B;
      data += val;
      } 
    }
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::EvaluateGradientAtPoint( PointType point, GradientType &gradient ) 
{
  for ( unsigned int i = 0; i < ImageDimension; i++)
    {
    point[i] -= this->m_Origin[i];
    point[i] /= ( static_cast<RealType>( this->m_Size[i]*this->m_Spacing[i] ) );
    }  
  this->EvaluateGradient( point, gradient );
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::EvaluateGradientAtIndex( IndexType idx, GradientType &gradient ) 
{
  PointType point;
  this->GetOutput()->TransformIndexToPhysicalPoint( idx, point );
  this->EvaluateGradientAtPoint( point, gradient );
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::EvaluateGradientAtContinuousIndex( ContinuousIndexType idx, GradientType &gradient ) 
{
  PointType point;
  this->GetOutput()->TransformContinuousIndexToPhysicalPoint( idx, gradient );
  this->EvaluateGradientAtPoint( point, gradient );
}

template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::EvaluateGradient( PointType params, GradientType &gradient ) 
{
  vnl_vector<RealType> p( ImageDimension );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( params[i] >= 1.0 )
      {
      params[i] = 0.99999;
      }
    if ( params[i] < 0.0 )
      {
      params[i] = 0.0;
      }
    p[i] = static_cast<RealType>( params[i] ) 
         * static_cast<RealType>( this->m_CurrentNumberOfControlPoints[i]
                          - this->m_SplineOrder[i] );
    }

  typename RealImageType::RegionType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    size[i] = this->m_SplineOrder[i]+1;
    } 
  typename RealImageType::Pointer w;
  w = RealImageType::New();  
  w->SetRegions( size );
  w->Allocate(); 

  ImageRegionIteratorWithIndex<RealImageType> 
     Itw( w, w->GetLargestPossibleRegion() );

  PointDataType val;
  gradient.set_size( val.Size(), ImageDimension );
  gradient.fill( 0.0 );
  
  for ( unsigned int j = 0; j < gradient.cols(); j++ )
    {
    for ( Itw.GoToBegin(); !Itw.IsAtEnd(); ++Itw )
      {
      RealType B = 1.0;
      typename RealImageType::IndexType idx = Itw.GetIndex();
      for ( unsigned int k = 0; k < ImageDimension; k++ )
        {  
        RealType u = p[k] - static_cast<RealType>( static_cast<unsigned>(p[k]) + idx[k] ) 
                          + 0.5*static_cast<RealType>( this->m_SplineOrder[k] - 1 );
        if ( j == k )
          { 
          B *= this->m_Kernel[k]->EvaluateDerivative( u );
          }
        else
          {
          B *= this->m_Kernel[k]->Evaluate( u );
          }
        }  
      for ( unsigned int k = 0; k < ImageDimension; k++ )
        {
        idx[k] += static_cast<unsigned int>( p[k] );
        if ( this->m_CloseDimension[k] )
          {
          idx[k] %= this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[k];
          }  
        }      
      if ( this->m_PhiLattice->GetLargestPossibleRegion().IsInside( idx ) )
        {
        val = this->m_PhiLattice->GetPixel( idx );  
        val *= B;
        for ( unsigned int k = 0; k < val.Size(); k++ )
          {
          gradient( k, j ) += val[k];
          }
        } 
      }
   }  
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputPointSet, class TOutputImage>
void
BSplineDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {    
    this->m_Kernel[i]->Print( os, indent );
    }
  os << indent << "B-spline order: " 
     << this->m_SplineOrder << std::endl;
  os << indent << "Number Of control points: " 
     << this->m_NumberOfControlPoints << std::endl;
  os << indent << "Close dimension: " 
     << this->m_CloseDimension << std::endl;
  os << indent << "Number of levels " 
     << this->m_NumberOfLevels << std::endl;
}

} // end namespace itk

#endif
