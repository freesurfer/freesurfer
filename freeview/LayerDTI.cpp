/**
 * @brief Layer class for DTI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#include "LayerDTI.h"
#include "LayerPropertyDTI.h"
#include "MyUtils.h"
#include "FSVolume.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageMapToColors.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include <QDebug>

LayerDTI::LayerDTI( LayerMRI* ref, QObject* parent ) : LayerMRI( ref, parent ),
  m_vectorSource( NULL),
  m_eigenvalueSource(NULL)
{
  m_strTypeNames.push_back( "DTI" );
  if ( mProperty )
  {
    delete mProperty;
  }

  m_vectorData = vtkFloatArray::New();

  mProperty = new LayerPropertyDTI( this );
  LayerMRI::ConnectProperty();

  SetEditable( false );
}

LayerDTI::~LayerDTI()
{
  if ( m_vectorSource )
    delete m_vectorSource;

  if (m_eigenvalueSource)
    delete m_eigenvalueSource;

  m_vectorData->Delete();
}

bool LayerDTI::LoadDTIFromFile( )
{
  if ( !LayerMRI::LoadVolumeFromFile() )
  {
    return false;
  }

  if ( m_vectorSource )
  {
    delete m_vectorSource;
  }

  m_vectorSource = new FSVolume( m_volumeRef );
  m_vectorSource->SetResampleToRAS( m_bResampleToRAS );

  if ( !m_vectorSource->MRIRead(  m_sVectorFileName,
                                  m_sRegFilename.isEmpty() ? NULL : m_sRegFilename ) )
  {
    return false;
  }

  if (!m_sEigenvalueFileName.isEmpty())
  {
    m_eigenvalueSource = new FSVolume(m_volumeRef);
    if (!m_eigenvalueSource->MRIRead(m_sEigenvalueFileName, m_sRegFilename.isEmpty() ? NULL : m_sRegFilename ) )
    {
      return false;
    }
    else
    {
      double scale = GetProperty()->GetVectorScale() / m_eigenvalueSource->GetMaxValue()*2;
      GetProperty()->SetVectorScale(scale);
    }
  }

  if ( m_vectorSource->GetNumberOfFrames() < 3 )
  {
    std::cerr << "Vector data file is not valid.";
    return false;
  }

  InitializeDTIColorMap();

  return true;
}

void LayerDTI::InitializeDTIColorMap()
{
  vtkImageData* rasDTI = m_vectorSource->GetImageOutput();
  int* dim = rasDTI->GetDimensions();
  long long nSize = ((long long)dim[0])*dim[1]*dim[2];
  double v[4] = { 0, 0, 0, 1 };
  int c[3];
  vtkDataArray* vectors = rasDTI->GetPointData()->GetScalars();
  //  qDebug() << vectors->GetDataTypeAsString();
  m_vectorData->DeepCopy(vectors);
  vtkFloatArray* fas = vtkFloatArray::New();
  fas->DeepCopy( m_imageData->GetPointData()->GetScalars() );
#if VTK_MAJOR_VERSION > 5
  m_imageData->AllocateScalars(VTK_FLOAT, 2);
#else
  m_imageData->SetNumberOfScalarComponents( 2 );
  m_imageData->AllocateScalars();
#endif
  vtkMatrix4x4* rotation_mat = vtkMatrix4x4::New();
  rotation_mat->Identity();
  MATRIX* reg = m_vectorSource->GetRegMatrix();
  if ( reg )
  {
    for ( int i = 0; i < 3; i++ )
    {
      for ( int j = 0; j < 3; j++ )
      {
        rotation_mat->SetElement( j, i, *MATRIX_RELT( reg, i+1, j+1 ) );
      }
    }
  }
  float* ptr = (float*)m_imageData->GetScalarPointer();
  for ( long long i = 0; i < nSize; i++ )
  {
    vectors->GetTuple( i, v );
    rotation_mat->MultiplyPoint( v, v );
    vtkMath::Normalize( v );
    m_vectorData->SetTuple(i, v);
    double fa = fas->GetComponent( i, 0 );
    for ( int j = 0; j < 3; j++ )
    {
      c[j] = (int)(fabs(v[j]) * fa * 64);
      if ( c[j] > 63 )
      {
        c[j] = 63;
      }
    }
    float scalar = c[0]*64*64 + c[1]*64 + c[2];
    //    int x = i%dim[0];
    //    int y = (i/dim[0])%dim[1];
    //    int z = i/(dim[0]*dim[1]);
    //    m_imageData->SetScalarComponentFromFloat( x, y, z, 0, fa );
    //    m_imageData->SetScalarComponentFromFloat( x, y, z, 1, scalar );
    *(ptr+i*2) = fa;
    *(ptr+i*2+1) = scalar;
  }
  rotation_mat->Delete();
  fas->Delete();
}

void LayerDTI::UpdateColorMap()
{
  if ( GetProperty()->GetColorMap() == LayerPropertyMRI::DirectionCoded )
  {
    for ( int i = 0; i < 3; i++ )
    {
      mColorMap[i]->SetLookupTable( GetProperty()->GetDirectionCodedTable() );
      mColorMap[i]->SetActiveComponent( 1 );
    }
    if (m_imageData.GetPointer() && m_imageData->GetNumberOfScalarComponents() > 1)
    {
      int* dim = m_imageData->GetDimensions();
      long long nSize = ((long long)dim[0])*dim[1]*dim[2];
      float* ptr = (float*)m_imageData->GetScalarPointer();
      double v[4] = {0, 0, 0, 1};
      int c[3];
      double dMin = GetProperty()->GetMinGenericThreshold(),
          dMax = GetProperty()->GetMaxGenericThreshold();
      if (dMax > dMin)
      {
        for ( long long i = 0; i < nSize; i++ )
        {
          m_vectorData->GetTuple( i, v );
          float fa = *(ptr+i*2);
          for ( int j = 0; j < 3; j++ )
          {
            c[j] = (int)(fabs(v[j]) * (fa-dMin)/(dMax-dMin) * 64);
            if ( c[j] > 63 )
            {
              c[j] = 63;
            }
          }
          float scalar = c[0]*64*64 + c[1]*64 + c[2];
          *(ptr+i*2+1) = scalar;
        }
      }
      m_imageData->Modified();
    }

    emit ActorUpdated();
  }
  else
  {
    LayerMRI::UpdateColorMap();
  }
}

bool LayerDTI::GetVectorValue( double* pos, double* v_out )
{
  vtkImageData* rasDTI = m_vectorSource->GetImageOutput();
  if ( rasDTI == NULL )
  {
    return 0;
  }

  double* orig = rasDTI->GetOrigin();
  double* vsize = rasDTI->GetSpacing();
  int* ext = rasDTI->GetExtent();

  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos[i] - orig[i] ) / vsize[i] + 0.5 );
  }

  if ( n[0] < ext[0] || n[0] > ext[1] ||
       n[1] < ext[2] || n[1] > ext[3] ||
       n[2] < ext[4] || n[2] > ext[5] )
  {
    return false;
  }
  else
  {
    v_out[0] = rasDTI->GetScalarComponentAsDouble( n[0], n[1], n[2], 0 );
    v_out[1] = rasDTI->GetScalarComponentAsDouble( n[0], n[1], n[2], 1 );
    v_out[2] = rasDTI->GetScalarComponentAsDouble( n[0], n[1], n[2], 2 );
    return true;
  }
}

bool LayerDTI::DoRotate( std::vector<RotationElement>& rotations )
{
  m_bResampleToRAS = false;
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  m_vectorSource->SetResampleToRAS( m_bResampleToRAS );

  bool ret = LayerMRI::DoRotate( rotations ) && m_vectorSource->Rotate( rotations );

  InitializeDTIColorMap();
  return ret;
}

void LayerDTI::DoRestore()
{
  std::vector<RotationElement> rotations;

  LayerMRI::DoRestore();
  m_vectorSource->Rotate( rotations );
}

void LayerDTI::UpdateVectorActor( int nPlane )
{
  LayerMRI::UpdateVectorActor( nPlane, m_vectorSource->GetImageOutput(),
                               m_eigenvalueSource ? m_eigenvalueSource->GetImageOutput() : NULL );
}
