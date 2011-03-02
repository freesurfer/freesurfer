/**
 * @file  LayerDTI.cpp
 * @brief Layer class for DTI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.13 $
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
 */

#include <wx/wx.h>
#include "LayerDTI.h"
#include "LayerPropertiesDTI.h"
#include "MyUtils.h"
#include "FSVolume.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageMapToColors.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"

LayerDTI::LayerDTI( LayerMRI* ref ) : LayerMRI( ref ),
    m_vectorSource( NULL)
{
  m_strTypeNames.push_back( "DTI" );
  if ( mProperties )
    delete mProperties;

  mProperties = new LayerPropertiesDTI();
  mProperties->AddListener( this );

  SetEditable( false );
}

LayerDTI::~LayerDTI()
{
  if ( m_vectorSource )
    delete m_vectorSource;
}

bool LayerDTI::LoadDTIFromFile( wxWindow* wnd, wxCommandEvent& event )
{
  if ( !LayerMRI::LoadVolumeFromFile( wnd, event ) )
    return false;

  if ( m_vectorSource )
    delete m_vectorSource;

  m_vectorSource = new FSVolume( m_volumeRef );
  m_vectorSource->SetResampleToRAS( m_bResampleToRAS );
  event.SetInt( 25 );

  if ( !m_vectorSource->MRIRead(  m_sVectorFileName.c_str(),
                                  m_sRegFilename.size() > 0 ? m_sRegFilename.c_str() : NULL,
                                  wnd,
                                  event ) )
    return false;

  if ( m_vectorSource->GetNumberOfFrames() < 3 )
  {
    SetErrorString( "Vector data file is not valid." );
    return false;
  }

  event.SetInt( 50 );
  InitializeDTIColorMap( wnd, event );

  event.SetInt( 100 );
  wxPostEvent( wnd, event );

  return true;
}

void LayerDTI::InitializeDTIColorMap( wxWindow* wnd, wxCommandEvent& event )
{
  vtkImageData* rasDTI = m_vectorSource->GetImageOutput();
  int* dim = rasDTI->GetDimensions();
  int nSize = dim[0]*dim[1]*dim[2];
  double v[4] = { 0, 0, 0, 1 };
  int c[3];
  vtkDataArray* vectors = rasDTI->GetPointData()->GetScalars();
  vtkFloatArray* fas = vtkFloatArray::New();
  fas->DeepCopy( m_imageData->GetPointData()->GetScalars() );
  m_imageData->SetNumberOfScalarComponents( 2 );
  m_imageData->AllocateScalars();
  int nProgressStep = ( 99-event.GetInt() ) / 5;
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
  for ( int i = 0; i < nSize; i++ )
  {
    vectors->GetTuple( i, v );
    rotation_mat->MultiplyPoint( v, v );
    vtkMath::Normalize( v );
    double fa = fas->GetComponent( i, 0 );
    for ( int j = 0; j < 3; j++ )
    {
      c[j] = (int)(fabs(v[j]) * fa * 64);
      if ( c[j] > 63 )
        c[j] = 63;
    }
    float scalar = c[0]*64*64 + c[1]*64 + c[2];
    int x = i%dim[0];
    int y = (i/dim[0])%dim[1];
    int z = i/(dim[0]*dim[1]);
    m_imageData->SetScalarComponentFromFloat( x, y, z, 0, fa );
    m_imageData->SetScalarComponentFromFloat( x, y, z, 1, scalar );
    if ( wnd && nSize >= 5 && i%(nSize/5) == 0 )
    {
      event.SetInt( event.GetInt() + nProgressStep );
      wxPostEvent( wnd, event );
    }
  }
  rotation_mat->Delete();
  fas->Delete();
}

void LayerDTI::UpdateColorMap()
{
  if ( GetProperties()->GetColorMap() == LayerPropertiesMRI::DirectionCoded )
  {
    for ( int i = 0; i < 3; i++ )
    {
      mColorMap[i]->SetLookupTable( GetProperties()->GetDirectionCodedTable() );
      mColorMap[i]->SetActiveComponent( 1 );
    }
  }
  else
    LayerMRI::UpdateColorMap();
}

bool LayerDTI::GetVectorValue( double* pos, double* v_out )
{
  vtkImageData* rasDTI = m_vectorSource->GetImageOutput();
  if ( rasDTI == NULL )
    return 0;

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
    return false;
  else
  {
    v_out[0] = rasDTI->GetScalarComponentAsDouble( n[0], n[1], n[2], 0 );
    v_out[1] = rasDTI->GetScalarComponentAsDouble( n[0], n[1], n[2], 1 );
    v_out[2] = rasDTI->GetScalarComponentAsDouble( n[0], n[1], n[2], 2 );
    return true;
  }
}

bool LayerDTI::DoRotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event )
{
  m_bResampleToRAS = false;
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  m_vectorSource->SetResampleToRAS( m_bResampleToRAS );

  bool ret = LayerMRI::DoRotate( rotations, wnd, event ) && m_vectorSource->Rotate( rotations, wnd, event );

  InitializeDTIColorMap( wnd, event );
  return ret;
}

void LayerDTI::DoRestore()
{
  std::vector<RotationElement> rotations;
  
  LayerMRI::DoRestore();
  wxCommandEvent e;
  m_vectorSource->Rotate( rotations, NULL, e );
}

void LayerDTI::UpdateVectorActor( int nPlane )
{
  LayerMRI::UpdateVectorActor( nPlane, m_vectorSource->GetImageOutput() );
}
