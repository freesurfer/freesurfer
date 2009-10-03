/**
 * @file  LayerPLabel.cpp
 * @brief Layer class for p-label volumes.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/10/03 01:18:34 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <wx/wx.h>
#include "LayerPLabel.h"
#include "LayerPropertiesDTI.h"
#include "MyUtils.h"
#include "FSVolume.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageMapToColors.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "MainWindow.h"
#include "LUTDataHolder.h"
#include <vtkImageActor.h>
#include <vtkImageReslice.h>
#include <wx/filename.h>

LayerPLabel::LayerPLabel( LayerMRI* ref ) : LayerMRI( ref ),
    m_volumeTemp( NULL)
{
  m_strTypeNames.push_back( "PLabel" );

  SetEditable( false );
}

LayerPLabel::~LayerPLabel()
{
  if ( m_volumeTemp )
    delete m_volumeTemp;
}

bool LayerPLabel::LoadVolumeFiles( wxWindow* wnd, wxCommandEvent& event )
{
  if ( m_sFilenames.size() == 0 )
    return false;
  
  m_sFilename = m_sFilenames[0].c_str();
  if ( !LayerMRI::LoadVolumeFromFile( wnd, event ) )
    return false;
  
  if ( m_volumeTemp )
    delete m_volumeTemp;
  m_volumeTemp = NULL;
  
  LUTDataHolder* luts = MainWindow::GetMainWindowPointer()->GetLUTData();
  COLOR_TABLE* ct = luts->GetColorTable( m_sLUT );
  if ( !ct )
  {
    cerr << "Can not find look up table." << endl;
  }
  
  event.SetInt( 0 );

  m_imageData = vtkSmartPointer<vtkImageData>::New();
  m_imageData->SetScalarTypeToUnsignedChar();
  m_imageData->SetNumberOfScalarComponents(4);
  for ( size_t i = 0; i < m_sFilenames.size(); i++ )
  {
    wxString fn = wxFileName( m_sFilenames[i] ).GetName();
    if ( !m_sFilenamePrefix.IsEmpty() )
    {
      fn = fn.Mid( m_sFilenamePrefix.Len() );
      SetName( m_sFilenamePrefix.char_str() );
    }
    
    int r, g, b;
    if ( CTABrgbAtIndexi( ct, CTABentryNameToIndex( fn.char_str(), ct ), &r, &g, &b ) != 0 )
    {
      cerr << "Can not find index for color name " << fn.c_str() << endl;
      return false;
    }
    
    int color[4] = { r, g, b, 255 };
    m_volumeTemp = new FSVolume( m_volumeRef );
    if ( !m_volumeTemp->MRIRead( m_sFilenames[i].c_str(),
                                   m_sRegFilename.size() > 0 ? m_sRegFilename.c_str() : NULL,
                                   wnd,
                                   event ) )
    {
      cerr << "Can not load volume file " << m_sFilenames[i].c_str() << endl;
      return false;
    }
    
    vtkImageData* imageData = m_volumeTemp->GetImageOutput();
    if ( i == 0 ) // initialize m_imageData
    {
      m_imageData->SetDimensions( imageData->GetDimensions() );
      m_imageData->SetOrigin( imageData->GetOrigin() );
      m_imageData->SetSpacing( imageData->GetSpacing() );
      m_imageData->SetExtent( imageData->GetExtent() );
      m_imageData->AllocateScalars();
    }
    
    int* dim = m_imageData->GetDimensions();
    for ( int ni = 0; ni < dim[0]; ni++ )
    {
      for ( int nj = 0; nj < dim[1]; nj++ )
      {
        for ( int nk = 0; nk < dim[2]; nk++ )
        {
          float pvalue = imageData->GetScalarComponentAsFloat( ni, nj, nk, 0 );
          for ( int m = 0; m < 4; m++ )
          { 
            float fvalue = 0;
            if ( i != 0 )
              fvalue = m_imageData->GetScalarComponentAsFloat( ni, nj, nk, m );
            
            fvalue += color[m]*pvalue/255;
            if ( fvalue > 255 )
              fvalue = 255;
            m_imageData->SetScalarComponentFromFloat( ni, nj, nk, m, fvalue );
          }
        }
      }
    }
    
    delete m_volumeTemp;
    m_volumeTemp = NULL;
  }
  
  InitializeActors();
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetInput( mReslice[i]->GetOutput() );
    m_sliceActor3D[i]->SetInput( mReslice[i]->GetOutput() );
  }
  
  event.SetInt( 90 );

  event.SetInt( 100 );
  wxPostEvent( wnd, event );

  return true;
}

/*
void LayerPLabel::InitializeDTIColorMap( wxWindow* wnd, wxCommandEvent& event )
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
    if ( nSize >= 5 && i%(nSize/5) == 0 )
    {
      event.SetInt( event.GetInt() + nProgressStep );
      wxPostEvent( wnd, event );
    }
  }
  rotation_mat->Delete();
  fas->Delete();
}
*/
void LayerPLabel::UpdateColorMap()
{
  // over-ride parent class, do nothing
}


bool LayerPLabel::Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event )
{
  /* 
  m_bResampleToRAS = false;
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  m_vectorSource->SetResampleToRAS( m_bResampleToRAS );

  bool ret = LayerMRI::Rotate( rotations, wnd, event ) && m_vectorSource->Rotate( rotations, wnd, event );

  InitializeDTIColorMap( wnd, event );
  return ret;*/
  return true;
}

