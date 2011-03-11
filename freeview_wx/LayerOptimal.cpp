/**
 * @file  LayerOptimal.cpp
 * @brief Layer class for computed optimal volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
 *    $Revision: 1.1 $
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
#include "LayerOptimal.h"
#include "LayerPropertiesDTI.h"
#include "MyUtils.h"
#include "FSVolume.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageMapToColors.h"
#include "vtkLookupTable.h"
#include "vtkImageData.h"

LayerOptimal::LayerOptimal( LayerMRI* ref ) : LayerMRI( ref )
{
  m_strTypeNames.push_back( "Optimal" );

  SetEditable( false );
}

LayerOptimal::~LayerOptimal()
{}


bool LayerOptimal::Create( LayerMRI* layer_label, std::vector<LayerMRI*> layers )
{
  m_layerLabel = layer_label;
  m_layers = layers;

  m_layerLabel->AddListener( this );
  for ( size_t i = 0; i < layers.size(); i++ )
    m_layers[i]->AddListener( this );

  bool ret = LayerMRI::Create( layers[0], false );

  UpdateVolume();

  return ret;
}

bool LayerOptimal::UpdateVolume()
{
  if ( m_layers.size() < 2 )
    return false;

  // find labeled voxels in label volume
  vtkImageData* label_vol = m_layerLabel->GetImageData();
  int* dim = label_vol->GetDimensions();
// cout << m_layerLabel->GetImageData()->GetScalarType() << "  " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
  int nsize = dim[0] * dim[1] * dim[2];
  int fClass1 = 0, fClass2 = 0;
  int* vox1 = new int[nsize], *vox2 = new int[nsize];
  int nsize1 = 0, nsize2 = 0;

  switch ( label_vol->GetScalarType() )
  {
  case VTK_UNSIGNED_CHAR:
  {
    unsigned char* p = (unsigned char*)label_vol->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      if ( p[i] != 0 )
      {
        if ( fClass1 == 0 )
        {
          fClass1 = (int)p[i];
          vox1[nsize1++] = i;
        }
        else if ( (int)p[i] == fClass1 )
        {
          vox1[nsize1++] = i;
        }
        else if ( fClass2 == 0 )
        {
          fClass2 = (int)p[i];
          vox2[nsize2++] = i;
        }
        else if ( (int)p[i] == fClass2 )
        {
          vox2[nsize2++] = i;
        }
      }
    }
  }
  break;
  case VTK_INT:
  {
    int* p = (int*)label_vol->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      if ( p[i] != 0 )
      {
        if ( fClass1 == 0 )
        {
          fClass1 = (int)p[i];
          vox1[nsize1++] = i;
        }
        else if ( (int)p[i] == fClass1 )
        {
          vox1[nsize1++] = i;
        }
        else if ( fClass2 == 0 )
        {
          fClass2 = (int)p[i];
          vox2[nsize2++] = i;
        }
        else if ( (int)p[i] == fClass2 )
        {
          vox2[nsize2++] = i;
        }
      }
    }
  }
  break;
  case VTK_FLOAT:
  {
    float* p = (float*)label_vol->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      if ( p[i] != 0 )
      {
        if ( fClass1 == 0 )
        {
          fClass1 = (int)p[i];
          vox1[nsize1++] = i;
        }
        else if ( (int)p[i] == fClass1 )
        {
          vox1[nsize1++] = i;
        }
        else if ( fClass2 == 0 )
        {
          fClass2 = (int)p[i];
          vox2[nsize2++] = i;
        }
        else if ( (int)p[i] == fClass2 )
        {
          vox2[nsize2++] = i;
        }
      }
    }
  }
  break;
  case VTK_SHORT:
  {
    short* p = (short*)label_vol->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      if ( p[i] != 0 )
      {
        if ( fClass1 == 0 )
        {
          fClass1 = (int)p[i];
          vox1[nsize1++] = i;
        }
        else if ( (int)p[i] == fClass1 )
        {
          vox1[nsize1++] = i;
        }
        else if ( fClass2 == 0 )
        {
          fClass2 = (int)p[i];
          vox2[nsize2++] = i;
        }
        else if ( (int)p[i] == fClass2 )
        {
          vox2[nsize2++] = i;
        }
      }
    }
  }
  break;
  case VTK_LONG:
  {
    long* p = (long*)label_vol->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      if ( p[i] != 0 )
      {
        if ( fClass1 == 0 )
        {
          fClass1 = (int)p[i];
          vox1[nsize1++] = i;
        }
        else if ( (int)p[i] == fClass1 )
        {
          vox1[nsize1++] = i;
        }
        else if ( fClass2 == 0 )
        {
          fClass2 = (int)p[i];
          vox2[nsize2++] = i;
        }
        else if ( (int)p[i] == fClass2 )
        {
          vox2[nsize2++] = i;
        }
      }
    }
  }
  break;
  }

  if ( nsize1 < 2 || nsize2 < 2 )
    return false;

  std::vector<void*> vols;
  for ( size_t i = 0; i < m_layers.size(); i++ )
  {
    vols.push_back( m_layers[i]->GetImageData()->GetScalarPointer() );
    // cout << m_layers[i]->GetImageData()->GetScalarType() << endl;
  }

  void* out_p = GetImageData()->GetScalarPointer();
  switch ( m_layers[0]->GetImageData()->GetScalarType() )
  {
  case VTK_FLOAT:
    if ( !MyUtils::CalculateOptimalVolume( vox1, nsize1, vox2, nsize2, vols, (float*)out_p, nsize ) )
      return false;
    break;
  case VTK_INT:
    if ( !MyUtils::CalculateOptimalVolume( vox1, nsize1, vox2, nsize2, vols, (int*)out_p, nsize ) )
      return false;
    break;
  case VTK_SHORT:
    if ( !MyUtils::CalculateOptimalVolume( vox1, nsize1, vox2, nsize2, vols, (short*)out_p, nsize ) )
      return false;
    break;
  case VTK_UNSIGNED_CHAR:
    if ( !MyUtils::CalculateOptimalVolume( vox1, nsize1, vox2, nsize2, vols, (unsigned char*)out_p, nsize ) )
      return false;
    break;
  case VTK_LONG:
    if ( !MyUtils::CalculateOptimalVolume( vox1, nsize1, vox2, nsize2, vols, (long*)out_p, nsize ) )
      return false;
    break;
  }

  delete[] vox1;
  delete[] vox2;

  SetModified();

  return true;
}


void LayerOptimal::DoListenToMessage( std::string const iMessage, void* iData, void* sender )
{
  if ( iData == m_layerLabel && iMessage == "LayerEdited" )
  {
    this->UpdateVolume();
    // this->SendBroadcast( "LayerActorUpdated", this );
  }
  else if ( iMessage == "LayerObjectDeleted" )
  {
    for ( size_t i = 0; i < m_layers.size(); i++ )
    {
      if ( iData == m_layers[i] )
      {
        m_layers.erase( m_layers.begin() + i );
      }
    }
  }

  LayerMRI::DoListenToMessage( iMessage, iData, sender );
}
