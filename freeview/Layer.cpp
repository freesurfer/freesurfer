/**
 * @file  Layer.cpp
 * @brief Base Layer class. A layer is an independent data object with 2D and 3D graphical representations.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.19 $
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

#include "wx/wx.h"
#include "Layer.h"
#include "LayerProperties.h"
#include <math.h>

#define CLOSE_DISTANCE 1e-6

Layer::Layer() : Listener( "Layer" ), Broadcaster( "Layer" )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dSlicePosition[i] = 0;
    m_dWorldOrigin[i] = 0;
    m_dWorldVoxelSize[i] = 1;
    m_dWorldSize[i] = 0;
    m_dTranslate[i] = 0;
    m_dScale[i] = 1;
  }
  m_bLocked = false;
  mProperties = NULL;
}

Layer::~Layer()
{
  if ( mProperties )
    delete mProperties;
  
  SendBroadcast( "LayerObjectDeleted", this );
}

void Layer::SetName( const char* name )
{
  if ( m_strName != name )
  {
    m_strName = name;
    SendBroadcast( "LayerNameChanged", this );
  }
}

bool Layer::IsTypeOf( std::string tname )
{
  for ( size_t i = 0; i < m_strTypeNames.size(); i++ )
  {
    if ( m_strTypeNames[i] == tname )
      return true;
  }
  return false;
}

std::string Layer::GetEndType()
{
  if ( m_strTypeNames.size() > 0 )
    return m_strTypeNames[ m_strTypeNames.size()-1 ];
  else
    return "";
}

double* Layer::GetWorldOrigin()
{
  return m_dWorldOrigin;
}

void Layer::GetWorldOrigin( double* origin )
{
  for ( int i = 0; i < 3; i++ )
    origin[i] = m_dWorldOrigin[i];
}

void Layer::SetWorldOrigin( double* origin )
{
  for ( int i = 0; i < 3; i++ )
    m_dWorldOrigin[i] = origin[i];
}

double* Layer::GetWorldSize()
{
  return m_dWorldSize;
}

void Layer::GetWorldSize( double* size )
{
  for ( int i = 0; i < 3; i++ )
    size[i] = m_dWorldSize[i];
}

void Layer::SetWorldSize( double* size )
{
  for ( int i = 0; i < 3; i++ )
    m_dWorldSize[i] = size[i];
}

double* Layer::GetWorldVoxelSize()
{
  return m_dWorldVoxelSize;
}

void Layer::GetWorldVoxelSize( double* vs )
{
  for ( int i = 0; i < 3; i++ )
    vs[i] = m_dWorldVoxelSize[i];
}

void Layer::SetWorldVoxelSize( double* vs )
{
  for ( int i = 0; i < 3; i++ )
    m_dWorldVoxelSize[i] = vs[i];
}

double* Layer::GetSlicePosition()
{
  return m_dSlicePosition;
}

void Layer::GetSlicePosition( double* slicePos )
{
  for ( int i = 0; i < 3; i++ )
  {
    slicePos[i] = m_dSlicePosition[i];
  }
}

void Layer::SetSlicePosition( double* slicePos )
{
  this->BlockBroadcast( true );
  for ( int i = 0; i < 3; i++ )
  {
    SetSlicePosition( i, slicePos[i] );
  }
  this->BlockBroadcast( false );
//  this->SendBroadcast( "SlicePositionChanged", this );
//  this->SendBroadcast( "LayerActorUpdated", this );
}

void Layer::SetSlicePosition( int nPlane, double slicePos )
{
  wxASSERT( nPlane >= 0 && nPlane <= 2 );

  if ( fabs( slicePos - m_dSlicePosition[ nPlane ] ) > CLOSE_DISTANCE )
  {
    m_dSlicePosition[nPlane] = slicePos;
//    this->SendBroadcast( "SlicePositionChanged", this );
    OnSlicePositionChanged( nPlane );
//    this->SendBroadcast( "LayerActorUpdated", this );
  }
}

void Layer::RASToVoxel( const double* pos, int* n )
{
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( pos[i] - m_dWorldOrigin[i] ) / m_dWorldVoxelSize[i] + 0.5 );
  }
}

void Layer::VoxelToRAS( const int* n, double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = m_dWorldOrigin[i] + m_dWorldVoxelSize[i] * n[i];
  }
}

void Layer::Lock( bool bLock )
{
  m_bLocked = bLock;
  this->SendBroadcast( "LayerLockChanged", this );
}

void Layer::DoListenToMessage( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == _( "ShowInfoChanged" ) )
    this->SendBroadcast( "LayerShowInfoChanged", iData, this );
}

void Layer::GetBounds( double* bounds )
{
  double* origin = GetWorldOrigin();
  double* size = GetWorldSize();
  for ( int i = 0; i < 3; i++ )
  {
    bounds[i*2] = origin[i];
    bounds[i*2+1] = origin[i] + size[i];
  }
}

void Layer::GetDisplayBounds( double* bounds )
{
  this->GetBounds( bounds );
}

bool Layer::Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event )
{
  bool ret = DoRotate( rotations, wnd, event ); 
  if ( ret )
  {
    ResetTranslate();
    ResetScale();   
    this->SendBroadcast( "LayerTransformed", this, this );
  }
  return ret;
}

bool Layer::Translate( double x, double y, double z )
{
  double pos[3] = { x, y, z };
  return Translate( pos );
}
  
bool Layer::Translate( double* dPos )
{
  double offset[3];
  for ( int i = 0; i < 3; i++ )
    offset[i] = dPos[i] - m_dTranslate[i];
  
  DoTranslate( offset );
  
  for ( int i = 0; i < 3; i++ )
    m_dTranslate[i] = dPos[i];
  
  ResetScale();
  this->SendBroadcast( "LayerTransformed", this, this );
  
  return true;
}

void Layer::Scale( double* scale, int nSampleMethod )
{
  double rscale[3];
  for ( int i = 0; i < 3; i++ )
    rscale[i] = scale[i] / m_dScale[i];
  
  DoScale( rscale, nSampleMethod );
  
  for ( int i = 0; i < 3; i++ )
    m_dScale[i] = scale[i];
  
  ResetTranslate();
  this->SendBroadcast( "LayerTransformed", this, this );
}

// reset transformations
void Layer::Restore()
{
  DoRestore();
  
  for ( int i = 0; i < 3; i++ )
  {
    m_dTranslate[i] = 0;
    m_dScale[i] = 1;
  }
  
  this->SendBroadcast( "LayerTransformed", this, this );
}
