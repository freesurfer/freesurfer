/**
 * @file  Layer.cpp
 * @brief Base Layer class. A layer is an independent data object with 2D and 3D graphical representations.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/05/03 19:50:01 $
 *    $Revision: 1.25 $
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


#include "Layer.h"
#include "LayerProperty.h"
#include <math.h>

#define CLOSE_DISTANCE 1e-6

int Layer::m_nLastID = 0;

Layer::Layer( QObject* parent ) : QObject( parent )
{
  // assign unique ID
  m_nID = m_nLastID + 1;
  m_nLastID++;

  for ( int i = 0; i < 3; i++ )
  {
    m_dSlicePosition[i] = 0;
    m_dWorldOrigin[i] = 0;
    m_dWorldVoxelSize[i] = 1;
    m_dWorldSize[i] = 0;
    m_dTranslate[i] = 0;
    m_dScale[i] = 1;
    m_dRotate[i] = 0;
  }
  m_bLocked = false;
  mProperty = NULL;
  connect(this, SIGNAL(VisibilityChanged(bool)), this, SIGNAL(ActorUpdated()));
}

Layer::~Layer()
{
}

void Layer::SetName( const QString& name )
{
  if ( m_strName != name )
  {
    m_strName = name;
    emit NameChanged( name );
  }
}

bool Layer::IsTypeOf( const QString& tname )
{
  return m_strTypeNames.contains( tname );
}

QString Layer::GetEndType() const
{
  if ( m_strTypeNames.size() > 0 )
  {
    return m_strTypeNames[ m_strTypeNames.size()-1 ];
  }
  else
  {
    return "";
  }
}

double* Layer::GetWorldOrigin()
{
  return m_dWorldOrigin;
}

void Layer::GetWorldOrigin( double* origin )
{
  for ( int i = 0; i < 3; i++ )
  {
    origin[i] = m_dWorldOrigin[i];
  }
}

void Layer::SetWorldOrigin( double* origin )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldOrigin[i] = origin[i];
  }
}

double* Layer::GetWorldSize()
{
  return m_dWorldSize;
}

void Layer::GetWorldSize( double* size )
{
  for ( int i = 0; i < 3; i++ )
  {
    size[i] = m_dWorldSize[i];
  }
}

void Layer::SetWorldSize( double* size )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldSize[i] = size[i];
  }
}

double* Layer::GetWorldVoxelSize()
{
  return m_dWorldVoxelSize;
}

void Layer::GetWorldVoxelSize( double* vs )
{
  for ( int i = 0; i < 3; i++ )
  {
    vs[i] = m_dWorldVoxelSize[i];
  }
}

void Layer::SetWorldVoxelSize( double* vs )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldVoxelSize[i] = vs[i];
  }
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
  this->blockSignals( true );
  for ( int i = 0; i < 3; i++ )
  {
    SetSlicePosition( i, slicePos[i] );
  }
  this->blockSignals( false );;
}

void Layer::SetSlicePosition( int nPlane, double slicePos )
{
  if ( fabs( slicePos - m_dSlicePosition[ nPlane ] ) > CLOSE_DISTANCE )
  {
    m_dSlicePosition[nPlane] = slicePos;
    OnSlicePositionChanged( nPlane );
  }
}

void Layer::Lock( bool bLock )
{
  m_bLocked = bLock;
  emit Locked( bLock );
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


bool Layer::Rotate( std::vector<RotationElement>& rotations )
{
  bool ret = DoRotate( rotations );
  if ( ret )
  {
    ResetTranslate();
    ResetScale();
    emit Transformed();
  }
  return ret;
}

bool Layer::Transform(double *mat, int sample_method)
{
  DoTransform(mat, sample_method);
  emit Transformed();

  return true;
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
  {
    offset[i] = dPos[i] - m_dTranslate[i];
  }

  DoTranslate( offset );

  for ( int i = 0; i < 3; i++ )
  {
    m_dTranslate[i] = dPos[i];
  }

  ResetScale();
  emit Transformed();

  return true;
}

void Layer::Scale( double* scale, int nSampleMethod )
{
  double rscale[3];
  for ( int i = 0; i < 3; i++ )
  {
    rscale[i] = scale[i] / m_dScale[i];
  }

  DoScale( rscale, nSampleMethod );

  for ( int i = 0; i < 3; i++ )
  {
    m_dScale[i] = scale[i];
  }

  ResetTranslate();
  emit Transformed();
}

// reset transformations
void Layer::Restore()
{
  DoRestore();

  for ( int i = 0; i < 3; i++ )
  {
    m_dTranslate[i] = 0;
    m_dScale[i] = 1;
    m_dRotate[i] = 0;
  }

  emit Transformed();
}

void Layer::SetRotate(double *rotate, bool bAroundCenter)
{
  m_dRotate[0] = rotate[0];
  m_dRotate[1] = rotate[1];
  m_dRotate[2] = rotate[2];
  m_bRotateAroundCenter = bAroundCenter;

  UpdateTransform();
}

void Layer::SetTranslate(double *offset)
{
  m_dTranslate[0] = offset[0];
  m_dTranslate[1] = offset[1];
  m_dTranslate[2] = offset[2];
  UpdateTransform();
}

void Layer::SetScale(double *scale)
{
  m_dScale[0] = scale[0];
  m_dScale[1] = scale[1];
  m_dScale[2] = scale[2];
  UpdateTransform();
}

void Layer::UpdateTransform(int sample_method)
{
  DoTransform(sample_method);
  emit Transformed();
}
