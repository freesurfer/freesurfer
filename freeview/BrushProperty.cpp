/**
 * @file  BrushProperty.cpp
 * @brief Class to hold brush properties for voxel editing
 *
 * Simple mix-in class for use with the Listener class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:45 $
 *    $Revision: 1.11 $
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


#include "BrushProperty.h"
#include "LayerVolumeBase.h"
#include <QSettings>

BrushProperty::BrushProperty (QObject* parent) : QObject(parent),
    m_nBrushSize( 1 ),
    m_nBrushTolerance( 0 ),
    m_bEnableDrawRange( false ),
    m_bEnableExcludeRange( false ),
    m_bDrawConnectedOnly( false ),
    m_layerRef( NULL )
{
    m_dDrawRange[0] = 0;
    m_dDrawRange[1] = 1000000;
    m_dExcludeRange[0] = 0;
    m_dExcludeRange[1] = 0;
    QSettings settings;
    m_nBrushSize = settings.value("/BrushProperty/Size", 1 ).toInt();
    // config->Read( _T("/BrushProperty/Tolerance"), &m_nBrushTolerance, 0L );
    m_nBrushTolerance = 0;
    m_bEnableDrawRange      = settings.value( "/BrushProperty/EnableDrawRange", false ).toBool();
    m_bEnableExcludeRange   = settings.value( "/BrushProperty/EnableExcludeRange", false ).toBool();
    m_bDrawConnectedOnly    = settings.value( "/BrushProperty/DrawConnected", false ).toBool();
    m_dDrawRange[0] = settings.value( "/BrushProperty/DrawRangeLow", 0 ).toDouble();
    m_dDrawRange[1] = settings.value( "/BrushProperty/DrawRangeHigh", 1000000 ).toDouble();
    m_dExcludeRange[0] = settings.value( "/BrushProperty/ExcludeRangeLow", 0 ).toDouble();
    m_dExcludeRange[1] = settings.value( "/BrushProperty/ExcludeRangeHigh", 0 ).toDouble();
}

BrushProperty::~BrushProperty()
{
    QSettings settings;
    settings.setValue( "/BrushProperty/Size", m_nBrushSize );
    settings.setValue( "/BrushProperty/Tolerance", m_nBrushTolerance );
    settings.setValue( "/BrushProperty/EnableDrawRange", m_bEnableDrawRange );
    settings.setValue( "/BrushProperty/EnableExcludeRange", m_bEnableExcludeRange );
    settings.setValue( "/BrushProperty/DrawConnected", m_bDrawConnectedOnly );
    settings.setValue( "/BrushProperty/DrawRangeLow", m_dDrawRange[0] );
    settings.setValue( "/BrushProperty/DrawRangeHigh", m_dDrawRange[1] );
    settings.setValue( "/BrushProperty/ExcludeRangeLow", m_dExcludeRange[0] );
    settings.setValue( "/BrushProperty/ExcludeRangeHigh", m_dExcludeRange[1] );
}

int BrushProperty::GetBrushSize()
{
  return m_nBrushSize;
}

void BrushProperty::SetBrushSize( int nSize )
{
  m_nBrushSize = nSize;
}

int BrushProperty::GetBrushTolerance()
{
  return m_nBrushTolerance;
}

void BrushProperty::SetBrushTolerance( int nTolerance )
{
  m_nBrushTolerance = nTolerance;
}

LayerVolumeBase* BrushProperty::GetReferenceLayer()
{
  return m_layerRef;
}

void BrushProperty::SetReferenceLayer( LayerVolumeBase* layer )
{
  m_layerRef = layer;
}

double* BrushProperty::GetDrawRange()
{
  return m_dDrawRange;
}

void BrushProperty::SetDrawRange( double* range )
{
  SetDrawRange( range[0], range[1] );
}

void BrushProperty::SetDrawRange( double low, double high )
{
  m_dDrawRange[0] = low;
  m_dDrawRange[1] = high;
}


bool BrushProperty::GetDrawRangeEnabled()
{
  return m_bEnableDrawRange;
}

void BrushProperty::SetDrawRangeEnabled( bool bEnable )
{
  m_bEnableDrawRange = bEnable;
}

double* BrushProperty::GetExcludeRange()
{
  return m_dExcludeRange;
}

void BrushProperty::SetExcludeRange( double* range )
{
  SetExcludeRange( range[0], range[1] );
}

void BrushProperty::SetExcludeRange( double low, double high )
{
  m_dExcludeRange[0] = low;
  m_dExcludeRange[1] = high;
}

bool BrushProperty::GetExcludeRangeEnabled()
{
  return m_bEnableExcludeRange;
}

void BrushProperty::SetExcludeRangeEnabled( bool bEnable )
{
  m_bEnableExcludeRange = bEnable;
}

bool BrushProperty::GetDrawConnectedOnly()
{
  return m_bDrawConnectedOnly;
}

void BrushProperty::SetDrawConnectedOnly( bool bEnable )
{
  m_bDrawConnectedOnly = bEnable;
}

void BrushProperty::OnLayerRemoved(Layer* layer)
{
    if (layer == m_layerRef)
        SetReferenceLayer(NULL);
}
