/**
 * @file  VolumeFilter.cpp
 * @brief Base VolumeFilter class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/07/20 19:34:09 $
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

#include "VolumeFilter.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>

VolumeFilter::VolumeFilter( LayerMRI* input, LayerMRI* output ) : 
    Listener( "VolumeFilter" ), 
    Broadcaster( "VolumeFilter" )
{
  SetVolumes( input, output );
}

VolumeFilter::~VolumeFilter()
{
}

void VolumeFilter::SetVolumes( LayerMRI* input, LayerMRI* output )
{
  m_MRIInput = input;
  m_MRIOutput = output;
  m_MRIInput->AddListener( this );
  m_MRIOutput->AddListener( this );
}

void VolumeFilter::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "LayerObjectDeleted" )
  {
    if ( sender == m_MRIInput || iData == m_MRIInput )
      m_MRIInput = NULL;
    else if ( sender == m_MRIOutput || iData == m_MRIOutput )
      m_MRIOutput = NULL;
  }
}

bool VolumeFilter::ReadyToUpdate()
{
  return ( m_MRIInput && m_MRIOutput );
}

