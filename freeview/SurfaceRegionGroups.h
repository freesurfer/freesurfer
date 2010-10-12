/**
 * @file  SurfaceRegionGroups.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/10/12 21:22:32 $
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

#ifndef SurfaceRegionGroups_h
#define SurfaceRegionGroups_h

#include "Broadcaster.h"
#include "Listener.h"
#include <wx/colour.h>

class LayerMRI;
class SurfaceRegion;

class SurfaceRegionGroups : public Broadcaster, public Listener
{
public:
  SurfaceRegionGroups( LayerMRI* owner );
  virtual ~SurfaceRegionGroups();

  wxColour GetGroupColor( int nGroup );
  void SetGroupColor( int nGroup, const wxColour& color );
  
  int GetGroupIdRange( SurfaceRegion* reg );

private:
  LayerMRI*     m_mri;
  std::vector<wxColour>    m_colors;
};

#endif


