/**
 * @file  SurfaceRegionGroups.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:42 $
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


