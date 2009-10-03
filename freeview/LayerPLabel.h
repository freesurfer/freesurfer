/**
 * @file  LayerPLabel.h
 * @brief Layer class for P-Label volumes.
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

#ifndef LayerPLabel_h
#define LayerPLabel_h

#include "LayerMRI.h"
#include "vtkSmartPointer.h"
#include <string>

class wxWindow;
class wxCommandEvent;
class FSVolume;

class LayerPLabel : public LayerMRI
{
public:
  LayerPLabel( LayerMRI* ref );
  virtual ~LayerPLabel();

  bool LoadVolumeFiles( wxWindow* wnd, wxCommandEvent& event );

  void SetVolumeFileNames( const wxArrayString& filenames )
  {
    m_sFilenames = filenames;
  }

  void SetFileNamePrefix( const wxString& prefix )
  {
    m_sFilenamePrefix = prefix;
  }

  void SetLUT( const wxString& lut )
  {
    m_sLUT = lut;
  }
  
  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );

protected:
  void UpdateColorMap();

  FSVolume*       m_volumeTemp;
  wxArrayString   m_sFilenames;
  wxString        m_sFilenamePrefix;
  wxString        m_sLUT;
};

#endif


