/**
 * @file  LayerPLabel.h
 * @brief Layer class for P-Label volumes.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/01/05 18:02:57 $
 *    $Revision: 1.3 $
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
class vtkImageData;

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

  std::string GetLabelName( double* pos );
  
  virtual double GetVoxelValue( double* pos );
  
protected:
  bool DoRotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );
  void UpdateColorMap();

  FSVolume*       m_volumeTemp;
  wxArrayString   m_sFilenames;
  wxString        m_sFilenamePrefix;
  wxString        m_sLUT;
  vtkSmartPointer<vtkImageData> m_imageIndex;
};

#endif


