/**
 * @file  LayerPLabel.h
 * @brief Layer class for P-Label volumes.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
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


