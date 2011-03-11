/**
 * @file  ControlPanel.h
 * @brief Main control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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
#ifndef ControlPanel_h
#define ControlPanel_h

#include <wx/wx.h>


class wxAuiNotebook;
class wxListBox;
class PanelVolume;
class PanelROI;
class PanelSurface;
class PanelWayPoints;
class PanelSceneSetting;

class ControlPanel : public wxPanel
{
public:
  ControlPanel(wxWindow* parent);
  virtual ~ControlPanel();

  void RaisePage( const wxString& title );

  wxString GetCurrentLayerCollectionName();

  void UpdateUI( bool bForce = false );

  wxAuiNotebook*      m_notebook;
  PanelVolume*        m_panelVolume;
  PanelROI*           m_panelROI;
  PanelSurface*       m_panelSurface;
  PanelWayPoints*     m_panelWayPoints;
  PanelSceneSetting*  m_panelSceneSetting;
  
private:

  DECLARE_EVENT_TABLE()
};

#endif

