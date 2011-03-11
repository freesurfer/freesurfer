/**
 * @file  PixelInfoPanel.h
 * @brief Main control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:40 $
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
#ifndef PixelInfoPanel_h
#define PixelInfoPanel_h

#include <wx/wx.h>
#include <wx/listctrl.h>
#include "Listener.h"
#include "Broadcaster.h"

#define PIXEL_INFO_LISTCTRL_ID        1000
#define PIXEL_INFO_INLINE_EDITOR_ID   1001

#define MENU_SHOW_VOXEL_COORDINATES   10000
#define MENU_SHOW_SHORT_NAME          10001
#define MENU_SHOW_SURFACE_COORDINATES 10002
#define MENU_SHOW_SURFACE_CURVATURE   10003

WX_DEFINE_ARRAY_LONG( long, ArrayLong );

class wxTextCtrl;
class wxContextMenuEvent;
class wxCommandEvent;
class LayerSurface;

class PixelInfoListCtrl : public wxListCtrl, public Listener
{
  friend class PixelInfoPanel;

public:
  PixelInfoListCtrl( wxWindow* parent, const wxString& name );
  virtual ~PixelInfoListCtrl();

  virtual wxString OnGetItemText(long item, long column) const;
  virtual wxListItemAttr *OnGetItemAttr(long item) const;

  void OnSize           ( wxSizeEvent& event );
  void OnContextMenu    ( wxContextMenuEvent& event );
  void OnDrag           ( wxListEvent& event );
  void OnButtonDown     ( wxMouseEvent& event );
  void OnFinishEditing  ( wxCommandEvent& event );
  void OnShowVoxelCoordinates   ( wxCommandEvent& event );
  void OnShowSurfaceCoordinates ( wxCommandEvent& event );
  void OnShowSurfaceCurvature   ( wxCommandEvent& event );
  void OnShowShortName          ( wxCommandEvent& event );

  void AddItem( const wxString& name, 
                const wxString& value = _(""),
                long userData = 0, 
                bool bUpdateCount = true );
  void AddSurfaceItem( LayerSurface* surf, double* pos, bool bUpdateCount = true );

  void SetRASPosition( double* pos );

  void HideEditor();

  bool GetShowVoxelCoordinates()
  {
    return m_bShowVoxelCoordinates;
  }

  void SetShowVoxelCoordinates( bool bShow );

  bool GetShowShortName()
  {
    return m_bShowShortName;
  }

  void SetShowShortName( bool bShow );
  
  bool GetShowSurfaceCoordinates()
  {
    return m_bShowSurfaceCoordinates;
  }
  
  void SetShowSurfaceCoordinates( bool bShow );
  
  bool GetShowSurfaceCurvature()
  {
    return m_bShowSurfaceCurvature;
  }
  void SetShowSurfaceCurvature( bool bShow );

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

protected:
  void OnInternalIdle();
  wxString GetShortName( wxString& name );

private:
  void Reset();
  void UpdateList();
  void DoUpdateList();

  void ResizeColumns();
  
  wxString AppendSpaceString( const wxString& input, int nDefaultLength = 0 );

  wxArrayString m_listName;
  wxArrayString m_listValue;
  ArrayLong     m_listPtr;

  wxListItemAttr  m_attr;

  double  m_dRASPosition[3];
  bool    m_bNeedUpdate;

  int     m_nRowEdited;

  bool    m_bShowVoxelCoordinates;
  bool    m_bShowSurfaceCoordinates;
  bool    m_bShowSurfaceCurvature;
  bool    m_bShowShortName;

  wxTextCtrl* m_textEditor;

  DECLARE_EVENT_TABLE()
};

class PixelInfoPanel : public wxPanel, public Broadcaster, public Listener
{
public:
  PixelInfoPanel( wxWindow* parent );
  virtual ~PixelInfoPanel();

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

  void ToggleShowVoxelCoordinates();
  void ToggleShortName();

private:
  PixelInfoListCtrl*  m_listCtrlCursor;
  PixelInfoListCtrl*  m_listCtrlMouse;

  DECLARE_EVENT_TABLE()
};

#endif

