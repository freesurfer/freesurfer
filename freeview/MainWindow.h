/**
 * @file  MainWindow.h
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/08/08 20:13:39 $
 *    $Revision: 1.13 $
 *
 * Copyright (C) 2002-2007,
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
 
#ifndef MainWindow_h
#define MainWindow_h

#include <wx/wx.h>
#include <wx/frame.h>
#include <wx/dynarray.h>
#include "Listener.h"
#include "Broadcaster.h"
#include "LayerCollectionManager.h"
#include <vector>

#define ID_WORKER_THREAD	wxID_HIGHEST + 1

class ControlPanel;
class wxPanel;
class wxSplitterWindow;
class wxFileHistory;
class wxSpinEvent;
class RenderView;
class RenderView2D;
class RenderView3D;
class PixelInfoPanel;
class WindowQuickReference;
class StatusBar;
class LayerCollectionManager;
class LUTDataHolder;
class wxToolBar;
class BrushProperty;
class ToolWindowEdit;

struct Settings2D
{
	bool	SyncZoomFactor;
};

struct SettingsScreenshot
{
	bool	HideCursor;
	bool	HideCoords;
	int		Magnification;
};

class MainWindow : public wxFrame, public Listener, public Broadcaster
{
public:
    MainWindow();
    virtual ~MainWindow();
	
	enum ViewLayout { VL_1X1 = 0, VL_2X2, VL_1N3 };
	enum MainView 	{ MV_Sagittal = 0, MV_Coronal, MV_Axial, MV_3D };
	
    // event handlers (these functions should _not_ be virtual)
    void OnFileOpen				( wxCommandEvent& event );
	void OnFileNew				( wxCommandEvent& event );
	void OnFileNewUpdateUI		( wxUpdateUIEvent& event );
	void OnFileSave				( wxCommandEvent& event );
	void OnFileSaveUpdateUI		( wxUpdateUIEvent& event );
	void OnFileSaveAs			( wxCommandEvent& event );
	void OnFileSaveAsUpdateUI	( wxUpdateUIEvent& event );
    void OnFileExit				( wxCommandEvent& event );
	void OnFileRecent			( wxCommandEvent& event );
	void OnFileNewROI			( wxCommandEvent& event );
	void OnFileNewROIUpdateUI	( wxUpdateUIEvent& event );
	void OnFileLoadROI			( wxCommandEvent& event );
	void OnFileLoadROIUpdateUI	( wxUpdateUIEvent& event );
	void OnFileSaveROI			( wxCommandEvent& event );
	void OnFileSaveROIUpdateUI	( wxUpdateUIEvent& event );
	void OnFileSaveROIAs		( wxCommandEvent& event );
	void OnFileSaveROIAsUpdateUI( wxUpdateUIEvent& event );
	void OnFileLoadDTI			( wxCommandEvent& event );
	void OnFileSaveScreenshot	( wxCommandEvent& event );
	void OnFileSaveScreenshotUpdateUI( wxUpdateUIEvent& event );
	void OnFileLoadSurface		( wxCommandEvent& event );
	
	void OnFileNewWayPoints			( wxCommandEvent& event );
	void OnFileNewWayPointsUpdateUI	( wxUpdateUIEvent& event );
	void OnFileLoadWayPoints		( wxCommandEvent& event );
	void OnFileLoadWayPointsUpdateUI	( wxUpdateUIEvent& event );
	void OnFileSaveWayPoints			( wxCommandEvent& event );
	void OnFileSaveWayPointsUpdateUI	( wxUpdateUIEvent& event );
	void OnFileSaveWayPointsAs			( wxCommandEvent& event );
	void OnFileSaveWayPointsAsUpdateUI	( wxUpdateUIEvent& event );
	
	void OnViewLayout1X1		( wxCommandEvent& event );
	void OnViewLayout1X1UpdateUI( wxUpdateUIEvent& event );
	void OnViewLayout2X2		( wxCommandEvent& event );
	void OnViewLayout2X2UpdateUI( wxUpdateUIEvent& event );
	void OnViewLayout1N3		( wxCommandEvent& event );
	void OnViewLayout1N3UpdateUI( wxUpdateUIEvent& event );
	void OnViewSagittal			( wxCommandEvent& event );
	void OnViewSagittalUpdateUI	( wxUpdateUIEvent& event );
	void OnViewCoronal			( wxCommandEvent& event );
	void OnViewCoronalUpdateUI	( wxUpdateUIEvent& event );
	void OnViewAxial			( wxCommandEvent& event );
	void OnViewAxialUpdateUI	( wxUpdateUIEvent& event );
	void OnView3D				( wxCommandEvent& event );
	void OnView3DUpdateUI		( wxUpdateUIEvent& event );
	void OnViewReset			( wxCommandEvent& event );
	void OnViewResetUpdateUI	( wxUpdateUIEvent& event );
	
	void OnViewCycleLayer					( wxCommandEvent& event );
	void OnViewCycleLayerUpdateUI			( wxUpdateUIEvent& event );
	void OnViewToggleVoxelCoordinates		( wxCommandEvent& event );
	void OnViewToggleVolumeVisibility		( wxCommandEvent& event );
	void OnViewToggleVolumeVisibilityUpdateUI( wxUpdateUIEvent& event );
	void OnViewToggleROIVisibility			( wxCommandEvent& event );
	void OnViewToggleROIVisibilityUpdateUI	( wxUpdateUIEvent& event );
	void OnViewToggleSurfaceVisibility		( wxCommandEvent& event );
	void OnViewToggleSurfaceVisibilityUpdateUI( wxUpdateUIEvent& event );
	void OnViewToggleWayPointsVisibility	( wxCommandEvent& event );
	void OnViewToggleWayPointsVisibilityUpdateUI( wxUpdateUIEvent& event );
	void OnViewToggleCursorVisibility			( wxCommandEvent& event );
	void OnViewToggleCursorVisibilityUpdateUI	( wxUpdateUIEvent& event );
	
	void OnViewSurfaceMain			( wxCommandEvent& event );
	void OnViewSurfaceMainUpdateUI	( wxUpdateUIEvent& event );
	void OnViewSurfaceInflated		( wxCommandEvent& event );
	void OnViewSurfaceInflatedUpdateUI( wxUpdateUIEvent& event );
	void OnViewSurfaceWhite			( wxCommandEvent& event );
	void OnViewSurfaceWhiteUpdateUI	( wxUpdateUIEvent& event );
	void OnViewSurfacePial			( wxCommandEvent& event );
	void OnViewSurfacePialUpdateUI	( wxUpdateUIEvent& event );
	void OnViewSurfaceOriginal		( wxCommandEvent& event );
	void OnViewSurfaceOriginalUpdateUI( wxUpdateUIEvent& event );
	
	void OnModeNavigate			( wxCommandEvent& event );
	void OnModeNavigateUpdateUI	( wxUpdateUIEvent& event);	
	void OnModeVoxelEdit		( wxCommandEvent& event );
	void OnModeVoxelEditUpdateUI( wxUpdateUIEvent& event);
	void OnModeROIEdit			( wxCommandEvent& event );
	void OnModeROIEditUpdateUI	( wxUpdateUIEvent& event);
	void OnModeWayPointsEdit			( wxCommandEvent& event );
	void OnModeWayPointsEditUpdateUI	( wxUpdateUIEvent& event);
	
	void OnActionVoxelFreehand	( wxCommandEvent& event );
	void OnActionVoxelFreehandUpdateUI( wxUpdateUIEvent& event );
	void OnActionVoxelFill		( wxCommandEvent& event );
	void OnActionVoxelFillUpdateUI( wxUpdateUIEvent& event );
	void OnActionVoxelPolyline	( wxCommandEvent& event );
	void OnActionVoxelPolylineUpdateUI( wxUpdateUIEvent& event );
	void OnActionROIFreehand	( wxCommandEvent& event );
	void OnActionROIFreehandUpdateUI( wxUpdateUIEvent& event );
	void OnActionROIFill		( wxCommandEvent& event );
	void OnActionROIFillUpdateUI( wxUpdateUIEvent& event );
	void OnActionROIPolyline	( wxCommandEvent& event );
	void OnActionROIPolylineUpdateUI( wxUpdateUIEvent& event );
	
	void OnEditCopy( wxCommandEvent& event );
	void OnEditCopyUpdateUI( wxUpdateUIEvent& event );
	void OnEditPaste( wxCommandEvent& event );
	void OnEditPasteUpdateUI( wxUpdateUIEvent& event );
	void OnEditUndo( wxCommandEvent& event );
	void OnEditUndoUpdateUI( wxUpdateUIEvent& event );
	void OnEditRedo( wxCommandEvent& event );
	void OnEditRedoUpdateUI( wxUpdateUIEvent& event );
	void OnEditPreferences( wxCommandEvent& event );
	
	void OnHelpQuickReference( wxCommandEvent& event );
	void OnHelpAbout( wxCommandEvent& event );
		
	void OnWorkerThreadResponse( wxCommandEvent& event );
	
	void OnSpinBrushSize( wxSpinEvent& event );
	void OnSpinBrushTolerance( wxSpinEvent& event );
	void OnCheckBrushTemplate( wxCommandEvent& event );
	void OnChoiceBrushTemplate( wxCommandEvent& event );
	
	void OnActivate	( wxActivateEvent& event );
	void OnIconize	( wxIconizeEvent& event );
	void OnClose	( wxCloseEvent &event );
	void OnKeyDown	( wxKeyEvent &event );
	
	void LoadVolume();
	void NewVolume();
	void SaveVolume();
	void SaveVolumeAs();
	
	void LoadROI();
	void NewROI();
	void SaveROI();
	void SaveROIAs();
	
	void LoadWayPoints();
	void NewWayPoints();
	void SaveWayPoints();
	void SaveWayPointsAs();
	
	void LoadSurface();

	void LoadVolumeFile( const wxString& fn, bool bResample = true, int nColorMap = 0 );
	void LoadDTIFile( const wxString& fn_fa, const wxString& fn_vector, bool Resample = true );
	void LoadROIFile( const wxString& fn );	
	void LoadSurfaceFile( const wxString& fn );
	void LoadWayPointsFile( const wxString& fn );
	
	bool IsSaving()
		{ return m_bSaving; }
	
	bool IsLoading()
		{ return m_bLoading; }
	
	LayerCollection* GetLayerCollection( std::string strType );
	LayerCollectionManager* GetLayerCollectionManager();
	
	LUTDataHolder* GetLUTData()
		{ return m_luts; }
	
	static MainWindow* GetMainWindowPointer();

	void NeedRedraw( int nCount = 2 );	
		
	void AddScript( const wxArrayString& script );
	void RunScript(); 
	
	int GetActiveViewId();
	
	RenderView* GetActiveView();
	
	RenderView* GetRenderView( int nIndex )
		{ return m_viewRender[nIndex]; }
	
	BrushProperty* GetBrushProperty()
		{ return m_propertyBrush; }
	
	Settings2D Get2DSettings()
		{ return m_settings2D; }
	
	SettingsScreenshot GetScreenshotSettings()
		{ return m_settingsScreenshot; }
	
	void SetAction( int nAction );
	
	void SetMode( int nMode );
	
protected:
	void OnInternalIdle(); 
	virtual void DoListenToMessage ( std::string const iMsg, void* const iData );
	
private:
	void SetViewLayout( int nLayout );
	void SetMainView( int nView );
	void UpdateToolbars();
	void DoUpdateToolbars();
	
	ControlPanel*		m_controlPanel;
	PixelInfoPanel*		m_pixelInfoPanel;
	wxSplitterWindow*	m_splitterMain;
	wxSplitterWindow*	m_splitterSub;
	wxPanel*			m_renderViewHolder;	
	WindowQuickReference*	m_wndQuickReference;
	StatusBar*			m_statusBar;
	wxToolBar*			m_toolbarVoxelEdit;
	wxToolBar*			m_toolbarROIEdit;
	wxToolBar*			m_toolbarBrush;
	wxPanel*			m_panelToolbarHolder;
	ToolWindowEdit*		m_toolWindowEdit;
	
	RenderView2D*		m_viewAxial;
	RenderView2D*		m_viewSagittal;
	RenderView2D*		m_viewCoronal;
	RenderView3D*		m_view3D;
	RenderView*			m_viewRender[4];
	int					m_nViewLayout;
	int					m_nMainView;
	int					m_nPrevActiveViewId;
	
	LayerCollectionManager*	m_layerCollectionManager;
	
	Settings2D			m_settings2D;
	SettingsScreenshot	m_settingsScreenshot;
	
	wxString			m_strLastDir;
	wxFileHistory*		m_fileHistory;
	int					m_nMaxRecentFiles;
	bool				m_bResampleToRAS;
	int					m_nScreenshotFilterIndex;
	
	LUTDataHolder*		m_luts;
	
	int				m_nRedrawCount;
	bool			m_bToUpdateToolbars;
	
	bool			m_bSaving;
	bool			m_bLoading;
    
	std::vector<wxArrayString>	m_scripts;
	
	BrushProperty*	m_propertyBrush;
	
    DECLARE_CLASS( MainWindow )
    // any class wishing to process wxWindows events must use this macro
    DECLARE_EVENT_TABLE()

};

#endif 


