/**
 * @file  MainWindow.h
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/03/27 20:39:00 $
 *    $Revision: 1.2 $
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

#include <wx/frame.h>
#include "Listener.h"
#include "Broadcaster.h"
#include "LayerCollectionManager.h"

#define ID_WORKER_THREAD	wxID_HIGHEST + 1

class ControlPanel;
class wxPanel;
class wxSplitterWindow;
class wxFileHistory;
class RenderView2D;
class RenderView3D;
class PixelInfoPanel;
class WindowQuickReference;
class StatusBar;
class LayerCollectionManager;
class LUTDataHolder;

class MainWindow : public wxFrame, public Listener, public Broadcaster
{
public:
    MainWindow();
    virtual ~MainWindow();
	
	enum ViewLayout { VL_1X1 = 0, VL_2X2, VL_1N3 };
	enum MainView 	{ MV_Sagittal = 0, MV_Coronal, MV_Axial, MV_3D };
	
    // event handlers (these functions should _not_ be virtual)
    void OnFileOpen( wxCommandEvent& event );
	void OnFileNew( wxCommandEvent& event );
	void OnFileNewUpdateUI( wxUpdateUIEvent& event );
	void OnFileSave( wxCommandEvent& event );
	void OnFileSaveUpdateUI( wxUpdateUIEvent& event );
	void OnFileSaveAs( wxCommandEvent& event );
	void OnFileSaveAsUpdateUI( wxUpdateUIEvent& event );
    void OnFileExit( wxCommandEvent& event );
	void OnFileRecent( wxCommandEvent& event );
	void OnFileNewROI( wxCommandEvent& event );
	void OnFileNewROIUpdateUI( wxUpdateUIEvent& event );
	void OnFileLoadROI( wxCommandEvent& event );
	void OnFileLoadROIUpdateUI( wxUpdateUIEvent& event );
	void OnFileSaveROI( wxCommandEvent& event );
	void OnFileSaveROIUpdateUI( wxUpdateUIEvent& event );
	void OnFileSaveROIAs( wxCommandEvent& event );
	void OnFileSaveROIAsUpdateUI( wxUpdateUIEvent& event );
	void OnFileLoadDTI( wxCommandEvent& event );
	
	void OnViewLayout1X1( wxCommandEvent& event );
	void OnViewLayout1X1UpdateUI( wxUpdateUIEvent& event );
	void OnViewLayout2X2( wxCommandEvent& event );
	void OnViewLayout2X2UpdateUI( wxUpdateUIEvent& event );
	void OnViewLayout1N3( wxCommandEvent& event );
	void OnViewLayout1N3UpdateUI( wxUpdateUIEvent& event );
	void OnViewSagittal( wxCommandEvent& event );
	void OnViewSagittalUpdateUI( wxUpdateUIEvent& event );
	void OnViewCoronal( wxCommandEvent& event );
	void OnViewCoronalUpdateUI( wxUpdateUIEvent& event );
	void OnViewAxial( wxCommandEvent& event );
	void OnViewAxialUpdateUI( wxUpdateUIEvent& event );
	void OnView3D( wxCommandEvent& event );
	void OnView3DUpdateUI( wxUpdateUIEvent& event );
	
	void OnActionNavigate( wxCommandEvent& event );
	void OnActionNavigateUpdateUI( wxUpdateUIEvent& event);
	void OnActionVoxelFreehand( wxCommandEvent& event );
	void OnActionVoxelFreehandUpdateUI( wxUpdateUIEvent& event );
	void OnActionVoxelFill( wxCommandEvent& event );
	void OnActionVoxelFillUpdateUI( wxUpdateUIEvent& event );
	void OnActionVoxelPolyline( wxCommandEvent& event );
	void OnActionVoxelPolylineUpdateUI( wxUpdateUIEvent& event );
	void OnActionROIFreehand( wxCommandEvent& event );
	void OnActionROIFreehandUpdateUI( wxUpdateUIEvent& event );
	void OnActionROIFill( wxCommandEvent& event );
	void OnActionROIFillUpdateUI( wxUpdateUIEvent& event );
	void OnActionROIPolyline( wxCommandEvent& event );
	void OnActionROIPolylineUpdateUI( wxUpdateUIEvent& event );
	
	void OnEditUndo( wxCommandEvent& event );
	void OnEditUndoUpdateUI( wxUpdateUIEvent& event );
	void OnEditRedo( wxCommandEvent& event );
	void OnEditRedoUpdateUI( wxUpdateUIEvent& event );
	void OnEditPreferences( wxCommandEvent& event );
	
	void OnHelpQuickReference( wxCommandEvent& event );
		
	void OnWorkerThreadResponse( wxCommandEvent& event );
	
	void OnActivate( wxActivateEvent& event );
	void OnClose( wxCloseEvent &event );
	void OnKeyDown( wxKeyEvent &event );
	
	void LoadVolume();
	void NewVolume();
	void SaveVolume();
	void SaveVolumeAs();
	
	void LoadROI();
	void NewROI();
	void SaveROI();
	void SaveROIAs();

	void LoadVolumeFile( const wxString& fn, bool bResample = true, int nColorMap = 0 );
	void LoadDTIFile( const wxString& fn_fa, const wxString& fn_vector, bool Resample = true );
	void LoadROIFile( const wxString& fn );	
	
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
	
protected:
	void OnInternalIdle(); 
	virtual void DoListenToMessage ( std::string const iMsg, void* const iData );
	
private:
	void SetViewLayout( int nLayout );
	void SetMainView( int nView );
	
	ControlPanel*		m_controlPanel;
	PixelInfoPanel*		m_pixelInfoPanel;
	wxSplitterWindow*	m_splitterMain;
	wxSplitterWindow*	m_splitterSub;
	wxPanel*			m_renderViewHolder;	
	WindowQuickReference*	m_wndQuickReference;
	StatusBar*			m_statusBar;
	
	RenderView2D*		m_viewAxial;
	RenderView2D*		m_viewSagittal;
	RenderView2D*		m_viewCoronal;
	RenderView3D*		m_view3D;
	int					m_nViewLayout;
	int					m_nMainView;
	
	LayerCollectionManager*	m_layerCollectionManager;
	
	wxString			m_strLastDir;
	wxFileHistory*		m_fileHistory;
	int					m_nMaxRecentFiles;
	bool				m_bResampleToRAS;
	
	LUTDataHolder*		m_luts;
	
	int				m_nRedrawCount;
	
	bool			m_bSaving;
	bool			m_bLoading;
    
    DECLARE_CLASS( MainWindow )
    // any class wishing to process wxWindows events must use this macro
    DECLARE_EVENT_TABLE()

};

#endif 


