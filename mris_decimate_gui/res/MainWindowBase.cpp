///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Dec 29 2008)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include "MainWindowBase.h"

///////////////////////////////////////////////////////////////////////////

MainWindowBase::MainWindowBase( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );
	
	MenuBar = new wxMenuBar( 0 );
	FileMenu = new wxMenu();
	wxMenuItem* FileOpen;
	FileOpen = new wxMenuItem( FileMenu, wxID_ANY, wxString( wxT("&Load Surface...") ) + wxT('\t') + wxT("CTRL+O"), wxEmptyString, wxITEM_NORMAL );
	FileMenu->Append( FileOpen );
	
	m_saveSurface = new wxMenuItem( FileMenu, ID_SAVE_SURFACE, wxString( wxT("&Save Surface") ) + wxT('\t') + wxT("CTRL+S"), wxEmptyString, wxITEM_NORMAL );
	FileMenu->Append( m_saveSurface );
	m_saveSurface->Enable( false );
	
	m_saveSurfaceAs = new wxMenuItem( FileMenu, ID_SAVE_SURFACE_AS, wxString( wxT("Save Surface &As") ) + wxT('\t') + wxT("CTRL+SHIFT+S"), wxEmptyString, wxITEM_NORMAL );
	FileMenu->Append( m_saveSurfaceAs );
	m_saveSurfaceAs->Enable( false );
	
	FileMenu->AppendSeparator();
	
	wxMenuItem* FileExit;
	FileExit = new wxMenuItem( FileMenu, wxID_ANY, wxString( wxT("E&xit") ) + wxT('\t') + wxT("ALT+F4"), wxEmptyString, wxITEM_NORMAL );
	FileMenu->Append( FileExit );
	
	MenuBar->Append( FileMenu, wxT("&File") );
	
	help = new wxMenu();
	wxMenuItem* About;
	About = new wxMenuItem( help, ID_ABOUT, wxString( wxT("About") ) , wxEmptyString, wxITEM_NORMAL );
	help->Append( About );
	
	MenuBar->Append( help, wxT("&Help") );
	
	this->SetMenuBar( MenuBar );
	
	
	// Connect Events
	this->Connect( wxEVT_CLOSE_WINDOW, wxCloseEventHandler( MainWindowBase::OnClose ) );
	this->Connect( FileOpen->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileOpen ) );
	this->Connect( m_saveSurface->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileSave ) );
	this->Connect( m_saveSurfaceAs->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileSaveAs ) );
	this->Connect( FileExit->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileExit ) );
	this->Connect( About->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnAbout ) );
}

MainWindowBase::~MainWindowBase()
{
	// Disconnect Events
	this->Disconnect( wxEVT_CLOSE_WINDOW, wxCloseEventHandler( MainWindowBase::OnClose ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileOpen ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileSave ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileSaveAs ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnFileExit ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindowBase::OnAbout ) );
}
