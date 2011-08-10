// Licensed under MIT license; see license.txt.

#include <sbl/gui/MainWindow.h>
#include <sbl/core/Init.h>
#include <sbl/core/Command.h>
#include <sbl/system/Timer.h> // for testProgressDialog
#include <sbl/gui/ImageSeqViewer.h>
#include <sbl/gui/ConfigViewer.h>
#include <sbl/gui/ConfigEditor.h>
namespace sbl {


//-------------------------------------------
// MAIN WINDOW CLASS
//-------------------------------------------


// current (single, global) instance of MainWindow
MainWindow *MainWindow::s_instance = NULL;


// basic constructor
MainWindow::MainWindow( const wxString &title, const wxPoint &pos, const wxSize &size )
                : wxFrame( (wxFrame *) NULL, -1, title, pos, size ) {

    // store pointer to self for external access
    s_instance = this;

    // no active timer
    m_timer = NULL;

    // this window will be the parent of any dialog windows
    setDialogParent( this );

    // add status bar
    CreateStatusBar();

    // add vertical sizer
    wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );

    // add tab control
    m_notebook = new wxNotebook( this, -1 );
    m_notebook->SetSizeHints( 700, 500 );
    sizer->Add( m_notebook, 1, wxEXPAND );

    // add command entry area
    m_commandInputText = new wxTextCtrl( this, -1 );
    m_commandInputText->SetMaxLength( 1000 );
    m_commandInputText->SetFocus();
    sizer->Add( m_commandInputText, 0, wxEXPAND );

    // add command output tab
    m_commandOutputList = new wxListBox( m_notebook, -1 );
    addTab( m_commandOutputList, "Command Output" );

    // add image sequence viewer
    ImageSeqViewer *imageSeqViewer = new ImageSeqViewer( m_notebook, -1, true );
    addTab( imageSeqViewer, "Image Viewer" );

    // finish layout
    SetSizerAndFit( sizer );

    // connect other event handlers
    m_commandInputText->Connect( wxEVT_KEY_DOWN, wxKeyEventHandler( MainWindow::onCommandInputKey ), NULL, this );
}


// basic destructor
MainWindow::~MainWindow() {
    if (m_timer)
        m_timer->Stop();
    runCleanUp();
    s_instance = NULL;
}


/// switch to a given tab
void MainWindow::switchTo( int tabIndex ) {
    m_notebook->SetSelection( tabIndex );
}


/// switch to a given tab
void MainWindow::switchTo( const wxWindow *window ) {
    for (unsigned int i = 0; i < m_notebook->GetPageCount(); i++) {
        if (m_notebook->GetPage( i ) == window) {
            m_notebook->SetSelection( i );
            break;
        }
    }
}


/// add a tab to the notebook
void MainWindow::addTab( wxWindow *window, const String &name ) {
    m_notebook->AddPage( window, name.c_str() );
}


/// display text in the command output viewer
void MainWindow::dispText( const char *text ) {
    if (m_commandOutputList->GetCount() > 5000)
        m_commandOutputList->Delete( 0 );
    m_commandOutputList->AppendString( text );
    m_commandOutputList->ScrollLines( 1 );
}


/// if command found for this code, execute it
void MainWindow::execKeyboardShortcut( int modifier, int code ) {
    KeyboardShortcut *ks = m_keyboardShortcuts.find( code );
    if (ks && ((ks->modifier & modifier) || ks->modifier == 0)) 
        execCommand( ks->command, false );
}


// handle keyboard input for command entry widget
void MainWindow::onCommandInputKey( wxKeyEvent &event ) {
    if (event.GetKeyCode() == WXK_RETURN || event.GetKeyCode() == WXK_NUMPAD_ENTER) {
        String text = m_commandInputText->GetValue().c_str();
        if (text.length())
            execCommand( text, true ); // fix(later): use a second thread?
        if (s_instance) // if exit command, don't clear the input text box because it no longer exists
            m_commandInputText->SetValue( "" );
    } else if (event.GetKeyCode() == WXK_TAB) {
        String prefix = m_commandInputText->GetValue().c_str();
        String completed = tabCompleteCommand( prefix );
        m_commandInputText->SetValue( completed.c_str() );
        m_commandInputText->SetInsertionPointEnd();
    } else if (event.GetKeyCode() == WXK_UP) {
        String history = nextHistoryCommand( -1 );
        m_commandInputText->SetValue( history.c_str() );
        m_commandInputText->SetInsertionPointEnd();    
    } else if (event.GetKeyCode() == WXK_DOWN) {
        String history = nextHistoryCommand( 1 );
        m_commandInputText->SetValue( history.c_str() );
        m_commandInputText->SetInsertionPointEnd();
    } else {
        execKeyboardShortcut( event.GetModifiers(), event.GetKeyCode() );
        event.Skip();
    }
}


// handle mouse wheel movement; redirect to image sequence viewer
void MainWindow::onMouseWheel( wxMouseEvent &event ) {
    if (ImageSeqViewer::instance()) 
        ImageSeqViewer::instance()->onMouseWheel( event );
}


// handle timer event to run register timed callbacks
// fix(later): pay attention to each TimerEvent's msecInterval
void MainWindow::onTimer( wxTimerEvent &event ) {
    for (int i = 0; i < m_timerEvents.count(); i++)
        m_timerEvents[ i ].callback();
}


/// register a function to be periodically called using wxWidgets timers
/// (currently restricted to msecInterval == 100)
void MainWindow::addTimerEvent( int msecInterval, void (*callback)() ) {
    assertAlways( msecInterval == 100 ); // fix(later): remove this constraint

    // if needed, add a timer
    if (m_timer == NULL) {
        m_timer = new wxTimer( this, -1 ); // fix(clean): do we need to deallocate this?
        m_timer->Connect( wxEVT_TIMER, wxTimerEventHandler( MainWindow::onTimer ), NULL, this );
        m_timer->Start( 100 );
    }

    // store timed event callback
    TimerEvent *e = new TimerEvent;
    e->msecInterval = msecInterval;
    e->callback = callback;
    m_timerEvents.append( e );
}


//-------------------------------------------
// DISPLAY WRAPPER FUNCTIONS
//-------------------------------------------


// display text in main window
void dispInfoMainWindow( const char *text ) {
    if (MainWindow::instance())
        MainWindow::instance()->dispText( text );
}


// display an error/warning dialog
void dispErrorMainWindow( const char *text ) {
    if (MainWindow::instance()) {
        MainWindow::instance()->dispText( text );
        wxMessageBox( text, MainApp::instance()->GetAppName().c_str(), wxOK | wxICON_EXCLAMATION, MainWindow::instance() );
    }
}


// display a message in the main window's status bar
void dispStatusMainWindow( const char *text ) {
    if (MainWindow::instance())
        MainWindow::instance()->SetStatusText( text );
}


// display/update a progress dialog
void dispProgressMainWindow( int index, int count ) {
    if (MainWindow::instance()) 
        ProgressDialog::dispProgress( index, count );
}


//-------------------------------------------
// USER INTERFACE COMMANDS
//-------------------------------------------


// switch to the specified a tab
void switchTab( Config &conf ) {
    int tabIndex = conf.readInt( "tabIndex" );
    MainWindow::instance()->switchTo( tabIndex );
}


// close the main window
void exitUserInterface( Config &conf ) {
    MainWindow::instance()->Close( true );
}


// test progress dialog
void testProgressDialog( Config &conf ) {
    bool mode = conf.readBool( "mode" );
    for (int i = 0; i < 1000; i++) {
        if (mode)
            progress( i );
        else
            progress( i, 1000 );
        delaySeconds( 0.01f );
        if (checkCommandEvents()) {
            disp( 1, "cancelled at %d", i );
            break;
        }
    }
    if (mode)
        progressDone();
}


//-------------------------------------------
// MAIN APP CLASS
//-------------------------------------------


// current (single, global) instance of MainApp
MainApp *MainApp::s_instance = NULL;


// allow wxWidgets to do processing
void runGUIEvents() {
    if (MainApp::instance()) 
        MainApp::instance()->Yield( true );
}


/// run all GUI SBL initialization 
bool MainApp::OnInit() {

    // store pointer to self for runGUIEvents
    s_instance = this;

    // create the main window
    MainWindow *mainWindow = new MainWindow( GetAppName().c_str(), wxPoint( 50, 50 ), wxSize( 900, 500 ) );

    // redirect displayed text to command output text area
    setDispCallback( dispInfoMainWindow );
    setErrorCallback( dispErrorMainWindow );
    setStatusCallback( dispStatusMainWindow );
    setProgressCallback( dispProgressMainWindow );

    // set other callbacks
    setCommandEventCallback( runGUIEvents );

    // register gui commands
    registerCommand( "x", exitUserInterface );
    registerCommand( "tab", switchTab );
#ifdef REGISTER_TEST_COMMANDS
    registerCommand( "testprogress", testProgressDialog );
#endif
    initConfigViewer();
    initImageSeqViewer();

    // add keyboard shortcuts
    mainWindow->addKeyboardShortcut( WXK_F1, "tab 0" );
    mainWindow->addKeyboardShortcut( WXK_F2, "tab 1" );
    mainWindow->addKeyboardShortcut( WXK_F3, "tab 2" );
    mainWindow->addKeyboardShortcut( WXK_F4, "tab 3" );
    mainWindow->addKeyboardShortcut( WXK_F5, "tab 4" );
    mainWindow->addKeyboardShortcut( WXK_F6, "tab 5" );
    mainWindow->addKeyboardShortcut( WXK_F7, "tab 6" );
    mainWindow->addKeyboardShortcut( WXK_F8, "tab 7" );
    mainWindow->addKeyboardShortcut( WXK_ESCAPE, "cancel" );

    // show main window
    mainWindow->Maximize( true );
    mainWindow->Show( TRUE );
    SetTopWindow( mainWindow );

    // init all SBL modules
    initModules();
    return TRUE;
} 


// override event filtering for keyboard and mouse wheel handling
// see here: http://wiki.wxwidgets.org/Catching_key_events_globally
// our approach: need to keep focus on the input text box; could not find good way to redirect keystrokes there (because native control?)
// better approach: do not allow other controls to accept focus?
int MainApp::FilterEvent( wxEvent &event ) { 
    if (event.GetEventType() == wxEVT_KEY_DOWN)
        MainWindow::instance()->setFocusCommandInput();
    if (event.GetEventType() == wxEVT_MOUSEWHEEL) {
        MainWindow::instance()->onMouseWheel( (wxMouseEvent &) event );
        return 1;
    }
    return -1;
}


} // end namespace sbl

