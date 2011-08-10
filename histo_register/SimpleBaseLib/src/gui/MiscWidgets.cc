// Licensed under MIT license; see license.txt.

#include <sbl/gui/MiscWidgets.h>
#include <sbl/core/Command.h>
namespace sbl {


//-------------------------------------------
// PROGRESS DIALOG CLASS
//-------------------------------------------


// basic constructor
ProgressDialog::ProgressDialog( wxWindow *parent, wxWindowID id ) : wxDialog( parent, id, wxString( "Progress" ), wxDefaultPosition, wxDefaultSize, wxCAPTION ) {

    // populate dialog contents
    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );
    m_progressLabel = new wxStaticText( this, -1, "", wxDefaultPosition, wxDefaultSize, wxALIGN_CENTER );
    m_progressLabel->SetMinSize( wxSize( 120, 18 ) );
    topSizer->Add( m_progressLabel, 1, wxEXPAND | wxALL, 5 );
    wxButton *cancelButton = new wxButton( this, -1, "Cancel" );
    topSizer->Add( cancelButton, 1, wxEXPAND | wxBOTTOM | wxRIGHT | wxLEFT, 5 );
    SetSizerAndFit( topSizer );
    SetDoubleBuffered( true );

    // connect event handlers
    cancelButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ProgressDialog::onCancelButton ), NULL, this );
}


/// update the displayed progress;
/// (assumes index in [0, count - 1]);
/// closes progress bar when index == count - 1; 
/// if count == -1, assumes unknown number if items
void ProgressDialog::updateProgress( int index, int count ) {
    if (count > 0) {
        m_progressLabel->SetLabel( sprintF( "%d of %d", index + 1, count ).c_str() );
    } else {
        m_progressLabel->SetLabel( sprintF( "%d", index ).c_str() );
    }
}


/// called when user presses cancel button
void ProgressDialog::onCancelButton( wxCommandEvent &event ) {
    MakeModal( false );
    Hide();
    wxEndBusyCursor();
    setCancelCommand( true );
}


// display or update a progress dialog
void ProgressDialog::dispProgress( int index, int count ) {
    static ProgressDialog *s_progressDialog = NULL;
    if (index + 1 < count || count <= 0) {
        if (s_progressDialog == NULL) {
            s_progressDialog = new ProgressDialog( dialogParent(), -1 );
        }
        if (s_progressDialog->IsVisible() == false) {
            s_progressDialog->Show();
            s_progressDialog->MakeModal( true ); // use this rather than ShowModal because want to keep running our code
            wxBeginBusyCursor();
        }
        s_progressDialog->updateProgress( index, count );
    } else {
        if (s_progressDialog) {
            s_progressDialog->MakeModal( false );
            s_progressDialog->Hide();
            delete s_progressDialog;
            s_progressDialog = NULL;
        }
        if (wxIsBusy())
            wxEndBusyCursor();
    }
}


//-------------------------------------------
// TOOL BAR CLASS
//-------------------------------------------


// basic constructor
ToolBar::ToolBar( wxWindow *parent, wxWindowID id ) : wxPanel( parent, id ) {

    // create horizontal sizer to contain buttons
    m_sizer = new wxBoxSizer( wxHORIZONTAL );
    SetSizer( m_sizer );
}


/// add a button that runs the given command
void ToolBar::addButton( const String &caption, const String &command, const String &tip ) {

    // add button
    wxButton *button = new wxButton( this, m_commands.count(), caption.c_str() );
    if (tip.length())
        button->SetToolTip( tip.c_str() );
    button->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ToolBar::onButton ), NULL, this );
    m_sizer->Add( button );

    // store associated command
    m_commands.appendCopy( command );
}


// called when user presses one of the buttons
void ToolBar::onButton( wxCommandEvent &event ) {
    int id = event.GetId();
    if (id >= 0 && id < m_commands.count()) {
        const String &command = m_commands[ id ];
        execCommand( command, false );
    }
}


//-------------------------------------------
// DIALOG PARENT
//-------------------------------------------


// the parent widget to use for dialog widgets
wxFrame *g_dialogParent = NULL;


/// set the parent widget to use for dialog widgets
void setDialogParent( wxFrame *parent ) {
    g_dialogParent = parent;
}


/// get the parent widget to use for dialog widgets
wxFrame *dialogParent() {
    assertAlways( g_dialogParent );
    return g_dialogParent;
}


} // end namespace sbl

