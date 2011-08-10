// Licensed under MIT license; see license.txt.

#include <sbl/gui/ConfigViewer.h>
#include <sbl/core/Command.h> // for test command
#include <sbl/gui/MiscWidgets.h> // for test command
namespace sbl {


//-------------------------------------------
// CONFIG VIEWER CLASS
//-------------------------------------------


// basic constructor
ConfigViewer::ConfigViewer( wxWindow *parent, wxWindowID id ) : wxPanel( parent, id ) {

    // create top-level sizer; will have two colums
    m_topSizer = new wxFlexGridSizer( 2, 0, 15 );
    SetSizer( m_topSizer );
}


/// display a Config object for editing
void ConfigViewer::init( Config &conf ) {

    // add each config entry
    for (int i = 0; i < conf.entryCount(); i++) {
        const ConfigEntry &entry = conf.entry( i );
        if (entry.type != CONFIG_ENTRY_BLANK)
            addEntry( entry );
    }
}


/// updates values of currently displayed config entries;
/// assumes order is unchanged; adds missing entries, assuming those are at end
void ConfigViewer::updateValues( Config &conf ) {

    // add any missing entries (assumes only added entries at end)
    if (conf.entryCount() > m_valueControls.count()) {
        while (m_valueControls.count() < conf.entryCount()) {
            addEntry( conf.entry( m_valueControls.count() ) );
        }
        SetSizerAndFit( m_topSizer );
    }

    // update values
    for (int i = 0; i < conf.entryCount(); i++) 
        m_valueControls[ i ].SetLabel( conf.entry( i ).value.c_str() );
}


/// add a single config entry
void ConfigViewer::addEntry( const ConfigEntry &entry ) {

    // create widgets for name and value
    wxStaticText *nameText = new wxStaticText( this, -1, entry.name.c_str() );
    wxStaticText *valueText = new wxStaticText( this, -1, entry.value.c_str() );
    m_valueControls.append( valueText ); // just adding a blank entry to keep numbering correct

    // set appearance based on whether section heading
    int topBorder = 0;
    if (entry.type == CONFIG_ENTRY_SECTION) {
        wxFont headingFont( 10, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD );
        nameText->SetFont( headingFont );
        if (m_valueControls.count() > 1)
            topBorder = 10;
    } else {
        wxFont boldFont( 10, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD );
        valueText->SetFont( boldFont );
    }

    // add to sizer
    m_topSizer->Add( nameText, 1, wxALIGN_LEFT | wxEXPAND | wxTOP, topBorder );
    m_topSizer->Add( valueText, 1, wxEXPAND );
}


//-------------------------------------------
// CONFIG VIEWER DIALOG CLASS
//-------------------------------------------


// basic constructor
ConfigViewerDialog::ConfigViewerDialog( wxWindow *parent, wxWindowID id, Config &conf, const String &title ) : wxDialog( parent, id, wxString( title.c_str() ) ) {
    m_sizer = new wxBoxSizer( wxVERTICAL );
    m_configViewer = new ConfigViewer( this, -1 );
    m_configViewer->init( conf );
    m_sizer->Add( m_configViewer, 1, wxALL, 10 );
    SetSizerAndFit( m_sizer );
}


//-------------------------------------------
// TEST COMMAND
//-------------------------------------------


// test the config viewer as a dialog box
void testConfigViewer( Config &conf ) {
    
    // create a config for testing the viewer
    Config testConf;
    testConf.writeSection( "First Section" );
    testConf.writeString( "aString", "some text" );
    testConf.writeInt( "anInt", 137 );
    testConf.writeSection( "Second Section" );
    testConf.writeDouble( "aDouble", 3.1415 );
    testConf.writeBool( "aBool", true );

    // show the dialog
    ConfigViewerDialog *configViewerDialog = new ConfigViewerDialog( dialogParent(), -1, testConf, "Config Viewer Test" );
    configViewerDialog->ShowModal();
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initConfigViewer() {
#ifdef REGISTER_TEST_COMMANDS
    registerCommand( "testconfview", testConfigViewer );
#endif
}


} // end namespace sbl

