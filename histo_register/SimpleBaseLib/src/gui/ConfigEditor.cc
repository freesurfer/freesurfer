// Licensed under MIT license; see license.txt.

#include <sbl/gui/ConfigEditor.h>
#include <sbl/core/Command.h>
#include <sbl/system/FileSystem.h>
#include <sbl/gui/MainWindow.h>
#include <sbl/gui/MiscWidgets.h>
namespace sbl {


//-------------------------------------------
// CONFIG ENTRY EDITOR CLASS
//-------------------------------------------


// basic constructor
ConfigEntryEditor::ConfigEntryEditor( wxWindow *parent, wxWindowID id, ConfigEntry &configEntry ) : wxPanel( parent, id ), m_configEntry( configEntry ) {
    m_checkBox = NULL;
    m_text = NULL;

    // create top-level sizer
    wxBoxSizer *sizer = new wxBoxSizer( wxHORIZONTAL );
    SetSizer( sizer );

    // add label
    String desc = configEntry.type == CONFIG_ENTRY_SECTION ? configEntry.name : descriptionFromName( configEntry.name );
    wxStaticText *label = new wxStaticText( this, -1, desc.c_str(), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    if (configEntry.type == CONFIG_ENTRY_SECTION) {
        wxFont boldFont( 10, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD );
        label->SetFont( boldFont );
        label->SetWindowStyle( wxALIGN_LEFT );
        label->SetMinSize( wxSize( 100, 22 ) );
    } else {
        label->SetMinSize( wxSize( 100, 18 ) );
    }
    sizer->Add( label, 1, wxEXPAND );
    sizer->AddSpacer( 5 );

    // add buttons for file/path
    if (configEntry.type == CONFIG_ENTRY_PATH || configEntry.type == CONFIG_ENTRY_FILE) {
        m_text = new wxTextCtrl( this, -1, configEntry.value.c_str() );
        sizer->Add( m_text, 3, wxEXPAND );
        wxButton *button = new wxButton( this, -1, "Browse" );
        sizer->Add( button, 1, wxEXPAND | wxLEFT, 5 );
        button->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ConfigEntryEditor::onBrowseButton ), NULL, this );

    // add checkbox for bool entries
    } else if (configEntry.type == CONFIG_ENTRY_BOOL) {
        m_checkBox = new wxCheckBox( this, -1, "" );
        m_checkBox->SetValue( configEntry.value.toBool() );
        sizer->Add( m_checkBox, 1, wxEXPAND );

    // otherwise, just use a text box
    } else if (configEntry.type != CONFIG_ENTRY_SECTION) {
        m_text = new wxTextCtrl( this, -1, configEntry.value.c_str() );
        sizer->Add( m_text, 4, wxEXPAND );
    }
}


/// store changes from user interface into config object
void ConfigEntryEditor::storeChanges() {
    if (m_configEntry.type == CONFIG_ENTRY_BOOL) {
        assertAlways( m_checkBox );
        if (m_checkBox->GetValue()) 
            m_configEntry.value = "1";
        else
            m_configEntry.value = "0";
    } else {
        assertAlways( m_text );
        m_configEntry.value = String( m_text->GetValue().mb_str() );
    }
}


// button callback for file/path config entries
void ConfigEntryEditor::onBrowseButton( wxCommandEvent &event ) {
    if (m_configEntry.type == CONFIG_ENTRY_PATH) {
        wxDirDialog *dirDialog = new wxDirDialog( this, "Select a path", m_configEntry.value.c_str() );
        if (dirDialog->ShowModal() == wxID_OK) {
            String path = dirDialog->GetPath();
            m_text->SetValue( path.c_str() );
        }
    } else {
        wxFileDialog *fileDialog = new wxFileDialog( this, "Select a file" );
        fileDialog->SetDirectory( pathPart( m_configEntry.value ).c_str() ); 
        if (fileDialog->ShowModal() == wxID_OK) {
            String fileName = fileDialog->GetPath();
            m_text->SetValue( fileName.c_str() );
        }
    }
}


//-------------------------------------------
// CONFIG EDITOR CLASS
//-------------------------------------------


// basic constructor
ConfigEditor::ConfigEditor( wxWindow *parent, wxWindowID id, Config &config ) : wxPanel( parent, id ), m_config( config ) {

    // create top-level sizer
    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );

    // add each config entry
    for (int i = 0; i < config.entryCount(); i++) {
        ConfigEntry &entry = config.entry( i );
        if (entry.type != CONFIG_ENTRY_BLANK) {
            ConfigEntryEditor *entryEditor = new ConfigEntryEditor( this, -1, entry );
            if (entry.type == CONFIG_ENTRY_SECTION && m_entryEditors.count() > 0) 
                topSizer->AddSpacer( 10 );
            topSizer->Add( entryEditor );
            m_entryEditors.append( entryEditor );
        }
    }
    SetSizerAndFit( topSizer );
}


/// store changes from user interface into config object
void ConfigEditor::storeChanges() {
    for (int i = 0; i < m_entryEditors.count(); i++) 
        m_entryEditors[ i ].storeChanges();
}


//-------------------------------------------
// CONFIG EDITOR DIALOG CLASS
//-------------------------------------------


// basic constructor
ConfigEditorDialog::ConfigEditorDialog( wxWindow *parent, wxWindowID id, Config &config, const String &title ) : wxDialog( parent, id, wxString( title.c_str() ), wxDefaultPosition, wxDefaultSize, wxCAPTION ) {

    // create main sizer
    wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );

    // add config editor
    m_configEditor = new ConfigEditor( this, -1, config );
    topSizer->Add( m_configEditor, 1, wxEXPAND | wxALL, 5 );

    // add ok/cancel buttons
    wxBoxSizer *buttonSizer = new wxBoxSizer( wxHORIZONTAL );
    topSizer->Add( buttonSizer, 0, wxEXPAND );
    wxButton *okButton = new wxButton( this, -1, "OK" );
    okButton->SetDefault();
    wxButton *cancelButton = new wxButton( this, -1, "Cancel" );
    buttonSizer->AddStretchSpacer( 2 );
    buttonSizer->Add( okButton, 1, wxEXPAND | wxBOTTOM | wxLEFT, 5 );
    buttonSizer->Add( cancelButton, 1, wxEXPAND | wxBOTTOM | wxRIGHT | wxLEFT, 5 );

    // finish layout
    SetSizerAndFit( topSizer );

    // connect event handlers
    okButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ConfigEditorDialog::onOkButton ), NULL, this );
    cancelButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ConfigEditorDialog::onCancelButton ), NULL, this );
}


// ok button callback; save changes
void ConfigEditorDialog::onOkButton( wxCommandEvent &event ) {
    m_configEditor->storeChanges();
    EndModal( 1 );
}


// cancel button callback; don't save changes
void ConfigEditorDialog::onCancelButton( wxCommandEvent &event ) {
    EndModal( 0 );
}


//-------------------------------------------
// EXTERNAL INTERFACE
//-------------------------------------------


/// open a window for editing the given config; returns false if user cancels
bool editConfig( Config &conf, const String &title ) {
    ConfigEditorDialog *configEditorDialog = new ConfigEditorDialog( dialogParent(), -1, conf, title );
    int ret = configEditorDialog->ShowModal();
    return ret ? true : false;
}


} // end namespace sbl

