#ifndef _SBL_CONFIG_EDITOR_H_
#define _SBL_CONFIG_EDITOR_H_
#include <sbl/core/Config.h>
#include <wx/wx.h>
namespace sbl {


//-------------------------------------------
// CONFIG ENTRY EDITOR CLASS
//-------------------------------------------


/// The ConfigEntryEditor widget provides an interface for editing a single Config entry.
class ConfigEntryEditor : public wxPanel {
public:

    // basic constructor
    ConfigEntryEditor( wxWindow *parent, wxWindowID id, ConfigEntry &configEntry );

    /// store changes from user interface into config object
    void storeChanges();

    // button callback for file/path config entries
    void onBrowseButton( wxCommandEvent &event );

private:

    // the entry being edited
    ConfigEntry &m_configEntry;

    // widgets holding the edited value
    wxTextCtrl *m_text;
    wxCheckBox *m_checkBox;

    // disable copy constructor and assignment operator
    ConfigEntryEditor( const ConfigEntryEditor &x );
    ConfigEntryEditor &operator=( const ConfigEntryEditor &x );
};


//-------------------------------------------
// CONFIG EDITOR CLASS
//-------------------------------------------


/// The ConfigEditor widget provides an interface for editing a Config object.
class ConfigEditor : public wxPanel {
public:

    /// basic constructor
    ConfigEditor( wxWindow *parent, wxWindowID id, Config &config );

    /// store changes from user interface into config object
    void storeChanges();

private:

    // the config being edited
    Config &m_config;

    // an editor widget for each (non-blank) entry
    PtrArray<ConfigEntryEditor> m_entryEditors;

    // disable copy constructor and assignment operator
    ConfigEditor( const ConfigEditor &x );
    ConfigEditor &operator=( const ConfigEditor &x );
};


//-------------------------------------------
// CONFIG EDITOR DIALOG CLASS
//-------------------------------------------


/// The ConfigEditorDialog class creates a pop-up window for editing a Config object.
class ConfigEditorDialog : public wxDialog {
public:

    // basic constructor
    ConfigEditorDialog( wxWindow *parent, wxWindowID id, Config &config, const String &title );

    // ok button callback; save changes
    void onOkButton( wxCommandEvent &event );

    // cancel button callback; don't save changes
    void onCancelButton( wxCommandEvent &event );

private:

    // the editor widget contained within this dialog
    ConfigEditor *m_configEditor;

    // disable copy constructor and assignment operator
    ConfigEditorDialog( const ConfigEditorDialog &x );
    ConfigEditorDialog &operator=( const ConfigEditorDialog &x );
};


//-------------------------------------------
// EXTERNAL INTERFACE
//-------------------------------------------


/// open a window for editing the given config; returns false if user cancels
bool editConfig( Config &conf, const String &title );


} // end namespace sbl
#endif // _SBL_CONFIG_EDITOR_H_

