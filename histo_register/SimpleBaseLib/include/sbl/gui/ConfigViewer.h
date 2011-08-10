#ifndef _SBL_CONFIG_VIEWER_H_
#define _SBL_CONFIG_VIEWER_H_
#include <sbl/core/Config.h>
#include <wx/wx.h>
namespace sbl {


// register commands, etc. defined in this module
void initConfigViewer();


//-------------------------------------------
// CONFIG VIEWER CLASS
//-------------------------------------------


/// The ConfigViewer displays a set of config entries, but does not allow editing.
class ConfigViewer : public wxPanel {
public:

    // basic constructor
    ConfigViewer( wxWindow *parent, wxWindowID id );

    /// display a Config object for editing
    void init( Config &conf );

    /// updates values of currently displayed config entries;
    /// assumes order is unchanged; adds missing entries, assuming those are at end
    void updateValues( Config &conf );

private:

    /// add a single config entry
    void addEntry( const ConfigEntry &entry );

    // the sizer containing the names and values
    wxFlexGridSizer *m_topSizer;

    // the widgets displaying the values
    PtrArray<wxStaticText> m_valueControls;
};


//-------------------------------------------
// CONFIG VIEWER DIALOG CLASS
//-------------------------------------------


/// The ConfigViewerDialog class wraps the ConfigViewer widget in a dialog box.
class ConfigViewerDialog : public wxDialog {
public:

    // basic constructor
    ConfigViewerDialog( wxWindow *parent, wxWindowID id, Config &conf, const String &title );

    /// the config viewer widget
    ConfigViewer &configViewer() { return *m_configViewer; }

    /// the top-level sizer containing the config viewer widget
    wxBoxSizer &sizer() { return *m_sizer; }

private:

    // internal data
    wxBoxSizer *m_sizer;
    ConfigViewer *m_configViewer;
};


} // end namespace sbl
#endif // _SBL_CONFIG_VIEWER_H_

