#ifndef _SBL_MISC_WIDGETS_H_
#define _SBL_MISC_WIDGETS_H_
#include <wx/wx.h>
#include <sbl/core/String.h>
namespace sbl {


//-------------------------------------------
// PROGRESS DIALOG CLASS
//-------------------------------------------


/// The ProgressDialog class is a dialog window that shows progress of a long-running process.
class ProgressDialog : public wxDialog {
public:

    // basic constructor
    ProgressDialog( wxWindow *parent, wxWindowID id );

    /// update the displayed progress;
    /// (assumes index in [0, count - 1]);
    /// closes progress bar when index == count - 1; 
    /// if count == -1, assumes unknown number if items
    void updateProgress( int index, int count );

    /// called when user presses cancel button
    void onCancelButton( wxCommandEvent &event );

    // display or update a progress dialog
    static void dispProgress( int index, int count );

private:

    /// the label showing the current progress
    wxStaticText *m_progressLabel;
};


//-------------------------------------------
// TOOL BAR CLASS
//-------------------------------------------


/// The ToolBar widget shows a set of buttons; each button executes a command.
class ToolBar : public wxPanel {
public:

    // basic constructor
    ToolBar( wxWindow *parent, wxWindowID id );

    /// add a button that runs the given command
    void addButton( const String &caption, const String &command, const String &tip );

    // called when user presses one of the buttons
    void onButton( wxCommandEvent &event );    

private:

    // the top-level sizer used for adding buttons
    wxBoxSizer *m_sizer;

    // the command associated with each button
    Array<String> m_commands;
};


//-------------------------------------------
// DIALOG PARENT
//-------------------------------------------


/// set the parent widget to use for dialog widgets
void setDialogParent( wxFrame *parent );


/// get the parent widget to use for dialog widgets
wxFrame *dialogParent();


} // end namespace sbl
#endif // _SBL_MISC_WIDGETS_H_

