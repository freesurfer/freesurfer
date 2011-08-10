#ifndef _SBL_MAIN_WINDOW_H_
#define _SBL_MAIN_WINDOW_H_
#include <wx/wx.h>
#include <wx/notebook.h>
#include <sbl/core/String.h>
#include <sbl/core/Dict.h>
#include <sbl/gui/MiscWidgets.h>
namespace sbl {


/// The KeyboardShortcut struct holds a command associated with a keyboard code.
struct KeyboardShortcut {
    KeyboardShortcut( const String &_command, int _modifier ) { command = _command; modifier = _modifier; }
    String command;
    int modifier;
};


/// The TimerEvent struct represents a callback to be run periodically by the wxWidgets timer code.
struct TimerEvent {
    int msecInterval;
    void (*callback)();
};


//-------------------------------------------
// MAIN WINDOW CLASS
//-------------------------------------------


/// The MainWindow class is the top-level interface window.  
/// It provides a command-entry control and tabs for assorted viewers.
class MainWindow : public wxFrame {
public:

    // basic constructor/destructor
    MainWindow( const wxString &title, const wxPoint &pos, const wxSize &size );
    ~MainWindow();

    //-------------------------------------------
    // TABS
    //-------------------------------------------

    /// the notebook control holds tabs for various viewers
    wxNotebook &notebook() { return *m_notebook; }

    /// switch to a given tab
    void switchTo( int tabIndex );
    void switchTo( const wxWindow *window );

    /// add a tab to the notebook
    void addTab( wxWindow *window, const String &name );

    //-------------------------------------------
    // DISPLAY 
    //-------------------------------------------

    /// display text in the command output viewer
    void dispText( const char *text );

    //-------------------------------------------
    // KEYBOARD SHORTCUTS
    //-------------------------------------------

    /// add a keyboard shortcut that runs a command (uses wxWidgets key codes);
    /// assumes no more than one shortcut per code
    inline void addKeyboardShortcut( int code, const String &command ) { m_keyboardShortcuts.add( code, new KeyboardShortcut( command, 0 ) ); }

    /// add a keyboard shortcut that runs a command (uses wxWidgets key codes and modifiers);
    /// assumes no more than one shortcut per code
    inline void addKeyboardShortcut( int modifier, int code, const String &command ) { m_keyboardShortcuts.add( code, new KeyboardShortcut( command, modifier ) ); }

    /// if command found for this code, execute it
    void execKeyboardShortcut( int modifier, int code );

    //-------------------------------------------
    // OTHER METHODS
    //-------------------------------------------

    /// set focus on command entry widget
    inline void setFocusCommandInput() { m_commandInputText->SetFocus(); }

    // handle keyboard input for command entry widget
    void onCommandInputKey( wxKeyEvent &event );

    // handle mouse wheel movement
    void onMouseWheel( wxMouseEvent &event );

    // handle timer event to run register timed callbacks
    void onTimer( wxTimerEvent &event );

    /// register a function to be periodically called using wxWidgets timers
    /// (currently restricted to msecInterval == 100)
    void addTimerEvent( int msecInterval, void (*callback)() );

    /// access current (single, global) MainWindow instance
    static MainWindow *instance() { return s_instance; }

private:

    // internal data
    wxNotebook *m_notebook;
    wxTextCtrl *m_commandInputText;
    wxListBox *m_commandOutputList;
    wxTimer *m_timer;
    Dict<KeyboardShortcut> m_keyboardShortcuts;
    Array<TimerEvent> m_timerEvents;

    // current (single, global) instance of MainWindow
    static MainWindow *s_instance;
};


//-------------------------------------------
// MAIN APP CLASS
//-------------------------------------------


/// The MainApp class is the top-level object for the user interface; create a sub-class to specify own start-up code.
class MainApp : public wxApp {
public:

    /// run all GUI SBL initialization 
    virtual bool OnInit();

    // override event filtering for keyboard handling
    int FilterEvent( wxEvent &event );

    /// access current (single, global) MainApp instance
    static MainApp *instance() { return s_instance; }

private:

    // current (single, global) instance of MainApp
    static MainApp *s_instance;
};


} // end namespace sbl
#endif // _SBL_MAIN_WINDOW_H_

