/**
 * @file  MainApp.h
 * @brief Main application.
 *
 */
/*
 * Original Author: Dan Ginsburg
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:30 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef MainApp_h
#define MainApp_h

#include <wx/wx.h>

class MainWindow;

// Define a new application type, each program should derive a class from wxApp
class MainApp : public wxApp
{
public:
    MainApp();
    virtual ~MainApp();
    void CreateMainWindow();

    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit();
    virtual int OnExit();


private:
    void InitializeGUI();

    MainWindow* m_wndMain;

    DECLARE_EVENT_TABLE()
};

DECLARE_APP( MainApp )

#endif // MainApp_H


