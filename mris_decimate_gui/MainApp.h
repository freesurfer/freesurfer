/**
 * @file  MainApp.h
 * @brief Main application.
 *
 */
/*
 * Original Author: Dan Ginsburg
 * CVS Revision Info:
 *    $Author: ginsburg $
 *    $Date: 2010/08/04 20:38:52 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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


