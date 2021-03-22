/**
 * @brief Abstract class used to post progress messages and status bar updates
 *
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef ProgressUpdateGUI_H
#define ProgressUpdateGUI_H

class ProgressUpdateGUI {
public:
  virtual ~ProgressUpdateGUI () {};
  virtual void BeginActionWithProgress ( const char * isTitle ) = 0;
  virtual void UpdateProgressMessage ( const char* isMessage ) = 0;
  virtual void UpdateProgressPercent ( float iPercent ) = 0; // 0 - 100
  virtual void EndActionWithProgress () = 0;
};

#endif // ProgressUpdateGUI_H
