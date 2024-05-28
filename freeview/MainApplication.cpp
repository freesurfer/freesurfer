/**
 * @brief Application object.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
 *
 */


#include "MainApplication.h"
#include <QKeyEvent>
#include <QDebug>

MainApplication::MainApplication( int & argc, char ** argv ) :
  QApplication(argc, argv)
{
}


void MainApplication::SetLargeFont(bool b)
{
  QString strg;
  if (b)
  {
#ifdef Q_OS_MAC
    strg = "QWidget {font-size: 14pt;}";
#else
    strg = "QWidget {font-size: 12pt;}";
#endif
  }
  setStyleSheet(strg);
}
