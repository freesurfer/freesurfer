/**
 * @file  DialogAbout.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:17 $
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

#include "DialogAbout.h"
#include "ui_DialogAbout.h"

DialogAbout::DialogAbout(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogAbout)
{
  ui->setupUi(this);
  QString strg = ui->labelVersion->text();
  strg.replace("xxx", QString(__DATE__) + " " + __TIME__);
  ui->labelVersion->setText(strg);
}

DialogAbout::~DialogAbout()
{
  delete ui;
}
