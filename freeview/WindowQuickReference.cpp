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
 */
#include "WindowQuickReference.h"
#include "ui_WindowQuickReference.h"
#include <QFile>
#include <QSettings>

WindowQuickReference::WindowQuickReference(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::WindowQuickReference)
{
  ui->setupUi(this);
  setWindowFlags( Qt::Tool );
  QFile file(":/resource/QuickRef.html");
  file.open(QIODevice::ReadOnly | QIODevice::Text);
  QByteArray ba = file.readAll();
#ifdef Q_OS_MAC
  ba.replace("Ctrl +", "Cmd +");
#endif
  ui->textBrowser->setHtml(ba);

  QSettings settings;
  restoreGeometry(settings.value("WindowQuickRef/Geometry").toByteArray());
}

WindowQuickReference::~WindowQuickReference()
{
  QSettings settings;
  settings.setValue("WindowQuickRef/Geometry", this->saveGeometry());
  delete ui;
}
