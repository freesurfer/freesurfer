/**
 * @file  WindowGroupPlot.cpp
 * @brief Tool window to plot group data
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/03/29 20:35:50 $
 *    $Revision: 1.1 $
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

#include "WindowGroupPlot.h"
#include "ui_WindowGroupPlot.h"
#include <QSettings>
#include "FSGroupDescriptor.h"

WindowGroupPlot::WindowGroupPlot(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WindowGroupPlot)
{
    ui->setupUi(this);
    this->setWindowFlags(Qt::Tool);
    this->setWindowTitle("Plot");

    QSettings s;
    QVariant v = s.value("WindowPlot/Geomerty");
    if (v.isValid())
      this->restoreGeometry(v.toByteArray());
}

WindowGroupPlot::~WindowGroupPlot()
{
  QSettings s;
  s.setValue("WindowPlot/Geomerty", this->saveGeometry());
  delete ui;
}

void WindowGroupPlot::SetFsgdData(FSGroupDescriptor *fsgd)
{
  ui->widgetPlot->SetFsgdData(fsgd);
  ui->comboBoxVariable->blockSignals(true);
  ui->comboBoxVariable->clear();
  for (int i = 0; i < fsgd->m_variables.size(); i++)
    ui->comboBoxVariable->addItem(fsgd->m_variables[i].label);
  ui->comboBoxVariable->setCurrentIndex(0);
  ui->comboBoxVariable->blockSignals(false);
}
