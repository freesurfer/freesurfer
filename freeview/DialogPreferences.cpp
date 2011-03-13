/**
 * @file  DialogPreferences.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:17 $
 *    $Revision: 1.15 $
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

#include "DialogPreferences.h"
#include "ui_DialogPreferences.h"
#include "MainWindow.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "Cursor2D.h"
#include "Cursor3D.h"
#include "TermWidget.h"

DialogPreferences::DialogPreferences(QWidget *parent) :
  QDialog(parent),
  UIUpdateHelper(),
  ui(new Ui::DialogPreferences)
{
  ui->setupUi(this);
  ui->buttonBox->button(QDialogButtonBox::Close)->setDefault(true);
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  for (int i = 0; i < 4; i++)
  {
    connect(ui->colorPickerBackground, SIGNAL(colorChanged(QColor)),
            mainwnd->GetRenderView(i), SLOT(SetBackgroundColor(QColor)));
  }
  for (int i = 0; i < 3; i++)
  {
    connect(ui->colorPickerCursor, SIGNAL(colorChanged(QColor)),
            ((RenderView2D*)mainwnd->GetRenderView(i))->GetCursor2D(), SLOT(SetColor(QColor)));
    connect(ui->comboBoxCursorStyle, SIGNAL(currentIndexChanged(int)),
            ((RenderView2D*)mainwnd->GetRenderView(i))->GetCursor2D(), SLOT(SetStyle(int)));
  }
  connect(ui->colorPickerCursor, SIGNAL(colorChanged(QColor)),
          ((RenderView3D*)mainwnd->GetRenderView(3))->GetCursor3D(), SLOT(SetColor(QColor)));
  connect(ui->checkBoxSyncZoom, SIGNAL(toggled(bool)),
          mainwnd, SLOT(SyncZoom(bool)));
  connect(ui->radioButtonThemeDark, SIGNAL(toggled(bool)),
          mainwnd->GetCommandConsole(), SLOT(SetDarkTheme(bool)));

#ifdef Q_WS_MAC
  ui->groupBoxMac->setEnabled(true);
  ui->groupBoxMac->show();
  connect(ui->checkBoxMacUnified, SIGNAL(toggled(bool)), mainwnd, SLOT(SetUnifiedTitleAndToolBar(bool)));
  connect(ui->checkBoxCommandKey, SIGNAL(toggled(bool)), mainwnd, SLOT(SetUseCommandControl(bool)));
#else
  ui->groupBoxMac->setEnabled(false);
  ui->groupBoxMac->hide();
#endif
}

DialogPreferences::~DialogPreferences()
{
  delete ui;
}

void DialogPreferences::SetSettings(const QVariantMap &map)
{
  BlockAllSignals(this, true);
  ui->colorPickerBackground->setCurrentColor(map["BackgroundColor"].value<QColor>());
  ui->colorPickerCursor->setCurrentColor(map["CursorColor"].value<QColor>());
  ui->comboBoxCursorStyle->setCurrentIndex(map["CursorStyle"].toInt());
  ui->checkBoxSaveCopy->setChecked(map["SaveCopy"].toBool());
  ui->checkBoxSyncZoom->setChecked(map["SyncZoom"].toBool());
  ui->checkBoxCommandKey->setChecked(map["MacUseCommand"].toBool());
  ui->checkBoxMacUnified->setChecked(map["MacUnifiedTitleBar"].toBool());
  ui->radioButtonThemeDark->setChecked(map["DarkConsole"].toBool());
  ui->radioButtonThemeLight->setChecked(!map["DarkConsole"].toBool());
  BlockAllSignals(this, false);
}

QVariantMap DialogPreferences::GetSettings()
{
  QVariantMap map;
  map["BackgroundColor"] = ui->colorPickerBackground->currentColor();
  map["CursorColor"] = ui->colorPickerCursor->currentColor();
  map["CursorStyle"] = ui->comboBoxCursorStyle->currentIndex();
  map["SaveCopy"] = ui->checkBoxSaveCopy->isChecked();
  map["SyncZoom"] = ui->checkBoxSyncZoom->isChecked();
  map["MacUseCommand"] = ui->checkBoxCommandKey->isChecked();
  map["MacUnifiedTitleBar"] = ui->checkBoxMacUnified->isChecked();
  map["DarkConsole"] = ui->radioButtonThemeDark->isChecked();
  return map;
}

void DialogPreferences::OnClicked(QAbstractButton* btn)
{
  if (ui->buttonBox->buttonRole(btn) == QDialogButtonBox::ResetRole)
  {
    ui->colorPickerBackground->setCurrentColor(Qt::black);
    ui->colorPickerCursor->setCurrentColor(Qt::red);
    ui->comboBoxCursorStyle->setCurrentIndex(0);
    ui->checkBoxSaveCopy->setChecked(true);
    ui->checkBoxSyncZoom->setChecked(true);
    ui->radioButtonThemeDark->setChecked(true);
#ifdef Q_WS_MAC
    ui->checkBoxCommandKey->setChecked(false);
    ui->checkBoxMacUnified->setChecked(false);
#endif
  }
}
