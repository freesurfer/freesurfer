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
#include "DialogPreferences.h"
#include "ui_DialogPreferences.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "Cursor2D.h"
#include "Cursor3D.h"
#include "Annotation2D.h"
#include "TermWidget.h"
#include <QMessageBox>

DialogPreferences::DialogPreferences(QWidget *parent) :
  QDialog(parent),
  UIUpdateHelper(),
  ui(new Ui::DialogPreferences)
{
  ui->setupUi(this);
  ui->buttonBox->button(QDialogButtonBox::Close)->setDefault(true);
  ui->groupBoxCommandConsole->hide();

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  ui->comboBoxShortcutCycleLayer->setProperty("action", QVariant::fromValue<QObject*>(mainwnd->ui->actionCycleLayer));
  ui->comboBoxShortcutToggleSurface->setProperty("action", QVariant::fromValue<QObject*>(mainwnd->ui->actionToggleSurfaceVisibility));
  ui->comboBoxShortcutToggleVolume->setProperty("action", QVariant::fromValue<QObject*>(mainwnd->ui->actionToggleVolumeVisibility));

  for (int i = 0; i < 4; i++)
  {
    connect(ui->colorPickerBackground, SIGNAL(colorChanged(QColor)),
            mainwnd->GetRenderView(i), SLOT(SetBackgroundColor(QColor)));
  }
  for (int i = 0; i < 3; i++)
  {
    connect(ui->colorPickerCursor, SIGNAL(colorChanged(QColor)),
            ((RenderView2D*)mainwnd->GetRenderView(i))->GetCursor2D(), SLOT(SetColor(QColor)));
    connect(ui->horizontalSliderSize2D, SIGNAL(valueChanged(int)),
            ((RenderView2D*)mainwnd->GetRenderView(i))->GetCursor2D(), SLOT(SetSize(int)));
    connect(ui->horizontalSliderThickness2D, SIGNAL(valueChanged(int)),
            ((RenderView2D*)mainwnd->GetRenderView(i))->GetCursor2D(), SLOT(SetThickness(int)));
    connect(ui->colorPickerAnnotation, SIGNAL(colorChanged(QColor)),
            ((RenderView2D*)mainwnd->GetRenderView(i))->GetAnnotation2D(), SLOT(SetColor(QColor)));
    connect(ui->spinBoxFontSize, SIGNAL(valueChanged(int)),
            ((RenderView2D*)mainwnd->GetRenderView(i)), SLOT(SetTextSize(int)));
    connect(ui->checkBoxAutoScaleFont, SIGNAL(toggled(bool)),
            ((RenderView2D*)mainwnd->GetRenderView(i)), SLOT(SetAutoScaleText(bool)));
  }
  connect(ui->colorPickerCursor, SIGNAL(colorChanged(QColor)),
          ((RenderView3D*)mainwnd->GetRenderView(3))->GetCursor3D(), SLOT(SetColor(QColor)));
  connect(ui->colorPickerCursor, SIGNAL(colorChanged(QColor)),
          ((RenderView3D*)mainwnd->GetRenderView(3))->GetInflatedSurfCursor(), SLOT(SetColor(QColor)));
  connect(ui->horizontalSliderSize3D, SIGNAL(valueChanged(int)),
          ((RenderView3D*)mainwnd->GetRenderView(3))->GetCursor3D(), SLOT(SetSize(int)));
  connect(ui->horizontalSliderSize3D, SIGNAL(valueChanged(int)),
          ((RenderView3D*)mainwnd->GetRenderView(3))->GetInflatedSurfCursor(), SLOT(SetSize(int)));
  connect(ui->horizontalSliderThickness3D, SIGNAL(valueChanged(int)),
          ((RenderView3D*)mainwnd->GetRenderView(3))->GetCursor3D(), SLOT(SetThickness(int)));
  connect(ui->horizontalSliderThickness3D, SIGNAL(valueChanged(int)),
          ((RenderView3D*)mainwnd->GetRenderView(3))->GetInflatedSurfCursor(), SLOT(SetThickness(int)));
  connect(ui->checkBoxSyncZoom, SIGNAL(toggled(bool)),
          mainwnd, SLOT(SyncZoom(bool)));
  connect(ui->radioButtonThemeDark, SIGNAL(toggled(bool)),
          mainwnd->GetCommandConsole(), SLOT(SetDarkTheme(bool)));

  connect(ui->comboBox3DScaleStyle, SIGNAL(currentIndexChanged(int)), ((RenderView3D*)mainwnd->GetRenderView(3)), SLOT(SetAxesFlyMode(int)));

  ui->checkBoxMacUnified->hide();

#ifdef Q_OS_MAC
  ui->groupBoxMac->setEnabled(true);
  ui->groupBoxMac->show();
//  connect(ui->checkBoxMacUnified, SIGNAL(toggled(bool)), mainwnd, SLOT(SetUnifiedTitleAndToolBar(bool)));
  connect(ui->checkBoxCommandKey, SIGNAL(toggled(bool)), mainwnd, SLOT(SetUseCommandControl(bool)));
#else
  ui->groupBoxMac->setEnabled(false);
  ui->groupBoxMac->hide();
#endif

  connect(ui->colorPickerAnnotation, SIGNAL(colorChanged(QColor)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->colorPickerBackground, SIGNAL(colorChanged(QColor)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->colorPickerCursor, SIGNAL(colorChanged(QColor)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->horizontalSliderSize2D, SIGNAL(valueChanged(int)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->horizontalSliderSize3D, SIGNAL(valueChanged(int)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxRightButtonErase, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxSaveCopy, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxSyncZoom, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->radioButtonThemeDark, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxAutoReorientView, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxDecimalVoxelCoord, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxAutoScaleFont, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->spinBoxFontSize, SIGNAL(valueChanged(int)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxAutoMidToMin, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->spinBoxPrecision, SIGNAL(valueChanged(int)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxComma, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxClickToLock, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));
  connect(ui->checkBoxAllowDeleteKey, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateSettings()));

  connect(ui->spinBoxPrecision, SIGNAL(valueChanged(int)), mainwnd, SLOT(UpdateInfoPanel()), Qt::QueuedConnection);
  connect(ui->checkBoxComma, SIGNAL(toggled(bool)), mainwnd, SLOT(UpdateInfoPanel()), Qt::QueuedConnection);

  QList<QComboBox*> list_combos;
  list_combos << ui->comboBoxShortcutCycleLayer << ui->comboBoxShortcutToggleSurface
              << ui->comboBoxShortcutToggleVolume;
  foreach (QComboBox* combo, list_combos)
  {
    combo->clear();
    QAction* act = qobject_cast<QAction*>(combo->property("action").value<QObject*>());
    if (act)
    {
      combo->addItem(act->shortcut().toString());
      for (int i = 1; i <= 12; i++)
        combo->addItem(tr("F%1").arg(i));
      combo->addItem(tr("Pause"));
      combo->setCurrentIndex(0);
    }
    connect(combo, SIGNAL(currentIndexChanged(QString)), this, SLOT(OnComboShortcutChanged(QString)));
  }
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
  ui->horizontalSliderSize2D->setValue(map["CursorSize"].toInt());
  ui->horizontalSliderSize3D->setValue(map["CursorSize3D"].toInt());
  ui->horizontalSliderThickness2D->setValue(map["CursorThickness"].toInt());
  ui->horizontalSliderThickness3D->setValue(map["CursorThickness3D"].toInt());
  ui->checkBoxSaveCopy->setChecked(map["SaveCopy"].toBool());
  ui->checkBoxSyncZoom->setChecked(map["SyncZoom"].toBool());
  ui->checkBoxCommandKey->setChecked(map["MacUseCommand"].toBool());
  ui->checkBoxMacUnified->setChecked(map["MacUnifiedTitleBar"].toBool());
  ui->radioButtonThemeDark->setChecked(map["DarkConsole"].toBool());
  ui->radioButtonThemeLight->setChecked(!map["DarkConsole"].toBool());
  ui->colorPickerAnnotation->setCurrentColor(map["AnnotationColor"].value<QColor>());
  ui->checkBoxRightButtonErase->setChecked(map["RightButtonErase"].toBool());
  ui->checkBoxAutoReorientView->setChecked(map["AutoReorientView"].toBool());
  ui->checkBoxAutoMidToMin->setChecked(map["AutoSetMidToMin"].toBool());
  ui->checkBoxDecimalVoxelCoord->setChecked(map["DecimalVoxelCoord"].toBool());
  ui->checkBoxAutoScaleFont->setChecked(map["AutoScaleText"].toBool());
  ui->spinBoxFontSize->setValue(map["TextSize"].toInt());
  ui->spinBoxPrecision->setValue(map["Precision"].toInt());
  ui->checkBoxComma->setChecked(map["UseComma"].toBool());
  ui->checkBoxClickToLock->setChecked(map["ClickToLock"].toBool());
  ui->comboBox3DScaleStyle->setCurrentIndex(map["3DAxesFlyMode"].toInt());
  ui->checkBoxAllowDeleteKey->setChecked(map["AllowDeleteKey"].toBool());

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QString val = map.value("ShortcutCycleLayer").toString();
  if (!val.isEmpty() && val != "Default")
  {
    SetActionShortcut(mainwnd->ui->actionCycleLayer, val);
    ui->comboBoxShortcutCycleLayer->setProperty("shortcut_text", val);
    SetCurrentComboText(ui->comboBoxShortcutCycleLayer, val);
  }
  val = map.value("ShortcutToggleVolume").toString();
  if (!val.isEmpty() && val != "Default")
  {
    SetActionShortcut(mainwnd->ui->actionToggleVolumeVisibility, val);
    ui->comboBoxShortcutToggleVolume->setProperty("shortcut_text", val);
    SetCurrentComboText(ui->comboBoxShortcutToggleVolume, val);
  }
  val = map.value("ShortcutToggleSurface").toString();
  if (!val.isEmpty() && val != "Default")
  {
    SetActionShortcut(mainwnd->ui->actionToggleSurfaceVisibility, val);
    ui->comboBoxShortcutToggleSurface->setProperty("shortcut_text", val);
    SetCurrentComboText(ui->comboBoxShortcutToggleSurface, val);
  }
  BlockAllSignals(this, false);
}

QVariantMap DialogPreferences::GetSettings()
{
  QVariantMap map;
  map["BackgroundColor"] = ui->colorPickerBackground->currentColor();
  map["CursorColor"] = ui->colorPickerCursor->currentColor();
  map["CursorSize"] = ui->horizontalSliderSize2D->value();
  map["CursorSize3D"] = ui->horizontalSliderSize3D->value();
  map["CursorThickness"] = ui->horizontalSliderThickness2D->value();
  map["CursorThickness3D"] = ui->horizontalSliderThickness3D->value();
  map["SaveCopy"] = ui->checkBoxSaveCopy->isChecked();
  map["SyncZoom"] = ui->checkBoxSyncZoom->isChecked();
  map["MacUseCommand"] = ui->checkBoxCommandKey->isChecked();
  map["MacUnifiedTitleBar"] = ui->checkBoxMacUnified->isChecked();
  map["DarkConsole"] = ui->radioButtonThemeDark->isChecked();
  map["AnnotationColor"] = ui->colorPickerAnnotation->currentColor();
  map["RightButtonErase"] = ui->checkBoxRightButtonErase->isChecked();
  map["AutoReorientView"] = ui->checkBoxAutoReorientView->isChecked();
  map["AutoSetMidToMin"] = ui->checkBoxAutoMidToMin->isChecked();
  map["DecimalVoxelCoord"] = ui->checkBoxDecimalVoxelCoord->isChecked();
  map["ShortcutCycleLayer"] = ui->comboBoxShortcutCycleLayer->currentText();
  map["ShortcutToggleVolume"] = ui->comboBoxShortcutToggleVolume->currentText();
  map["ShortcutToggleSurface"] = ui->comboBoxShortcutToggleSurface->currentText();
  map["TextSize"] = ui->spinBoxFontSize->value();
  map["AutoScaleText"] = ui->checkBoxAutoScaleFont->isChecked();
  map["Precision"] = ui->spinBoxPrecision->value();
  map["UseComma"] = ui->checkBoxComma->isChecked();
  map["ClickToLock"] = ui->checkBoxClickToLock->isChecked();
  map["3DAxesFlyMode"] = ui->comboBox3DScaleStyle->currentIndex();
  map["AllowDeleteKey"] = ui->checkBoxAllowDeleteKey->isChecked();
  return map;
}

void DialogPreferences::SetCurrentComboText(QComboBox *combo, const QString &text)
{
#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
  combo->setCurrentText(text);
#else
  int n = combo->findText(text);
  if (n < 0)
    n = 0;
  combo->setCurrentIndex(n);
#endif
}

void DialogPreferences::OnClicked(QAbstractButton* btn)
{
  if (ui->buttonBox->buttonRole(btn) == QDialogButtonBox::ResetRole)
  {
    ui->colorPickerBackground->setCurrentColor(Qt::black);
    ui->colorPickerCursor->setCurrentColor(Qt::red);
    ui->colorPickerAnnotation->setCurrentColor(Qt::white);
    ui->horizontalSliderSize2D->setValue(5);
    ui->horizontalSliderSize3D->setValue(1);
    ui->horizontalSliderThickness2D->setValue(1);
    ui->horizontalSliderThickness3D->setValue(1);
    ui->checkBoxSaveCopy->setChecked(true);
    ui->checkBoxRightButtonErase->setChecked(false);
    ui->checkBoxSyncZoom->setChecked(true);
    ui->radioButtonThemeDark->setChecked(true);
    ui->comboBoxShortcutCycleLayer->setCurrentIndex(0);
    ui->comboBoxShortcutToggleVolume->setCurrentIndex(0);
    ui->comboBoxShortcutToggleSurface->setCurrentIndex(0);
#ifdef Q_OS_MAC
    ui->checkBoxCommandKey->setChecked(false);
    ui->checkBoxMacUnified->setChecked(false);
#endif
  }
}

void DialogPreferences::OnComboShortcutChanged(const QString& text)
{
  QStringList list;
  list << ui->comboBoxShortcutCycleLayer->currentText() <<
          ui->comboBoxShortcutToggleSurface->currentText() <<
          ui->comboBoxShortcutToggleVolume->currentText();
  list.removeOne(text);
  if (list.contains(text))
  {
    QMessageBox::warning(this, "Shortcut", "Conflict shortcut. Please select another one.");
    SetCurrentComboText(qobject_cast<QComboBox*>(sender()), sender()->property("shortcut_text").toString());
  }
  else
  {
    sender()->setProperty("shortcut_text", text);
    SetActionShortcut(qobject_cast<QAction*>(sender()->property("action").value<QObject*>()), text);
    MainWindow::GetMainWindow()->UpdateSettings();
  }
}

void DialogPreferences::SetActionShortcut(QAction *act, const QString &text)
{
  QList<QKeySequence> list = act->shortcuts();
  for (int i = 1; i <= 12; i++)
    list.removeAll(QKeySequence(tr("F%1").arg(i)));
  list.removeAll(QKeySequence(tr("Pause")));

  if (!list.contains(QKeySequence(text)))
    list << QKeySequence(text);

  act->setShortcuts(list);
}
