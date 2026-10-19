#include "DialogSegmentationSettings.h"
#include "AppTheme.h"
#include "SegmentedControl.h"

#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFrame>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QProcessEnvironment>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <QStandardPaths>
#include <QStyle>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>

namespace
{
// The settings green, on the pills, Save and focused fields.
const QColor kAccent(0x00, 0x7a, 0x00);
const int kControlHeight = 28;

const QStringList kFoldValues = { "auto", "1", "5" };
const QStringList kDeviceValues = { "auto", "gpu", "cpu" };

SegmentedControl* MakePill(const QStringList& items)
{
  SegmentedControl* s = new SegmentedControl;
  s->SetItems(items);
  s->SetAccent(kAccent);
  s->SetTrackColors(QColor("#dadada"), QColor("#1b1b1b"));
  s->setFixedHeight(kControlHeight + 2);
  s->setCursor(Qt::PointingHandCursor);
  return s;
}
}

bool DialogSegmentationSettings::GpuPresent()
{
#if defined(Q_OS_MACOS) && defined(Q_PROCESSOR_ARM)
  return true; // the daemon uses MPS
#elif defined(Q_OS_LINUX)
  return QFileInfo::exists("/proc/driver/nvidia/version") || !QStandardPaths::findExecutable("nvidia-smi").isEmpty();
#else
  return !QStandardPaths::findExecutable("nvidia-smi").isEmpty();
#endif
}

bool DialogSegmentationSettings::Enabled()
{
  return QSettings().value(KeyEnabled(), true).toBool();
}

QString DialogSegmentationSettings::Folds()
{
  const QString v = QSettings().value(KeyFolds()).toString();
  return kFoldValues.contains(v) ? v : QString("auto");
}

QString DialogSegmentationSettings::Device()
{
  const QString v = QSettings().value(KeyDevice()).toString();
  return kDeviceValues.contains(v) ? v : QString("auto");
}

int DialogSegmentationSettings::MaxThreads()
{
  bool bOk = false;
  const int n = QSettings().value(KeyMaxThreads()).toInt(&bOk);
  return (bOk && n >= 1) ? n : 8;
}

QString DialogSegmentationSettings::ModelDir()
{
  const QString v = QSettings().value(KeyModelDir()).toString();
  return v.isEmpty() ? DefaultModelDir() : v;
}

QString DialogSegmentationSettings::DefaultModelDir()
{
  // NNUNET_MODEL_DIR, else the model installed with FreeSurfer's Python.
  const QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  const QString strEnv = env.value("NNUNET_MODEL_DIR");
  const QString strInstalled = env.value("FREESURFER_HOME_FSPYTHON") + "/models/nnUNetTrainer__nnUNetPlans__2d";
  if ((strEnv.isEmpty() || !QDir(strEnv).exists()) && QDir(strInstalled).exists())
    return strInstalled;
  return strEnv;
}

bool DialogSegmentationSettings::ModelFound(const QString& strDir)
{
  return !strDir.isEmpty() && QFileInfo(QDir(strDir).filePath("plans.json")).isFile();
}

DialogSegmentationSettings::DialogSegmentationSettings(QWidget* parent) : QDialog(parent)
{
  setWindowTitle("Automatic Segmentation Settings");
  setMinimumWidth(680);

  QVBoxLayout* v = new QVBoxLayout(this);
  v->setContentsMargins(28, 24, 28, 20);
  v->setSpacing(6);

  m_pTitle = new QLabel("Automatic segmentation settings");
  m_pSubtitle = new QLabel("Changes apply the next time photos are segmented.");
  m_pSubtitle->setWordWrap(true);
  v->addWidget(m_pTitle);
  v->addWidget(m_pSubtitle);
  v->addSpacing(14);

  QGridLayout* grid = new QGridLayout;
  grid->setHorizontalSpacing(18);
  grid->setVerticalSpacing(4);
  grid->setColumnStretch(1, 1);
  v->addLayout(grid);

  m_pEnabled = MakePill(QStringList() << "On" << "Off");
  m_pEnabled->SetCurrentIndex(Enabled() ? 0 : 1, false);
  m_pEnabledCaption = AddRow(grid, 0, "Automatic segmentation", m_pEnabled, QString());
  connect(m_pEnabled, &SegmentedControl::currentIndexChanged, this, &DialogSegmentationSettings::OnEnabledChanged);

  m_pFolds = MakePill(QStringList() << "Automatic" << "Single fold" << "Five folds");
  m_pFolds->SetCurrentIndex(kFoldValues.indexOf(Folds()), false);
  AddRow(grid, 2, "Model folds", m_pFolds,
         QString::fromUtf8("Automatic: 5 on a GPU, 1 on the CPU. Five is steadier "
                           "but 4\xe2\x80\x93" "5\xc3\x97 slower."));

  m_pDevice = MakePill(QStringList() << "Automatic" << "GPU" << "CPU");
  m_pDevice->SetCurrentIndex(kDeviceValues.indexOf(Device()), false);
  AddRow(grid, 4, "Run on", m_pDevice,
         "Automatic uses the GPU if there is one.");

  {
    QWidget* w = new QWidget;
    QHBoxLayout* hl = new QHBoxLayout(w);
    hl->setContentsMargins(0, 0, 0, 0);
    m_pThreads = new QSpinBox;
    m_pThreads->setObjectName("segSettingsSpin");
    m_pThreads->setRange(1, 256);
    m_pThreads->setValue(MaxThreads());
    m_pThreads->setFixedSize(90, kControlHeight);
    hl->addWidget(m_pThreads);
    hl->addStretch(1);
    m_pThreadsCaption = AddRow(grid, 6, "Maximum threads", w, QString());
    connect(m_pThreads, QOverload<int>::of(&QSpinBox::valueChanged), this, &DialogSegmentationSettings::UpdateThreadsCaption);
  }

  QWidget* pPath = new QWidget;
  QHBoxLayout* h = new QHBoxLayout(pPath);
  h->setContentsMargins(0, 0, 0, 0);
  h->setSpacing(8);
  m_pModelDir = new QLineEdit(ModelDir());
  m_pModelDir->setObjectName("segSettingsPath");
  m_pModelDir->setPlaceholderText("Folder of the nnUNet model");
  m_pModelDir->setFixedHeight(kControlHeight);
  m_pChoose = new QPushButton(QString::fromUtf8("Choose\xe2\x80\xa6"));
  m_pChoose->setObjectName("segSettingsGhost");
  m_pChoose->setCursor(Qt::PointingHandCursor);
  h->addWidget(m_pModelDir, 1);
  h->addWidget(m_pChoose);
  m_pModelCaption = AddRow(grid, 8, "Model folder", pPath, QString());
  connect(m_pChoose, &QPushButton::clicked, this, &DialogSegmentationSettings::OnChooseFolder);
  connect(m_pModelDir, &QLineEdit::textChanged, this, &DialogSegmentationSettings::UpdateModelCaption);
  connect(m_pFolds, &SegmentedControl::currentIndexChanged, this, &DialogSegmentationSettings::UpdateModelCaption);
  connect(m_pDevice, &SegmentedControl::currentIndexChanged, this, &DialogSegmentationSettings::UpdateModelCaption);

  v->addStretch(1);
  v->addSpacing(14);
  // The model's citation, as a footnote above the buttons.
  QFrame* rule = new QFrame;
  rule->setObjectName("segSettingsRule");
  rule->setFixedHeight(1);
  v->addWidget(rule);
  v->addSpacing(6);
  QLabel* pCitation = new QLabel(QString::fromUtf8(
      "Automatic segmentation uses the nnU-Net software package:<br>"
      "Isensee, F., Jaeger, P. F., Kohl, S. A., Petersen, J., &amp; Maier-Hein, K. H. (2021).<br>"
      "nnU-Net: a self-configuring method for deep learning-based biomedical image segmentation.<br>"
      "<i>Nature Methods, 18</i>(2), 203\xe2\x80\x93" "211."));
  pCitation->setObjectName("segSettingsCitation");
  pCitation->setTextFormat(Qt::RichText);
  pCitation->setWordWrap(true);
  pCitation->setTextInteractionFlags(Qt::TextSelectableByMouse);
  v->addWidget(pCitation);
  v->addSpacing(14);
  QHBoxLayout* foot = new QHBoxLayout;
  foot->addStretch(1);
  QPushButton* pCancel = new QPushButton("Cancel");
  pCancel->setObjectName("segSettingsGhost");
  pCancel->setCursor(Qt::PointingHandCursor);
  QPushButton* pSave = new QPushButton("Save");
  pSave->setObjectName("segSettingsPrimary");
  pSave->setCursor(Qt::PointingHandCursor);
  pSave->setDefault(true);
  foot->addWidget(pCancel);
  foot->addWidget(pSave);
  v->addLayout(foot);
  connect(pCancel, &QPushButton::clicked, this, &QDialog::reject);
  connect(pSave, &QPushButton::clicked, this, [this]() { Save(); accept(); });

  Retheme(AppTheme::EffectiveDarkMode());
  UpdateModelCaption();
  UpdateThreadsCaption();
  UpdateEnabled();
  // Start on Save, not with a focus ring on the first control.
  pSave->setFocus();
}

namespace
{
// The warning before turning automatic segmentation off; true to turn it off.
bool ConfirmTurnOff(QWidget* parent)
{
  const bool bDark = AppTheme::EffectiveDarkMode();
  QDialog dlg(parent);
  dlg.setWindowTitle("Turn Off Automatic Segmentation?");
  QVBoxLayout* v = new QVBoxLayout(&dlg);
  v->setContentsMargins(24, 20, 24, 18);
  v->setSpacing(8);
  QLabel* title = new QLabel("Turn off automatic segmentation?");
  title->setStyleSheet(QString("font-weight: bold; font-size: 18px; color: %1;").arg(bDark ? "#f2f2f2" : "#262626"));
  v->addWidget(title);
  QLabel* body = new QLabel("Only turn this off if the automatic segmentation is failing on your photographs. Without it, "
                            "every mask has to be drawn by hand from scratch, which can be very time-consuming.");
  body->setWordWrap(true);
  body->setFixedWidth(420);
  body->setStyleSheet(QString("font-size: 14px; color: %1;").arg(bDark ? "#d8d8d8" : "#3a3a3a"));
  v->addWidget(body);
  v->addSpacing(8);
  QHBoxLayout* h = new QHBoxLayout;
  h->addStretch(1);
  QPushButton* keep = new QPushButton("Keep it on");
  QPushButton* off = new QPushButton("Turn off");
  keep->setStyleSheet("QPushButton {background: #007a00; padding: 0 22px; font-size: 13px; font-weight: bold; border: none;"
                      " border-radius: 14px; color: white; min-height: 28px; max-height: 28px;}"
                      "QPushButton:hover {background: #006b00;}");
  off->setStyleSheet(QString("QPushButton {background: transparent; padding: 0 14px; font-size: 13px; border: 1px solid %1;"
                             " border-radius: 14px; color: %1; min-height: 28px; max-height: 28px;}"
                             "QPushButton:hover {background: %2;}")
                     .arg(bDark ? "#e07b8a" : "#a2293c", bDark ? "#4a3035" : "#f6e3e6"));
  keep->setDefault(true);
  h->addWidget(off);
  h->addWidget(keep);
  v->addLayout(h);
  QObject::connect(keep, &QPushButton::clicked, &dlg, &QDialog::reject);
  QObject::connect(off, &QPushButton::clicked, &dlg, &QDialog::accept);
  return dlg.exec() == QDialog::Accepted;
}
}

void DialogSegmentationSettings::OnEnabledChanged(int nIndex)
{
  if (nIndex == 1 && !ConfirmTurnOff(this))
  {
    m_pEnabled->SetCurrentIndex(0, false);
    return;
  }
  UpdateEnabled();
}

void DialogSegmentationSettings::UpdateEnabled()
{
  const bool bOn = m_pEnabled->CurrentIndex() == 0;
  m_pEnabledCaption->setText(bOn ? QString("Masks are made while photos are corrected; you touch them up.")
                                 : QString("Every mask is drawn by hand from scratch."));
  m_pEnabledCaption->setProperty("level", bOn ? "" : "warning");
  m_pEnabledCaption->style()->unpolish(m_pEnabledCaption);
  m_pEnabledCaption->style()->polish(m_pEnabledCaption);
  if (isVisible())
    QTimer::singleShot(0, this, &DialogSegmentationSettings::FitCaptions);
  // The rest only matters when the model runs. A named list: a range-for over
  // "QList() << a << b" walks a temporary that is already gone.
  const QList<QWidget*> controls = { m_pFolds, m_pDevice, m_pThreads, m_pModelDir, m_pChoose };
  for (QWidget* w : controls)
    w->setEnabled(bOn);
  for (int i = 1; i < m_listLabels.size(); i++)
    m_listLabels[i]->setEnabled(bOn);
  for (int i = 1; i < m_listCaptions.size(); i++)
    m_listCaptions[i]->setEnabled(bOn);
}

void DialogSegmentationSettings::showEvent(QShowEvent* e)
{
  QDialog::showEvent(e);
  QTimer::singleShot(0, this, &DialogSegmentationSettings::FitCaptions);
}

void DialogSegmentationSettings::FitCaptions()
{
  // A wrapped label in a grid row gets the height of its first layout, which
  // clipped the last line of longer captions. Once the widths are known, each
  // gets the height its text needs, and the window grows to fit.
  for (QLabel* c : m_listCaptions)
    if (c->width() > 0)
      c->setMinimumHeight(c->heightForWidth(c->width()));
  adjustSize();
}

void DialogSegmentationSettings::UpdateThreadsCaption()
{
  const int nAvailable = QThread::idealThreadCount();
  const bool bOver = m_pThreads->value() > nAvailable;
  m_pThreadsCaption->setText(QString("This computer has %1; a higher number uses %1.").arg(nAvailable));
  m_pThreadsCaption->setProperty("level", bOver ? "warning" : "");
  m_pThreadsCaption->style()->unpolish(m_pThreadsCaption);
  m_pThreadsCaption->style()->polish(m_pThreadsCaption);
}

QLabel* DialogSegmentationSettings::AddRow(QGridLayout* grid, int nRow, const QString& strLabel, QWidget* pControl,
                                           const QString& strCaption)
{
  QLabel* l = new QLabel(strLabel);
  l->setObjectName("segSettingsLabel");
  l->setFixedHeight(kControlHeight + 2);
  grid->addWidget(l, nRow, 0, Qt::AlignLeft | Qt::AlignTop);
  grid->addWidget(pControl, nRow, 1);
  QLabel* c = new QLabel(strCaption);
  c->setObjectName("segSettingsCaption");
  c->setWordWrap(true);
  grid->addWidget(c, nRow + 1, 1);
  grid->setRowMinimumHeight(nRow + 1, 0);
  grid->addItem(new QSpacerItem(0, 14), nRow + 1, 0);
  m_listLabels << l;
  m_listCaptions << c;
  return c;
}

void DialogSegmentationSettings::OnChooseFolder()
{
  const QString str = QFileDialog::getExistingDirectory(this, "Choose the nnUNet model folder", m_pModelDir->text());
  if (!str.isEmpty())
    m_pModelDir->setText(QDir::toNativeSeparators(str));
}

void DialogSegmentationSettings::UpdateModelCaption()
{
  // What the daemon loads: plans.json, dataset.json and fold_N/ folders.
  const QDir d(m_pModelDir->text().trimmed());
  QString str, level;
  if (m_pModelDir->text().trimmed().isEmpty())
  {
    str = "Choose the model folder (nnUNetTrainer__nnUNetPlans__2d).";
    level = "error";
  }
  else if (!d.exists())
  {
    str = "Folder not found.";
    level = "error";
  }
  else if (!QFileInfo(d.filePath("plans.json")).isFile() || !QFileInfo(d.filePath("dataset.json")).isFile())
  {
    str = "No nnUNet model here (plans.json, dataset.json, fold_0 to fold_4).";
    level = "error";
  }
  else
  {
    int nFolds = 0;
    for (int i = 0; i < 5; i++)
      if (QFileInfo(d.filePath(QString("fold_%1").arg(i))).isDir())
        nFolds++;
    // Automatic means five folds on a graphics card, as in the daemon.
    const QString strFolds = kFoldValues.value(m_pFolds->CurrentIndex());
    const bool bFive = strFolds == "5"
        || (strFolds == "auto" && GpuPresent() && kDeviceValues.value(m_pDevice->CurrentIndex()) != "cpu");
    if (nFolds == 0)
    {
      str = "No folds here (fold_0 to fold_4).";
      level = "error";
    }
    else if (bFive && nFolds < 5)
    {
      str = QString("Only %1 of 5 folds here: five folds will not work.").arg(nFolds);
      level = "warning";
    }
    else
      str = QString("Model found (%1 of 5 folds).").arg(nFolds);
  }
  m_pModelCaption->setText(str);
  m_pModelCaption->setProperty("level", level);
  m_pModelCaption->style()->unpolish(m_pModelCaption);
  m_pModelCaption->style()->polish(m_pModelCaption);
}

void DialogSegmentationSettings::Save()
{
  QSettings s;
  s.setValue(KeyEnabled(), m_pEnabled->CurrentIndex() == 0);
  s.setValue(KeyFolds(), kFoldValues.value(m_pFolds->CurrentIndex(), "auto"));
  // The model folder only when not the default, so the default keeps following
  // the environment.
  s.setValue(KeyDevice(), kDeviceValues.value(m_pDevice->CurrentIndex(), "auto"));
  if (m_pModelDir->text().trimmed() == DefaultModelDir())
    s.remove(KeyModelDir());
  else
    s.setValue(KeyModelDir(), m_pModelDir->text().trimmed());
  s.setValue(KeyMaxThreads(), m_pThreads->value());
}

void DialogSegmentationSettings::Retheme(bool bDark)
{
  const QString strBody    = bDark ? "#f2f2f2" : "#262626";
  const QString strCaption = bDark ? "#9a9a9a" : "#666666";
  const QString strFieldBg = bDark ? "#333333" : "#ffffff";
  const QString strFieldBd = bDark ? "#5f5f5f" : "#c0c0c0";
  const QString strError   = bDark ? "#e07b8a" : "#a2293c";
  const QString strWarning = bDark ? "#e0b252" : "#8a6100";
  const QString strGhostBd = bDark ? "#777777" : "#999999";
  const QString strGhostTx = bDark ? "#dddddd" : "#555555";
  const QString strGhostHv = bDark ? "#4a4a4a" : "#e4e4e4";
  setStyleSheet(QString(
      "QDialog {background: %1;}"
      "QLabel {color: %2; font-size: 13px;}"
      "SegmentedControl {font-size: 13px;}"
      "#segSettingsLabel {font-size: 14px; font-weight: bold;}"
      "#segSettingsCaption {font-size: 12px; color: %3;}"
      "#segSettingsCitation {font-size: 11px; color: %3;}"
      "#segSettingsCaption[level=\"error\"] {color: %4;}"
      "#segSettingsCaption[level=\"warning\"] {color: %5;}"
      "#segSettingsLabel:disabled, #segSettingsCaption:disabled {color: %11;}"
      "#segSettingsPath {background: %6; border: 1px solid %7; border-radius: 4px; padding: 2px 6px; color: %2; font-size: 13px;}"
      "#segSettingsPath:focus {border: 1px solid #007a00;}"
      "#segSettingsSpin {background: %6; border: 1px solid %7; border-radius: 4px; padding: 2px 6px; color: %2; font-size: 13px;}"
      "#segSettingsSpin:focus {border: 1px solid #007a00;}"
      "#segSettingsGhost {background: transparent; padding: 0 14px; font-size: 13px; border: 1px solid %8;"
      "  border-radius: 14px; color: %9; min-height: 28px; max-height: 28px;}"
      "#segSettingsGhost:hover {background: %10;}"
      "#segSettingsPrimary {background: #007a00; padding: 0 22px; font-size: 13px; font-weight: bold; border: none;"
      "  border-radius: 14px; color: white; min-height: 28px; max-height: 28px;}"
      "#segSettingsPrimary:hover {background: #006b00;}"
      "#segSettingsPrimary:pressed {background: #005c00;}")
      .arg(bDark ? "#2d2d2d" : "#f0f0f0", strBody, strCaption, strError, strWarning, strFieldBg, strFieldBd, strGhostBd, strGhostTx)
      .arg(strGhostHv, bDark ? "#6a6a6a" : "#aaaaaa"));
  m_pTitle->setStyleSheet(QString("font-weight: bold; font-size: 20px; color: %1;").arg(strBody));
  m_pSubtitle->setStyleSheet(QString("font-size: 13px; color: %1;").arg(strCaption));
  findChild<QFrame*>("segSettingsRule")->setStyleSheet(QString("#segSettingsRule {background: %1; border: none;}")
                                                       .arg(bDark ? "#444444" : "#d5d5d5"));
}
