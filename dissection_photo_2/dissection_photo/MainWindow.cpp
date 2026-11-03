#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "CommonDef.h"
#include <QFileDialog>
#include <QSettings>
#include <QMessageBox>
#include <QDialogButtonBox>
#include <QStandardPaths>
#include <QBuffer>
#include <QScreen>
#include <QTemporaryDir>
#include <QTimer>
#include <QCryptographicHash>
#include <QTextStream>
#include <QDebug>
#include <QDateTime>
#include <QMimeData>
#include <QDropEvent>
#include <QDragEnterEvent>
#include <QSet>
#include <QApplication>
#include <QPainter>
#include <QIcon>
#include <QImageReader>
#include <QAbstractButton>
#include <QDialog>
#include <QLabel>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include "cnpy.h"
#include "MyUtils.h"
#include "ProgressWindow.h"
#include "SegmentedControl.h"
#include "CompactSlider.h"
#include "ImageGridView.h"
#include "AppMotion.h"
#include <QAbstractAnimation>
#include <QPointer>
#include "StepperBar.h"
#include "ReconPage.h"
#include "DialogSegmentationSettings.h"
#include <QCloseEvent>
#include <QSaveFile>
#include <QJsonArray>
#include <QJsonObject>
#include <QJsonDocument>
#include <QTextEdit>
#include <QPlainTextEdit>
#include <QAbstractSpinBox>
#include <QLineEdit>
#include <QToolButton>

static int ShowBox(QWidget* parent, const QString& strTitle, const QString& strBody,
                   const QString& strCaption, const QStringList& listButtons,
                   int nDanger = -1);
static void Warn(QWidget* parent, const QString& strTitle, const QString& strBody);
static QPixmap MarkerPictogram(int nLogicalSize, bool bTargetCentre,
                               const QColor& target, const QColor& rest);
#include <QScrollArea>
#include <QUrl>
#include "AppTheme.h"
#include <cstdio>

static bool CCSavedIsStale(const QString& strSavedNpy, const QString& strMask);

#define SCRIPT_RETROSPECTIVE "func_retrospective_correction.py"

// True once a pipeline output is finished, not merely present. The Python steps
// write straight onto the final path, so a half-written file exists; the header
// check rejects it cheaply (a full decode would cost ~100 ms per image).
static bool IsCompleteImageOutput(const QString& path)
{
  QFileInfo fi(path);
  if (!fi.exists() || fi.size() <= 0)
    return false;
  QImageReader reader(path);
  if (reader.size().isValid())
    return true;
  // A format this build cannot read (TIFF is optional on some Linux builds):
  // fall back to exists-and-non-empty rather than incomplete forever.
  return reader.error() == QImageReader::UnsupportedFormatError;
}

// For outputs with no header to check (CC <base>_mask.npy files):
// exists-and-non-empty.
static bool IsCompleteDataOutput(const QString& path)
{
  QFileInfo fi(path);
  return fi.exists() && fi.size() > 0;
}

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow), m_nIndex(0)
{
  ui->setupUi(this);
  addAction(ui->actionToggleMask);
  addAction(ui->actionUndoMaskEdit);
  addAction(ui->actionRedoMaskEdit);
  addAction(ui->actionToggleDarkMode);
  SetUpNavigationShortcuts();
  SetUpStepper();
  SetUpReconPage();
  SetUpStartCaption();
  NormaliseSegPanel();
  SetUpCalibrationInstructions();
  // Both pages are taller than an 800px window: scroll the content, keep the
  // footer pinned.
  ScrollPageContent(ui->pageCalibration, ui->pushButtonSaveCalibration);
  ScrollPageContent(ui->pageSegEdit, ui->frameCCProceedDivider);

  GoToStage(ui->pageStartUp);

  m_pImageGridView = new ImageGridView(ui->widgetMain);
  ui->verticalLayout_10->insertWidget(1, m_pImageGridView);
  ui->verticalLayout_10->setStretchFactor(m_pImageGridView, 1);
  m_pImageGridView->hide();
  connect(m_pImageGridView, &ImageGridView::selectionChanged, this, [this](int nSelected, int) {
    if (m_pImageGridView->IsSelectable())
      UpdateStartCaption(nSelected);
  });
  connect(ui->stackedWidget, SIGNAL(currentChanged(int)), SLOT(UpdateStartupPreview()));
  // Connected in code: uic's function-pointer connections are typed against
  // QMainWindow and cannot see this class's slots.
  connect(ui->actionSwitchProject, SIGNAL(triggered()), SLOT(OnActionSwitchProject()));
  // Same reason as actionSwitchProject above.
  connect(ui->pushButtonUndoCC, SIGNAL(clicked()), SLOT(OnButtonUndoCC()));
  connect(ui->pushButtonRedoCC, SIGNAL(clicked()), SLOT(OnButtonRedoCC()));
  // Ctrl + Z / Ctrl + Shift + Z: the mask undo on Segmentation (MainWindow.ui), the
  // box undo here, through the buttons so a greyed-out one does nothing.
  connect(ui->actionUndoMaskEdit, &QAction::triggered, this, [this]()
  {
    if (ui->stackedWidget->currentWidget() == ui->pageCC && ui->pushButtonUndoCC->isEnabled())
      ui->pushButtonUndoCC->click();
  });
  connect(ui->actionRedoMaskEdit, &QAction::triggered, this, [this]()
  {
    if (ui->stackedWidget->currentWidget() == ui->pageCC && ui->pushButtonRedoCC->isEnabled())
      ui->pushButtonRedoCC->click();
  });
  ui->pushButtonUndoCC->setToolTip("Undo the last box (Ctrl+Z)");
  ui->pushButtonRedoCC->setToolTip("Redo (Ctrl+Shift+Z)");
  // The CC page's primary action goes forward to Reconstruction.
  connect(ui->pushButtonQuit, SIGNAL(clicked()), SLOT(OnButtonProceedToRecon()));
  connect(ui->toolButtonDownloadMarkers, SIGNAL(clicked()), SLOT(OnButtonDownloadMarkers()));
  // The line under the printout reads as its caption: left-aligned and centred
  // under the button (sized in FitMarkerCaption()).
  ui->labelDownloadMarkers->setAlignment(Qt::AlignLeft | Qt::AlignTop);
  ui->verticalLayout_12->setAlignment(ui->labelDownloadMarkers, Qt::AlignHCenter);
  // Connect after SetItems()/SetCurrentIndex(): SetItems() resets to index 0, so
  // setting the default would otherwise emit a spurious change.
  ui->segControlPointMode->SetItems(QStringList() << "2 Point" << "3 Point" << "4 Point");
  ui->segControlPointMode->SetCurrentIndex(2, false); // 4 Point default, no animation on startup
  connect(ui->segControlPointMode, SIGNAL(currentIndexChanged(int)), SLOT(OnPointModeChanged(int)));

  // "Show Original | Corrected". "Corrected" is disabled until the photo has a
  // complete corrected output (UpdatePreviewAvailability()). Connect last so the
  // initial state does not trigger a load.
  ui->segControlPreview->SetItems(QStringList() << "Original" << "Corrected");
  ui->segControlPreview->SetItemToolTips(QStringList() << "Show the photo as taken"
                                         << "Show the corrected photo (correct this photo first)");
  ui->segControlPreview->SetItemEnabled(1, false);
  // Left-column track tones, as on the Reconstruction page's pills.
  ui->segControlPreview->SetTrackColors(QColor("#dadada"), QColor("#1b1b1b"));
  ui->segControlPreview->SetCurrentIndex(0, false);
  connect(ui->segControlPreview, &SegmentedControl::currentIndexChanged, this,
          [this](int n) { OnButtonTogglePreview(n == 1); });

  // Sliders: wired in code for the same uic reason as actionSwitchProject.
  // Corrected images are ~700-4000 px, so a 1-200 px brush covers fine and broad
  // strokes.
  ui->sliderBrushSize->SetRange(1, 200);
  ui->sliderBrushSize->SetValue(10);
  ui->sliderBrushSize->SetLabel("Brush size");
  ui->sliderBrushSize->SetSuffix("px");
  ui->sliderBrushSize->SetSwatchColor(QColor(80, 140, 180));
  ui->sliderBrushSize->SetSwatchIcon(":/resource/icons/brush.png");
  connect(ui->sliderBrushSize, SIGNAL(valueChanged(int)), ui->widgetImageView, SLOT(SetBrushSize(int)));
  // Clicking a slider's swatch selects that tool.
  connect(ui->sliderBrushSize, &CompactSlider::swatchClicked, this,
          [this]() { if (m_pToolBrush) m_pToolBrush->setChecked(true); });

  ui->sliderEraserSize->SetRange(1, 200);
  ui->sliderEraserSize->SetValue(10);
  ui->sliderEraserSize->SetLabel("Eraser size");
  ui->sliderEraserSize->SetSuffix("px");
  ui->sliderEraserSize->SetSwatchColor(QColor(180, 100, 80));
  ui->sliderEraserSize->SetSwatchIcon(":/resource/icons/eraser.png");
  // Separate from the brush size; see WidgetImageView::CurrentToolSize().
  connect(ui->sliderEraserSize, SIGNAL(valueChanged(int)), ui->widgetImageView, SLOT(SetEraserSize(int)));
  connect(ui->sliderEraserSize, &CompactSlider::swatchClicked, this,
          [this]() { if (m_pToolEraser) m_pToolEraser->setChecked(true); });

  ui->sliderSegOpacity->SetRange(0, 100);
  ui->sliderSegOpacity->SetValue(70);
  ui->sliderSegOpacity->SetLabel("Mask opacity");
  ui->sliderSegOpacity->SetSuffix("%");
  // No swatch on opacity sliders: they control no tool.
  ui->sliderSegOpacity->SetSwatchVisible(false);
  connect(ui->sliderSegOpacity, SIGNAL(valueChanged(int)), SLOT(OnSliderSegOpacity(int)));

  ui->horizontalSliderLabelOpacity->SetRange(0, 100);
  ui->horizontalSliderLabelOpacity->SetValue(80);
  ui->horizontalSliderLabelOpacity->SetLabel("Slab overlay");
  ui->horizontalSliderLabelOpacity->SetSuffix(" %");
  ui->horizontalSliderLabelOpacity->SetSwatchVisible(false);
  connect(ui->horizontalSliderLabelOpacity, SIGNAL(valueChanged(int)), SLOT(OnSliderLabelOpacity(int)));
  OnSliderLabelOpacity(ui->horizontalSliderLabelOpacity->Value());

  m_strPythonCmd = PY_COMMAND;

  SetupScriptPath();

  QSettings s;
  QRect rc = s.value("MainWindow/Geometry").toRect();
  if (rc.isValid() && QGuiApplication::primaryScreen()->geometry().contains(rc))
    setGeometry(rc);
  ui->splitter->restoreState(s.value("MainWindow/SplitterState").toByteArray());
  // Window resizes go to the image pane; only the splitter handle changes the
  // control column's width, so its labels don't reflow.
  ui->splitter->setStretchFactor(0, 0);
  ui->splitter->setStretchFactor(1, 1);

  SetDarkMode(AppTheme::EffectiveDarkMode());
  // Icons after SetDarkMode(), so later theme rules cannot land on top of them.
  ApplyButtonIcons();
  // Size the counter for its widest text once (after SetDarkMode(), which sets
  // the font), so Previous/Next never shift.
  for (QLabel* l : {ui->labelIndex, ui->labelIndexSeg, ui->labelIndexCC})
  {
    l->ensurePolished();
    l->setMinimumWidth(l->fontMetrics().horizontalAdvance("999 of 999") + 8);
  }

  // Ctrl+Shift+D is this window's only light/dark toggle; the Hub has the visible
  // control.
  connect(ui->actionToggleDarkMode, &QAction::triggered, this,
          [this]() { OnToggleDarkMode(!m_bDarkMode); });

  // AnimateHeight() clears the minimum height when done, so remember the .ui
  // floor.
  m_nCalibrationBoxFloor = ui->widgetDropCalibration->minimumHeight();
  connect(ui->toolButtonToggleCalibration, &QAbstractButton::toggled,
          this, &MainWindow::OnToggleCalibrationDisclosure);

  ui->widgetImageView->SetNumberOfExpectedPoints(4);
  connect(ui->widgetImageView, SIGNAL(PointEdited()), SLOT(OnPointEdited()), Qt::QueuedConnection);
  // The calibration marker count follows the same signal, so it is right
  // whether a point was added, moved or cleared.
  connect(ui->widgetImageView, SIGNAL(PointEdited()), SLOT(UpdateCalibrationProgress()), Qt::QueuedConnection);
  // Said on the photo, like Connected components' "Mask created and saved";
  // switching to Original reloads the view, which clears it.
  connect(ui->widgetImageView, &WidgetImageView::PointRejected, this, [this]()
  {
    ui->widgetImageView->ShowMessage("To redo the correction, toggle to the original photograph.");
  });
  // A mask drawn for a photo that had none can be the last one Connected components waits for.
  connect(ui->widgetImageView, &WidgetImageView::MaskEdited, this, [this]()
  {
    if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
      UpdateSegReadiness();
  }, Qt::QueuedConnection);

  m_wndProgress = new ProgressWindow(this);
  m_wndProgress->hide();

  qRegisterMetaType<QProcess::ExitStatus>("QProcess::ExitStatus");
  m_proc = new QProcess(this);
  connect(m_proc, SIGNAL(readyReadStandardOutput()), SLOT(OnProcessOutputMessage()));
  connect(m_proc, SIGNAL(readyReadStandardError()), SLOT(OnProcessErrorMessage()));
  connect(m_proc, SIGNAL(started()), SLOT(OnProcessStarted()));
  connect(m_proc, SIGNAL(finished(int, QProcess::ExitStatus)), SLOT(OnProcessFinished()), Qt::QueuedConnection);
  connect(m_proc, SIGNAL(errorOccurred(QProcess::ProcessError)), SLOT(OnProcessError(QProcess::ProcessError)));

  connect(&m_fileWatcher, SIGNAL(directoryChanged(QString)), SLOT(OnFileChanged(QString)), Qt::QueuedConnection);

  m_procSeg = new QProcess(this);
  m_procSeg->setStandardOutputFile(QProcess::nullDevice());
  connect(m_procSeg, &QProcess::readyReadStandardError, this, [this]()
  {
    const QString str = m_procSeg->readAllStandardError().trimmed();
    qDebug() << str;
    m_strNNUnetStderr = (m_strNNUnetStderr + str + "\n").right(20000);
  });
  // Queued, as for m_proc: the handler may show a dialog. A crash also emits
  // finished(); a failed start emits only errorOccurred().
  connect(m_procSeg, static_cast<void (QProcess::*)(int, QProcess::ExitStatus)>(&QProcess::finished), this,
          [this](int nCode, QProcess::ExitStatus status) { OnSegDaemonFinished(status, nCode); },
          Qt::QueuedConnection);
  connect(m_procSeg, &QProcess::errorOccurred, this, [this](QProcess::ProcessError er)
  {
    if (er == QProcess::FailedToStart)
    {
      m_strNNUnetStderr += m_procSeg->errorString();
      OnSegDaemonFinished(QProcess::CrashExit, -1);
    }
  }, Qt::QueuedConnection);

  ui->widgetDropCalibration->setAcceptDrops(true);
  ui->widgetDropCalibration->setMouseTracking(true);
  ui->widgetDropImages->setAcceptDrops(true);
  ui->widgetDropCalibrationPhoto->setAcceptDrops(true);
  ui->widgetDropCalibration->installEventFilter(this);
  qApp->installEventFilter(this);
  ui->widgetDropImages->installEventFilter(this);
  ui->widgetDropCalibrationPhoto->installEventFilter(this);
  connect(ui->listWidgetThumbnail, SIGNAL(CurrentIndexChanged(int)), SLOT(OnThumbnailIndexChanged(int)));

  // Floor the left column (plus its 32px padding) at the widest button it
  // holds, measured since font substitution changes widths. QAbstractButton
  // includes toolButtonToggleCalibration.
  int nWidestButton = 0;
  foreach (QAbstractButton* btn, ui->widgetControl->findChildren<QAbstractButton*>())
    nWidestButton = qMax(nWidestButton, btn->sizeHint().width());
  if (nWidestButton + 32 > ui->widgetControl->minimumWidth())
    ui->widgetControl->setMinimumWidth(nWidestButton + 32);

  EqualizePrimaryActionButtonWidths();
}

void MainWindow::EqualizePrimaryActionButtonWidths()
{
  // The Continue buttons of Correction, Segmentation and CC (pushButtonQuit).
  // pushButtonProcess and pushButtonCreateMask are excluded: they share a row
  // and must fit the 380px column. Hints are valid: setupUi() applied the sheets.
  const QList<QAbstractButton*> listButtons = {
    ui->pushButtonSegmentation, ui->pushButtonCC, ui->pushButtonQuit
  };

  int nWidest = 0;
  foreach (QAbstractButton* pBtn, listButtons)
    nWidest = qMax(nWidest, pBtn->sizeHint().width());
  foreach (QAbstractButton* pBtn, listButtons)
  {
    if (pBtn->minimumWidth() < nWidest)
      pBtn->setMinimumWidth(nWidest);
  }
}

MainWindow::~MainWindow()
{
  QSettings s;
  s.setValue("MainWindow/Geometry", geometry());
  s.setValue("MainWindow/SplitterState",  ui->splitter->saveState());

  // Kill a running step on close, with a bounded wait, so no partially written
  // output is left looking finished. Signals blocked: no error dialogs now.
  if (m_proc && m_proc->state() != QProcess::NotRunning)
  {
    m_proc->blockSignals(true);
    m_proc->kill();
    m_proc->waitForFinished(3000);
  }
  if (m_procCCWarmup && m_procCCWarmup->state() != QProcess::NotRunning)
  {
    m_procCCWarmup->blockSignals(true);
    m_procCCWarmup->kill();
    m_procCCWarmup->waitForFinished(1000);
  }
  StopSegDaemon();

  delete ui;
}

bool MainWindow::eventFilter(QObject *watched, QEvent *event)
{
  // Installed on the application: a click anywhere in the window outside a
  // text box ends typing in it (clicks on a button already did).
  if (event->type() == QEvent::MouseButtonPress && watched->isWidgetType())
  {
    QWidget* pFocus = QApplication::focusWidget();
    QWidget* w = static_cast<QWidget*>(watched);
    if (pFocus && pFocus->window() == this && w->window() == this && w != pFocus && !pFocus->isAncestorOf(w)
        && (qobject_cast<QLineEdit*>(pFocus) || qobject_cast<QAbstractSpinBox*>(pFocus)))
      pFocus->clearFocus();
  }
  if (watched == ui->widgetDropImages || watched == ui->widgetDropCalibration || watched == ui->widgetDropCalibrationPhoto)
  {
    if (event->type() == QEvent::DragEnter)
    {
      QDragEnterEvent *e = static_cast<QDragEnterEvent *>(event);
      QList<QUrl> urls = e->mimeData()->urls();
      if (!urls.isEmpty() && urls.last().scheme().toLower() == "file")
        e->acceptProposedAction();
    }
    else if (event->type() == QEvent::Drop)
    {
      QDropEvent *e = static_cast<QDropEvent *>(event);
      QList<QUrl> urls = e->mimeData()->urls();
      if (!urls.isEmpty() && urls.last().scheme().toLower() != "file")
        Warn(this, "Nothing to add", "Drop a folder of photos, or the photo files themselves.");
      else if (!urls.isEmpty())
      {
        e->acceptProposedAction();
        QFileInfo fi(urls.last().toLocalFile());
        if (watched == ui->widgetDropImages)
          SetInputFolder(fi.isDir()?fi.absoluteFilePath():fi.absolutePath());
        else if (watched == ui->widgetDropCalibration)
          SetCalibrationFile(fi.absoluteFilePath());
        else if (watched == ui->widgetDropCalibrationPhoto)
          SetCalibrationPhoto(fi.absoluteFilePath());
      }
    }
  }

  return QMainWindow::eventFilter(watched, event);
}

void MainWindow::resizeEvent(QResizeEvent* e)
{
  RepositionProgressWindow();
  QMainWindow::resizeEvent(e);
}

void MainWindow::showEvent(QShowEvent* e)
{
  QMainWindow::showEvent(e);
  // A project loaded before the first show() opens the calibration box at zero
  // width; re-measure its minimum once laid out.
  QTimer::singleShot(0, this, [this]()
  {
    if (ui->toolButtonToggleCalibration->isChecked())
      OnToggleCalibrationDisclosure(true);
  });
}

// The top and bottom bars follow the page:
//
//   page             top bar   bottom bar (carousel)
//   Start            hidden    hidden
//   Calibration      hidden*   hidden
//   Correction       shown     shown
//   Segmentation     hidden    shown
//   Connected comp.  hidden    shown
//   Reconstruction   hidden    hidden
//
// * OnButtonCalibration() shows it once a photo is loaded.
// The step bar spans the whole window: centralwidget's layout moves onto a body
// widget under a new vertical [strip, body] layout, so the .ui stays as Designer
// writes it.
void MainWindow::SetUpStepper()
{
  m_pStepper = new StepperBar(ui->centralwidget);
  m_pStepper->SetSteps(QStringList() << "Add photos" << "Correction" << "Segmentation"
                                     << "Connected components" << "Reconstruction");
  // Home, at the bar's left end: back to the Projects Hub, as File > Switch Project.
  m_pHomeButton = new QToolButton;
  m_pHomeButton->setObjectName("stepperHome");
  m_pHomeButton->setToolTip("Back to Projects");
  m_pHomeButton->setAccessibleName("Back to Projects");
  m_pHomeButton->setCursor(Qt::PointingHandCursor);
  m_pHomeButton->setIconSize(QSize(18, 18));
  m_pHomeButton->setFixedSize(34, 24);
  connect(m_pHomeButton, &QToolButton::clicked, this, &MainWindow::OnActionSwitchProject);
  m_pStepper->SetLeadingWidget(m_pHomeButton);

  QWidget* pBody = new QWidget(ui->centralwidget);
  pBody->setLayout(ui->centralwidget->layout());   // reparents the splitter's layout

  QVBoxLayout* pOuter = new QVBoxLayout(ui->centralwidget);
  pOuter->setContentsMargins(0, 0, 0, 0);
  pOuter->setSpacing(0);
  pOuter->addWidget(m_pStepper);
  pOuter->addWidget(pBody, 1);

  // A finished step is reopened through the same slots its Back button uses, not
  // GoToStage(), which lacks the per-page setup (edit mode, accent, index, masks).
  connect(m_pStepper, &StepperBar::StepActivated, this, [this](int n)
  {
    // The step cannot be left while a reconstruction runs.
    if (m_pReconPage && m_pReconPage->IsRunning())
      return;
    switch (n)
    {
    // Always backward to a finished step: bBack skips forward-only work.
    case 0: OnButtonBackToStart(); break;
    case 1: OnButtonContinue(true); break;
    case 2: OnButtonProceedToSeg(true); break;
    case 3: OnButtonProceedToCC(true); break;
    default: break;
    }
  });
}

// Steps in order. Calibration is a detour from Add photos and reports as step 0.
// A step is Done when its work is finished, not when the user has passed it.
// Done steps are clickable, so an unfinished step gets no tick and no click; its
// page's Back button still reaches it.
void MainWindow::UpdateStepper(QWidget* pPage)
{
  if (!m_pStepper)
    return;
  int nCurrent = 0;
  if (pPage == ui->pageCorrection)         nCurrent = 1;
  else if (pPage == ui->pageSegEdit)       nCurrent = 2;
  else if (pPage == ui->pageCC)            nCurrent = 3;
  else if (pPage == ui->pagePhotoRecon)    nCurrent = 4;

  bool bDone[5] = { false, false, false, false, false };
  bDone[0] = !m_listInputFiles.isEmpty();
  // Re-list the original photos: from Segmentation on, m_listInputFiles holds the
  // corrected images. Only when a tick can be shown (this reads every header).
  if (nCurrent > 1)
  {
    // Guarded: QDir("") is the working directory, not "no folder".
    const QFileInfoList photos = m_strInputFolder.isEmpty() ? QFileInfoList()
        : QDir(m_strInputFolder).entryInfoList(ImageNameFilters(), QDir::Files, QDir::Name);
    bDone[1] = AllImagesProcessed(photos);
    // Segmentation is done when every photo has a mask on disk. Counted from the
    // folder: masks can still arrive after leaving the page.
    bDone[2] = !photos.isEmpty() && !m_strMaskFolder.isEmpty()
               && QDir(m_strMaskFolder).entryInfoList(QDir::Files).size() == photos.size();
  }
  // The repointed list is right here: CC files are named after the corrected
  // image (<base>_deformed_mask.npy).
  bDone[3] = CCReadyForRecon();

  QVector<StepperBar::State> states;
  for (int i = 0; i < 5; i++)
  {
    if (i == nCurrent)
      states << StepperBar::Current;
    else if (i < nCurrent && bDone[i])
      states << StepperBar::Done;
    else
      states << StepperBar::Future;
  }
  m_pStepper->SetStates(states);
}

// ReconPage runs mri_3d_photo_recon. Its command pane goes in widgetMain beside
// the image view and grid; UpdateStartupPreview() picks which is shown.
void MainWindow::SetUpReconPage()
{
  m_pReconPage = new ReconPage(ui->pagePhotoRecon);
  ui->label_6->setText("3D reconstruction");
  ui->verticalLayout_3->insertWidget(1, m_pReconPage, 1);

  // Index 2 is just before frameBottomBar; the Start page's grid is inserted
  // at 1 later, which keeps this ahead of the bottom bar.
  ui->verticalLayout_10->insertWidget(2, m_pReconPage->CommandPane(), 1);
  m_pReconPage->CommandPane()->hide();

  connect(m_pReconPage, &ReconPage::BackRequested, this, [this]() { OnButtonProceedToCC(true); });
}

void MainWindow::UpdateReconPage()
{
  if (!m_pReconPage)
    return;
  m_pReconPage->SetProject(m_strRootOutputFolder);
}

void MainWindow::SetUpStartCaption()
{
  // Insert into the page's vertical layout, on the row above the button's row:
  // inside the button's own row the stretch squeezes it to nothing.
  QBoxLayout* pPage = qobject_cast<QBoxLayout*>(ui->pageStartUp->layout());
  if (!pPage)
    return;
  int nRow = -1;
  for (int i = 0; i < pPage->count() && nRow < 0; i++)
  {
    QLayoutItem* it = pPage->itemAt(i);
    if (it->widget() && it->widget()->isAncestorOf(ui->pushButtonGoNew))
      nRow = i;
    else if (it->layout() && it->layout()->indexOf(ui->pushButtonGoNew) >= 0)
      nRow = i;
  }
  if (nRow < 0)
    return;
  m_pStartCaption = new QLabel(ui->pageStartUp);
  m_pStartCaption->setObjectName("startCaption");
  m_pStartCaption->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  pPage->insertWidget(nRow, m_pStartCaption);

  UpdateStartCaption(0);
}

void MainWindow::UpdateStartCaption(int nPhotos)
{
  if (!m_pStartCaption)
    return;
  QString str;
  if (nPhotos <= 0)
    str = "No photographs chosen yet";
  else if (nPhotos == 1)
    str = "1 photograph ready";
  else
    str = QString("%1 photographs ready").arg(nPhotos);
  if (!m_strCalibrationFile.isEmpty())
    str += QString::fromUtf8(" \xc2\xb7 calibrated");
  m_pStartCaption->setText(str);
}

// Segmentation's shortcuts card moves to the bottom, above the divider, so the
// page reads: title, description, Back, controls, card.
void MainWindow::NormaliseSegPanel()
{
  QBoxLayout* pPage = qobject_cast<QBoxLayout*>(ui->pageSegEdit->layout());
  if (!pPage)
    return;

  const int nHotkeys = pPage->indexOf(ui->labelSegHotkeys);
  const int nDivider = pPage->indexOf(ui->frameCCProceedDivider);
  if (nHotkeys < 0 || nDivider < 0)
    return;

  // Tool row: a segmented control that shows which tool is active.
  {
    QBoxLayout* pCtrls = qobject_cast<QBoxLayout*>(ui->widgetSegCtrls->layout());
    if (pCtrls)
    {
      QWidget* pRow = new QWidget(ui->widgetSegCtrls);
      QHBoxLayout* h = new QHBoxLayout(pRow);
      h->setContentsMargins(0, 0, 0, 0);
      h->setSpacing(6);
      h->addStretch(1);
      struct Tool { const char* icon; const char* text; const char* tip; QToolButton** target; };
      const Tool tools[] = {
        { ":/resource/icons/brush.png",  "Brush",  "Brush (Q). Ctrl + drag on the photo to paint.",      &m_pToolBrush  },
        { ":/resource/icons/eraser.png", "Eraser", "Eraser (W). Ctrl + drag on the photo to erase.",     &m_pToolEraser },
        { ":/resource/icons/bucket.png", "Fill",   "Fill (E). Ctrl + click on the photo to flood fill.", &m_pToolFill   },
      };
      for (int i = 0; i < 3; i++)
      {
        QToolButton* b = new QToolButton(pRow);
        b->setObjectName("segTool");
        b->setText(tools[i].text);
        b->setToolTip(tools[i].tip);
        b->setIcon(QIcon(QString::fromLatin1(tools[i].icon)));
        b->setIconSize(QSize(16, 16));
        b->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
        b->setCheckable(true);
        b->setAutoExclusive(true);
        b->setCursor(Qt::PointingHandCursor);
        // Function-pointer connect: SetMaskToolFromRow() is not a slot. Only act on the
        // button turning on; exclusive buttons also emit toggled(false).
        connect(b, &QToolButton::toggled, this, [this](bool bOn) { if (bOn) SetMaskToolFromRow(); });
        h->addWidget(b);
        *tools[i].target = b;
      }
      h->addStretch(1);
      m_pToolBrush->setChecked(true);
      pCtrls->insertWidget(0, pRow);

      // Q / W / E check the tool buttons, so the row, canvas and cursor badge always
      // agree. Window-level actions so they work with focus on the photo; only on
      // Segmentation. A focused text field keeps its keys.
      struct Key { const char* seq; QToolButton* btn; };
      const Key keys[] = { { "Q", m_pToolBrush }, { "W", m_pToolEraser }, { "E", m_pToolFill } };
      for (int i = 0; i < 3; i++)
      {
        QAction* a = new QAction(this);
        a->setShortcut(QKeySequence(QString::fromLatin1(keys[i].seq)));
        QToolButton* btn = keys[i].btn;
        connect(a, &QAction::triggered, this, [this, btn]()
        {
          if (ui->stackedWidget->currentWidget() == ui->pageSegEdit && btn->isEnabled())
            btn->setChecked(true);
        });
        addAction(a);
      }
    }
  }

  // Names the state rather than the action.
  ui->pushButtonToggleMaskVisibility->setText("Show mask");

  // What Connected components waits for, under the description.
  {
    const int nDesc = pPage->indexOf(ui->label_10);
    if (nDesc >= 0)
    {
      m_pSegPreviewNote = new QLabel(ui->pageSegEdit);
      m_pSegPreviewNote->setObjectName("segPreviewNote");
      m_pSegPreviewNote->setTextFormat(Qt::RichText);
      m_pSegPreviewNote->setWordWrap(true);
      pPage->insertWidget(nDesc + 1, m_pSegPreviewNote);
    }
  }

  pPage->removeWidget(ui->labelSegHotkeys);
  pPage->insertWidget(pPage->indexOf(ui->frameCCProceedDivider), ui->labelSegHotkeys);

  // Readiness line above the primary action, as on the Start page.
  int nBtnRow = -1;
  for (int i = 0; i < pPage->count() && nBtnRow < 0; i++)
  {
    QLayoutItem* it = pPage->itemAt(i);
    if (it->widget() && it->widget()->isAncestorOf(ui->pushButtonCC))
      nBtnRow = i;
    else if (it->layout() && it->layout()->indexOf(ui->pushButtonCC) >= 0)
      nBtnRow = i;
    else if (it->widget() == ui->pushButtonCC)
      nBtnRow = i;
  }
  if (nBtnRow >= 0)
  {
    m_pSegContinueCaption = new QLabel(ui->pageSegEdit);
    m_pSegContinueCaption->setObjectName("segContinueCaption");
    m_pSegContinueCaption->setAlignment(Qt::AlignCenter);
    pPage->insertWidget(nBtnRow, m_pSegContinueCaption);
  }
  UpdateSegContinueCaption();
}

// The tool buttons are auto-exclusive, so exactly one is checked.
void MainWindow::SetMaskToolFromRow()
{
  if (!m_pToolBrush)
    return;
  int nTool = WidgetImageView::MT_BRUSH;
  if (m_pToolEraser->isChecked())
    nTool = WidgetImageView::MT_ERASER;
  else if (m_pToolFill->isChecked())
    nTool = WidgetImageView::MT_BUCKET;
  ui->widgetImageView->SetMaskTool(nTool);
}

void MainWindow::UpdateSegContinueCaption()
{
  if (!m_pSegContinueCaption)
    return;
  const int nTotal = m_listInputFiles.size();
  const int nReady = m_listMaskFiles.size();
  if (nTotal <= 0)
    m_pSegContinueCaption->setText(QString());
  else if (nReady < nTotal)
    m_pSegContinueCaption->setText(QString("%1 of %2 masks ready").arg(nReady).arg(nTotal));
  else
    m_pSegContinueCaption->setText(QString("All %1 masks ready").arg(nTotal));
}

// Caption each pictogram with the click it stands for, and count placed
// markers.
void MainWindow::SetUpCalibrationInstructions()
{
  // horizontalLayout_5 is nested, so find the row that holds label_7 before
  // building anything.
  QBoxLayout* pRow = NULL;
  foreach (QBoxLayout* l, ui->widgetCalibrationControl->findChildren<QBoxLayout*>())
  {
    if (l->indexOf(ui->label_7) >= 0 && l->indexOf(ui->label_8) >= 0)
    {
      pRow = l;
      break;
    }
  }
  if (!pRow)
    return;

  m_pCalibCap1 = new QLabel(QString::fromUtf8("1 \xc2\xb7 Centre"));
  m_pCalibCap2 = new QLabel(QString::fromUtf8("2 \xc2\xb7 Outer ring"));
  m_pCalibCap1->setObjectName("calibCaption");
  m_pCalibCap2->setObjectName("calibCaption");
  // Insert the later one first: inserting at n7 would otherwise shift n8.
  pRow->insertWidget(pRow->indexOf(ui->label_8) + 1, m_pCalibCap2);
  pRow->insertWidget(pRow->indexOf(ui->label_7) + 1, m_pCalibCap1);

  // The count goes above the pictograms, under the instruction sentence.
  QBoxLayout* pOuter = qobject_cast<QBoxLayout*>(ui->widgetCalibrationControl->layout());
  if (pOuter)
  {
    const int nInstr = pOuter->indexOf(ui->label_18);
    if (nInstr >= 0)
    {
      m_pCalibProgress = new QLabel();
      m_pCalibProgress->setObjectName("calibProgress");
      pOuter->insertWidget(nInstr + 1, m_pCalibProgress);
    }
  }
  UpdateCalibrationProgress();
}

void MainWindow::UpdateCalibrationProgress()
{
  if (!m_pCalibProgress)
    return;
  // Two clicks per marker: centre, then a point on the outer ring.
  const int nPoints = ui->widgetImageView->GetEditedPoints().size();
  m_pCalibProgress->setText(QString("%1 of 4 markers placed").arg(qMin(4, nPoints / 2)));
}

static bool LayoutContains(QLayout* pLayout, QWidget* pWidget)
{
  if (!pLayout)
    return false;
  for (int i = 0; i < pLayout->count(); i++)
  {
    QLayoutItem* it = pLayout->itemAt(i);
    if (it->widget() == pWidget)
      return true;
    if (it->layout() && LayoutContains(it->layout(), pWidget))
      return true;
  }
  return false;
}

// Nested layouts are wrapped in a widget and setLayout() on it: addLayout()
// does not reparent their widgets, which then draw at stale positions.
void MainWindow::ScrollPageContent(QWidget* pPage, QWidget* pFooter)
{
  QBoxLayout* pCol = qobject_cast<QBoxLayout*>(pPage->layout());
  if (!pCol || !pFooter)
    return;

  // indexOf() only sees direct children; the footer can be nested one level down.
  int nFooter = -1;
  for (int i = 0; i < pCol->count() && nFooter < 0; i++)
  {
    QLayoutItem* it = pCol->itemAt(i);
    if (it->widget() == pFooter ||
        (it->widget() && it->widget()->isAncestorOf(pFooter)) ||
        (it->layout() && LayoutContains(it->layout(), pFooter)))
      nFooter = i;
  }
  if (nFooter <= 0)
    return;

  QWidget* pContent = new QWidget;
  QVBoxLayout* pInner = new QVBoxLayout(pContent);
  pInner->setContentsMargins(0, 0, 0, 0);
  pInner->setSpacing(pCol->spacing());
  // Always take the front item, which preserves the order.
  for (int i = 0; i < nFooter; i++)
  {
    QLayoutItem* it = pCol->takeAt(0);
    if (it->widget())
    {
      pInner->addWidget(it->widget());   // addWidget reparents; this is safe
    }
    else if (it->layout())
    {
      QWidget* pWrap = new QWidget;
      // Keep the child layout's own margins; as a top-level layout it would get the
      // style default.
      const QMargins marg = it->layout()->contentsMargins();
      pWrap->setLayout(it->layout());    // reparents the layout AND its widgets
      it->layout()->setContentsMargins(marg);
      pInner->addWidget(pWrap);
    }
    else if (it->spacerItem())
    {
      pInner->addSpacing(qMax(0, it->spacerItem()->sizeHint().height()));
    }
    // Free widget/spacer wrapper items only. For a nested layout the item is the
    // layout itself, now owned by pWrap.
    if (!it->layout())
      delete it;
  }
  pInner->addStretch(1);

  QScrollArea* pScroll = new QScrollArea(pPage);
  pScroll->setObjectName("panelScroll");
  pScroll->setWidgetResizable(true);
  pScroll->setFrameShape(QFrame::NoFrame);
  pScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  pScroll->setWidget(pContent);
  // style.qss fills every QWidget, the scroll area's viewport included, which
  // would paint over the panel colour behind it.
  pScroll->setStyleSheet("#panelScroll, #panelScroll > QWidget > QWidget {background: transparent;}");
  pCol->insertWidget(0, pScroll, 1);
}

void MainWindow::GoToStage(QWidget* pPage)
{
  const bool bCanvas = (pPage == ui->pageCorrection ||
                        pPage == ui->pageSegEdit ||
                        pPage == ui->pageCC);
  // Set the bars before switching so the crossfade snapshot shows the outgoing
  // page as it was.
  ui->frameTopBar->setVisible(pPage == ui->pageCorrection);
  ui->frameBottomBar->setVisible(bCanvas);
  // Calibration's recentre button takes its olive accent here; the other pages
  // set it in their handlers.
  if (pPage == ui->pageCalibration)
    SetRecenterAccent("#a6a233", "#928e2c", "#7e7b26");
  AppMotion::CrossfadeStackAfterWork(ui->stackedWidget, pPage);
  UpdateStepper(pPage);
  if (pPage == ui->pagePhotoRecon)
    UpdateReconPage();
  if (pPage == ui->pageSegEdit)
    UpdateSegContinueCaption();
}

void MainWindow::SetRecenterAccent(const QString& strBase, const QString& strHover,
                                   const QString& strPressed)
{
  // The .ui rule with only the colours substituted.
  ui->toolButton->setStyleSheet(
        QString("QToolButton {background: %1; padding: 6 8 6 8;font-size: 14px;"
                " border: none; border-radius: 16; color:white}"
                "QToolButton:disabled {background: #999999;}"
                "QToolButton:hover {background: %2;}"
                "QToolButton:pressed {background: %3;}")
        .arg(strBase).arg(strHover).arg(strPressed));
}

void MainWindow::ApplyButtonIcons()
{
  // These sit on a filled accent with white text, so white glyphs work in both
  // themes.
  const QIcon iconUndo = AppTheme::TintedIcon("undo", Qt::white);
  const QIcon iconRedo = AppTheme::TintedIcon("undo", Qt::white, true);
  QPushButton* undos[] = {ui->pushButtonMaskEditUndo, ui->pushButtonUndoCC,
                          ui->pushButtonUndoCalibration};
  QPushButton* redos[] = {ui->pushButtonMaskEditRedo, ui->pushButtonRedoCC,
                          ui->pushButtonRedoCallibration};
  for (int i = 0; i < 3; i++)
  {
    undos[i]->setIcon(iconUndo);
    undos[i]->setIconSize(QSize(16, 16));
    redos[i]->setIcon(iconRedo);
    redos[i]->setIconSize(QSize(16, 16));
  }

  ui->toolButton->setIcon(AppTheme::TintedIcon("recenter", Qt::white));
  // QToolButton defaults to icon-only; keep the "Re-Center" text beside the icon.
  ui->toolButton->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);

  // Toggle Mask's text colour flips when checked, so give the icon On and Off
  // tints.
  QIcon iconMask;
  const QIcon iconMaskOff = AppTheme::TintedIcon("hide", QColor(0x55, 0x55, 0x55));
  const QIcon iconMaskOn  = AppTheme::TintedIcon("hide", Qt::white);
  if (!iconMaskOff.isNull() && !iconMaskOn.isNull())
  {
    iconMask.addPixmap(iconMaskOff.pixmap(24, 24), QIcon::Normal, QIcon::Off);
    iconMask.addPixmap(iconMaskOn.pixmap(24, 24),  QIcon::Normal, QIcon::On);
    ui->pushButtonToggleMaskVisibility->setIcon(iconMask);
    ui->pushButtonToggleMaskVisibility->setIconSize(QSize(16, 16));
  }
}

// Previous/Next and Undo/Redo are semibold. Appended after
// RethemeGhostPillButtons(), which rebuilds some of these sheets; the
// font-weight check avoids stacking rules. 600 because DejaVu Sans has no 500.
static void EmboldenNavButtons(std::initializer_list<QPushButton*> buttons)
{
  for (QPushButton* b : buttons)
  {
    const QString strCss = b->styleSheet();
    if (!strCss.contains("font-weight"))
      b->setStyleSheet(strCss + "\nQPushButton {font-weight: 600;}");
  }
}

void MainWindow::RethemeGhostPillButtons(bool bDark)
{
  // Dark tones match ImageGridView's lightbox pills.
  const QString strBorder  = bDark ? "#777777" : "#999999";
  const QString strText    = bDark ? "#dddddd" : "#555555";
  const QString strHoverBg = bDark ? "#4a4a4a" : "#f0f0f0";
  const QString strHoverBd = bDark ? "#999999" : "#888888";
  const QString strHoverTx = bDark ? "#ffffff" : "#333333";
  const QString strPressBg = bDark ? "#2e2e2e" : "#e2e2e2";
  const QString strPressBd = bDark ? "#aaaaaa" : "#777777";
  const QString strPressTx = bDark ? "#ffffff" : "#333333";
  const QString strDisTx   = bDark ? "#6a6a6a" : "#bbbbbb";
  const QString strDisBd   = bDark ? "#3f3f3f" : "#dddddd";

  // Padding, font size and radius copied from each button's .ui rule, so only
  // colours change.
  auto Pill = [&](const QString& strPadding, int nFont, int nRadius) {
    return QString(
          "QPushButton {background: transparent; padding: %1; font-size: %2px;"
          " border: 1px solid %4; border-radius: %3; color:%5}"
          "QPushButton:disabled {background: transparent; color: %6;"
          " border: 1px solid %7;}"
          "QPushButton:hover {background: %8; border: 1px solid %9;"
          " border-radius: %3; color:%10}"
          "QPushButton:pressed {background: %11; border: 1px solid %12;"
          " border-radius: %3; color:%13}")
        .arg(strPadding).arg(nFont).arg(nRadius)
        .arg(strBorder, strText, strDisTx, strDisBd,
             strHoverBg, strHoverBd, strHoverTx)
        .arg(strPressBg, strPressBd, strPressTx);
  };

  const QString strSmallBack = Pill("4 12 4 12", 14, 14);
  // Calibration pictograms.
  {
    // Darker green in light mode for contrast on this panel.
    const QColor cTarget(bDark ? "#8bd48b" : "#55a355");
    const QColor cRest(bDark ? "#a8a8a8" : "#6a6a6a");
    const int nPict = 64;
    ui->label_7->setPixmap(MarkerPictogram(nPict, true, cTarget, cRest));
    ui->label_8->setPixmap(MarkerPictogram(nPict, false, cTarget, cRest));
    // The pixmap carries its own device pixel ratio: don't rescale it, and fix the
    // label size so the layout does not squeeze it.
    ui->label_7->setScaledContents(false);
    ui->label_8->setScaledContents(false);
    ui->label_7->setFixedSize(nPict, nPict);
    ui->label_8->setFixedSize(nPict, nPict);
    ui->label_7->setAlignment(Qt::AlignCenter);
    ui->label_8->setAlignment(Qt::AlignCenter);
  }
  if (m_pToolBrush)
  {
    // The armed tool uses Segmentation's accent.
    const QString strToolCss = QString(
        "QToolButton {background: %1; border: 1px solid %2; border-radius: 13px;"
        "  padding: 4px 12px; font-size: 13px; color: %3;}"
        "QToolButton:hover {background: %4;}"
        "QToolButton:checked {background: #a6a233; border: 1px solid #a6a233; color: white;}"
        "QToolButton:checked:hover {background: #928e2c; border: 1px solid #928e2c;}")
        .arg(bDark ? "#3a3a3a" : "#e6e6e6",
             bDark ? "#4f4f4f" : "#cfcfcf",
             bDark ? "#d8d8d8" : "#444444",
             bDark ? "#464646" : "#dcdcdc");
    m_pToolBrush->setStyleSheet(strToolCss);
    m_pToolEraser->setStyleSheet(strToolCss);
    m_pToolFill->setStyleSheet(strToolCss);
  }
  ui->corrContinueCaption->setStyleSheet(
        QString("#corrContinueCaption {background: transparent; font-size: 13px; color: %1;}")
        .arg(bDark ? "#a8a8a8" : "#6a6a6a"));
  if (m_pSegContinueCaption)
  {
    m_pSegContinueCaption->setStyleSheet(
          QString("#segContinueCaption {background: transparent; font-size: 13px; color: %1;}")
          .arg(bDark ? "#a8a8a8" : "#6a6a6a"));
  }
  if (m_pSegPreviewNote)
  {
    // The shortcuts card's surface, marked with Segmentation's accent.
    m_pSegPreviewNote->setStyleSheet(
          QString("#segPreviewNote {background: %1; color: %2; font-size: 14px; border-left: 3px solid #a6a233;"
                  " border-radius: 4px; padding: 8px 10px;}")
          .arg(bDark ? "#3a3a3a" : "#e4e4e4", bDark ? "#e3e3e3" : "#262626"));
  }
  if (m_pCalibCap1)
  {
    const QString strCalibCss =
        QString("#calibCaption {background: transparent; font-size: 13px; color: %1;}"
                "#calibProgress {background: transparent; font-size: 13px; color: %1;}")
        .arg(bDark ? "#a8a8a8" : "#6a6a6a");
    m_pCalibCap1->setStyleSheet(strCalibCss);
    m_pCalibCap2->setStyleSheet(strCalibCss);
    if (m_pCalibProgress)
      m_pCalibProgress->setStyleSheet(strCalibCss);
  }
  if (m_pStartCaption)
  {
    // Transparent: style.qss fills every QWidget, QLabel included.
    m_pStartCaption->setStyleSheet(
          QString("#startCaption {background: transparent; font-size: 13px; color: %1;}")
          .arg(bDark ? "#a8a8a8" : "#6a6a6a"));
  }
  if (m_pReconPage)
    m_pReconPage->SetDarkMode(bDark);
  ui->pushButtonBackToCorrection->setStyleSheet(strSmallBack);
  ui->pushButtonBackToSegmentation->setStyleSheet(strSmallBack);
  ui->pushButtonBackToStartFromCorrection->setStyleSheet(strSmallBack);

  const QString strCalibration = Pill("5 14 5 14", 14, 14);
  // Calibration's Undo/Redo are filled with the page accent; Clear stays a ghost
  // pill. 1px more padding than Clear makes up for "border: none".
  const QString strCalibrationAccent =
      "QPushButton {background: #a6a233; padding: 6 15 6 15; font-size: 14px;"
      " border: none; border-radius: 14; color:white}"
      "QPushButton:disabled {background: #999999;}"
      "QPushButton:hover {background: #928e2c;}"
      "QPushButton:pressed {background: #7e7b26;}";
  ui->pushButtonUndoCalibration->setStyleSheet(strCalibrationAccent);
  // "Callibration": the .ui's object name is misspelled.
  ui->pushButtonRedoCallibration->setStyleSheet(strCalibrationAccent);
  ui->pushButtonClearCalibration->setStyleSheet(strCalibration);

  // CC Undo/Redo are filled with the page accent, 40px tall like the Save row
  // (the extra 1px padding makes up for "border: none").
  const QString strUndoRedoCC =
      "QPushButton {background: #457c92; padding: 9 15 9 15; font-size: 14px;"
      " border: none; border-radius: 18; color:white}"
      "QPushButton:disabled {background: #999999;}"
      "QPushButton:hover {background: #3c6d80;}"
      "QPushButton:pressed {background: #345e6e;}";
  ui->pushButtonUndoCC->setStyleSheet(strUndoRedoCC);
  ui->pushButtonRedoCC->setStyleSheet(strUndoRedoCC);

  const QString strWidePill = Pill("8 18 8 18", 17, 18);
  ui->pushButtonClearCC->setStyleSheet(strWidePill);
  // Correction's Clear points sits beside the 14px "Correct this photo", so
  // it takes that size (both 38px tall) rather than the wide 17px pill.
  ui->pushButtonClear->setStyleSheet(Pill("7 17 7 17", 14, 18));
  ui->pushButtonBackToStart->setStyleSheet(strWidePill);

  ui->pushButtonChangeOutputFolder->setStyleSheet(Pill("3 12 3 12", 12, 12));
}

// Rich text cannot be tinted by QSS, so each hotkey glyph is a data-URL <img>
// pre-coloured for the theme, rendered at 32px for a 16px slot (crisp at DPR 2).
static QString HotkeyGlyphImg(const QString& strName, const QColor& color)
{
  const QPixmap pix = AppTheme::TintedPixmap(strName, color);
  if (pix.isNull())
    return QString();   // resource missing: the row degrades to text-only
  const QImage img = pix.toImage().scaled(32, 32, Qt::IgnoreAspectRatio,
                                          Qt::SmoothTransformation);
  QByteArray ba;
  QBuffer buf(&ba);
  buf.open(QIODevice::WriteOnly);
  img.save(&buf, "PNG");
  return QString("<img src=\"data:image/png;base64,%1\" width=\"16\" height=\"16\">")
      .arg(QString::fromLatin1(ba.toBase64()));
}

// The window's message box: bold headline, body, optional dim caption (a path
// or variable), buttons named by outcome. Button 0 is the default. Escape and
// the close box return the last button, except with nDanger: then button 0 is
// also the escape button, so any dismissal is safe and returns -1.
// The title is repeated as a headline because macOS hides dialog titles; a
// spacer row stops QMessageBox from collapsing to ~240px.
static int ShowBox(QWidget* parent, const QString& strTitle, const QString& strBody,
                   const QString& strCaption, const QStringList& listButtons,
                   int nDanger)
{
  QMessageBox box(parent);
  box.setIcon(QMessageBox::NoIcon);
  box.setWindowTitle(strTitle);
  box.setTextFormat(Qt::RichText);
  box.setText("<b>" + strTitle.toHtmlEscaped() + "</b>");
  QString strHtml = "<p style=\"margin-top:0\">" + strBody.toHtmlEscaped().replace("\n", "<br/>") + "</p>";
  if (!strCaption.isEmpty())
  {
    const QString strDim = AppTheme::EffectiveDarkMode() ? "#a8a8a8" : "#6a6a6a";
    strHtml += QString("<p style=\"color:%1\">%2</p>").arg(strDim, strCaption.toHtmlEscaped());
  }
  box.setInformativeText(strHtml);
  // Headline and body on the type ladder (16 / 14 px); the box's own sheet,
  // never qApp's.
  box.setStyleSheet("QLabel#qt_msgbox_label { font-size: 16px; }"
                    "QLabel#qt_msgbox_informativelabel { font-size: 14px; }"
                    "QPushButton { font-size: 14px; padding: 7px 16px; }");

  QList<QPushButton*> buttons;
  for (int i = 0; i < listButtons.size(); i++)
  {
    QMessageBox::ButtonRole role = QMessageBox::ActionRole;
    if (i == 0)
      role = QMessageBox::AcceptRole;
    else if (i == listButtons.size() - 1)
      role = QMessageBox::RejectRole;
    buttons << box.addButton(listButtons[i], role);
    buttons.last()->setCursor(Qt::PointingHandCursor);
  }
  box.setDefaultButton(buttons.first());
  // nDanger 0 is ignored: button 0 is the safe default and escape button,
  // so it can never be the destructive one.
  const bool bDanger = (nDanger > 0 && nDanger < buttons.size());
  box.setEscapeButton((buttons.size() > 1 && !bDanger) ? buttons.last() : buttons.first());
  if (bDanger)
  {
    // Safe choice: filled teal on the right. Destructive: red text button to its
    // left.
    const bool bDark = AppTheme::EffectiveDarkMode();
    const QString strShape = "border: 1px solid transparent; border-radius: 8px; padding: 7px 16px; font-size: 14px; font-weight: bold;";
    buttons.first()->setStyleSheet(
        "QPushButton { background: #457c92; color: white; " + strShape + " }"
        "QPushButton:hover { background: #3c6d80; } QPushButton:pressed { background: #34606f; }");
    buttons[nDanger]->setStyleSheet(
        QString("QPushButton { background: transparent; color: %1; %2 }"
                "QPushButton:hover { background: %3; }"
                "QPushButton:pressed { background: #a2293c; color: white; }")
            .arg(bDark ? "#e08a8d" : "#a2293c", strShape, bDark ? "#3a2528" : "#f6e4e7"));
    if (QDialogButtonBox* pButtonBox = box.findChild<QDialogButtonBox*>())
      pButtonBox->setStyleSheet("QDialogButtonBox { button-layout: 3; }");
  }
  if (QGridLayout* grid = qobject_cast<QGridLayout*>(box.layout()))
    grid->addItem(new QSpacerItem(400, 0, QSizePolicy::Minimum, QSizePolicy::Minimum),
                  grid->rowCount(), 0, 1, grid->columnCount());
  box.exec();
  return buttons.indexOf(qobject_cast<QPushButton*>(box.clickedButton()));
}

// Single-button form. Text after the first blank line becomes the caption.
static void Warn(QWidget* parent, const QString& strTitle, const QString& strBody)
{
  const int nSplit = strBody.indexOf("\n\n");
  const QString strMain = (nSplit < 0) ? strBody : strBody.left(nSplit);
  const QString strCaption = (nSplit < 0) ? QString() : strBody.mid(nSplit + 2).trimmed();
  ShowBox(parent, strTitle, strMain, strCaption, QStringList() << "Close");
}

// Calibration click pictograms: outer ring plus centre dot, the part to click
// in the target colour. Painted at 2x.
static QPixmap MarkerPictogram(int nLogicalSize, bool bTargetCentre,
                               const QColor& target, const QColor& rest)
{
  const int nDev = nLogicalSize * 2;
  QPixmap out(nDev, nDev);
  out.fill(Qt::transparent);
  QPainter p(&out);
  p.setRenderHint(QPainter::Antialiasing, true);
  const QPointF c(nDev / 2.0, nDev / 2.0);
  const qreal rRing = nDev * 0.40, rDot = nDev * 0.11, wRing = 4.0;
  p.setBrush(Qt::NoBrush);
  p.setPen(QPen(bTargetCentre ? rest : target, wRing));
  p.drawEllipse(c, rRing, rRing);
  p.setPen(Qt::NoPen);
  p.setBrush(bTargetCentre ? target : rest);
  p.drawEllipse(c, rDot, rDot);
  p.end();
  out.setDevicePixelRatio(2.0);
  return out;
}

namespace
{
struct HotkeyRow { const char* action; const char* keys; };
}

// Shortcuts card table (glyph | action | keys), shared by the editing pages. An
// empty action is a group break; rows 0 and 1 are Pan and Zoom and carry
// glyphs. white-space:pre keeps each chord on one line; the font size is set
// because rich text ignores the widget font.
static QString HotkeysHtml(const HotkeyRow* rows, size_t nRows, const QColor& colorGlyph)
{
  const QString strPan  = HotkeyGlyphImg("pan",  colorGlyph);
  const QString strZoom = HotkeyGlyphImg("zoom", colorGlyph);
  QString strHtml = "<html><head/><body style=\" font-size:12px;\">"
                    "<p style=\" letter-spacing:1px; margin-bottom:4px;\"><b>SHORTCUTS</b></p>"
                    "<table cellspacing=\"0\" cellpadding=\"1\">";
  for (size_t i = 0; i < nRows; ++i)
  {
    if (!*rows[i].action)
    {
      strHtml += "<tr><td colspan=\"3\" style=\" font-size:4px;\">&nbsp;</td></tr>";
      continue;
    }
    const QString strGlyph = (i == 0) ? strPan : (i == 1) ? strZoom : QString();
    strHtml += QString("<tr>"
                       "<td style=\" padding-right:6px;\">%1</td>"
                       "<td style=\" padding-right:16px; white-space:pre;\">%2</td>"
                       "<td style=\" white-space:pre;\">%3</td>"
                       "</tr>")
                   .arg(strGlyph, QString::fromLatin1(rows[i].action),
                        QString::fromLatin1(rows[i].keys));
  }
  strHtml += "</table></body></html>";
  return strHtml;
}

static QString SegHotkeysHtml(const QColor& colorGlyph)
{
  // Every Segmentation binding, taken from MaskActionFor(), the tool actions and
  // the window actions. No wheel row: wheelEvent() does nothing.
  const HotkeyRow rows[] = {
    {"Pan",                   "Left drag"},
    {"Zoom",                  "Right drag"},
    {"",                      ""},                          // group break
    {"Brush / Eraser / Fill", "Q / W / E"},
    {"Use selected tool",     "Ctrl + Left drag"},
    {"",                      ""},
    {"Flood fill",            "Ctrl + Right click"},
    {"Flood erase",           "Ctrl + Shift + Right click"},
    {"",                      ""},
    {"Undo / Redo",           "Ctrl + Z / Ctrl + Shift + Z"},
    {"Show mask",             "Ctrl + T"},
    {"",                      ""},
    {"Previous / Next photo", "Ctrl + A / Ctrl + D"},
    {"Continue",              "Ctrl + S"},
  };
  return HotkeysHtml(rows, sizeof(rows) / sizeof(rows[0]), colorGlyph);
}

static QString CorrectionHotkeysHtml(const QColor& colorGlyph)
{
  // No keyboard shortcuts act on points.
  const HotkeyRow rows[] = {
    {"Pan",                   "Left drag"},
    {"Zoom",                  "Right drag"},
    {"",                      ""},
    {"Place point",           "Ctrl + Left click"},
    {"Move last point",       "Ctrl + click once all placed"},
    {"Correct and next",      "Ctrl + W"},
    {"",                      ""},
    {"Previous / Next photo", "Ctrl + A / Ctrl + D"},
    {"Continue",              "Ctrl + S"},
  };
  return HotkeysHtml(rows, sizeof(rows) / sizeof(rows[0]), colorGlyph);
}

// Each shortcut presses the page's own button, and only when it is enabled, so
// every greyed-out rule applies unchanged. (Ctrl is Command on macOS.)
void MainWindow::SetUpNavigationShortcuts()
{
  // 0 previous, 1 next, 2 continue, 3 do this photo and go to the next one
  struct Key { const char* seq; int nWhich; };
  const Key keys[] = { { "Ctrl+A", 0 }, { "Ctrl+D", 1 }, { "Ctrl+S", 2 }, { "Ctrl+W", 3 } };
  for (const Key& k : keys)
  {
    QAction* a = new QAction(this);
    a->setShortcut(QKeySequence(QString::fromLatin1(k.seq)));
    const int nWhich = k.nWhich;
    connect(a, &QAction::triggered, this, [this, nWhich]()
    {
      QWidget* pPage = ui->stackedWidget->currentWidget();
      auto press = [](QAbstractButton* b)
      {
        // isEnabled() also covers a disabled page or panel around the button.
        if (!b || !b->isEnabled() || !b->isVisible())
          return false;
        b->click();
        return true;
      };
      if (nWhich == 3)
      {
        if (pPage == ui->pageCorrection)
        {
          // Next only once this press started a correction; a refused one (e.g. too
          // few points) stays put.
          const int n = m_nIndex;
          const bool bWasRunning = m_proc->state() != QProcess::NotRunning;
          if (press(ui->pushButtonProcess) && m_nIndex == n && !bWasRunning
              && m_proc->state() != QProcess::NotRunning && m_setProcessStarted.contains(n))
            press(ui->pushButtonNext);
        }
        else if (pPage == ui->pageCC)
        {
          if (press(ui->pushButtonCreateMask))
            press(ui->pushButtonNextCC);
        }
        return;
      }
      // Typing in a text field keeps Ctrl + A / Ctrl + D (select all, delete).
      QWidget* pFocus = QApplication::focusWidget();
      if (nWhich != 2 && (qobject_cast<QLineEdit*>(pFocus) || qobject_cast<QAbstractSpinBox*>(pFocus)
                          || qobject_cast<QPlainTextEdit*>(pFocus) || qobject_cast<QTextEdit*>(pFocus)))
        return;
      QAbstractButton* pButtons[3] = { nullptr, nullptr, nullptr };
      if (pPage == ui->pageCorrection)
      {
        pButtons[0] = ui->pushButtonPrevious; pButtons[1] = ui->pushButtonNext; pButtons[2] = ui->pushButtonSegmentation;
      }
      else if (pPage == ui->pageSegEdit)
      {
        pButtons[0] = ui->pushButtonPreviousSeg; pButtons[1] = ui->pushButtonNextSeg; pButtons[2] = ui->pushButtonCC;
      }
      else if (pPage == ui->pageCC)
      {
        pButtons[0] = ui->pushButtonPreviousCC; pButtons[1] = ui->pushButtonNextCC; pButtons[2] = ui->pushButtonQuit;
      }
      press(pButtons[nWhich]);
    });
    addAction(a);
  }
}

static QString CCHotkeysHtml(const QColor& colorGlyph)
{
  const HotkeyRow rows[] = {
    {"Pan",                   "Left drag"},
    {"Zoom",                  "Right drag"},
    {"",                      ""},
    {"Draw a box",            "Ctrl + Left drag"},
    {"Undo / Redo",           "Ctrl + Z / Ctrl + Shift + Z"},
    {"Save components and next", "Ctrl + W"},
    {"",                      ""},
    {"Previous / Next photo", "Ctrl + A / Ctrl + D"},
    {"Continue",              "Ctrl + S"},
  };
  return HotkeysHtml(rows, sizeof(rows) / sizeof(rows[0]), colorGlyph);
}

void MainWindow::SetDarkMode(bool bDark)
{
  if (m_pHomeButton)
  {
    // Flat on the bar (#fafafa / #383838), with the step pills' hover tone and
    // the current step's blue.
    m_pHomeButton->setIcon(AppTheme::TintedIcon("home", bDark ? QColor("#6fb3c9") : QColor("#35637a")));
    m_pHomeButton->setStyleSheet(QString("#stepperHome {background: transparent; border: none; border-radius: 12px;}"
                                         "#stepperHome:hover {background: %1;}"
                                         "#stepperHome:pressed {background: %2;}")
                                 .arg(bDark ? "#474747" : "#ececec", bDark ? "#525252" : "#e0e0e0"));
  }
  m_bDarkMode = bDark;
  qApp->setPalette(AppTheme::BuildPalette(bDark));

  // Input fields get per-theme fills and a fake inset bevel (darker top/left,
  // lighter bottom/right). Same selector as style.qss's input rule and appended
  // after it, so these win per property. Keep the tones in sync with style.qss.
  QString strInputCss;
  if (bDark)
  {
    strInputCss = "QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {"
                  "  background: #333333; color: #f0f0f0;"
                  "  border-top: 1px solid #262626; border-left: 1px solid #262626;"
                  "  border-bottom: 1px solid #444444; border-right: 1px solid #444444;"
                  "  border-radius: 4px; padding: 3px 6px;}"
                  "QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {"
                  "  border: 1px solid #6fb3c9;}"
                  // Top-bar fields take the bar's tone; id + type outranks the rule above.
                  "#frameTopBar QLineEdit {"
                  "  background: #383838; color: #f0f0f0;"
                  "  border-top: 1px solid #2a2a2a; border-left: 1px solid #2a2a2a;"
                  "  border-bottom: 1px solid #4c4c4c; border-right: 1px solid #4c4c4c;"
                  "  border-radius: 4px; padding: 3px 6px;}"
                  "#frameTopBar QLineEdit:focus {border: 1px solid #6fb3c9;}";
  }
  else
  {
    strInputCss = "QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {"
                  "  background: #f7f7f7; color: #262626;"
                  "  border-top: 1px solid #d9d9d9; border-left: 1px solid #d9d9d9;"
                  "  border-bottom: 1px solid #fbfbfb; border-right: 1px solid #fbfbfb;"
                  "  border-radius: 4px; padding: 3px 6px;}"
                  "QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {"
                  "  border: 1px solid #457c92;}"
                  "#frameTopBar QLineEdit {"
                  "  background: #fafafa; color: #262626;"
                  "  border-top: 1px solid #dcdcdc; border-left: 1px solid #dcdcdc;"
                  "  border-bottom: 1px solid #ffffff; border-right: 1px solid #ffffff;"
                  "  border-radius: 4px; padding: 3px 6px;}"
                  "#frameTopBar QLineEdit:focus {border: 1px solid #457c92;}";
  }

  // Bar captions must not take the app-wide QWidget background.
  strInputCss += "#frameTopBar QLabel {background: transparent;}";
  strInputCss += "#frameBottomBar QLabel {background: transparent;}";

  QColor colorPreview;
  QString strPanelCss;
  if (bDark)
  {
    colorPreview = QColor(0x10, 0x10, 0x10);
    strPanelCss = "#widgetControl {background: #2a2a2a;} "
                          "#widgetMain {background: #101010;} "
                          "#frameTopBar, #frameBottomBar {background: #383838;} "
                          // In-window only on Linux; darker than the column and bars, as chrome.
                          "QMenuBar {background: #1e1e1e; color: #f0f0f0;} "
                          // Required: once QMenuBar takes a stylesheet background,
                          // Fusion otherwise draws the items from the palette.
                          "QMenuBar::item {background: transparent;} "
                          "QMenuBar::item:selected {background: #457c92; color: #ffffff;} "
                          "#labelTitle, #labelTitle_2, #labelTitle_3, #label_9, #label_11 "
                          // Page titles: bottom rule only; radius 0 so the edge does not curl.
                          "{border: none; border-bottom: 1px solid #555555; border-radius: 0; padding: 6px 2px;} "
                          "#frameStartupDivider, #frameSegmentationProceedDivider, "
                          "#frameCCProceedDivider, #frameFinishCaseDivider {background: #3a3a3a;} "
                          // On #frameBottomBar, so it needs its own tone to be visible.
                          "#frameCarouselDivider {background: #4f4f4f;} "
                          "#labelSegHotkeys, #labelCorrectionHotkeys, #labelCCHotkeys {background: #3a3a3a;"
                          " border-radius: 6px; padding: 8px 10px;} "
                          // pageSegEdit is disabled while nnUNet runs; keep its labels at full text
                          // colour.
                          "#pageSegEdit QLabel:disabled {color: white;} "
                          // Transparent, listed by id: "#widgetMain QWidget" would outrank the top-bar input rules.
                          "#widgetThumbnail, #widgetIndex, #widgetIndexSeg, #widgetIndexCC, "
                          "#widgetPointMode {background: transparent;} "
                          + strInputCss;
    ui->labelSegHotkeys->setText(SegHotkeysHtml(QColor(0xe3, 0xe3, 0xe3)));
    ui->labelCorrectionHotkeys->setText(CorrectionHotkeysHtml(QColor(0xe3, 0xe3, 0xe3)));
    ui->labelCCHotkeys->setText(CCHotkeysHtml(QColor(0xe3, 0xe3, 0xe3)));
    // The inline progress captions carry a light-theme grey in the .ui.
    ui->labelSegProgress->setStyleSheet("color: #a8a8a8; font-size: 13px;");
    ui->labelCCProgress->setStyleSheet("color: #a8a8a8; font-size: 13px;");
  }
  else
  {
    colorPreview = QColor(0xc6, 0xc6, 0xc6);
    strPanelCss = "#widgetControl {background: #f0f0f0;} "
                          "#widgetMain {background: #c6c6c6;} "
                          "#frameTopBar, #frameBottomBar {background: #fafafa;} "
                          "QMenuBar {background: #dcdcdc; color: #262626;} "
                          "QMenuBar::item {background: transparent;} "
                          "QMenuBar::item:selected {background: #457c92; color: #ffffff;} "
                          "#labelTitle, #labelTitle_2, #labelTitle_3, #label_9, #label_11 "
                          "{border: none; border-bottom: 1px solid #d0d0d0; border-radius: 0; padding: 6px 2px;} "
                          "#frameStartupDivider, #frameSegmentationProceedDivider, "
                          "#frameCCProceedDivider, #frameFinishCaseDivider {background: #dcdcdc;} "
                          "#frameCarouselDivider {background: #d8d8d8;} "
                          "#labelSegHotkeys, #labelCorrectionHotkeys, #labelCCHotkeys {background: #e4e4e4;"
                          " border-radius: 6px; padding: 8px 10px;} "
                          "#pageSegEdit QLabel:disabled {color: #262626;} "
                          "#widgetThumbnail, #widgetIndex, #widgetIndexSeg, #widgetIndexCC, "
                          "#widgetPointMode {background: transparent;} "
                          + strInputCss;
    ui->labelSegHotkeys->setText(SegHotkeysHtml(QColor(0x44, 0x44, 0x44)));
    ui->labelCorrectionHotkeys->setText(CorrectionHotkeysHtml(QColor(0x44, 0x44, 0x44)));
    ui->labelCCHotkeys->setText(CCHotkeysHtml(QColor(0x44, 0x44, 0x44)));
    ui->labelSegProgress->setStyleSheet("color: #666666; font-size: 13px;");
    ui->labelCCProgress->setStyleSheet("color: #666666; font-size: 13px;");
  }
  // WidgetImageView paints its own background: keep it the #widgetMain colour.
  ui->widgetImageView->SetBackgroundColor(colorPreview);
  m_pImageGridView->SetDarkMode(bDark);
  RethemeGhostPillButtons(bDark);
  EmboldenNavButtons({ui->pushButtonPrevious, ui->pushButtonNext,
                      ui->pushButtonPreviousSeg, ui->pushButtonNextSeg,
                      ui->pushButtonPreviousCC, ui->pushButtonNextCC,
                      ui->pushButtonMaskEditUndo, ui->pushButtonMaskEditRedo,
                      ui->pushButtonUndoCC, ui->pushButtonRedoCC,
                      ui->pushButtonUndoCalibration, ui->pushButtonRedoCallibration});
  qApp->setStyleSheet(AppTheme::BaseStyleSheet(bDark) + "\n" + strPanelCss);
  // Size the top-bar fields for four digits in the applied font (+ padding and
  // border). Measured after setStyleSheet() re-polishes.
  {
    const int nField = ui->lineEditRectWidth->fontMetrics().horizontalAdvance("0000") + 14;
    ui->lineEditRectWidth->setFixedWidth(nField);
    ui->lineEditRectHeight->setFixedWidth(nField);
  }
}

void MainWindow::OnToggleDarkMode(bool bDark)
{
  // An explicit choice: persist it for the next start and the Hub.
  AppTheme::SetMode(bDark ? "dark" : "light");
  SetDarkMode(bDark);
}

void MainWindow::moveEvent(QMoveEvent* e)
{
  RepositionProgressWindow();
  QMainWindow::moveEvent(e);
}

void MainWindow::RepositionProgressWindow()
{
  QRect rc = m_wndProgress->rect();
  rc.moveCenter(rect().center());
  rc.moveTopLeft(rc.topLeft()+this->pos());
  m_wndProgress->setGeometry(rc);
}

void MainWindow::LoadProject(const QString& projectFolder)
{
  if (projectFolder.isEmpty())
    return;

  ProjectInfo info;
  if (!ProjectManager::LoadProject(projectFolder, info))
    return; // Folder no longer exists; land on the normal empty StartUp page.

  // No build saves "reconstruction"; treat it as connected_components. Whether
  // Reconstruction is offered is decided from the files below.
  if (info.currentStage == "reconstruction")
    info.currentStage = "connected_components";

  m_strCurrentProjectFolder = projectFolder;
  m_currentProjectInfo = info;
  setWindowTitle(info.name.isEmpty() ? "Dissection Photo" : QString("Dissection Photo — %1").arg(info.name));

  SetInputFolder(info.inputFolder, info.rootOutputFolder);
  if (!info.calibrationFile.isEmpty())
    SetCalibrationFile(info.calibrationFile);

  // Ask which step to resume at instead of jumping to the saved stage. Not while
  // this project's reconstruction runs (left running by Switch Project): back to it.
  QString strStage = info.currentStage;
  ReconRunner* pRun = ReconRunner::ForRoot(info.rootOutputFolder);
  if (pRun && pRun->IsActive())
    strStage = "reconstruction";
  if (!strStage.isEmpty() && strStage != "correction" && strStage != "reconstruction")
  {
    auto StageDisplayName = [](const QString& s) -> QString {
      if (s == "connected_components")
        return "Connected Components";
      if (s == "segmentation")
        return "Segmentation";
      return "Correction";
    };

    // A plain QDialog, not QMessageBox: QMessageBox caps its width, which clipped
    // these buttons, and orders buttons by role, not insertion.
    QDialog dlg(this);
    dlg.setWindowTitle("Resume Project");

    // Whether the case is finished, tested against the folder:
    // AllCCMasksSaved() needs state the replay has not filled yet. Requires a
    // corrected photo, a mask and a non-empty CC .npy for every photo; a mask
    // mismatch would make the Segmentation replay rerun nnUNet.
    const QString strRoot = info.rootOutputFolder;
    bool bFinished = (strStage == "connected_components");
    if (bFinished)
    {
      const QFileInfoList listDeformed = QDir(strRoot + "/deformed").entryInfoList(ImageNameFilters(), QDir::Files, QDir::Name);
      bFinished = !listDeformed.isEmpty()
          && listDeformed.size() == QDir(info.inputFolder).entryInfoList(ImageNameFilters(), QDir::Files).size()
          && listDeformed.size() == QDir(strRoot + "/mask").entryInfoList(QDir::Files).size();
      foreach (const QFileInfo& fi, listDeformed)
      {
        if (!bFinished)
          break;
        const QString strSaved = QFileInfo(strRoot + "/connected_components", fi.completeBaseName() + "_mask.npy").absoluteFilePath();
        // Components saved before the mask was last edited are not finished.
        bFinished = IsCompleteDataOutput(strSaved)
            && !CCSavedIsStale(strSaved, strRoot + "/mask/" + fi.completeBaseName() + ".png");
      }
    }

    const QString strName = info.name.isEmpty() ? QDir(projectFolder).dirName() : info.name;
    QLabel* pLabel = new QLabel(bFinished
          ? QString("\"%1\" is finished: every photo has its connected components saved.\nWhere would you like to continue?").arg(strName)
          : QString("\"%1\" was last at the %2 step.\nWhere would you like to continue?").arg(strName, StageDisplayName(strStage)),
          &dlg);
    pLabel->setStyleSheet("font-size: 14px;");
    pLabel->setWordWrap(true);

    // The dark tones match RethemeGhostPillButtons().
    const QString strBackCss = m_bDarkMode
        ? "QPushButton {background: transparent; padding: 6 16 6 16; font-size: 14px; font-weight: bold; border: 1px solid #777777; border-radius: 14; color:#dddddd}"
          "QPushButton:hover {background: #4a4a4a; border: 1px solid #999999; color:#ffffff}"
          "QPushButton:pressed {background: #2e2e2e; border: 1px solid #aaaaaa; color:#ffffff}"
        : "QPushButton {background: transparent; padding: 6 16 6 16; font-size: 14px; font-weight: bold; border: 1px solid #999999; border-radius: 14; color:#555555}"
          "QPushButton:hover {background: #f0f0f0; border: 1px solid #888888; color:#333333}"
          "QPushButton:pressed {background: #e2e2e2; border: 1px solid #777777; color:#333333}";
    auto StepCss = [](const QString& color, const QString& hover, const QString& pressed) {
      return QString("QPushButton {background: %1; padding: 6 20 6 20; font-size: 14px; font-weight: bold; border: none; border-radius: 14; color:white}"
                     "QPushButton:hover {background: %2;}"
                     "QPushButton:pressed {background: %3;}").arg(color, hover, pressed);
    };

    // Only a button click sets it; Escape or closing leaves it empty.
    QString strChosen;
    auto AddButton = [&](const QString& strText, const QString& strCss,
                         const QString& strStageValue) -> QPushButton* {
      QPushButton* pBtn = new QPushButton(strText, &dlg);
      pBtn->setStyleSheet(strCss);
      pBtn->setCursor(Qt::PointingHandCursor);
      // Never squeeze a button below its text width.
      pBtn->setMinimumWidth(pBtn->sizeHint().width());
      connect(pBtn, &QPushButton::clicked, &dlg, [&strChosen, strStageValue, &dlg]() {
        strChosen = strStageValue;
        dlg.accept();
      });
      return pBtn;
    };

    QHBoxLayout* pButtonRow = new QHBoxLayout;
    pButtonRow->setSpacing(10);
    // Back first, then the steps in pipeline order.
    QPushButton* pBackBtn = AddButton("Back to Projects", strBackCss, "__back__");
    pButtonRow->addWidget(pBackBtn);
    pButtonRow->addStretch(1);
    // The furthest step offered is the default, never "Back to Projects".
    QPushButton* pLastStepBtn = AddButton("Correction", StepCss("#a2293c", "#8e2434", "#7b1f2d"), "correction");
    pButtonRow->addWidget(pLastStepBtn);
    if (strStage == "segmentation" || strStage == "connected_components")
    {
      pLastStepBtn = AddButton("Segmentation", StepCss("#a6a233", "#928e2c", "#7e7b26"), "segmentation");
      pButtonRow->addWidget(pLastStepBtn);
    }
    if (strStage == "connected_components")
    {
      pLastStepBtn = AddButton("Connected Components", StepCss("#457c92", "#3c6d80", "#34606f"), "connected_components");
      pButtonRow->addWidget(pLastStepBtn);

      if (bFinished)
      {
        pLastStepBtn = AddButton("Reconstruction", StepCss("#89448e", "#783b7c", "#68336b"), "reconstruction");
        pButtonRow->addWidget(pLastStepBtn);
      }
    }

    QVBoxLayout* pDlgLayout = new QVBoxLayout(&dlg);
    pDlgLayout->setContentsMargins(20, 18, 20, 16);
    pDlgLayout->setSpacing(18);
    pDlgLayout->addWidget(pLabel);
    pDlgLayout->addLayout(pButtonRow);

    dlg.setSizeGripEnabled(false);
    pLastStepBtn->setDefault(true);

    dlg.exec();

    if (strChosen.isEmpty())
      return;
    if (strChosen == "__back__")
    {
      emit RequestSwitchProject();
      return;
    }
    strStage = strChosen;
  }

  // Replay forward navigation up to the chosen stage; each handler re-validates
  // against disk as on a fresh visit.
  m_bResumingProject = true;
  const bool bToRecon = (strStage == "reconstruction");
  if (strStage == "correction" || strStage == "segmentation" || strStage == "connected_components" || bToRecon)
    OnButtonContinue();
  if (strStage == "segmentation" || strStage == "connected_components" || bToRecon)
    OnButtonProceedToSeg();
  if (strStage == "connected_components" || bToRecon)
    OnButtonProceedToCC();
  // Only if the replay reached Connected components and the case is still
  // finished (the same test as the Reconstruction gate).
  if (bToRecon && ui->stackedWidget->currentWidget() == ui->pageCC && CCReadyForRecon())
    OnButtonProceedToRecon();
  m_bResumingProject = false;
  // Deferred from the replay's OnButtonContinue(): a later stage starts the
  // daemon itself, with only the photos it needs.
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
    SegmentCorrectedPhotos();
}

void MainWindow::SaveProjectState(const QString& stage)
{
  if (m_strCurrentProjectFolder.isEmpty())
    return;

  m_currentProjectInfo.rootOutputFolder = m_strCurrentProjectFolder;
  m_currentProjectInfo.inputFolder = m_strInputFolder;
  m_currentProjectInfo.calibrationFile = m_strCalibrationFile;
  m_currentProjectInfo.calibratedMode = m_bCalibratedMode;
  // currentStage never moves back: LoadProject()'s replay looks like forward
  // navigation, and resuming at Correction would otherwise erase later stages.
  auto StageRank = [](const QString& s) -> int {
    if (s == "connected_components") return 3;
    if (s == "segmentation")         return 2;
    if (s == "correction")           return 1;
    return 0;
  };
  if (StageRank(stage) >= StageRank(m_currentProjectInfo.currentStage))
    m_currentProjectInfo.currentStage = stage;
  if (m_currentProjectInfo.inputImageCountAtCreation <= 0)
    m_currentProjectInfo.inputImageCountAtCreation = m_listInputFiles.size();
  ProjectManager::SaveProject(m_currentProjectInfo);
}

void MainWindow::OnActionSwitchProject()
{
  // Keep edits not yet saved: a mask on Segmentation, boxes on Connected components.
  if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
    ui->widgetImageView->SaveMaskIfEdited();
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
    SaveCCPhoto();
  StopSegDaemon();
  emit RequestSwitchProject();
  // A running reconstruction keeps going; closeEvent() asks only on quit.
  m_bSwitchingProject = true;
  close();
}

void MainWindow::closeEvent(QCloseEvent* e)
{
  if (!m_bSwitchingProject && !ReconRunner::ConfirmQuit(this))
  {
    e->ignore();
    return;
  }
  QMainWindow::closeEvent(e);
}

void MainWindow::SetInputFolder(const QString &path, const QString& root_output)
{
  QString strActualPath = path;
  QString strActualRootOutput = root_output;

  // A newly picked folder with a project open: the project folder is the output
  // root, and the photos are copied into <project>/photos on Get started
  // (CopyPendingPhotosIntoProject()), so they can be unticked first. Resume
  // passes root_output explicitly; picking photos/ itself copies nothing.
  if (!path.isEmpty() && root_output.isEmpty() && !m_strCurrentProjectFolder.isEmpty())
  {
    QString strPhotosDir = QDir(m_strCurrentProjectFolder).filePath("photos");
    if (QDir(path).absolutePath() != QDir(strPhotosDir).absolutePath())
    {
      m_strPendingSourceFolder = path;
    }
    else
    {
      m_strPendingSourceFolder.clear();
      strActualPath = strPhotosDir;
    }
    strActualRootOutput = m_strCurrentProjectFolder;
  }
  else if (path != m_strPendingSourceFolder)
  {
    // Resume, Clear, or a different folder. Changing only the output root
    // (OnButtonChangeOutputFolder() passes the same path) keeps the pick.
    m_strPendingSourceFolder.clear();
  }

  m_strInputFolder = strActualPath;
  if (strActualPath.isEmpty())
  {
    ui->labelImageFolder->clear();
    ui->labelOutputFolder->clear();
    ui->pushButtonChangeOutputFolder->setVisible(false);
  }
  else
  {
    if (strActualRootOutput.isEmpty())
      m_strRootOutputFolder = strActualPath;
    else
      m_strRootOutputFolder = strActualRootOutput;

    // Font size repeated in the markup: rich text ignores the widget stylesheet.
    ui->labelImageFolder->setText("<html><head/><body style=\" font-size:13px;\"><b>" +
                                  strActualPath + "</b><br /><a href='clear'>Clear</a></body></html>");
    ui->labelOutputFolder->setText("Output: " + m_strRootOutputFolder);
    ui->pushButtonChangeOutputFolder->setVisible(true);
  }
  UpdateStartupPreview();
}

QStringList MainWindow::ImageNameFilters()
{
  static const QStringList filters = QStringList()
      << "*.jpg" << "*.jpeg" << "*.png" << "*.tif" << "*.tiff" << "*.bmp"
      << "*.JPG" << "*.JPEG" << "*.PNG" << "*.TIF" << "*.TIFF" << "*.BMP";
  return filters;
}

bool MainWindow::CopyPendingPhotosIntoProject()
{
  const QString strSource = m_strPendingSourceFolder;
  if (!QDir(strSource).exists())
  {
    Warn(this, "Photo folder is missing",
         "It is no longer where it was. Choose the folder that holds the photos again.");
    return false;
  }

  // Use the grid's picks only if it shows this folder's current photos; otherwise
  // take every photo.
  QStringList listAll;
  foreach (const QFileInfo& fi, QDir(strSource).entryInfoList(ImageNameFilters(), QDir::Files, QDir::Name))
    listAll << fi.absoluteFilePath();
  QStringList listCopy = listAll;
  if (m_pImageGridView->IsSelectable() && m_pImageGridView->Files() == listAll)
    listCopy = m_pImageGridView->SelectedFiles();

  if (listAll.isEmpty())
  {
    Warn(this, "No photos in that folder",
         "Choose a folder that holds photos.\n\nLooked in " + strSource + ".");
    return false;
  }
  if (listCopy.isEmpty())
  {
    Warn(this, "No photographs selected",
         "Tick at least one photograph in the grid, then continue.");
    return false;
  }

  const QString strPhotosDir = QDir(m_strCurrentProjectFolder).filePath("photos");
  QDir().mkpath(strPhotosDir);
  QStringList listFailed;
  foreach (const QString& strFile, listCopy)
  {
    const QFileInfo fi(strFile);
    const QString strDest = QDir(strPhotosDir).filePath(fi.fileName());
    if (!QFile::exists(strDest) && !QFile::copy(fi.absoluteFilePath(), strDest))
      listFailed << fi.fileName();
  }
  if (!listFailed.isEmpty())
  {
    Warn(this, "Some photos could not be copied",
         QString("%1 of %2 photographs could not be copied into %3:\n\n%4")
         .arg(listFailed.size()).arg(listCopy.size()).arg(strPhotosDir, listFailed.join("\n")));
    return false;
  }

  m_strPendingSourceFolder.clear();
  SetInputFolder(strPhotosDir, m_strRootOutputFolder);
  return true;
}

void MainWindow::UpdateStartupPreview()
{
  bool bShowGrid = false;
  QStringList listFiles;
  if (ui->stackedWidget->currentWidget() == ui->pageStartUp &&
      !m_strInputFolder.isEmpty() && QDir(m_strInputFolder).exists())
  {
    QFileInfoList fis = QDir(m_strInputFolder).entryInfoList(ImageNameFilters(), QDir::Files, QDir::Name);
    foreach (const QFileInfo& fi, fis)
      listFiles << fi.absoluteFilePath();
    bShowGrid = !listFiles.isEmpty();
  }
  // Photos can be ticked or left out only while they are still waiting to be
  // copied; the project's own photos folder is shown as it is.
  const bool bPicking = !m_strPendingSourceFolder.isEmpty() &&
      QDir(m_strPendingSourceFolder).absolutePath() == QDir(m_strInputFolder).absolutePath();
  m_pImageGridView->SetSelectable(bPicking);
  if (!bShowGrid || listFiles != m_listStartupPreviewFiles || !bPicking)
    UpdateStartCaption(listFiles.size());
  else
    UpdateStartCaption(m_pImageGridView->SelectedCount());

  // Reconstruction shows its command and log where other steps show a photo.
  // Decided from the page, so every way off the step restores the image view.
  const bool bRecon = (ui->stackedWidget->currentWidget() == ui->pagePhotoRecon);
  if (m_pReconPage)
    m_pReconPage->CommandPane()->setVisible(bRecon);

  if (bShowGrid)
  {
    // Only reload when the file list changed, so the lightbox isn't reset.
    if (listFiles != m_listStartupPreviewFiles)
    {
      m_listStartupPreviewFiles = listFiles;
      m_pImageGridView->SetImages(listFiles);
    }
    ui->widgetImageView->hide();
    m_pImageGridView->show();
  }
  else
  {
    m_listStartupPreviewFiles.clear();
    m_pImageGridView->Clear();
    m_pImageGridView->hide();
    ui->widgetImageView->setVisible(!bRecon);
  }
}

void MainWindow::SetCalibrationFile(const QString &path)
{
  m_strCalibrationFile = path;
  if (path.isEmpty())
  {
    ui->labelCalibrationFile->clear();
  }
  else
  {
    // Bold path like the input folder; 13px in the markup (rich text ignores QSS
    // sizes); Clear link in the box's olive.
    ui->labelCalibrationFile->setText("<html><head/><body style=\" font-size:13px;\"><b>" + path.toHtmlEscaped() +
                                      "</b><br /><a style='color:#a6a233' href='clear'>Clear</a></body></html>");
    // Open the calibration box whenever a file is set (e.g. dropped), so the toggle
    // matches what is shown.
    if (ui->toolButtonToggleCalibration->isChecked())
      OnToggleCalibrationDisclosure(true); // re-target a box still opening
    else
      ui->toolButtonToggleCalibration->setChecked(true);
  }
  if (path.isEmpty() && ui->toolButtonToggleCalibration->isChecked())
    OnToggleCalibrationDisclosure(true);
}

// Eases the calibration box open or closed, starting from its current height so
// re-toggling mid-way reverses smoothly. Before the window is shown it just
// shows or hides. Reduced motion is handled in AnimateHeight().
void MainWindow::OnToggleCalibrationDisclosure(bool bOpen)
{
  QWidget* zone = ui->widgetDropCalibration;
  const bool bAnimating = zone->isVisible() && zone->maximumHeight() < QWIDGETSIZE_MAX;
  const int nFrom = !zone->isVisible() ? 0 : (bAnimating ? zone->maximumHeight() : zone->height());

  // Open height measured at the column's real width (wrapped labels make
  // sizeHint() wrong). The minimum is that height at the narrowest column width,
  // never below the .ui floor, so a short window cannot squeeze the wrapped text.
  int nTo = 0;
  int nFloor = m_nCalibrationBoxFloor;
  QLayout* lay = zone->layout();
  const int nWidth = ui->widgetDropImages->width() > 0 ? ui->widgetDropImages->width() : zone->width();
  // Only once shown: before that the column has no width. showEvent() re-runs
  // this.
  if (isVisible() && lay && lay->hasHeightForWidth() && nWidth > 0)
  {
    const int nSlack = qMax(0, ui->widgetControl->width() - ui->widgetControl->minimumWidth());
    nFloor = qMax(nFloor, lay->totalHeightForWidth(qMax(1, nWidth - nSlack)));
    nTo = lay->totalHeightForWidth(nWidth);
  }
  else
    nTo = zone->sizeHint().height();
  nTo = bOpen ? qMax(nTo, nFloor) : 0;

  // Already fully open (e.g. a file was chosen with the box open): the
  // layout grows or shrinks it on its own.
  if (bOpen && zone->isVisible() && !bAnimating)
  {
    zone->setMinimumHeight(nFloor);
    return;
  }
  if (!bOpen && !zone->isVisible())
    return;

  if (!ui->pageStartUp->isVisible())
  {
    // Off screen: no animation, but still stop one that may be running.
    AppMotion::AnimateHeight(zone, nTo, nTo, AppMotion::Fast, true);
    zone->setVisible(bOpen);
    if (bOpen)
      zone->setMinimumHeight(nFloor);
    return;
  }

  AppMotion::AnimateHeight(zone, nFrom, nTo, AppMotion::Slow, true);
  if (!bOpen)
    return;
  // AnimateHeight() leaves the minimum height at 0 when it finishes; put the
  // floor back then (or now, if it finished at once under reduced motion).
  QAbstractAnimation* pAnim = zone->findChild<QAbstractAnimation*>(QString("AppMotionHeight"), Qt::FindDirectChildrenOnly);
  if (pAnim && pAnim->state() != QAbstractAnimation::Stopped) // Paused: driven by the wall clock
  {
    QPointer<QWidget> guard(zone);
    connect(pAnim, &QAbstractAnimation::finished, this, [guard, nFloor]()
    {
      if (guard && guard->isVisible())
        guard->setMinimumHeight(nFloor);
    });
  }
  else
    zone->setMinimumHeight(nFloor);
}

void MainWindow::SetCalibrationPhoto(const QString &path)
{
  m_strCalibrationPhoto = path;
  if (path.isEmpty())
  {
    ui->labelCalibrationPhoto->clear();
  }
  else
  {
    ui->labelCalibrationPhoto->setText(path);
    OnButtonCalibration();
  }
}

void MainWindow::OnButtonInputFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", m_strCurrentProjectFolder);
  if (!path.isEmpty())
  {
    SetInputFolder(path);
  }
}

void MainWindow::OnButtonOutputFolderCC()
{
  QString path = m_strRootOutputFolder + "/connected_components";
  QDir().mkpath(path);
  QStringList fns;
  QList<int> listSaved, listStale;
  foreach (QFileInfo fi, m_listInputFiles)
  {
    // Saved components show the tick, unless their mask was edited since.
    const QString strSaved = QFileInfo(path, fi.completeBaseName() + "_mask.npy").absoluteFilePath();
    if (IsCompleteDataOutput(strSaved))
    {
      if (CCSavedIsStale(strSaved, MaskPathFor(fi.completeBaseName())))
        listStale << fns.size();
      else
        listSaved << fns.size();
    }
    fns << fi.absoluteFilePath();
  }
  foreach (int i, listStale)
    m_setCCUnsaved.insert(i);
  ui->listWidgetThumbnail->LoadImageFiles(fns, QStringList(), listSaved);
  // Set before CCSavedPath() below.
  m_strFinalOutputFolder = path.trimmed();
  // Saved photos show their components (built off the GUI thread, after the
  // plain thumbnails).
  foreach (int i, listSaved)
    ui->listWidgetThumbnail->SetItemComponents(i, fns[i], CCMaskPath(i), CCSavedPath(i), m_listStockColors,
                                               ui->horizontalSliderLabelOpacity->Value()/100.0);
  // Boxes edited on an earlier visit and left unsaved: still unsaved, so no tick.
  foreach (int i, listSaved)
  {
    if (m_listRegionData.contains(i) && m_listRegionData[i] != LoadCCBoxes(i))
    {
      m_setCCUnsaved.insert(i);
      ui->listWidgetThumbnail->SetItemChecked(i, false);
    }
  }

  // No "folder is not empty" prompt: the pass only writes scratch files, and only
  // for masks saved since their last pass. LoadImage() restores the boxes.
  RunOrWaitForMaskToCC();
}

void MainWindow::RunOrWaitForMaskToCC()
{
  const bool bOnCC = (ui->stackedWidget->currentWidget() == ui->pageCC);
  const bool bCanLoad = bOnCC && m_nIndex >= 0 && m_nIndex < m_listInputFiles.size();
  auto ShowCurrent = [&]() { if (bCanLoad) LoadImage(m_nIndex); UpdateIndex(); };
  // A pass is running, or just exited with its queued finish pending.
  // FinishCCPass() starts another pass for masks saved meanwhile.
  if (m_bCCPassStarted)
  {
    UpdateCCProgress();
    ShowCurrent();
    return;
  }
  // Nothing changed since the last pass (usually the Segmentation warm-up):
  // use its output straight away.
  const bool bWarmupRunning = m_procCCWarmup && m_procCCWarmup->state() != QProcess::NotRunning;
  if (!m_bCCPassActive && !bWarmupRunning && StaleCCMasks().isEmpty())
  {
    ui->widgetCCProgress->setVisible(false);
    ShowCurrent();
    return;
  }
  m_bCCPassActive = true;
  m_nCCNewest = -1;
  UpdateCCProgress();
  ShowCurrent();
  // Drop leftover wait-connections (Qt::SingleShotConnection is Qt 6 only).
  if (m_procCCWarmup)
    disconnect(m_procCCWarmup, SIGNAL(finished(int,QProcess::ExitStatus)),
               this, SLOT(RunOrWaitForMaskToCC()));
  if (bWarmupRunning)
  {
    // The Segmentation warm-up writes the same scratch folder; wait for it instead
    // of racing. What it finishes counts as fresh.
    ui->labelCCProgress->setText("Finishing the connected-components pass already running...");
    connect(m_procCCWarmup, SIGNAL(finished(int,QProcess::ExitStatus)),
            this, SLOT(RunOrWaitForMaskToCC()));
    return;
  }
  // m_proc runs one task at a time (a resumed project may still be running
  // nnUNet). Wait, queued, so OnProcessFinished() handles that task first.
  disconnect(m_proc, SIGNAL(finished(int,QProcess::ExitStatus)),
             this, SLOT(RunOrWaitForMaskToCC()));
  // A wait that ends after the user has left the page must not start the
  // pass there.
  if (!bOnCC)
  {
    m_bCCPassActive = false;
    ui->widgetCCProgress->setVisible(false);
    return;
  }
  if (m_proc->state() != QProcess::NotRunning)
  {
    ui->labelCCProgress->setText("Waiting for the step still running to finish...");
    connect(m_proc, SIGNAL(finished(int,QProcess::ExitStatus)),
            this, SLOT(RunOrWaitForMaskToCC()), Qt::QueuedConnection);
    return;
  }

  // Reprocess only stale masks, from a scratch copy of just those. Their old .npz
  // files are removed first so readiness can be read from which files exist; the
  // script writes them in sorted path order.
  QStringList listStale = StaleCCMasks();
  if (listStale.isEmpty())
  {
    // The warm-up that was running covered everything.
    m_bCCPassActive = false;
    ui->widgetCCProgress->setVisible(false);
    ShowCurrent();
    return;
  }
  listStale.sort(Qt::CaseSensitive);
  m_listCCOrder = listStale;
  const QString strPassIn = StageCCMasks(m_listCCOrder, "cc_in_pass", m_hashCCPassSig);
  m_nCCNewest = -1;
  m_bCCPassStarted = true;

  m_proc->setProperty("task", "mask_to_cc");
  m_proc->start(m_strPythonCmd, QStringList() << m_strPyScriptMaskToCC
                << "--in_dir" << strPassIn << "--out_dir" << m_sTempDir);
  // A start that failed straight away has already been through
  // OnProcessError() -> FinishCCPass(false).
  if (!m_bCCPassActive)
    return;
  if (!m_timerCCPoll)
  {
    m_timerCCPoll = new QTimer(this);
    m_timerCCPoll->setInterval(250);
    connect(m_timerCCPoll, &QTimer::timeout, this, &MainWindow::PollCCPass);
  }
  m_timerCCPoll->start();
  UpdateCCProgress();
  UpdateIndex();
}

// Masks are paired with photos by name, not position, so a missing or stray
// file cannot shift later photos onto the wrong mask.
QString MainWindow::CCMaskPath(int n) const
{
  if (n < 0 || n >= m_listInputFiles.size() || m_strMaskFolder.isEmpty())
    return QString();
  return QFileInfo(m_strMaskFolder, m_listInputFiles[n].completeBaseName() + ".png").absoluteFilePath();
}

QString MainWindow::CCNpzPath(const QString& strMask) const
{
  return QFileInfo(m_sTempDir, QFileInfo(strMask).completeBaseName() + ".npz").absoluteFilePath();
}

// A mask's identity as the pass read it: modification time and size, compared
// exactly (never against this machine's clock), plus an MD5 when the file system
// stores whole seconds.
MainWindow::CCMaskSig MainWindow::ReadCCMaskSig(const QString& strMask, bool bWithDigest)
{
  CCMaskSig sig;
  QFileInfo fi(strMask);
  if (!fi.isFile())
    return sig;
  sig.mtime = fi.lastModified().toMSecsSinceEpoch();
  sig.size = fi.size();
  if (bWithDigest)
  {
    QFile f(strMask);
    if (f.open(QIODevice::ReadOnly))
      sig.md5 = QCryptographicHash::hash(f.readAll(), QCryptographicHash::Md5);
  }
  return sig;
}

bool MainWindow::SameCCMaskSig(const CCMaskSig& recorded, const QString& strMask)
{
  if (recorded.size < 0)
    return false;
  const CCMaskSig now = ReadCCMaskSig(strMask, !recorded.md5.isEmpty());
  return now.size == recorded.size && now.mtime == recorded.mtime && now.md5 == recorded.md5;
}

bool MainWindow::IsCCReady(int n) const
{
  const QString strMask = CCMaskPath(n);
  if (strMask.isEmpty() || !QFileInfo(strMask).isFile())
    return false;
  if (m_bCCPassStarted)
  {
    const int nPos = m_listCCOrder.indexOf(strMask);
    if (nPos >= 0)
    {
      // In this pass: written, not the file being written now (a later one
      // exists), and the mask not saved again since it was staged.
      return nPos < m_nCCNewest && QFileInfo(CCNpzPath(strMask)).isFile()
          && SameCCMaskSig(m_hashCCPassSig.value(strMask), strMask);
    }
  }
  return IsCCMaskFresh(strMask);
}

bool MainWindow::IsCCMaskFresh(const QString& strMask) const
{
  if (!m_hashCCNpzFrom.contains(strMask) || !QFileInfo(CCNpzPath(strMask)).isFile())
    return false;
  return SameCCMaskSig(m_hashCCNpzFrom.value(strMask), strMask);
}

QStringList MainWindow::StaleCCMasks() const
{
  QStringList list;
  for (int i = 0; i < m_listInputFiles.size(); i++)
  {
    const QString strMask = CCMaskPath(i);
    if (!strMask.isEmpty() && QFileInfo(strMask).isFile() && !IsCCMaskFresh(strMask))
      list << strMask;
  }
  return list;
}

QString MainWindow::StageCCMasks(const QStringList& listMasks, const QString& strDirName,
                                 QHash<QString, CCMaskSig>& hashSig)
{
  const QString strDir = QFileInfo(m_sTempDir, strDirName).absoluteFilePath();
  QDir().mkpath(strDir);
  ClearFolder(strDir);
  hashSig.clear();
  foreach (const QString& strMask, listMasks)
  {
    // Identity read BEFORE the copy: a save between the two then makes the
    // result look stale (redone), never look fresh with old content.
    CCMaskSig sig = ReadCCMaskSig(strMask, false);
    if (sig.mtime % 1000 == 0)
      sig = ReadCCMaskSig(strMask, true);
    hashSig[strMask] = sig;
    // The old result goes, so an .npz after the run can only be this run's.
    QFile::remove(CCNpzPath(strMask));
    m_hashCCNpzFrom.remove(strMask);
    QFile::copy(strMask, QFileInfo(strDir, QFileInfo(strMask).fileName()).absoluteFilePath());
  }
  return strDir;
}

void MainWindow::PollCCPass()
{
  if (!m_bCCPassStarted)
  {
    if (m_timerCCPoll)
      m_timerCCPoll->stop();
    return;
  }
  int nNewest = -1;
  for (int i = 0; i < m_listCCOrder.size(); i++)
  {
    if (QFileInfo(CCNpzPath(m_listCCOrder[i])).isFile())
      nNewest = i;
  }
  if (nNewest == m_nCCNewest)
    return;
  m_nCCNewest = nNewest;
  UpdateCCProgress();
  if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    if (m_bCCViewPending && IsCCReady(m_nIndex) && m_nIndex < m_listInputFiles.size())
      LoadImage(m_nIndex);
    UpdateIndex();
  }
}

void MainWindow::UpdateCCProgress()
{
  if (!m_bCCPassActive)
  {
    ui->widgetCCProgress->setVisible(false);
    return;
  }
  ui->widgetCCProgress->setVisible(true);
  if (!m_bCCPassStarted)
  {
    ui->progressBarCC->setRange(0, 0);
    ui->labelCCProgress->setText("Waiting to find connected components...");
    return;
  }
  const int nTotal = m_listCCOrder.size();
  const int nDone = qBound(0, m_nCCNewest, nTotal);
  ui->progressBarCC->setRange(0, qMax(nTotal, 1));
  ui->progressBarCC->setValue(nDone);
  ui->labelCCProgress->setText(QString("Finding connected components in the background: %1 of %2 ready")
                               .arg(nDone).arg(nTotal));
}

void MainWindow::FinishCCPass(bool bSucceeded)
{
  if (!m_bCCPassActive)
    return;
  const bool bWasStarted = m_bCCPassStarted;
  // Masks this pass tried but got no result for, unchanged since: running
  // again would only fail again, so they are not retried automatically.
  QSet<QString> setNoResult;
  if (bWasStarted)
  {
    foreach (const QString& strMask, m_listCCOrder)
    {
      const bool bHasNpz = QFileInfo(CCNpzPath(strMask)).isFile();
      if (bSucceeded && bHasNpz)
        m_hashCCNpzFrom[strMask] = m_hashCCPassSig.value(strMask);
      else if (SameCCMaskSig(m_hashCCPassSig.value(strMask), strMask))
        setNoResult.insert(strMask);
    }
  }
  m_bCCPassActive = false;
  m_bCCPassStarted = false;
  if (m_timerCCPoll)
    m_timerCCPoll->stop();
  ui->widgetCCProgress->setVisible(false);
  if (ui->stackedWidget->currentWidget() != ui->pageCC)
    return;
  // Masks saved on Segmentation while this pass ran: go again for those.
  if (bWasStarted && bSucceeded)
  {
    bool bMore = false;
    foreach (const QString& strMask, StaleCCMasks())
      bMore = bMore || !setNoResult.contains(strMask);
    if (bMore)
    {
      QTimer::singleShot(0, this, SLOT(RunOrWaitForMaskToCC()));
      return;
    }
  }
  // Stay on the photo the user is on; only a photo still waiting is reloaded.
  if (m_bCCViewPending && m_nIndex >= 0 && m_nIndex < m_listInputFiles.size())
    LoadImage(m_nIndex);
  UpdateIndex();
}

void MainWindow::ShowProgressWindow(bool bShow)
{
  m_wndProgress->setVisible(bShow);
  if (bShow)
    RepositionProgressWindow();
}

void MainWindow::OnLinkActivated(const QString &link)
{
  if (sender() == ui->labelCalibrationFile)
  {
    SetCalibrationFile("");
  }
  else if (sender() == ui->labelImageFolder)
  {
    if (link == "clear")
      SetInputFolder("");
  }
}

void MainWindow::OnButtonChangeOutputFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", m_strRootOutputFolder);
  if (!path.isEmpty())
  {
    SetInputFolder(m_strInputFolder, path);
  }
}

void MainWindow::OnButtonCalibrationFile()
{
  QString path = QFileDialog::getOpenFileName(this, "Select Calibration File", m_strCurrentProjectFolder);
  if (!path.isEmpty())
  {
    SetCalibrationFile(path);
  }
}

void MainWindow::OnButtonCalibrationPhoto()
{
  QString path = QFileDialog::getOpenFileName(this, "Select Image File", m_strCurrentProjectFolder);
  if (!path.isEmpty())
  {
    SetCalibrationPhoto(path);
  }
}

void MainWindow::OnButtonBackToStart()
{
  StopSegDaemon();
  GoToStage(ui->pageStartUp);
  SetCalibrationPhoto("");
  ui->widgetImageView->Clear();
}

void MainWindow::OnButtonSaveCalibration()
{
  double w = 1, h = 1;
  bool bOK;
  w = ui->lineEditRectWidth->text().trimmed().toDouble(&bOK);
  if (bOK)
    h = ui->lineEditRectHeight->text().trimmed().toDouble(&bOK);
  if (!bOK)
  {
    Warn(this, "Marker size is missing",
                         "Enter the width and height of the printed marker square, in millimetres.");
    return;
  }
  QList<QPoint> pts = ui->widgetImageView->GetEditedPoints();
  if (pts.size() < 8)
  {
    Warn(this, "Marker outline is incomplete",
                         "Place all 4 markers and give each one a radius before saving.");
    return;
  }

  QString strDefaultDir = !m_strCurrentProjectFolder.isEmpty() ? m_strCurrentProjectFolder
                         : (!m_strRootOutputFolder.isEmpty() ? m_strRootOutputFolder
                         : QFileInfo(m_strCalibrationPhoto).absolutePath());
  QString strProjectName = !m_currentProjectInfo.name.isEmpty() ? m_currentProjectInfo.name
                          : (!m_strCurrentProjectFolder.isEmpty() ? QDir(m_strCurrentProjectFolder).dirName()
                          : QFileInfo(m_strCalibrationPhoto).completeBaseName());
  QString strDefaultName = strProjectName + "_calibration.npz";
  QString fn = QFileDialog::getSaveFileName(this, "Specify Output npz File",
                                             QDir(strDefaultDir).filePath(strDefaultName), "*.npz");
  if (!fn.isEmpty())
  {
    QStringList list;
    foreach (QPoint pt, pts)
      list << QString::number(pt.x()) << QString::number(pt.y());

    m_proc->start(m_strPythonCmd, QStringList() << m_strPyScriptFiducialsCalibration
                  << "--in_img" << ui->widgetImageView->GetFilename()
                  << "--points" << list
                  << "--width" << QString::number(w)
                  << "--height" << QString::number(h)
                  << "--out_file" << fn);
    m_proc->setProperty("task", "fiducials_calibration");
    ui->widgetImageView->ShowMessage("Processing...");
  }
}

void MainWindow::OnPointEdited()
{
  QList<QPoint> pts = ui->widgetImageView->GetEditedPoints();
  if (ui->pageCalibration->isVisible())
  {
    ui->pushButtonSaveCalibration->setEnabled(pts.size() >= 8);
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
    UpdateCorrClearEnabled();
}

void MainWindow::OnButtonLoadMask()
{
  QString path = m_strRootOutputFolder + "/mask";
  QDir().mkpath(path);
  if (!path.isEmpty())
  {
    m_strMaskFolder = path;
    m_listMaskFiles = QDir(path).entryInfoList(QDir::Files, QDir::Name);
    UpdateSegContinueCaption();
    m_listInputFiles = CorrectedPhotos();
    ui->listWidgetThumbnail->Clear();
    QStringList fns, mask_fns;
    foreach (QFileInfo fi, m_listInputFiles)
    {
      fns << fi.absoluteFilePath();
      mask_fns << QFileInfo(path, fi.completeBaseName() + ".png").absoluteFilePath();
    }
    // Missing masks are skipped; they are added as they arrive.
    ui->listWidgetThumbnail->LoadImageFiles(fns, mask_fns);

    if (!m_listInputFiles.isEmpty())
    {
      // Off in Automatic segmentation settings: no model, no questions; masks that
      // are missing are drawn by hand.
      bool bRunNNUnet = DialogSegmentationSettings::Enabled();
      bool bResumeRemainingOnly = false;
      // When resuming straight to Connected components, this page is only passed
      // through: keep the masks without asking. A resume that stops here still asks.
      bool bSkipPromptForResume = m_bResumingProject && m_currentProjectInfo.currentStage == "connected_components";
      if (bRunNNUnet && m_listInputFiles.size() == m_listMaskFiles.size())
      {
        if (bSkipPromptForResume)
        {
          bRunNNUnet = false;
        }
        else
        {
          // "Keep masks" is button 0 (default and escape). Only a click on "Start again"
          // clears the masks.
          if (ShowBox(this, "Masks already exist",
                      "This project already has a mask for every photo. Keep them, or start again with "
                      "fresh automatic masks. Starting again replaces every mask, including any you edited.",
                      QString(),
                      QStringList() << "Keep masks" << "Start again", 1) != 1)
            bRunNNUnet = false;
        }
      }
      else if (bRunNNUnet && !m_listMaskFiles.isEmpty())
      {
        // Some masks exist, usually from an interrupted run.
        if (bSkipPromptForResume)
        {
          bResumeRemainingOnly = true;
        }
        else
        {
          // Same rule: anything but a click on "Start again" keeps them.
          int nRemaining = m_listInputFiles.size() - m_listMaskFiles.size();
          bResumeRemainingOnly =
              (ShowBox(this, "Some masks exist",
                       QString("%1 of %2 masks were made in an earlier session. Finish the remaining %3, "
                               "or start again and replace all of them.")
                           .arg(m_listMaskFiles.size()).arg(m_listInputFiles.size()).arg(nRemaining),
                       QString(),
                       QStringList() << "Finish remaining" << "Start again", 1) != 1);
        }
      }
      if (bRunNNUnet && (m_strNNUnetScriptFolder.isEmpty() || m_strNNUnetModelFolder.isEmpty()))
      {
        Warn(this, "Segmentation model was not found",
                               "Masks cannot be generated automatically. You can still draw them by hand.\n\n"
                               "Choose the model folder in Automatic segmentation settings (the gear on the "
                               "Projects screen), or set NNUNET_SCRIPT_DIR and NNUNET_MODEL_DIR, then try again.");
        // On to drawing by hand.
        bRunNNUnet = false;
      }
      if (bRunNNUnet)
      {
        // A running CC pass is moot once masks are regenerated: stop it without its
        // finish being handled. Anything else must finish first.
        if (m_proc->state() != QProcess::NotRunning)
        {
          if (m_proc->property("task").toString() != "mask_to_cc")
          {
            Warn(this, "Another step is still running",
                 "Wait for it to finish, then try again.");
            return;
          }
          m_proc->blockSignals(true);
          m_proc->kill();
          m_proc->waitForFinished(3000);
          m_proc->blockSignals(false);
          FinishCCPass(false);
        }

        QSet<QString> existingMaskBaseNames;
        if (bResumeRemainingOnly)
        {
          foreach (const QFileInfo& fi, m_listMaskFiles)
            existingMaskBaseNames.insert(fi.completeBaseName());
        }
        else
        {
          ClearFolder(path);
          // Every mask is new: none has been previewed, and no thumbnail shows one.
          m_setSegPreviewed.clear();
          ui->listWidgetThumbnail->LoadImageFiles(fns);
        }

        QList<int> listNeeded;
        for (int i = 0; i < m_listInputFiles.size(); i++)
        {
          if (!existingMaskBaseNames.contains(m_listInputFiles[i].completeBaseName()))
            listNeeded << i;
        }
        if (listNeeded.isEmpty())
          bRunNNUnet = false;
        else
        {
          // Usually the daemon started on Correction is still running and has
          // segmented most photos already; otherwise start one now.
          StartSegDaemon();
          QStringList watch_list, listUnqueued;
          m_mapInputImageSize.clear();
          foreach (int i, listNeeded)
          {
            const QFileInfo& fi = m_listInputFiles[i];
            const QString strNpz = QFileInfo(m_strSegOutDir, fi.completeBaseName()+".npz").absoluteFilePath();
            // Neither done nor waiting (never queued, or failed: retried). Not watched if it
            // could not be queued: nothing would come, and the results behind it would wait
            // forever.
            if (!QFileInfo::exists(strNpz)
                && !QFileInfo::exists(QFileInfo(m_strSegInDir, fi.fileName()+".txt").absoluteFilePath())
                && !QueueForSegmentation(fi.absoluteFilePath(), false))
            {
              listUnqueued << fi.completeBaseName();
              continue;
            }
            watch_list << strNpz;
            // Read the size from the header instead of decoding each photo on the GUI
            // thread. autoTransform is off, so this is the raw pre-EXIF size NpyToImage()
            // expects. Fall back to a full decode if the header gives none.
            QImageReader szReader(fi.absoluteFilePath());
            QSize szInput = szReader.size();
            if (!szInput.isValid())
              szInput = QImage(fi.absoluteFilePath()).size();
            m_mapInputImageSize[watch_list.last()] = szInput;
          }
          if (watch_list.isEmpty())
          {
            // None could be queued: they are drawn by hand.
            bRunNNUnet = false;
            Warn(this, "Some photos were not segmented",
                 QString("%1 photo(s) could not be segmented automatically: %2.\n\n"
                         "Draw their masks by hand, or try again.").arg(listUnqueued.size()).arg(listUnqueued.join(", ")));
          }
          else
          {
            FinishSegQueue();
            ui->pageSegEdit->setEnabled(false);
            m_fileWatcher.addPath(m_strSegOutDir);
            m_listQueuedFiles = watch_list;
            m_nTotalSegFiles = watch_list.size();
            m_bSegFirstMaskShown = false;
            // Reported with the failures when the daemon exits.
            m_listSegFailed = listUnqueued;
            m_nIndex = 0;
            // The overlay covers the model load until the first mask arrives
            // (OnFileChanged() hides it); then the inline bar shows progress.
            m_wndProgress->SetIndeterminate("Loading segmentation model...");
            ShowProgressWindow(true);
            ui->widgetSegProgress->setVisible(true);
            ui->progressBarSeg->setMaximum(qMax(m_nTotalSegFiles, 1));
            ui->progressBarSeg->setValue(0);
            ui->labelSegProgress->setText(QString("Segmenting images in the background: 0 of %1 ready").arg(m_nTotalSegFiles));
            // Results may be there already.
            OnFileChanged(m_strSegOutDir);
          }
        }
      }
      if (!bRunNNUnet)
      {
        // Nothing to wait for: every photo keeps its mask, or none could be queued.
        StopSegDaemon();
        m_bSegFirstMaskShown = true;
        LoadImage(0);
        ui->widgetSegCtrls->setEnabled(true);
        UpdateIndex();
        // All masks exist: start the CC warm-up now.
        MaybeWarmUpConnectedComponents();
      }
    }
  }
}

bool MainWindow::StartSegDaemon()
{
  if (m_procSeg->state() != QProcess::NotRunning)
    return true;
  // Off in Automatic segmentation settings: every mask is drawn by hand.
  if (!DialogSegmentationSettings::Enabled() || m_strNNUnetScriptFolder.isEmpty() || m_strNNUnetModelFolder.isEmpty())
    return false;
  StopSegDaemon();
  // Not m_sTempDir itself: the CC pass writes its .npz files there.
  const qint64 nMs = QDateTime::currentMSecsSinceEpoch();
  m_strSegInDir = QString("%1/seg_in_%2").arg(m_sTempDir).arg(nMs);
  m_strSegOutDir = QString("%1/seg_out_%2").arg(m_sTempDir).arg(nMs);
  QDir().mkpath(m_strSegInDir);
  QDir().mkpath(m_strSegOutDir);
  m_strNNUnetStderr.clear();
  m_elapsedTimer.start();
  m_procSeg->start(m_strPythonCmd, QStringList()
                   << m_strNNUnetScriptFolder + "/process_directory_no_upsample_daemon.py"
                   << "--in_dir" << m_strSegInDir
                   << "--out_dir" << m_strSegOutDir
                   << "--model_path" << m_strNNUnetModelFolder
                   << "--device" << DialogSegmentationSettings::Device()
                   << "--folds" << DialogSegmentationSettings::Folds()
                   << "--max_threads" << QString::number(DialogSegmentationSettings::MaxThreads()));
  return true;
}

void MainWindow::SegmentCorrectedPhotos()
{
  // A daemon kept from Segmentation must not exit while photos may be
  // re-corrected.
  if (m_procSeg->state() != QProcess::NotRunning)
    QFile::remove(m_strSegInDir + "/ALL_DONE.txt");
  else if (StartSegDaemon())
  {
    foreach (const QFileInfo& fi, m_listInputFiles)
    {
      const QString strImage = GetDeformedImagePath(fi);
      if (!QFileInfo::exists(MaskPathFor(QFileInfo(strImage).completeBaseName())))
        QueueForSegmentation(strImage, false);
    }
  }
}

QString MainWindow::MaskPathFor(const QString& strBase) const
{
  return m_strRootOutputFolder + "/mask/" + strBase + ".png";
}

bool MainWindow::QueueForSegmentation(const QString& strImage, bool bRecorrected)
{
  if (!IsCompleteImageOutput(strImage))
    return false;
  const QFileInfo fi(strImage);
  const QString strBase = fi.completeBaseName();
  if (bRecorrected)
  {
    QFile::remove(MaskPathFor(strBase));
    QFile::remove(m_strRootOutputFolder + "/connected_components/" + strBase + "_mask.npy");
  }
  const QString strNpz = QFileInfo(m_strSegOutDir, strBase + ".npz").absoluteFilePath();
  const QString strFailed = QFileInfo(m_strSegOutDir, strBase + ".failed").absoluteFilePath();
  if (m_procSeg->state() == QProcess::NotRunning)
  {
    if (!m_strSegOutDir.isEmpty())
    {
      QFile::remove(strNpz);
      QFile::remove(strFailed);
    }
    return false;
  }
  // Copy, then rename into place (atomic), so the daemon never reads a partial
  // image. The old result is removed only after the rename: a job still
  // running on the old image checks it before publishing, so it cannot leave
  // a stale result behind.
  const QString strDest = QFileInfo(m_strSegInDir, fi.fileName()).absoluteFilePath();
  QFile::remove(strDest + ".part");
  if (!QFile::copy(strImage, strDest + ".part")
      || ::rename(QFile::encodeName(strDest + ".part").constData(), QFile::encodeName(strDest).constData()) != 0)
  {
    qDebug() << "Could not queue" << strImage << "for segmentation";
    return false;
  }
  QFile::remove(strNpz);
  QFile::remove(strFailed);
  QFile call(strDest + ".txt");
  if (!call.open(QIODevice::WriteOnly))
  {
    qDebug() << "Could not write" << call.fileName();
    return false;
  }
  return true;
}

void MainWindow::FinishSegQueue()
{
  QFile call(m_strSegInDir + "/ALL_DONE.txt");
  if (!call.open(QIODevice::WriteOnly))
    qDebug() << "Could not write" << call.fileName();
}

void MainWindow::ForgetSegProgress()
{
  if (m_listQueuedFiles.isEmpty())
    return;
  m_listQueuedFiles.clear();
  ui->widgetSegProgress->setVisible(false);
  ShowProgressWindow(false);
  ui->pageSegEdit->setEnabled(true);
}

bool MainWindow::SegResultReady() const
{
  if (m_listQueuedFiles.isEmpty())
    return false;
  const QFileInfo fiNpz(m_listQueuedFiles.first());
  return fiNpz.exists() || QFileInfo::exists(fiNpz.path() + "/" + fiNpz.completeBaseName() + ".failed");
}

void MainWindow::StopSegDaemon()
{
  if (m_procSeg->state() != QProcess::NotRunning)
  {
    m_procSeg->blockSignals(true);
    m_procSeg->kill();
    m_procSeg->waitForFinished(3000);
    m_procSeg->blockSignals(false);
  }
  if (!m_strSegOutDir.isEmpty())
  {
    m_fileWatcher.removePath(m_strSegOutDir);
    QDir(m_strSegInDir).removeRecursively();
    QDir(m_strSegOutDir).removeRecursively();
  }
  m_strSegInDir.clear();
  m_strSegOutDir.clear();
  m_bSegDaemonExited = false;
  // Masks still expected from it can no longer arrive.
  ForgetSegProgress();
}

void MainWindow::OnSegDaemonFinished(QProcess::ExitStatus status, int nCode)
{
  // A queued report from a daemon already stopped or replaced.
  if (m_procSeg->state() != QProcess::NotRunning || m_strSegOutDir.isEmpty())
    return;
  qDebug() << "Segmentation daemon finished after" << m_elapsedTimer.elapsed()/1000 << "s, exit code" << nCode;
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    // Nothing waits for it here; Segmentation starts a new one. After a normal
    // exit (it may have seen ALL_DONE just before returning here removed it),
    // replace it; never after a crash, to avoid a restart loop.
    qDebug() << m_strNNUnetStderr.right(2000);
    StopSegDaemon();
    if (status == QProcess::NormalExit)
      SegmentCorrectedPhotos();
    return;
  }
  m_bSegDaemonExited = true;
  m_bSegDaemonOk = (status == QProcess::NormalExit && nCode == 0);
  // Takes what is left one per turn, then calls ReportSegDaemonExit().
  OnFileChanged(m_strSegOutDir);
}

void MainWindow::ReportSegDaemonExit()
{
  // Photos with no mask: failed, or never reached because it stopped early.
  QStringList listMissing = m_listSegFailed;
  foreach (const QString& strNpz, m_listQueuedFiles)
    listMissing << QFileInfo(strNpz).completeBaseName();
  const bool bModelFailed = !m_bSegDaemonOk && m_listSegFailed.isEmpty();
  const QString strDetail = m_strNNUnetStderr.trimmed().right(1500);
  // Clean up before the dialog, whose event loop can deliver more changes.
  StopSegDaemon();
  UpdateIndex();
  UpdateStepper(ui->stackedWidget->currentWidget());
  if (bModelFailed)
    Warn(this, "Segmentation did not finish",
         QString("The automatic segmentation stopped, for example because the model could not load. "
                 "You can draw the masks by hand, or try again.") +
         (strDetail.isEmpty() ? QString() : "\n\n" + strDetail));
  else if (!listMissing.isEmpty())
    Warn(this, "Some photos were not segmented",
         QString("%1 photo(s) could not be segmented automatically: %2.\n\n"
                 "Draw their masks by hand, or try again.").arg(listMissing.size()).arg(listMissing.join(", ")));
}

void MainWindow::OnButtonContinue(bool bBack)
{
  // Serves Continue from the Start page and Back from Segmentation. Going back
  // skips the "already corrected" tip and the stage write.
  const bool bNavigatingBack = bBack || (sender() == ui->pushButtonBackToCorrection);

  if (!m_strPendingSourceFolder.isEmpty() && !CopyPendingPhotosIntoProject())
    return;

  m_strOutputFolder = m_strRootOutputFolder + "/deformed";
  if (!QFileInfo::exists(m_strOutputFolder))
    QDir().mkpath(m_strOutputFolder);

  if (m_strInputFolder.isEmpty() || !QDir(m_strInputFolder).exists())
  {
    Warn(this, "Photo folder is missing",
                         "Choose the folder that holds the photos, then continue.");
    return;
  }
  if (m_strOutputFolder.isEmpty() || !QDir(m_strOutputFolder).exists())
  {
    Warn(this, "Output folder is missing",
                         "Choose where the corrected photos should be written, then continue.");
    return;
  }

  if (!m_strCalibrationFile.isEmpty() && !QFile(m_strCalibrationFile).exists())
  {
    Warn(this, "Calibration file is missing",
                         "It is no longer where it was. Choose it again, or continue without one.");
    return;
  }
  m_listInputFiles = QDir(m_strInputFolder).entryInfoList(ImageNameFilters(), QDir::Files, QDir::Name);
  if (m_listInputFiles.isEmpty())
  {
    Warn(this, "No photos in that folder",
                         "Choose a folder that holds photos.\n\nLooked in " + m_strInputFolder + ".");
    return;
  }

  GoToStage(ui->pageCorrection);
  SetRecenterAccent("#a2293c", "#8e2434", "#7b1f2d");
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_POINT);
  ui->widgetImageView->SetNumberOfExpectedPoints(4);
  ui->widgetImageView->SetMaskOpacity(ui->sliderSegOpacity->Value()/100.0);
  // Loading re-shows the mask, so Show mask must be checked too.
  ui->pushButtonToggleMaskVisibility->setChecked(true);

  m_bCalibratedMode = !m_strCalibrationFile.isEmpty();
  ui->labelTitle->setText(m_bCalibratedMode?"Calibrated":"Retrospective Correction");
  // Same 16px markup as the .ui text, worded for the mode.
  ui->label_2->setText(QString("<html><head/><body><p style=\" font-size:16px;\">%1</p></body></html>").arg(
      m_bCalibratedMode
      ? "Get photographs in the same pixel scale. Review the pixel corrected images, and re-correct images "
        "when necessary by selecting the center of each fiducial marker in the original image."
      : "Get photographs in the same pixel scale. Choose the correction mode above, provide the real "
        "dimensions, and select the key points in your image."));
  ui->widgetPointMode->setVisible(!m_bCalibratedMode);

  ui->widgetIndex->show();
  ui->widgetIndexSeg->hide();
  ui->widgetIndexCC->hide();

  // Back from Segmentation: masks still wanted are collected again there.
  ForgetSegProgress();
  // LoadProject() starts it after its replay, if it stops here.
  if (!m_bResumingProject)
    SegmentCorrectedPhotos();

  QStringList fns;
  QList<int> checkedIndices;
  for (int i = 0; i < m_listInputFiles.size(); i++)
  {
    fns << m_listInputFiles[i].absoluteFilePath();
    // Green "corrected" border once the corrected output is complete.
    if (IsCompleteImageOutput(GetDeformedImagePath(m_listInputFiles[i])))
      checkedIndices << i;
  }
  ui->listWidgetThumbnail->LoadImageFiles(fns, QStringList(), checkedIndices);

  // Before UpdateIndex(), whose caption update would otherwise find the
  // cached count stale and read every header a second time.
  RecountCorrected();
  m_bCalibratedRunTried = false;
  if (!m_listInputFiles.isEmpty())
  {
    UpdatePreviewAvailability();
    LoadImage(m_nIndex, false, PreviewChecked());
  }
  UpdateIndex();

  ui->pushButtonNext->setEnabled(false);

  if (AllImagesProcessed())
  {
    ui->pushButtonSegmentation->setEnabled(true);
    UpdateCorrContinueCaption();
    if (!m_bResumingProject && !bNavigatingBack)
      Warn(this->ui->widgetControl, "These photos are already corrected",
                               "Every photo in this folder already has a corrected version. Go straight to "
                               "Segmentation, or redo any photo you want to change.");
  }

  // Stage is a high-water mark: don't rewind it when stepping back.
  if (!bNavigatingBack)
    SaveProjectState("correction");
}

void MainWindow::OnPointModeChanged(int index)
{
  bool b2Points = (index == 0);
  ui->labelRectWidth->setText(b2Points?"Ruler width (mm)":"Box width (mm)");

  // Hiding "Box height" shrinks the top bar's size hint, which makes the splitter
  // re-divide its panes. Keep the user's divider position across the toggle;
  // QSplitter clamps restored sizes to the new minimums.
  const QList<int> listSplitterSizes = ui->splitter->sizes();
  ui->labelRectHeight->setVisible(!b2Points);
  ui->lineEditRectHeight->setVisible(!b2Points);
  ui->splitter->setSizes(listSplitterSizes);
  ui->widgetImageView->SetNumberOfExpectedPoints(b2Points?2:(index == 1?3:4));
}

void MainWindow::SetupScriptPath()
{
  m_strNNUnetScriptFolder = QProcessEnvironment::systemEnvironment().value("NNUNET_SCRIPT_DIR");
  // From Automatic segmentation settings (default: NNUNET_MODEL_DIR, else FreeSurfer's).
  m_strNNUnetModelFolder = DialogSegmentationSettings::ModelDir();
  if (!DialogSegmentationSettings::ModelFound(m_strNNUnetModelFolder))
  {
    qDebug() << "No nnUNet model in" << m_strNNUnetModelFolder;
    m_strNNUnetModelFolder.clear();
  }
  if (m_strNNUnetScriptFolder.isEmpty())
  {
    m_strNNUnetScriptFolder = QProcessEnvironment::systemEnvironment().value("FREESURFER_HOME_FSPYTHON") + "/python/packages/nnUNet_v2";
    if (!QFileInfo(m_strNNUnetScriptFolder).exists())
    {
      m_strNNUnetScriptFolder.clear();
      qDebug() << "NNUNET_SCRIPT_DIR is not set!";
    }
  }

  static QTemporaryDir dir;
  m_sTempDir = dir.path();
  {
    QFile file(QString(":/")+SCRIPT_RETROSPECTIVE);
    if (!file.open(QFile::ReadOnly | QFile::Text))
      qDebug() << "Could not read" << file.fileName();
    QTextStream in(&file);
    QString str = in.readAll();
    str.replace("./resources/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
    str.replace("./resources/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
    file.close();
    file.setFileName(QFileInfo(dir.path(),SCRIPT_RETROSPECTIVE).absoluteFilePath());
    if (file.open(QFile::ReadWrite))
    {
      QTextStream out(&file);
      out << str;
    }
    file.close();
    QFile::copy(":/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
    QFile::copy(":/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
    m_strPyScriptRetrospective = file.fileName();

    if (!QFile::exists(m_strPyScriptRetrospective))
    {
      Warn(this, "Correction step is not installed",
                            tr("Dissection Photo could not find the script that corrects the photos.\n\n"
                               "Looked for %1. Check that FreeSurfer is set up in the terminal this "
                               "was started from.").arg(SCRIPT_RETROSPECTIVE));
      QTimer::singleShot(0, qApp, SLOT(quit()));
    }
  }
  {
    QFile::copy(":/registration.py", QFileInfo(dir.path(),"registration.py").absoluteFilePath());
    m_strPyScriptFiducialsDetection = QFileInfo(dir.path(),"func_fiducials_detection.py").absoluteFilePath();
    QFile::copy(":/func_fiducials_detection.py", m_strPyScriptFiducialsDetection);
    m_strPyScriptFiducialsCalibration = QFileInfo(dir.path(),"func_fiducials_calibration.py").absoluteFilePath();
    QFile::copy(":/func_fiducials_calibration.py", m_strPyScriptFiducialsCalibration);
  }
  {
    m_listStockColors << QColor(255,100,100) << QColor(255,255,100) << QColor(100,255,100)
                      << QColor(110,245,255) << QColor(75,100,255) << QColor(255,128,0)
                      << QColor(100,150,170) << QColor(120,60,220);
    m_strPyScriptMaskToCC = QFileInfo(dir.path(),"func_mask_to_cc.py").absoluteFilePath();
    QFile::copy(":/func_mask_to_cc.py", m_strPyScriptMaskToCC);
  }
}

void MainWindow::MaybeWarmUpConnectedComponents()
{
  if (m_strPyScriptMaskToCC.isEmpty() || m_strMaskFolder.isEmpty())
    return;
  // Never on the CC page: the pass overwrites the scratch .npz files that page
  // reads.
  if (ui->stackedWidget->currentWidget() != ui->pageSegEdit)
    return;
  // Skip if a previous background run is still going -- it'll pick up
  // this mask (and any others changed since) the next time this is called.
  if (m_procCCWarmup && m_procCCWarmup->state() != QProcess::NotRunning)
    return;
  // A page pass still writing the same scratch files.
  if (m_bCCPassActive)
    return;
  // Only masks saved since their last pass; usually just the one edited.
  const QStringList listStale = StaleCCMasks();
  if (listStale.isEmpty())
    return;

  if (!m_procCCWarmup)
  {
    m_procCCWarmup = new QProcess(this);
    // Connected before any RunOrWaitForMaskToCC() wait on this process, so
    // the results are recorded before that re-entry looks for stale masks.
    connect(m_procCCWarmup, static_cast<void (QProcess::*)(int, QProcess::ExitStatus)>(&QProcess::finished),
            this, [this](int nCode, QProcess::ExitStatus status) {
      if (status != QProcess::NormalExit || nCode != 0)
        return;
      // StageCCMasks() removed the old .npz, so one present now is this run's.
      for (QHash<QString, CCMaskSig>::const_iterator it = m_hashCCWarmSig.constBegin();
           it != m_hashCCWarmSig.constEnd(); ++it)
      {
        if (QFileInfo(CCNpzPath(it.key())).isFile())
          m_hashCCNpzFrom[it.key()] = it.value();
      }
    });
  }
  const QString strWarmIn = StageCCMasks(listStale, "cc_in_warm", m_hashCCWarmSig);
  m_procCCWarmup->start(m_strPythonCmd, QStringList() << m_strPyScriptMaskToCC
                        << "--in_dir" << strWarmIn << "--out_dir" << m_sTempDir);
}

void MainWindow::UpdateThumbnailWithMask(int nIndex)
{
  if (nIndex < 0 || nIndex >= m_listInputFiles.size())
    return;
  // Built off the GUI thread from the saved mask, like the first thumbnails.
  ui->listWidgetThumbnail->SetItemFiles(nIndex, m_listInputFiles[nIndex].absoluteFilePath(), CCMaskPath(nIndex));
}

bool MainWindow::SegMaskPending(int n) const
{
  return n >= 0 && n < m_listInputFiles.size() && !m_strSegOutDir.isEmpty()
      && m_listQueuedFiles.contains(QFileInfo(m_strSegOutDir, m_listInputFiles[n].completeBaseName() + ".npz").absoluteFilePath());
}

void MainWindow::UpdateSegReadiness()
{
  // Photos are numbered from 1, as on the carousel badges.
  QStringList listNotPreviewed, listToDraw;
  int nPending = 0;
  for (int i = 0; i < m_listInputFiles.size(); i++)
  {
    // The mask being drawn on screen counts: continuing saves it.
    const bool bMask = QFileInfo(CCMaskPath(i)).isFile()
        || (i == m_nIndex && ui->widgetImageView->IsMaskEdited());
    if (SegMaskPending(i))
      nPending++;
    else if (!bMask)
      listToDraw << QString::number(i+1);
    else if (!m_setSegPreviewed.contains(i))
      listNotPreviewed << QString::number(i+1);
  }
  const bool bReady = !m_listInputFiles.isEmpty() && nPending == 0 && listToDraw.isEmpty() && listNotPreviewed.isEmpty();
  ui->pushButtonCC->setEnabled(bReady);
  if (!m_pSegPreviewNote)
    return;

  auto numbers = [](const QStringList& l)
  {
    const int nShown = 8;
    return l.size() <= nShown ? l.join(", ") : QStringList(l.mid(0, nShown)).join(", ") + QString(" and %1 more").arg(l.size() - nShown);
  };
  QString str = "<b>Preview every photo's mask before moving on to Connected components.</b>";
  if (bReady)
    str += "<br>All masks have been previewed.";
  if (!listNotPreviewed.isEmpty())
    str += QString("<br>Not previewed yet: %1.").arg(numbers(listNotPreviewed));
  if (!listToDraw.isEmpty())
    str += DialogSegmentationSettings::Enabled()
        ? QString("<br>Draw these masks by hand (no automatic mask): %1.").arg(numbers(listToDraw))
        : QString("<br>Automatic segmentation is off: draw these masks by hand: %1.").arg(numbers(listToDraw));
  if (nPending > 0)
    str += QString("<br>%1 mask(s) still being segmented.").arg(nPending);
  // Font size repeated in the markup: rich text ignores the widget stylesheet.
  m_pSegPreviewNote->setText("<span style=\"font-size:14px;\">" + str + "</span>");
}

bool MainWindow::AllImagesProcessed() const
{
  return AllImagesProcessed(m_listInputFiles);
}

bool MainWindow::AllImagesProcessed(const QFileInfoList& originals) const
{
  if (originals.isEmpty())
    return false;
  // Check each expected output, not a file count (stray or half-written files
  // would pass). Cheap count pre-filter first, since this runs on every
  // navigation and network home folders are slow.
  if (QDir(m_strOutputFolder).entryInfoList(ImageNameFilters(), QDir::Files).size()
      < originals.size())
    return false;
  foreach (const QFileInfo& fi, originals)
  {
    if (!IsCompleteImageOutput(GetDeformedImagePath(fi)))
      return false;
  }
  return true;
}

bool MainWindow::AllCCMasksSaved() const
{
  if (m_listInputFiles.isEmpty() || m_strFinalOutputFolder.isEmpty())
    return false;
  foreach (const QFileInfo& fi, m_listInputFiles)
  {
    if (!IsCompleteDataOutput(QFileInfo(m_strFinalOutputFolder,
                                        fi.completeBaseName() + "_mask.npy").absoluteFilePath()))
      return false;
  }
  return true;
}

QString MainWindow::GetDeformedImagePath(const QFileInfo& originalFile) const
{
  // Same naming as func_retrospective_correction.py: split on the last dot.
  return QFileInfo(m_strOutputFolder, originalFile.completeBaseName() + "_deformed." + originalFile.suffix()).absoluteFilePath();
}

QFileInfoList MainWindow::CorrectedPhotos() const
{
  // Not a sort of the corrected names: "_deformed" sorts differently from the
  // extension ("a_1_deformed" comes after "a_10_deformed", "a_1.jpg" before "a_10.jpg").
  QFileInfoList list;
  foreach (const QFileInfo& fi, QDir(m_strInputFolder).entryInfoList(ImageNameFilters(), QDir::Files, QDir::Name))
  {
    QFileInfo corrected(GetDeformedImagePath(fi));
    if (corrected.isFile())
      list << corrected;
  }
  return list;
}

void MainWindow::UpdateIndex()
{
  QLabel* label = ui->labelIndex;
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    ui->pushButtonPrevious->setEnabled(m_nIndex > 0);
    // Next unlocks once "Correct this photo" has been started for this image, or it
    // is already corrected. Calibrated mode corrects every photo in one run.
    bool bCanAdvance = m_bCalibratedMode || m_setProcessStarted.contains(m_nIndex)
                     || (m_nIndex >= 0 && m_nIndex < m_listInputFiles.size()
                         && IsCompleteImageOutput(GetDeformedImagePath(m_listInputFiles[m_nIndex])));
    ui->pushButtonNext->setEnabled(bCanAdvance && m_nIndex < m_listInputFiles.size()-1);
    ui->pushButtonSegmentation->setEnabled(AllImagesProcessed());
    UpdateCorrClearEnabled();
    UpdateCorrContinueCaption();
  }
  if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    label = ui->labelIndexSeg;
    m_listMaskFiles = QDir(m_strMaskFolder).entryInfoList(QDir::Files, QDir::Name);
    UpdateSegContinueCaption();
    ui->pushButtonPreviousSeg->setEnabled(m_nIndex > 0);
    ui->pushButtonNextSeg->setEnabled(m_nIndex < m_listInputFiles.size()-1);
    // The photo on screen counts as previewed once its mask is there.
    if (m_nIndex >= 0 && m_nIndex < m_listInputFiles.size() && !SegMaskPending(m_nIndex))
      m_setSegPreviewed.insert(m_nIndex);
    UpdateSegReadiness();
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    label = ui->labelIndexCC;
    ui->pushButtonPreviousCC->setEnabled(m_nIndex > 0);
    // Next unlocks once this photo's components are saved, or it has boxes on
    // screen (leaving the photo saves them).
    const bool bBoxes = !m_bCCViewPending && IsCCReady(m_nIndex) && !ui->widgetImageView->GetEditedRegions().isEmpty();
    bool bCCCanAdvance = m_setCCProcessed.contains(m_nIndex) || bBoxes
                       || (m_nIndex >= 0 && m_nIndex < m_listInputFiles.size()
                           && IsCompleteDataOutput(QFileInfo(m_strFinalOutputFolder,
                               m_listInputFiles[m_nIndex].completeBaseName() + "_mask.npy").absoluteFilePath()));
    ui->pushButtonNextCC->setEnabled(bCCCanAdvance && m_nIndex < m_listInputFiles.size()-1);
    ui->pushButtonQuit->setEnabled(CCReadyForRecon());
    ui->pushButtonCreateMask->setEnabled(bBoxes);
  }
  ui->listWidgetThumbnail->blockSignals(true);
  ui->listWidgetThumbnail->SetCurrentIndex(m_nIndex);
  ui->listWidgetThumbnail->blockSignals(false);
  label->setText(tr("%1 of %2").arg(m_nIndex+1).arg(m_listInputFiles.size()));
}

void MainWindow::LoadImage(int n, bool bNoProcessing, bool bPreview)
{
  QList<QPoint> pts;
  QList<RECT_REGION> rects;
  bool bNoPreviewYet = false;
  QString fn = m_listInputFiles[n].absoluteFilePath(), mask_fn;
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    if (m_listPointData.contains(n))
      pts = m_listPointData[n];
    if (bPreview)
    {
      QString fn_p = GetDeformedImagePath(QFileInfo(fn));
      if (QFile::exists(fn_p))
      {
        fn = fn_p;
        pts.clear();
        ui->widgetImageView->HideMessage();
      }
      else
      {
        // Said after the load below, which clears messages.
        bNoPreviewYet = true;
      }
    }
    else
    {
      ui->widgetImageView->HideMessage();
    }
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    m_listMaskFiles = QDir(m_strMaskFolder).entryInfoList(QDir::Files, QDir::Name);
    UpdateSegContinueCaption();
    // By name, not position. A photo with no mask gets an empty one to draw
    // on, unless its mask is still coming from the daemon.
    mask_fn = CCMaskPath(n);
    if (!QFileInfo(mask_fn).isFile() && SegMaskPending(n))
      mask_fn.clear();
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    if (QFileInfo(CCMaskPath(n)).isFile())
      mask_fn = CCMaskPath(n);
    // Boxes saved with the components (earlier visit or session) come back as drawn.
    if (!m_listRegionData.contains(n))
    {
      const QList<RECT_REGION> saved = LoadCCBoxes(n);
      if (!saved.isEmpty())
        m_listRegionData[n] = saved;
    }
    if (m_listRegionData.contains(n))
      rects = m_listRegionData[n];
  }
  ui->widgetImageView->LoadImage(fn, mask_fn, pts, rects);
  ui->widgetImageView->SetPointsLocked(fn != m_listInputFiles[n].absoluteFilePath()
                                       && ui->stackedWidget->currentWidget() == ui->pageCorrection);
  // Loading shows or hides this photo's points (e.g. switching Original/Corrected)
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
    UpdateCorrClearEnabled();
  if (bNoPreviewYet)
    ui->widgetImageView->ShowMessage("No corrected preview available yet for this image.\nClick \"Correct this photo\" to generate it.");
  else if (m_setCalibratedWorking.contains(n) && ui->stackedWidget->currentWidget() == ui->pageCorrection)
    ui->widgetImageView->ShowMessage("Processing...");
  if (bNoProcessing)
    return;

  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    // Calibrated mode corrects every photo once, when Correction opens. Later
    // loads are for review.
    if (!m_strCalibrationFile.isEmpty() && !m_bCalibratedRunTried && m_proc->state() == QProcess::NotRunning)
      StartCalibratedCorrection();
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    m_bCCViewPending = !IsCCReady(n);
    if (m_bCCViewPending)
    {
      const bool bWorking = m_bCCPassActive
          || (m_procCCWarmup && m_procCCWarmup->state() != QProcess::NotRunning);
      if (!QFileInfo(CCMaskPath(n)).isFile())
        ui->widgetImageView->ShowMessage("This photo has no mask yet.\nGo back to Segmentation to make one.");
      else if (bWorking)
        ui->widgetImageView->ShowMessage("Finding connected components for this photo...");
      else
        ui->widgetImageView->ShowMessage("Connected components could not be found for this photo.\nGo back to Segmentation and check its mask.");
      return;
    }
    QString fn = CCNpzPath(CCMaskPath(n));
    if (QFile::exists(fn) && !m_maskProcessor.Load(fn))
    {
      Warn(this, "Corrected photos could not be read",
                             "The correction step finished, but its output could not be read back.\n\n"
                             "The .npz files in the output folder may be incomplete. Try running the correction again.");
      return;
    }

    m_bCCShowingSavedOnly = false;
    if (!rects.isEmpty())
    {
      QList<QPoint> pts;
      for (int i = 0; i < rects.size(); i++)
        pts << rects[i].first << rects[i].second;
      m_maskProcessor.LoadSelections(pts);
      ui->widgetImageView->SetOverlay(m_maskProcessor.GetMaskImage(m_listStockColors));
    }
    // Saved before boxes were stored: show the saved components as they are. New
    // boxes start the photo over (OnLastRegionEdited()), so slab numbers never mix.
    // Not after Clear, which leaves an empty box list.
    else if (!m_listRegionData.contains(n) && m_maskProcessor.LoadSaved(CCSavedPath(n)))
    {
      m_bCCShowingSavedOnly = true;
      ui->widgetImageView->SetOverlay(m_maskProcessor.GetMaskImage(m_listStockColors));
      ui->widgetImageView->ShowMessage("Saved components. Drawing a box starts this photo over.");
    }
    if (n == m_nIndex)
      UpdateCCUnsaved();
  }
}

QList<QPoint> MainWindow::GetCalibrationPointsList(const QVariantMap& info)
{
  QList<QPoint> pts;
  QPoint pt;
  for (int i = 0; i < 4; i++)
  {
    QString id_str = QString("corner%1").arg(i);
    pt.setX(info[id_str].toList().first().toDouble());
    pt.setY(info[id_str].toList().last().toDouble());
    pts << pt;
  }
  return pts;
}

void MainWindow::UpdatePreviewAvailability()
{
  if (ui->stackedWidget->currentWidget() != ui->pageCorrection
      || m_nIndex < 0 || m_nIndex >= m_listInputFiles.size())
    return;

  QString fn_p = GetDeformedImagePath(m_listInputFiles[m_nIndex]);
  bool bHasCorrected = IsCompleteImageOutput(fn_p);

  // Enable before selecting: SetCurrentIndex() ignores a disabled segment.
  ui->segControlPreview->SetItemEnabled(1, bHasCorrected);
  SetPreviewChecked(bHasCorrected);
}

bool MainWindow::PreviewChecked() const
{
  return ui->segControlPreview->CurrentIndex() == 1;
}

void MainWindow::SetPreviewChecked(bool bChecked)
{
  // Silent: callers load the matching image themselves, so this must not also
  // reach OnButtonTogglePreview() (a second load + crossfade).
  ui->segControlPreview->blockSignals(true);
  ui->segControlPreview->SetCurrentIndex(bChecked ? 1 : 0, false);
  ui->segControlPreview->blockSignals(false);
}

// Clear points is only meaningful with a point to clear, and not while a
// correction is running (the same moments "Correct this photo" is greyed).
void MainWindow::UpdateCorrClearEnabled()
{
  ui->pushButtonClear->setEnabled(ui->pushButtonProcess->isEnabled()
                                  && !ui->widgetImageView->GetEditedPoints().isEmpty());
}

// Reads one header per photo; called only when the answer can change. Other
// updates reuse this count.
void MainWindow::RecountCorrected()
{
  m_nCorrectedCount = 0;
  foreach (const QFileInfo& fi, m_listInputFiles)
  {
    if (IsCompleteImageOutput(GetDeformedImagePath(fi)))
      m_nCorrectedCount++;
  }
}

// The caption starts with "All" exactly when the button is enabled.
void MainWindow::UpdateCorrContinueCaption()
{
  const int nTotal = m_listInputFiles.size();
  const bool bAll = ui->pushButtonSegmentation->isEnabled() && nTotal > 0;
  // If the cached count disagrees with the button's fresh check, recount
  // (e.g. files changed outside the app).
  if (bAll != (nTotal > 0 && m_nCorrectedCount >= nTotal))
    RecountCorrected();
  QString str;
  if (bAll)
    str = (nTotal == 1) ? QString("All photos corrected") : QString("All %1 photos corrected").arg(nTotal);
  else
  {
    const QString strDone = QString(" \u00b7 %1 of %2 done").arg(qMin(m_nCorrectedCount, qMax(0, nTotal - 1))).arg(nTotal);
    if (m_proc && m_proc->state() != QProcess::NotRunning)
      str = (m_bCalibratedMode ? QString("Correcting automatically") : QString("Correcting\u2026")) + strDone;
    else
      str = QString("Correct all photos to continue") + (nTotal > 0 ? strDone : QString());
  }
  ui->corrContinueCaption->setText(str);
}

void MainWindow::StepPhoto(int nDelta)
{
  if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    if (ui->widgetImageView->SaveMaskIfEdited())
      ui->listWidgetThumbnail->SetItemImage(m_nIndex, ui->widgetImageView->GetImageInView());
    // Always, not only after an edit: continuing needs only a preview of each mask.
    MaybeWarmUpConnectedComponents();
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
    SaveCCPhoto();
  m_nIndex += nDelta;
  UpdatePreviewAvailability();
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection && PreviewChecked())
    LoadImage(m_nIndex, true, true);
  else
    LoadImage(m_nIndex);
  UpdateIndex();
}

void MainWindow::OnButtonNext()
{
  if (m_nIndex < m_listInputFiles.size()-1)
    StepPhoto(1);
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
    ui->pushButtonQuit->setEnabled(CCReadyForRecon());
  else if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    ui->pushButtonSegmentation->setEnabled(AllImagesProcessed());
    UpdateCorrContinueCaption();
  }
}

void MainWindow::OnButtonPrevious()
{
  StepPhoto(-1);
}

void MainWindow::OnButtonProcess()
{
  bool bOK;
  double dWidth = ui->lineEditRectWidth->text().trimmed().toDouble(&bOK);
  double dHeight = 1;
  int nExpectedPoints = 2;
  if (ui->segControlPointMode->CurrentIndex() == 2 || ui->segControlPointMode->CurrentIndex() == 1)
  {
    nExpectedPoints = (ui->segControlPointMode->CurrentIndex() == 2?4:3);
    dHeight = ui->lineEditRectHeight->text().trimmed().toDouble(&bOK);
    if (!bOK)
      dHeight = -1;
  }
  if (ui->widgetImageView->GetFilename() != m_listInputFiles[m_nIndex].absoluteFilePath())
    Warn(this, "Showing the corrected photo", "Switch to Original and place the points on the photo as taken.");
  else if (ui->widgetImageView->GetEditedPoints().size() < nExpectedPoints)
    Warn(this, "One more point needed", "Click the photo to place the next point.");
  else if (!bOK || dWidth <= 0 || dHeight <= 0)
    Warn(this, "Length is missing", "Enter the length of the rectangle, in millimetres.");
  else
  {
    QList<QPoint> pts = ui->widgetImageView->GetEditedPoints();
    m_listPointData[m_nIndex] = pts;
    // Remember which photo this run corrects: the user may move on before it
    // finishes.
    m_proc->setProperty("processed_index", m_nIndex);
    // Refresh now so Next unlocks as soon as correction starts.
    m_setProcessStarted.insert(m_nIndex);
    UpdateIndex();
    RunRetrospectiveCorrection(m_listInputFiles[m_nIndex].absoluteFilePath(), pts, dWidth, dHeight);
  }
}

void MainWindow::RunRetrospectiveCorrection(const QString &infile, const QList<QPoint> &pts, double w, double h)
{
  QStringList strList;
  foreach (QPoint pt, pts)
    strList << QString::number(pt.x()) << QString::number(pt.y());
  m_proc->start(m_strPythonCmd, QStringList() << m_strPyScriptRetrospective
                << "--in_img" << infile
                << "--points" << strList
                << "--width" << QString::number(w)
                << "--height" << QString::number(h)
                << "--out_dir" << m_strOutputFolder);
  m_proc->setProperty("task", "retrospective_correction");
  // Queued for segmentation when done, even if the user has left Correction.
  m_proc->setProperty("output_image", GetDeformedImagePath(QFileInfo(infile)));
}

void MainWindow::OnProcessError(QProcess::ProcessError er)
{
  QString str = "Unknown error occurred during process.";
  switch (er)
  {
  case QProcess::FailedToStart:
    str = "Failed to start python script";
    break;
  case QProcess::Crashed:
    str = "Process crashed";
    break;
  case QProcess::Timedout:
    str = "Process timed out";
    break;
  case QProcess::ReadError:
    str = "Failed to read";
    break;
  case QProcess::WriteError:
    str = "Failed to write";
    break;
  default:
    break;
  }

  ShowProgressWindow(false);
  ui->widgetImageView->HideMessage();
  ui->pushButtonProcess->setEnabled(true);
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    UpdateCorrClearEnabled();
    UpdateCorrContinueCaption();
  }
  // A pass that never started emits no finished(); a crash emits both, and
  // FinishCCPass() ignores the second call.
  if (sender() == m_proc && m_proc->property("task").toString() == "mask_to_cc")
    FinishCCPass(false);
  qDebug() << str;
  Warn(this, "That step did not finish", str);
}

void MainWindow::OnProcessStarted()
{
  ui->pushButtonProcess->setEnabled(false);
  // The CC pass runs behind the page with its own progress, and a calibrated
  // run marks the photos it is working on as it reports them.
  const QString strTask = m_proc->property("task").toString();
  if (strTask != "mask_to_cc" && strTask != "calibrated_correction")
    ui->widgetImageView->ShowMessage("Processing...");
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    UpdateCorrClearEnabled();
    UpdateCorrContinueCaption();
  }
}

void MainWindow::OnProcessFinished()
{
  ui->widgetImageView->HideMessage();
  if (m_proc->exitStatus() == QProcess::NormalExit)
    ui->pushButtonNext->setEnabled(true);
  ui->pushButtonProcess->setEnabled(true);
  if (m_nIndex >= m_listInputFiles.size()-1)
    ui->pushButtonSegmentation->setEnabled(AllImagesProcessed());
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
    UpdateCorrClearEnabled();

  ShowProgressWindow(false);

  QString task;
  if (sender())
    task = sender()->property("task").toString();
  if (task == "calibrated_correction")
  {
    OnProcessOutputMessage();   // lines still in the pipe
    m_setCalibratedWorking.clear();
    UpdateIndex();
    // A crash has already been reported by OnProcessError().
    if (m_proc->exitStatus() == QProcess::NormalExit && m_proc->exitCode() != 0)
      Warn(this, "Automatic correction could not run",
           "The photos could not be corrected with this calibration file. Check that it is the "
           "calibration saved for these photos, or correct each photo by selecting its markers.");
    else if (!m_listCalibratedFailed.isEmpty())
      Warn(this, "Markers not found on some photos",
           QString("The calibration markers could not be found on %1. On each of these, select the "
                   "center of the four markers and click Correct this photo.\n\n%2")
               .arg(m_listCalibratedFailed.size() == 1 ? QString("1 photo") : QString("%1 photos").arg(m_listCalibratedFailed.size()),
                    m_listCalibratedFailed.join("\n")));
  }
  else if (task == "mask_to_cc")
  {
    const bool bOk = (m_proc->exitStatus() == QProcess::NormalExit && m_proc->exitCode() == 0);
    FinishCCPass(bOk);
    // A crash has already been reported by OnProcessError().
    if (!bOk && m_proc->exitStatus() == QProcess::NormalExit)
      Warn(this, "Connected components could not be found",
           "The connected-components step stopped with an error. Go back to Segmentation, "
           "check the masks, and try again.");
  }
  else if (task == "retrospective_correction")
  {
    if (m_proc->exitStatus() == QProcess::NormalExit && m_proc->exitCode() == 0)
      QueueForSegmentation(m_proc->property("output_image").toString(), true);
    // Tick the photo that was processed, not the one on screen.
    int nProcessedIndex = m_proc->property("processed_index").isValid()
        ? m_proc->property("processed_index").toInt() : m_nIndex;
    ui->listWidgetThumbnail->SetItemChecked(nProcessedIndex, true);
    RecountCorrected();
    UpdateIndex();
    // Refresh the preview toggle now that the output exists (the last photo has no
    // Next to trigger it).
    UpdatePreviewAvailability();
    if (ui->stackedWidget->currentWidget() == ui->pageCorrection && PreviewChecked())
      LoadImage(m_nIndex, true, true);
  }
}

void MainWindow::OnProcessOutputMessage()
{
  if (m_proc->property("task").toString() != "calibrated_correction")
  {
    qDebug() << m_proc->readAllStandardOutput();
    return;
  }
  while (m_proc->canReadLine())
    OnCalibratedLine(QString::fromUtf8(m_proc->readLine()).trimmed());
}

void MainWindow::StartCalibratedCorrection()
{
  m_bCalibratedRunTried = true;
  // Only photos without a complete output, e.g. after an interrupted run.
  QStringList pending;
  foreach (const QFileInfo& fi, m_listInputFiles)
    if (!IsCompleteImageOutput(GetDeformedImagePath(fi)))
      pending << fi.absoluteFilePath();
  if (pending.isEmpty())
    return;
  m_setCalibratedWorking.clear();
  m_listCalibratedFailed.clear();
  m_proc->setProperty("task", "calibrated_correction");
  m_proc->start(m_strPythonCmd, QStringList() << m_strPyScriptFiducialsDetection
                << "--in_images" << pending
                << "--calibration_file" << m_strCalibrationFile
                << "--out_dir" << m_strOutputFolder);
}

// One progress line of func_fiducials_detection.py: SIZE w h, START photo,
// DONE photo "x y x y x y x y", FAILED photo reason.
void MainWindow::OnCalibratedLine(const QString& strLine)
{
  const QStringList f = strLine.split('\t');
  if (f.size() < 2)
  {
    if (!strLine.isEmpty())
      qDebug() << strLine;
    return;
  }
  if (f[0] == "SIZE" && f.size() >= 3)
  {
    // For correcting a photo by hand afterwards.
    ui->lineEditRectWidth->setText(QString::number(f[1].toDouble()));
    ui->lineEditRectHeight->setText(QString::number(f[2].toDouble()));
    return;
  }
  int n = -1;
  for (int i = 0; i < m_listInputFiles.size() && n < 0; i++)
    if (m_listInputFiles[i].absoluteFilePath() == f[1])
      n = i;
  if (n < 0)
    return;
  const bool bOnScreen = (n == m_nIndex && ui->stackedWidget->currentWidget() == ui->pageCorrection);
  if (f[0] == "START")
  {
    m_setCalibratedWorking.insert(n);
    if (bOnScreen)
      ui->widgetImageView->ShowMessage("Processing...");
    return;
  }
  m_setCalibratedWorking.remove(n);
  if (bOnScreen)
    ui->widgetImageView->HideMessage();
  if (f[0] == "DONE" && f.size() >= 3)
  {
    // This photo's own markers, shown on its Original.
    const QStringList c = f[2].split(' ', Qt::SkipEmptyParts);
    QList<QPoint> pts;
    for (int i = 0; i + 1 < c.size(); i += 2)
      pts << QPoint(qRound(c[i].toDouble()), qRound(c[i+1].toDouble()));
    m_listPointData[n] = pts;
    QueueForSegmentation(GetDeformedImagePath(m_listInputFiles[n]), true);
    ui->listWidgetThumbnail->SetItemChecked(n, true);
    m_nCorrectedCount++;
    ui->pushButtonSegmentation->setEnabled(m_nCorrectedCount >= m_listInputFiles.size() && AllImagesProcessed());
    if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
    {
      UpdateCorrContinueCaption();
      // Show the photo on screen corrected, as after "Correct this photo".
      if (n == m_nIndex)
      {
        UpdatePreviewAvailability();
        LoadImage(m_nIndex, true, PreviewChecked());
      }
    }
  }
  else if (f[0] == "FAILED")
    m_listCalibratedFailed << m_listInputFiles[n].fileName();
}

void MainWindow::OnProcessErrorMessage()
{
  QString str = m_proc->readAllStandardError().trimmed();
  if (!str.isEmpty())
    qDebug() << str;
}

void MainWindow::OnButtonClear()
{
  ui->widgetImageView->ClearEdits();
  m_maskProcessor.ClearBuffer();
  // CC page only (this slot also serves Correction's Clear points): drop the
  // saved regions too, or navigating back would restore them.
  if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    // An empty list, not none: coming back shows the cleared photo, not the boxes
    // stored with its saved components.
    m_listRegionData[m_nIndex] = QList<RECT_REGION>();
    UpdateCCUnsaved();
  }
  UpdateIndex();
}

void MainWindow::RebuildCCOverlayFromCurrentRegions()
{
  // Replay the current region list on a reset buffer: fill values follow list
  // order (MaskProcessor::LoadSelections()).
  m_maskProcessor.ClearBuffer();
  QList<RECT_REGION> regions = ui->widgetImageView->GetEditedRegions();
  QList<QPoint> pts;
  // Range-for: the comma in RECT_REGION (QPair<QPoint,QPoint>) breaks Q_FOREACH.
  for (const RECT_REGION& rc : regions)
    pts << rc.first << rc.second;
  if (!pts.isEmpty())
    m_maskProcessor.LoadSelections(pts);
  ui->widgetImageView->SetOverlay(m_maskProcessor.GetMaskImage(m_listStockColors));
  m_listRegionData[m_nIndex] = regions;
  UpdateCCUnsaved();
  UpdateIndex();
}

void MainWindow::OnButtonUndoCC()
{
  if (m_bCCViewPending)
    return;
  ui->widgetImageView->UndoLastRegion();
  RebuildCCOverlayFromCurrentRegions();
}

void MainWindow::OnButtonRedoCC()
{
  if (m_bCCViewPending)
    return;
  ui->widgetImageView->RedoLastRegion();
  RebuildCCOverlayFromCurrentRegions();
}

void MainWindow::OnButtonProceedToSeg(bool bBack)
{
  // Serves Continue from Correction and Back from Connected components.
  const bool bNavigatingBack = bBack || (sender() == ui->pushButtonBackToSegmentation);

  // Never block Back.
  if (!bNavigatingBack && !AllImagesProcessed())
  {
    Warn(this, "Some photos are not corrected",
                         "Every photo needs to be corrected before you can continue. "
                         "Correct the ones still marked not yet corrected in the strip below.");
    return;
  }
  // Refuse before the page changes and the masks prompt appears. A CC pass is
  // stopped by OnButtonLoadMask() if masks are regenerated.
  if (m_proc->state() != QProcess::NotRunning && m_proc->property("task").toString() != "mask_to_cc")
  {
    Warn(this, "Another step is still running",
         "Wait for it to finish, then try again.");
    return;
  }
  GoToStage(ui->pageSegEdit);
  SetRecenterAccent("#a6a233", "#928e2c", "#7e7b26");
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_EDIT_MASK);
  m_nIndex = 0;
  ui->widgetIndex->hide();
  ui->widgetIndexCC->hide();
  ui->widgetIndexSeg->show();
  OnButtonLoadMask();
  // Not when going back: the stage is a high-water mark for the Resume dialog.
  if (!bNavigatingBack)
    SaveProjectState("segmentation");
}

void MainWindow::OnButtonProceedToCC(bool bBack)
{
  // Entering Connected components does the same work in either direction.
  Q_UNUSED(bBack);
  // Gate Continue to Reconstruction on saved components right away; the .ui
  // starts it enabled.
  ui->pushButtonQuit->setEnabled(CCReadyForRecon());

  ui->widgetImageView->SaveMaskIfEdited();
  // No mask is waiting for the daemon any more (it is about to exit, or was
  // never told to); don't leave it running.
  if (m_listQueuedFiles.isEmpty())
    StopSegDaemon();

  GoToStage(ui->pageCC);
  SetRecenterAccent("#457c92", "#3c6d80", "#34606f");
  m_nIndex = 0;
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_REGION);
  // UniqueConnection: the page can be entered many times per session.
  connect(ui->widgetImageView, SIGNAL(LastRegionEdited(int)),
          SLOT(OnLastRegionEdited(int)), (Qt::ConnectionType)(Qt::QueuedConnection | Qt::UniqueConnection));

  m_listInputFiles = CorrectedPhotos();
  // Masks may have changed since: leaving a photo with boxes saves it again once.
  m_mapCCSavedRegions.clear();
  m_setCCUnsaved.clear();
  ui->widgetIndex->hide();
  ui->widgetIndexCC->show();
  ui->widgetIndexSeg->hide();
  OnButtonOutputFolderCC();
  SaveProjectState("connected_components");
}

void MainWindow::OnButtonProceedToRecon()
{
  GoToStage(ui->pagePhotoRecon);
  m_nIndex = 0;
}

void MainWindow::OnButtonProceedToCalibration()
{
  GoToStage(ui->pageCalibration);
  ui->pushButtonSaveCalibration->hide();
  ui->pushButtonSaveCalibration->setEnabled(false);
  ui->widgetCalibrationControl->hide();
  ui->widgetDropCalibrationPhoto->show();
  ui->widgetDownloadMarkers->show();
  // After layout, when the button has its final size.
  QTimer::singleShot(0, this, [this]() { FitMarkerCaption(); });
}

void MainWindow::FitMarkerCaption()
{
  // At most 30% wider than the printout button. Fixed, not maximum: a centred
  // word-wrapped label otherwise shrinks to its own narrow size hint. The height
  // is set too: this page's layout does not grow a label for its wrapped lines.
  QLabel* l = ui->labelDownloadMarkers;
  const int nButton = ui->toolButtonDownloadMarkers->width();
  if (nButton <= 0)
    return;
  const int w = qRound(nButton * 1.3);
  l->setFixedWidth(w);
  l->setMinimumHeight(l->fontMetrics().boundingRect(QRect(0, 0, w, 10000), Qt::TextWordWrap | Qt::AlignLeft, l->text()).height());
  // Its panel grows to match, so the printout button keeps its full height.
  QWidget* pPanel = ui->widgetDownloadMarkers;
  pPanel->layout()->invalidate();
  pPanel->setMinimumHeight(pPanel->layout()->minimumSize().height());
}

void MainWindow::OnButtonDownloadMarkers()
{
  QString strDir = QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);
  if (strDir.isEmpty())
    strDir = QDir::homePath();
  QString fn = QFileDialog::getSaveFileName(this, "Save Marker Printout",
                                            QDir(strDir).filePath("fiducial_markers_printout.png"),
                                            "PNG image (*.png)");
  if (fn.isEmpty())
    return;
  if (!fn.endsWith(".png", Qt::CaseInsensitive))
    fn += ".png";

  // QFile::copy won't overwrite, and a copy from resources is read-only: remove
  // first, then make it writable.
  if (QFile::exists(fn) && !QFile::remove(fn))
  {
    Warn(this, "Marker sheet was not saved",
                         "The file already there could not be replaced.\n\n" + fn);
    return;
  }
  if (!QFile::copy(":/resource/fiducial_markers_printout.png", fn))
  {
    Warn(this, "Marker sheet was not saved",
                         "The file could not be written.\n\n" + fn);
    return;
  }
  QFile::setPermissions(fn, QFile::ReadOwner | QFile::WriteOwner | QFile::ReadGroup | QFile::ReadOther);
}

void MainWindow::OnButtonCalibration()
{
  if (m_strCalibrationPhoto.isEmpty())
  {
    Warn(this, "No calibration photo chosen",
                         "Choose the photo of the printed marker sheet, then calibrate.");
    return;
  }
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_CALIBRATION);
  ui->widgetImageView->SetPointsLocked(false);
  ui->widgetImageView->ClearEdits();
  if (ui->widgetImageView->LoadImage(m_strCalibrationPhoto))
  {
    ui->pushButtonSaveCalibration->show();
    ui->widgetDropCalibrationPhoto->hide();
    ui->widgetDownloadMarkers->hide();
    ui->widgetCalibrationControl->show();
    ui->frameTopBar->show();
    ui->widgetPointMode->hide();
  }
  else
    Warn(this, "No photo loaded", "Choose a photo before continuing.");
}

void MainWindow::OnSliderSegOpacity(int n)
{
  // Only while the mask is shown; OnToggleMask() applies the value when it is
  // shown again.
  if (ui->pushButtonToggleMaskVisibility->isChecked())
    ui->widgetImageView->SetMaskOpacity(n/100.0);
}

void MainWindow::OnSliderLabelOpacity(int n)
{
  ui->widgetImageView->SetCombinedMaskOpacity(n/100.0);
}

void MainWindow::OnLastRegionEdited(int n)
{
  if (ui->stackedWidget->currentWidget() == ui->pageCC && (m_bCCViewPending || !IsCCReady(m_nIndex)))
  {
    ui->widgetImageView->UndoLastRegion();
    ui->widgetImageView->ShowMessage("Finding connected components for this photo...");
    return;
  }
  if (m_bCCShowingSavedOnly)
  {
    // The first new box on a photo shown from its saved components alone.
    m_maskProcessor.ClearBuffer();
    m_bCCShowingSavedOnly = false;
  }
  ui->widgetImageView->HideMessage();
  RECT_REGION rc = ui->widgetImageView->GetEditedRegions().last();
  if (m_maskProcessor.ProcessSelection(rc.first, rc.second, n+1))
  {
    ui->widgetImageView->SetOverlay(m_maskProcessor.GetMaskImage(m_listStockColors));
    // Save regions as they are drawn so navigating away keeps them.
    m_listRegionData[m_nIndex] = ui->widgetImageView->GetEditedRegions();
    UpdateCCUnsaved();
    UpdateIndex();
  }
  else
  {
    Warn(this, "That region is empty",
                         "No new slab was found inside it. Draw the region over a piece of "
                         "tissue that is not already numbered.");
  }
}

void MainWindow::OnButtonCreateMask()
{
  if (ui->widgetImageView->GetEditedRegions().isEmpty())
  {
    Warn(this, "No boxes on this photo",
         "Draw a box around each slab (Ctrl + Left drag) before saving its components.");
    return;
  }
  if (!SaveCCPhoto(true))
    return;
  // Stays on the photo; Next, the carousel or Ctrl + W move on (and save too).
  ui->widgetImageView->ShowMessage("Mask created and saved");
}

QString MainWindow::CCSavedPath(int n) const
{
  if (n < 0 || n >= m_listInputFiles.size() || m_strFinalOutputFolder.isEmpty())
    return QString();
  return QFileInfo(m_strFinalOutputFolder, m_listInputFiles[n].completeBaseName() + "_mask.npy").absoluteFilePath();
}

// Saved components made from an older mask: the mask was edited (or re-made) on
// Segmentation after they were saved, so they must be saved again.
static bool CCSavedIsStale(const QString& strSavedNpy, const QString& strMask)
{
  const QFileInfo npy(strSavedNpy), mask(strMask);
  return npy.isFile() && mask.isFile() && mask.lastModified() > npy.lastModified();
}

// One file for the project, not one per photo: connected_components/ holds only
// what the reconstruction reads. {"<photo base>": [[x1, y1, x2, y2], ...], ...}
static QString CCBoxesFile(const QString& strRoot)
{
  return strRoot + "/cc_boxes.json";
}

void MainWindow::SaveCCBoxes(int n, const QList<RECT_REGION>& list)
{
  if (n < 0 || n >= m_listInputFiles.size() || m_strRootOutputFolder.isEmpty())
    return;
  QFile f(CCBoxesFile(m_strRootOutputFolder));
  QJsonObject all;
  if (f.open(QIODevice::ReadOnly))
  {
    all = QJsonDocument::fromJson(f.readAll()).object();
    f.close();
  }
  QJsonArray boxes;
  for (const RECT_REGION& rc : list)
    boxes.append(QJsonArray({rc.first.x(), rc.first.y(), rc.second.x(), rc.second.y()}));
  all[m_listInputFiles[n].completeBaseName()] = boxes;
  QSaveFile out(CCBoxesFile(m_strRootOutputFolder));
  if (out.open(QIODevice::WriteOnly))
  {
    out.write(QJsonDocument(all).toJson(QJsonDocument::Indented));
    out.commit();
  }
}

QList<RECT_REGION> MainWindow::LoadCCBoxes(int n) const
{
  QList<RECT_REGION> list;
  // Only alongside the saved components they made: without them the boxes are stale.
  if (n < 0 || n >= m_listInputFiles.size() || !IsCompleteDataOutput(CCSavedPath(n)))
    return list;
  QFile f(CCBoxesFile(m_strRootOutputFolder));
  if (!f.open(QIODevice::ReadOnly))
    return list;
  const QJsonArray boxes = QJsonDocument::fromJson(f.readAll()).object().value(m_listInputFiles[n].completeBaseName()).toArray();
  for (const QJsonValue& v : boxes)
  {
    const QJsonArray a = v.toArray();
    if (a.size() == 4)
      list << RECT_REGION(QPoint(a[0].toInt(), a[1].toInt()), QPoint(a[2].toInt(), a[3].toInt()));
  }
  return list;
}

bool MainWindow::CCReadyForRecon() const
{
  return AllCCMasksSaved() && m_setCCUnsaved.isEmpty();
}

void MainWindow::UpdateCCUnsaved()
{
  if (ui->stackedWidget->currentWidget() != ui->pageCC || m_nIndex < 0 || m_nIndex >= m_listInputFiles.size())
    return;
  // Only photos with saved components can have unsaved changes; the others have
  // no tick and already hold back Continue to Reconstruction.
  bool bUnsaved = false;
  if (IsCompleteDataOutput(CCSavedPath(m_nIndex)))
  {
    const QList<RECT_REGION> saved = m_mapCCSavedRegions.contains(m_nIndex) ? m_mapCCSavedRegions[m_nIndex] : LoadCCBoxes(m_nIndex);
    // Edited boxes, or a mask edited on Segmentation since the save.
    bUnsaved = ui->widgetImageView->GetEditedRegions() != saved
        || CCSavedIsStale(CCSavedPath(m_nIndex), CCMaskPath(m_nIndex));
  }
  if (bUnsaved == m_setCCUnsaved.contains(m_nIndex))
    return;
  if (bUnsaved)
    m_setCCUnsaved.insert(m_nIndex);
  else
    m_setCCUnsaved.remove(m_nIndex);
  // No tick while the edits are unsaved; back on save (or after undoing them).
  ui->listWidgetThumbnail->SetItemChecked(m_nIndex, !bUnsaved);
  ui->pushButtonQuit->setEnabled(CCReadyForRecon());
  UpdateStepper(ui->pageCC);
}

bool MainWindow::SaveCCPhoto(bool bForce)
{
  if (ui->stackedWidget->currentWidget() != ui->pageCC || m_nIndex < 0 || m_nIndex >= m_listInputFiles.size()
      || m_bCCViewPending || !IsCCReady(m_nIndex))
    return false;
  // Never a photo without boxes: its components would be empty.
  const QList<RECT_REGION> list = ui->widgetImageView->GetEditedRegions();
  if (list.isEmpty())
    return false;
  if (!bForce && m_mapCCSavedRegions.contains(m_nIndex) && m_mapCCSavedRegions[m_nIndex] == list)
    return false;
  m_listRegionData[m_nIndex] = list;

  // completeBaseName(), not an extension replace, which would mangle "sample.png.png".
  QString fn = QFileInfo(ui->widgetImageView->GetFilename()).completeBaseName() + "_mask.npy";
  m_maskProcessor.SaveToNpy(QFileInfo(m_strFinalOutputFolder, fn).absoluteFilePath());
  SaveCCBoxes(m_nIndex, list);
  m_setCCProcessed.insert(m_nIndex);
  m_mapCCSavedRegions[m_nIndex] = list;
  m_setCCUnsaved.remove(m_nIndex);
  ui->pushButtonNextCC->setEnabled(true);
  ui->listWidgetThumbnail->SetItemImage(m_nIndex, ui->widgetImageView->GetImageInView());
  // The same green outline and tick as a corrected photo on Correction.
  ui->listWidgetThumbnail->SetItemChecked(m_nIndex, true);
  ui->pushButtonQuit->setEnabled(CCReadyForRecon());
  UpdateStepper(ui->pageCC);
  return true;
}

void MainWindow::OnToggleMask()
{
  // Decided from the button, not the opacity (a shown mask at 0% has opacity 0).
  // A click has already flipped the button; Ctrl+T has not.
  QPushButton* pShow = ui->pushButtonToggleMaskVisibility;
  const bool bShow = (sender() == pShow) ? pShow->isChecked() : !pShow->isChecked();
  if (bShow)
    ui->widgetImageView->SetMaskOpacity(ui->sliderSegOpacity->Value()/100.0);
  else
    ui->widgetImageView->SetMaskOpacity(0);

  // Sync the button for every trigger, with signals blocked.
  ui->pushButtonToggleMaskVisibility->blockSignals(true);
  ui->pushButtonToggleMaskVisibility->setChecked(bShow);
  ui->pushButtonToggleMaskVisibility->blockSignals(false);
}

void MainWindow::OnFileChanged(const QString& path)
{
  // The queue holds full paths. Results may already exist or arrive in any
  // order; each .npz appears complete (the daemon renames it into place).
  Q_UNUSED(path);
  if (SegResultReady())
  {
    const QFileInfo fiNpz(m_listQueuedFiles.takeFirst());
    const QString strMaskBaseName = fiNpz.completeBaseName();
    int nPhoto = -1;
    for (int i = 0; i < m_listInputFiles.size() && nPhoto < 0; i++)
    {
      if (m_listInputFiles[i].completeBaseName() == strMaskBaseName)
        nPhoto = i;
    }
    if (fiNpz.exists())
    {
      const QString fn = fiNpz.absoluteFilePath();
      QImage image = NpyToImage(fn, m_mapInputImageSize[fn].toSize());
      image.save(QFileInfo(m_strMaskFolder, strMaskBaseName+".png").absoluteFilePath());
      // This mask may complete Segmentation while the user is on a later page, where
      // the tick is shown.
      if (ui->stackedWidget->currentWidget() == ui->pageCC
          || ui->stackedWidget->currentWidget() == ui->pagePhotoRecon)
        UpdateStepper(ui->stackedWidget->currentWidget());

      // A new mask has not been previewed yet.
      UpdateThumbnailWithMask(nPhoto);
      m_setSegPreviewed.remove(nPhoto);
    }
    else
    {
      // The daemon gave up on this photo: it is drawn by hand, on the empty
      // mask LoadImage() gives it.
      m_listSegFailed << strMaskBaseName;
    }
    // On screen without its mask: show it now (or the empty one to draw on), unless
    // something was drawn meanwhile. UpdateIndex() below counts it as previewed.
    if (nPhoto >= 0 && nPhoto == m_nIndex && m_bSegFirstMaskShown
        && ui->stackedWidget->currentWidget() == ui->pageSegEdit
        && !ui->widgetImageView->IsMaskEdited())
      LoadImage(nPhoto);
    // A result exists, so the model has loaded: hide the loading overlay
    // (idempotent).
    ui->pageSegEdit->setEnabled(true);
    ShowProgressWindow(false);
    ui->widgetSegCtrls->setEnabled(true);

    int nDone = m_nTotalSegFiles - m_listQueuedFiles.size();
    if (m_listQueuedFiles.isEmpty())
    {
      // All masks exist: start the CC warm-up now so the CC page opens without a wait.
      MaybeWarmUpConnectedComponents();
      ui->widgetSegProgress->setVisible(false);
    }
    else
    {
      ui->progressBarSeg->setValue(nDone);
      ui->labelSegProgress->setText(QString("Segmenting images in the background: %1 of %2 ready").arg(nDone).arg(m_nTotalSegFiles));
    }
    UpdateIndex();
  }

  // Show the first photo with a mask as soon as there is one, once per visit.
  if (!m_bSegFirstMaskShown && ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    for (int i = 0; i < m_listInputFiles.size(); i++)
    {
      if (!QFileInfo(CCMaskPath(i)).isFile())
        continue;
      m_bSegFirstMaskShown = true;
      ui->pageSegEdit->setEnabled(true);
      ShowProgressWindow(false);
      ui->widgetSegCtrls->setEnabled(true);
      m_nIndex = i;
      LoadImage(i);
      UpdateIndex();
      if (!SegResultReady() && !m_listQueuedFiles.isEmpty() && m_procSeg->state() != QProcess::NotRunning)
        Warn(this, "Masks are still being generated",
                                 "You can start editing the one that is ready while the rest finish.");
      break;
    }
  }

  // One result per event-loop turn, so a backlog does not freeze the page. Not
  // a 0 ms timer: a chain of those starves other timers on macOS.
  if (SegResultReady())
    QTimer::singleShot(1, this, [this]() { OnFileChanged(m_strSegOutDir); });
  else if (m_bSegDaemonExited)
    ReportSegDaemonExit();
}

QImage MainWindow::NpyToImage(const QString& npy_in, const QSize& sz_in)
{
  QImage image;
  if (!QFile::exists(npy_in))
  {
    qDebug() << "File does not exist:" << npy_in;
    return image;
  }
  cnpy::NpyArray ar = cnpy::npz_load(qPrintable(npy_in), "probabilities");
  if (ar.shape.size() != 4)
  {
    qDebug() << "Could not load numpy file " << npy_in;
    return image;
  }

  float* ptr = ar.data<float>();
  int w = ar.shape[3];
  int h = ar.shape[2];
  image = QImage(QSize(w, h), QImage::Format_ARGB32);
  image.fill(QColor(0,0,0,255));

  for (int y = 0; y < h; y++)
  {
    QRgb* p = (QRgb*)image.scanLine(y);
    for (int x = 0; x < w; x++)
    {
      if (ptr[w*h+y*w+x] >= 0.5)
      {
        p[x] = qRgba(255,255,255,255);
      }
    }
  }

  QSize sz = sz_in;
  if (!sz.isValid())
    sz = QImage(m_listInputFiles.first().absoluteFilePath()).size();
  image = image.scaled(sz, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
  w = image.width();
  h = image.height();
  for (int y = 0; y < h; y++)
  {
    QRgb* p = (QRgb*)image.scanLine(y);
    for (int x = 0; x < w; x++)
    {
      if (qRed(p[x]) >= 128)
        p[x] = qRgba(255,255,255,255);
      else
        p[x] = qRgba(0,0,0,255);
    }
  }

  return image;
}

void MainWindow::ClearFolder(const QString& path)
{
  QFileInfoList flist = QDir(path).entryInfoList(QDir::Files, QDir::Name);
  foreach (QFileInfo fi, flist)
    QFile::remove(fi.absoluteFilePath());
}

void MainWindow::OnThumbnailIndexChanged(int n)
{
  if (n >= 0)
  {
    if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
    {
      if (ui->widgetImageView->SaveMaskIfEdited())
      {
        QImage img = ui->widgetImageView->GetImageInView();
        if (!img.isNull() && m_nIndex < ui->listWidgetThumbnail->count())
          ui->listWidgetThumbnail->SetItemImage(m_nIndex, img);
      }
      // Always, as in StepPhoto().
      MaybeWarmUpConnectedComponents();
    }
    // Leaving a photo saves its components, as "Save components" would.
    else if (ui->stackedWidget->currentWidget() == ui->pageCC && n != m_nIndex)
      SaveCCPhoto();
    m_nIndex = n;
    UpdatePreviewAvailability();
    if (ui->stackedWidget->currentWidget() == ui->pageCorrection && PreviewChecked())
      LoadImage(n, true, true);
    else
      LoadImage(n);
    UpdateIndex();
  }
}

void MainWindow::OnButtonTogglePreview(bool b)
{
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection && m_nIndex >= 0)
  {
    // Crossfade the toggle only. Armed here, not in LoadImage(), because bPreview is
    // also true when navigating with the toggle on. Keep it inside this guard,
    // right before its load: the arm is consumed by the next LoadImage().
    ui->widgetImageView->BeginPreviewCrossfade();
    LoadImage(m_nIndex, true, b);
  }
}
