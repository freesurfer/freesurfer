#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileInfoList>
#include <QProcess>
#include <QVariantMap>
#include <QElapsedTimer>
#include "WidgetImageView.h"
#include "MaskProcessor.h"
#include "ProjectManager.h"
#include <QFileSystemWatcher>
#include <QTimer>
#include <QHash>
#include <QSet>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class QLabel;
class QToolButton;
class ProgressWindow;
class ImageGridView;
class StepperBar;
class ReconPage;

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

  void resizeEvent(QResizeEvent* e);
  void moveEvent(QMoveEvent* e);
  void showEvent(QShowEvent* e);
  // Quitting with a reconstruction running asks first; Switch Project does not.
  void closeEvent(QCloseEvent* e);

  bool eventFilter(QObject *watched, QEvent *event);

  // Opens a project from the Projects Hub, asks where to resume, and replays the
  // normal stage-entry handlers up to that stage, so validation and prompts match
  // a fresh visit. No-op for an empty folder.
  void LoadProject(const QString& projectFolder);

signals:
  // Emitted by File > Switch Project; main() shows the Hub again.
  void RequestSwitchProject();

public slots:
  void OnButtonInputFolder();
  void OnButtonOutputFolderCC();
  // Runs the connected-components pass, or waits for a running warm-up or step to
  // finish first.
  void RunOrWaitForMaskToCC();
  // bBack: the caller is navigating backward. Not inferable from sender(): the
  // step bar calls these directly.
  void OnButtonContinue(bool bBack = false);
  void OnButtonCalibrationFile();
  void OnPointModeChanged(int index);
  void OnButtonPrevious();
  void OnButtonNext();
  void OnButtonProcess();
  void OnButtonClear();
  // Undo/redo a CC box and rebuild the overlay to match.
  void OnButtonUndoCC();
  void OnButtonRedoCC();
  void OnButtonLoadMask();
  void OnButtonCreateMask();
  void OnLinkActivated(const QString& link);
  void OnButtonChangeOutputFolder();
  void OnButtonTogglePreview(bool b);
  void OnButtonCalibrationPhoto();
  void OnButtonBackToStart();
  void OnButtonSaveCalibration();

private slots:
  void OnProcessOutputMessage();
  void OnProcessErrorMessage();
  void OnProcessStarted();
  void OnProcessFinished();
  void OnProcessError(QProcess::ProcessError);
  void OnButtonProceedToSeg(bool bBack = false);
  void OnButtonProceedToCC(bool bBack = false);
  void OnButtonProceedToRecon();
  void OnButtonProceedToCalibration();
  void OnButtonDownloadMarkers();
  void OnSliderSegOpacity(int n);
  void OnSliderLabelOpacity(int n);
  void OnLastRegionEdited(int n);
  void OnToggleMask();
  void OnFileChanged(const QString& path);
  void ShowProgressWindow(bool bShow = true);
  void OnThumbnailIndexChanged(int n);
  void OnButtonCalibration();
  void OnPointEdited();
  void OnActionSwitchProject();

  // Shows the Start page grid instead of widgetImageView while pageStartUp is
  // current and has photos. Runs on every page change and after a folder is
  // picked.
  void UpdateStartupPreview();
  void OnToggleDarkMode(bool bDark);

private:
  void SetupScriptPath();

  // The segmentation daemon (process_directory_no_upsample_daemon.py) segments
  // corrected photos in the background from Correction on. StartSegDaemon()
  // returns false when nnUNet is not configured.
  bool StartSegDaemon();
  // On Correction: keeps the daemon running, or starts one with the corrected
  // photos that have no mask yet.
  void SegmentCorrectedPhotos();
  // Hands a corrected image to the daemon (only if one runs and the image is
  // complete) and drops its old result. bRecorrected also deletes its mask and
  // connected components, so they are made again.
  // False when the photo could not be handed to the daemon (it is not running,
  // or the copy failed): no result will ever come for it.
  bool QueueForSegmentation(const QString& strImage, bool bRecorrected);
  // The saved mask for a corrected photo's completeBaseName().
  QString MaskPathFor(const QString& strBase) const;
  // Tells the daemon to exit once its queue is done.
  void FinishSegQueue();
  // Kills the daemon, if running, and forgets its folders and queue.
  void StopSegDaemon();
  // Off Correction, the results still waiting are taken one per event-loop
  // turn first; then ReportSegDaemonExit() warns about missing masks.
  void OnSegDaemonFinished(QProcess::ExitStatus status, int nCode);
  void ReportSegDaemonExit();
  // Stops waiting for masks: clears m_listQueuedFiles and the progress UI.
  void ForgetSegProgress();
  // The head of m_listQueuedFiles has a result (.npz or .failed).
  bool SegResultReady() const;

  // Rebuilds one carousel thumbnail with its newly saved mask.
  void UpdateThumbnailWithMask(int nIndex);
  // Photo n's mask is still expected from the daemon.
  bool SegMaskPending(int n) const;
  // Enables "Continue to Connected Components" once every photo has a mask that
  // has been on screen, and says in the note at the top what is still missing.
  void UpdateSegReadiness();
  void SetUpNavigationShortcuts();
  // Sizes the caption under "Download Marker Printout" to the button.
  void FitMarkerCaption();
  // Connected components: saves the photo on screen when it has boxes that
  // were not saved yet (bForce: even if unchanged). False when nothing was saved.
  bool SaveCCPhoto(bool bForce = false);
  // Every photo saved, and no saved photo edited since: Continue to Reconstruction.
  bool CCReadyForRecon() const;
  // Marks the photo on screen unsaved (no tick, no Continue) when its boxes
  // differ from the ones it was saved with.
  void UpdateCCUnsaved();
  // The saved components of photo n (<base>_mask.npy), and the boxes they were
  // made from, kept in <project>/cc_boxes.json.
  QString CCSavedPath(int n) const;
  void SaveCCBoxes(int n, const QList<RECT_REGION>& list);
  QList<RECT_REGION> LoadCCBoxes(int n) const;

  // On Segmentation navigation (and when all masks arrive), runs the CC pass in
  // the background on masks saved since their last pass, so the CC page usually
  // opens with nothing to do.
  void MaybeWarmUpConnectedComponents();
  void UpdateIndex();
  // Next/Previous: saves the edits on screen (mask, components), moves nDelta photos, loads it.
  void StepPhoto(int nDelta);

  // Photo formats used when scanning folders, so other files (e.g. a calibration
  // .npz) are never treated as photos.
  static QStringList ImageNameFilters();

  // Every input photo has a complete corrected output (gates Continue).
  bool AllImagesProcessed() const;
  // Same test over an explicit list of original photos; the no-arg form uses
  // m_listInputFiles, which holds originals only on Start and Correction.
  bool AllImagesProcessed(const QFileInfoList& originals) const;
  // Every photo has saved components in m_strFinalOutputFolder.
  bool AllCCMasksSaved() const;
  void RebuildCCOverlayFromCurrentRegions();
  // The corrected output path for an input photo, matching
  // func_retrospective_correction.py's naming.
  QString GetDeformedImagePath(const QFileInfo& originalFile) const;
  // The corrected photos, in the same order as the originals on the Correction page.
  QFileInfoList CorrectedPhotos() const;
  void LoadImage(int n, bool no_processing = false, bool bPreview = false);
  void ClearFolder(const QString& path);
  QImage NpyToImage(const QString& fn, const QSize& sz);
  QList<QPoint> GetCalibrationPointsList(const QVariantMap& info);
  void RepositionProgressWindow();

  // Enables "Corrected" when the current photo has a complete corrected output,
  // and selects it if so. Call after m_nIndex changes, before LoadImage().
  void UpdatePreviewAvailability();
  void SetInputFolder(const QString& path, const QString& root_output = "");
  void SetCalibrationFile(const QString& path);
  void SetCalibrationPhoto(const QString& path);
  void RunRetrospectiveCorrection(const QString& infile, const QList<QPoint>& pts, double w, double h);

  // Saves paths, mode and the stage (as a high-water mark) to project.json. No-op
  // without a project.
  void SaveProjectState(const QString& stage);

  // Applies the light or dark palette and the panel stylesheet (on top of
  // AppTheme::BaseStyleSheet()).
  void SetDarkMode(bool bDark);

  // The single place the stage changes: switches the page, sets the canvas bars
  // and updates the step bar. Per-page setup (edit mode, accent, index widgets)
  // stays with each caller.
  void GoToStage(QWidget* pPage);
  // Builds the pipeline strip across the top of the window and puts the
  // existing splitter underneath it.
  void SetUpStepper();
  // Builds the body of the Reconstruction page, which the .ui leaves as a
  // title and two spacers, and fills in the folders the pipeline produced.
  void SetUpReconPage();
  void SetUpStartCaption();
  // Puts the Segmentation panel in the same order as Correction and
  // Connected components.
  void NormaliseSegPanel();
  void UpdateStartCaption(int nPhotos);
  // Connected components background pass; see m_bCCPassActive.
  bool IsCCReady(int n) const;
  void PollCCPass();
  void UpdateCCProgress();
  void FinishCCPass(bool bSucceeded = true);
  // Per-mask freshness of the scratch .npz: it exists and came from a pass
  // that started after the mask was last saved.
  bool IsCCMaskFresh(const QString& strMask) const;
  QStringList StaleCCMasks() const;
  // The mask for photo n, paired by name; and the scratch result for a mask.
  QString CCMaskPath(int n) const;
  QString CCNpzPath(const QString& strMask) const;
  struct CCMaskSig { qint64 mtime = -1; qint64 size = -1; QByteArray md5; };
  static CCMaskSig ReadCCMaskSig(const QString& strMask, bool bWithDigest);
  static bool SameCCMaskSig(const CCMaskSig& recorded, const QString& strMask);
  // Copies these masks into a scratch folder in m_sTempDir, removing their old
  // results and recording each mask's identity first. Returns the folder.
  QString StageCCMasks(const QStringList& listMasks, const QString& strDirName,
                       QHash<QString, CCMaskSig>& hashSig);
  // Copies the photos ticked in the grid from m_strPendingSourceFolder into
  // <project>/photos and points the input folder there. Returns false (after a
  // warning) on failure.
  bool CopyPendingPhotosIntoProject();
  void UpdateReconPage();
  // Recomputes which steps are done / current / still ahead for pPage.
  void UpdateStepper(QWidget* pPage);

  // Colours the shared top-bar Re-Center button with the current step's accent.
  void SetRecenterAccent(const QString& strBase, const QString& strHover,
                         const QString& strPressed);

  // Re-styles the ghost pill buttons (Undo/Redo/Clear, Back, Change) for the
  // theme. Their .ui stylesheets beat the app sheet, so they must be restyled per
  // widget; geometry is kept, only colours change.
  void RethemeGhostPillButtons(bool bDark);

  // Sets the icon glyphs on buttons. Called once: the icons don't depend on the
  // theme.
  void ApplyButtonIcons();

  // Gives the primary action buttons of the editing pages one shared minimum
  // width, measured because substituted fonts render wider on Linux.
  void EqualizePrimaryActionButtonWidths();

  bool m_bDarkMode = false;

  Ui::MainWindow *ui;

  // The project's root folder as opened from the Hub. Stays fixed even if the
  // output folder is changed, so project.json stays where it was registered.
  // Empty without a project.
  QString m_strCurrentProjectFolder;
  ProjectInfo m_currentProjectInfo;
  // True while LoadProject() replays navigation; suppresses one-time notices.
  bool m_bResumingProject = false;
  QString  m_strInputFolder;
  QString  m_strRootOutputFolder;
  QString  m_strOutputFolder;
  QString  m_strCalibrationFile;
  QString  m_strCalibrationPhoto;
  QString  m_strFinalOutputFolder;
  QString  m_strMaskFolder;

  QString  m_strNNUnetScriptFolder;
  QString  m_strNNUnetModelFolder;

  QString  m_strPythonCmd;

  QFileInfoList  m_listInputFiles;
  QFileInfoList  m_listMaskFiles;
  // Each photo's last region list, restored by LoadImage() on the CC page. Kept
  // in sync with every edit, undo, redo and clear.
  QMap<int, QList<RECT_REGION> > m_listRegionData;
  QList<QColor> m_listStockColors;

  int m_nIndex;
  QMap<int, QList<QPoint> > m_listPointData;
  bool m_bCalibratedMode;
  QString m_strPyScriptRetrospective;
  QString m_strPyScriptFiducialsDetection;
  QString m_strPyScriptFiducialsCalibration;
  QString m_strPyScriptMaskToCC;

  QProcess* m_proc;

  // Also separate from m_proc -- see MaybeWarmUpConnectedComponents().
  QProcess* m_procCCWarmup = nullptr;
  // The CC pass (func_mask_to_cc.py on m_proc) runs in the background. It writes
  // one .npz per mask in sorted path order; a photo is ready once a later one
  // exists or the pass has finished. m_bCCPassActive covers the whole request,
  // including waiting for another process.
  bool m_bCCPassActive = false;
  bool m_bCCPassStarted = false;
  bool m_bCCViewPending = false;       // the photo on screen is waiting for its .npz
  QStringList m_listCCOrder;           // mask paths, in the script's order
  int m_nCCNewest = -1;                // highest position in m_listCCOrder written so far
  QTimer* m_timerCCPoll = nullptr;
  // Mask path -> the identity the mask had when the pass that wrote its .npz read
  // it. Recorded when the pass finishes cleanly; fresh while it still matches.
  QHash<QString, CCMaskSig> m_hashCCNpzFrom;
  QHash<QString, CCMaskSig> m_hashCCPassSig;   // staged by the page pass
  QHash<QString, CCMaskSig> m_hashCCWarmSig;   // staged by the warm-up

  QString m_sTempDir;

  MaskProcessor  m_maskProcessor;
  QElapsedTimer  m_elapsedTimer;
  QFileSystemWatcher  m_fileWatcher;
  QStringList    m_listQueuedFiles;
  int            m_nTotalSegFiles = 0;
  QVariantMap    m_mapInputImageSize;

  QProcess*      m_procSeg = nullptr;
  // Its input and output folders, fresh per daemon.
  QString        m_strSegInDir;
  QString        m_strSegOutDir;
  // A mask has been shown (and the wait noted) since Segmentation was entered.
  bool           m_bSegFirstMaskShown = false;
  // Photos the daemon gave up on this Segmentation visit.
  QStringList    m_listSegFailed;
  // The daemon exited (and how) while Segmentation still takes its results.
  bool           m_bSegDaemonExited = false;
  bool           m_bSegDaemonOk = false;
  // Its stderr (the tail), shown if it fails. Cleared per daemon.
  QString        m_strNNUnetStderr;

  // Calibrated mode: one run detects each photo's markers and corrects it.
  // Tried once per visit to Correction; photos it is working on, and photos
  // whose markers it could not find (left to correct by hand).
  bool           m_bCalibratedRunTried = false;
  QSet<int>      m_setCalibratedWorking;
  QStringList    m_listCalibratedFailed;
  void StartCalibratedCorrection();
  void OnCalibratedLine(const QString& strLine);

  // Photos for which "Correct this photo" was started (retrospective mode).
  // Gates Next on Correction.
  QSet<int>      m_setProcessStarted;

  // Photos whose components were saved this visit. Gates Next on Connected
  // components.
  QSet<int>      m_setCCProcessed;
  // The boxes each photo's components were last saved with this visit.
  QMap<int, QList<RECT_REGION> > m_mapCCSavedRegions;
  // The photo on screen shows saved components that have no stored boxes.
  bool           m_bCCShowingSavedOnly = false;
  // Photos with saved components whose boxes were changed and not saved yet.
  QSet<int>      m_setCCUnsaved;

  // Segmentation photos viewed at least once this visit. Continue to connected
  // components requires all of them. Cleared on entry.
  QSet<int>      m_setSegPreviewed;

  ProgressWindow* m_wndProgress;

  // Start page grid preview in widgetMain (see UpdateStartupPreview()). Parented
  // to widgetMain.
  ImageGridView* m_pImageGridView = nullptr;
  StepperBar*    m_pStepper = nullptr;
  // Readiness caption above Get started.
  QLabel*        m_pStartCaption = nullptr;
  // The calibration pictogram captions and marker count.
  QLabel*        m_pCalibCap1 = nullptr;
  QLabel*        m_pCalibCap2 = nullptr;
  QLabel*        m_pCalibProgress = nullptr;
  int            m_nCalibrationBoxFloor = 0; // widgetDropCalibration's .ui minimum height
  void SetUpCalibrationInstructions();
  // Puts everything above pFooter into a scroll area, leaving the footer
  // pinned. See the implementation for why this cannot use addLayout().
  void ScrollPageContent(QWidget* pPage, QWidget* pFooter);
private slots:
  void UpdateCalibrationProgress();
  void OnToggleCalibrationDisclosure(bool bOpen);
private:
  QLabel*        m_pSegContinueCaption = nullptr;
  bool           m_bSwitchingProject = false;
  // Home button at the left of the step bar.
  QToolButton*   m_pHomeButton = nullptr;
  // Under the Segmentation description: every mask must be previewed first.
  QLabel*        m_pSegPreviewNote = nullptr;
  // Segmentation tool row (brush / eraser / fill).
  QToolButton*   m_pToolBrush = nullptr;
  QToolButton*   m_pToolEraser = nullptr;
  QToolButton*   m_pToolFill = nullptr;
  void SetMaskToolFromRow();
  void UpdateSegContinueCaption();
  // Correction page: the Show pill's state, the Clear points gate, and the
  // readiness caption above "Continue to Segmentation" (corrContinueCaption).
  bool PreviewChecked() const;
  void SetPreviewChecked(bool bChecked);
  void UpdateCorrClearEnabled();
  void RecountCorrected();
  void UpdateCorrContinueCaption();
  int            m_nCorrectedCount = 0;  // cached by RecountCorrected()
  ReconPage*     m_pReconPage = nullptr;
  // Files last passed to SetImages(), so redundant calls don't reload or reset the
  // lightbox.
  QStringList m_listStartupPreviewFiles;
  // A folder picked on the Start page whose photos are not copied yet. The ticked
  // photos are copied on Get started. Empty otherwise.
  QString m_strPendingSourceFolder;
};
#endif // MAINWINDOW_H
