#ifndef RECONPAGE_H
#define RECONPAGE_H

// The pipeline's last step: run FreeSurfer's mri_3d_photo_recon on the earlier
// steps' output. A guided form on the left; the command and its log on the
// right, where other steps show a photo.
//
// Two roots, one class
// --------------------
//  * this            the form in the left column (pagePhotoRecon): a scrolling
//                    form above a footer with the readiness caption and the
//                    primary button.
//  * CommandPane()   the Command card over the Log card. MainWindow inserts it
//                    into widgetMain, shown only on pagePhotoRecon; it is then
//                    no longer our child, so SetDarkMode() styles both roots.
//  One state machine drives both halves, and one BuildArguments() feeds both
//  the preview and the run.
//
// Contract
// --------
//  * BuildArguments() emits the six arguments utils.py get_arguments requires
//    (photo resolution fixed at 0.1 mm; an empty slab thickness blocks Run)
//    and --threads; other empty fields add no flag. No GPU field on macOS.
//  * SetProject() pre-fills paths from the project root and restores its saved
//    form. It never runs anything or touches the form during a run.
//  * Availability is checked up front (FREESURFER_HOME, mri_3d_photo_recon,
//    fspython, imputation model). DISSECTION_PHOTO_RECON_CMD overrides the
//    program and skips these checks (used for testing).
//  * A run belongs to its project, not to the page: ReconRunner keeps it (and
//    its log) for the life of the app, so Switch Project leaves it running and
//    the project's page shows it again when reopened. Quitting stops it
//    (ReconRunner::ConfirmQuit() asks first).

#include <QWidget>
#include <QProcess>
#include <QStringList>
#include <QElapsedTimer>
#include <QMap>
#include <QList>
#include <QString>
#include <QVariantMap>
#include <QPointer>

class QLabel;
class QLineEdit;
class QCheckBox;
class QRadioButton;
class QPushButton;
class QPlainTextEdit;
class QProgressBar;
class QScrollArea;
class QVBoxLayout;
class QTimer;
class SegmentedControl;
class ReconFormRow;

// One reconstruction: the process, its log and its outcome. Kept per project
// root until the app quits, after it finishes too, so the result is still shown
// when the project is reopened.
class ReconRunner : public QObject
{
  Q_OBJECT
public:
  enum Outcome { Active, Succeeded, Failed, Stopped };

  // The run for this project (running or finished this session), or null.
  static ReconRunner* ForRoot(const QString& strRoot);
  // Replaces any finished run of the project.
  static ReconRunner* Start(const QString& strRoot, const QString& strProgram, const QStringList& args,
                            const QString& strCommandLine);
  static QStringList RunningRoots();
  // With runs going: asks whether to stop them. False: the user chose to keep them.
  static bool ConfirmQuit(QWidget* parent);
  // Terminates every run, then kills what is left (at quit).
  static void StopAll();

  bool IsActive() const { return m_outcome == Active; }
  bool IsStopping() const { return m_bStopRequested && m_outcome == Active; }
  Outcome GetOutcome() const { return m_outcome; }
  QStringList Lines() const { return m_listLines; }
  QString Phase() const { return m_strPhase; }
  qint64 ElapsedMs() const { return m_outcome == Active ? m_timer.elapsed() : m_nFinalMs; }
  int ExitCode() const { return m_nExitCode; }
  bool Crashed() const { return m_bCrashed; }
  bool FailedToStart() const { return m_bFailedToStart; }
  // SIGTERM, then SIGKILL after 5 s.
  void Stop();

signals:
  void LineAdded(const QString& strLine);
  void PhaseChanged(const QString& strPhase);
  void Finished();

private:
  ReconRunner(const QString& strRoot, QObject* parent);
  void AppendChunk(const QString& strChunk, QString& strBuffer);
  void NoteLine(const QString& strLine);
  void OnProcFinished(int nExitCode, QProcess::ExitStatus status);
  void OnProcError(QProcess::ProcessError err);
  void Finish(Outcome o);

  QString m_strRoot;
  QString m_strProgram;
  QProcess* m_proc = nullptr;
  QTimer* m_pKill = nullptr;
  QElapsedTimer m_timer;
  qint64 m_nFinalMs = 0;
  Outcome m_outcome = Active;
  bool m_bStopRequested = false;
  bool m_bCrashed = false;
  bool m_bFailedToStart = false;
  int m_nExitCode = 0;
  QString m_strBufOut;
  QString m_strBufErr;
  QString m_strPhase;
  QStringList m_listLines;
};

class ReconPage : public QWidget
{
  Q_OBJECT
public:
  enum State { Idle, Running, Stopping, Succeeded, Failed, Stopped, Unavailable };

  struct Issue
  {
    QString key;
    QString message;
    bool blocking;
  };

  explicit ReconPage(QWidget* parent = nullptr);
  ~ReconPage();

  QWidget* CommandPane() const { return m_pPane; }

  void SetProject(const QString& strRoot);
  bool IsRunning() const;

  // Pure: form -> argv, without the program.
  QStringList BuildArguments() const;
  // multiline: one flag (and its value) per line, as the Command card shows
  // it; otherwise the single shell line Copy command copies.
  QString CommandPreview(bool bMultiline = false) const;

  // Pure: the form's problems, blocking ones first per field.
  QList<Issue> Validate() const;

  // The resolved program path, or empty when FreeSurfer is not set up.
  QString Program() const;

  QVariantMap FormValues() const;
  // Where a project's form is saved (QSettings group), from its folder. Renaming
  // the folder moves the group (ProjectManager::RenameProject()).
  static QString SettingsKeyFor(const QString& strRoot);
  void SetFormValues(const QVariantMap& v);
  void SetAdvancedOpen(bool bOpen);
  void SetExpertOpen(bool bOpen);

  // Stops a running reconstruction without asking (the Stop button asks).
  void StopWithoutConfirm();

  // Re-applies the page's own colours to both roots.
  void SetDarkMode(bool bDark);

signals:
  void BackRequested();

protected:
  bool eventFilter(QObject* o, QEvent* e) override;
  void showEvent(QShowEvent* e) override;

private:
  // construction
  QWidget* BuildForm();
  QWidget* BuildFooter();
  QWidget* BuildPane();
  void AddSection(QVBoxLayout* into, const QString& strText, bool bFirst = false);
  void AddRow(QVBoxLayout* into, const QString& strKey, const QString& strLabel,
              QWidget* pControl, const QString& strHelp);
  QLineEdit* MakeNumberEdit(const QString& strKey, const QString& strPlaceholder, bool bInteger);
  QLineEdit* MakeExpertEdit(const QString& strKey);
  QWidget* MakeUnitRow(QLineEdit* pEdit, const QString& strUnit, int nWidth);
  QWidget* MakePathRow(const QString& strKey, const QString& strPlaceholder, bool bDirectory);
  QCheckBox* MakeCheck(const QString& strKey, const QString& strText);
  SegmentedControl* MakePill(const QStringList& items, int nDefault);

  // form logic
  void OnFormEdited();
  void UpdateVisibility();
  void UpdatePreview();
  void OnValidate();
  void ApplyIssues(const QList<Issue>& issues);
  void RefreshInventory();
  void SetNarrow(bool bNarrow);
  QString Hemisphere() const;
  QString ReferenceMode() const;
  QString Text(const char* key) const;
  // The run's folders: the project's own inputs, and the named output folder
  // inside it.
  QString PhotoDir() const;
  QString MaskDir() const;
  QString OutputDir() const;
  QStringList ListPhotos(const QString& strDir, QString* pExt = nullptr, int* pOther = nullptr) const;
  QStringList ListMasks(const QString& strDir) const;
  void ResetForm();
  void FillDefaultReferencePaths();

  // settings
  QString SettingsKey() const;
  void SaveSettings();
  void LoadSettings();

  // availability and state
  void CheckAvailability();
  void SetState(State s);
  void SetCaptionText(QLabel* l, const QString& strText, const char* role);

  // running
  void OnPrimary();
  void OnCopyCommand();
  void OnShowOutput();
  // Shows this run (its log so far, and its state) and follows it.
  void Attach(ReconRunner* r);
  void OnRunFinished();
  void OnTick();
  QString ElapsedText(qint64 ms) const;

  State m_state = Idle;
  QString m_strRoot;
  QPointer<ReconRunner> m_pRunner;
  qint64 m_nFinalMs = 0;
  bool m_bLoading = false;
  bool m_bNarrow = false;
  QString m_strUnavailCopy;
  // The output folder the default SynthSR / SynthSeg paths point into.
  QString m_strRefDefaultsDir;
  // The stored settings had fields this version no longer has.
  bool m_bSettingsOutdated = false;

  QTimer* m_pTick = nullptr;
  QTimer* m_pValidate = nullptr;
  QTimer* m_pSave = nullptr;
  QTimer* m_pCopied = nullptr;

  QMap<QString, QLineEdit*> m_edits;
  QMap<QString, QCheckBox*> m_checks;
  QMap<QString, QLabel*> m_captions;
  QMap<QString, QString> m_help;
  QMap<QString, ReconFormRow*> m_rows;
  QList<ReconFormRow*> m_listRows;

  // left root
  QScrollArea* m_pScroll = nullptr;
  QWidget* m_pForm = nullptr;
  QLabel* m_pInventory = nullptr;
  SegmentedControl* m_pHemisphere = nullptr;
  SegmentedControl* m_pFace = nullptr;
  SegmentedControl* m_pOrder = nullptr;
  SegmentedControl* m_pTissue = nullptr;
  QRadioButton* m_pRefNone = nullptr;
  QRadioButton* m_pRefMesh = nullptr;
  QRadioButton* m_pRefMri = nullptr;
  QRadioButton* m_pRefBoth = nullptr;
  QWidget* m_pMeshBlock = nullptr;
  QWidget* m_pMriBlock = nullptr;
  QPushButton* m_pAdvancedToggle = nullptr;
  QWidget* m_pAdvanced = nullptr;
  QPushButton* m_pExpertToggle = nullptr;
  QWidget* m_pExpert = nullptr;
  QLabel* m_pReady = nullptr;
  QPushButton* m_pBack = nullptr;
  QPushButton* m_pRun = nullptr;

  // right root
  QPointer<QWidget> m_pPane;
  QPlainTextEdit* m_pPreview = nullptr;
  QLabel* m_pPreviewNote = nullptr;
  QPushButton* m_pCopy = nullptr;
  QLabel* m_pPhase = nullptr;
  QLabel* m_pElapsed = nullptr;
  QProgressBar* m_pProgress = nullptr;
  QLabel* m_pStatusLead = nullptr;
  QLabel* m_pStatus = nullptr;
  QLabel* m_pStatusDetail = nullptr;
  QPushButton* m_pShowOutput = nullptr;
  QPushButton* m_pUnavailCopy = nullptr;
  QPushButton* m_pCheckAgain = nullptr;
  QPlainTextEdit* m_pLog = nullptr;
  QVBoxLayout* m_pLogLayout = nullptr;
  QVBoxLayout* m_pStatusLayout = nullptr;
};

#endif // RECONPAGE_H
