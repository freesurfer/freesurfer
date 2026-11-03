#ifndef DIALOGSEGMENTATIONSETTINGS_H
#define DIALOGSEGMENTATIONSETTINGS_H

// Settings for the automatic segmentation (the nnUNet daemon), opened from the
// gear on the Projects Hub.
//
// Values are kept in QSettings under "Segmentation/"; the Key*() names and the
// stored values ("auto"/"1"/"5", "auto"/"gpu"/"cpu", a thread count) match the
// daemon's --folds, --device and --max_threads options.

#include <QDialog>

class QLabel;
class QLineEdit;
class QPushButton;
class SegmentedControl;

class DialogSegmentationSettings : public QDialog
{
  Q_OBJECT
public:
  explicit DialogSegmentationSettings(QWidget* parent = nullptr);

  static const char* KeyEnabled()  { return "Segmentation/Enabled"; }
  static const char* KeyFolds()    { return "Segmentation/Folds"; }
  static const char* KeyDevice()   { return "Segmentation/Device"; }
  static const char* KeyModelDir() { return "Segmentation/ModelDir"; }
  static const char* KeyMaxThreads() { return "Segmentation/MaxThreads"; }

  // The saved values, or the defaults: on, automatic folds (five on a graphics
  // card, one on the processor), automatic device, and DefaultModelDir(). Off
  // means every mask is drawn by hand.
  static bool Enabled();
  static QString Folds();
  // NNUNET_MODEL_DIR, else the model installed with FreeSurfer.
  static QString DefaultModelDir();
  // The folder holds an nnUNet model (plans.json).
  static bool ModelFound(const QString& strDir);
  static QString Device();
  static QString ModelDir();
  // The most processor threads segmentation may use (default 8, the daemon's cap).
  static int MaxThreads();
  // A GPU the daemon can use: Apple silicon, or an NVIDIA driver.
  static bool GpuPresent();

protected:
  void showEvent(class QShowEvent* e) override;

private:
  void FitCaptions();
  void Save();
  void OnChooseFolder();
  void UpdateModelCaption();
  void UpdateThreadsCaption();
  // Asks before turning automatic segmentation off; updates the rows below it.
  void OnEnabledChanged(int nIndex);
  void UpdateEnabled();
  void Retheme(bool bDark);
  QLabel* AddRow(class QGridLayout* grid, int nRow, const QString& strLabel, QWidget* pControl, const QString& strCaption);

  SegmentedControl* m_pEnabled = nullptr;
  QLabel* m_pEnabledCaption = nullptr;
  SegmentedControl* m_pFolds = nullptr;
  SegmentedControl* m_pDevice = nullptr;
  QLineEdit* m_pModelDir = nullptr;
  QPushButton* m_pChoose = nullptr;
  QLabel* m_pTitle = nullptr;
  QLabel* m_pSubtitle = nullptr;
  QLabel* m_pModelCaption = nullptr;
  class QSpinBox* m_pThreads = nullptr;
  QLabel* m_pThreadsCaption = nullptr;
  QList<QLabel*> m_listLabels;
  QList<QLabel*> m_listCaptions;
};

#endif // DIALOGSEGMENTATIONSETTINGS_H
