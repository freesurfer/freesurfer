#ifndef PROJECTMANAGER_H
#define PROJECTMANAGER_H

#include <QString>
#include <QStringList>

// One project = one photo folder taken through the pipeline. Identity is the
// root output folder (rootOutputFolder), which the app owns; the input folder
// may be shared or read-only.
struct ProjectInfo
{
  QString name;
  QString inputFolder;
  QString rootOutputFolder;
  QString calibrationFile;
  bool    calibratedMode = false;
  QString createdAt;      // ISO 8601, QDateTime::toString(Qt::ISODate)
  QString lastOpenedAt;
  // "correction" | "segmentation" | "connected_components" | "" (not started)
  QString currentStage;
  int     inputImageCountAtCreation = 0;
};

// Progress derived from file counts, so projects changed outside the app don't
// show a stale cached stage.
struct ProjectProgress
{
  int inputCount = 0;
  int correctedCount = 0;
  int maskCount = 0;
  int ccCount = 0;
  bool reconDone = false;   // a photo_recon*.mgz exists in <root>/recon
  bool inputFolderMissing = false;
  bool outputFolderMissing = false;
};

namespace ProjectManager
{
  // Registry of known projects (root output folder paths) in QSettings, so the
  // Hub needn't scan the filesystem. Re-registering is safe.
  QStringList GetRegisteredProjectPaths();
  void RegisterProject(const QString& rootOutputFolder);
  void UnregisterProject(const QString& rootOutputFolder);

  // Per-project state lives at <rootOutputFolder>/.dissection_photo/project.json.
  QString GetProjectFilePath(const QString& rootOutputFolder);

  // Loads project state. If project.json is missing or unreadable, infers it
  // from the folder (name from basename, stage from populated outputs).
  // Returns false (out untouched) only if the folder no longer exists.
  bool LoadProject(const QString& rootOutputFolder, ProjectInfo& out);

  // Writes project.json, creating the hidden directory as needed.
  bool SaveProject(const ProjectInfo& info);

  // Counts files in the project's known output subfolders against its
  // input count.
  ProjectProgress DeriveProgress(const ProjectInfo& info);

  // The folder name a project name gives: characters a folder name can't hold
  // become "_" (as New Project does).
  QString FolderNameFor(const QString& strName);

  // Renames the project and its folder (same parent): the folder on disk, the
  // paths kept in project.json, the Hub's list (same place) and the
  // Reconstruction form's saved settings. Returns false with *pError set, and
  // nothing changed, when it can't: no name, a folder of that name already
  // there, a reconstruction running, or the rename refused by the system.
  bool RenameProject(const QString& strOldRoot, const QString& strNewName, QString* pNewRoot, QString* pError);

  // Stage summary and 0..1 completion from DeriveProgress(), falling back to
  // info.currentStage when counts are ambiguous (e.g. a stage with no files yet).
  QString StageLabel(const ProjectInfo& info, const ProjectProgress& progress);
  double  StageFraction(const ProjectInfo& info, const ProjectProgress& progress);
}

#endif // PROJECTMANAGER_H
