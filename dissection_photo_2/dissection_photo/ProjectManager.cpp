#include "ProjectManager.h"

#include <QSettings>
#include "ReconPage.h"
#include <QVariantMap>
#include <QRegularExpression>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QDateTime>

namespace
{
const char* REGISTRY_KEY = "ProjectsHub/Paths";
const char* PROJECT_SUBDIR = ".dissection_photo";
const char* PROJECT_FILE = "project.json";
}

QStringList ProjectManager::GetRegisteredProjectPaths()
{
  QSettings s;
  QStringList list = s.value(REGISTRY_KEY).toStringList();
  list.removeDuplicates();
  return list;
}

void ProjectManager::RegisterProject(const QString& rootOutputFolder)
{
  if (rootOutputFolder.isEmpty())
    return;
  QSettings s;
  QStringList list = s.value(REGISTRY_KEY).toStringList();
  if (!list.contains(rootOutputFolder))
  {
    list.prepend(rootOutputFolder);
    s.setValue(REGISTRY_KEY, list);
  }
}

void ProjectManager::UnregisterProject(const QString& rootOutputFolder)
{
  QSettings s;
  QStringList list = s.value(REGISTRY_KEY).toStringList();
  if (list.removeAll(rootOutputFolder) > 0)
    s.setValue(REGISTRY_KEY, list);
}

QString ProjectManager::GetProjectFilePath(const QString& rootOutputFolder)
{
  return QDir(rootOutputFolder).absoluteFilePath(QString(PROJECT_SUBDIR) + "/" + PROJECT_FILE);
}

namespace
{
// Used when project.json is missing or unreadable (e.g. older app versions).
ProjectInfo InferProjectFromFolder(const QString& rootOutputFolder)
{
  ProjectInfo info;
  info.rootOutputFolder = rootOutputFolder;
  info.inputFolder = rootOutputFolder;
  info.name = QDir(rootOutputFolder).dirName();
  if (info.name.isEmpty())
    info.name = "Untitled Project";

  int nCC = QDir(rootOutputFolder + "/connected_components").entryInfoList(QDir::Files).size();
  int nMask = QDir(rootOutputFolder + "/mask").entryInfoList(QDir::Files).size();
  int nCorrected = QDir(rootOutputFolder + "/deformed").entryInfoList(QDir::Files).size();
  if (nCC > 0)
    info.currentStage = "connected_components";
  else if (nMask > 0)
    info.currentStage = "segmentation";
  else if (nCorrected > 0)
    info.currentStage = "correction";

  QFileInfo fi(rootOutputFolder);
  info.createdAt = fi.birthTime().isValid() ? fi.birthTime().toString(Qt::ISODate) : fi.lastModified().toString(Qt::ISODate);
  info.lastOpenedAt = info.createdAt;
  return info;
}
}

bool ProjectManager::LoadProject(const QString& rootOutputFolder, ProjectInfo& out)
{
  if (!QDir(rootOutputFolder).exists())
    return false;

  QFile file(GetProjectFilePath(rootOutputFolder));
  if (file.open(QFile::ReadOnly))
  {
    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(file.readAll(), &err);
    file.close();
    if (err.error == QJsonParseError::NoError && doc.isObject())
    {
      QJsonObject o = doc.object();
      out.name = o.value("project_name").toString();
      out.inputFolder = o.value("input_folder").toString();
      out.rootOutputFolder = rootOutputFolder;
      out.calibrationFile = o.value("calibration_file").toString();
      out.calibratedMode = o.value("calibrated_mode").toBool();
      out.createdAt = o.value("created_at").toString();
      out.lastOpenedAt = o.value("last_opened_at").toString();
      out.currentStage = o.value("current_stage").toString();
      out.inputImageCountAtCreation = o.value("input_image_count_at_creation").toInt();
      if (out.name.isEmpty())
        out.name = QDir(rootOutputFolder).dirName();
      return true;
    }
  }

  // Missing or corrupt -- fall back to inference and silently rewrite a
  // fresh project.json so subsequent loads are fast and consistent.
  out = InferProjectFromFolder(rootOutputFolder);
  SaveProject(out);
  return true;
}

QString ProjectManager::FolderNameFor(const QString& strName)
{
  QString str = strName.trimmed();
  str.replace(QRegularExpression("[\\\\/:*?\"<>|]"), "_");
  // A trailing dot or space is dropped by some file systems, and a leading dot
  // hides the folder.
  while (str.endsWith('.') || str.endsWith(' '))
    str.chop(1);
  while (str.startsWith('.'))
    str.remove(0, 1);
  return str.trimmed();
}

bool ProjectManager::RenameProject(const QString& strOldRoot, const QString& strNewName, QString* pNewRoot, QString* pError)
{
  auto fail = [pError](const QString& s) { if (pError) *pError = s; return false; };
  const QString strName = strNewName.trimmed();
  const QString strFolder = FolderNameFor(strName);
  if (strName.isEmpty() || strFolder.isEmpty())
    return fail("Enter a name for the project.");

  ProjectInfo info;
  if (!LoadProject(strOldRoot, info))
    return fail("The project folder could not be found.");
  if (ReconRunner* r = ReconRunner::ForRoot(strOldRoot))
    if (r->IsActive())
      return fail("A reconstruction is running for this project. Wait for it to finish, or stop it, before renaming.");

  const QString strOld = QDir::cleanPath(QFileInfo(strOldRoot).absoluteFilePath());
  const QString strParent = QFileInfo(strOld).absolutePath();
  const QString strNew = QDir(strParent).filePath(strFolder);

  if (strNew != strOld)
  {
    const bool bCaseOnly = strNew.compare(strOld, Qt::CaseInsensitive) == 0;
    // A case-only change finds the folder itself on a case-insensitive disk.
    if (!bCaseOnly && QFileInfo::exists(strNew))
      return fail(QString::fromUtf8("A folder named \u201c%1\u201d already exists in %2. Choose a different name for this project.")
                  .arg(strFolder, QDir::toNativeSeparators(strParent)));
    bool bOk;
    if (bCaseOnly)
    {
      // Through a temporary name, so the change also happens where case is ignored.
      const QString strTmp = strOld + ".renaming";
      bOk = QDir().rename(strOld, strTmp) && QDir().rename(strTmp, strNew);
      if (!bOk && QFileInfo::exists(strTmp))
        QDir().rename(strTmp, strOld);
    }
    else
      bOk = QDir().rename(strOld, strNew);
    if (!bOk)
      return fail(QString::fromUtf8("The folder could not be renamed. It may be open in another app, or you may not have "
                                    "permission to change %1.").arg(QDir::toNativeSeparators(strParent)));
  }

  // Paths inside the project move with it.
  auto moved = [&strOld, &strNew](const QString& p)
  {
    const QString c = QDir::cleanPath(p);
    if (c == strOld)
      return strNew;
    if (c.startsWith(strOld + "/"))
      return strNew + c.mid(strOld.length());
    return p;
  };
  info.name = strName;
  info.rootOutputFolder = strNew;
  info.inputFolder = moved(info.inputFolder);
  info.calibrationFile = moved(info.calibrationFile);

  if (strNew != strOld)
  {
    // Same place in the Hub's list, under the new path.
    QSettings st;
    QStringList list = st.value(REGISTRY_KEY).toStringList();
    bool bFound = false;
    for (QString& p : list)
    {
      if (QDir::cleanPath(QFileInfo(p).absoluteFilePath()) == strOld)
      {
        p = strNew;
        bFound = true;
      }
    }
    if (!bFound)
      list.prepend(strNew);
    list.removeDuplicates();
    st.setValue(REGISTRY_KEY, list);
    // The Reconstruction form's saved settings are filed under the folder path.
    const QString strOldKey = ReconPage::SettingsKeyFor(strOld), strNewKey = ReconPage::SettingsKeyFor(strNew);
    st.beginGroup(strOldKey);
    const QStringList keys = st.allKeys();
    QVariantMap values;
    for (const QString& k : keys)
      values[k] = st.value(k);
    st.endGroup();
    for (auto it = values.constBegin(); it != values.constEnd(); ++it)
      st.setValue(strNewKey + "/" + it.key(), it.value());
    st.remove(strOldKey);
  }
  SaveProject(info);
  if (pNewRoot)
    *pNewRoot = strNew;
  return true;
}

bool ProjectManager::SaveProject(const ProjectInfo& info)
{
  if (info.rootOutputFolder.isEmpty())
    return false;

  // mkpath() also creates a missing project folder.
  QDir().mkpath(QDir(info.rootOutputFolder).absoluteFilePath(PROJECT_SUBDIR));

  QJsonObject o;
  o["schema_version"] = 1;
  o["project_name"] = info.name;
  o["input_folder"] = info.inputFolder;
  o["root_output_folder"] = info.rootOutputFolder;
  o["calibration_file"] = info.calibrationFile;
  o["calibrated_mode"] = info.calibratedMode;
  o["created_at"] = info.createdAt.isEmpty() ? QDateTime::currentDateTime().toString(Qt::ISODate) : info.createdAt;
  o["last_opened_at"] = QDateTime::currentDateTime().toString(Qt::ISODate);
  o["current_stage"] = info.currentStage;
  o["input_image_count_at_creation"] = info.inputImageCountAtCreation;

  QFile file(GetProjectFilePath(info.rootOutputFolder));
  if (!file.open(QFile::WriteOnly | QFile::Truncate))
    return false;
  file.write(QJsonDocument(o).toJson(QJsonDocument::Indented));
  file.close();

  RegisterProject(info.rootOutputFolder);
  return true;
}

ProjectProgress ProjectManager::DeriveProgress(const ProjectInfo& info)
{
  ProjectProgress p;
  // An empty inputFolder is "Not started"; only a set path missing from disk is
  // "Missing".
  if (!info.inputFolder.isEmpty())
  {
    if (QDir(info.inputFolder).exists())
      p.inputCount = QDir(info.inputFolder).entryInfoList(QDir::Files).size();
    else
      p.inputFolderMissing = true;
  }
  if (p.inputCount == 0)
    p.inputCount = info.inputImageCountAtCreation;

  if (!QDir(info.rootOutputFolder).exists())
  {
    p.outputFolderMissing = true;
    return p;
  }

  p.correctedCount = QDir(info.rootOutputFolder + "/deformed").entryInfoList(QDir::Files).size();
  p.maskCount = QDir(info.rootOutputFolder + "/mask").entryList(QStringList("*.png"), QDir::Files).size();
  p.ccCount = QDir(info.rootOutputFolder + "/connected_components").entryInfoList(QDir::Files).size();
  p.reconDone = !QDir(info.rootOutputFolder + "/recon").entryList(QStringList("photo_recon*.mgz"), QDir::Files).isEmpty();
  return p;
}

QString ProjectManager::StageLabel(const ProjectInfo& info, const ProjectProgress& p)
{
  if (p.outputFolderMissing || (p.inputFolderMissing && p.inputCount == 0))
    return "Missing";

  int total = p.inputCount;
  if (p.reconDone)
    return "Complete";
  if (total > 0 && p.ccCount >= total)
    return "Reconstruction";
  if (info.currentStage == "connected_components" || p.ccCount > 0)
    return total > 0 ? QString("Connected Components — %1 of %2").arg(p.ccCount).arg(total) : "Connected Components";
  if (info.currentStage == "segmentation" || p.maskCount > 0)
    return total > 0 ? QString("Segmentation — %1 of %2").arg(p.maskCount).arg(total) : "Segmentation";
  if (info.currentStage == "correction" || p.correctedCount > 0)
    return total > 0 ? QString("Correction — %1 of %2").arg(p.correctedCount).arg(total) : "Correction";
  return "Not started";
}

double ProjectManager::StageFraction(const ProjectInfo& info, const ProjectProgress& p)
{
  // Whole-case progress: Correction, Segmentation, Connected components and
  // Reconstruction are a quarter each; earlier steps count as done and the
  // current one adds its share of photos.
  const int total = p.inputCount;
  if (p.reconDone)
    return 1.0;
  if (total > 0 && p.ccCount >= total)
    return 0.75;
  if (total <= 0)
    return 0.0;
  if (info.currentStage == "connected_components" || p.ccCount > 0)
    return 0.50 + 0.25 * qBound(0.0, p.ccCount / (double)total, 1.0);
  if (info.currentStage == "segmentation" || p.maskCount > 0)
    return 0.25 + 0.25 * qBound(0.0, p.maskCount / (double)total, 1.0);
  if (info.currentStage == "correction" || p.correctedCount > 0)
    return 0.25 * qBound(0.0, p.correctedCount / (double)total, 1.0);
  return 0.0;
}
