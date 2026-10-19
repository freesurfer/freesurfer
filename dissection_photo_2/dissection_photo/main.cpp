#include "MainWindow.h"
#include "DialogWelcome.h"
#include "DialogProjectsHub.h"
#include "DialogSegmentationSettings.h"
#include "ReconPage.h"
#include "AppTheme.h"

#include <QIcon>
#include <QApplication>
#include <QSettings>
#include <QStyleFactory>
#include <QMessageBox>
#include <QDir>
#include <QProcessEnvironment>
#include <QDebug>

namespace
{
// Same checks as dissection_photo_gui.sh / mac_wrapper, so a directly launched
// binary also fails fast with a clear message. Messages are also printed to the
// terminal.
bool CheckEnvironmentOrExplain()
{
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();

  QString strFreesurferHome = env.value("FREESURFER_HOME");
  if (strFreesurferHome.isEmpty() || !QDir(strFreesurferHome).exists())
  {
    qCritical().noquote() << "error: FREESURFER_HOME has not been properly set";
    QMessageBox::critical(nullptr, "FreeSurfer Not Sourced",
        "FREESURFER_HOME has not been properly set.\n\n"
        "This application needs FreeSurfer to be sourced first (e.g. via "
        "SetUpFreeSurfer.sh) -- without it, the required Python files "
        "cannot be found.");
    return false;
  }

  // The model folder from Automatic segmentation settings (default: NNUNET_MODEL_DIR,
  // else FreeSurfer's). Not checked when automatic segmentation is off.
  const QString strModelDir = DialogSegmentationSettings::ModelDir();
  if (DialogSegmentationSettings::Enabled() && !DialogSegmentationSettings::ModelFound(strModelDir))
  {
    QString strDetails = QString(
        "You can choose the model folder, or turn automatic segmentation off, in Automatic "
        "segmentation settings (the gear on the Projects screen).\n"
        "\n"
        "Unable to detect the required files in NNUNET_MODEL_DIR: %1\n"
        "If this is your first time using this util, follow the instructions\n"
        "below to download, unpack, install, and configure the model files.\n"
        "\n"
        "Download the model file from the FTP server. Copying the following link\n"
        "into your browser will begin the download: https://ftp.nmr.mgh.harvard.edu/pub/dist/lcnpublic/dist/dissection_photo_model/nnUNetTrainer__nnUNetPlans__2d.tar.gz\n"
        "\n"
        "Unpack the the downloaded tarball:\n"
        "  tar -xzf nnUNetTrainer__nnUNetPlans__2d.tar.gz\n"
        "\n"
        "Move the extracted directory to your models/ directory in your FS\n"
        "install. If FS was built with FSPYTHON_INSTALL_TREE=OFF, the models/\n"
        "dir will be located at $FREESURFER_HOME/models\n"
        "If FSPYTHON_INSTALL_TREE=ON, then models/ will be located at\n"
        "$FREESURFER_HOME/../fspython/python/models\n"
        "\n"
        "Once the model files are unpacked and installed into the models/dir,\n"
        "set NNUNET_MODEL_DIR to be the root of the directory you just moved\n"
        "into the models/ dir.\n"
        "\n"
        "Then set NNUNET_SCRIPTS_DIR to be the path to the nnUNet_v2 dir under\n"
        "the packages/ directory.\n"
        "If FSPYTHON_INSTALL_TREE=OFF was set during the build, this dir can\n"
        "be found under $FREESURFER_HOME/python/packages.\n"
        "If FSPYTHON_INSTALL_TREE=ON was set during the build, this dir can be\n"
        "under $FREESURFER_HOME/../fspython/python/packages").arg(strModelDir);

    // A warning, not a stop: the settings that fix it are on the Projects screen.
    qWarning().noquote() << strDetails;
    QMessageBox::warning(nullptr, "nnUNet Model Files Not Found", strDetails);
  }

  return true;
}
}

int main(int argc, char *argv[])
{
  qputenv("QT_AUTO_SCREEN_SCALE_FACTOR", "1");
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
  // No-op (and deprecated) on Qt 6; needed on Qt 5, where scaling is off.
  QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif
  QApplication a(argc, argv);
  // The Dock icon while running (and every window's), bundled or not.
  a.setWindowIcon(QIcon(":/resource/dissection_photo_icon.png"));

  // Fail before any window shows if FreeSurfer or the nnUNet model is missing.
  if (!CheckEnvironmentOrExplain())
    return 1;

  // Fusion on every platform, so unstyled widgets look the same on macOS and
  // Linux. Light/dark still follows AppTheme.
  QApplication::setStyle(QStyleFactory::create("Fusion"));

  // Query the OS setting before BuildPalette() overrides the app palette (the
  // Qt < 6.5 fallback reads it). Mode is "system" by default, or the user's
  // light/dark override. Apply() sets palette and stylesheet together: the
  // stylesheet would otherwise override the dark palette.
  AppTheme::Apply(AppTheme::EffectiveDarkMode());

  QSettings s;
  if (!s.value("WelcomeWindow/DoNotShow").toBool())
  {
    DialogWelcome dlg;
    dlg.exec();
    if (dlg.GetDoNotShow())
      s.setValue("WelcomeWindow/DoNotShow", true);
  }

  // Hub -> MainWindow repeats: File > Switch Project returns to the Hub instead
  // of quitting. Sequential exec() calls are safe.
  int ret = 0;
  bool bRestart = true;
  while (bRestart)
  {
    bRestart = false;

    DialogProjectsHub hub;
    if (hub.exec() != QDialog::Accepted)
    {
      // User closed the Hub itself -- quit the app, unless a reconstruction left
      // running by Switch Project should keep going (back to the Hub).
      if (ReconRunner::ConfirmQuit(nullptr))
        break;
      bRestart = true;
      continue;
    }

    MainWindow w;
    // Connected before LoadProject(), which can synchronously ask to return to
    // the Hub (Back in the Resume dialog).
    QObject::connect(&w, &MainWindow::RequestSwitchProject, [&bRestart]() { bRestart = true; });
    // An accepted Hub always means a chosen project; LoadProject() also
    // ignores an empty path.
    w.LoadProject(hub.GetSelectedProjectFolder());
    // Don't show the window if LoadProject() already asked to return to the
    // Hub.
    if (!bRestart)
    {
      w.show();
      ret = a.exec();
    }
  }
  ReconRunner::StopAll();
  return ret;
}
