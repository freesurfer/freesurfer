#include "DialogProjectsHub.h"
#include "ProjectCardWidget.h"
#include "DialogWelcome.h"
#include "AppTheme.h"
#include <QFileInfo>
#include "DialogSegmentationSettings.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QScrollArea>
#include <QStackedWidget>
#include <QLabel>
#include <QPushButton>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QDesktopServices>
#include <QUrl>
#include <QDir>
#include <QDateTime>
#include <QLineEdit>
#include <QRegularExpression>
#include <algorithm>

namespace
{
const int CARD_W = 200;
const int GRID_COLUMNS = 4;

QString NewProjectTileStyle(bool bDark)
{
  if (bDark)
    return "QPushButton {background: #333333; border: 3px dashed #57a2c0; border-radius: 12px; color: #6fb3c9; font-size: 14px; font-weight: bold;}"
           "QPushButton:hover {background: #3b4448; border: 3px dashed #6fb3c9;}"
           "QPushButton:pressed {background: #2b3336;}";
  return "QPushButton {background: #fafafa; border: 3px dashed #57a2c0; border-radius: 12px; color: #457c92; font-size: 14px; font-weight: bold;}"
         "QPushButton:hover {background: #f0f7fa; border: 3px dashed #457c92;}"
         "QPushButton:pressed {background: #e6f0f4;}";
}

// Header link buttons: #457c92 is muddy on dark, so dark mode uses the lighter
// link colour from BuildPalette().
QString LinkButtonStyle(bool bDark)
{
  QString strColor = bDark ? "#6fb3c9" : "#457c92";
  return QString("QPushButton {background: transparent; border: none; color: %1; font-size: 13px;}"
                 "QPushButton:hover {text-decoration: underline;}").arg(strColor);
}
}

DialogProjectsHub::DialogProjectsHub(QWidget* parent) : QDialog(parent)
{
  setWindowTitle("Photo Reconstruction Projects");
  resize(900, 640);
  // Width 0: the layout's own minimum (real font metrics) applies. A fixed width
  // clipped the title with wider Linux fonts; the grid scrolls, so cards can't
  // raise this floor.
  setMinimumSize(0, 420);

  QVBoxLayout* mainLayout = new QVBoxLayout(this);
  mainLayout->setContentsMargins(24, 20, 24, 16);
  mainLayout->setSpacing(12);

  QHBoxLayout* header = new QHBoxLayout();
  m_labelTitle = new QLabel("Photo Reconstruction Projects", this);
  header->addWidget(m_labelTitle);
  header->addStretch();
  // Segmentation settings (folds, device, threads, model folder): a gear in the
  // same pill as the mode toggle beside it.
  m_btnSettings = new QPushButton(this);
  m_btnSettings->setCursor(Qt::PointingHandCursor);
  m_btnSettings->setToolTip("Automatic segmentation settings");
  m_btnSettings->setAccessibleName("Automatic segmentation settings");
  m_btnSettings->setIconSize(QSize(18, 18));
  connect(m_btnSettings, &QPushButton::clicked, this, [this]()
  {
    DialogSegmentationSettings dlg(this);
    dlg.exec();
  });
  header->addWidget(m_btnSettings);
  // The first screen on every launch, so the appearance override (AppTheme.h)
  // lives here. Defaults to the OS setting.
  m_btnDarkMode = new QPushButton(this);
  m_btnDarkMode->setCheckable(true);
  m_btnDarkMode->setChecked(AppTheme::EffectiveDarkMode());
  m_btnDarkMode->setCursor(Qt::PointingHandCursor);
  m_btnDarkMode->setText(m_btnDarkMode->isChecked() ? "☀  Light Mode" : "☽  Dark Mode");
  connect(m_btnDarkMode, &QPushButton::toggled, this, [this](bool bDark) {
    m_btnDarkMode->setText(bDark ? "☀  Light Mode" : "☽  Dark Mode");
    OnToggleDarkMode(bDark);
  });
  header->addWidget(m_btnDarkMode);
  m_btnHelp = new QPushButton("How This Works", this);
  m_btnHelp->setCursor(Qt::PointingHandCursor);
  connect(m_btnHelp, &QPushButton::clicked, this, &DialogProjectsHub::OnShowWelcome);
  header->addWidget(m_btnHelp);
  mainLayout->addLayout(header);

  m_stack = new QStackedWidget(this);
  mainLayout->addWidget(m_stack, 1);

  m_pageEmpty = new QWidget(this);
  QVBoxLayout* emptyLayout = new QVBoxLayout(m_pageEmpty);
  emptyLayout->addStretch();
  QLabel* emptyIcon = new QLabel(m_pageEmpty);
  QPixmap iconPix(":/resource/slabs@2x.png");
  if (!iconPix.isNull())
    emptyIcon->setPixmap(iconPix.scaled(80, 80, Qt::KeepAspectRatio, Qt::SmoothTransformation));
  emptyIcon->setAlignment(Qt::AlignCenter);
  emptyLayout->addWidget(emptyIcon);
  QLabel* emptyTitle = new QLabel("No projects yet", m_pageEmpty);
  emptyTitle->setAlignment(Qt::AlignCenter);
  emptyTitle->setStyleSheet("font-weight: bold; font-size: 16px; margin-top: 12px;");
  emptyLayout->addWidget(emptyTitle);
  m_labelEmptyBody = new QLabel("Start your first project to begin spatial correction, segmentation,\nand connected-components export.", m_pageEmpty);
  m_labelEmptyBody->setAlignment(Qt::AlignCenter);
  emptyLayout->addWidget(m_labelEmptyBody);
  m_btnEmptyNew = new QPushButton("+  Start Your First Project", m_pageEmpty);
  m_btnEmptyNew->setFixedSize(260, 60);
  m_btnEmptyNew->setCursor(Qt::PointingHandCursor);
  connect(m_btnEmptyNew, &QPushButton::clicked, this, &DialogProjectsHub::OnNewProject);
  QHBoxLayout* emptyBtnRow = new QHBoxLayout();
  emptyBtnRow->addStretch();
  emptyBtnRow->addWidget(m_btnEmptyNew);
  emptyBtnRow->addStretch();
  emptyLayout->addLayout(emptyBtnRow);
  emptyLayout->addStretch();
  m_stack->addWidget(m_pageEmpty);

  m_pageGrid = new QWidget(this);
  QVBoxLayout* gridPageLayout = new QVBoxLayout(m_pageGrid);
  gridPageLayout->setContentsMargins(0, 0, 0, 0);
  QScrollArea* scroll = new QScrollArea(m_pageGrid);
  scroll->setWidgetResizable(true);
  scroll->setFrameShape(QFrame::NoFrame);
  scroll->setStyleSheet("QScrollArea {background: transparent;}");
  QWidget* gridContainer = new QWidget();
  gridContainer->setStyleSheet("background: transparent;");
  m_gridLayout = new QGridLayout(gridContainer);
  m_gridLayout->setSpacing(16);
  m_gridLayout->setContentsMargins(4, 4, 4, 4);
  m_gridLayout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
  scroll->setWidget(gridContainer);
  gridPageLayout->addWidget(scroll);
  m_stack->addWidget(m_pageGrid);

  Retheme(AppTheme::EffectiveDarkMode());
  RefreshGrid();
}

void DialogProjectsHub::Retheme(bool bDark)
{
  // #262626 in light mode, the app's softened text colour.
  m_labelTitle->setStyleSheet(QString("font-weight: bold; font-size: 24px; color: %1;")
                              .arg(bDark ? "#f2f2f2" : "#262626"));
  m_labelEmptyBody->setStyleSheet(QString("font-size: 14px; color: %1; margin-bottom: 16px;")
                                  .arg(bDark ? "#9a9a9a" : "#666666"));
  m_btnEmptyNew->setStyleSheet(NewProjectTileStyle(bDark));
  m_btnDarkMode->setStyleSheet(AppTheme::ToggleStyleSheet(bDark));
  // The toggle's pill, narrower: just the gear, in the toggle's text colour.
  m_btnSettings->setStyleSheet(AppTheme::ToggleStyleSheet(bDark) + "QPushButton {padding: 5px 9px;}");
  m_btnSettings->setIcon(AppTheme::TintedIcon("settings", bDark ? QColor("#d8d8d8") : QColor("#4a4a4a")));
  m_btnHelp->setStyleSheet(LinkButtonStyle(bDark));
}

void DialogProjectsHub::RefreshGrid()
{
  QLayoutItem* item;
  while ((item = m_gridLayout->takeAt(0)) != nullptr)
  {
    if (item->widget())
      item->widget()->deleteLater();
    delete item;
  }

  QStringList paths = ProjectManager::GetRegisteredProjectPaths();

  QList<ProjectInfo> infos;
  foreach (const QString& path, paths)
  {
    ProjectInfo info;
    if (!ProjectManager::LoadProject(path, info))
    {
      // Folder gone: still listed (as Missing) so it can be located or removed.
      info.rootOutputFolder = path;
      info.name = QDir(path).dirName();
    }
    infos << info;
  }
  // Most-recently-opened first.
  std::sort(infos.begin(), infos.end(), [](const ProjectInfo& a, const ProjectInfo& b) {
    return QDateTime::fromString(a.lastOpenedAt, Qt::ISODate) > QDateTime::fromString(b.lastOpenedAt, Qt::ISODate);
  });

  if (infos.isEmpty())
  {
    m_stack->setCurrentWidget(m_pageEmpty);
    return;
  }
  m_stack->setCurrentWidget(m_pageGrid);

  QPushButton* newTile = new QPushButton("+\nNew Project");
  newTile->setFixedSize(CARD_W, 240);
  newTile->setCursor(Qt::PointingHandCursor);
  newTile->setStyleSheet(NewProjectTileStyle(AppTheme::EffectiveDarkMode()));
  connect(newTile, &QPushButton::clicked, this, &DialogProjectsHub::OnNewProject);
  m_gridLayout->addWidget(newTile, 0, 0);

  int row = 0, col = 1;
  foreach (const ProjectInfo& info, infos)
  {
    ProjectCardWidget* card = new ProjectCardWidget(info);
    connect(card, &ProjectCardWidget::Opened, this, &DialogProjectsHub::OnCardOpened);
    connect(card, &ProjectCardWidget::RenameRequested, this, &DialogProjectsHub::OnCardRenameRequested);
    connect(card, &ProjectCardWidget::RemoveRequested, this, &DialogProjectsHub::OnCardRemoveRequested);
    connect(card, &ProjectCardWidget::RevealRequested, this, &DialogProjectsHub::OnCardRevealRequested);
    connect(card, &ProjectCardWidget::LocateRequested, this, &DialogProjectsHub::OnCardLocateRequested);
    m_gridLayout->addWidget(card, row, col);
    col++;
    if (col >= GRID_COLUMNS)
    {
      col = 0;
      row++;
    }
  }
}

void DialogProjectsHub::OnCardOpened(const QString& projectFolder)
{
  ProjectInfo info;
  if (ProjectManager::LoadProject(projectFolder, info))
  {
    info.lastOpenedAt = QDateTime::currentDateTime().toString(Qt::ISODate);
    ProjectManager::SaveProject(info);
  }
  m_strSelectedFolder = projectFolder;
  accept();
}

namespace
{
// Asks for the new name, showing the folder it will become. Returns false on
// Cancel; strName holds what was typed either way.
bool AskProjectName(QWidget* parent, const QString& strFolderParent, QString& strName)
{
  const bool bDark = AppTheme::EffectiveDarkMode();
  QDialog dlg(parent);
  dlg.setWindowTitle("Rename Project");
  dlg.setMinimumWidth(460);
  QVBoxLayout* v = new QVBoxLayout(&dlg);
  v->setContentsMargins(24, 20, 24, 18);
  v->setSpacing(8);
  QLabel* title = new QLabel("Rename project");
  title->setStyleSheet(QString("font-weight: bold; font-size: 18px; color: %1;").arg(bDark ? "#f2f2f2" : "#262626"));
  v->addWidget(title);
  QLineEdit* edit = new QLineEdit(strName);
  edit->setObjectName("renameEdit");
  edit->setStyleSheet(QString("#renameEdit {background: %1; border: 1px solid %2; border-radius: 4px; padding: 4px 8px;"
                              " font-size: 14px; color: %3;} #renameEdit:focus {border: 1px solid #457c92;}")
                      .arg(bDark ? "#333333" : "#ffffff", bDark ? "#5f5f5f" : "#c0c0c0", bDark ? "#f2f2f2" : "#262626"));
  edit->selectAll();
  v->addWidget(edit);
  QLabel* caption = new QLabel;
  caption->setWordWrap(true);
  caption->setTextInteractionFlags(Qt::TextSelectableByMouse);
  caption->setStyleSheet(QString("font-size: 12px; color: %1;").arg(bDark ? "#9a9a9a" : "#666666"));
  v->addWidget(caption);
  v->addSpacing(8);
  QHBoxLayout* h = new QHBoxLayout;
  h->addStretch(1);
  QPushButton* cancel = new QPushButton("Cancel");
  QPushButton* ok = new QPushButton("Rename");
  const QString strGhost = QString("QPushButton {background: transparent; padding: 0 14px; font-size: 13px; border: 1px solid %1;"
                                   " border-radius: 14px; color: %2; min-height: 28px; max-height: 28px;}"
                                   "QPushButton:hover {background: %3;}")
                           .arg(bDark ? "#777777" : "#999999", bDark ? "#dddddd" : "#555555", bDark ? "#4a4a4a" : "#e4e4e4");
  cancel->setStyleSheet(strGhost);
  ok->setStyleSheet("QPushButton {background: #457c92; padding: 0 22px; font-size: 13px; font-weight: bold; border: none;"
                    " border-radius: 14px; color: white; min-height: 28px; max-height: 28px;}"
                    "QPushButton:hover {background: #3c6d80;} QPushButton:disabled {background: #9fb7c1;}");
  ok->setDefault(true);
  h->addWidget(cancel);
  h->addWidget(ok);
  v->addLayout(h);
  QObject::connect(cancel, &QPushButton::clicked, &dlg, &QDialog::reject);
  QObject::connect(ok, &QPushButton::clicked, &dlg, &QDialog::accept);
  auto update = [&]()
  {
    const QString strFolder = ProjectManager::FolderNameFor(edit->text());
    ok->setEnabled(!strFolder.isEmpty());
    caption->setText(strFolder.isEmpty() ? QString("Enter a name for the project.")
        : QString::fromUtf8("Its folder becomes \u201c%1\u201d in %2.").arg(strFolder, QDir::toNativeSeparators(strFolderParent)));
  };
  QObject::connect(edit, &QLineEdit::textChanged, &dlg, update);
  update();
  const bool bAccepted = dlg.exec() == QDialog::Accepted;
  strName = edit->text().trimmed();
  return bAccepted;
}
}

// Why a rename didn't happen, styled like the Rename window (a QMessageBox wraps
// a long path into a narrow column).
void ShowRenameProblem(QWidget* parent, const QString& strText)
{
  const bool bDark = AppTheme::EffectiveDarkMode();
  QDialog dlg(parent);
  dlg.setWindowTitle("Can't Rename the Project");
  QVBoxLayout* v = new QVBoxLayout(&dlg);
  v->setContentsMargins(24, 20, 24, 18);
  v->setSpacing(8);
  QLabel* title = new QLabel("Can't rename the project");
  title->setStyleSheet(QString("font-weight: bold; font-size: 18px; color: %1;").arg(bDark ? "#f2f2f2" : "#262626"));
  v->addWidget(title);
  QLabel* body = new QLabel(strText);
  body->setWordWrap(true);
  body->setFixedWidth(420);
  body->setTextInteractionFlags(Qt::TextSelectableByMouse);
  body->setStyleSheet(QString("font-size: 14px; color: %1;").arg(bDark ? "#d8d8d8" : "#3a3a3a"));
  v->addWidget(body);
  v->addSpacing(8);
  QHBoxLayout* h = new QHBoxLayout;
  h->addStretch(1);
  QPushButton* ok = new QPushButton("OK");
  ok->setStyleSheet("QPushButton {background: #457c92; padding: 0 22px; font-size: 13px; font-weight: bold; border: none;"
                    " border-radius: 14px; color: white; min-height: 28px; max-height: 28px;}"
                    "QPushButton:hover {background: #3c6d80;}");
  ok->setDefault(true);
  h->addWidget(ok);
  v->addLayout(h);
  QObject::connect(ok, &QPushButton::clicked, &dlg, &QDialog::accept);
  dlg.exec();
}

void DialogProjectsHub::OnCardRenameRequested(const QString& projectFolder)
{
  ProjectInfo info;
  if (!ProjectManager::LoadProject(projectFolder, info))
    return;
  const QString strParent = QFileInfo(QDir::cleanPath(projectFolder)).absolutePath();
  QString strName = info.name;
  // Until it works or is cancelled: a refusal (e.g. the folder name is taken)
  // explains itself and reopens with what was typed.
  while (AskProjectName(this, strParent, strName))
  {
    if (strName == info.name && ProjectManager::FolderNameFor(strName) == QDir(projectFolder).dirName())
      return;
    QString strError;
    if (ProjectManager::RenameProject(projectFolder, strName, nullptr, &strError))
    {
      RefreshGrid();
      return;
    }
    ShowRenameProblem(this, strError);
  }
}

void DialogProjectsHub::OnCardRemoveRequested(const QString& projectFolder)
{
  ProjectInfo info;
  ProjectManager::LoadProject(projectFolder, info);
  QString name = info.name.isEmpty() ? QDir(projectFolder).dirName() : info.name;
  if (QMessageBox::question(this, "Remove Project",
        QString("Remove \"%1\" from Projects? The photos and files on disk will not be touched.").arg(name),
        QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
  {
    ProjectManager::UnregisterProject(projectFolder);
    RefreshGrid();
  }
}

void DialogProjectsHub::OnCardRevealRequested(const QString& projectFolder)
{
  QDesktopServices::openUrl(QUrl::fromLocalFile(projectFolder));
}

void DialogProjectsHub::OnCardLocateRequested(const QString& projectFolder)
{
  QString newPath = QFileDialog::getExistingDirectory(this, "Locate Project Folder", "");
  if (newPath.isEmpty())
    return;
  ProjectInfo info;
  info.rootOutputFolder = projectFolder; // registry identity stays the same
  info.name = QDir(projectFolder).dirName();
  ProjectManager::LoadProject(projectFolder, info); // fill in what we can
  info.inputFolder = newPath;
  ProjectManager::SaveProject(info);
  RefreshGrid();
}

void DialogProjectsHub::OnNewProject()
{
  // Only name and folder are asked here; photos are added afterward on the
  // Start page (MainWindow::SetInputFolder) and copied into <project
  // folder>/photos, so the project folder is always a location the user chose.
  bool ok = false;
  QString name = QInputDialog::getText(this, "Name Your Project", "Project name:", QLineEdit::Normal, "", &ok).trimmed();
  if (!ok)
    return;
  if (name.isEmpty())
    name = "Untitled Project";

  QString safeName = name;
  safeName.replace(QRegularExpression("[\\\\/:*?\"<>|]"), "_");

  // A Save-style panel (location plus folder name, pre-filled with the project
  // name): this creates a new folder, and Save panels can do that natively.
  QFileDialog dlg(this, "Save Project Folder");
  dlg.setAcceptMode(QFileDialog::AcceptSave);
  dlg.setFileMode(QFileDialog::Directory);
  dlg.setOption(QFileDialog::ShowDirsOnly, false);
  dlg.selectFile(safeName);
  if (dlg.exec() != QDialog::Accepted || dlg.selectedFiles().isEmpty())
    return;
  QString projectFolder = dlg.selectedFiles().first();

  if (QFileInfo::exists(projectFolder) && !QDir(projectFolder).isEmpty())
  {
    QMessageBox::warning(this, "Error",
        QString("\"%1\" already exists and isn't empty. Choose a different name or location.").arg(projectFolder));
    return;
  }
  if (!QDir().mkpath(projectFolder))
  {
    QMessageBox::warning(this, "Error", QString("Could not create the project folder \"%1\".").arg(projectFolder));
    return;
  }

  // Disambiguate a name collision with an existing project -- cosmetic
  // only, the registry key is the folder path, not the name.
  QStringList existingNames;
  foreach (const QString& p, ProjectManager::GetRegisteredProjectPaths())
  {
    ProjectInfo existing;
    if (ProjectManager::LoadProject(p, existing))
      existingNames << existing.name;
  }
  QString finalName = name;
  int suffix = 2;
  while (existingNames.contains(finalName))
    finalName = QString("%1 (%2)").arg(name).arg(suffix++);

  ProjectInfo info;
  info.name = finalName;
  info.rootOutputFolder = projectFolder;
  info.createdAt = QDateTime::currentDateTime().toString(Qt::ISODate);
  info.lastOpenedAt = info.createdAt;
  ProjectManager::SaveProject(info);

  m_strSelectedFolder = info.rootOutputFolder;
  accept();
}

void DialogProjectsHub::OnShowWelcome()
{
  DialogWelcome dlg(this);
  dlg.exec();

  // The Welcome dialog has its own mode toggle, so re-theme on return (widget
  // sheets aren't re-applied otherwise). Signals are blocked so syncing the
  // button doesn't persist the mode as an explicit choice.
  const bool bDark = AppTheme::EffectiveDarkMode();
  m_btnDarkMode->blockSignals(true);
  m_btnDarkMode->setChecked(bDark);
  m_btnDarkMode->setText(bDark ? "☀  Light Mode" : "☽  Dark Mode");
  m_btnDarkMode->blockSignals(false);
  Retheme(bDark);
  RefreshGrid();
}

void DialogProjectsHub::OnToggleDarkMode(bool bDark)
{
  // Persisted as an explicit override, so MainWindow picks it up via
  // EffectiveDarkMode(). Apply(), not setPalette(): the app-wide sheet must be
  // swapped too (see AppTheme::BaseStyleSheet()). RefreshGrid() re-themes cards.
  AppTheme::SetMode(bDark ? "dark" : "light");
  AppTheme::Apply(bDark);
  Retheme(bDark);
  RefreshGrid();
}
