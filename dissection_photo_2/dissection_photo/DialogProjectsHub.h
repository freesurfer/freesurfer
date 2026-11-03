#ifndef DIALOGPROJECTSHUB_H
#define DIALOGPROJECTSHUB_H

#include <QDialog>
#include "ProjectManager.h"

class QGridLayout;
class QLabel;
class QPushButton;
class QStackedWidget;
class QWidget;

// Launcher shown before MainWindow: lists known projects, resumes one, or
// creates a new named one. Persistence: see ProjectManager.
class DialogProjectsHub : public QDialog
{
  Q_OBJECT
public:
  explicit DialogProjectsHub(QWidget* parent = nullptr);

  // Valid after exec() returns Accepted: the root folder of the project to
  // load. Accepting always means a project was chosen or created.
  QString GetSelectedProjectFolder() const { return m_strSelectedFolder; }

private slots:
  void OnCardOpened(const QString& projectFolder);
  void OnCardRenameRequested(const QString& projectFolder);
  void OnCardRemoveRequested(const QString& projectFolder);
  void OnCardRevealRequested(const QString& projectFolder);
  void OnCardLocateRequested(const QString& projectFolder);
  void OnNewProject();
  void OnShowWelcome();
  void OnToggleDarkMode(bool bDark);

private:
  void RefreshGrid();
  // Re-applies this dialog's widget-level colours for the mode; they beat the
  // palette and the app-wide sheet.
  void Retheme(bool bDark);

  QStackedWidget* m_stack;
  QWidget* m_pageEmpty;
  QWidget* m_pageGrid;
  QGridLayout* m_gridLayout;
  QLabel* m_labelTitle;
  QLabel* m_labelEmptyBody;
  QPushButton* m_btnEmptyNew;
  QPushButton* m_btnSettings;
  QPushButton* m_btnDarkMode;
  QPushButton* m_btnHelp;

  QString m_strSelectedFolder;
};

#endif // DIALOGPROJECTSHUB_H
