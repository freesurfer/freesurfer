#ifndef PROJECTCARDWIDGET_H
#define PROJECTCARDWIDGET_H

#include <QFrame>
#include "ProjectManager.h"

class QLabel;
class QProgressBar;

// One project's tile in the Projects Hub grid: thumbnail, name, stage +
// progress, and a metadata line. Click to open; right-click for
// rename/remove/reveal.
class ProjectCardWidget : public QFrame
{
  Q_OBJECT
public:
  explicit ProjectCardWidget(const ProjectInfo& info, QWidget* parent = nullptr);

signals:
  void Opened(const QString& projectFolder);
  void RenameRequested(const QString& projectFolder);
  void RemoveRequested(const QString& projectFolder);
  void RevealRequested(const QString& projectFolder);
  void LocateRequested(const QString& projectFolder);

protected:
  void mouseReleaseEvent(QMouseEvent* e) override;
  void mouseDoubleClickEvent(QMouseEvent* e) override;
  void contextMenuEvent(QContextMenuEvent* e) override;

private:
  void LoadThumbnail();

  ProjectInfo m_info;
  QLabel* m_labelThumbnail;
  QLabel* m_labelName;
  QLabel* m_labelStage;
  QProgressBar* m_progressBar;
  QLabel* m_labelMeta;
  bool m_bMissing = false;
};

#endif // PROJECTCARDWIDGET_H
