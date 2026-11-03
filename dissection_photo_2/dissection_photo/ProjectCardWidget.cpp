#include "ProjectCardWidget.h"
#include "AppTheme.h"
#include "ReconPage.h"

#include <QLabel>
#include <QProgressBar>
#include <QVBoxLayout>
#include <QDir>
#include <QImageReader>
#include <QMouseEvent>
#include <QMenu>
#include <QDateTime>
#include <QFontMetrics>
#include <QTimer>

ProjectCardWidget::ProjectCardWidget(const ProjectInfo& info, QWidget* parent)
  : QFrame(parent), m_info(info)
{
  setFixedSize(200, 240);
  setCursor(Qt::PointingHandCursor);
  setObjectName("ProjectCard");

  // Colours are per mode (a card's own sheet beats the palette). Cards are
  // rebuilt on every toggle, so reading the mode once is enough.
  const bool bDark = AppTheme::EffectiveDarkMode();
  const QString strCardBg     = bDark ? "#353535" : "#ffffff";
  const QString strCardBorder = bDark ? "#4a4a4a" : "#cccccc";
  const QString strCardHover  = bDark ? "#3d3d3d" : "#fafafa";
  const QString strHoverEdge  = bDark ? "#6fb3c9" : "#457c92";
  const QString strThumbBg    = bDark ? "#262626" : "#f0f0f0";
  const QString strNameColor  = bDark ? "#f2f2f2" : "#262626";
  const QString strAccent     = bDark ? "#6fb3c9" : "#457c92";
  const QString strTrackBg    = bDark ? "#2a2a2a" : "#e6e6e6";
  const QString strMetaColor  = bDark ? "#9a9a9a" : "#666666";

  setStyleSheet(QString("QFrame#ProjectCard {background: %1; border: 1px solid %2; border-radius: 10px;}"
                        "QFrame#ProjectCard:hover {background: %3; border: 1px solid %4;}")
                .arg(strCardBg, strCardBorder, strCardHover, strHoverEdge));

  QVBoxLayout* layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  layout->setSpacing(0);

  m_labelThumbnail = new QLabel(this);
  m_labelThumbnail->setFixedSize(198, 120);
  m_labelThumbnail->setAlignment(Qt::AlignCenter);
  m_labelThumbnail->setStyleSheet(QString("background: %1; border-top-left-radius: 9px; border-top-right-radius: 9px;")
                                  .arg(strThumbBg));
  layout->addWidget(m_labelThumbnail);

  QWidget* body = new QWidget(this);
  QVBoxLayout* bodyLayout = new QVBoxLayout(body);
  bodyLayout->setContentsMargins(12, 10, 12, 10);
  bodyLayout->setSpacing(4);

  m_labelName = new QLabel(body);
  // #262626 in light mode, the app's softened text colour.
  m_labelName->setStyleSheet(QString("font-weight: bold; font-size: 14px; color: %1;").arg(strNameColor));
  bodyLayout->addWidget(m_labelName);

  m_labelStage = new QLabel(body);
  m_labelStage->setStyleSheet(QString("font-size: 12px; color: %1;").arg(strAccent));
  bodyLayout->addWidget(m_labelStage);

  m_progressBar = new QProgressBar(body);
  m_progressBar->setFixedHeight(6);
  m_progressBar->setTextVisible(false);
  m_progressBar->setStyleSheet(QString("QProgressBar {border: none; background: %1; border-radius: 3px;}"
                                       "QProgressBar::chunk {background: %2; border-radius: 3px;}")
                               .arg(strTrackBg, strAccent));
  bodyLayout->addWidget(m_progressBar);

  bodyLayout->addStretch();

  m_labelMeta = new QLabel(body);
  m_labelMeta->setStyleSheet(QString("font-size: 12px; color: %1;").arg(strMetaColor));
  bodyLayout->addWidget(m_labelMeta);

  layout->addWidget(body);

  ProjectProgress progress = ProjectManager::DeriveProgress(m_info);
  m_bMissing = progress.outputFolderMissing || (progress.inputFolderMissing && progress.inputCount == 0);

  QFontMetrics fm(m_labelName->font());
  m_labelName->setText(fm.elidedText(m_info.name.isEmpty() ? "Untitled Project" : m_info.name, Qt::ElideRight, 176));
  m_labelName->setToolTip(m_info.name);

  if (m_bMissing)
  {
    // Brand red (#a2293c) is too dark to read on the dark card; lightened
    // toward the same hue instead of swapped for a different one.
    m_labelStage->setText("⚠ Folder not found");
    m_labelStage->setStyleSheet(QString("font-size: 12px; color: %1;").arg(bDark ? "#e07b8a" : "#a2293c"));
    m_progressBar->setVisible(false);
    // A shade off the normal card fill, marking it as inert -- lighter in
    // light mode, darker in dark mode (both directions read as "recessed").
    setStyleSheet(styleSheet() + QString("QFrame#ProjectCard {background: %1;}").arg(bDark ? "#2e2e2e" : "#f7f7f7"));
  }
  else
  {
    QString label = ProjectManager::StageLabel(m_info, progress);
    m_labelStage->setText(label);
    if (label == "Complete")
      m_labelStage->setStyleSheet(QString("font-size: 12px; color: %1; font-weight: bold;").arg(strAccent));
    m_progressBar->setValue((int)(ProjectManager::StageFraction(m_info, progress) * 100));

    // A reconstruction left running by Switch Project, or finished since: kept
    // current while the Hub is open.
    if (ReconRunner::ForRoot(m_info.rootOutputFolder))
    {
      const QString strError = bDark ? "#e07b8a" : "#a2293c";
      auto update = [this, strAccent, strError, label]()
      {
        ReconRunner* r = ReconRunner::ForRoot(m_info.rootOutputFolder);
        if (!r)
          return;
        const qint64 s = r->ElapsedMs() / 1000;
        const QString strFor = s < 60 ? QString("%1 s").arg(s) : QString("%1 min").arg(s / 60);
        QString strText = label, strColor = strAccent;
        switch (r->GetOutcome())
        {
        case ReconRunner::Active:    strText = QString::fromUtf8("Reconstructing\xe2\x80\xa6 %1").arg(strFor); break;
        case ReconRunner::Succeeded: strText = "Reconstruction finished"; break;
        case ReconRunner::Failed:    strText = "Reconstruction did not finish"; strColor = strError; break;
        case ReconRunner::Stopped:   strText = "Reconstruction stopped"; break;
        }
        m_labelStage->setText(strText);
        m_labelStage->setStyleSheet(QString("font-size: 12px; color: %1; font-weight: bold;").arg(strColor));
        // Moving while it runs.
        if (r->IsActive())
          m_progressBar->setRange(0, 0);
        else
        {
          m_progressBar->setRange(0, 100);
          m_progressBar->setValue(100);
        }
      };
      update();
      QTimer* pTimer = new QTimer(this);
      connect(pTimer, &QTimer::timeout, this, update);
      pTimer->start(1000);
    }
  }

  QString metaText;
  if (progress.inputCount > 0)
    metaText = QString("%1 image%2").arg(progress.inputCount).arg(progress.inputCount == 1 ? "" : "s");
  QDateTime lastOpened = QDateTime::fromString(m_info.lastOpenedAt, Qt::ISODate);
  if (lastOpened.isValid())
  {
    qint64 days = lastOpened.daysTo(QDateTime::currentDateTime());
    QString when = days <= 0 ? "today" : (days == 1 ? "yesterday" : QString("%1 days ago").arg(days));
    metaText += (metaText.isEmpty() ? "" : " · ") + QString("opened %1").arg(when);
  }
  m_labelMeta->setText(metaText);

  LoadThumbnail();
}

void ProjectCardWidget::LoadThumbnail()
{
  QString firstImage;
  if (!m_bMissing && !m_info.inputFolder.isEmpty() && QDir(m_info.inputFolder).exists())
  {
    QFileInfoList files = QDir(m_info.inputFolder).entryInfoList(QDir::Files, QDir::Name);
    if (!files.isEmpty())
      firstImage = files.first().absoluteFilePath();
  }

  QPixmap pix;
  if (!firstImage.isEmpty())
  {
    QImageReader reader(firstImage);
    QSize sz = reader.size();
    if (!sz.isEmpty())
    {
      QSize target(198, 120);
      QSize scaled = sz.scaled(target, Qt::KeepAspectRatioByExpanding);
      reader.setScaledSize(scaled);
      QImage img = reader.read();
      if (!img.isNull())
      {
        // Center-crop to exactly fill the thumbnail rect.
        int x = qMax(0, (img.width() - target.width()) / 2);
        int y = qMax(0, (img.height() - target.height()) / 2);
        pix = QPixmap::fromImage(img.copy(x, y, target.width(), target.height()));
      }
    }
  }

  if (!pix.isNull())
  {
    m_labelThumbnail->setPixmap(pix);
    return;
  }

  // Fallback: a bundled stage-appropriate icon.
  QString iconPath = ":/resource/upload@2x.png";
  if (m_info.currentStage == "connected_components")
    iconPath = ":/resource/cc@2x.png";
  else if (m_info.currentStage == "segmentation")
    iconPath = ":/resource/slabs@2x.png";
  else if (m_info.currentStage == "correction")
    iconPath = ":/resource/ruler@2x.png";
  QPixmap icon(iconPath);
  if (!icon.isNull())
    m_labelThumbnail->setPixmap(icon.scaled(56, 56, Qt::KeepAspectRatio, Qt::SmoothTransformation));
}

void ProjectCardWidget::mouseReleaseEvent(QMouseEvent* e)
{
  if (e->button() == Qt::LeftButton && rect().contains(e->pos()))
    emit Opened(m_info.rootOutputFolder);
  QFrame::mouseReleaseEvent(e);
}

void ProjectCardWidget::mouseDoubleClickEvent(QMouseEvent* e)
{
  emit Opened(m_info.rootOutputFolder);
  QFrame::mouseDoubleClickEvent(e);
}

void ProjectCardWidget::contextMenuEvent(QContextMenuEvent* e)
{
  // No parent: a child menu would inherit the Hub grid's unscoped "background:
  // transparent" (DialogProjectsHub) and draw see-through. On its own it takes
  // the app stylesheet's QMenu colours for the current theme.
  QMenu menu;
  if (m_bMissing)
  {
    QAction* locate = menu.addAction("Locate Folder…");
    QAction* remove = menu.addAction("Remove from List");
    QAction* chosen = menu.exec(e->globalPos());
    if (chosen == locate)
      emit LocateRequested(m_info.rootOutputFolder);
    else if (chosen == remove)
      emit RemoveRequested(m_info.rootOutputFolder);
    return;
  }
  QAction* open = menu.addAction("Open");
  QAction* rename = menu.addAction("Rename…");
  QAction* reveal = menu.addAction("Reveal in Finder");
  menu.addSeparator();
  QAction* remove = menu.addAction("Remove from List");
  QAction* chosen = menu.exec(e->globalPos());
  if (chosen == open)
    emit Opened(m_info.rootOutputFolder);
  else if (chosen == rename)
    emit RenameRequested(m_info.rootOutputFolder);
  else if (chosen == reveal)
    emit RevealRequested(m_info.rootOutputFolder);
  else if (chosen == remove)
    emit RemoveRequested(m_info.rootOutputFolder);
}
