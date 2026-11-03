#include "ImageGridView.h"
#include "QThumbnailWidget.h"
#include "MyUtils.h"

#include <QStackedWidget>
#include <QListWidget>
#include <QListWidgetItem>
#include <QAbstractItemView>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QImage>
#include <QPainter>
#include <QApplication>
#include <QEventLoop>
#include <QScrollBar>
#include <QStyle>
#include <QEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QStyledItemDelegate>
#include <QPainterPath>

namespace
{
// Grid cells are cropped-to-fit squares (MakeThumbnailThread center-crops
// to whatever size it's given), similar to Lightroom's square grid mode.
const int kGridCellSize = 160;
const int kCarouselCellSize = 68;

// Smallest horizontal gap between cells; UpdateGridMetrics() spreads leftover
// width evenly, so the real gap is usually larger.
const int kGridMinGap = 8;
// Fixed row gap, independent of the width-dependent column gap.
const int kGridRowGap = 16;

// Back pill for the lightbox top bar. Themed here, not by
// RethemeGhostPillButtons(): it sits on a different surface, so its tones
// differ from the left-column Back buttons.
const char* kBackButtonCssLight =
    "QPushButton {background: transparent; padding: 4 12 4 12; font-size: 14px; border: 1px solid #b0b0b0; border-radius: 14; color:#555555}"
    "QPushButton:hover {background: #ececec; border: 1px solid #8a8a8a; color:#333333}"
    "QPushButton:pressed {background: #dcdcdc; border: 1px solid #777777; color:#333333}";

// Dark-mode version with inverted tones for the dark top bar.
const char* kBackButtonCssDark =
    "QPushButton {background: transparent; padding: 4 12 4 12; font-size: 14px; border: 1px solid #5c5c5c; border-radius: 14; color:#dddddd}"
    "QPushButton:hover {background: #464646; border: 1px solid #8a8a8a; color:#ffffff}"
    "QPushButton:pressed {background: #565656; border: 1px solid #a0a0a0; color:#ffffff}";

// Space between the photo and its cell's item rect: the 3px item border plus
// Qt's 3px icon-mode margin (QThumbnailWidget::ItemPadding() is both sides).
const int kIconInset = 6;
const int kBadgeSize = 24;
const int kBadgeInset = 6;
// The carousel's thumbnails are 68px, so its badge is smaller.
const int kSmallBadgeSize = 18;
const int kSmallBadgeInset = 4;
// Teal, the app's neutral accent, so a picked photo does not read as the
// green "already corrected" tick the Correction strip uses.
const char* kBadgeOn = "#457c92";
}

// Paints the pick state over each cell: a tick badge top right and a veil over
// photos left out. Clicks are handled in eventFilter().
class GridSelectionDelegate : public QStyledItemDelegate
{
public:
  GridSelectionDelegate(ImageGridView* pView, QObject* parent, bool bSmall = false)
    : QStyledItemDelegate(parent), m_pView(pView), m_bSmall(bSmall) {}

  void paint(QPainter* p, const QStyleOptionViewItem& option, const QModelIndex& index) const override
  {
    QStyledItemDelegate::paint(p, option, index);
    if (!m_pView->IsSelectable())
      return;
    const bool bOn = m_pView->IsSelected(index.row());
    p->save();
    p->setRenderHint(QPainter::Antialiasing, true);
    if (!bOn)
    {
      QColor veil = m_pView->m_bDark ? QColor(16, 16, 16, 170) : QColor(224, 224, 224, 170);
      p->fillRect(option.rect.adjusted(kIconInset, kIconInset, -kIconInset, -kIconInset), veil);
    }
    QRectF rc = ImageGridView::BadgeRect(option.rect, m_bSmall);
    p->setPen(QPen(Qt::white, m_bSmall ? 1.5 : 2));
    p->setBrush(bOn ? QColor(kBadgeOn) : QColor(0, 0, 0, 90));
    p->drawEllipse(rc.adjusted(1, 1, -1, -1));
    if (bOn)
    {
      QPainterPath tick;
      tick.moveTo(rc.left() + rc.width()*0.28, rc.top() + rc.height()*0.52);
      tick.lineTo(rc.left() + rc.width()*0.44, rc.top() + rc.height()*0.68);
      tick.lineTo(rc.left() + rc.width()*0.73, rc.top() + rc.height()*0.36);
      QPen pen(Qt::white, m_bSmall ? 1.8 : 2.4);
      pen.setCapStyle(Qt::RoundCap);
      pen.setJoinStyle(Qt::RoundJoin);
      p->setPen(pen);
      p->setBrush(Qt::NoBrush);
      p->drawPath(tick);
    }
    p->restore();
  }

private:
  ImageGridView* m_pView;
  bool m_bSmall;
};

ImageGridView::ImageGridView(QWidget* parent)
  : QWidget(parent)
  , m_nCurrentIndex(-1)
  , m_bUpdatingGridMetrics(false)
  , m_bSelectable(false)
  , m_bDark(false)
{
  QVBoxLayout* pOuterLayout = new QVBoxLayout(this);
  pOuterLayout->setContentsMargins(0, 0, 0, 0);
  pOuterLayout->setSpacing(0);

  m_pStack = new QStackedWidget(this);
  // style.qss's app-wide QWidget background also matches this stack (a type
  // selector matches subclasses) and would paint over the panel colour from
  // paintEvent(). A rule in the widget's own stylesheet wins.
  m_pStack->setStyleSheet("QStackedWidget {background: transparent;}");
  pOuterLayout->addWidget(m_pStack);

  // Grid plus the selection bar above it. Transparent for the same reason as
  // m_pStack.
  m_pGridPage = new QWidget(m_pStack);
  m_pGridPage->setObjectName("gridPage");
  QVBoxLayout* pGridPageLayout = new QVBoxLayout(m_pGridPage);
  pGridPageLayout->setContentsMargins(0, 0, 0, 0);
  pGridPageLayout->setSpacing(0);

  m_pSelectionBar = new QWidget(m_pGridPage);
  m_pSelectionBar->setObjectName("gridSelectionBar");
  QHBoxLayout* pSelLayout = new QHBoxLayout(m_pSelectionBar);
  pSelLayout->setContentsMargins(16, 8, 10, 8);
  pSelLayout->setSpacing(8);
  m_pSelectionLabel = new QLabel(m_pSelectionBar);
  m_pSelectionLabel->setWordWrap(true);
  pSelLayout->addWidget(m_pSelectionLabel, 1);
  m_pSelectAll = new QPushButton("Select all", m_pSelectionBar);
  m_pSelectNone = new QPushButton("Select none", m_pSelectionBar);
  foreach (QPushButton* pBtn, QList<QPushButton*>() << m_pSelectAll << m_pSelectNone)
  {
    pBtn->setCursor(Qt::PointingHandCursor);
    pBtn->setFocusPolicy(Qt::NoFocus);
    pSelLayout->addWidget(pBtn);
  }
  connect(m_pSelectAll, &QPushButton::clicked, this, [this]() { SetAllSelected(true); });
  connect(m_pSelectNone, &QPushButton::clicked, this, [this]() { SetAllSelected(false); });
  pGridPageLayout->addWidget(m_pSelectionBar);
  m_pSelectionBar->hide();

  m_pGrid = new QThumbnailWidget(m_pGridPage);
  m_pGrid->setViewMode(QListWidget::IconMode);
  m_pGrid->setFlow(QListWidget::LeftToRight);
  m_pGrid->setWrapping(true);
  m_pGrid->setResizeMode(QListWidget::Adjust);   // reflows on resize -> responsive grid
  m_pGrid->setMovement(QListWidget::Static);
  // Fallback before the first layout; once UpdateGridMetrics() sets a grid size,
  // QListView ignores spacing.
  m_pGrid->setSpacing(8);
  m_pGrid->setIconSize(QSize(kGridCellSize, kGridCellSize));
  m_pGrid->setUniformItemSizes(true);
  m_pGrid->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  m_pGrid->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
  m_pGrid->setSelectionMode(QAbstractItemView::SingleSelection);
  m_pGrid->setFrameShape(QFrame::NoFrame);
  m_pGrid->setMouseTracking(true);
  // Don't let the viewport fill its own background over paintEvent()'s panel
  // colour.
  m_pGrid->viewport()->setAutoFillBackground(false);
  // Recompute cell size on the viewport's resize: its width is final by then,
  // unlike in this widget's resizeEvent().
  m_pGrid->viewport()->installEventFilter(this);
  connect(m_pGrid, SIGNAL(itemDoubleClicked(QListWidgetItem*)), SLOT(OnGridItemDoubleClicked(QListWidgetItem*)));
  m_pGrid->setItemDelegate(new GridSelectionDelegate(this, m_pGrid));
  m_pGrid->installEventFilter(this); // Space picks / leaves out the current photo

  pGridPageLayout->addWidget(m_pGrid, 1);
  m_pStack->addWidget(m_pGridPage);

  // --- Lightbox page ---------------------------------------------------
  m_pLightboxPage = new QWidget(m_pStack);
  m_pLightboxPage->setObjectName("lightboxPage");
  QVBoxLayout* pLightboxLayout = new QVBoxLayout(m_pLightboxPage);
  // No margins or spacing: the top bar and carousel band run edge to edge; each
  // has its own padding.
  pLightboxLayout->setContentsMargins(0, 0, 0, 0);
  pLightboxLayout->setSpacing(0);

  m_pTopBar = new QWidget(m_pLightboxPage);
  m_pTopBar->setObjectName("lightboxTopBar");
  QHBoxLayout* pTopBarLayout = new QHBoxLayout(m_pTopBar);
  pTopBarLayout->setContentsMargins(10, 8, 10, 8);
  m_pBackButton = new QPushButton("← Back", m_pTopBar);
  m_pBackButton->setCursor(Qt::PointingHandCursor);
  connect(m_pBackButton, SIGNAL(clicked()), SLOT(OnBackClicked()));
  pTopBarLayout->addWidget(m_pBackButton);
  pTopBarLayout->addStretch(1);
  m_pLightboxCount = new QLabel(m_pTopBar);
  pTopBarLayout->addWidget(m_pLightboxCount);
  pTopBarLayout->addSpacing(8);
  m_pLightboxSelect = new QPushButton(m_pTopBar);
  m_pLightboxSelect->setCheckable(true);
  m_pLightboxSelect->setCursor(Qt::PointingHandCursor);
  m_pLightboxSelect->setFocusPolicy(Qt::NoFocus);
  connect(m_pLightboxSelect, &QPushButton::clicked, this, [this]() {
    if (m_nCurrentIndex >= 0)
      SetSelected(m_nCurrentIndex, !IsSelected(m_nCurrentIndex));
  });
  pTopBarLayout->addWidget(m_pLightboxSelect);
  m_pLightboxCount->hide();
  m_pLightboxSelect->hide();
  pLightboxLayout->addWidget(m_pTopBar);

  // The image's 10px padding belongs to the dark image area, not the page, so no
  // dark gutter shows beside the bands.
  m_pImageArea = new QWidget(m_pLightboxPage);
  m_pImageArea->setObjectName("lightboxImageArea");
  QVBoxLayout* pImageAreaLayout = new QVBoxLayout(m_pImageArea);
  pImageAreaLayout->setContentsMargins(10, 10, 10, 10);
  pImageAreaLayout->setSpacing(0);

  m_pImageLabel = new QLabel(m_pImageArea);
  m_pImageLabel->setAlignment(Qt::AlignCenter);
  // Ignored size policy so the label doesn't hold the panel at the last pixmap's
  // size.
  m_pImageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
  // No extra bottom margin: m_pImageArea's padding separates it from the strip.
  m_pImageLabel->setMinimumSize(1, 1);
  pImageAreaLayout->addWidget(m_pImageLabel, 1);
  pLightboxLayout->addWidget(m_pImageArea, 1);

  m_pCarouselBand = new QWidget(m_pLightboxPage);
  m_pCarouselBand->setObjectName("lightboxCarouselBand");
  QVBoxLayout* pCarouselBandLayout = new QVBoxLayout(m_pCarouselBand);
  pCarouselBandLayout->setContentsMargins(10, 6, 10, 6);
  pCarouselBandLayout->setSpacing(0);

  m_pCarousel = new QThumbnailWidget(m_pCarouselBand);
  // No frame -- matches the grid's own m_pGrid->setFrameShape(NoFrame), and
  // one less thing eating into the fixed height below.
  m_pCarousel->setFrameShape(QFrame::NoFrame);
  // A single horizontally scrolling row: IconMode wraps by default, and the
  // vertical scrollbar is off, so wrapped thumbnails would be unreachable.
  m_pCarousel->setFlow(QListWidget::LeftToRight);
  m_pCarousel->setWrapping(false);
  // 80px item (68 + 12 padding) + 12px horizontal scrollbar, rounded up to 94.
  m_pCarousel->setFixedHeight(94);
  m_pCarousel->setIconSize(QSize(kCarouselCellSize, kCarouselCellSize));
  m_pCarousel->viewport()->setAutoFillBackground(false); // see the identical call on m_pGrid above
  // The same pick badges as the grid, smaller, and clickable here too.
  m_pCarousel->setItemDelegate(new GridSelectionDelegate(this, m_pCarousel, true));
  m_pCarousel->setMouseTracking(true);
  m_pCarousel->viewport()->installEventFilter(this);
  m_pCarousel->installEventFilter(this);
  connect(m_pCarousel, SIGNAL(CurrentIndexChanged(int)), SLOT(OnCarouselIndexChanged(int)));
  pCarouselBandLayout->addWidget(m_pCarousel);
  pLightboxLayout->addWidget(m_pCarouselBand);

  m_pStack->addWidget(m_pLightboxPage);

  m_pStack->setCurrentWidget(m_pGridPage);

  SetDarkMode(false); // establishes an initial (light) theme; MainWindow re-applies the real one right after construction
}

void ImageGridView::SetDarkMode(bool bDark)
{
  // Surface tones, deepest first (matching MainWindow::SetDarkMode()):
  //
  //   surface                              light      dark
  //   lightbox page / image area           #b0b0b0    #0a0a0a
  //   grid page background (paintEvent)    #c6c6c6    #101010   = #widgetMain
  //   grid cell background                 #e0e0e0    #242424
  //   lightbox top bar + carousel band     #fafafa    #383838   = #frameTopBar
  //
  // Accents (neutral, not the brand red):
  //
  //   accent                               light      dark
  //   hover outline                        #8a8a8a    #7d7d7d
  //   selected outline                     #4f4f4f    #d0d0d0
  //
  // "outline: none" hides Fusion's dotted focus rectangle.
  QString strHover    = bDark ? "#7d7d7d" : "#8a8a8a";
  QString strSelected = bDark ? "#d0d0d0" : "#4f4f4f";
  QString strPanelBg = bDark ? "#101010" : "#c6c6c6";
  // Cell tone, separated from near-white photos.
  QString strCellBg = bDark ? "#242424" : "#e0e0e0";

  // Painted in paintEvent(); stored here where the tone is computed.
  m_colorPanelBg = QColor(strPanelBg);
  update();

  m_pGrid->setStyleSheet(
        QString("QListWidget {background: transparent; border: none; outline: none;}"
                "QListWidget::item {background-color: %2; border: 3px solid transparent; border-radius: 4px; outline: none;}"
                "QListWidget::item:hover {border: 3px solid %1; outline: none;}"
                "QListWidget::item:selected {border: 3px solid %3; background: %2; outline: none;}")
        .arg(strHover, strCellBg, strSelected));

  // In the lightbox the image area is the darkest surface and the carousel band
  // lighter (the reverse of the grid page, where cells sit lighter on the panel).
  QString strLightboxBg = bDark ? "#0a0a0a" : "#b0b0b0";
  // The lightbox bands use #frameTopBar/#frameBottomBar's tones, so no new bar
  // colours are introduced.
  QString strBandBg = bDark ? "#383838" : "#fafafa";
  QString strCarouselBg = strBandBg; // no seam between the strip and its band
  // Object-name selectors: each band needs its own tone, and unnamed containers
  // would take style.qss's app-wide QWidget background.
  m_pLightboxPage->setStyleSheet(
        QString("#lightboxPage, #lightboxImageArea {background: %1;}"
                "#lightboxTopBar, #lightboxCarouselBand {background: %2;}")
        .arg(strLightboxBg, strBandBg));
  m_pBackButton->setStyleSheet(bDark ? kBackButtonCssDark : kBackButtonCssLight);
  m_bDark = bDark;
  m_pGridPage->setStyleSheet(
        QString("#gridPage {background: transparent;}"
                "#gridSelectionBar {background: %1;}"
                "#gridSelectionBar QLabel {background: transparent; font-size: 13px; color: %2;}")
        .arg(strBandBg, bDark ? "#dddddd" : "#333333"));
  const QString strPickCss = QString(bDark ? kBackButtonCssDark : kBackButtonCssLight) +
      (bDark ? "QPushButton:disabled {border: 1px solid #444444; color:#6a6a6a}"
             : "QPushButton:disabled {border: 1px solid #d6d6d6; color:#b0b0b0}");
  m_pSelectAll->setStyleSheet(strPickCss);
  m_pSelectNone->setStyleSheet(strPickCss);
  m_pLightboxCount->setStyleSheet(QString("QLabel {background: transparent; font-size: 13px; color: %1;}")
                                  .arg(bDark ? "#dddddd" : "#333333"));
  // Selected: the teal badge colour, filled. Left out: the ghost pill.
  m_pLightboxSelect->setStyleSheet(QString(bDark ? kBackButtonCssDark : kBackButtonCssLight) +
      "QPushButton:checked {background: #457c92; border: 1px solid #457c92; color: white; font-weight: bold;}"
      "QPushButton:checked:hover {background: #3c6d80; border: 1px solid #3c6d80;}"
      "QPushButton:checked:pressed {background: #345e6e; border: 1px solid #345e6e;}");
  m_pGrid->viewport()->update();
  // The image border has its own quiet tone, distinct from the selection
  // outline.
  QString strImageEdge = bDark ? "#303030" : "#9c9c9c";
  m_pImageLabel->setStyleSheet(QString("QLabel {background: %1; border: 1px solid %2;}").arg(strLightboxBg, strImageEdge));
  m_pCarousel->setStyleSheet(
        QString("QListWidget {background: %2; outline: none;}"
                "QListWidget::item {background-color: transparent; border: 3px solid transparent; outline: none;}"
                "QListWidget::item:hover {border: 3px solid %1; outline: none;}"
                "QListWidget::item:selected {background-color: transparent; border: 3px solid %3; outline: none;}")
        .arg(strHover, strCarouselBg, strSelected));
}

void ImageGridView::SetImages(const QStringList& files)
{
  m_listFiles = files;
  m_vSelected = QVector<bool>(files.size(), true);
  m_nCurrentIndex = -1;
  m_fullPixmap = QPixmap();
  UpdateSelectionBar();
  emit selectionChanged(SelectedCount(), m_listFiles.size());
  m_pGrid->LoadImageFiles(files);
  m_pCarousel->LoadImageFiles(files);
  m_pStack->setCurrentWidget(m_pGridPage);
}

void ImageGridView::Clear()
{
  m_listFiles.clear();
  m_vSelected.clear();
  UpdateSelectionBar();
  m_nCurrentIndex = -1;
  m_fullPixmap = QPixmap();
  m_pGrid->Clear();
  m_pCarousel->Clear();
  m_pStack->setCurrentWidget(m_pGridPage);
}

void ImageGridView::OnGridItemDoubleClicked(QListWidgetItem* item)
{
  ShowLightbox(m_pGrid->row(item));
}

void ImageGridView::ShowLightbox(int n)
{
  if (n < 0 || n >= m_listFiles.size())
    return;

  m_nCurrentIndex = n;
  m_pCarousel->SetCurrentIndex(n);
  m_pStack->setCurrentWidget(m_pLightboxPage);
  UpdateSelectionBar();
  // The label has not been laid out on first show (setCurrentWidget() only queues
  // it), so process pending events before measuring, or the image is scaled to a
  // stale size.
  QApplication::processEvents(QEventLoop::ExcludeUserInputEvents | QEventLoop::ExcludeSocketNotifiers);
  UpdateLightboxPixmap();
}

void ImageGridView::OnCarouselIndexChanged(int n)
{
  // QThumbnailWidget selects row 0 when its first thumbnail loads; ignore that
  // while the grid is showing, to avoid a full-resolution decode.
  if (m_pStack->currentWidget() != m_pLightboxPage)
    return;
  if (n < 0 || n >= m_listFiles.size() || n == m_nCurrentIndex)
    return;

  m_nCurrentIndex = n;
  UpdateSelectionBar();
  UpdateLightboxPixmap();
}

void ImageGridView::UpdateLightboxPixmap()
{
  if (m_nCurrentIndex < 0 || m_nCurrentIndex >= m_listFiles.size())
    return;

  QImage img = MyUtils::ReadImageWithExifAwareness(m_listFiles[m_nCurrentIndex]);
  m_fullPixmap = QPixmap::fromImage(img);
  QSize szTarget = m_pImageLabel->size();
  if (!m_fullPixmap.isNull() && szTarget.width() > 0 && szTarget.height() > 0)
    m_pImageLabel->setPixmap(m_fullPixmap.scaled(szTarget, Qt::KeepAspectRatio, Qt::SmoothTransformation));
  else
    m_pImageLabel->setPixmap(m_fullPixmap);
}

void ImageGridView::OnBackClicked()
{
  m_pStack->setCurrentWidget(m_pGridPage);
}

void ImageGridView::resizeEvent(QResizeEvent* e)
{
  QWidget::resizeEvent(e);
  if (m_pStack->currentWidget() == m_pLightboxPage && !m_fullPixmap.isNull())
  {
    QSize szTarget = m_pImageLabel->size();
    if (szTarget.width() > 0 && szTarget.height() > 0)
      m_pImageLabel->setPixmap(m_fullPixmap.scaled(szTarget, Qt::KeepAspectRatio, Qt::SmoothTransformation));
  }
}

void ImageGridView::SetSelectable(bool bSelectable)
{
  if (m_bSelectable == bSelectable)
    return;
  m_bSelectable = bSelectable;
  UpdateSelectionBar();
  if (!bSelectable)
  {
    m_pGrid->viewport()->unsetCursor();
    m_pCarousel->viewport()->unsetCursor();
  }
  emit selectionChanged(SelectedCount(), m_listFiles.size());
}

bool ImageGridView::IsSelected(int n) const
{
  if (!m_bSelectable)
    return true;
  return n >= 0 && n < m_vSelected.size() && m_vSelected[n];
}

void ImageGridView::SetSelected(int n, bool bSelected)
{
  if (n < 0 || n >= m_vSelected.size() || m_vSelected[n] == bSelected)
    return;
  m_vSelected[n] = bSelected;
  UpdateSelectionBar();
  emit selectionChanged(SelectedCount(), m_listFiles.size());
}

void ImageGridView::SetAllSelected(bool bSelected)
{
  m_vSelected.fill(bSelected);
  UpdateSelectionBar();
  emit selectionChanged(SelectedCount(), m_listFiles.size());
}

QStringList ImageGridView::SelectedFiles() const
{
  QStringList list;
  for (int i = 0; i < m_listFiles.size(); i++)
  {
    if (IsSelected(i))
      list << m_listFiles[i];
  }
  return list;
}

int ImageGridView::SelectedCount() const
{
  if (!m_bSelectable)
    return m_listFiles.size();
  return m_vSelected.count(true);
}

QRect ImageGridView::BadgeRect(const QRect& rcItem, bool bSmall)
{
  const int nSize = bSmall ? kSmallBadgeSize : kBadgeSize;
  const int nInset = bSmall ? kSmallBadgeInset : kBadgeInset;
  const QRect rcIcon = rcItem.adjusted(kIconInset, kIconInset, -kIconInset, -kIconInset);
  return QRect(rcIcon.right() + 1 - nInset - nSize, rcIcon.top() + nInset, nSize, nSize);
}

void ImageGridView::UpdateSelectionBar()
{
  // The bar is there whenever there are photos: it always says how to open
  // the carousel, and while picking it also counts the selection.
  m_pSelectionBar->setVisible(!m_listFiles.isEmpty());
  const int nTotal = m_listFiles.size();
  const int nOn = SelectedCount();
  const QString strHint = "Double-click a photo to open the carousel view.";
  if (m_bSelectable)
    m_pSelectionLabel->setText(QString("%1<br/><b>%2 of %3</b> photographs selected. Click a tick to leave a photo out.")
                               .arg(strHint).arg(nOn).arg(nTotal));
  else
    m_pSelectionLabel->setText(strHint);
  m_pSelectAll->setVisible(m_bSelectable);
  m_pSelectNone->setVisible(m_bSelectable);
  m_pSelectAll->setEnabled(nOn < nTotal);
  m_pSelectNone->setEnabled(nOn > 0);

  const bool bPick = m_bSelectable && m_nCurrentIndex >= 0 && m_nCurrentIndex < nTotal;
  m_pLightboxCount->setVisible(bPick);
  m_pLightboxSelect->setVisible(bPick);
  if (bPick)
  {
    const bool bOn = IsSelected(m_nCurrentIndex);
    m_pLightboxCount->setText(QString("%1 of %2 selected").arg(nOn).arg(nTotal));
    m_pLightboxSelect->setChecked(bOn);
    m_pLightboxSelect->setText(bOn ? QString::fromUtf8("\xe2\x9c\x93 Selected") : QString("Left out"));
    m_pLightboxSelect->setToolTip(bOn ? "Click to leave this photo out" : "Click to include this photo");
  }
  m_pGrid->viewport()->update();
  m_pCarousel->viewport()->update();
}

bool ImageGridView::HandleBadgeMouse(QListWidget* pView, QEvent* e, bool bSmall)
{
  // Badge clicks toggle the photo and are consumed, so they don't move the
  // current cell or open the lightbox. A little slack makes the badge easy to hit.
  const QEvent::Type t = e->type();
  if (t != QEvent::MouseButtonPress && t != QEvent::MouseButtonDblClick &&
      t != QEvent::MouseButtonRelease && t != QEvent::MouseMove)
    return false;
  QMouseEvent* me = static_cast<QMouseEvent*>(e);
  const QModelIndex idx = pView->indexAt(me->pos());
  const int nSlack = bSmall ? 3 : 4;
  const bool bOnBadge = idx.isValid() &&
      BadgeRect(pView->visualRect(idx), bSmall).adjusted(-nSlack, -nSlack, nSlack, nSlack).contains(me->pos());
  if (t == QEvent::MouseMove)
  {
    if (bOnBadge)
      pView->viewport()->setCursor(Qt::PointingHandCursor);
    else
      pView->viewport()->unsetCursor();
    return false;
  }
  if (bOnBadge && me->button() == Qt::LeftButton)
  {
    if (t != QEvent::MouseButtonRelease)
      SetSelected(idx.row(), !IsSelected(idx.row()));
    return true;
  }
  return false;
}

bool ImageGridView::eventFilter(QObject* obj, QEvent* e)
{
  if (obj == m_pGrid->viewport() && e->type() == QEvent::Resize)
    UpdateGridMetrics();

  if (m_bSelectable && obj == m_pGrid->viewport())
  {
    if (HandleBadgeMouse(m_pGrid, e, false))
      return true;
  }
  else if (m_bSelectable && obj == m_pCarousel->viewport())
  {
    if (HandleBadgeMouse(m_pCarousel, e, true))
      return true;
  }
  else if (m_bSelectable && (obj == m_pGrid || obj == m_pCarousel) && e->type() == QEvent::KeyPress)
  {
    // Space picks or leaves out the current photo, in either view.
    QListWidget* pView = (obj == m_pGrid) ? static_cast<QListWidget*>(m_pGrid) : static_cast<QListWidget*>(m_pCarousel);
    QKeyEvent* ke = static_cast<QKeyEvent*>(e);
    if (ke->key() == Qt::Key_Space && pView->currentRow() >= 0)
    {
      SetSelected(pView->currentRow(), !IsSelected(pView->currentRow()));
      return true;
    }
  }
  return QWidget::eventFilter(obj, e);
}

void ImageGridView::UpdateGridMetrics()
{
  // Spread leftover width evenly between columns instead of leaving a gutter on
  // the right. setGridSize() controls layout (spacing is ignored once a grid size
  // is set), and QListView centres each item in its cell, so a cell width of
  // viewport / columns gives equal outer margins.
  if (m_bUpdatingGridMetrics)
    return;

  const int nItem = kGridCellSize + QThumbnailWidget::ItemPadding();

  // Measure against the widget width minus a permanently reserved scrollbar, not
  // the live viewport width. The viewport lags a scrollbar appearing or
  // disappearing, which lost a column; reserving it keeps the result stable.
  const int nFrame = m_pGrid->frameWidth() * 2;
  int nScrollBar = m_pGrid->verticalScrollBar() ? m_pGrid->verticalScrollBar()->sizeHint().width() : 0;
  if (nScrollBar <= 0)
    nScrollBar = m_pGrid->style()->pixelMetric(QStyle::PM_ScrollBarExtent, nullptr, m_pGrid);
  const int nAvail = m_pGrid->width() - nFrame - nScrollBar;
  if (nAvail <= 0)
    return;

  int nCols = nAvail / (nItem + kGridMinGap);
  if (nCols < 1)
    nCols = 1;
  // nAvail - 1: QListView wraps when cells reach the right edge, so an exact fit
  // would lose a column.
  int nCellW = (nAvail - 1) / nCols;
  // Floor at nItem, not nItem + kGridMinGap: that floor undid the (nAvail - 1)
  // guard at exact multiples and lost a column. The clamp only applies when a
  // single item is wider than the space.
  if (nCellW < nItem)
    nCellW = nItem;

  const QSize szCell(nCellW, nItem + kGridRowGap);
  // Only update when the size changed: setGridSize() can toggle the scrollbar,
  // which resizes the viewport and re-enters here. The flag guards the same loop.
  if (m_pGrid->gridSize() == szCell)
    return;

  m_bUpdatingGridMetrics = true;
  m_pGrid->setGridSize(szCell);
  m_bUpdatingGridMetrics = false;
}

void ImageGridView::paintEvent(QPaintEvent* /*e*/)
{
  if (m_colorPanelBg.isValid())
  {
    QPainter p(this);
    p.fillRect(rect(), m_colorPanelBg);
  }
}
