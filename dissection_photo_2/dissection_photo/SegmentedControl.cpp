#include "SegmentedControl.h"

#include <QApplication>
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <QPropertyAnimation>
#include "AppMotion.h"
#include <QFontMetrics>
#include <QKeyEvent>
#include <QHelpEvent>
#include <QToolTip>

namespace
{
const char* kBrandRed = "#a2293c";
const int kMargin = 3;         // gap between track edge and indicator
const int kAnimDurationMs = 220;

// Fixed, not palette-derived: on #frameTopBar a palette track was invisible.
// Both are darker than the bar so the track reads as recessed.
const char* kTrackLight = "#e2e2e2";
const char* kTrackDark  = "#252525";

// Over QPalette::WindowText; keeps unselected labels above 4.5:1 on both tracks.
const qreal kUnselectedTextAlpha = 0.72;
}

SegmentedControl::SegmentedControl(QWidget* parent)
  : QWidget(parent)
  , m_nCurrentIndex(0)
  , m_dIndicatorPos(0.0)
  , m_pAnimation(new QPropertyAnimation(this, "indicatorPos", this))
  , m_colorAccent(kBrandRed)
  , m_colorTrackLight(kTrackLight)
  , m_colorTrackDark(kTrackDark)
{
  m_pAnimation->setDuration(AppMotion::Ms(kAnimDurationMs));
  m_pAnimation->setEasingCurve(AppMotion::Enter());
  setMinimumHeight(30);
  // Not StrongFocus, so a click leaves no focus ring.
  setFocusPolicy(Qt::TabFocus);
}

void SegmentedControl::SetItemEnabled(int index, bool bEnabled)
{
  if (index < 0 || index >= m_items.size())
    return;
  while (m_enabled.size() < m_items.size())
    m_enabled << true;
  if (m_enabled[index] == bEnabled)
    return;
  m_enabled[index] = bEnabled;
  update();
}

bool SegmentedControl::IsItemEnabled(int index) const
{
  if (index < 0 || index >= m_items.size())
    return false;
  return index >= m_enabled.size() || m_enabled[index];
}

void SegmentedControl::SetAccent(const QColor& color)
{
  m_colorAccent = color;
  update();
}

void SegmentedControl::SetTrackColors(const QColor& light, const QColor& dark)
{
  m_colorTrackLight = light;
  m_colorTrackDark = dark;
  update();
}

void SegmentedControl::SetItems(const QStringList& items)
{
  m_items = items;
  m_enabled.clear();
  m_pAnimation->stop();
  m_nCurrentIndex = 0;
  m_dIndicatorPos = 0.0;
  updateGeometry();
  update();
}

void SegmentedControl::SetCurrentIndex(int index, bool animate)
{
  if (m_items.isEmpty())
    return;

  index = qBound(0, index, m_items.size() - 1);
  // A disabled segment can never become current, whoever asks.
  if (!IsItemEnabled(index))
    return;
  bool bChanged = (index != m_nCurrentIndex);
  m_nCurrentIndex = index;

  m_pAnimation->stop();
  // Zero-length animations don't reliably land on the end value, so reduced
  // motion sets it directly.
  if (animate && AppMotion::Ms(kAnimDurationMs) > 0)
  {
    m_pAnimation->setStartValue(m_dIndicatorPos);
    m_pAnimation->setEndValue(static_cast<qreal>(index));
    m_pAnimation->start();
  }
  else
  {
    SetIndicatorPos(static_cast<qreal>(index));
  }

  if (bChanged)
    emit currentIndexChanged(index);
}

void SegmentedControl::SetIndicatorPos(qreal pos)
{
  m_dIndicatorPos = pos;
  update();
}

int SegmentedControl::SegmentMinWidth(const QFontMetrics& fmBold, const QString& str) const
{
  // Advance, not boundingRect(): drawText() lays out on advance widths. +1:
  // elidedText() works in qreal and otherwise elides at exactly this width.
  return fmBold.horizontalAdvance(str) + 2 * kMargin + 1;
}

QSize SegmentedControl::sizeHint() const
{
  // Measured bold: the selected segment is always drawn bold.
  QFont boldFont = font();
  boldFont.setBold(true);
  QFontMetrics fmBold(boldFont);

  int nWidth = 0;
  for (const QString& item : m_items)
    nWidth += fmBold.horizontalAdvance(item) + 32;
  if (nWidth < 1)
    nWidth = 200;
  // The minimum (widest item * N) can exceed this sum.
  return QSize(qMax(nWidth, minimumSizeHint().width()), 32);
}

QSize SegmentedControl::minimumSizeHint() const
{
  QFont boldFont = font();
  boldFont.setBold(true);
  QFontMetrics fmBold(boldFont);

  // Widest * N, not the sum: every segment gets total/N.
  int nMaxSegment = 0;
  for (const QString& item : m_items)
    nMaxSegment = qMax(nMaxSegment, SegmentMinWidth(fmBold, item));
  int nWidth = nMaxSegment * m_items.size();
  if (nWidth < 1)
    nWidth = 200;

  return QSize(nWidth, qMax(32, minimumHeight()));
}

qreal SegmentedControl::SegmentWidth() const
{
  if (m_items.isEmpty())
    return 0.0;
  return static_cast<qreal>(width()) / m_items.size();
}

QRectF SegmentedControl::SegmentTextRect(int index) const
{
  qreal dSegW = SegmentWidth();
  return QRectF(index * dSegW + kMargin, 0, qMax(0.0, dSegW - 2 * kMargin), height());
}

QRectF SegmentedControl::IndicatorRect() const
{
  qreal dSegW = SegmentWidth();
  return QRectF(m_dIndicatorPos * dSegW + kMargin, kMargin, qMax(0.0, dSegW - 2 * kMargin),
                height() - 2 * kMargin);
}

void SegmentedControl::paintEvent(QPaintEvent*)
{
  if (m_items.isEmpty())
    return;

  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing, true);

  QRectF rcTrack(0, 0, width(), height());
  qreal dTrackRadius = rcTrack.height() / 2.0;

  // The application palette, not ours: a "background: transparent" page sheet
  // makes this widget's own Window transparent, which reads as dark.
  bool bDark = QApplication::palette().color(QPalette::Window).lightness() < 128;
  QColor trackColor(bDark ? m_colorTrackDark : m_colorTrackLight);
  QColor unselectedTextColor = palette().color(QPalette::WindowText);
  unselectedTextColor.setAlphaF(kUnselectedTextAlpha);
  QColor disabledTextColor = palette().color(QPalette::WindowText);
  disabledTextColor.setAlphaF(0.35);
  const bool bWidgetEnabled = isEnabled();

  QPainterPath trackPath;
  trackPath.addRoundedRect(rcTrack, dTrackRadius, dTrackRadius);
  painter.fillPath(trackPath, trackColor);

  QPainterPath indicatorPath;
  QRectF rcIndicator = IndicatorRect();
  qreal dIndicatorRadius = rcIndicator.height() / 2.0;
  indicatorPath.addRoundedRect(rcIndicator, dIndicatorRadius, dIndicatorRadius);
  painter.fillPath(indicatorPath, bWidgetEnabled ? m_colorAccent : QColor(bDark ? "#5a5a5a" : "#b4b4b4"));

  QFont boldFont = font();
  boldFont.setBold(true);
  QFont normalFont = font();
  normalFont.setBold(false);

  for (int i = 0; i < m_items.size(); i++)
  {
    QRectF rcText = SegmentTextRect(i);
    // Keyed to the animated position, so the label changes as the pill crosses.
    bool bSelected = (qAbs(m_dIndicatorPos - i) < 0.5);
    painter.setFont(bSelected ? boldFont : normalFont);
    if (bSelected)
      painter.setPen(bWidgetEnabled ? QColor(Qt::white) : QColor(bDark ? "#c8c8c8" : "#f4f4f4"));
    else
      painter.setPen((bWidgetEnabled && IsItemEnabled(i)) ? unselectedTextColor : disabledTextColor);

    // The .ui's minimumSize overrides minimumSizeHint(), so elide to the pill.
    QFontMetrics fmDraw(bSelected ? boldFont : normalFont);
    QString strLabel = fmDraw.elidedText(m_items.at(i), Qt::ElideRight,
                                         static_cast<int>(rcText.width()));
    painter.drawText(rcText, Qt::AlignCenter, strLabel);
  }

  // Keyboard focus ring (Tab only).
  if (hasFocus())
  {
    QRectF rcRing = rcTrack.adjusted(0.5, 0.5, -0.5, -0.5);
    painter.setPen(QPen(palette().color(QPalette::Highlight), 1.0));
    painter.setBrush(Qt::NoBrush);
    painter.drawRoundedRect(rcRing, rcRing.height() / 2.0, rcRing.height() / 2.0);
  }
}

void SegmentedControl::mousePressEvent(QMouseEvent* e)
{
  if (e->button() != Qt::LeftButton)
  {
    QWidget::mousePressEvent(e);
    return;
  }

  // SegmentWidth() is 0 without items; SetCurrentIndex() ignores disabled ones.
  qreal dSegW = SegmentWidth();
  if (dSegW <= 0)
    return;
  SetCurrentIndex(static_cast<int>(e->pos().x() / dSegW), true);
}

void SegmentedControl::keyPressEvent(QKeyEvent* e)
{
  int nStep = 0;
  if (e->key() == Qt::Key_Left)
    nStep = -1;
  else if (e->key() == Qt::Key_Right)
    nStep = 1;
  if (nStep == 0 || m_items.isEmpty())
  {
    QWidget::keyPressEvent(e);
    return;
  }
  for (int n = m_nCurrentIndex + nStep; n >= 0 && n < m_items.size(); n += nStep)
  {
    if (IsItemEnabled(n))
    {
      SetCurrentIndex(n, true);
      break;
    }
  }
  e->accept();
}

bool SegmentedControl::event(QEvent* e)
{
  if (e->type() == QEvent::ToolTip && !m_items.isEmpty())
  {
    QHelpEvent* h = static_cast<QHelpEvent*>(e);
    qreal dSegW = SegmentWidth();
    int n = dSegW > 0 ? qBound(0, static_cast<int>(h->pos().x() / dSegW), m_items.size() - 1) : -1;
    if (n >= 0 && n < m_tips.size() && !m_tips.at(n).isEmpty())
    {
      QToolTip::showText(h->globalPos(), m_tips.at(n), this);
      return true;
    }
  }
  return QWidget::event(e);
}
