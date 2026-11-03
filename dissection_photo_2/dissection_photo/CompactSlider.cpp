#include "CompactSlider.h"

#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <QFileInfo>
#include <QFontMetrics>
#include <algorithm>

namespace
{
// Unfilled track tone. A single #2d2d2d track vanished on the dark QWidget
// background (#2d2d2d), so dark mode uses a recessed #151515. That the light
// track equals the dark background is a coincidence: keep two constants.
const QColor kTrackLight(45, 45, 45);
const QColor kTrackDark(21, 21, 21);

// Drop shadow under swatch and track, painted by hand: a
// QGraphicsDropShadowEffect would also shadow the text, can't shadow the two
// shapes separately, and re-renders offscreen on every scrub. Falloff is
// kShadowPasses 1px-larger draws. One alpha for both themes: 8-bit rounding
// makes nearby alphas identical on the dark surface.
const int kShadowOffsetY  = 2;
const int kShadowPasses   = 3;
const int kShadowAlpha    = 11;   // per pass; composites to ~12% at the core

// Room below the shapes for the shadow tail. 8, not the 4 the shadow uses: the
// .ui pins these sliders to 34px, so less room would grow the shapes past 26px
// and out of step with the buttons around them.
const int kShadowRoom = 8;

// Minimum gap between label and value before the label elides.
const int kTextGap = 8;

// Glyph inset inside the 22px swatch; the 24px assets already have ~2px margin.
const int kIconInset = 1;
}

CompactSlider::CompactSlider(QWidget* parent)
  : QWidget(parent),
    m_nMin(0), m_nMax(100), m_nValue(0),
    m_colorSwatch(140, 140, 140),
    m_bDragging(false),
    m_bSwatchVisible(true),
    m_dIconDpr(0)
{
  setCursor(Qt::PointingHandCursor);
  setFocusPolicy(Qt::ClickFocus);
}

void CompactSlider::SetRange(int min, int max)
{
  m_nMin = min;
  m_nMax = std::max(min, max);
  m_nValue = qBound(m_nMin, m_nValue, m_nMax);
  update();
}

void CompactSlider::SetValue(int val)
{
  val = qBound(m_nMin, val, m_nMax);
  if (val == m_nValue)
    return;

  m_nValue = val;
  update();
}

void CompactSlider::SetLabel(const QString& label)
{
  m_sLabel = label;
  update();
}

void CompactSlider::SetSuffix(const QString& suffix)
{
  m_sSuffix = suffix;
  update();
}

void CompactSlider::SetSwatchVisible(bool bVisible)
{
  if (m_bSwatchVisible != bVisible)
  {
    m_bSwatchVisible = bVisible;
    update();
  }
}

void CompactSlider::SetSwatchColor(const QColor& color)
{
  m_colorSwatch = color;
  update();
}

void CompactSlider::SetSwatchIcon(const QString& strResourcePath)
{
  if (m_strIconPath == strResourcePath)
    return;

  m_strIconPath = strResourcePath;
  m_pixmapIcon = QPixmap();
  m_dIconDpr = 0;   // reload on the next paint
  update();
}

// The @2x path is chosen explicitly: Qt's automatic @2x lookup did not apply
// here on Qt 5 or Qt 6.
void CompactSlider::EnsureIconPixmap()
{
  qreal dDpr = devicePixelRatioF();
  if (m_strIconPath.isEmpty())
  {
    m_pixmapIcon = QPixmap();
    m_dIconDpr = dDpr;
    return;
  }
  // Keyed on DPR, so a missing file isn't retried on every repaint.
  if (qFuzzyCompare(m_dIconDpr, dDpr))
    return;

  m_dIconDpr = dDpr;
  m_pixmapIcon = QPixmap();

  if (dDpr > 1.5)
  {
    QFileInfo fi(m_strIconPath);
    QString str2x = fi.path() + "/" + fi.completeBaseName() + "@2x." + fi.suffix();
    if (m_pixmapIcon.load(str2x))
      m_pixmapIcon.setDevicePixelRatio(2.0);
  }
  if (m_pixmapIcon.isNull())
    m_pixmapIcon.load(m_strIconPath);
}

QSize CompactSlider::sizeHint() const
{
  return QSize(200, 26 + kShadowRoom);
}

QSize CompactSlider::minimumSizeHint() const
{
  return QSize(120, 24 + kShadowRoom);
}

QRect CompactSlider::ContentRect() const
{
  // Below 34px give up shadow room, never the shapes; floor at SWATCH_SIZE so
  // the text stays readable.
  int nRoom = std::min(kShadowRoom, std::max(0, height() - SWATCH_SIZE));
  return QRect(0, 0, width(), std::max(0, height() - nRoom));
}

QRect CompactSlider::TrackRect() const
{
  QRect rcContent = ContentRect();
  int x = m_bSwatchVisible ? SWATCH_SIZE + SWATCH_SPACING : 0;
  return QRect(x, rcContent.top(), std::max(0, width() - x), rcContent.height());
}

QRect CompactSlider::TrackHitRect() const
{
  QRect rcTrack = TrackRect();
  return QRect(rcTrack.left(), 0, rcTrack.width(), height());
}

QRect CompactSlider::SwatchRect() const
{
  QRect rcContent = ContentRect();
  return QRect(0, rcContent.top() + (rcContent.height() - SWATCH_SIZE) / 2,
               SWATCH_SIZE, SWATCH_SIZE);
}

void CompactSlider::PaintDropShadow(QPainter& p, const QRect& rc, int nRadius) const
{
  p.save();
  p.setPen(Qt::NoPen);
  p.setBrush(QColor(0, 0, 0, kShadowAlpha));
  for (int i = 0; i < kShadowPasses; i++)
  {
    QRect rcPass = rc.translated(0, kShadowOffsetY).adjusted(-i, -i, i, i);
    p.drawRoundedRect(rcPass, nRadius + i, nRadius + i);
  }
  p.restore();
}

void CompactSlider::PaintSwatch(QPainter& p, const QRect& rcSwatch, bool bDark)
{
  EnsureIconPixmap();

  if (m_pixmapIcon.isNull())
  {
    p.setPen(QColor(60, 60, 60));
    p.setBrush(m_colorSwatch);
    p.drawRoundedRect(rcSwatch, 4, 4);
    return;
  }

  // The glyphs are light grey, invisible on the light panel, so they sit on a
  // plate in the track's dark tone.
  p.setPen(Qt::NoPen);
  p.setBrush(bDark ? kTrackDark : kTrackLight);
  p.drawRoundedRect(rcSwatch, 4, 4);

  QRect rcGlyph = rcSwatch.adjusted(kIconInset, kIconInset, -kIconInset, -kIconInset);
  // scaled() takes device pixels: scale by the pixmap's DPR, then tag it so it
  // draws at logical size.
  qreal dDpr = m_pixmapIcon.devicePixelRatio();
  QPixmap pmDraw = m_pixmapIcon.scaled(QSize(qRound(rcGlyph.width() * dDpr),
                                             qRound(rcGlyph.height() * dDpr)),
                                       Qt::KeepAspectRatio, Qt::SmoothTransformation);
  pmDraw.setDevicePixelRatio(dDpr);
  QSize szLogical(qRound(pmDraw.width() / dDpr), qRound(pmDraw.height() / dDpr));
  QPoint ptTopLeft(rcGlyph.left() + (rcGlyph.width() - szLogical.width()) / 2,
                   rcGlyph.top() + (rcGlyph.height() - szLogical.height()) / 2);
  p.drawPixmap(ptTopLeft, pmDraw);
}

void CompactSlider::paintEvent(QPaintEvent* /*e*/)
{
  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing);

  // Read at paint time: MainWindow re-themes live.
  bool bDark = palette().color(QPalette::Window).lightness() < 128;

  QRect rcSwatch = SwatchRect();
  QRect rcTrack = TrackRect();
  if (rcTrack.width() <= 0)
    return;

  // Both shadows first, then both shapes, so the swatch shadow doesn't tint the
  // track's left edge.
  if (m_bSwatchVisible)
    PaintDropShadow(p, rcSwatch, 4);
  PaintDropShadow(p, rcTrack, 5);

  if (m_bSwatchVisible)
    PaintSwatch(p, rcSwatch, bDark);

  p.setPen(Qt::NoPen);
  p.setBrush(bDark ? kTrackDark : kTrackLight);
  p.drawRoundedRect(rcTrack, 5, 5);

  double frac = 0;
  if (m_nMax > m_nMin)
    frac = double(m_nValue - m_nMin) / double(m_nMax - m_nMin);
  frac = qBound(0.0, frac, 1.0);

  int fillWidth = int(rcTrack.width() * frac);
  if (fillWidth > 0)
  {
    QRect rcFill(rcTrack.left(), rcTrack.top(), fillWidth, rcTrack.height());
    QPainterPath path;
    path.addRoundedRect(rcTrack, 5, 5);
    p.save();
    p.setClipPath(path);
    p.setBrush(QColor(150, 150, 150));
    p.drawRect(rcFill);
    p.restore();
  }

  QString strValue = QString::number(m_nValue) + m_sSuffix;

  p.setPen(QColor(235, 235, 235));
  QFont f = p.font();
  // The stylesheet sets pixel size, so pointSizeF() is -1; scale whichever is set.
  if (f.pointSizeF() > 0)
    f.setPointSizeF(f.pointSizeF() * 0.9);
  else if (f.pixelSize() > 0)
    f.setPixelSize(qMax(1, qRound(f.pixelSize() * 0.9)));
  p.setFont(f);

  // 8px inset clears the track's 5px corner radius.
  QRect rcText = rcTrack.adjusted(8, 0, -8, 0);
  QFontMetrics fm(f);
  int nValueWidth = fm.horizontalAdvance(strValue);

  // If label and value collide, the label elides, never the number.
  if (!m_sLabel.isEmpty())
  {
    int nLabelRoom = rcText.width() - nValueWidth - kTextGap;
    // No room: drop the label; elidedText() at zero width returns an ellipsis
    // that would draw through the value.
    if (nLabelRoom > 0)
    {
      QString strLabel = m_sLabel;
      if (nLabelRoom < fm.horizontalAdvance(strLabel))
        strLabel = fm.elidedText(strLabel, Qt::ElideRight, nLabelRoom);
      p.drawText(rcText, Qt::AlignLeft | Qt::AlignVCenter, strLabel);
    }
  }
  p.drawText(rcText, Qt::AlignRight | Qt::AlignVCenter, strValue);
}

void CompactSlider::UpdateValueFromPos(int x)
{
  QRect rcTrack = TrackRect();
  if (rcTrack.width() <= 0)
    return;

  double frac = double(x - rcTrack.left()) / double(rcTrack.width());
  frac = qBound(0.0, frac, 1.0);

  int val = m_nMin + int(frac * (m_nMax - m_nMin) + 0.5);
  if (val != m_nValue)
  {
    m_nValue = val;
    update();
    emit valueChanged(m_nValue);
  }
}

void CompactSlider::mousePressEvent(QMouseEvent* e)
{
  // Full widget height, like the track, so shadow room never shrinks a target.
  QRect rcSwatchHit(0, 0, SWATCH_SIZE, height());
  if (m_bSwatchVisible && e->button() == Qt::LeftButton && rcSwatchHit.contains(e->pos()))
  {
    emit swatchClicked();
  }
  else if (e->button() == Qt::LeftButton && TrackHitRect().contains(e->pos()))
  {
    m_bDragging = true;
    UpdateValueFromPos(e->pos().x());
  }
  QWidget::mousePressEvent(e);
}

void CompactSlider::mouseMoveEvent(QMouseEvent* e)
{
  if (m_bDragging)
    UpdateValueFromPos(e->pos().x());
  QWidget::mouseMoveEvent(e);
}

void CompactSlider::mouseReleaseEvent(QMouseEvent* e)
{
  m_bDragging = false;
  QWidget::mouseReleaseEvent(e);
}
