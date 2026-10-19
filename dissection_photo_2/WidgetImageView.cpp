#include "WidgetImageView.h"
#include <QPainter>
#include <QFile>
#include <QMouseEvent>
#include "math.h"
#include "MyUtils.h"
#include <QGuiApplication>
#include <QStyleHints>
#include <QPropertyAnimation>

// Default fit leaves a thin margin of background around the image.
static const double DEFAULT_FIT_SCALE = 0.92;

// Preview crossfade: AppMotion::Base (180 ms) with OutCubic, like the other
// animations.
static const int PREVIEW_FADE_DURATION_MS = 180;

static int s_nMotionReduced = -1;   // -1 unread, 0 off, 1 on
bool WidgetImageView::MotionReduced()
{
  if (s_nMotionReduced < 0)
  {
    const QByteArray env = qgetenv("DISSECTION_PHOTO_REDUCED_MOTION");
    s_nMotionReduced = (env == "1") ? 1 : (env == "0") ? 0
                     : (QGuiApplication::platformName() == QLatin1String("offscreen")) ? 1 : 0;
  }
  return s_nMotionReduced == 1;
}

WidgetImageView::WidgetImageView(QWidget *parent)
  : QWidget(parent), m_dScale(DEFAULT_FIT_SCALE), m_ptOffset(QPoint(0,0)), m_bPanning(false), m_bZooming(false), m_bDrawing(false),
    m_nNumberOfExpectedPoints(2), m_nEditMode(0), m_dMaskOpacity(0.7), m_dCombinedMaskOpacity(1.0),
    m_nBrushSize(10), m_bMaskEdited(false), m_bModifierDown(false)
{
  setMouseTracking(true);
  m_colorPen = QColor(50,255,50);

  if (parent)
    parent->installEventFilter(this);
}

bool WidgetImageView::eventFilter(QObject *watched, QEvent *event)
{
  if (watched == parentWidget() && event->type() == QEvent::KeyRelease
      && static_cast<QKeyEvent*>(event)->key() == Qt::Key_Control)
    unsetCursor();
  return QWidget::eventFilter(watched, event);
}

bool WidgetImageView::LoadImage(const QString& filename, const QString& mask, const QList<QPoint>& points, const QList<RECT_REGION>& regions)
{
  // Consume the one-shot arm from BeginPreviewCrossfade(). Every other load
  // cancels any fade, so Next mid-fade never blends the old photo over the new.
  const bool bCrossfade = m_bPreviewFadePending;
  m_bPreviewFadePending = false;
  if (!bCrossfade)
    CancelPreviewCrossfade();

  m_sMessage.clear();
  m_listMaskUndoBuffer.clear();
  m_listMaskRedoBuffer.clear();
  m_sFilename = filename;
  m_sMaskFilename = mask;
  m_imageOriginal = m_imageMask = QImage();
  // Callers save before loading; unsaved strokes go with the old mask.
  m_bMaskEdited = false;
  if (!mask.isEmpty() && !QFile::exists(mask))
  {
    // No mask yet: start from an empty one, written only once drawn on.
    m_imageOriginal = MyUtils::ReadImageWithExifAwareness(filename);
    m_imageMask = QImage(m_imageOriginal.size(), QImage::Format_ARGB32);
    m_imageMask.fill(m_nEditMode == EM_EDIT_MASK ? QColor(0,0,0,0) : QColor(0,0,0,255));
  }
  else if (!mask.isEmpty())
  {
    // The mask has no EXIF tag of its own; rotate it by the photo's orientation so
    // it stays aligned with the rotated photo.
    m_imageMask = MyUtils::ApplyExifOrientation(QImage(m_sMaskFilename), MyUtils::GetExifOrientation(filename))
                      .convertToFormat(QImage::Format_ARGB32);
    if (m_nEditMode == EM_EDIT_MASK)  // remove background
      SetAlphaByMask(m_imageMask);
  }
  m_imageOverlay = QImage();
  PrepareImage();
  if (!m_image.isNull())
  {
    m_dScale = DEFAULT_FIT_SCALE;
    m_ptOffset = QPoint(0,0);
    m_listPoints = points;
    m_listRedoPoints.clear();
    m_listRegions = regions;
    m_listRedoRegions.clear();
    UpdateScaledImage();
    // After UpdateScaledImage(), so the incoming frame the outgoing one is
    // fading out over is already the one paintEvent() will draw.
    if (bCrossfade)
      StartPreviewCrossfade();
    return true;
  }
  else
  {
    // Nothing loaded: release the retained frame.
    CancelPreviewCrossfade();
    return false;
  }
}

void WidgetImageView::SetAlphaByMask(QImage& image)
{
  for (int y = 0; y < image.height(); y++)
  {
    QRgb* p_dest = (QRgb*)image.scanLine(y);
    for (int x = 0; x < image.width(); x++)
    {
      if (qRed(p_dest[x]) == 0)
      {
        p_dest[x] = qRgba(0,0,0,0);
      }
    }
  }
}

void WidgetImageView::PrepareImage()
{
  // Nothing loaded yet: don't try to open an empty filename.
  if (m_imageOriginal.isNull() && !m_sFilename.isEmpty())
    m_imageOriginal = MyUtils::ReadImageWithExifAwareness(m_sFilename);
  if (!m_imageOriginal.isNull())
  {
    m_image = m_imageOriginal;
    QImage image = m_imageMask;
    if (!m_sMaskFilename.isEmpty())
    {
      QPainter p(&m_image);
      if (m_nEditMode == EM_EDIT_MASK)
      {
        p.setOpacity(m_dMaskOpacity);
        QPainter p2(&image);
        p2.setCompositionMode(QPainter::CompositionMode_SourceIn);
        p2.fillRect(image.rect(), QColor(0,255,0));
      }
      else
      {
        p.setCompositionMode(QPainter::CompositionMode_Multiply);
      }
      p.drawImage(0, 0, image);
      p.end();
    }
    if (!m_imageOverlay.isNull())
    {
      QPainter p(&m_image);
      p.setOpacity(m_dCombinedMaskOpacity);
      p.drawImage(0, 0, m_imageOverlay);
      p.end();
    }
  }
}

void WidgetImageView::SetEditMode(int n)
{
  m_nEditMode = n;
  if (n == EM_CALIBRATION)
    m_nNumberOfExpectedPoints = 8;
}

void WidgetImageView::SetOverlay(const QImage& overlay_image)
{
  m_imageOverlay = overlay_image;
  PrepareImage();
  UpdateScaledImage();
}

void WidgetImageView::paintEvent(QPaintEvent *e)
{
  QPainter p(this);
  p.fillRect(rect(), m_colorBackground);
  if (m_sFilename.isEmpty())
    return;

  QRect target = m_imageScaled.rect();
  target.moveCenter(rect().center()+m_ptOffset);
  p.drawImage(target, m_imageScaled);

  // Preview crossfade: the outgoing frame is drawn over the incoming image at the
  // complementary opacity, in its own rect (the corrected image can have a
  // different size and aspect ratio, so the overhang blends over background
  // during the fade). Opacity only; pan and zoom apply to both. Drawn before the
  // overlays so points, regions and messages stay crisp.
  if (!m_imageFadeFrom.isNull() && m_dPreviewFadeOpacity < 1.0)
  {
    QRect rcFade = m_imageFadeFrom.rect();
    rcFade.moveCenter(rect().center()+m_ptOffset);
    const qreal dPrevOpacity = p.opacity();
    p.setOpacity(1.0 - m_dPreviewFadeOpacity);
    p.drawImage(rcFade, m_imageFadeFrom);
    p.setOpacity(dPrevOpacity);
  }

  if (m_nEditMode == EM_POINT)
  {
    p.setPen(m_colorPen);
    p.setBrush(m_colorPen);
    QList<QPoint> pts;
    foreach (QPoint pt, m_listPoints)
    {
      pt = pt*m_imageScaled.width()/m_image.width() + target.topLeft();
      p.drawEllipse(pt, 2, 2);
      pts << pt;
    }
    if (pts.size() > 1 && m_nNumberOfExpectedPoints == 2)
    {
      p.setPen(QPen(m_colorPen, 2));
      p.drawLine(pts[0], pts[1]);
    }
  }
  else if (m_nEditMode == EM_REGION)
  {
    p.setPen(QPen(QColor(255,255,255), 2));
    p.setBrush(Qt::NoBrush);
    for (int i = 0; i < m_listRegions.size(); i++)
    {
      if (i == m_listRegions.size()-1 && m_bDrawing)
        p.setPen(QPen(QColor(255,30,30), 2));

      QRect rc(m_listRegions[i].first*m_imageScaled.width()/m_image.width() + target.topLeft(),
               m_listRegions[i].second*m_imageScaled.width()/m_image.width() + target.topLeft());
      p.drawRect(rc);
    }
  }
  else if (m_nEditMode == EM_CALIBRATION)
  {
    for (int i = 0; i < m_listPoints.size(); i++)
    {
      QPoint pt = m_listPoints[i]*m_imageScaled.width()/m_image.width() + target.topLeft();
      if (i%2 == 0)
      {
        p.setBrush(m_colorPen);
        p.setPen(m_colorPen);
        p.drawEllipse(pt, 2, 2);
      }
      else
      {
        QPoint prev_pt = m_listPoints[i-1]*m_imageScaled.width()/m_image.width() + target.topLeft();
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(m_colorPen,2));
        double r = sqrt((pt.x()-prev_pt.x())*(pt.x()-prev_pt.x())+(pt.y()-prev_pt.y())*(pt.y()-prev_pt.y()));
        p.drawEllipse(QPointF(prev_pt), r, r);
      }
    }
  }
  else if (m_nEditMode == EM_EDIT_MASK)
  {
    // Hover preview of what a click would do: a size circle for brush/eraser, a
    // marker for flood actions, darker while a modifier is held. Follows
    // m_nActiveAction, so flood fill and flood erase are told apart.
    QPointF pt = this->mapFromGlobal(QCursor::pos());
    int nAlpha = m_bModifierDown ? 140 : 80;
    // Radius of the circle drawn below; the tool glyph is placed off it.
    double dCursorRadius = 5.0;
    const bool bFlood = (m_nActiveAction == MA_FLOOD_FILL || m_nActiveAction == MA_FLOOD_ERASE);
    // Both erases take the warm tone (flood erase too, so the two floods are not
    // one indistinguishable marker); flood fill is green, the brush blue.
    QColor color = (m_nActiveAction == MA_ERASE || m_nActiveAction == MA_FLOOD_ERASE) ? QColor(180, 100, 80)
                 : bFlood ? QColor(120, 160, 90) : QColor(80, 140, 180);
    p.setBrush(QColor(color.red(), color.green(), color.blue(), nAlpha));
    p.setPen(QPen(color, 1));
    if (!bFlood)
      dCursorRadius = CurrentToolSize()/2.0*m_imageScaled.width()/m_image.width();
    p.drawEllipse(pt, dCursorRadius, dCursorRadius);
    DrawToolCursorIcon(p, pt, m_nActiveAction, dCursorRadius);
  }
  if (!m_sMessage.isEmpty())
  {
    QFont f = p.font();
    f.setPixelSize(16);
    p.setFont(f);
    QRect rc = p.fontMetrics().boundingRect(m_sMessage);
    rc.adjust(-20,-9, 20, 9);
    rc.moveCenter(rect().center());
    rc.moveBottom(rect().bottom()-18);
    p.setPen(Qt::NoPen);
    p.setBrush(QColor(0,0,0,200));
    p.drawRoundedRect(rc, rc.height()/2, rc.height()/2);
    p.setPen(Qt::white);
    p.drawText(rc, m_sMessage, Qt::AlignHCenter|Qt::AlignVCenter);
  }
}

// 18 logical px glyph plus 3 px of scrim per side: a 24 px badge, one toolbar icon.
static const double TOOL_CURSOR_ICON_SIZE = 18.0;
static const double TOOL_CURSOR_ICON_PAD = 3.0;

const QPixmap& WidgetImageView::ToolCursorIcon(int nAction)
{
  if (!m_bToolIconsLoaded)
  {
    m_bToolIconsLoaded = true;

    // Load the @2x file explicitly: Qt's automatic @2x lookup did not apply here,
    // and a 48px source downsamples cleanly at any DPR.
    // MA_FLOOD_ERASE uses the bucket glyph: colour means add vs remove, glyph means
    // flood vs stroke. With the eraser glyph it would look identical to stroke-erase.
    const char* aszNames[4] = { "brush", "eraser", "bucket", "bucket" };
    for (int i = 0; i < 4; i++)
    {
      QPixmap pix(QString(":/resource/icons/%1@2x.png").arg(aszNames[i]));
      if (pix.isNull())
        continue;   // resource absent -- leave this slot null, see below

      // An "@2x" filename makes Qt tag the pixmap DPR 2; untag it so the tint copy
      // below is 1:1. The draw rect sets the on-screen size.
      pix.setDevicePixelRatio(1.0);

      // Tint white (only the asset's alpha matters), like the carousel number badge.
      QPixmap out(pix.size());
      out.fill(Qt::transparent);
      QPainter pp(&out);
      pp.drawPixmap(0, 0, pix);
      pp.setCompositionMode(QPainter::CompositionMode_SourceIn);
      pp.fillRect(out.rect(), Qt::white);
      pp.end();
      m_pixToolCursorIcons[i] = out;
    }
  }

  // Out of range (MA_NONE should not get here): no badge.
  if (nAction < MA_DRAW || nAction > MA_FLOOD_ERASE)
    return m_pixNullToolIcon;
  return m_pixToolCursorIcons[nAction];
}

void WidgetImageView::DrawToolCursorIcon(QPainter& p, const QPointF& pt, int nAction, double dCircleRadius)
{
  const QPixmap& pix = ToolCursorIcon(nAction);
  if (pix.isNull())
    return;   // no icon resource (e.g. renamed): the circle only

  // Offset lower-right, off the hot spot. The radius ranges from 0.5 to 100+ px,
  // so the offset is clamped to 8..34 px: a tiny brush doesn't cover its hot spot,
  // and a huge one doesn't push the badge away from the cursor (it sits inside
  // the circle, where the scrim keeps it legible).
  const double dOffset = qBound(8.0, dCircleRadius*0.7071 + 4.0, 34.0);
  const double dBadge = TOOL_CURSOR_ICON_SIZE + 2*TOOL_CURSOR_ICON_PAD;
  const QRectF rcBadge(pt.x() + dOffset, pt.y() + dOffset, dBadge, dBadge);

  p.save();
  // Antialiasing for the rounded scrim; smooth transform to downsample the 48px
  // glyph.
  p.setRenderHint(QPainter::Antialiasing, true);
  p.setRenderHint(QPainter::SmoothPixmapTransform, true);

  // Scrim values from the carousel number badge: a dark translucent body with a
  // light hairline reads on light and dark photos. Fixed alphas, independent of
  // the modifier state.
  const double dPen = 1.0;
  p.setPen(Qt::NoPen);
  p.setBrush(QColor(0, 0, 0, 70));
  p.drawRoundedRect(rcBadge.translated(0, dPen), 6.0, 6.0);

  const QRectF rcBody = rcBadge.adjusted(dPen/2.0, dPen/2.0, -dPen/2.0, -dPen/2.0);
  p.setBrush(QColor(0, 0, 0, 170));
  p.setPen(QPen(QColor(255, 255, 255, 115), dPen));
  p.drawRoundedRect(rcBody, 6.0, 6.0);

  // Rect overload: fixes the glyph at 18 logical px whatever the pixmap's DPR.
  QRectF rcGlyph(0, 0, TOOL_CURSOR_ICON_SIZE, TOOL_CURSOR_ICON_SIZE);
  rcGlyph.moveCenter(rcBadge.center());
  p.setBrush(Qt::NoBrush);
  p.drawPixmap(rcGlyph, pix, QRectF(pix.rect()));
  p.restore();
}

int WidgetImageView::MaskActionFor(Qt::MouseButton button, Qt::KeyboardModifiers mods) const
{
  const bool bCtrl = mods & Qt::ControlModifier;
  const bool bShift = mods & Qt::ShiftModifier;

  // What the toolbar selection alone means, used both for Ctrl+Left (the
  // "apply the selected tool" binding) and for the idle hover preview.
  int nSelected = MA_DRAW;
  if (m_nMaskTool == MT_ERASER)
    nSelected = MA_ERASE;
  else if (m_nMaskTool == MT_BUCKET)
    nSelected = MA_FLOOD_FILL;

  if (button == Qt::RightButton)
  {
    // Plain right-drag is zoom in every mode; only the Ctrl bindings are
    // mask operations here.
    if (bCtrl && bShift)
      return MA_FLOOD_ERASE;
    if (bCtrl)
      return MA_FLOOD_FILL;
    return MA_NONE;
  }

  // Q / W / E pick the tool and Ctrl + Left applies it (Command on macOS, where
  // Qt maps ControlModifier to Command).
  if (bCtrl)
    return nSelected;

  // No modifier (or Shift): a click pans; the hover preview (NoButton) shows what
  // Ctrl would do.
  return (button == Qt::NoButton) ? nSelected : MA_NONE;
}

void WidgetImageView::RefreshMaskAction(Qt::KeyboardModifiers mods)
{
  if (m_nEditMode != EM_EDIT_MASK)
    return;
  // Never mid-stroke: a stroke keeps the action (and size) it started with.
  if (m_bDrawing)
    return;
  int n = MaskActionFor(Qt::NoButton, mods);
  if (n != MA_NONE)
    m_nActiveAction = n;
}

bool WidgetImageView::BeginMaskAction(QMouseEvent* e)
{
  int nAction = MaskActionFor(e->button(), e->modifiers());
  if (nAction == MA_NONE)
    return false;

  m_nActiveAction = nAction;
  // m_nMaskTool is not changed: flood overrides are momentary.
  m_bErasing = (nAction == MA_ERASE);
  m_bModifierDown = true;
  m_listMaskUndoBuffer << m_imageMask;
  m_listMaskRedoBuffer.clear(); // a new edit invalidates whatever was available to redo

  if (nAction == MA_FLOOD_FILL || nAction == MA_FLOOD_ERASE)
  {
    // A flood is one click: m_bDrawing stays false so a following drag doesn't
    // paint.
    m_bDrawing = false;
    FloodFillMaskImage(m_ptPress, nAction == MA_FLOOD_ERASE);
  }
  else
  {
    m_bDrawing = true;
    FreeHandOnMaskImage(m_ptPrev, m_ptPrev);
  }
  return true;
}

void WidgetImageView::mousePressEvent(QMouseEvent *e)
{
  m_ptPress = e->pos();
  m_ptPrev = e->pos();
  m_ptOldOffset = m_ptOffset;
  m_dOldScale = m_dScale;

  // Mask editing has its own bindings (see MaskAction).
  if (m_nEditMode == EM_EDIT_MASK && BeginMaskAction(e))
    return;
  // Any other left press on Segmentation (plain or Shift) pans, instead of
  // reaching the region-drawing code below.
  if (m_nEditMode == EM_EDIT_MASK && e->button() == Qt::LeftButton)
  {
    m_bPanning = true;
    return;
  }

  if (e->button() == Qt::LeftButton)
  {
    if (e->modifiers() & Qt::ControlModifier)
    {
      m_bDrawing = true;
      QPoint pt = ScreenToImage(m_ptPress);
      m_listRegions << RECT_REGION(pt, pt);
      m_listRedoRegions.clear(); // a new region invalidates whatever was available to redo
      setCursor(Qt::CrossCursor);
    }
    else if (!(e->modifiers() & Qt::ShiftModifier))
      m_bPanning = true;
  }
  else if (e->button() == Qt::RightButton)
  {
    if (e->modifiers() & Qt::ControlModifier)
      m_bDrawing = true;
    else
      m_bZooming = true;
  }
}

QPoint WidgetImageView::ScreenToImage(const QPoint &pt_in)
{
  QRect target = m_imageScaled.rect();
  target.moveCenter(rect().center()+m_ptOffset);
  return (pt_in - target.topLeft())*m_image.width()/m_imageScaled.width();
}

void WidgetImageView::mouseReleaseEvent(QMouseEvent *e)
{
  // A click is anything within the platform's drag distance (normally 10 px; the
  // old 2 px dropped clicks made while the pointer was still moving). The point
  // goes where it was pressed.
  const bool bClick = (e->pos() - m_ptPress).manhattanLength() < QGuiApplication::styleHints()->startDragDistance();
  if (m_bDrawing && !m_imageScaled.isNull())
  {
    if ((m_nEditMode == EM_POINT || m_nEditMode == EM_CALIBRATION) && !m_bPointsLocked && bClick)
    {
      QPoint pt = ScreenToImage(m_ptPress);
      if (m_listPoints.size() < m_nNumberOfExpectedPoints)
        m_listPoints << pt;
      else
        m_listPoints[m_nNumberOfExpectedPoints-1] = pt;
      m_listRedoPoints.clear();
      emit PointEdited();
      update();
    }
    else if ((m_nEditMode == EM_POINT || m_nEditMode == EM_CALIBRATION) && m_bPointsLocked
             && bClick)
      emit PointRejected();
    // Ctrl+Right also sets m_bDrawing. Test the released button (left) so a
    // Ctrl+Right drag cannot touch an earlier region, and isEmpty() before last().
    else if (m_nEditMode == EM_REGION && e->button() == Qt::LeftButton
             && !m_listRegions.isEmpty())
    {
      update();
      QPoint pt = m_listRegions.last().first - m_listRegions.last().second;
      if (qAbs(pt.x()) < 3 || qAbs(pt.y()) < 3)
        m_listRegions.removeLast();
      else
        emit LastRegionEdited(m_listRegions.size()-1);
    }
  }
  else if (m_bZooming)
  {
    UpdateScaledImage(true);
  }
  m_bPanning = false;
  m_bZooming = false;
  m_bDrawing = false;
  m_bErasing = false;
  // The stroke is over, so any momentary Ctrl override is over with it --
  // fall back to whatever the toolbar still says is selected.
  RefreshMaskAction(e->modifiers());
  unsetCursor();
  update();
}

void WidgetImageView::mouseMoveEvent(QMouseEvent *e)
{
  QPoint dpt = e->pos() - m_ptPress;
  m_bModifierDown = false;
  if (m_bPanning)
  {
    m_ptOffset = m_ptOldOffset + dpt;
    update();
  }
  else if (m_bZooming)
  {
    int dn = e->y()-m_ptPress.y();
    if (dn > 0)
      m_dScale = m_dOldScale/pow(1.002, dn);
    else
      m_dScale = m_dOldScale*pow(1.002, -dn);
    if (m_dScale < 0.25)
      m_dScale = 0.25;
    else if (m_dScale > 10)
      m_dScale = 8;
    UpdateScaledImage();
  }
  else if (m_bDrawing && m_nEditMode == EM_REGION && !m_listRegions.isEmpty())
  {
    const QPoint pt_0 = ScreenToImage(m_ptPress), pt_1 = ScreenToImage(e->pos());
    m_listRegions.last() = RECT_REGION(QPoint(qMin(pt_0.x(), pt_1.x()), qMin(pt_0.y(), pt_1.y())),
                                       QPoint(qMax(pt_0.x(), pt_1.x()), qMax(pt_0.y(), pt_1.y())));
    update();
  }
  else if (m_bDrawing && m_nEditMode == EM_EDIT_MASK)
  {
    FreeHandOnMaskImage(m_ptPrev, e->pos());
    m_ptPrev = e->pos();
  }
  // Follow the held modifiers so the preview updates on release.
  RefreshMaskAction(e->modifiers());
  // On Segmentation only Ctrl arms a mask action (Shift alone does nothing there).
  if (m_nEditMode == EM_EDIT_MASK && (e->modifiers() & Qt::ControlModifier))
  {
    m_bModifierDown = true;
    update();
  }
  else if (m_nEditMode != EM_EDIT_MASK && (e->modifiers() & Qt::ControlModifier || e->modifiers() & Qt::ShiftModifier))
    setCursor(Qt::CrossCursor);
  else if (m_nEditMode == EM_EDIT_MASK)
  {
    // Repaint while hovering so the tool preview follows the cursor.
    update();
  }

  QWidget::mouseMoveEvent(e);
}

// On macOS Qt maps both Key_Control and ControlModifier to Command, so
// Key_Control is correct everywhere. Key_Meta (the Mac Control key) is not
// bound here; accepting it flickered the cursor.
static bool IsMaskModifierKey(int key)
{
  return key == Qt::Key_Control;
}

void WidgetImageView::keyPressEvent(QKeyEvent* e)
{
  if (m_nEditMode == EM_EDIT_MASK && IsMaskModifierKey(e->key()))
  {
    // queryKeyboardModifiers(): during a modifier's own press, e->modifiers() may
    // not include it yet.
    Qt::KeyboardModifiers mods = QGuiApplication::queryKeyboardModifiers();
    m_bModifierDown = true;
    RefreshMaskAction(mods);
    update();
  }
}

void WidgetImageView::keyReleaseEvent(QKeyEvent* e)
{
  if (IsMaskModifierKey(e->key()))
  {
    Qt::KeyboardModifiers mods = QGuiApplication::queryKeyboardModifiers();
    // Only up once no relevant modifier is held.
    m_bModifierDown = (mods & Qt::ControlModifier) != 0;
    RefreshMaskAction(mods);
    update();
  }
}

// Empty on purpose: QWidget's version ignores the event, passing it to the parent.
void WidgetImageView::wheelEvent(QWheelEvent* e)
{
}

void WidgetImageView::resizeEvent(QResizeEvent *e)
{
  UpdateScaledImage();
}

void WidgetImageView::UpdateScaledImage(bool bSmooth)
{
  if (!m_image.isNull() && height() > 0)
  {
    if (1.0*m_image.width()/m_image.height() > 1.0*width()/height())
      m_imageScaled = m_image.scaledToWidth(width()*m_dScale, bSmooth?Qt::SmoothTransformation:Qt::FastTransformation);
    else
      m_imageScaled = m_image.scaledToHeight(height()*m_dScale, bSmooth?Qt::SmoothTransformation:Qt::FastTransformation);
  }
  update();
}

void WidgetImageView::ClearEdits()
{
  m_listPoints.clear();
  m_listRedoPoints.clear();
  m_listRegions.clear();
  m_listRedoRegions.clear();
  m_imageOverlay = QImage();
  HideMessage();
  UpdateAll();
  emit PointEdited();
}

void WidgetImageView::SetMaskOpacity(double val)
{
  m_dMaskOpacity = val;
  UpdateAll();
}

void WidgetImageView::SetCombinedMaskOpacity(double val)
{
  m_dCombinedMaskOpacity = val;
  UpdateAll();
}

void WidgetImageView::SetEditedPoints(const QList<QPoint>& pts)
{
  m_listPoints = pts;
  m_imageOverlay = QImage();
  UpdateAll();
  HideMessage();
}

void WidgetImageView::FreeHandOnMaskImage(const QPoint& scr_pt1, const QPoint& scr_pt2)
{
  QPoint pt1 = ScreenToImage(scr_pt1), pt2 = ScreenToImage(scr_pt2);
  if (pt1 == pt2)
    UpdatePointOnMaskImage(pt1);
  else
  {
    int x0 = pt1.x(), x1 = pt2.x(), y0 = pt1.y(), y1 = pt2.y();
    int dx = x1 - x0;
    int dy = y1 - y0;
    double t = 0.5;
    if ( abs( dx ) > abs( dy ) )
    {
      double m = (double) dy / (double) dx;
      t += y0;
      dx = ( dx < 0 ? -1 : 1 );
      m *= dx;
      while ( x0 != x1 )
      {
        x0 += dx;
        t += m;
        UpdatePointOnMaskImage(QPoint(x0, t));
      }
    }
    else
    {
      double m = (double) dx / (double) dy;
      t += x0;
      dy = ( dy < 0 ? -1 : 1 );
      m *= dy;
      while ( y0 != y1 )
      {
        y0 += dy;
        t += m;
        UpdatePointOnMaskImage(QPoint(t, y0));
      }
    }
  }
  UpdateAll();
}

void WidgetImageView::UpdatePointOnMaskImage(const QPoint& pt_in)
{
  int nSize = CurrentToolSize();
  int nr = nSize/2+1;
  int w = m_imageMask.width(), h = m_imageMask.height();
  QRgb fill_val = QColor(m_bErasing?QColor(0,0,0,0):Qt::white).rgba();
  for (int i = -nr+1; i < nr; i++)
  {
    for (int j = -nr+1; j < nr; j++)
    {
      QPoint pt(pt_in.x()+i, pt_in.y()+j);
      if (pt.x() >= 0 && pt.x() < w && pt.y() >= 0 && pt.y() < h &&
          sqrt(i*i+j*j) <= nSize/2.0)
      {
        QRgb* p = (QRgb*)m_imageMask.scanLine(pt.y());
        p[pt.x()] = fill_val;
      }
    }
  }
  SetMaskEdited();
}

void WidgetImageView::SetMaskEdited()
{
  if (m_bMaskEdited)
    return;
  m_bMaskEdited = true;
  emit MaskEdited();
}

bool WidgetImageView::SaveMaskIfEdited()
{
  if (m_bMaskEdited && !m_imageMask.isNull())
  {
    QImage image = m_imageMask;
    // restore black background
    for (int j = 0; j < image.height(); j++)
    {
      QRgb* p = (QRgb*)image.scanLine(j);
      for (int i = 0; i < image.width(); i++)
      {
        if (qRed(p[i]) == 0)
          p[i] = qRgba(0,0,0,255);
      }
    }
    image.save(m_sMaskFilename);
    m_bMaskEdited = false;
    return true;
  }
  return false;
}

void WidgetImageView::FloodFillMaskImage(const QPoint& pt_in, bool bErase)
{
  QPoint pt = ScreenToImage(pt_in);
  int w = m_imageMask.width(), h = m_imageMask.height();
  if (pt.x() < 0 || pt.x() >= w || pt.y() < 0 || pt.y() >= h)
    return;

  // Fill and erase are the same scan with mask and background swapped.
  // FILL: the mask is the wall; a click inside a hole fills that hole (a click on
  // the mask is a no-op).
  // ERASE: the background is the wall and the flood marks one connected blob. The
  // wall value (2) must differ from the mark (1), or the second pass would clear
  // the background too.
  const char cWall = bErase ? 2 : 1;
  const char cMark = 1;
  char* buf = new char[w*h];
  for (int j = 0; j < h; j++)
  {
    QRgb* p = (QRgb*)m_imageMask.scanLine(j);
    for (int i = 0; i < w; i++)
    {
      const bool bMask = (qRed(p[i]) > 0);
      buf[j*w+i] = (bMask != bErase) ? cWall : 0;
    }
  }
  MyUtils::FloodFill(buf, pt.x(), pt.y(), w, h, cMark, cWall);
  const QRgb valPaint = bErase ? qRgba(0, 0, 0, 0) : QColor(Qt::white).rgb();
  for (int j = 0; j < h; j++)
  {
    QRgb* p = (QRgb*)m_imageMask.scanLine(j);
    for (int i = 0; i < w; i++)
    {
      if (buf[j*w+i] == cMark)
        p[i] = valPaint;
    }
  }
  delete[] buf;
  SetMaskEdited();
  UpdateAll();
}

void WidgetImageView::UpdateAll()
{
  PrepareImage();
  UpdateScaledImage();   // also schedules the repaint
}

void WidgetImageView::UndoLastMaskEdit()
{
  if (!m_listMaskUndoBuffer.isEmpty())
  {
    m_listMaskRedoBuffer << m_imageMask;
    m_imageMask = m_listMaskUndoBuffer.takeLast();
    UpdateAll();
  }
}

void WidgetImageView::RedoLastMaskEdit()
{
  if (!m_listMaskRedoBuffer.isEmpty())
  {
    m_listMaskUndoBuffer << m_imageMask;
    m_imageMask = m_listMaskRedoBuffer.takeLast();
    UpdateAll();
  }
}

void WidgetImageView::BeginPreviewCrossfade()
{
  // Nothing on screen to fade from: don't arm.
  if (m_sFilename.isEmpty() || m_imageScaled.isNull())
    return;

  if (!m_pPreviewFadeAnim)
  {
    m_pPreviewFadeAnim = new QPropertyAnimation(this, "previewFadeOpacity", this);
    m_pPreviewFadeAnim->setDuration(MotionReduced() ? 0 : PREVIEW_FADE_DURATION_MS);
    m_pPreviewFadeAnim->setEasingCurve(QEasingCurve::OutCubic);
    // Connected once: the animation object is reused for every fade.
    connect(m_pPreviewFadeAnim, &QAbstractAnimation::finished, this, [this]()
    {
      m_imageFadeFrom = QImage();
      update();
    });
  }

  // Stop, then snapshot, so rapid toggles restart the single animation.
  m_pPreviewFadeAnim->stop();

  // A standalone copy of the last scaled image (not the on-screen blend), so a
  // very fast double toggle may pop slightly; grabbing the blend would re-enter
  // paintEvent.
  m_imageFadeFrom = m_imageScaled.copy();
  m_bPreviewFadePending = true;
}

void WidgetImageView::StartPreviewCrossfade()
{
  if (m_imageFadeFrom.isNull() || !m_pPreviewFadeAnim)
    return;

  m_pPreviewFadeAnim->stop();
  m_pPreviewFadeAnim->setStartValue(0.0);
  m_pPreviewFadeAnim->setEndValue(1.0);
  m_pPreviewFadeAnim->start();
}

void WidgetImageView::CancelPreviewCrossfade()
{
  if (m_pPreviewFadeAnim)
    m_pPreviewFadeAnim->stop();
  m_imageFadeFrom = QImage();
  // Back to "incoming image fully opaque", so paintEvent()'s guard is false
  // even if a retained frame ever reappears without a fade being started.
  m_dPreviewFadeOpacity = 1.0;
}

void WidgetImageView::Reset()
{
  m_dScale = DEFAULT_FIT_SCALE;
  m_ptOffset = QPoint(0,0);
  UpdateScaledImage();
}

void WidgetImageView::Clear()
{
  // Nothing is displayed: stop the fade.
  m_bPreviewFadePending = false;
  CancelPreviewCrossfade();
  m_sMessage.clear();
  m_sFilename.clear();
  update();
}

void WidgetImageView::UndoLastPoint()
{
  if (!m_listPoints.isEmpty())
  {
    m_listRedoPoints << m_listPoints.takeLast();
    update();
    emit PointEdited();
  }
}

void WidgetImageView::RedoLastPoint()
{
  if (!m_listRedoPoints.isEmpty() && !m_bPointsLocked)
  {
    m_listPoints << m_listRedoPoints.takeLast();
    update();
    emit PointEdited();
  }
}

void WidgetImageView::UndoLastRegion()
{
  if (!m_listRegions.isEmpty())
  {
    m_listRedoRegions << m_listRegions.takeLast();
    update();
  }
}

void WidgetImageView::RedoLastRegion()
{
  if (!m_listRedoRegions.isEmpty())
  {
    m_listRegions << m_listRedoRegions.takeLast();
    update();
  }
}
