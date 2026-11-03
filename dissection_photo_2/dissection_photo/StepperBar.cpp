#include "StepperBar.h"
#include "CommonDef.h"
#include "AppMotion.h"
#include "AppTheme.h"

#include <QFontMetrics>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>

namespace
{
// The app's existing colours, so the bar matches the window; source noted.
struct Pal
{
  QColor bar, accent, onAccent, accentText, success;
  QColor textBody, textSecondary, textDisabled;
  QColor nHover, nPress, focusRing;
  int radius, fsLabel, fsNumeral;
};

Pal CurrentPal()
{
  const bool d = AppTheme::EffectiveDarkMode();
  Pal t;
  t.bar           = QColor(d ? "#383838" : "#fafafa");   // #frameTopBar
  t.accent        = QColor("#457c92");                   // primary QPushButton
  t.onAccent      = QColor("#ffffff");                   // its text
  t.accentText    = QColor(d ? "#6fb3c9" : "#35637a");   // button:pressed / QLineEdit:focus
  t.success       = QColor(d ? "#8bd48b" : "#55a355");   // "processed" outline, darker in light for contrast
  t.textBody      = QColor(d ? "#f0f0f0" : "#262626");   // #frameTopBar QLineEdit
  t.textSecondary = QColor(d ? "#b0b0b0" : "#5f5f5f");   // derived: between body and disabled
  t.textDisabled  = QColor(d ? "#787878" : "#a0a0a0");   // QMenu::item:disabled
  t.nHover        = QColor(d ? "#464646" : "#ececec");   // lightbox back-button pill
  t.nPress        = QColor(d ? "#565656" : "#dcdcdc");   // same pill, pressed
  t.focusRing     = QColor("#6fb3c9");                   // QLineEdit:focus
  t.radius        = 6;
  t.fsLabel       = 13;
  t.fsNumeral     = 12;
  return t;
}

const int kHeight = 32;
const int kGutter = 16;
const int kDisc = 20;
const int kGap = 8;         // disc -> label
const int kStepGap = 32;    // between steps (the connector lives here)
const int kPillInset = 8;   // hover pill grows this much past disc/label
const int kPillHeight = 24;
const int kLeadingMargin = 12;  // bar edge -> Home button
const int kLeadingGap = 12;     // Home button -> first step's hover pill
}

StepperBar::StepperBar(QWidget* parent)
  : QWidget(parent)
{
  setFocusPolicy(Qt::TabFocus);
  setMouseTracking(true);
  setAttribute(Qt::WA_Hover, true);
  setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
}

void StepperBar::SetSteps(const QStringList& names)
{
  m_items.clear();
  for (const QString& s : names)
  {
    Item it;
    it.name = s;
    m_items << it;
  }
  m_nHover = m_nPressed = m_nFocus = -1;
  Layout();
  updateGeometry();
  update();
}

void StepperBar::SetStates(const QVector<State>& states)
{
  bool bChanged = false;
  for (int i = 0; i < m_items.size(); i++)
  {
    const State s = i < states.size() ? states[i] : Future;
    if (m_items[i].state != s)
      bChanged = true;
  }
  if (!bChanged)
    return;
  auto apply = [this, states]()
  {
    for (int i = 0; i < m_items.size(); i++)
      m_items[i].state = i < states.size() ? states[i] : Future;
    if (m_nFocus >= 0 && !Clickable(m_nFocus))
      m_nFocus = NextClickable(-1, 1);
    update();
  };
  if (isVisible())
    AppMotion::CrossfadeInPlace(this, apply, AppMotion::Fast);
  else
    apply();
}

QSize StepperBar::sizeHint() const
{
  // Wide enough for every label at its Semibold width, so a label never
  // moves when it becomes current.
  int w = StepsLeft() + kGutter;
  QFontMetrics fm(LabelFont(true));
  for (int i = 0; i < m_items.size(); i++)
    w += kDisc + kGap + fm.horizontalAdvance(m_items[i].name) + (i ? kStepGap : 0);
  return QSize(w, kHeight);
}

QSize StepperBar::minimumSizeHint() const
{
  return QSize(0, kHeight);
}

QFont StepperBar::LabelFont(bool bCurrent) const
{
  const Pal t = CurrentPal();
  QFont f = font();
  f.setPixelSize(t.fsLabel);
  f.setWeight(bCurrent ? QFont::DemiBold : QFont::Normal);
  return f;
}

void StepperBar::Layout()
{
  // Measured Semibold for every step so geometry doesn't change with state;
  // labels elide in paintEvent when too narrow.
  QFontMetrics fm(LabelFont(true));
  const int y = (height() - kDisc) / 2;
  int x = StepsLeft();
  for (int i = 0; i < m_items.size(); i++)
  {
    Item& it = m_items[i];
    // +2: horizontalAdvance() rounds down, so eliding at exactly this width
    // could clip the string it was measured from.
    const int wLabel = fm.horizontalAdvance(it.name) + 2;
    it.rcDisc = QRect(x, y, kDisc, kDisc);
    it.rcLabel = QRect(x + kDisc + kGap, 0, wLabel, height());
    it.rcPill = QRect(x - kPillInset, (height() - kPillHeight) / 2,
                      kDisc + kGap + wLabel + kPillInset * 2, kPillHeight);
    x += kDisc + kGap + wLabel + kStepGap;
  }
}

int StepperBar::IndexAt(const QPoint& pt) const
{
  for (int i = 0; i < m_items.size(); i++)
    if (m_items[i].rcPill.contains(pt))
      return i;
  return -1;
}

bool StepperBar::Clickable(int i) const
{
  return isEnabled() && i >= 0 && i < m_items.size() && m_items[i].state == Done;
}

int StepperBar::NextClickable(int from, int step) const
{
  for (int i = from + step; i >= 0 && i < m_items.size(); i += step)
    if (Clickable(i))
      return i;
  return -1;
}

void StepperBar::Activate(int i)
{
  if (Clickable(i))
  {
    m_nFocus = i;   // a click that lands moves the keyboard position with it
    emit StepActivated(i);
  }
}

void StepperBar::paintEvent(QPaintEvent*)
{
  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing, true);
  const Pal t = CurrentPal();
  p.fillRect(rect(), t.bar);
  if (m_items.isEmpty())
    return;

  const bool bEnabled = isEnabled();
  const bool bKeyFocus = hasFocus() || property("focused").toString() == "true";
  const QFont fontLabel = LabelFont(false), fontCurrent = LabelFont(true);
  QFont fontNumeral = font();
  // Same mono font list as the canvas badges.
  FsSetFontFamilies(fontNumeral, QStringList() << "Menlo" << "DejaVu Sans Mono" << "Noto Sans Mono" << "Liberation Mono");
  fontNumeral.setPixelSize(t.fsNumeral);

  for (int i = 0; i < m_items.size(); i++)
  {
    const Item& it = m_items[i];
    if (i + 1 < m_items.size())
    {
      const int x0 = it.rcPill.right() + 1 + kPillInset / 2;
      const int x1 = m_items[i + 1].rcPill.left() - kPillInset / 2;
      const int yc = height() / 2;
      p.setPen(QPen(t.nPress, 1.0));
      p.drawLine(QPoint(x0, yc), QPoint(x1, yc));
    }

    if (Clickable(i) && (m_nHover == i || m_nPressed == i))
    {
      QPainterPath pill;
      pill.addRoundedRect(it.rcPill, t.radius, t.radius);
      p.fillPath(pill, m_nPressed == i ? t.nPress : t.nHover);
    }

    const QRectF rcDisc(it.rcDisc);
    p.setPen(Qt::NoPen);
    if (it.state == Current)
    {
      p.setBrush(t.accent);
      p.drawEllipse(rcDisc);
      p.setFont(fontNumeral);
      p.setPen(t.onAccent);
      p.drawText(rcDisc, Qt::AlignCenter, QString::number(i + 1));
    }
    else if (it.state == Done)
    {
      // Ring plus tick: a ring alone reads as not done.
      const QColor cTick = bEnabled ? t.success : t.textDisabled;
      p.setPen(QPen(cTick, 1.5));
      p.setBrush(Qt::NoBrush);
      p.drawEllipse(rcDisc.adjusted(1, 1, -1, -1));
      QPainterPath tick;
      tick.moveTo(rcDisc.left() + rcDisc.width() * 0.30, rcDisc.top() + rcDisc.height() * 0.52);
      tick.lineTo(rcDisc.left() + rcDisc.width() * 0.44, rcDisc.top() + rcDisc.height() * 0.68);
      tick.lineTo(rcDisc.left() + rcDisc.width() * 0.72, rcDisc.top() + rcDisc.height() * 0.34);
      p.setPen(QPen(cTick, 2.0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
      p.drawPath(tick);
    }
    else
    {
      p.setPen(QPen(t.nPress, 1.0));
      p.setBrush(Qt::NoBrush);
      p.drawEllipse(rcDisc.adjusted(0.5, 0.5, -0.5, -0.5));
      p.setFont(fontNumeral);
      p.setPen(t.textDisabled);
      p.drawText(rcDisc, Qt::AlignCenter, QString::number(i + 1));
    }

    const bool bCurrent = it.state == Current;
    p.setFont(bCurrent ? fontCurrent : fontLabel);
    QColor cLabel = t.textDisabled;
    if (bCurrent)
      cLabel = t.accentText;
    else if (it.state == Done && bEnabled)
      cLabel = (m_nHover == i || m_nPressed == i) ? t.textBody : t.textSecondary;
    p.setPen(cLabel);
    const int nRight = width() - kGutter;
    const int wAvail = qMax(0, qMin(it.rcLabel.width(), nRight - it.rcLabel.left()));
    const QFontMetrics fmLabel(p.font());
    const QString strLabel = fmLabel.horizontalAdvance(it.name) <= wAvail
                               ? it.name
                               : fmLabel.elidedText(it.name, Qt::ElideRight, wAvail);
    p.drawText(QRect(it.rcLabel.left(), 0, wAvail, height()), Qt::AlignLeft | Qt::AlignVCenter, strLabel);

    if (bKeyFocus && m_nFocus == i)
    {
      QPainterPath ring;
      ring.addRoundedRect(QRectF(it.rcPill).adjusted(1, 1, -1, -1), t.radius, t.radius);
      p.setPen(QPen(t.focusRing, 2.0));
      p.setBrush(Qt::NoBrush);
      p.drawPath(ring);
    }
  }
}

void StepperBar::mouseMoveEvent(QMouseEvent* e)
{
  const int n = IndexAt(e->pos());
  if (n != m_nHover)
  {
    m_nHover = n;
    setCursor(Clickable(n) ? Qt::PointingHandCursor : Qt::ArrowCursor);
    update();
  }
  QWidget::mouseMoveEvent(e);
}

void StepperBar::mousePressEvent(QMouseEvent* e)
{
  if (e->button() != Qt::LeftButton)
  {
    QWidget::mousePressEvent(e);
    return;
  }
  const int n = IndexAt(e->pos());
  if (Clickable(n))
  {
    m_nPressed = n;
    // Not m_nFocus: the ring tracks the keyboard; Activate() moves it when a
    // click lands.
    update();
  }
}

void StepperBar::mouseReleaseEvent(QMouseEvent* e)
{
  const int n = m_nPressed;
  m_nPressed = -1;
  update();
  if (n >= 0 && IndexAt(e->pos()) == n)
    Activate(n);
}

void StepperBar::leaveEvent(QEvent* e)
{
  m_nHover = -1;
  m_nPressed = -1;
  setCursor(Qt::ArrowCursor);
  update();
  QWidget::leaveEvent(e);
}

void StepperBar::keyPressEvent(QKeyEvent* e)
{
  int step = 0;
  if (e->key() == Qt::Key_Left || e->key() == Qt::Key_Up)
    step = -1;
  else if (e->key() == Qt::Key_Right || e->key() == Qt::Key_Down)
    step = 1;
  if (step != 0)
  {
    const int n = NextClickable(m_nFocus < 0 ? (step > 0 ? -1 : m_items.size()) : m_nFocus, step);
    if (n >= 0)
    {
      m_nFocus = n;
      update();
    }
    e->accept();
    return;
  }
  if (e->key() == Qt::Key_Space || e->key() == Qt::Key_Return || e->key() == Qt::Key_Enter)
  {
    if (Clickable(m_nFocus))
    {
      Activate(m_nFocus);
      e->accept();
      return;
    }
  }
  QWidget::keyPressEvent(e);
}

void StepperBar::focusInEvent(QFocusEvent* e)
{
  if (!Clickable(m_nFocus))
    m_nFocus = e->reason() == Qt::BacktabFocusReason ? NextClickable(m_items.size(), -1) : NextClickable(-1, 1);
  update();
  QWidget::focusInEvent(e);
}

void StepperBar::focusOutEvent(QFocusEvent* e)
{
  update();
  QWidget::focusOutEvent(e);
}

void StepperBar::resizeEvent(QResizeEvent* e)
{
  Layout();
  PlaceLeading();
  QWidget::resizeEvent(e);
}

void StepperBar::SetLeadingWidget(QWidget* w)
{
  m_pLeading = w;
  if (w)
    w->setParent(this);
  Layout();
  PlaceLeading();
  updateGeometry();
  if (w)
    w->show();
}

QSize StepperBar::LeadingSize() const
{
  // Its real size: the size hint alone is narrower than a fixed-width button.
  return m_pLeading ? m_pLeading->sizeHint().expandedTo(m_pLeading->minimumSize()).boundedTo(m_pLeading->maximumSize())
                    : QSize();
}

int StepperBar::StepsLeft() const
{
  return m_pLeading ? kLeadingMargin + LeadingSize().width() + kLeadingGap + kPillInset : kGutter;
}

void StepperBar::PlaceLeading()
{
  if (!m_pLeading)
    return;
  const QSize sz = LeadingSize();
  m_pLeading->setGeometry(kLeadingMargin, (height() - sz.height()) / 2, sz.width(), sz.height());
  update();
}

void StepperBar::changeEvent(QEvent* e)
{
  // QWidget::changeEvent() repaints for all three.
  if (e->type() == QEvent::EnabledChange || e->type() == QEvent::FontChange || e->type() == QEvent::StyleChange)
    Layout();
  QWidget::changeEvent(e);
}
