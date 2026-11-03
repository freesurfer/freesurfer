#include "AppMotion.h"

#include <QAbstractAnimation>
#include <QElapsedTimer>
#include <QEvent>
#include <QGraphicsOpacityEffect>
#include <QGuiApplication>
#include <QLabel>
#include <QPixmap>
#include <QPointer>
#include <QPropertyAnimation>
#include <QStackedWidget>
#include <QTimer>
#include <QVariantAnimation>
#include <QWidget>
#include <QList>

namespace AppMotion
{

QEasingCurve Enter() { return QEasingCurve(QEasingCurve::OutCubic); }
QEasingCurve Exit()  { return QEasingCurve(QEasingCurve::OutQuad); }
QEasingCurve Leave() { return QEasingCurve(QEasingCurve::InCubic); }

namespace
{
int g_nReduced = -1;   // -1 = not yet read

bool ReadReduced()
{
  const QByteArray env = qgetenv("DISSECTION_PHOTO_REDUCED_MOTION").trimmed();
  if (env == "1" || env.compare("true", Qt::CaseInsensitive) == 0)
    return true;
  if (env == "0" || env.compare("false", Qt::CaseInsensitive) == 0)
    return false;
  // Not set: an offscreen platform (the QA harness) wants final states.
  if (QGuiApplication::instance() && QGuiApplication::platformName() == QLatin1String("offscreen"))
    return true;
  return qgetenv("QT_QPA_PLATFORM").trimmed() == "offscreen";
}
}

bool Reduced()
{
  if (g_nReduced < 0)
    g_nReduced = ReadReduced() ? 1 : 0;
  return g_nReduced == 1;
}

int Ms(Duration d) { return Reduced() ? 0 : int(d); }
int Ms(int ms)     { return Reduced() ? 0 : ms; }

qint64 Now()
{
  static QElapsedTimer s_timer;
  if (!s_timer.isValid())
    s_timer.start();
  return s_timer.elapsed();
}

namespace
{
// A snapshot label over a host. Dies with its fade, with a host resize, or
// when the next veil replaces it.
class Veil : public QLabel
{
public:
  explicit Veil(QWidget* host) : QLabel(host), m_host(host)
  {
    setObjectName(QStringLiteral("AppMotionVeil"));
    setAttribute(Qt::WA_TransparentForMouseEvents);
    setAttribute(Qt::WA_NoSystemBackground);
    host->installEventFilter(this);
  }
  ~Veil() override
  {
    if (m_host)
      m_host->removeEventFilter(this);
  }
  bool eventFilter(QObject* watched, QEvent* e) override
  {
    if (watched == m_host && e->type() == QEvent::Resize)
    {
      // A resized snapshot looks wrong; finish now.
      hide();
      deleteLater();
    }
    return QLabel::eventFilter(watched, e);
  }
private:
  QPointer<QWidget> m_host;
};

// Deletes any veil on host now: the snapshot must not linger for even a frame
// under the new one.
void FinishVeils(QWidget* host)
{
  // By name: Veil has no Q_OBJECT, so findChildren<Veil*> is unavailable.
  const QList<QLabel*> veils = host->findChildren<QLabel*>(QStringLiteral("AppMotionVeil"), Qt::FindDirectChildrenOnly);
  for (QLabel* v : veils)
  {
    v->hide();
    delete v;
  }
}

// Qt's animation driver starts on a queued call, so the first visible step
// trails a click by ~36 ms. Drive() starts the animation (so it owns its
// values, emits finished() and honours DeleteWhenStopped), then pauses it to
// leave Qt's driver; the Ticker sets currentTime = elapsed + kLead on each tick.
// A one-shot half-frame step covers a paint just before the first tick. The
// drive is the animation's child, so ending the animation ends it.
class WallClockDrive : public QObject
{
public:
  static const int kLead = 16;
  explicit WallClockDrive(QAbstractAnimation* a) : QObject(a), m_a(a), m_t0(Now())
  {
    QObject::connect(Ticker::Instance(), &Ticker::tick, this, [this](qint64 now) { Step(now); });
    Ticker::Subscribe(this);
    QTimer::singleShot(kLead / 2, this, [this]() { Step(Now()); });
    Step(m_t0);
  }
private:
  void Step(qint64 now)
  {
    if (!m_a || m_a->state() == QAbstractAnimation::Stopped)
      return;
    // Reaching totalDuration stops the animation and emits finished().
    m_a->setCurrentTime(qMin(int(now - m_t0) + kLead, m_a->totalDuration()));
  }
  QPointer<QAbstractAnimation> m_a;
  qint64 m_t0;
};

void Drive(QAbstractAnimation* a)
{
  a->start(QAbstractAnimation::DeleteWhenStopped);
  if (a->state() != QAbstractAnimation::Running)   // zero duration: finished in start()
    return;
  a->pause();
  new WallClockDrive(a);
}

void FadeVeilOut(Veil* veil, Duration d)
{
  auto* fx = new QGraphicsOpacityEffect(veil);
  fx->setOpacity(1.0);
  veil->setGraphicsEffect(fx);
  auto* a = new QPropertyAnimation(fx, "opacity", veil);
  a->setStartValue(1.0);
  a->setEndValue(0.0);
  a->setDuration(Ms(d));
  a->setEasingCurve(Exit());
  QPointer<Veil> guard(veil);
  QObject::connect(a, &QAbstractAnimation::finished, veil, [guard]()
  {
    if (guard)
    {
      guard->hide();
      guard->deleteLater();
    }
  });
  Drive(a);
}

Veil* CoverWithSnapshot(QWidget* host)
{
  QPixmap snap = host->grab();
  Veil* veil = new Veil(host);
  veil->setPixmap(snap);
  veil->setGeometry(host->rect());
  veil->raise();
  veil->show();
  return veil;
}

const char* kHeightAnimName = "AppMotionHeight";

void StopAnim(QWidget* w, const char* name)
{
  if (QAbstractAnimation* a = w->findChild<QAbstractAnimation*>(QLatin1String(name), Qt::FindDirectChildrenOnly))
  {
    // An interrupted animation must not run its completion handler.
    QObject::disconnect(a, &QAbstractAnimation::finished, nullptr, nullptr);
    a->stop();     // DeleteWhenStopped posts a deleteLater...
    delete a;      // ...delete now so the name is free; Qt drops the posted event
  }
}
} // namespace

void CrossfadeStackAfterWork(QStackedWidget* stack, QWidget* next, Duration d)
{
  if (!stack || !next)
    return;
  FinishVeils(stack);
  if (stack->currentWidget() == next || Reduced() || !stack->isVisible() || !stack->currentWidget())
  {
    stack->setCurrentWidget(next);
    return;
  }
  Veil* veil = CoverWithSnapshot(stack);
  stack->setCurrentWidget(next);     // live underneath from frame 0
  veil->raise();
  QPointer<Veil> guard(veil);
  QTimer::singleShot(0, veil, [guard, d]()
  {
    if (guard)
      FadeVeilOut(guard, d);
  });
}

void CrossfadeInPlace(QWidget* w, const std::function<void()>& mutate, Duration d)
{
  if (!w)
    return;
  FinishVeils(w);
  if (Reduced() || !w->isVisible())
  {
    if (mutate)
      mutate();
    return;
  }
  Veil* veil = CoverWithSnapshot(w);
  if (mutate)
    mutate();
  veil->raise();
  FadeVeilOut(veil, d);
}

void AnimateHeight(QWidget* w, int from, int to, Duration d, bool hideAtEnd)
{
  if (!w)
    return;
  StopAnim(w, kHeightAnimName);
  from = qMax(0, from);
  to = qMax(0, to);
  auto finish = [](QWidget* w, int to, bool hideAtEnd)
  {
    if (to == 0 && hideAtEnd)
      w->hide();
    w->setMinimumHeight(0);
    w->setMaximumHeight(QWIDGETSIZE_MAX);
  };
  if (Reduced() || from == to)
  {
    if (to > 0)
      w->show();
    finish(w, to, hideAtEnd);
    return;
  }
  w->setMinimumHeight(from);
  w->setMaximumHeight(from);
  if (from == 0)
    w->show();
  auto* a = new QVariantAnimation(w);
  a->setObjectName(QLatin1String(kHeightAnimName));
  a->setStartValue(from);
  a->setEndValue(to);
  a->setDuration(Ms(d));
  a->setEasingCurve(to > from ? Enter() : Leave());
  QPointer<QWidget> guard(w);
  QObject::connect(a, &QVariantAnimation::valueChanged, w, [guard](const QVariant& v)
  {
    if (!guard)
      return;
    const int h = v.toInt();
    guard->setMinimumHeight(h);
    guard->setMaximumHeight(h);
  });
  QObject::connect(a, &QAbstractAnimation::finished, w, [guard, to, hideAtEnd, finish]()
  {
    if (guard)
      finish(guard, to, hideAtEnd);
  });
  Drive(a);
}

struct Ticker::Private
{
  QTimer timer;
  QList<QPointer<QObject>> subs;
};

Ticker::Ticker() : QObject(nullptr), d(new Private)
{
  d->timer.setInterval(16);
  d->timer.setTimerType(Qt::PreciseTimer);
  connect(&d->timer, &QTimer::timeout, this, &Ticker::OnTimeout);
}

Ticker* Ticker::Instance()
{
  static Ticker* s_instance = new Ticker();   // lives for the process; no parent needed
  return s_instance;
}

void Ticker::Prune()
{
  for (int i = d->subs.size() - 1; i >= 0; --i)
    if (d->subs[i].isNull())
      d->subs.removeAt(i);
  if (d->subs.isEmpty())
    d->timer.stop();
}

void Ticker::OnTimeout()
{
  Prune();
  if (!d->subs.isEmpty())
    emit tick(Now());
}

void Ticker::Subscribe(QObject* who)
{
  if (!who)
    return;
  Ticker* t = Instance();
  t->Prune();
  for (const QPointer<QObject>& s : t->d->subs)
    if (s.data() == who)
      return;
  t->d->subs.append(QPointer<QObject>(who));
  if (!t->d->timer.isActive())
    t->d->timer.start();
}

} // namespace AppMotion
