#ifndef APPMOTION_H
#define APPMOTION_H

// Animation durations, curves, the reduced-motion switch and shared
// primitives. GUI thread only.
//  * Under Reduced() every primitive completes within the call (deterministic
//    screenshots; a kill switch for users).
//  * Crossfade*/AnimateHeight run from the wall clock (the shared Ticker), one
//    frame ahead, not Qt's animation driver, so the click's own paint shows the
//    first step; the fade ends at duration - 16 ms.
//  * QGraphicsOpacityEffect only on short-lived Crossfade* veils.

#include <QEasingCurve>
#include <QObject>
#include <functional>

class QStackedWidget;
class QWidget;

namespace AppMotion
{

enum Duration { Fast = 120, Base = 180, Slow = 240 };   // ms

QEasingCurve Enter();   // appearing, moving in, growing
QEasingCurve Exit();    // fading out
QEasingCurve Leave();   // moving off, collapsing

// DISSECTION_PHOTO_REDUCED_MOTION=1, or the "offscreen" platform unless the
// variable is "0". Read once.
bool Reduced();
int  Ms(Duration d);              // d, or 0 when Reduced()
int  Ms(int ms);                  // ms, or 0 when Reduced()

// Monotonic milliseconds since the first call.
qint64 Now();

// Stacked-page cross-fade for a caller about to block: snapshots the current
// page into a veil, switches the stack at once, and fades the veil out after
// one event-loop turn, so the fade doesn't elapse during the blocking work
// (which holds only while the caller doesn't call processEvents()).
void CrossfadeStackAfterWork(QStackedWidget* stack, QWidget* next, Duration d = Base);

// Snapshot w, run `mutate`, fade the snapshot out.
void CrossfadeInPlace(QWidget* w, const std::function<void()>& mutate, Duration d = Fast);

// Expand/collapse by animating min and max height together. Shows w first if
// `from` is 0; hides at the end if `to` is 0 and hideAtEnd; then restores min 0
// / max QWIDGETSIZE_MAX.
void AnimateHeight(QWidget* w, int from, int to, Duration d = Slow, bool hideAtEnd = true);

// Shared 16 ms timer, running only while it has (live) subscribers.
class Ticker : public QObject
{
  Q_OBJECT
public:
  static Ticker* Instance();
  static void Subscribe(QObject* who);

signals:
  void tick(qint64 nowMs);

private:
  Ticker();
  void OnTimeout();
  void Prune();
  struct Private;
  Private* d;
};

} // namespace AppMotion

#endif // APPMOTION_H
