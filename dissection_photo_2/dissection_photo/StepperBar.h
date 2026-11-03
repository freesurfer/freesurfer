#ifndef STEPPERBAR_H
#define STEPPERBAR_H

// Pipeline stepper: a 32px custom-painted strip atop the main window showing
// the five steps. Done steps are clickable (StepActivated); current shows the
// accent disc; future steps are inert. Tab focus, Left/Right between clickable
// steps, Return/Space activate; the 2px focus ring sits inside the hover pill
// so nothing moves.

#include <QVector>
#include <QWidget>

class StepperBar : public QWidget
{
  Q_OBJECT
public:
  enum State { Future = 0, Done = 1, Current = 2 };

  explicit StepperBar(QWidget* parent = nullptr);

  void SetSteps(const QStringList& names);
  // One state per step; a shorter vector leaves the rest Future. Crossfades
  // when anything changed.
  void SetStates(const QVector<State>& states);
  // Reparents w (the Home button) to the bar's left end; the steps start after it.
  void SetLeadingWidget(QWidget* w);

  QSize sizeHint() const override;
  QSize minimumSizeHint() const override;

signals:
  void StepActivated(int index);

protected:
  void paintEvent(QPaintEvent*) override;
  void mouseMoveEvent(QMouseEvent*) override;
  void mousePressEvent(QMouseEvent*) override;
  void mouseReleaseEvent(QMouseEvent*) override;
  void leaveEvent(QEvent*) override;
  void keyPressEvent(QKeyEvent*) override;
  void focusInEvent(QFocusEvent*) override;
  void focusOutEvent(QFocusEvent*) override;
  void resizeEvent(QResizeEvent*) override;
  void changeEvent(QEvent*) override;

private:
  struct Item
  {
    QString name;
    State state = Future;
    QRect rcPill;   // hover / focus / hit rect
    QRect rcDisc;
    QRect rcLabel;
  };
  void Layout();
  int  IndexAt(const QPoint& pt) const;
  bool Clickable(int i) const;
  int  NextClickable(int from, int step) const;
  QFont LabelFont(bool bCurrent) const;
  void Activate(int i);

  QVector<Item> m_items;
  int m_nHover = -1;
  int m_nPressed = -1;
  int m_nFocus = -1;   // keyboard cursor; -1 = none
  QWidget* m_pLeading = nullptr;
  QSize LeadingSize() const;
  int  StepsLeft() const;
  void PlaceLeading();
};

#endif // STEPPERBAR_H
