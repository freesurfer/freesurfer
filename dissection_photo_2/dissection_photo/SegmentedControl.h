#ifndef SEGMENTEDCONTROL_H
#define SEGMENTEDCONTROL_H

#include <QWidget>
#include <QStringList>
#include <QColor>

class QFontMetrics;
class QPropertyAnimation;

// Pill track split into N equal segments; the selection pill slides between
// them (QPropertyAnimation on "indicatorPos").
class SegmentedControl : public QWidget
{
  Q_OBJECT
  Q_PROPERTY(qreal indicatorPos READ IndicatorPos WRITE SetIndicatorPos)

public:
  explicit SegmentedControl(QWidget* parent = nullptr);

  // Sets the segment labels. Resets the current index to 0 (no animation).
  void SetItems(const QStringList& items);

  int CurrentIndex() const { return m_nCurrentIndex; }

  // Emits currentIndexChanged() only on a change.
  void SetCurrentIndex(int index, bool animate = true);

  // Per-segment enabled state. A disabled segment is dimmed, ignores clicks and
  // arrow keys, and can't become current. SetItems() re-enables all.
  void SetItemEnabled(int index, bool bEnabled);
  bool IsItemEnabled(int index) const;

  // One tooltip per segment.
  void SetItemToolTips(const QStringList& tips) { m_tips = tips; }

  // Selection pill fill; defaults to brand red.
  void SetAccent(const QColor& color);

  // Track tones for light and dark themes; the defaults suit #frameTopBar.
  void SetTrackColors(const QColor& light, const QColor& dark);

  QSize sizeHint() const override;

  // Width at which the bold selected label still fits its indicator.
  QSize minimumSizeHint() const override;

  // In segment units, so geometry stays correct on resize mid-animation.
  qreal IndicatorPos() const { return m_dIndicatorPos; }
  void SetIndicatorPos(qreal pos);

signals:
  void currentIndexChanged(int index);

protected:
  void paintEvent(QPaintEvent*) override;
  void mousePressEvent(QMouseEvent*) override;
  void keyPressEvent(QKeyEvent*) override;
  bool event(QEvent*) override;

private:
  QRectF IndicatorRect() const;
  qreal SegmentWidth() const;

  // Segment 'index' inset by kMargin per side (as wide as the pill); never
  // negative width.
  QRectF SegmentTextRect(int index) const;

  // Width a segment needs for 'str' in bold plus the 2*kMargin inset.
  int SegmentMinWidth(const QFontMetrics& fmBold, const QString& str) const;

  QStringList m_items;
  QStringList m_tips;
  QList<bool> m_enabled;  // may be shorter than m_items; missing = enabled
  int m_nCurrentIndex;
  qreal m_dIndicatorPos;
  QPropertyAnimation* m_pAnimation;
  QColor m_colorAccent;
  QColor m_colorTrackLight;
  QColor m_colorTrackDark;
};

#endif // SEGMENTEDCONTROL_H
