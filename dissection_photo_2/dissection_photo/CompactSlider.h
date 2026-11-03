#ifndef COMPACTSLIDER_H
#define COMPACTSLIDER_H

#include <QWidget>
#include <QString>
#include <QColor>
#include <QPixmap>

// Compact scrubber slider: an optional swatch, then a dark rounded track with a
// lighter value fill, label left and value right. Click or drag on the track
// to set the value. The track stays dark with light text in both themes.
class CompactSlider : public QWidget
{
  Q_OBJECT
public:
  explicit CompactSlider(QWidget* parent = nullptr);

  void SetRange(int min, int max);
  void SetValue(int val);
  int Value() const { return m_nValue; }

  void SetLabel(const QString& label);
  void SetSuffix(const QString& suffix);
  void SetSwatchColor(const QColor& color);
  // Hides the swatch and gives its space to the track.
  void SetSwatchVisible(bool bVisible);

  // Draws a tool glyph (a resource path) in the swatch instead of the flat
  // colour; an empty string restores the colour.
  void SetSwatchIcon(const QString& strResourcePath);

  QSize sizeHint() const override;
  QSize minimumSizeHint() const override;

signals:
  void valueChanged(int val);
  void swatchClicked();

protected:
  void paintEvent(QPaintEvent* e) override;
  void mousePressEvent(QMouseEvent* e) override;
  void mouseMoveEvent(QMouseEvent* e) override;
  void mouseReleaseEvent(QMouseEvent* e) override;

private:
  // The widget minus the bottom strip reserved for the drop shadows.
  QRect ContentRect() const;
  QRect TrackRect() const;
  // Full widget height over the track's x-range, so the shadow strip stays
  // clickable.
  QRect TrackHitRect() const;
  QRect SwatchRect() const;
  void UpdateValueFromPos(int x);
  void PaintDropShadow(QPainter& p, const QRect& rc, int nRadius) const;
  void PaintSwatch(QPainter& p, const QRect& rcSwatch, bool bDark);
  void EnsureIconPixmap();

  int     m_nMin;
  int     m_nMax;
  int     m_nValue;
  QString m_sLabel;
  QString m_sSuffix;
  QColor  m_colorSwatch;
  bool    m_bDragging;
  bool    m_bSwatchVisible;

  QString m_strIconPath;
  QPixmap m_pixmapIcon;
  // DPR the cached pixmap was loaded for; reloads when it changes.
  qreal   m_dIconDpr;

  static const int SWATCH_SIZE = 22;
  static const int SWATCH_SPACING = 6;
};

#endif // COMPACTSLIDER_H
