#ifndef TOOLWINDOWLESIONPOPUP_H
#define TOOLWINDOWLESIONPOPUP_H

#include <QWidget>
#include "UIUpdateHelper.h"
#include <QPointer>

class LayerPointSet;
struct ControlPoint;

namespace Ui {
class ToolWindowLesionPopup;
}

class ToolWindowLesionPopup : public QWidget, public UIUpdateHelper
{
  Q_OBJECT

  public:
  explicit ToolWindowLesionPopup(QWidget *parent = nullptr);
  ~ToolWindowLesionPopup();

  signals:
  void GoToPointChanged(int);
  void GoToPointTriggered();
  void PointColorChanged(const QColor& c);
  void RadiusTextChanged(const QString& text);
  void OpacityChanged(double val);

  public:
  void UpdateUI(LayerPointSet* ps);
  void UpdatePointInfo(int nIndex, ControlPoint* p);
  void SetCurrentPoint(int nIndex);

  public slots:
  void OnCheckBoxShowRegion(bool b);
  void OnCheckBoxSetRegion();
  void OnSpinBoxGoToPoint(int n);
  void OnTextAdditionalComment();
  void OnCheckBoxNonLesionReason(bool b);
  void OnConfidenceScoreChanged(int);
  void OnSliderOpacity(int n);

  private:
  void BlockAllSignals(bool bBlock);

  private:
  Ui::ToolWindowLesionPopup *ui;
  QPointer<LayerPointSet>  m_curLayer;
};

#endif // TOOLWINDOWLESIONPOPUP_H
