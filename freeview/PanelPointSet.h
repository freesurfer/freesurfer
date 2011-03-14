#ifndef PANELPOINTSET_H
#define PANELPOINTSET_H

#include "PanelLayer.h"

namespace Ui {
    class PanelPointSet;
}

class PanelPointSet : public PanelLayer
{
    Q_OBJECT

public:
    explicit PanelPointSet(QWidget *parent = 0);
    ~PanelPointSet();

protected:
    void DoUpdateWidgets();
    void DoIdle();
    virtual void ConnectLayer( Layer* layer );
    void LoadScalarValues();

protected slots:
    void OnSliderOpacity( int nVal );
    void OnSliderMin(int nVal);
    void OnSliderMid(int nVal);
    void OnSliderMax(int nVal);
    void OnSliderOffset(int nVal);
    void OnLineEditMin(const QString& text);
    void OnLineEditMid(const QString& text);
    void OnLineEditMax(const QString& text);
    void OnLineEditOffset(const QString& text);
    void OnLineEditRadius(const QString& text);
    void OnLineEditSplineRadius(const QString& text);
    void OnComboScalarMap(int nSel);

private:
    Ui::PanelPointSet *ui;
    QList<QWidget*> m_widgetlistSolidColor;
    QList<QWidget*> m_widgetlistHeatScale;
    QList<QWidget*> m_widgetlistSpline;
};

#endif // PANELPOINTSET_H
