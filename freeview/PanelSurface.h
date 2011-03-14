#ifndef PANELSURFACE_H
#define PANELSURFACE_H

#include "PanelLayer.h"
#include <QList>

namespace Ui {
    class PanelSurface;
}

class QAction;
class WindowConfigureOverlay;

class PanelSurface : public PanelLayer
{
    Q_OBJECT

public:
    explicit PanelSurface(QWidget *parent = 0);
    ~PanelSurface();

protected:
    void DoUpdateWidgets();
    void DoIdle();
    virtual void ConnectLayer( Layer* layer );

protected slots:
    void OnChangeSurfaceType( QAction* act );
    void OnSliderOpacity( int nVal );
    void OnSliderMidPoint( int nVal );
    void OnSliderSlope( int nVal );
    void OnComboCurvature( int nSel );
    void OnLineEditMidPoint( const QString& text );
    void OnLineEditSlope( const QString& text );
    void OnComboOverlay( int nSel );
    void OnComboAnnotation( int nSel );
    void OnComboLabel( int nSel );
    void OnComboVector( int nSel );
    void OnButtonConfigureOverlay();
    void OnEditPositionOffset();

private:
    Ui::PanelSurface *ui;

    QList<QWidget*>  m_widgetsMidPoint;
    QList<QWidget*>  m_widgetsSlope;
    QList<QWidget*>  m_widgetsVector;
    QList<QWidget*>  m_widgetsVertex;
    QList<QWidget*>  m_widgetsMesh;
    QList<QWidget*>  m_widgetsLabel;

    WindowConfigureOverlay* m_wndConfigureOverlay;
};

#endif // PANELSURFACE_H
