#ifndef PANELVOLUME_H
#define PANELVOLUME_H

#include "PanelLayer.h"
#include <QList>

extern "C"
{
#include "colortab.h"
}

namespace Ui {
  class PanelVolume;
}

class LUTDataHolder;

class PanelVolume : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelVolume(QWidget *parent = 0);
  ~PanelVolume();

protected slots:
    void OnCheckShowContour( bool bShow );
    void OnSliderOpacity( int nVal );
    void OnComboColorMap( int nSel );
    void OnComboLookupTable( int nSel );
    void OnColorTableCurrentItemChanged( QTreeWidgetItem* item );
    void OnLineEditBrushValue( const QString& strg );

    void OnSliderWindow( int );
    void OnSliderLevel( int );
    void OnSliderMin( int );
    void OnSliderMid( int );
    void OnSliderMax( int );
    void OnSliderOffset( int );
    void OnLineEditWindow( const QString& text );
    void OnLineEditLevel( const QString& text );
    void OnLineEditMin( const QString& text );
    void OnLineEditMid( const QString& text );
    void OnLineEditMax( const QString& text );
    void OnLineEditOffset( const QString& text );
    void OnSliderContourMin(int);
    void OnSliderContourMax(int);
    void OnSliderContourSmooth(int);
    void OnContourValueChanged();
    void OnContourSave();
    void OnCopySettings();
    void OnPasteSettings();
    void OnPasteSettingsToAll();

    void UpdateColorLabel();

protected:
    void PopulateColorTable( COLOR_TABLE* ctab );
    void DoUpdateWidgets();
    void DoIdle();
    virtual void ConnectLayer( Layer* layer );

private:
  Ui::PanelVolume *ui;

  QList<QWidget*> m_widgetlistGrayScale;
  QList<QWidget*> m_widgetlistHeatScale;
  QList<QWidget*> m_widgetlistGenericColorMap;
  QList<QWidget*> m_widgetlistLUT;
  QList<QWidget*> m_widgetlistDirectionCode;
  QList<QWidget*> m_widgetlistFrame;
  QList<QWidget*> m_widgetlistVector;
  QList<QWidget*> m_widgetlistContour;
  QList<QWidget*> m_widgetlistNormalDisplay;
  QList<QWidget*> m_widgetlistEditable;

  LUTDataHolder* m_luts;

  COLOR_TABLE*  m_curCTAB;
};

#endif // PANELVOLUME_H
