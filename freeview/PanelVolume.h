/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#ifndef PANELVOLUME_H
#define PANELVOLUME_H

#include "PanelLayer.h"
#include <QList>
#include <QVector>

#include "colortab.h"

class LayerMRI;

namespace Ui
{
class PanelVolume;
}

class LUTDataHolder;
class QTreeWidget;

class ColorTableItem : public QTreeWidgetItem
{
public:
  explicit ColorTableItem(int type = Type) : QTreeWidgetItem(type) {}
  explicit ColorTableItem(QTreeWidget* tree) : QTreeWidgetItem(tree) {}

  enum SORT_TYPE  { ST_VALUE = 0, ST_NAME };

  virtual bool operator < ( const QTreeWidgetItem& other ) const;

  static int  SortType;
  static bool SortAscending;
};

class PanelVolume : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelVolume(QWidget *parent = 0);
  ~PanelVolume();

  bool eventFilter(QObject *watched, QEvent *event);

  QList<LayerMRI*> GetLinkedVolumes();

protected slots:
  void OnCheckShowContour( bool bShow );
  void OnCheckShowLabelContour( bool bShow );
  void OnSliderOpacity( int nVal );
  void OnComboColorMap( int nSel );
  void OnComboLookupTable( int nSel );
  void OnColorTableCurrentItemChanged( QTreeWidgetItem* item );
  void OnColorTableItemClicked( QTreeWidgetItem* item);
  void OnColorTableItemDoubleClicked( QTreeWidgetItem* item = NULL );
  void OnLineEditBrushValue( const QString& strg = NULL );
  void OnCheckBoxSelectAllLabels(int nState);
  void OnColorTableItemChanged( QTreeWidgetItem* item );
  void OnColorTableSortingChanged();

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
  void OnCopySettings();
  void OnPasteSettings();
  void OnPasteSettingsToAll();
  void OnSliderTrackVolumeMin(int);
  void OnTrackVolumeThresholdChanged();
  void OnLockLayer(bool);

  void UpdateColorLabel();
  void UpdateTrackVolumeThreshold();

  void OnActiveFrameChanged(int nFrame);

  void OnShowExistingLabelsOnly(bool b);

  void ShowAllLabels()
  {
      OnShowExistingLabelsOnly(false);
  }

  void OnComboMask( int nSel );

  void OnComboCorrelationSurface(int nSel);

  void OnCheckUsePercentile(bool b);

  void OnLineEditVectorDisplayScale(const QString& strg);

  void OnLineEditVectorLineWidth(const QString& strg);

  void OnLineEditVectorNormThreshold(const QString &strg);

  void OnLineEditProjectionMapRangeChanged();

  void OnComboProjectionMapType(int nType);

  void OnLineEditMaskThreshold(const QString& text);

  void OnCustomContextMenu(const QPoint& pt);

  void OnButtonResetWindowLevel();

  void OnCheckBoxSetDisplayVector(bool b);

  void OnCheckBoxSetDisplayTensor(bool b);

  void OnCheckBoxSetNormalizeVector(bool b);

  void OnCheckBoxSetDisplayRGB(bool b);

  void OnGoToFirstPoint();

  void OnGoToNextPoint();

  void OnColorTableChangeColor();

  void OnCheckVoxelizedContour( bool bVoxelize );

  void OnCheckBoxSetAutoMid(bool b);

  void UpdateOpacity(double val);

  void OnLineEditClearBackgroundValue(const QString& text);

  void RefreshColorTable();

protected:
  void PopulateColorTable( COLOR_TABLE* ctab, bool bForce = false );
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
  QList<QWidget*> m_widgetlistContourNormal;
  QList<QWidget*> m_widgetlistNormalDisplay;
  QList<QWidget*> m_widgetlistEditable;
  QList<QWidget*> m_widgetlistVolumeTrack;
  QList<QWidget*> m_widgetlistVolumeTrackSpecs;
  QList<QWidget*> m_widgetlistNonVolumeTrack;

  LUTDataHolder* m_luts;

  COLOR_TABLE*  m_curCTAB;
  bool          m_bShowExistingLabelsOnly;

  QVector<double> m_voxelList;
  int           m_nCurrentVoxelIndex;
};

#endif // PANELVOLUME_H
