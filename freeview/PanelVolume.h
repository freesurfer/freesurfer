/**
 * @file  PanelVolume.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/04 18:06:50 $
 *    $Revision: 1.42 $
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

extern "C"
{
#include "colortab.h"
}

namespace Ui
{
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
  void OnCheckShowLabelContour( bool bShow );
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
  void OnSliderTrackVolumeMin(int);
  void OnTrackVolumeThresholdChanged();

  void UpdateColorLabel();
  void UpdateTrackVolumeThreshold();

  void OnActiveFrameChanged(int nFrame);

  void OnShowExistingLabelsOnly(bool b);

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
  QList<QWidget*> m_widgetlistContourNormal;
  QList<QWidget*> m_widgetlistContourLabel;
  QList<QWidget*> m_widgetlistNormalDisplay;
  QList<QWidget*> m_widgetlistEditable;
  QList<QWidget*> m_widgetlistVolumeTrack;
  QList<QWidget*> m_widgetlistVolumeTrackSpecs;
  QList<QWidget*> m_widgetlistNonVolumeTrack;

  LUTDataHolder* m_luts;

  COLOR_TABLE*  m_curCTAB;
  bool          m_bShowExistingLabelsOnly;
};

#endif // PANELVOLUME_H
