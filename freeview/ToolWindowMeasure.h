/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#ifndef TOOLWINDOWMEASURE_H
#define TOOLWINDOWMEASURE_H

#include "UIUpdateHelper.h"
#include <QWidget>
#include <QList>

namespace Ui
{
class ToolWindowMeasure;
}

class Region2D;
class SurfaceRegion;
class Region3D;

class ToolWindowMeasure : public QWidget, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit ToolWindowMeasure(QWidget *parent = 0);
  ~ToolWindowMeasure();

  QString GetLabelStats();

public slots:
  void SetRegion(Region2D* reg = 0);
  void SetSurfaceRegion( SurfaceRegion* reg = 0 );
  void Set3DRegion( Region3D* reg = 0 );

protected slots:
  void OnIdle();
  void OnAction(QAction* act);
  void UpdateWidgets();

  void OnLoad();
  void OnSave();
  void OnSaveAll();
  void OnUpdate();
  void OnCopy();
  void OnExport();
  void OnSpinBoxId(int val);
  void OnSpinBoxGroup(int val);
  void OnColorGroup( const QColor& color );
  void OnButtonDelete3DRegion();
  void OnButtonDeleteAll3DRegions();

protected:
  virtual void showEvent(QShowEvent *);

private:
  Ui::ToolWindowMeasure *ui;

  QList<QWidget*> m_widgets2D;
  QList<QWidget*> m_widgets3D;
  QList<QWidget*> m_widgets3DDraw;
  Region2D*       m_region;
  SurfaceRegion*  m_surfaceRegion;
  Region3D*       m_3DRegion;
  bool            m_bToUpdateWidgets;
};

#endif // TOOLWINDOWMEASURE_H
