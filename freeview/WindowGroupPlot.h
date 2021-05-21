/**
 * @brief Tool window to plot group data
 *
 */
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

#ifndef WINDOWGROUPPLOT_H
#define WINDOWGROUPPLOT_H

#include <QWidget>
#include <QPixmap>
#include <QVariantMap>

namespace Ui {
class WindowGroupPlot;
}

class FSGroupDescriptor;
class QListWidgetItem;

class WindowGroupPlot : public QWidget
{
  Q_OBJECT

public:
  explicit WindowGroupPlot(QWidget *parent = 0);
  ~WindowGroupPlot();

  void SetFsgdData(FSGroupDescriptor* fsgd);

  void resizeEvent(QResizeEvent *e);

public slots:
  void SetCurrentVertex(int nVertex);
  void OnComboViewBy(int nIndex);
  void OnComboConfigClass(int nIndex);
  void OnComboConfigShape(const QString& strg);
  void OnConfigColor(const QColor& c);
  void OnCurrentItemChanged(QListWidgetItem* item);
  void OnCurrentDataIndexChanged(int nIndex);

private:
  void UpdateStockPixmaps();
  void UpdateCurrentConfig(const QString& shape, const QColor& c);

  Ui::WindowGroupPlot *ui;
  FSGroupDescriptor*  m_fsgd;
  QList<QPixmap>    m_listMarkerPixmaps;
};

#endif // WINDOWGROUPPLOT_H
