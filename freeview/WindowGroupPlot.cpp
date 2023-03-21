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

#include "WindowGroupPlot.h"
#include "ui_WindowGroupPlot.h"
#include <QSettings>
#include "FSGroupDescriptor.h"
#include <QListWidgetItem>
#include <QPainter>

WindowGroupPlot::WindowGroupPlot(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::WindowGroupPlot)
{
  ui->setupUi(this);
  this->setWindowFlags(Qt::Tool);
  this->setWindowTitle("Plot");
  connect(ui->comboBoxViewBy, SIGNAL(currentIndexChanged(int)), SLOT(OnComboViewBy(int)));
  connect(ui->comboBoxConfigClass, SIGNAL(currentIndexChanged(int)), SLOT(OnComboConfigClass(int)));
  connect(ui->comboBoxConfigShape, SIGNAL(currentTextChanged(QString)), SLOT(OnComboConfigShape(QString)));
  connect(ui->widgetConfigColor, SIGNAL(colorChanged(QColor)), SLOT(OnConfigColor(QColor)));
  connect(ui->listWidget, SIGNAL(currentItemChanged(QListWidgetItem*,QListWidgetItem*)),
          SLOT(OnCurrentItemChanged(QListWidgetItem*)));
  connect(ui->widgetPlot, SIGNAL(CurrentDataIndexChanged(int)), SLOT(OnCurrentDataIndexChanged(int)));

  QSettings s;
  QVariant v = s.value("WindowPlot/Geomerty");
  if (v.isValid())
    this->restoreGeometry(v.toByteArray());
}

WindowGroupPlot::~WindowGroupPlot()
{
  QSettings s;
  s.setValue("WindowPlot/Geomerty", this->saveGeometry());
  delete ui;
}

void WindowGroupPlot::SetFsgdData(FSGroupDescriptor *fsgd)
{
  QWidgetList widgets = findChildren<QWidget*>();
  foreach (QWidget* w, widgets)
    w->blockSignals(true);

  m_fsgd = fsgd;
  ui->widgetPlot->SetFsgdData(fsgd);
  setWindowTitle(fsgd->m_title);
  UpdateStockPixmaps();
  ui->comboBoxVariable->clear();
  for (int i = 0; i < fsgd->m_variables.size(); i++)
    ui->comboBoxVariable->addItem(fsgd->m_variables[i].label);
  ui->comboBoxVariable->setCurrentIndex(0);

  ui->comboBoxConfigClass->clear();
  for (int i = 0; i < fsgd->m_classes.size(); i++)
    ui->comboBoxConfigClass->addItem(fsgd->m_classes[i].label);
  ui->comboBoxConfigClass->setCurrentIndex(0);

  UpdateCurrentConfig(fsgd->m_classes[0].marker, fsgd->m_classes[0].color);

  OnComboViewBy(0);

  foreach (QWidget* w, widgets)
    w->blockSignals(false);
}

void WindowGroupPlot::UpdateCurrentConfig(const QString& shape, const QColor& c)
{
  ui->comboBoxConfigShape->blockSignals(true);
  ui->widgetConfigColor->blockSignals(true);
  for (int i = 0; i < ui->comboBoxConfigShape->count(); i++)
  {
    if (shape == ui->comboBoxConfigShape->itemText(i))
    {
      ui->comboBoxConfigShape->setCurrentIndex(i);
      break;
    }
  }
  ui->widgetConfigColor->setCurrentColor(c);
  ui->widgetConfigColor->blockSignals(false);
  ui->comboBoxConfigShape->blockSignals(false);
}

void WindowGroupPlot::SetCurrentVertex(int nVertex)
{
  if (isVisible())
  {
    ui->widgetPlot->SetCurrentVertex(nVertex);
    ui->labelVertexNumber->setText(QString("Vertex # %1").arg(nVertex));
  }
}

void WindowGroupPlot::OnComboViewBy(int nIndex)
{
  ui->listWidget->clear();
  if (nIndex == 0)
  {
    ui->listWidget->setFixedWidth(100);
    ui->listWidget->setWrapping(false);
    QList<FSGDClass> classes = m_fsgd->m_classes;
    for (int i = 0; i < classes.size(); i++)
    {
      QListWidgetItem* item = new QListWidgetItem;
      item->setText(classes[i].label);
      item->setIcon(QIcon(m_listMarkerPixmaps[i]));
      item->setData(Qt::UserRole, classes[i].label);
      ui->listWidget->addItem(item);
    }
  }
  else
  {
    ui->listWidget->setFixedWidth(180);
    ui->listWidget->setWrapping(true);
    QList<FSGDDataItem> data = m_fsgd->m_data;
    for (int i = 0; i < data.size(); i++)
    {
      QListWidgetItem* item = new QListWidgetItem;
      item->setIcon(QIcon(m_listMarkerPixmaps[data[i].class_id]));
      item->setText(data[i].subject_id);
      item->setData(Qt::UserRole, data[i].subject_id);
      ui->listWidget->addItem(item);
    }
  }
}

void WindowGroupPlot::resizeEvent(QResizeEvent *e)
{
  QWidget::resizeEvent(e);
  ui->listWidget->setWrapping(ui->listWidget->isWrapping());
}

void WindowGroupPlot::UpdateStockPixmaps()
{
  m_listMarkerPixmaps.clear();
  foreach (FSGDClass c, m_fsgd->m_classes)
  {
    QPixmap pixmap(16,16);
    pixmap.fill(QColor(0,0,0,0));
    QPainter p(&pixmap);
    WidgetGroupPlot::DrawMarker(&p, QPoint(8,8), c.marker, c.color, 6.0);
    p.end();
    m_listMarkerPixmaps << pixmap;
  }
}

void WindowGroupPlot::OnComboConfigClass(int nIndex)
{
  UpdateCurrentConfig(m_fsgd->m_classes[nIndex].marker, m_fsgd->m_classes[nIndex].color);
}

void WindowGroupPlot::OnComboConfigShape(const QString &strg)
{
  int n = ui->comboBoxConfigClass->currentIndex();
  m_fsgd->m_classes[n].marker = strg;
  ui->widgetPlot->update();
  UpdateStockPixmaps();
  OnComboViewBy(ui->comboBoxViewBy->currentIndex());
}

void WindowGroupPlot::OnConfigColor(const QColor &c)
{
  int n = ui->comboBoxConfigClass->currentIndex();
  m_fsgd->m_classes[n].color = c;
  ui->widgetPlot->update();
  UpdateStockPixmaps();
  OnComboViewBy(ui->comboBoxViewBy->currentIndex());
}

void WindowGroupPlot::OnCurrentItemChanged(QListWidgetItem *item)
{
  if (ui->comboBoxViewBy->currentIndex() > 0)
    ui->widgetPlot->SetCurrentDataIndex(ui->listWidget->row(item));
}

void WindowGroupPlot::OnCurrentDataIndexChanged(int nIndex)
{
  if (nIndex >= 0 && ui->comboBoxViewBy->currentIndex() > 0)
    ui->listWidget->setCurrentItem(ui->listWidget->item(nIndex));
}
