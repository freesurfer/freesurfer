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
#include "PanelPointSet.h"
#include "ui_PanelPointSet.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QToolBar>
#include "LayerPointSet.h"
#include "LayerPropertyPointSet.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include <QFileDialog>
#include "MyUtils.h"
#include <QLabel>
#include <QDateTime>
#include <QInputDialog>
#include <QScrollBar>
#include <QDebug>
#include <QTimer>
#include <QTreeWidgetItem>
#include <QMessageBox>
#include "DialogAddPointSetStat.h"
#include "DialogControlPointComment.h"
#include "ToolWindowLesionPopup.h"
#ifdef Q_OS_MAC
#include "MacHelper.h"
#endif

PanelPointSet::PanelPointSet(QWidget *parent) :
  PanelLayer("PointSet", parent),
  ui(new Ui::PanelPointSet)
{
  ui->setupUi(this);
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if (mainwnd)
  {
    ui->toolbar->addAction(mainwnd->ui->actionNewPointSet);
    ui->toolbar->addAction(mainwnd->ui->actionLoadPointSet);
    ui->toolbar->addAction(mainwnd->ui->actionSavePointSet);
    ui->toolbar->addAction(mainwnd->ui->actionClosePointSet);
  }

  m_widgetlistSolidColor << ui->colorpickerSplineColor;

  m_widgetlistHeatScale << ui->comboBoxScalarMap
                        << ui->labelScalarMap
                        << ui->sliderMax
                        << ui->sliderMid
                        << ui->sliderMin
                        << ui->sliderOffset
                        << ui->lineEditMax
                        << ui->lineEditMid
                        << ui->lineEditMin
                        << ui->lineEditOffset
                        << ui->labelMax
                        << ui->labelMid
                        << ui->labelMin
                        << ui->labelOffset;

  m_widgetlistSpline << ui->labelSplineColor
                     << ui->comboBoxSplineColor
                     << ui->lineEditSplineRadius
                     << ui->labelSplineRadius
                     << ui->checkBoxClosedSpline;

  m_self = qgetenv("USER");
  if (m_self.isEmpty())
    m_self = qgetenv("USERNAME");

#ifdef Q_OS_MAC
  if (MacHelper::IsDarkMode())
    ui->commentsContentWidget->setStyleSheet(QString("#commentsContentWidget {background-color:#1E1E1E;}"));
#endif

  m_toolLesionPopup = new ToolWindowLesionPopup(this);
  m_toolLesionPopup->hide();
  connect(m_toolLesionPopup, SIGNAL(GoToPointChanged(int)), SLOT(OnSpinBoxGoToPoint(int)), Qt::QueuedConnection);
  connect(m_toolLesionPopup, SIGNAL(GoToPointTriggered()), SLOT(OnButtonGoToPoint()), Qt::QueuedConnection);
  connect(m_toolLesionPopup, SIGNAL(RadiusTextChanged(QString)), SLOT(OnLineEditRadius(QString)));
  connect(m_toolLesionPopup, SIGNAL(PointColorChanged(QColor)), ui->colorpickerPointColor, SLOT(setCurrentColor(QColor)), Qt::QueuedConnection);
  connect(m_toolLesionPopup, SIGNAL(OpacityChanged(double)), ui->doubleSpinBoxOpacity, SLOT(setValue(double)));

  connect(ui->spinBoxOverallScore, SIGNAL(valueChanged(int)), SLOT(OnSpinBoxOverallScore(int)));
  connect(ui->spinBoxSecondQA, SIGNAL(valueChanged(int)), SLOT(OnSpinBoxSecondQA(int)));
  connect(ui->textEditOverallQuality, SIGNAL(textChanged()), SLOT(OnTextOverallQualityChanged()));
  connect(ui->checkBoxFixed, SIGNAL(toggled(bool)), SLOT(OnCheckBoxFixed(bool)));
  connect(ui->pushButtonLesionPopup, SIGNAL(clicked(bool)), SLOT(OnButtonLesionPopup()));
}

PanelPointSet::~PanelPointSet()
{
  delete ui;
}

void PanelPointSet::ConnectLayer( Layer* layer_in )
{
  PanelLayer::ConnectLayer( layer_in );

  LayerPointSet* layer = qobject_cast<LayerPointSet*>(layer_in);
  if ( !layer )
  {
    return;
  }

  LayerPropertyPointSet* p = layer->GetProperty();
  connect( p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::QueuedConnection );
  connect( layer, SIGNAL(PointAdded(int)), this, SLOT(UpdateWidgets()), Qt::QueuedConnection);
  connect( layer, SIGNAL(PointRemoved(int)), this, SLOT(UpdateWidgets()), Qt::QueuedConnection);
  connect( layer, SIGNAL(PointAdded(int)), this, SLOT(SetCurrentPoint(int)), Qt::QueuedConnection);
  connect( layer, SIGNAL(PointRemoved(int)), this, SLOT(SetCurrentPoint(int)), Qt::QueuedConnection);
  connect( layer, SIGNAL(PointSelected(int)), this, SLOT(SetCurrentPoint(int)), Qt::QueuedConnection);
  connect( ui->doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetOpacity(double)) );
  connect( ui->checkBoxShowSpline, SIGNAL(toggled(bool)), p, SLOT(SetShowSpline(bool)) );
  connect( ui->checkBoxSnapToCenter, SIGNAL(toggled(bool)), p, SLOT(SetSnapToVoxelCenter(bool)));
  connect( ui->colorpickerPointColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetColor(QColor)));
  connect( ui->colorpickerSplineColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetSplineColor(QColor)));
  connect( ui->comboBoxSplineColor, SIGNAL(currentIndexChanged(int)), p, SLOT(SetColorMap(int)));
  connect( ui->checkBoxClosedSpline, SIGNAL(toggled(bool)), p, SLOT(SetClosedSpline(bool)));
  connect( ui->checkBoxShowUnfixedOnly, SIGNAL(toggled(bool)), p, SLOT(SetShowUnfixedOnly(bool)));

  if (m_mapCurrentPoint.contains(layer))
    SetCurrentPoint(m_mapCurrentPoint[layer]);
}

void PanelPointSet::DoIdle()
{
  // update action status
  BlockAllSignals( true );

  BlockAllSignals( false );
}

void PanelPointSet::DoUpdateWidgets()
{
  BlockAllSignals( true );
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar )
    {
      allWidgets[i]->setEnabled(layer);
    }
  }
  int nColorMap = 0;
  bool bShowSpline = false;
  ui->lineEditFileName->clear();
  if ( layer )
  {
    ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );
    double* rgb = layer->GetProperty()->GetColor();
    ui->colorpickerPointColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    rgb = layer->GetProperty()->GetSplineColor();
    ui->colorpickerSplineColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    ui->lineEditFileName->setText( MyUtils::Win32PathProof(layer->GetFileName()) );
    ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );
    ChangeLineEditNumber( ui->lineEditRadius, layer->GetProperty()->GetRadius() );
    ChangeLineEditNumber( ui->lineEditSplineRadius, layer->GetProperty()->GetSplineRadius() );
    ui->labelNumberTotal->setText(QString("%1").arg(layer->GetNumberOfPoints()));
    ui->spinBoxGoToPoint->setRange(1, layer->GetNumberOfPoints());

    nColorMap = layer->GetProperty()->GetColorMap();
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    ui->sliderMin->setValue( (int)( ( layer->GetProperty()->GetHeatScaleMin() - fMin ) / ( fMax - fMin ) * 100 ) );
    ui->sliderMid->setValue( (int)( ( layer->GetProperty()->GetHeatScaleMid() - fMin ) / ( fMax - fMin ) * 100 ) );
    ui->sliderMax->setValue( (int)( ( layer->GetProperty()->GetHeatScaleMax() - fMin ) / ( fMax - fMin ) * 100 ) );
    ui->sliderOffset->setValue( (int)( ( layer->GetProperty()->GetHeatScaleOffset() + fMax ) / ( fMax + fMax ) * 100 ) );
    ChangeLineEditNumber( ui->lineEditMin, layer->GetProperty()->GetHeatScaleMin() );
    ChangeLineEditNumber( ui->lineEditMid, layer->GetProperty()->GetHeatScaleMid() );
    ChangeLineEditNumber( ui->lineEditMax, layer->GetProperty()->GetHeatScaleMax() );
    ChangeLineEditNumber( ui->lineEditOffset, layer->GetProperty()->GetHeatScaleOffset() );

    ui->comboBoxSplineColor->setCurrentIndex( nColorMap );

    ui->comboBoxScalarMap->clear();
    ui->comboBoxScalarMap->addItem("stat");
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" )->GetLayers();
    int nSel = 0;
    for ( int i = 0; i < layers.size(); i++ )
    {
      ui->comboBoxScalarMap->addItem( layers[i]->GetName(), QVariant::fromValue((QObject*)layers[i]) );
      if ( layer->GetProperty()->GetScalarType() == LayerPropertyPointSet::ScalarLayer &&
           layer->GetProperty()->GetScalarLayer() == layers[i] )
      {
        nSel = i+1;
      }
    }
    std::vector<ScalarValues> svs = layer->GetProperty()->GetScalarSets();
    for ( int i = 0; i < (int)svs.size(); i++ )
    {
      ui->comboBoxScalarMap->addItem( svs[i].strName );
      if ( layer->GetProperty()->GetScalarType() == LayerPropertyPointSet::ScalarSet &&
           layer->GetProperty()->GetScalarSet() == i )
      {
        nSel = i+1 + layers.size();
      }
    }
    ui->comboBoxScalarMap->addItem( "Load..." );
    if ( nSel >= 0 )
    {
      ui->comboBoxScalarMap->setCurrentIndex( nSel );
    }

    bShowSpline = layer->GetProperty()->GetShowSpline();
    ui->checkBoxShowSpline->setChecked( bShowSpline );
    ui->checkBoxSnapToCenter->setChecked( layer->GetProperty()->GetSnapToVoxelCenter() );
    ui->labelEndPointDistance->setText(QString("%1 mm").arg(layer->GetEndPointDistance(), 0, 'f', 3));
    ui->checkBoxClosedSpline->setChecked(layer->GetProperty()->GetClosedSpline());

    if (!layer->GetEnhancedData("overall_score").isNull())
      ui->spinBoxOverallScore->setValue(layer->GetEnhancedData("overall_score").toInt());
    else
      ui->spinBoxOverallScore->setValue(1);

    if (!layer->GetEnhancedData("qa_level").isNull())
      ui->spinBoxSecondQA->setValue(layer->GetEnhancedData("qa_level").toInt());
    else
      ui->spinBoxSecondQA->setValue(-1);

    ui->textEditOverallQuality->setPlainText(layer->GetEnhancedData("overall_quality").toString());
  }

  // MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
  ui->colorpickerPointColor->setEnabled( layer );
  ui->comboBoxSplineColor->setEnabled( layer );
  ui->pushButtonCommentAdd->setEnabled(layer && layer->GetNumberOfPoints() > 0);

  ShowWidgets( m_widgetlistSpline, bShowSpline );
  ShowWidgets( m_widgetlistSolidColor, bShowSpline && layer && nColorMap == LayerPropertyPointSet::SolidColor );
  ShowWidgets( m_widgetlistHeatScale, bShowSpline && layer && nColorMap == LayerPropertyPointSet::HeatScale );

  UpdatePointInfo();

  m_toolLesionPopup->UpdateUI(layer);

  BlockAllSignals( false );
}

void PanelPointSet::OnSliderOpacity( int nVal )
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    layer->GetProperty()->SetOpacity( nVal / 100.0 );
  }
}

void PanelPointSet::OnSliderMin(int nVal)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleMin( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }

}

void PanelPointSet::OnSliderMid(int nVal)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleMid( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelPointSet::OnSliderMax(int nVal)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleMax( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelPointSet::OnSliderOffset(int nVal)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleOffset( nVal / 100.0 * ( fMax + fMax ) - fMax );
  }
}

void PanelPointSet::OnLineEditMin(const QString& text)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && layer->GetProperty()->GetHeatScaleMin() != dVal )
    {
      layer->GetProperty()->SetHeatScaleMin( dVal );
    }
  }
}

void PanelPointSet::OnLineEditMid(const QString& text)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && layer->GetProperty()->GetHeatScaleMid() != dVal )
    {
      layer->GetProperty()->SetHeatScaleMid( dVal );
    }
  }
}

void PanelPointSet::OnLineEditMax(const QString& text)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && layer->GetProperty()->GetHeatScaleMax() != dVal )
    {
      layer->GetProperty()->SetHeatScaleMax( dVal );
    }
  }
}

void PanelPointSet::OnLineEditOffset(const QString& text)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && layer->GetProperty()->GetHeatScaleOffset() != dVal )
    {
      layer->GetProperty()->SetHeatScaleOffset( dVal );
    }
  }
}

void PanelPointSet::OnLineEditRadius(const QString& text)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && dVal >= 0 && layer->GetProperty()->GetRadius() != dVal )
    {
      layer->GetProperty()->SetRadius( dVal );
    }
  }
}

void PanelPointSet::OnLineEditSplineRadius(const QString& text)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && dVal > 0 && layer->GetProperty()->GetSplineRadius() != dVal )
    {
      layer->GetProperty()->SetSplineRadius( dVal );
    }
  }
}

void PanelPointSet::OnComboScalarMap(int nSel)
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(ui->comboBoxScalarMap->itemData( nSel ).value<QObject*>());
    if ( mri )
    {
      layer->GetProperty()->SetScalarLayer( mri );
    }
    else if (nSel == 0)
    {
      layer->GetProperty()->SetScalarToStat();
    }
    else if ( nSel == ui->comboBoxScalarMap->count() - 1 )
    {
      LoadScalarValues();
    }
    else
    {
      int offset = ui->comboBoxScalarMap->count()-1 - layer->GetProperty()->GetNumberOfScalarSets();
      layer->GetProperty()->SetScalarSet( nSel - offset );
    }
  }
}

void PanelPointSet::LoadScalarValues()
{
  QList<LayerPointSet*> layers = GetSelectedLayers<LayerPointSet*>();
  foreach (LayerPointSet* layer, layers)
  {
    QString fn = QFileDialog::getOpenFileName( this, "Select scalar file",
                                               "",
                                               "All files (*)");
    if ( !fn.isEmpty() )
    {
      if ( !layer->GetProperty()->LoadScalarsFromFile( fn ) )
      {
        cout << "Load scalar values failed.\n";
      }
    }
    UpdateWidgets();
  }
}

void PanelPointSet::SetCurrentPoint(int nIndex)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if (layer)
  {
    if (nIndex >= layer->GetNumberOfPoints())
      nIndex = layer->GetNumberOfPoints()-1;
    if (nIndex+1 > ui->spinBoxGoToPoint->maximum())
      ui->spinBoxGoToPoint->setMaximum(nIndex+1);
    if (ui->spinBoxGoToPoint->value() != nIndex+1)
    {
      ui->spinBoxGoToPoint->blockSignals(true);
      ui->spinBoxGoToPoint->setValue(nIndex+1);
      ui->spinBoxGoToPoint->blockSignals(false);

      m_toolLesionPopup->SetCurrentPoint(nIndex);

      DoUpdateWidgets();
    }
  }
}

void PanelPointSet::OnSpinBoxGoToPoint(int val)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if (layer && layer->GetNumberOfPoints() > 0)
  {
    double pt[3];
    layer->GetPoint(val-1, pt);
    MainWindow::GetMainWindow()->SetSlicePosition(pt);
    UpdatePointInfo();
  }
}

void PanelPointSet::OnButtonGoToPoint()
{
  OnSpinBoxGoToPoint(ui->spinBoxGoToPoint->value());
}

QLabel* PanelPointSet::MakeCommentItem(const QVariantMap& map, QLabel* label_in)
{
  QLabel* label = label_in;
  if (!label)
    label = new QLabel();
  label->setWordWrap(true);
  label->setTextInteractionFlags(label->textInteractionFlags() | Qt::TextSelectableByMouse);
  bool bDarkMode = false;
#ifdef Q_OS_MAC
  bDarkMode = MacHelper::IsDarkMode();
#endif
  QString prefilled = map.value("prefilled").toStringList().join(" / ");
  if (!map["text"].toString().isEmpty() && !prefilled.isEmpty())
    prefilled = " / " + prefilled;
  QString text = QString("<span style=\"color:rgba(%4,%4,%4,150);font-size:10px;\">[%1] (%2)</span><br />%3 %5")
      .arg(map["timestamp"].toDateTime().toString("yyyy-MM-dd hh:mm:ss"))
      .arg(map["user"].toString()).arg(map["text"].toString()).arg(bDarkMode?255:0).arg(prefilled);
  text += QString(" (<a href=\"edit\" style=\"font-size:11px;color:%1\">edit</a>)").arg(bDarkMode?"#00A6FF":"blue");
  if (map["user"].toString() == m_self)
    text += QString(" (<a href=\"delete\" style=\"font-size:11px;color:%1\">delete</a>)").arg(bDarkMode?"#00A6FF":"blue");
  label->setText(text);
  label->setStyleSheet("QLabel{font-size:12px;padding:2px;padding-top:3px;padding-bottom:3px;}");
  connect(label, SIGNAL(linkActivated(QString)), SLOT(OnCommentLabelClicked(QString)), Qt::UniqueConnection);
  label->setProperty("comment", map);
  return label;
}

void PanelPointSet::ScrollCommentsToBottom()
{
  ui->scrollAreaComments->verticalScrollBar()->setValue(ui->scrollAreaComments->verticalScrollBar()->maximum());
}

QTreeWidgetItem* PanelPointSet::AddStatItem(const QString &name, double value)
{
  QTreeWidgetItem* item = new QTreeWidgetItem(ui->treeWidgetStats);
  item->setText(0, name);
  item->setData(0, Qt::UserRole, name);
  item->setText(1, QString::number(value));
  item->setData(1, Qt::UserRole, value);
  item->setFlags( item->flags() | Qt::ItemIsEditable);
  return item;
}

void PanelPointSet::UpdatePointInfo()
{
  BlockAllSignals(true);
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if ( layer && layer->GetNumberOfPoints() > nIndex )
  {
    ui->labelMoreInfo->setText(tr("Information on Point #%1").arg(nIndex+1));
    // Update comments
    {
      QLayoutItem* item;
      while ( ( item = ui->layoutComments->takeAt( 0 ) ) != NULL )
      {
        delete item->widget();
        delete item;
      }
    }
    ControlPoint p = layer->GetPoint(nIndex);
    ui->checkBoxFixed->setChecked(p.info.value("fixed").toBool());
    QVariantList comments = p.info.value("comments").toList();
    foreach (QVariant v, comments)
    {
      ui->layoutComments->addWidget(MakeCommentItem(v.toMap()));
    }
    ui->layoutComments->addStretch(1);
    ui->scrollAreaComments->widget()->adjustSize();
    QTimer::singleShot(0, this, SLOT(ScrollCommentsToBottom()));

    // Update stats
    ui->treeWidgetStats->clear();
    QTreeWidgetItem* item = AddStatItem("legacy", p.value);
    item->setForeground(0, Qt::gray);
    QVariantMap stats = p.info.value("statistics").toMap();
    QStringList keys = stats.keys();
    foreach (QString key, keys)
    {
      AddStatItem(key, stats[key].toDouble());
    }

    m_mapCurrentPoint[layer] = nIndex;

    m_toolLesionPopup->UpdatePointInfo(nIndex, &p);
  }
  BlockAllSignals(false);
}

void PanelPointSet::OnButtonCommentAdd()
{
  DialogControlPointComment dlg(this);
  if (dlg.exec() == QDialog::Accepted)
  {
    QVariantMap map;
    map["text"] = dlg.GetComment();
    map["prefilled"] = dlg.GetPrefilledItems();
    // workaround for a QDateTime bug
    QDateTime local = QDateTime::currentDateTime();
    QDateTime utc = local.toUTC();
    utc.setTimeSpec(Qt::LocalTime);
    local.setOffsetFromUtc(utc.secsTo(local));
    map["timestamp"] = local;
    map["user"] = m_self;
    QLabel* label = MakeCommentItem(map);
    if (ui->layoutComments->count() > 0)
      ui->layoutComments->insertWidget(ui->layoutComments->count()-1, label);
    else
      ui->layoutComments->addWidget(label);
    ui->scrollAreaComments->widget()->adjustSize();
    QTimer::singleShot(0, this, SLOT(ScrollCommentsToBottom()));

    LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
    int nIndex = ui->spinBoxGoToPoint->value()-1;
    if (layer && nIndex < layer->GetNumberOfPoints())
    {
      ControlPoint p = layer->GetPoint(nIndex);
      QVariantList comments = p.info.value("comments").toList();
      comments << map;
      layer->UpdatePoint(nIndex, "comments", comments);
    }
  }
}

void PanelPointSet::OnCommentLabelClicked(const QString &link)
{
  QLabel* l = qobject_cast<QLabel*>(sender());
  if (!l)
    return;

  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (layer && nIndex < layer->GetNumberOfPoints())
  {
    ControlPoint p = layer->GetPoint(nIndex);
    QVariantList comments = p.info.value("comments").toList();
    if (link == "delete")
    {
      if (QMessageBox::question(this, "Delete Comments", "Are you sure you want to delete this comment?") != QMessageBox::Yes)
        return;

      l->hide();
      for (int i = 0; i < comments.size(); i++)
      {
        if (comments[i].toMap() == l->property("comment").toMap())
        {
          comments.removeAt(i);
          break;
        }
      }
    }
    else if (link == "edit")
    {
      for (int i = 0; i < comments.size(); i++)
      {
        if (comments[i].toMap() == l->property("comment").toMap())
        {
          QVariantMap map = comments[i].toMap();
          DialogControlPointComment dlg(this);
          dlg.SetComment(map["text"].toString(), map["prefilled"].toStringList());
          if (dlg.exec() != QDialog::Accepted)
            return;

          map["text"] = dlg.GetComment();
          map["prefilled"] = dlg.GetPrefilledItems();
          map["user"] = m_self;
          map["edited"] = true;
          QDateTime local = QDateTime::currentDateTime();
          QDateTime utc = local.toUTC();
          utc.setTimeSpec(Qt::LocalTime);
          local.setOffsetFromUtc(utc.secsTo(local));
          map["timestamp"] = local;
          comments[i] = map;
          MakeCommentItem(map, l);
          break;
        }
      }
    }
    layer->UpdatePoint(nIndex, "comments", comments);
  }
}

void PanelPointSet::OnStatItemChanged(QTreeWidgetItem *item, int col)
{
  if (item->text(col) == item->data(col, Qt::UserRole).toString())  // no text change
    return;
  else if (item->text(col).trimmed().isEmpty())
  {
    item->setText(col, item->data(col, Qt::UserRole).toString());
    return;
  }

  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (!layer || nIndex >= layer->GetNumberOfPoints())
    return;

  ControlPoint p = layer->GetPoint(nIndex);
  QVariantMap stats = p.info.value("statistics").toMap();
  QString old_name = item->data(0, Qt::UserRole).toString();
  if (col == 0) // change name
  {
    if (old_name == "legacy")
      item->setText(0, "legacy");
    else if (stats.contains(item->text(0)))
    {
      QMessageBox::information(this, "Change Stat", tr("\'%1\' already exists. Click on the inline text to edit existing stat.").arg(item->text(col)));
    }
    else
    {
      item->setData(0, Qt::UserRole, item->text(0));
      stats.remove(old_name);
      stats[item->text(0)] = item->data(1, Qt::UserRole);
      layer->UpdatePoint(nIndex, "statistics", stats);
    }
  }
  else  // change value
  {
    bool bOk;
    double new_val = item->text(col).toDouble(&bOk);
    if (!bOk)
    {
      item->setText(col, item->data(col, Qt::UserRole).toString());
      return;
    }

    item->setData(col, Qt::UserRole, new_val);
    if (old_name == "legacy")
    {
      layer->UpdatePoint(nIndex, "legacy_stat", new_val);
    }
    else
    {
      stats[item->text(0)] = new_val;
      layer->UpdatePoint(nIndex, "statistics", stats);
    }
  }
}

void PanelPointSet::OnButtonStatAdd()
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (layer && nIndex < layer->GetNumberOfPoints())
  {
    ControlPoint p = layer->GetPoint(nIndex);
    QVariantMap stats = p.info.value("statistics").toMap();
    DialogAddPointSetStat dlg(this);
    if (dlg.exec() == QDialog::Accepted)
    {
      if (stats.contains(dlg.GetStatName()))
      {
        QMessageBox::information(this, "Add Statistic", tr("\'%1\' already exists. Click on the inline text to edit existing stat.").arg(dlg.GetStatName()));
      }
      else
      {
        AddStatItem(dlg.GetStatName(), dlg.GetStatValue());
        stats[dlg.GetStatName()] = dlg.GetStatValue();
        layer->UpdatePoint(nIndex, "statistics", stats);
      }
    }
  }
}

void PanelPointSet::OnButtonStatDelete()
{
  QTreeWidgetItem* item = ui->treeWidgetStats->currentItem();
  if (item && ui->treeWidgetStats->indexOfTopLevelItem(item) != 0)
  {
    LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
    int nIndex = ui->spinBoxGoToPoint->value()-1;
    if (layer && nIndex < layer->GetNumberOfPoints())
    {
      ControlPoint p = layer->GetPoint(nIndex);
      QVariantMap stats = p.info.value("statistics").toMap();
      stats.remove(item->text(0));
      layer->UpdatePoint(nIndex, "statistics", stats);
      ui->treeWidgetStats->takeTopLevelItem(ui->treeWidgetStats->indexOfTopLevelItem(item));
      delete item;
    }
  }
}

void PanelPointSet::OnCheckBoxFixed(bool b)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (layer && nIndex < layer->GetNumberOfPoints())
  {
    layer->UpdatePoint(nIndex, "fixed", b);
  }
}

void PanelPointSet::OnCurrentStatItemChanged(QTreeWidgetItem *cur, QTreeWidgetItem *old)
{
  Q_UNUSED(cur);
  Q_UNUSED(old);
  ui->pushButtonStatDelete->setEnabled(cur && ui->treeWidgetStats->indexOfTopLevelItem(cur) != 0);
}

void PanelPointSet::OnTextOverallQualityChanged()
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if (layer)
    layer->SetEnhancedData("overall_quality", ui->textEditOverallQuality->toPlainText());
}


void PanelPointSet::OnSpinBoxOverallScore(int val)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if (layer)
    layer->SetEnhancedData("overall_score", val);
}

void PanelPointSet::OnSpinBoxSecondQA(int val)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if (layer)
    layer->SetEnhancedData("qa_level", val);
}

void PanelPointSet::OnButtonLesionPopup()
{
  m_toolLesionPopup->show();
  m_toolLesionPopup->raise();
}
