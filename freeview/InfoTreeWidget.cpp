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
#include "InfoTreeWidget.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPLabel.h"
#include "LayerPropertyMRI.h"
#include "LayerPropertySurface.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include "FSVolume.h"
#include "SurfaceOverlay.h"
#include "SurfaceAnnotation.h"
#include "MyUtils.h"
#include "LayerProperty.h"
#include <QTreeWidgetItem>
#include <QLineEdit>
#include <QHeaderView>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QShowEvent>
#include <QDebug>
#include <QMenu>
#include "RenderView3D.h"
#include "SurfacePath.h"
#include "MigrationDefs.h"


InfoTreeWidget::InfoTreeWidget(QWidget* parent) :
  QTreeWidget(parent),
  m_bShowSurfaceCurvature(false),
  m_bShowSurfaceNormal(false),
  m_bShowTkRegRAS(true),
  m_bForCursor(false),
  m_bShowSelectedLayersOnly(false)
{
  this->setAlternatingRowColors(true);
  this->setTextElideMode(Qt::ElideMiddle);
  m_editor = new QLineEdit(this);
  m_editor->hide();
  connect(this, SIGNAL(itemClicked(QTreeWidgetItem*,int)),
          this, SLOT(OnItemClicked(QTreeWidgetItem*,int)), Qt::QueuedConnection);
  connect(m_editor, SIGNAL(returnPressed()), this, SLOT(OnEditFinished()));
  connect(MainWindow::GetMainWindow()->GetLayerCollection("MRI"), SIGNAL(ActiveLayerChanged(Layer*)),
          this, SLOT(UpdateAll()), Qt::QueuedConnection);
}

void InfoTreeWidget::OnCursorPositionChanged()
{
  MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetSlicePosition(m_dRAS);
  UpdateAll();
}

void InfoTreeWidget::OnMousePositionChanged()
{
  MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetCurrentRASPosition(m_dRAS);
  UpdateAll();
}

void InfoTreeWidget::showEvent(QShowEvent * e)
{
  // hack to fix a qdesigner bug
  headerItem()->setText(1, "");
  QTreeWidget::showEvent(e);
}

void InfoTreeWidget::ShowHeaderText()
{
  headerItem()->setText(1, "Hold Shift to update");
  QFont fnt = headerItem()->font(1);
  fnt.setPointSize(8);
  headerItem()->setFont(1, fnt);
}

void InfoTreeWidget::ClearHeaderText()
{
  headerItem()->setText(1, "");
}

void InfoTreeWidget::SetShowSelectedLayersOnly(bool b)
{
  m_bShowSelectedLayersOnly = b;
  UpdateAll();
}

void InfoTreeWidget::OnLayerSelectionChanged()
{
  if (m_bShowSelectedLayersOnly)
    UpdateAll();
}

void InfoTreeWidget::UpdateAll()
{
  this->clear();
  m_editor->hide();
  MainWindow* mainWnd =  MainWindow::GetMainWindow();
  LayerCollection* lc_mri = mainWnd->GetLayerCollection( "MRI" );
  LayerCollection* lc_surf = mainWnd->GetLayerCollection( "Surface" );

  QList<Layer*> sel_mris = mainWnd->GetSelectedLayers("MRI");

  int nPrecision = MainWindow::GetMainWindow()->GetSetting("Precision").toInt();
  bool bComma = MainWindow::GetMainWindow()->GetSetting("UseComma").toBool();
  if ( lc_mri->IsEmpty() && lc_surf->IsEmpty())
  {
    return;
  }

  QTreeWidgetItem* item = new QTreeWidgetItem(this);
  item->setText(0, "RAS");

  double ras[3] = {m_dRAS[0], m_dRAS[1], m_dRAS[2]};
  if (!lc_mri->IsEmpty())
  {
    qobject_cast<LayerMRI*>(lc_mri->GetActiveLayer())->RemapPositionToRealRAS(m_dRAS, ras);
  }
  QVariantMap map;
  item->setText(1, QString("%1, %2, %3")
                .arg(ras[0], 0, 'f', 2)
      .arg(ras[1], 0, 'f', 2)
      .arg(ras[2], 0, 'f', 2));
  map["Type"] = "RAS";
  map["EditableText"] = item->text(1);
  item->setData(1, Qt::UserRole, map);

  if (!lc_mri->IsEmpty())
  {
    double tkRegRAS[3];
    LayerMRI* mri = qobject_cast<LayerMRI*>(lc_mri->GetActiveLayer());
    if (m_bShowTkRegRAS)
    {
      mri->NativeRASToTkReg(ras, tkRegRAS);
      item = new QTreeWidgetItem(this);
      item->setText(0, QString("TkReg RAS (%1)").arg(mri->GetName()));
      map.clear();
      item->setText(1, QString("%1, %2, %3")
                    .arg(tkRegRAS[0], 0, 'f', 2)
          .arg(tkRegRAS[1], 0, 'f', 2)
          .arg(tkRegRAS[2], 0, 'f', 2));
      map["Type"] = "TkRegRAS";
      map["EditableText"] = item->text(1);
      item->setData(1, Qt::UserRole, map);
    }
    FSVolume*vol = mri->GetSourceVolume();
    double tpos[3];
    if (vol->RASToTalairach(ras, tpos))
    {
      item = new QTreeWidgetItem(this);
      item->setText(0, QString("MNI305 (%1)").arg(mri->GetName()));
      map.clear();
      item->setText(1, QString("%1, %2, %3")
                    .arg(tpos[0], 0, 'f', 2)
          .arg(tpos[1], 0, 'f', 2)
          .arg(tpos[2], 0, 'f', 2));
      map["Type"] = "Talairach";
      map["EditableText"] = item->text(1);
      item->setData(1, Qt::UserRole, map);
    }
  }

  bool bDecimalIndex = MainWindow::GetMainWindow()->GetSetting("DecimalVoxelCoord").toBool();
  for (int i = 0; i < lc_mri->GetNumberOfLayers(); i++)
  {
    LayerMRI* layer = (LayerMRI*)lc_mri->GetLayer(i);
    if (m_bShowSelectedLayersOnly && !sel_mris.contains(layer))
      continue;
    double fIndex[3];
    if ( layer->GetProperty()->GetShowInfo() )
    {
      QTreeWidgetItem* item = new QTreeWidgetItem(this);
      item->setText(0, layer->GetName());
      layer->RASToOriginalIndex( ras, fIndex );
      double dvalue;
      //      if (layer->IsModified() || layer->GetCorrelationSurface())
      dvalue = layer->GetVoxelValue( m_dRAS );
      //      else
      //        dvalue = layer->GetVoxelValueByOriginalIndex(nIndex[0]+0.5, nIndex[1]+0.5, nIndex[2]+0.5);
      QString valueStrg = MyUtils::RealToNumber(dvalue, nPrecision);
      if (layer->GetNumberOfFrames() > 1 && layer->GetNumberOfFrames() <= 6)
      {
        QList<double> values = layer->GetVoxelValueByOriginalIndexAllFrames((int)(fIndex[0]+0.5), (int)(fIndex[1]+0.5), (int)(fIndex[2]+0.5));
        if (layer->GetDataType() == MRI_RGB)
        {
          int nval = (int)values[0];
          values.clear();
          values << (nval & 0x00ff) << ((nval >> 8) & 0x00ff) << ((nval >> 16) & 0x00ff);
        }
        QStringList strgs;
        for (int n = 0; n < values.size(); n++)
        {
          if (n == layer->GetActiveFrame())
            strgs << QString("*%1*").arg(MyUtils::RealToNumber(values[n], nPrecision));
          else
            strgs << MyUtils::RealToNumber(values[n], nPrecision);
        }
        valueStrg = strgs.join(bComma?", ":" ");
        if (values.size() == 6)
        {
          valueStrg = "("+strgs.mid(0, 3).join(bComma?", ":" ") + ") (" + strgs.mid(3,3).join(bComma?", ":" ")+")";
        }
      }
      for (int j = 0; j < 3; j++)
        fIndex[j] = ((int)(fIndex[j]*100+0.5))/100.0;
      QString editable;
      if (bDecimalIndex)
        editable = QString("%1, %2, %3").arg(fIndex[0]).arg(fIndex[1]).arg(fIndex[2]);
      else
      {
        int nIndex[3];
        layer->RASToOriginalIndex( ras, nIndex );
        editable = QString("%1, %2, %3").arg(nIndex[0]).arg(nIndex[1]).arg(nIndex[2]);
      }
      QString strg = QString("%1 \t[%2]").arg(valueStrg).arg(editable);
      QString labelStrg;
      if (layer->IsTypeOf("PLabel"))
      {
        labelStrg = ((LayerPLabel*)layer)->GetLabelName(m_dRAS);
      }
      else
      {
        labelStrg = layer->GetLabelName( dvalue );
      }
      if (!labelStrg.isEmpty())
      {
        strg += "  " + labelStrg;
      }
      item->setText(1, strg);
      item->setToolTip(1, strg);
      map.clear();
      map["Type"] = "MRI";
      map["EditableText"] = editable;
      map["Object"] = QVariant::fromValue((QObject*)layer);
      item->setData(1, Qt::UserRole, map);
    }
  }

  for (int i = 0; i < lc_surf->GetNumberOfLayers(); i++)
  {
    LayerSurface* surf = (LayerSurface*)lc_surf->GetLayer(i);
    if ( surf->GetProperty()->GetShowInfo() )
    {
      QTreeWidgetItem* item = new QTreeWidgetItem(this);
      item->setText(0, surf->GetName());
      int nVertex = -1;
      bool bMappingVertex = (surf->IsInflated() && surf->GetSourceSurface()->IsSurfaceLoaded(FSSurface::SurfaceWhite));
      if (bMappingVertex)
        nVertex = (m_bForCursor ? surf->GetCurrentVertex() : surf->GetMouseVertex());

      double sf_pos[3];
      if (bMappingVertex && nVertex >= 0)
        surf->GetSourceSurface()->GetSurfaceRASAtVertex(nVertex, sf_pos, FSSurface::SurfaceWhite);
      else
      {
        surf->GetSurfaceRASAtTarget( m_dRAS, sf_pos );
      }
      QString editable = QString("%1, %2, %3")
          .arg(sf_pos[0], 0, 'f', 2)
          .arg(sf_pos[1], 0, 'f', 2)
          .arg(sf_pos[2], 0, 'f', 2);
      item->setText(1, QString("SurfaceRAS\t[%1]").arg(editable));
      map.clear();
      map["Type"] = "SurfaceRAS";
      map["EditableText"] = editable;
      map["Object"] = QVariant::fromValue((QObject*)surf);
      item->setData(1, Qt::UserRole, map);

      if (nVertex < 0)
        nVertex = surf->GetVertexIndexAtTarget( m_dRAS, NULL );
      if ( nVertex >= 0 )
      {
        if (bMappingVertex)
          surf->GetSourceSurface()->GetSurfaceRASAtVertex(nVertex, sf_pos, FSSurface::SurfaceWhite);
        else
          surf->GetSurfaceRASAtVertex( nVertex, sf_pos );
        QTreeWidgetItem* item = new QTreeWidgetItem(this);
        item->setText(1, QString("Vertex \t%1  [%2, %3, %4]")
                      .arg(nVertex)
                      .arg(sf_pos[0], 0, 'f', 2)
            .arg(sf_pos[1], 0, 'f', 2)
            .arg(sf_pos[2], 0, 'f', 2));
        map.clear();
        map["Type"] = "SurfaceVertex";
        map["EditableText"] = QString::number(nVertex);
        map["Object"] = QVariant::fromValue((QObject*)surf);
        item->setData(1, Qt::UserRole, map);

        double vec[3];
        if (m_bShowSurfaceNormal)
        {
          surf->GetNormalAtVertex( nVertex, vec );
          item = new QTreeWidgetItem(this);
          item->setText(1, QString("Normal \t[%1, %2, %3]")
                        .arg(vec[0], 0, 'f', 2)
              .arg(vec[1], 0, 'f', 2)
              .arg(vec[2], 0, 'f', 2));
        }

        if ( surf->GetActiveVector() >= 0 )
        {
          surf->GetVectorAtVertex( nVertex, vec );
          item = new QTreeWidgetItem(this);
          item->setText(1, QString("Vector \t[%1, %2, %3]")
                        .arg(vec[0], 0, 'f', 2)
              .arg(vec[1], 0, 'f', 2)
              .arg(vec[2], 0, 'f', 2));
        }

        if ( surf->HasCurvature() && m_bShowSurfaceCurvature)
        {
          item = new QTreeWidgetItem(this);
          item->setText(1, QString("Curvature \t%1").arg(MyUtils::RealToNumber(surf->GetCurvatureValue(nVertex), nPrecision)));
        }

        int nOverlays = surf->GetNumberOfOverlays();
        for ( int i = 0; i < nOverlays; i++ )
        {
          SurfaceOverlay* overlay = surf->GetOverlay( i );
          item = new QTreeWidgetItem(this);
          item->setText(1, QString("%1 \t%2").arg(overlay->GetName()).arg(MyUtils::RealToNumber(overlay->GetDataAtVertex( nVertex ), nPrecision)));
          item->setToolTip(1, item->text(1));
        }

        int nAnnotations = surf->GetNumberOfAnnotations();
        for ( int i = 0; i < nAnnotations; i++ )
        {
          SurfaceAnnotation* annot = surf->GetAnnotation( i );
          item = new QTreeWidgetItem(this);
          item->setText(1, QString("%1 \t%2").arg(annot->GetName()).arg(annot->GetAnnotationNameAtVertex( nVertex )));
          item->setToolTip(1, item->text(1));
        }

        int nPath = surf->FindPathAt(nVertex);
        if (nPath >= 0)
        {
            SurfacePath* path = surf->GetMadePath(nPath);
            if (path)
            {
                item = new QTreeWidgetItem(this);
                item->setText(1, QString("Path %1 \t%2 mm").arg(nPath).arg(path->GetLength(), 0, 'f', 2));
            }
        }
      }
      else
      {
        QTreeWidgetItem* item = new QTreeWidgetItem(this);
        item->setText(1, "Vertex \tN/A");
        map.clear();
        map["Type"] = "SurfaceVertex";
        map["EditableText"] = "N/A";
        map["Object"] = QVariant::fromValue((QObject*)surf);
        item->setData(1, Qt::UserRole, map);
      }
    }
  }
}

void InfoTreeWidget::OnItemClicked(QTreeWidgetItem *item, int column)
{
  if ( !item || item != m_itemEdited)
  {
    m_editor->hide();
    m_itemEdited = item;
    return;
  }
  QVariantMap map = item->data(column, Qt::UserRole).toMap();
  if (map.contains("EditableText"))
  {
    m_editor->setText(map["EditableText"].toString());
    QRect rect = this->visualRect(indexFromItem(item, column));
    m_editor->show();
    m_editor->move(rect.topLeft() + QPoint(0,header()->height()));
    m_editor->resize(rect.size());
    m_itemEdited = item;
  }
  else
  {
    m_editor->hide();
  }
}

void InfoTreeWidget::OnEditFinished()
{
  if (!m_itemEdited)
  {
    return;
  }

  QStringList list = m_editor->text().trimmed().split(",", MD_SkipEmptyParts);
  QVariantMap map = m_itemEdited->data(1, Qt::UserRole).toMap();
  if ( list.size() < 3)
  {
    list = m_editor->text().trimmed().split(" ", MD_SkipEmptyParts);
  }
  QString type = map["Type"].toString();
  double ras[3];
  bool bSuccess = false;
  QObject* layer = map["Object"].value<QObject*>();
  LayerSurface* surf = NULL;
  if ( type == "SurfaceVertex")
  {
    bool bOK;
    int nVertex = list[0].toInt(&bOK);
    if (bOK && qobject_cast<LayerSurface*>(layer)->GetTargetAtVertex(nVertex, ras))
    {
      bSuccess = true;
      surf = qobject_cast<LayerSurface*>(layer);
      emit VertexChangeTriggered(nVertex);
    }
    else
    {
      std::cerr << "Error: Invalid input";
    }
  }
  else
  {
    if ( list.size() < 3 )
    {
      std::cerr << "Error: Need to enter 3 numbers.";
    }
    else
    {
      bool bOK;
      ras[0] = list[0].toDouble(&bOK);
      ras[1] = list[1].toDouble(&bOK);
      ras[2] = list[2].toDouble(&bOK);
      if (bOK)
      {
        if (type == "RAS")
        {
          LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetActiveLayer();
          if ( mri )
          {
            mri->RASToTarget( ras, ras );
          }
        }
        else if (type == "TkRegRAS")
        {
          LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetActiveLayer();
          if ( mri )
          {
            mri->TkRegToNativeRAS(ras, ras);
            mri->RASToTarget( ras, ras );
          }
        }
        else if (type == "Talairach")
        {
          LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetActiveLayer();
          if ( mri )
          {
            FSVolume* vol = mri->GetSourceVolume();
            vol->TalairachToRAS(ras, ras);
            mri->RASToTarget( ras, ras );
          }
        }
        else if (type == "MRI")
        {
          LayerMRI* mri = qobject_cast<LayerMRI*>(layer);
          //    int nv[3] = {(int)ras[0], (int)ras[1], (int)ras[2]};
          mri->OriginalVoxelToRAS( ras, ras );
          mri->RASToTarget( ras, ras );
        }
        else if (type == "SurfaceRAS")
        {
          qobject_cast<LayerSurface*>(layer)->GetTargetAtSurfaceRAS( ras, ras );
        }
        bSuccess = true;
      }
      else
      {
        std::cerr << "Error: Invalid input";
      }
    }
  }
  if (bSuccess)
  {
    m_editor->hide();
    if (surf)
    {
      RenderView3D* view = qobject_cast<RenderView3D*>(MainWindow::GetMainWindow()->GetRenderView(3));
      if (surf->IsInflated())
        view->MapInflatedCoords(surf, ras, ras, false);
      else
        view->MapToInflatedCoords(ras);
    }
    emit RASChangeTriggered(ras[0], ras[1], ras[2]);
  }
}

void InfoTreeWidget::keyPressEvent(QKeyEvent *event)
{
  if (event->key() == Qt::Key_Escape)
  {
    m_editor->hide();
  }

  QTreeWidget::keyPressEvent(event);
}

void InfoTreeWidget::mousePressEvent(QMouseEvent *event)
{
  if (!itemAt(event->pos()))
  {
    m_editor->hide();
  }

  QTreeWidget::mousePressEvent(event);
}

void InfoTreeWidget::UpdateTrackVolumeAnnotation(Layer *layer, const QVariantMap &info)
{
  for (int i = 0; i < this->topLevelItemCount(); i++)
  {
    QTreeWidgetItem* item = this->topLevelItem(i);
    if (item)
    {
      QVariantMap map = item->data(1, Qt::UserRole).toMap();
      if (map.contains("Object") && map["Object"].value<QObject*>() == layer)
      {
        item->setText(1, QString("%1 \t%2").arg(info["label"].toInt()).arg(info["name"].toString()));
      }
    }
  }
}

void InfoTreeWidget::contextMenuEvent(QContextMenuEvent * e)
{
  MainWindow* mainWnd = MainWindow::GetMainWindow();
  QList<Layer*> layers = mainWnd->GetLayerCollection( "MRI" )->GetLayers();
  QList<Layer*> surfs = mainWnd->GetLayerCollection( "Surface" )->GetLayers();
  layers += surfs;

  if ( layers.isEmpty())
    return;

  QMenu* menu = new QMenu;
  if (!layers.isEmpty())
  {
    QAction* act = new QAction("Show TkReg RAS", this);
    act->setCheckable(true);
    act->setChecked(m_bShowTkRegRAS);
    connect(act, SIGNAL(toggled(bool)), this, SLOT(OnToggleShowTkRegRAS(bool)));
    menu->addAction(act);
    menu->addSeparator();
  }
  foreach (Layer* layer, layers)
  {
    QAction* act = new QAction(layer->GetName(), this);
    act->setCheckable(true);
    act->setChecked(layer->GetProperty()->GetShowInfo());
    act->setData(QVariant::fromValue(qobject_cast<QObject*>(layer)));
    connect(act, SIGNAL(toggled(bool)), this, SLOT(OnToggleShowInfo(bool)));
    menu->addAction(act);
  }
  if (!surfs.isEmpty())
  {
    menu->addSeparator();
    QAction* act = new QAction("Show Surface Curvature", this);
    act->setCheckable(true);
    act->setChecked(m_bShowSurfaceCurvature);
    connect(act, SIGNAL(toggled(bool)), this, SLOT(OnToggleSurfaceCurvature(bool)));
    menu->addAction(act);
    act = new QAction("Show Surface Normal", this);
    act->setCheckable(true);
    act->setChecked(m_bShowSurfaceNormal);
    connect(act, SIGNAL(toggled(bool)), this, SLOT(OnToggleSurfaceNormal(bool)));
    menu->addAction(act);
  }

  QAction* act = new QAction("Show Highlighted Volumes Only", this);
  act->setCheckable(true);
  act->setChecked(m_bShowSelectedLayersOnly);
  connect(act, SIGNAL(toggled(bool)), SLOT(SetShowSelectedLayersOnly(bool)));
  menu->addSeparator();
  menu->addAction(act);
  menu->exec(e->globalPos());
}

void InfoTreeWidget::OnToggleShowInfo(bool bShow)
{
  QAction* act = qobject_cast<QAction*>(sender());
  if (act)
  {
    Layer* layer = qobject_cast<Layer*>(act->data().value<QObject*>());
    if (layer)
      layer->GetProperty()->SetShowInfo(bShow);
  }
}

void InfoTreeWidget::OnToggleSurfaceCurvature(bool show)
{
  m_bShowSurfaceCurvature = show;
  UpdateAll();
}

void InfoTreeWidget::OnToggleSurfaceNormal(bool show)
{
  m_bShowSurfaceNormal = show;
  UpdateAll();
}

void InfoTreeWidget::OnToggleShowTkRegRAS(bool bShow)
{
  m_bShowTkRegRAS = bShow;
  UpdateAll();
}
