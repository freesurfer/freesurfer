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
#include "PanelLayer.h"
#include "Layer.h"
#include "LayerProperty.h"
#include <QTimer>
#include <QApplication>
#include <QTreeWidget>
#include <QAction>
#include <QLineEdit>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QDebug>
#include "LayerTreeWidget.h"
#include "MainWindow.h"

PanelLayer::PanelLayer(const QString& layerType, QWidget *parent) :
  QScrollArea(parent),
  UIUpdateHelper(),
  m_bToUpdate( false ),
  m_currentLayer( NULL ),
  treeWidgetLayers( NULL ),
  m_layerType(layerType)
{
  m_layerCollection = MainWindow::GetMainWindow()->GetLayerCollection(layerType);
  connect(m_layerCollection, SIGNAL(destroyed()), this, SLOT(ResetLayerCollection()));
  QTimer* idleTimer = new QTimer(this);
  connect( idleTimer, SIGNAL(timeout()), this, SLOT(OnIdle()), Qt::QueuedConnection);
  idleTimer->start(500);
  QTimer* updateTimer = new QTimer(this);
  connect( updateTimer, SIGNAL(timeout()), this, SLOT(OnUpdate()), Qt::QueuedConnection);
  updateTimer->start(50);
}

void PanelLayer::UpdateWidgets()
{
  m_bToUpdate = true;
}

void PanelLayer::OnIdle()
{
  DoIdle();
}

void PanelLayer::OnUpdate()
{
  if ( m_bToUpdate )
  {
    DoUpdateWidgets();
    m_bToUpdate = false;
  }
}

void PanelLayer::SetCurrentLayer(Layer *layer)
{
  m_currentLayer = layer;
  UpdateWidgets();
}

void PanelLayer::InitializeLayerTreeWidget( QTreeWidget* treeWidget)
{
  allWidgets = this->widget()->findChildren<QWidget*>();
  allActions = this->findChildren<QAction*>();
  treeWidgetLayers = treeWidget;
  /*
  connect( cl, SIGNAL(LayerAdded(Layer*)), this, SLOT(OnLayerAdded(Layer*)), Qt::UniqueConnection );
  connect( cl, SIGNAL(LayerRemoved(Layer*)), this, SLOT(OnLayerRemoved(Layer*)), Qt::UniqueConnection );
  connect( cl, SIGNAL(LayerMoved(Layer*)), this, SLOT(OnLayerMoved(Layer*)), Qt::UniqueConnection );
  connect( cl, SIGNAL(ActiveLayerChanged(Layer*)), this, SLOT(OnActiveLayerChanged(Layer*)), Qt::UniqueConnection );
  connect( cl, SIGNAL(LayerCycled(Layer*)), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( cl, SIGNAL(LayerNameChanged()), this, SLOT(OnLayerNameChanged()), Qt::UniqueConnection );
  connect( cl, SIGNAL(LayerVisibilityChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( treeWidget, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)),
           this, SLOT(OnCurrentItemChanged(QTreeWidgetItem*)), Qt::UniqueConnection );
  connect( treeWidget, SIGNAL(itemSelectionChanged()), this,
           SLOT(OnItemSelectionChanged()));
  connect( treeWidget, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)),
           this, SLOT(UpdateWidgets()));
  connect( treeWidget, SIGNAL(itemChanged(QTreeWidgetItem*,int)), this, SLOT(OnItemChanged(QTreeWidgetItem*)) );
  connect( treeWidget, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(OnItemDoubleClicked(QModelIndex)));
*/

  //  UpdateWidgets();
}

void PanelLayer::ConnectLayer( Layer* layer )
{
  connect( layer, SIGNAL(Locked(bool)), treeWidgetLayers, SLOT(ForceUpdate()) );
}

void PanelLayer::DisconnectAllLayers()
{
  LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection(m_layerType);
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    for ( int j = 0; j < allWidgets.size(); j++ )
    {
      lc->GetLayer( i )->disconnect( this );
      lc->GetLayer( i )->GetProperty()->disconnect( this );
      allWidgets[j]->disconnect( lc->GetLayer( i ) );
      allWidgets[j]->disconnect( lc->GetLayer( i )->GetProperty() );
    }
  }
}

void PanelLayer::BlockAllSignals( bool bBlock )
{
  for ( int i = 0; i < allWidgets.size(); i++ )
  {
    allWidgets[i]->blockSignals( bBlock );
  }
  for ( int i = 0; i < allActions.size(); i++ )
  {
    allActions[i]->blockSignals( bBlock );
  }
}

