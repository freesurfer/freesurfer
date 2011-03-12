#include "PanelROI.h"
#include "ui_PanelROI.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QToolBar>
#include "LayerROI.h"
#include "LayerPropertyROI.h"
#include "MyUtils.h"

PanelROI::PanelROI(QWidget *parent) :
    PanelLayer(parent),
    ui(new Ui::PanelROI)
{
    ui->setupUi(this);
    MainWindow* mainwnd = MainWindow::GetMainWindow();
    if (mainwnd)
    {
      ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionNewROI);
      ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionLoadROI);
      ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionCloseROI);
      ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionSaveROI);
      ui->toolbar->insertSeparator(ui->actionMoveLayerUp);
    }

    LayerCollection* lc = mainwnd->GetLayerCollection("ROI");
    PanelLayer::InitializeLayerList( ui->treeWidgetLayers, lc );
}

PanelROI::~PanelROI()
{
    delete ui;
}

void PanelROI::ConnectLayer( Layer* layer_in )
{
  PanelLayer::ConnectLayer( layer_in );

  LayerROI* layer = qobject_cast<LayerROI*>(layer_in);
  if ( !layer )
    return;

  LayerPropertyROI* p = layer->GetProperty();
  connect( p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( ui->doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetOpacity(double)) );
  connect( ui->colorPickerColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetColor(QColor)) );
}

void PanelROI::DoIdle()
{
  // update action status
  BlockAllSignals( true );
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  LayerROI* layer = NULL;
  if ( item )
      layer = qobject_cast<LayerROI*>( item->data(0, Qt::UserRole).value<QObject*>() );
  int nItemIndex = ui->treeWidgetLayers->indexOfTopLevelItem(item);
  ui->actionMoveLayerUp->setEnabled( item && !layer->IsLocked() && ui->treeWidgetLayers->topLevelItemCount() > 1 &&
                                     nItemIndex != 0 );
  ui->actionMoveLayerDown->setEnabled( item && !layer->IsLocked() && ui->treeWidgetLayers->topLevelItemCount() > 1 &&
                                       nItemIndex < ui->treeWidgetLayers->topLevelItemCount()-1 );

  BlockAllSignals( false );
}

void PanelROI::OnSliderOpacity(int nVal)
{
    LayerROI* layer = GetCurrentLayer<LayerROI*>();
    if ( layer )
      layer->GetProperty()->SetOpacity( nVal / 100.0 );
}

void PanelROI::DoUpdateWidgets()
{
    BlockAllSignals( true );
    for ( int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++ )
    {
        QTreeWidgetItem* item = ui->treeWidgetLayers->topLevelItem( i );
        Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole).value<QObject*>() );
        if ( layer )
        {
            item->setCheckState( 0, (layer->IsVisible() ? Qt::Checked : Qt::Unchecked) );
        }
    }

    LayerROI* layer = GetCurrentLayer<LayerROI*>();
    for ( int i = 0; i < this->allWidgets.size(); i++ )
    {
        if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar )
            allWidgets[i]->setEnabled(layer);
    }

    ui->lineEditFileName->clear();
    if ( layer )
    {
      ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
      ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );

      double* rgb = layer->GetProperty()->GetColor();
      ui->colorPickerColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );

      ui->lineEditFileName->setText( MyUtils::Win32PathProof(layer->GetFileName()) );
      ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );
    }
    ui->doubleSpinBoxOpacity->setEnabled( layer );
    ui->colorPickerColor->setEnabled( layer );
    ui->sliderOpacity->setEnabled( layer );
    ui->labelColor->setEnabled( layer );
    ui->labelFileName->setEnabled( layer );
    ui->lineEditFileName->setEnabled( layer );
    ui->labelOpacity->setEnabled( layer );

    BlockAllSignals( false );
}
