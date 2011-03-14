#include "DialogNewROI.h"
#include "ui_DialogNewROI.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "MainWindow.h"
#include <QMessageBox>

DialogNewROI::DialogNewROI(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogNewROI)
{
    ui->setupUi(this);

    LayerCollection* col_mri = MainWindow::GetMainWindow()->GetLayerCollection("MRI");
    QList<Layer*> layers = col_mri->GetLayers();
    int nSel = 0;
    for ( int i = 0; i < layers.size(); i++ )
    {
        ui->comboBoxTemplate->addItem(layers[i]->GetName(), QVariant::fromValue((QObject*)layers[i]));
        if ( layers[i] == col_mri->GetActiveLayer() )
            nSel = i;
    }
    ui->comboBoxTemplate->setCurrentIndex(nSel);
    ui->lineEditName->setFocus();
}

DialogNewROI::~DialogNewROI()
{
    delete ui;
}

QString DialogNewROI::GetROIName()
{
  return ui->lineEditName->text().trimmed();
}

void DialogNewROI::SetROIName( const QString& name )
{
  ui->lineEditName->setText( name );
}

LayerMRI* DialogNewROI::GetTemplate()
{
  return qobject_cast<LayerMRI*>(
          ui->comboBoxTemplate->itemData(ui->comboBoxTemplate->currentIndex()).value<QObject*>());
}

void DialogNewROI::OnOK()
{
    if ( GetROIName().isEmpty())
    {
        QMessageBox::warning( this, "Error", "ROI name can not be empty." );
        return;
    }
    accept();
}
