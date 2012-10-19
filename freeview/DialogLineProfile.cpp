#include "DialogLineProfile.h"
#include "ui_DialogLineProfile.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPointSet.h"
#include "LayerCollection.h"
#include "LayerLineProfile.h"
#include <QFileDialog>
#include <QDebug>

DialogLineProfile::DialogLineProfile(QWidget *parent) :
    QDialog(parent),
    m_lineProfile(NULL),
    ui(new Ui::DialogLineProfile)
{
    ui->setupUi(this);
    ui->labelError->hide();
}

DialogLineProfile::~DialogLineProfile()
{
    delete ui;
}

void DialogLineProfile::OnCompute()
{
  qDebug() << "to compute";
  QVariant v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine1->currentIndex());
  LayerPointSet* layer1 = qobject_cast<LayerPointSet*>(v.value<QObject*>());
  v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine2->currentIndex());
  LayerPointSet* layer2 = qobject_cast<LayerPointSet*>(v.value<QObject*>());
  if (!layer1 || !layer2 || layer1 == layer2)
  {
    ui->labelError->setText(tr("Please select two point sets"));
    ui->labelError->show();
    return;
  }

  ui->labelError->hide();

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  int nViewId = mainwnd->GetActiveViewId();
  LayerCollection* col = mainwnd->GetLayerCollection("Supplement");
  if (!col)
    return;
  if (nViewId > 2)
  {
    ui->labelError->setText(tr("Please select a 2D slice view"));
    ui->labelError->show();
    return;
  }

  QList<Layer*> lineLayers = col->GetLayers("LineProfile");
  if (lineLayers.isEmpty())
  {
    m_lineProfile = new LayerLineProfile(nViewId, NULL, layer1, layer2);
    col->AddLayer(m_lineProfile);
  }
  else
  {
    m_lineProfile = NULL;
    foreach (Layer* layer, lineLayers)
    {
      LayerLineProfile* line = qobject_cast<LayerLineProfile*>(layer);
      if (line->GetPlane() == nViewId)
      {
        m_lineProfile = line;
        break;
      }
    }
    if (!m_lineProfile)
    {
      m_lineProfile = new LayerLineProfile(nViewId, NULL, layer1, layer2);
      col->AddLayer(m_lineProfile);
    }
    else
      m_lineProfile->SetSourceLayers(layer1, layer2);
  }

  m_lineProfile->Solve(GetResolution());
}

double DialogLineProfile::GetResolution()
{
  return ui->lineEditResolution->text().toDouble();
}

int DialogLineProfile::GetNumberOfSamples()
{
  return ui->lineEditSamplePoints->text().toInt();
}

void DialogLineProfile::OnExport()
{
  ui->labelError->hide();
  if (GetNumberOfSamples() < 2)
  {
    ui->labelError->setText("Sample points must be more than 1");
    ui->labelError->show();
    return;
  }

  QString fn = QFileDialog::getSaveFileName(this, "Select File To Save",
               "",
               "CSV files (*.csv);;All Files (*)");
  if (!fn.isEmpty())
  {
    MainWindow* mainwnd = MainWindow::GetMainWindow();
    LayerMRI* mri = qobject_cast<LayerMRI*>(mainwnd->GetActiveLayer("MRI"));
    if (mri)
    {
      if (!m_lineProfile->Export(fn, mri, GetNumberOfSamples()))
      {
        ui->labelError->setText("Failed to export");
        ui->labelError->show();
      }
    }
  }
}

void DialogLineProfile::OnComboIsoLine(int sel)
{

}

void DialogLineProfile::UpdatePointSetList()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QList<Layer*> layers = mainwnd->GetLayers("PointSet");
  QVariant v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine1->currentIndex());
  Layer* layer1 = qobject_cast<Layer*>(v.value<QObject*>());
  v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine2->currentIndex());
  Layer* layer2 = qobject_cast<Layer*>(v.value<QObject*>());
  ui->comboBoxIsoLine1->clear();
  ui->comboBoxIsoLine2->clear();
  foreach (Layer* layer, layers)
  {
    ui->comboBoxIsoLine1->addItem(layer->GetName(), qVariantFromValue<QObject*>(layer));
    ui->comboBoxIsoLine2->addItem(layer->GetName(), qVariantFromValue<QObject*>(layer));
  }

  if (layers.contains(layer1))
  {
    ui->comboBoxIsoLine1->setCurrentIndex(layers.indexOf(layer1));
  }
  else if (!layers.isEmpty())
    ui->comboBoxIsoLine1->setCurrentIndex(0);

  if (layers.contains(layer2))
  {
    ui->comboBoxIsoLine2->setCurrentIndex(layers.indexOf(layer2));
  }
  else if (layers.size() > 1)
    ui->comboBoxIsoLine2->setCurrentIndex(1);
}
