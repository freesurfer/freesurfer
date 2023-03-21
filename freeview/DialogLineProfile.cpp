#include "DialogLineProfile.h"
#include "ui_DialogLineProfile.h"
#include "MainWindow.h"
#include "LayerLineProfile.h"
#include "LayerMRI.h"
#include "LayerPointSet.h"
#include "LayerCollection.h"
#include "RenderView.h"
#include <QFileDialog>
#include <QDebug>
#include "LayerPropertyLineProfile.h"
#include "DialogSelectSplines.h"

DialogLineProfile::DialogLineProfile(QWidget *parent) :
  QDialog(parent),
  m_lineProfile(NULL),
  ui(new Ui::DialogLineProfile)
{
  ui->setupUi(this);
  ui->labelError->hide();
  UpdateWidgets();
}

DialogLineProfile::~DialogLineProfile()
{
  delete ui;
}

bool DialogLineProfile::Validate(LayerPointSet*& spline0_out, LayerPointSet*& spline1_out)
{
  QVariant v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine1->currentIndex());
  LayerPointSet* layer1 = qobject_cast<LayerPointSet*>(v.value<QObject*>());
  v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine2->currentIndex());
  LayerPointSet* layer2 = qobject_cast<LayerPointSet*>(v.value<QObject*>());
  if (!layer1 || !layer2 || layer1 == layer2)
  {
    ui->labelError->setText(tr("Please select two point sets"));
    ui->labelError->show();
    return false;
  }
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  int nViewId = mainwnd->GetActiveViewId();
  if (nViewId > 2)
  {
    ui->labelError->setText(tr("Please select a 2D slice view"));
    ui->labelError->show();
    return false;
  }
  spline0_out = layer1;
  spline1_out = layer2;
  return true;
}

void DialogLineProfile::OnCompute()
{
  ui->labelError->hide();
  LayerPointSet* layer1 = NULL, *layer2 = NULL;
  if (!Validate(layer1, layer2))
    return;

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  int nViewId = mainwnd->GetActiveViewId();

  LayerCollection* col = mainwnd->GetLayerCollection("Supplement");
  if (!col)
    return;

  QList<Layer*> lineLayers = col->GetLayers("LineProfile");
  if (!lineLayers.isEmpty())
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
    if (m_lineProfile)
    {
      col->RemoveLayer(m_lineProfile);
      m_lineProfile = NULL;
    }
  }
  m_lineProfile = new LayerLineProfile(nViewId, NULL, layer1, layer2);
  double r = ui->lineEditRadius->text().toDouble();
  if (r <= 0)
  {
    double vs[3];
    col->GetWorldVoxelSize(vs);
    r = qMin(vs[0], qMin(vs[1], vs[2]))*0.1;
    if (r <= 0)
      r = 0.1;
  }
  m_lineProfile->GetProperty()->SetRadius(r);
  col->AddLayer(m_lineProfile);

  double dVoxelSize = 1.0;
  LayerMRI* mri = qobject_cast<LayerMRI*>(mainwnd->GetActiveLayer("MRI"));
  if (mri)
  {
    double vs[3];
    mri->GetWorldVoxelSize(vs);
    dVoxelSize = qMin(vs[0], qMin(vs[1], vs[2]));
  }
  m_lineProfile->Solve(GetSpacing(), dVoxelSize, GetResolution(), GetOffset());

  mainwnd->SetMode(RenderView::IM_Navigate);

  UpdateWidgets();
}

void DialogLineProfile::OnClear()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  int nViewId = mainwnd->GetActiveViewId();

  LayerCollection* col = mainwnd->GetLayerCollection("Supplement");
  if (!col)
    return;

  QList<Layer*> lineLayers = col->GetLayers("LineProfile");
  foreach (Layer* layer, lineLayers)
  {
    LayerLineProfile* line = qobject_cast<LayerLineProfile*>(layer);
    if (line->GetPlane() == nViewId)
    {
      col->RemoveLayer(layer);
    }
  }
  m_lineProfile = NULL;
  UpdateWidgets();
}

double DialogLineProfile::GetResolution()
{
  return ui->lineEditResolution->text().toDouble();
}

double DialogLineProfile::GetSpacing()
{
  return ui->lineEditSpacing->text().toDouble();
}

double DialogLineProfile::GetOffset()
{
  return ui->lineEditOffset->text().toDouble();
}

int DialogLineProfile::GetNumberOfSamples()
{
  return ui->lineEditSamplePoints->text().toInt();
}

void DialogLineProfile::UpdateWidgets()
{
  ui->pushButtonExport->setEnabled(m_lineProfile);
  ui->pushButtonSave->setEnabled(m_lineProfile);
  ui->pushButtonExportThickness->setEnabled(m_lineProfile);
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

  QString fn = QFileDialog::getSaveFileName(this, "Export Line Profiles to File",
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

void DialogLineProfile::OnSave()
{
  ui->labelError->hide();
  QString fn = QFileDialog::getSaveFileName(this, "Select File to Save",
                                            "",
                                            "All Files (*)");
  if (!fn.isEmpty())
  {
    LayerLineProfile* lp = m_lineProfile;
    if (!lp)
    {
      LayerPointSet* ptset0, *ptset1;
      if (Validate(ptset0, ptset1))
      {
        MainWindow* mainwnd = MainWindow::GetMainWindow();
        lp = new LayerLineProfile(mainwnd->GetActiveViewId(), NULL, ptset0, ptset1);
      }
    }
    if (lp && !lp->Save(fn))
    {
      ui->labelError->setText("Failed to save splines");
      ui->labelError->show();
    }
    if (lp != m_lineProfile)
      delete lp;
  }
}

void DialogLineProfile::OnLoad()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  LayerMRI* mri = qobject_cast<LayerMRI*>(mainwnd->GetActiveLayer("MRI"));
  if (!mri)
  {
    return;
  }
  ui->labelError->hide();
  QString fn = QFileDialog::getOpenFileName(this, "Select File to Open",
                                            "",
                                            "All Files (*)");
  if (!fn.isEmpty())
  {
    LayerLineProfile* lp = LayerLineProfile::Load(fn, mri);
    if (!lp)
    {
      ui->labelError->setText("Failed to load splines");
      ui->labelError->show();
    }
    else
    {
      LayerCollection* col_wp = mainwnd->GetLayerCollection("PointSet");
      LayerCollection* col_sup = mainwnd->GetLayerCollection("Supplement");
      LayerCollection* col_mri = mainwnd->GetLayerCollection("MRI");
      if ( col_wp->IsEmpty() )
      {
        col_wp->SetWorldOrigin( col_mri->GetWorldOrigin() );
        col_wp->SetWorldSize( col_mri->GetWorldSize() );
        col_wp->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
        col_wp->SetSlicePosition( col_mri->GetSlicePosition() );
      }
      col_wp->AddLayer(lp->GetSpline0());
      col_wp->AddLayer(lp->GetSpline1());
      col_sup->AddLayer(lp);
      this->m_lineProfile = lp;
      ui->lineEditResolution->setText(QString::number(lp->GetResultion()));
      ui->lineEditSamplePoints->setText(QString::number(lp->GetNumberOfSamples()));
      ui->lineEditSpacing->setText(QString::number(lp->GetSpacing()));
      //   m_lineProfile->Solve(GetResolution());
    }
  }
}

void DialogLineProfile::OnComboIsoLine(int sel)
{
  Q_UNUSED(sel);
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
    ui->comboBoxIsoLine1->addItem(layer->GetName(), QVariant::fromValue<QObject*>(layer));
    ui->comboBoxIsoLine2->addItem(layer->GetName(), QVariant::fromValue<QObject*>(layer));
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

void DialogLineProfile::OnColorPicker(const QColor &color)
{
  if (m_lineProfile)
    m_lineProfile->GetProperty()->SetColor(color);
}

void DialogLineProfile::OnEditRadius(const QString &strg)
{
  if (m_lineProfile)
  {
    bool ok = false;
    double val = strg.toDouble(&ok);
    if (ok && val > 0)
      m_lineProfile->GetProperty()->SetRadius(val);
  }
}

void DialogLineProfile::OnSliderOpacity(int val)
{
  if (m_lineProfile)
  {
    m_lineProfile->GetProperty()->SetOpacity(val/100.0);
  }
}

void DialogLineProfile::OnLineProfileIdPicked(LayerLineProfile *lp, int nId)
{
  if (lp == m_lineProfile)
  {
    ui->labelActiveId->setText(QString::number(nId));
  }
}

void DialogLineProfile::OnExportThickness()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QList<Layer*> layers = mainwnd->GetLayers("PointSet");
  QVariant v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine1->currentIndex());
  LayerPointSet* layer1 = qobject_cast<LayerPointSet*>(v.value<QObject*>());
  v = ui->comboBoxIsoLine1->itemData(ui->comboBoxIsoLine2->currentIndex());
  LayerPointSet* layer2 = qobject_cast<LayerPointSet*>(v.value<QObject*>());
  layers.removeAll(layer1);
  layers.removeAll(layer2);
  DialogSelectSplines dlg(this);
  dlg.SetPointSets(layers);
  if (dlg.exec() == QDialog::Accepted)
  {
    ui->labelError->hide();
    if (layers.isEmpty())
    {
      ui->labelError->setText("No layer point sets to choose from.");
      ui->labelError->show();
      return;
    }

    QString fn = QFileDialog::getSaveFileName(this, "Export Thickness to File",
                                              "",
                                              "CSV files (*.csv);;All Files (*)");
    if (!fn.isEmpty())
    {
      if (!m_lineProfile->ExportThickness(fn, dlg.GetSelectedPointSets(), GetNumberOfSamples()))
      {
        ui->labelError->setText("Failed to export");
        ui->labelError->show();
      }
    }
  }
}
