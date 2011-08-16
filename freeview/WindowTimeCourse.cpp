#include "WindowTimeCourse.h"
#include "ui_WindowTimeCourse.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include "FSVolume.h"
#include <QSettings>

WindowTimeCourse::WindowTimeCourse(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WindowTimeCourse)
{
    ui->setupUi(this);
    this->setWindowFlags(Qt::Tool);
    this->setWindowTitle("Time Course");
    connect(ui->widgetPlot, SIGNAL(FrameChanged(int)), this, SLOT(OnFrameChanged(int)));

    QSettings s;
    QVariant v = s.value("WindowTimeCourse/Geomerty");
    if (v.isValid())
      this->restoreGeometry(v.toByteArray());
    ui->checkBoxAutoScale->setChecked(s.value("WindowTimeCourse/AutoScale", true).toBool());
}

WindowTimeCourse::~WindowTimeCourse()
{
  QSettings s;
  s.setValue("WindowTimeCourse/Geomerty", this->saveGeometry());
  s.setValue("WindowTimeCourse/AutoScale", ui->checkBoxAutoScale->isChecked());
  delete ui;
}

void WindowTimeCourse::UpdateData()
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer("MRI"));
  if (layer && layer->GetNumberOfFrames() > 1)
  {
    double ras[3];
    int n[3];
    MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetSlicePosition(ras);
    layer->RemapPositionToRealRAS(ras, ras);
    layer->RASToOriginalIndex(ras, n);
    QList<double> data;
    for (int i = 0; i < layer->GetNumberOfFrames(); i++)
      data <<  layer->GetVoxelValueByOriginalIndex(n[0], n[1], n[2], i);
    FSVolume* vol = layer->GetSourceVolume();
    ui->widgetPlot->SetTimeCourseData(data, vol->GetMinValue(), vol->GetMaxValue(), layer->GetTR());
    ui->widgetPlot->SetCurrentFrame(layer->GetActiveFrame());
  }
}

void WindowTimeCourse::OnFrameChanged(int frame)
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer("MRI"));
  if (layer && frame != layer->GetActiveFrame() && frame < layer->GetNumberOfFrames())
  {
    layer->SetActiveFrame(frame);
  }
}

void WindowTimeCourse::SetCurrentFrame(int n)
{
  ui->widgetPlot->SetCurrentFrame(n);
}
