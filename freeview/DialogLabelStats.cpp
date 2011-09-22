#include "DialogLabelStats.h"
#include "ui_DialogLabelStats.h"
#include "MainWindow.h"
#include "LayerMRI.h"

DialogLabelStats::DialogLabelStats(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DialogLabelStats)
{
    ui->setupUi(this);
    setWindowFlags(Qt::Dialog);
}

DialogLabelStats::~DialogLabelStats()
{
    delete ui;
}

void DialogLabelStats::OnSlicePositionChanged()
{
  if (!isVisible())
    return;

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  LayerMRI* mri = (LayerMRI*)mainwnd->GetActiveLayer("MRI");
  if (mri)
  {
    float fLabel, fArea = 0;
    int nCount = 0;
    mri->GetCurrentLabelStats( mainwnd->GetMainViewId(), &fLabel, &nCount, &fArea );
    ui->labelLabel->setText(QString::number((int)fLabel));
    ui->labelCount->setText(QString::number(nCount));
    ui->labelArea->setText(QString("%3 mm2").arg(fArea));
  }
}

void DialogLabelStats::showEvent(QShowEvent *e)
{
  OnSlicePositionChanged();
  QWidget::showEvent(e);
}
