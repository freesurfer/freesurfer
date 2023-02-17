#include "DialogSaveAllVolumes.h"
#include "ui_DialogSaveAllVolumes.h"
#include <QFileDialog>
#include "LayerMRI.h"
#include "MainWindow.h"
#include <QFileInfo>

ThreadSaveAll::ThreadSaveAll(QObject *parent) :
  QThread(parent)
{
}

void ThreadSaveAll::SaveAllVolumes(const QList<LayerMRI *> &layers)
{
  m_listMRIs = layers;
  start();
}

void ThreadSaveAll::run()
{
  for (int i = 0; i < m_listMRIs.size(); i++)
  {
    QString error_msg;
    if (!m_listMRIs[i]->SaveVolume())
      error_msg = tr("Failed to save volume to %").arg(m_listMRIs[i]->GetFileName());
    emit Progress((i+1)*100/m_listMRIs.size());
    emit Saved(m_listMRIs[i], error_msg);
  }
}

DialogSaveAllVolumes::DialogSaveAllVolumes(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSaveAllVolumes)
{
  ui->setupUi(this);
  ui->progressBar->hide();
  ui->labelErrorMessage->hide();

  connect(ui->toolButtonOpen, SIGNAL(clicked(bool)), SLOT(OnButtonOpen()));
  connect(ui->pushButtonSave, SIGNAL(clicked(bool)), SLOT(OnButtonSave()));
  m_thread = new ThreadSaveAll( this );
  connect(m_thread, SIGNAL(Progress(int)), ui->progressBar, SLOT(setValue(int)));
  connect(m_thread, SIGNAL(finished()), SLOT(OnFinished()));
}

DialogSaveAllVolumes::~DialogSaveAllVolumes()
{
  delete ui;
}

void DialogSaveAllVolumes::OnButtonOpen()
{
  QString fn = QFileDialog::getExistingDirectory(this, "Select output directory",
                                                 m_strLastDir);
  if ( !fn.isEmpty() )
  {
    ui->lineEditOutputFolder->setText(fn);
    ui->lineEditOutputFolder->setCursorPosition( fn.size() );
    m_strLastDir = fn;
  }
}

void DialogSaveAllVolumes::OnButtonSave()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  LayerMRI* cur_mri = qobject_cast<LayerMRI*>(mainwnd->GetActiveLayer("MRI"));
  if (!cur_mri)
    return;
  QList<Layer*> list = mainwnd->GetLayers("MRI");
  QList<LayerMRI*> layers;
  QString suffix = ui->lineEditSuffix->text().trimmed();
  foreach (Layer* layer, list)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(layer);
    if (mri != cur_mri)
      mri->CopyTransformation(cur_mri);
    mri->SetWriteResampled(!ui->checkBoxNoResample->isChecked());
    mri->SetCropToOriginal(ui->checkBoxCrop->isChecked());
    QString fn = mri->GetFileName();
    QString old_suffix = QFileInfo(fn).suffix();
    if (fn.endsWith("nii.gz", Qt::CaseInsensitive))
      old_suffix = "nii.gz";
    fn.resize(fn.size()-old_suffix.size());
    fn += suffix + "." + old_suffix;
    if (!ui->checkBoxSaveToSourceDir->isChecked())
    {
      QFileInfo fi(QDir(ui->lineEditOutputFolder->text().trimmed()), QFileInfo(fn).fileName());
      fn = fi.absoluteFilePath();
    }
    mri->setProperty("original_filename", mri->GetFileName());
    mri->SetFileName(fn);
    layers << mri;
  }
  m_thread->SaveAllVolumes(layers);
  ui->progressBar->show();
  ui->pushButtonSave->setDisabled(true);
}

void DialogSaveAllVolumes::OnVolumeSaved(LayerMRI *mri, const QString &error_msg)
{
  mri->Restore();
  mri->SetFileName(mri->property("original_filename").toString());
  if (!error_msg.isEmpty())
  {
    ui->labelErrorMessage->show();
    ui->labelErrorMessage->setText(error_msg);
  }
}

void DialogSaveAllVolumes::OnFinished()
{
  ui->pushButtonSave->setEnabled(true);
  if (!ui->labelErrorMessage->isVisible())
    accept();
}
