#include "DialogWriteMovieFrames.h"
#include "ui_DialogWriteMovieFrames.h"
#include "MainWindow.h"
#include "RenderView.h"
#include "RenderView2D.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QCloseEvent>
#include "Layer.h"

DialogWriteMovieFrames::DialogWriteMovieFrames(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogWriteMovieFrames)
{
    ui->setupUi(this);
    connect(MainWindow::GetMainWindow(), SIGNAL(MainViewChanged(int)),
            this, SLOT(UpdateUI()));
    m_timer.setInterval(1000);
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(OnTimeOut()));
}

DialogWriteMovieFrames::~DialogWriteMovieFrames()
{
    delete ui;
}

void DialogWriteMovieFrames::showEvent(QShowEvent *e)
{
    UpdateUI();

    QDialog::showEvent(e);
}

void DialogWriteMovieFrames::closeEvent(QCloseEvent *e)
{
    if (m_timer.isActive())
    {
        if ( QMessageBox::question(this, "Close", "Still writing movie frames. Do you want to abort?") != QMessageBox::Yes)
        {
            e->ignore();
            return;
        }
    }

    QDialog::closeEvent(e);
}

void DialogWriteMovieFrames::UpdateUI()
{
    m_view = MainWindow::GetMainWindow()->GetMainView();
    bool b3D = m_b3D = (MainWindow::GetMainWindow()->GetMainViewId() == 3);
    ui->labelAngleStep->setVisible(b3D);
    ui->doubleSpinBoxAngleStep->setVisible(b3D);
    ui->labelSliceStartNumber->setVisible(!b3D);
    ui->labelSliceStep->setVisible(!b3D);
    ui->spinBoxSliceStart->setVisible(!b3D);
    ui->spinBoxSliceStep->setVisible(!b3D);
    ui->pushButtonAbort->setEnabled(m_timer.isActive());
    ui->pushButtonClose->setEnabled(!m_timer.isActive());
    ui->pushButtonWrite->setEnabled(!m_timer.isActive());
}

void DialogWriteMovieFrames::OnOpen()
{
    QString strg = QFileDialog::getExistingDirectory(this, "Select Output Directory");
    if (!strg.isEmpty())
    {
        ui->lineEditOutputLocation->setText(strg);
        ui->lineEditOutputLocation->setCursorPosition(ui->lineEditOutputLocation->text().length());
    }
}

void DialogWriteMovieFrames::OnWrite()
{
    m_strOutputDir = ui->lineEditOutputLocation->text().trimmed();
    if (m_strOutputDir.isEmpty())
    {
        QMessageBox::warning(this, "Error", "Please enter output directory to save movie frames.");
        return;
    }
    QDir dir(m_strOutputDir);
    if (!dir.exists())
    {
        QMessageBox::warning(this, "Error", "Output directory does not exist.");
        return;
    }
    if (ui->doubleSpinBoxAngleStep->value() == 0)
    {
        QMessageBox::warning(this, "Error", "Angel step size can not be 0.");
        return;
    }
    if (m_strOutputDir[m_strOutputDir.length()-1] != '/' &&
        m_strOutputDir[m_strOutputDir.length()-1] != '\\')
        m_strOutputDir += QDir::separator();
    m_dAngleStepSize = ui->doubleSpinBoxAngleStep->value();
    m_nStartSlice = ui->spinBoxSliceStart->value();
    m_nSliceStepSize = ui->spinBoxSliceStep->value();
    m_nStepCount = 0;
    m_nTotalSteps = 1;
    if (m_b3D)
        m_nTotalSteps = (int)qAbs(360.0/m_dAngleStepSize+0.5);
    else
    {
        MainWindow* mwnd = MainWindow::GetMainWindow();
        int nView = mwnd->GetMainViewId();
        Layer* layer = mwnd->GetActiveLayer( "MRI" );
        if (!layer)
            layer = mwnd->GetActiveLayer("Surface");
        if (layer)
        {
            m_nTotalSteps = (int)((layer->GetWorldSize())[nView] / (layer->GetWorldVoxelSize())[nView]+0.5);
        }
        ((RenderView2D*)m_view)->SetSliceNumber( m_nStartSlice );
    }
    if (m_nTotalSteps < 1)
        m_nTotalSteps = 1;
    m_timer.start();
    UpdateUI();
    emit Started();
}

void DialogWriteMovieFrames::OnAbort()
{
    m_timer.stop();
    UpdateUI();
    emit Stopped();
}

void DialogWriteMovieFrames::OnTimeOut()
{
    QString fn;
    SettingsScreenshot settings = MainWindow::GetMainWindow()->GetScreenShotSettings();
    if ( m_b3D ) // 3D view
    {
        fn = QString("%1frame%2.%3").arg(m_strOutputDir)
             .arg(m_nStepCount, 3, 10, QChar('0'))
             .arg(ui->comboBoxExtension->currentText());
      m_view->SaveScreenShot( fn, settings.AntiAliasing, settings.Magnification );
      CameraOperations ops;
      ops << CameraOperation("azimuth", m_dAngleStepSize);
      m_view->SetCameraOperations(ops);
    }
    else
    {
        int nSlice = m_nStartSlice+m_nSliceStepSize*m_nStepCount;
        fn = QString("%1frame%2.%3").arg(m_strOutputDir)
             .arg(nSlice, 3, 10, QChar('0'))
             .arg(ui->comboBoxExtension->currentText());
      m_view->SaveScreenShot( fn,settings.AntiAliasing, settings.Magnification );
      if ( !((RenderView2D*)m_view)->SetSliceNumber( nSlice + m_nSliceStepSize ) )
        OnAbort();
  }

  m_view->RequestRedraw(true);
  m_nStepCount++;
  if (m_nStepCount >= m_nTotalSteps)
      OnAbort();
  else
      emit Progress(m_nStepCount*100/m_nTotalSteps);
}
