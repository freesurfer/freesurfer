/**
 * @brief Dialog window to apply volume crop
 *
 */
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
#include "DialogCropVolume.h"
#include "ui_DialogCropVolume.h"
#include "MainWindow.h"
#include "RenderView3D.h"
#include "VolumeCropper.h"
#include "LayerMRI.h"
#include <QShowEvent>
#include <QHideEvent>
#include <QSettings>

DialogCropVolume::DialogCropVolume(QWidget *parent, LayerMRI *mri) :
  QDialog(parent),
  ui(new Ui::DialogCropVolume)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint |
                        Qt::WindowCloseButtonHint| Qt::CustomizeWindowHint );
  m_spinRange[0] = ui->spinBoxMinX;
  m_spinRange[1] = ui->spinBoxMaxX;
  m_spinRange[2] = ui->spinBoxMinY;
  m_spinRange[3] = ui->spinBoxMaxY;
  m_spinRange[4] = ui->spinBoxMinZ;
  m_spinRange[5] = ui->spinBoxMaxZ;

  SetVolume( mri );

  VolumeCropper* vc = MainWindow::GetMainWindow()->GetVolumeCropper();
  connect(ui->pushButtonReset, SIGNAL(clicked()), vc, SLOT(Reset()));
  connect(ui->pushButtonApply, SIGNAL(clicked()), vc, SLOT(Apply()));
  connect(ui->pushButtonSaveAs, SIGNAL(clicked()),
          MainWindow::GetMainWindow(), SLOT(SaveVolumeAs()));

  QSettings s;
  QByteArray val = s.value("VolumeCropper/Geometry").toByteArray();
  if (!val.isEmpty())
    restoreGeometry(val);
  else
  {
    QWidget* p = parentWidget();
    if (p)
      move(p->rect().center());
  }
}

DialogCropVolume::~DialogCropVolume()
{
  QSettings s;
  s.setValue("VolumeCropper/Geometry", saveGeometry());
  delete ui;
}

void DialogCropVolume::SetVolume( LayerMRI* mri )
{
  m_mri = mri;
  if ( m_mri )
  {
    int* dim = m_mri->GetImageData()->GetDimensions();
    for ( int i = 0; i < 6; i++ )
    {
      m_spinRange[i]->setRange( 0, dim[i/2]-1 );
    }
  }
}

void DialogCropVolume::OnSpinRange(int nVal)
{
  for ( int i = 0; i < 6; i++ )
  {
    if ( qobject_cast<QSpinBox*>(sender()) == m_spinRange[i] )
    {
      if ( (i%2 == 0 && nVal < m_spinRange[i+1]->value()) ||
           (i%2 == 1 && nVal > m_spinRange[i-1]->value()) )
        MainWindow::GetMainWindow()->GetVolumeCropper()->SetExtent( i, nVal );
      break;
    }
  }
}

void DialogCropVolume::OnCropBoundChanged(LayerMRI *mri)
{
  if (mri == m_mri)
  {
    int* ext = MainWindow::GetMainWindow()->GetVolumeCropper()->GetExtent();
    for ( int i = 0; i < 6; i++ )
    {
      m_spinRange[i]->setValue(ext[i]);
    }
  }
}

void DialogCropVolume::OnLayerRemoved(Layer *layer)
{
  if (m_mri == layer )
  {
    MainWindow::GetMainWindow()->GetVolumeCropper()->SetEnabled( false );
    hide();
  }
}

void DialogCropVolume::showEvent(QShowEvent *)
{
  RenderView3D* view = (RenderView3D*)MainWindow::GetMainWindow()->GetRenderView( 3 );
  if ( view )
  {
    m_bShowSliceFrame = view->GetShowSliceFrames();
    view->SetShowSliceFrames( false );
  }
}

void DialogCropVolume::hideEvent(QHideEvent *)
{
  RenderView3D* view = (RenderView3D*)MainWindow::GetMainWindow()->GetRenderView( 3 );
  if ( view && m_bShowSliceFrame )
  {
    view->SetShowSliceFrames( true );
  }
}
