#include "DialogTransformVolume.h"
#include "ui_DialogTransformVolume.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>
extern "C"
{
#include "mri.h"
}

DialogTransformVolume::DialogTransformVolume(QWidget *parent) :
    QDialog(parent),
    UIUpdateHelper(),
    ui(new Ui::DialogTransformVolume)
{
    ui->setupUi(this);
    m_checkRotate[0] = ui->checkBoxRotateX;
    m_checkRotate[1] = ui->checkBoxRotateY;
    m_checkRotate[2] = ui->checkBoxRotateZ;
    m_comboRotate[0] = ui->comboBoxRotateX;
    m_comboRotate[1] = ui->comboBoxRotateY;
    m_comboRotate[2] = ui->comboBoxRotateZ;
    m_textAngle[0] = ui->lineEditRotateX;
    m_textAngle[1] = ui->lineEditRotateY;
    m_textAngle[2] = ui->lineEditRotateZ;
    m_scrollTranslate[0] = ui->scrollBarTranslateX;
    m_scrollTranslate[1] = ui->scrollBarTranslateY;
    m_scrollTranslate[2] = ui->scrollBarTranslateZ;
    m_textTranslate[0] = ui->lineEditTranslateX;
    m_textTranslate[1] = ui->lineEditTranslateY;
    m_textTranslate[2] = ui->lineEditTranslateZ;
    m_scrollScale[0] = ui->scrollBarScaleX;
    m_scrollScale[1] = ui->scrollBarScaleY;
    m_scrollScale[2] = ui->scrollBarScaleZ;
    m_textScale[0] = ui->lineEditScaleX;
    m_textScale[1] = ui->lineEditScaleY;
    m_textScale[2] = ui->lineEditScaleZ;

    connect(MainWindow::GetMainWindow()->GetLayerCollection("MRI"), SIGNAL(ActiveLayerChanged(Layer*)),
            this, SLOT(OnActiveLayerChanged()));
    connect(ui->pushButtonSaveVolumeAs, SIGNAL(clicked()),
            MainWindow::GetMainWindow(), SLOT(SaveVolumeAs()));
}

DialogTransformVolume::~DialogTransformVolume()
{
    delete ui;
}

// scope: 0 => translate related, 1 => scale related, 2 => both
void DialogTransformVolume::UpdateUI( int scope )
{
  LayerMRI* layer = (LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( layer )
  {
      QList<QWidget*> allwidgets = this->findChildren<QWidget*>();
      for ( int i = 0; i < allwidgets.size(); i++ )
        allwidgets[i]->blockSignals( true );
    if ( scope == 0 || scope == 2 )
    {
      double* vs = layer->GetWorldVoxelSize();
      double* ws = layer->GetWorldSize();
      double pos[3];
      layer->GetTranslate( pos );
      for ( int i = 0; i < 3; i++ )
      {
        int range = (int)( ws[i] / vs[i] + 0.5 ) * 2;
        int npos = (int)(pos[i] / vs[i]) + range/2;
        m_scrollTranslate[i]->setRange(0, range);
        m_scrollTranslate[i]->setValue(npos);
        ChangeLineEditNumber(m_textTranslate[i], pos[i]);
      }
    }
    if ( scope == 1 || scope == 2 )
    {
      double scale[3];
      layer->GetScale( scale );
      for ( int i = 0; i < 3; i++ )
      {
        if ( scale[i] >= 1 )
          m_scrollScale[i]->setValue( 50 + (int)( (scale[i]-1.0)*50 ) );
        else
          m_scrollScale[i]->setValue( 50 - (int)( (1.0-scale[i])*100 ) );

        ChangeLineEditNumber(m_textScale[i], scale[i]);
      }
    }

    ui->pushButtonRestore->setEnabled( layer->IsTransformed() );
    ui->pushButtonSaveReg->setEnabled( layer->IsTransformed() );
    ui->pushButtonSaveVolumeAs->setEnabled( layer->IsTransformed() );
    for ( int i = 0; i < allwidgets.size(); i++ )
      allwidgets[i]->blockSignals( false );
  }
}

bool DialogTransformVolume::GetRotation( int nIndex_in,
                      int& plane_out,
                      double& angle_out )
{
  if ( nIndex_in < 0 ||
       nIndex_in > 2 ||
       !m_checkRotate[ nIndex_in ]->isChecked() )
    return false;

  plane_out = m_comboRotate[ nIndex_in ]->currentIndex();
  bool bOK;
  double dVal = m_textAngle[ nIndex_in ]->text().toDouble(&bOK);
  if ( bOK)
      angle_out = dVal;
  return bOK;
}

void DialogTransformVolume::OnApply()
{
  int plane;
  double angle;
  if ( !m_checkRotate[0]->isChecked() &&
       !m_checkRotate[1]->isChecked() &&
       !m_checkRotate[2]->isChecked() )
  {
      QMessageBox::warning( this, "Error",
             "Must at least select one rotation.");
    return;
  }
  else if ( ( m_checkRotate[0]->isChecked() && !GetRotation( 0, plane, angle ) ) ||
            ( m_checkRotate[1]->isChecked() && !GetRotation( 1, plane, angle ) ) ||
            ( m_checkRotate[2]->isChecked() && !GetRotation( 2, plane, angle ) ) )
  {
      QMessageBox::warning( this, "Error",
             "Please enter correct rotation angle.");
    return;
  }

  DoRotate();

  UpdateUI();
}

void DialogTransformVolume::OnSaveReg()
{
    LayerMRI* layer_mri = ( LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
    if ( !layer_mri)
    {
      return;
    }

    QString filename = QFileDialog::getSaveFileName(this, "Save Registration",
                      QFileInfo( layer_mri->GetFileName() ).absolutePath(),
                      "LTA files (*.lta);;All files (*.*)");
    if ( !filename.isEmpty() )
    {
      layer_mri->SaveRegistration( filename );
    }
}

void DialogTransformVolume::OnRestore()
{
  LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( layer )
  {
    layer->Restore();
    UpdateUI();
  }
}

void DialogTransformVolume::DoRotate()
{
  LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( layer )
  {
    RotationElement re;
    if ( ui->radioButtonAroundCursor->isChecked() )
    {
      MainWindow::GetMainWindow()->GetLayerCollection( "MRI" )->
        GetSlicePosition( re.Point );
      layer->RemapPositionToRealRAS( re.Point, re.Point );
    }
    else
    {
    // use center of the volume to rotate
      layer->GetRASCenter( re.Point );
    }
    re.SampleMethod = SAMPLE_TRILINEAR;
    if ( ui->radioButtonNearestNeighbor->isChecked() )
      re.SampleMethod = SAMPLE_NEAREST;
//    else if ( m_radioSinc->GetValue() )
//      re.SampleMethod = SAMPLE_SINC;

    std::vector<RotationElement> rotations;
    for ( int i = 0; i < 3; i++ )
    {
      if ( GetRotation( i, re.Plane, re.Angle ) )
      {
        rotations.push_back( re );
      }
    }
    MainWindow::GetMainWindow()->RotateVolume( rotations, false );
  }
}

void DialogTransformVolume::OnActiveLayerChanged()
{
    if ( isVisible() )
        UpdateUI();
}


void DialogTransformVolume::OnLineEditTranslateX(const QString& text)
{
  RespondTextTranslate( 0 );
}

void DialogTransformVolume::OnLineEditTranslateY(const QString& text)
{
  RespondTextTranslate( 1 );
}

void DialogTransformVolume::OnLineEditTranslateZ(const QString& text)
{
  RespondTextTranslate( 2 );
}

void DialogTransformVolume::OnScrollBarTranslateX(int nVal)
{
  RespondScrollTranslate( 0 );
}

void DialogTransformVolume::OnScrollBarTranslateY(int nVal)
{
  RespondScrollTranslate( 1 );
}

void DialogTransformVolume::OnScrollBarTranslateZ(int nVal)
{
  RespondScrollTranslate( 2 );
}

void DialogTransformVolume::RespondTextTranslate( int n )
{
  if ( isVisible() )
  {
      bool bOK;
    double dvalue =m_textTranslate[n]->text().toDouble(&bOK);
    if ( bOK )
    {
      LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
      if ( layer )
      {
        double pos[3];
        layer->GetTranslate( pos );
        pos[n] = dvalue;
        layer->Translate( pos );
        MainWindow::GetMainWindow()->RequestRedraw();

        double* vs = layer->GetWorldVoxelSize();
        int range = m_scrollTranslate[n]->maximum();
        m_scrollTranslate[n]->blockSignals(true);
        m_scrollTranslate[n]->setValue(range/2 + (int)( pos[n] / vs[n] ) );
        m_scrollTranslate[n]->blockSignals(false);
        UpdateUI( 1 );
      }
    }
  }
}

void DialogTransformVolume::RespondScrollTranslate( int n )
{
  if ( isVisible() )
  {
    LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
    if ( layer )
    {
      double pos[3];
      layer->GetTranslate( pos );
      int range = m_scrollTranslate[n]->maximum();
      int npos = m_scrollTranslate[n]->value();
      double* vs = layer->GetWorldVoxelSize();
      pos[n] = ( npos - range/2 ) * vs[n];
      layer->Translate( pos );
      MainWindow::GetMainWindow()->RequestRedraw();
      ChangeLineEditNumber(m_textTranslate[n], pos[n] );
      UpdateUI( 1 );
    }
  }
}


void DialogTransformVolume::OnLineEditScaleX(const QString& text)
{
  RespondTextScale( 0 );
}

void DialogTransformVolume::OnLineEditScaleY(const QString& text)
{
  RespondTextScale( 1 );
}

void DialogTransformVolume::OnLineEditScaleZ(const QString& text)
{
  RespondTextScale( 2 );
}

void DialogTransformVolume::OnScrollBarScaleX(int nVal)
{
  RespondScrollScale( 0 );
}

void DialogTransformVolume::OnScrollBarScaleY(int nVal)
{
  RespondScrollScale( 1 );
}

void DialogTransformVolume::OnScrollBarScaleZ(int nVal)
{
  RespondScrollScale( 2 );
}

void DialogTransformVolume::RespondTextScale( int n )
{
  if ( isVisible() )
  {
      bool bOK;
    double dvalue = m_textScale[n]->text().toDouble(&bOK);
    if ( bOK && dvalue > 0 )
    {
      LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
      if ( layer )
      {
        double scale[3];
        layer->GetScale( scale );
        scale[n] = dvalue;
        layer->Scale( scale );
        MainWindow::GetMainWindow()->RequestRedraw();

        m_scrollScale[n]->blockSignals(true);
        if ( dvalue >= 1 )
          m_scrollScale[n]->setValue( 50 + (int)( (dvalue-1.0)*50 ) );
        else
          m_scrollScale[n]->setValue( 50 - (int)( (1.0-dvalue)*100 ) );
        m_scrollScale[n]->blockSignals(false);
        UpdateUI( 0 );
      }
    }
  }
}

void DialogTransformVolume::RespondScrollScale( int n )
{
  LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( layer )
  {
    double scale[3];
    layer->GetScale( scale );
    int npos = m_scrollScale[n]->value();
    if ( npos >= 50 )
      scale[n] = ( npos - 50 ) / 50.0 + 1.0;
    else
      scale[n] = ( npos - 50 ) / 100.0 + 1.0;
    layer->Scale( scale );
    MainWindow::GetMainWindow()->RequestRedraw();

    ChangeLineEditNumber( m_textScale[n], scale[n] );
    UpdateUI( 0 );
  }
}
