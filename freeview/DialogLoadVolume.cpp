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
#include "DialogLoadVolume.h"
#include "ui_DialogLoadVolume.h"
#include "LayerPropertyMRI.h"
#include "MainWindow.h"
#include "LUTDataHolder.h"
#include "MyUtils.h"
#include <QFileDialog>
#include <QFileInfo>
#include <QFileDialog>
#include <QMessageBox>



#include "mri.h"


DialogLoadVolume::DialogLoadVolume(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadVolume)
{
  ui->setupUi(this);
#ifdef Q_OS_MAC
  setWindowFlags( Qt::Sheet );
#endif

  UpdateLUT();
}

DialogLoadVolume::~DialogLoadVolume()
{
  delete ui;
}

void DialogLoadVolume::UpdateLUT()
{
  ui->comboBoxLUT->blockSignals( true );
  LUTDataHolder* luts = MainWindow::GetMainWindow()->GetLUTData();
  ui->comboBoxLUT->clear();
  for ( int i = 0; i < luts->GetCount(); i++ )
  {
    ui->comboBoxLUT->addItem( luts->GetName( i ) );
  }
  ui->comboBoxLUT->addItem( "Load lookup table..." );
  ui->comboBoxLUT->setCurrentIndex( 0 );
  ui->comboBoxLUT->blockSignals( false );
}

void DialogLoadVolume::OnOpen()
{
  QString fn = ui->comboBoxFilenames->currentText().trimmed();
  if (fn == "current folder")
    fn = QDir::currentPath();
  else
    fn =  QFileInfo(fn).absolutePath();
  QStringList filenames = QFileDialog::getOpenFileNames( this, "Select volume files",
                                                         fn, // MainWindow::AutoSelectLastDir( m_strLastDir, "mri" ),
                                                         "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !filenames.isEmpty() )
  {
    m_strLastDir = QFileInfo( filenames[0] ).canonicalPath();
    for (int i = 0; i < filenames.size(); i++)
    {
      filenames[i] = MyUtils::Win32PathProof(filenames[i]);
    }
    QString strg = filenames.join( ";" );
    ui->comboBoxFilenames->setEditText( strg );
    ui->comboBoxFilenames->lineEdit()->setCursorPosition( strg.size() );
  }
}

void DialogLoadVolume::OnOpenRegistration()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select registration file",
                                                   m_strLastDir,
                                                   "Registration files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditRegistration->setText( MyUtils::Win32PathProof(filename) );
    ui->lineEditRegistration->setCursorPosition( ui->lineEditRegistration->text().size() );
  }
}

void DialogLoadVolume::SetRecentFiles( const QStringList& filenames )
{
  QStringList fns = filenames;
  for (int i = 0; i < fns.size(); i++)
  {
    fns[i] = MyUtils::Win32PathProof(fns[i]);
  }
  fns.insert(0, "current folder");
  ui->comboBoxFilenames->clear();
  ui->comboBoxFilenames->addItems( fns );
  if ( !filenames.isEmpty() )
  {
    ui->comboBoxFilenames->setCurrentIndex( 0 );
    ui->comboBoxFilenames->lineEdit()->setCursorPosition( ui->comboBoxFilenames->currentText().size() );
  }
}

void DialogLoadVolume::OnColorMap( int nSel )
{
  ui->comboBoxLUT->setEnabled( nSel == LayerPropertyMRI::LUT );
  ui->labelLUT->setEnabled( nSel == LayerPropertyMRI::LUT );
}

void DialogLoadVolume::OnLUT( int nSel )
{
  LUTDataHolder* luts = MainWindow::GetMainWindow()->GetLUTData();
  if ( nSel >= luts->GetCount() )
  {
    QString filename = QFileDialog::getOpenFileName( this, "Load lookup table file",
                                                     m_strLastDir,
                                                     "LUT files (*)" );
    if ( !filename.isEmpty() && luts->LoadColorTable( filename ) )
    {
      UpdateLUT();
      ui->comboBoxLUT->setCurrentIndex( luts->GetCount() - 1 );
    }
    else
      ui->comboBoxLUT->setCurrentIndex(0);
  }
}

QStringList DialogLoadVolume::GetVolumeFileNames()
{
  QStringList fns = ui->comboBoxFilenames->currentText().trimmed().split( QRegExp( "[;]" ), QString::SkipEmptyParts );
  for (int i = 0; i < fns.size(); i++)
  {
    fns[i] = MyUtils::CygwinPathProof(fns[i]);
  }
  return fns;
}

QString DialogLoadVolume::GetRegFileName()
{
  if ( ui->checkBoxRegistration->isChecked() )
  {
    return MyUtils::CygwinPathProof(ui->lineEditRegistration->text().trimmed());
  }
  else
  {
    return "";
  }
}

bool DialogLoadVolume::IsToResample()
{
  return ui->checkBoxResampleToRAS->isChecked();
}

int DialogLoadVolume::GetSampleMethod()
{
  if ( ui->radioNearest->isChecked() )
  {
    return SAMPLE_NEAREST;
  }
  else if (ui->radioTrilinear->isChecked())
  {
    return SAMPLE_TRILINEAR;
  }
  else
    return SAMPLE_CUBIC_BSPLINE;
}

QString DialogLoadVolume::GetColorMap()
{
  QStringList names;
  names << "grayscale" << "lut" << "heat" << "jet" << "gecolor" << "nih" << "pet";
  return names[ui->comboBoxColorMap->currentIndex()];
}

QString DialogLoadVolume::GetLUT()
{
  return ui->comboBoxLUT->currentText();
}

void DialogLoadVolume::OnOK()
{
  if ( GetVolumeFileNames().isEmpty() || ui->comboBoxFilenames->currentText().trimmed() == "current folder")
  {
    QMessageBox::warning( this, "Error", "Please specify volume file to load.");
    return;
  }
  if ( ui->checkBoxRegistration->isChecked() && GetRegFileName().isEmpty() )
  {
    QMessageBox::warning( this, "Error", "Please specify registration file to use.");
    return;
  }
  accept();
}

bool DialogLoadVolume::GetLoadAsVector()
{
  return ui->checkBoxAsVector->isChecked();
}
