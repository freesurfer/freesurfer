#include "DialogLoadSurface.h"
#include "ui_DialogLoadSurface.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>
#include "MainWindow.h"
#include <QSettings>
#include <QLineEdit>
#include "MigrationDefs.h"

DialogLoadSurface::DialogLoadSurface(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadSurface)
{
  ui->setupUi(this);

  QSettings s;

  UpdateStatus();
}

DialogLoadSurface::~DialogLoadSurface()
{
  QSettings s;
  s.setValue("DialogLoadSurface/Filename", ui->comboBoxFilename->currentText().trimmed());

  delete ui;
}

void DialogLoadSurface::accept()
{
  if (ui->comboBoxFilename->currentText().trimmed().isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please enter the filename of the surface file to load");
  }
  else
    QDialog::accept();
}

void DialogLoadSurface::OnOpen()
{
  QString lastdir = QFileInfo(ui->comboBoxFilename->currentText().trimmed()).absoluteFilePath();
  QString filename = QFileDialog::getOpenFileName( this, "Select surface file",
                                                   MainWindow::AutoSelectLastDir(lastdir, "surf" ),
                                                   "Surface files (*)");
  if (!filename.isEmpty())
  {
    ui->comboBoxFilename->setCurrentText(filename);
    ui->comboBoxFilename->lineEdit()->setCursorPosition(filename.size());
  }

  UpdateStatus();
}

void DialogLoadSurface::SetRecentFiles( const QStringList& filenames )
{
  QStringList fns = filenames;
  fns.insert(0, "current folder");
  ui->comboBoxFilename->clear();
  ui->comboBoxFilename->addItems( fns );
  if ( !filenames.isEmpty() )
  {
    ui->comboBoxFilename->setCurrentIndex( 0 );
    ui->comboBoxFilename->setCurrentText("");
    ui->comboBoxFilename->lineEdit()->setCursorPosition( ui->comboBoxFilename->currentText().size() );
  }
}

void DialogLoadSurface::UpdateStatus()
{
}

QStringList DialogLoadSurface::GetFilenames()
{
  QStringList fns = ui->comboBoxFilename->currentText().trimmed().split( QRegularExpression( "[;]" ), MD_SkipEmptyParts );
  return fns;
}

bool DialogLoadSurface::GetLoadAll()
{
  return ui->checkBoxAll->isChecked();
}
