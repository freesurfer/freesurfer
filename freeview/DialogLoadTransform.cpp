#include "DialogLoadTransform.h"
#include "ui_DialogLoadTransform.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>



#include "mri.h"



DialogLoadTransform::DialogLoadTransform(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadTransform)
{
  ui->setupUi(this);
  connect(ui->toolButtonOpenRegistration, SIGNAL(clicked(bool)), SLOT(OnButtonOpen()));
}

DialogLoadTransform::~DialogLoadTransform()
{
  delete ui;
}

void DialogLoadTransform::OnOK()
{
  if (ui->lineEditRegistration->text().trimmed().isEmpty())
  {
    QMessageBox::warning(this, tr("Error"), tr("Please select transformation file to Load"));
  }
  else
    accept();
}

QString DialogLoadTransform::GetFilename()
{
  return ui->lineEditRegistration->text().trimmed();
}

int DialogLoadTransform::GetSampleMethod()
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

void DialogLoadTransform::OnButtonOpen()
{
  QString fn = QFileDialog::getOpenFileName( this, "Select registration file",
                                             m_strLastDir,
                                             "Registration files (*)");
  if ( !fn.isEmpty() )
  {
    ui->lineEditRegistration->setText(fn);
    ui->lineEditRegistration->setCursorPosition( fn.size() );
    m_strLastDir = QFileInfo(fn).absolutePath();
  }
}
