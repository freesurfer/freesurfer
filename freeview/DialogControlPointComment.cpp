#include "DialogControlPointComment.h"
#include "ui_DialogControlPointComment.h"
#include <QMessageBox>

DialogControlPointComment::DialogControlPointComment(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogControlPointComment)
{
  ui->setupUi(this);

  m_listCheckBoxes << ui->checkBoxMain << ui->checkBoxMain_2 << ui->checkBoxMain_3
    << ui->checkBoxMain_4 << ui->checkBoxMain_5 << ui->checkBoxMain_6
    << ui->checkBoxMain_7
    << ui->checkBoxSolution << ui->checkBoxSolution_2 << ui->checkBoxSolution_3
    << ui->checkBoxSolution_4
    << ui->checkBoxProgress << ui->checkBoxProgress_2 << ui->checkBoxProgress_3;

  foreach (QCheckBox* box, m_listCheckBoxes)
    connect(box, SIGNAL(toggled(bool)), SLOT(OnCheckBoxToggled(bool)), Qt::QueuedConnection);
}

DialogControlPointComment::~DialogControlPointComment()
{
  delete ui;
}

void DialogControlPointComment::SetLabel(const QString &text)
{
  ui->label->setText(text);
}

QStringList DialogControlPointComment::GetPrefilledItems()
{
  QStringList list;
  foreach (QCheckBox* box, m_listCheckBoxes)
  {
    if (box->isChecked())
    {
      list << box->text();
    }
  }
  return list;
}

void DialogControlPointComment::OnCheckBoxToggled(bool bChecked)
{
}

QString DialogControlPointComment::GetComment()
{
  return ui->textEdit->toPlainText().trimmed();
}

void DialogControlPointComment::SetComment(const QString &text, const QStringList& prefilled_items)
{
  ui->textEdit->clear();
  ui->textEdit->appendPlainText(text);

  foreach (QCheckBox* box, m_listCheckBoxes)
  {
    if (prefilled_items.contains(box->text(), Qt::CaseInsensitive))
    {
      box->blockSignals(true);
      box->setChecked(true);
      box->blockSignals(false);
    }
  }
}

void DialogControlPointComment::OnOk()
{
  if (GetComment().isEmpty())
    QMessageBox::warning(this, "Error", "Please enter your comment");
  else
    accept();
}
