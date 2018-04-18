#include "DialogSelectSplines.h"
#include "ui_DialogSelectSplines.h"
#include "LayerPointSet.h"
#include <QListWidgetItem>
#include <QMessageBox>

DialogSelectSplines::DialogSelectSplines(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSelectSplines)
{
  ui->setupUi(this);
//  connect(ui->listWidget, SIGNAL(itemChanged(QListWidgetItem*)), this, SLOT(OnItemChanged(QListWidgetItem*)));
  connect(ui->pushButtonOk, SIGNAL(clicked(bool)), SLOT(OnButtonOK()));
  connect(ui->toolButtonDown, SIGNAL(clicked(bool)), SLOT(OnButtonDown()));
  connect(ui->toolButtonUp, SIGNAL(clicked(bool)), SLOT(OnButtonUp()));
}

DialogSelectSplines::~DialogSelectSplines()
{
  delete ui;
}

void DialogSelectSplines::SetPointSets(const QList<Layer *> &list)
{
  ui->listWidget->clear();
  foreach (Layer* ps, list)
  {
    QListWidgetItem* item = new QListWidgetItem;
    item->setText(ps->GetName());
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    item->setCheckState(ps->IsVisible()?Qt::Checked:Qt::Unchecked);
    item->setData(Qt::UserRole, QVariant::fromValue<QObject*>(ps));
    ui->listWidget->insertItem(0, item);
  }
}

QList<LayerPointSet*> DialogSelectSplines::GetSelectedPointSets()
{
  QList<LayerPointSet*> list;
  for (int i = 0; i < ui->listWidget->count(); i++)
  {
    QListWidgetItem* item = ui->listWidget->item(i);
    if (item->checkState() == Qt::Checked)
    {
      LayerPointSet* ps = qobject_cast<LayerPointSet*>(item->data(Qt::UserRole).value<QObject*>());
      if (ps)
        list << ps;
    }
  }
  return list;
}

void DialogSelectSplines::OnButtonUp()
{
  int n = ui->listWidget->currentRow();
  if (n >= 0)
  {
    QListWidgetItem* item = ui->listWidget->takeItem(n);
    if (item)
    {
      n--;
      if (n < 0)
        n = ui->listWidget->count();
      ui->listWidget->insertItem(n, item);
      ui->listWidget->setCurrentRow(n);
    }
  }
}

void DialogSelectSplines::OnButtonDown()
{
  int n = ui->listWidget->currentRow();
  if (n >= 0)
  {
    QListWidgetItem* item = ui->listWidget->takeItem(n);
    if (item)
    {
      n++;
      if (n > ui->listWidget->count())
        n = 0;
      ui->listWidget->insertItem(n, item);
      ui->listWidget->setCurrentRow(n);
    }
  }
}

void DialogSelectSplines::OnButtonOK()
{
  if (GetSelectedPointSets().isEmpty())
    QMessageBox::warning(this, "Error", "No layer was selected.");
  else
    accept();
}
