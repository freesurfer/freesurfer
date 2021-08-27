#include "DialogTransformSurface.h"
#include "ui_DialogTransformSurface.h"
#include "vtkSmartPointer.h"
#include "vtkTransform.h"
#include "MainWindow.h"
#include "LayerSurface.h"
#include <QFileDialog>
#include "LayerCollection.h"
#include <QMessageBox>
#include <QDebug>

DialogTransformSurface::DialogTransformSurface(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogTransformSurface)
{
  ui->setupUi(this);
  ui->pushButtonRedo->hide();
  ui->pushButtonApply->hide();
  ui->tabWidget->removeTab(3);
  ui->tabFlip->hide();

  m_scrollRotate[0] = ui->scrollBarRotateX;
  m_scrollRotate[1] = ui->scrollBarRotateY;
  m_scrollRotate[2] = ui->scrollBarRotateZ;
  m_textRotate[0] = ui->lineEditRotateX;
  m_textRotate[1] = ui->lineEditRotateY;
  m_textRotate[2] = ui->lineEditRotateZ;
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

  for (int i = 0; i < 3; i++)
  {
    connect(m_scrollRotate[i], SIGNAL(valueChanged(int)), SLOT(OnScrollBarValueChanged(int)));
    connect(m_scrollTranslate[i], SIGNAL(valueChanged(int)), SLOT(OnScrollBarValueChanged(int)));
    connect(m_scrollScale[i], SIGNAL(valueChanged(int)), SLOT(OnScrollBarValueChanged(int)));
    connect(m_textRotate[i], SIGNAL(textChanged(QString)), SLOT(OnTextEdit(QString)));
    connect(m_textTranslate[i], SIGNAL(textChanged(QString)), SLOT(OnTextEdit(QString)));
    connect(m_textScale[i], SIGNAL(textChanged(QString)), SLOT(OnTextEdit(QString)));
  }

  connect(ui->radioButtonAroundCenter, SIGNAL(toggled(bool)), SLOT(ResetUI()));
  connect(ui->radioButtonAroundCursor, SIGNAL(toggled(bool)), SLOT(ResetUI()));

  connect(ui->pushButtonRestore, SIGNAL(clicked(bool)), SLOT(OnButtonRestore()));
  connect(ui->pushButtonUndo, SIGNAL(clicked(bool)), SLOT(OnButtonUndo()));
  connect(ui->pushButtonSaveReg, SIGNAL(clicked(bool)), SLOT(OnButtonSaveTransform()));

  UpdateUI();
}

DialogTransformSurface::~DialogTransformSurface()
{
  delete ui;
}


void DialogTransformSurface::showEvent(QShowEvent *)
{
  ResetUI();
}

void DialogTransformSurface::UpdateUI()
{
  LayerSurface* layer = ( LayerSurface* )MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  ui->pushButtonUndo->setEnabled(layer && layer->HasTransformUndo());
}

void DialogTransformSurface::ResetUI(const QList<QWidget *> &excluded)
{
  QList<QWidget*> allwidgets = this->findChildren<QWidget*>();
  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( true );
  }

  for (int i = 0; i < 3; i++)
  {
    if (!excluded.contains(m_scrollRotate[i]))
      m_scrollRotate[i]->setValue(0);
    if (!excluded.contains(m_scrollTranslate[i]))
      m_scrollTranslate[i]->setValue(500);
    if (!excluded.contains(m_scrollScale[i]))
      m_scrollScale[i]->setValue(200);
    if (!excluded.contains(m_textRotate[i]))
      m_textRotate[i]->setText("0");
    if (!excluded.contains(m_textTranslate[i]))
      m_textTranslate[i]->setText("0");
    if (!excluded.contains(m_textScale[i]))
      m_textScale[i]->setText("1");
  }

  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( false );
  }
}

void DialogTransformSurface::OnScrollBarValueChanged(int n_in)
{
  static QObject* sw = NULL;
  static double cb_pt[3] = {0,0,0};

  LayerSurface* layer = ( LayerSurface* )MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  QList<QWidget*> excluded;
  excluded << qobject_cast<QWidget*>(sender());
  vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
  if (layer && sw != sender())
    layer->GetCenterOfActor(cb_pt);
  double dIncrementTranslate = 0.1;
  for (int i = 0; i < 3; i++)
  {
    if (m_scrollRotate[i] == sender())
    {
      double cpt[3];
      if (ui->radioButtonAroundCenter->isChecked())
      {
        cpt[0] = cb_pt[0];
        cpt[1] = cb_pt[1];
        cpt[2] = cb_pt[2];
      }
      else
      {
        MainWindow::GetMainWindow()->GetLayerCollection( "Surface" )->
            GetSlicePosition(cpt);
      }
      t->Translate(cpt[0], cpt[1], cpt[2]);
      double dval = n_in/2.0;
      if (i == 0)
        t->RotateX(dval);
      else if (i == 1)
        t->RotateY(dval);
      else
        t->RotateZ(dval);
      t->Translate(-cpt[0], -cpt[1], -cpt[2]);

      ChangeLineEditNumber(m_textRotate[i], dval, 2, true );
      excluded << m_textRotate[i];
      break;
    }
    else if (m_scrollTranslate[i] == sender())
    {
      double pos[3] = {0,0,0};
      int range = m_scrollTranslate[i]->maximum();
      int npos = n_in;
      pos[i] = ( npos - range/2 ) * dIncrementTranslate;
      t->Translate( pos );
      ChangeLineEditNumber(m_textTranslate[i], pos[i], 2, true);
      excluded << m_textTranslate[i];
      break;
    }
    else if (m_scrollScale[i] == sender())
    {
      double scale[3] = {1,1,1};
      int npos = n_in;
      double nmax = m_scrollScale[i]->maximum();
      if ( npos >= nmax/2 )
      {
        scale[i] = ( npos - nmax/2 ) / (nmax/2) + 1.0;
      }
      else
      {
        scale[i] = ( npos - nmax/2 ) / nmax + 1.0;
      }
      t->Translate(cb_pt[0], cb_pt[1], cb_pt[2]);
      t->Scale( scale );
      t->Translate(-cb_pt[0], -cb_pt[1], -cb_pt[2]);
      ChangeLineEditNumber( m_textScale[i], scale[i], 4, true );
      excluded << m_textScale[i];
      break;
    }
  }

  if (sw != sender())
  {
    ResetUI(excluded);
    sw = sender();
    if (layer)
      layer->AppendIdentityTransform();
  }
  if (layer)
    layer->UpdateLastTransform(t);

  UpdateUI();
}

void DialogTransformSurface::OnTextEdit(const QString &text)
{
  static QObject* sw = NULL;
  static double cb_pt[3] = {0,0,0};

  LayerSurface* layer = ( LayerSurface* )MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  QList<QWidget*> excluded;
  excluded << qobject_cast<QWidget*>(sender());

  bool bOK;
  double dvalue = text.toDouble(&bOK);
  if ( !bOK )
    return;

  if (layer && sw != sender())
    layer->GetCenterOfActor(cb_pt);
  vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
  double dIncrementTranslate = 0.1;
  for (int i = 0; i < 3; i++)
  {
    if (sender() == m_textRotate[i])
    {
      double cpt[3];
      if (ui->radioButtonAroundCenter->isChecked())
      {
        cpt[0] = cb_pt[0];
        cpt[1] = cb_pt[1];
        cpt[2] = cb_pt[2];
      }
      else
      {
        MainWindow::GetMainWindow()->GetLayerCollection( "Surface" )->
            GetSlicePosition(cpt);
      }
      t->Translate(cpt[0], cpt[1], cpt[2]);
      if (i == 0)
        t->RotateX(dvalue);
      else if (i == 1)
        t->RotateY(dvalue);
      else
        t->RotateZ(dvalue);
      t->Translate(-cpt[0], -cpt[1], -cpt[2]);

      while (dvalue > 180)
        dvalue -= 360;
      while (dvalue < -180)
        dvalue += 360;

      m_scrollRotate[i]->blockSignals(true);
      m_scrollRotate[i]->setValue( (int)(dvalue*2) );
      m_scrollRotate[i]->blockSignals(false);
      excluded << m_scrollRotate[i];
      break;
    }
    else if (sender() == m_textTranslate[i])
    {
      double pos[3] = {0, 0, 0};
      pos[i] = dvalue;
      t->Translate( pos );
      int range = m_scrollTranslate[i]->maximum();
      m_scrollTranslate[i]->blockSignals(true);
      m_scrollTranslate[i]->setValue(range/2 + (int)( pos[i] / dIncrementTranslate ) );
      m_scrollTranslate[i]->blockSignals(false);
      excluded << m_scrollTranslate[i];
      break;
    }
    else if (sender() == m_textScale[i] && dvalue > 0)
    {
      double scale[3] = {1,1,1};
      scale[i] = dvalue;
      t->Translate(cb_pt[0], cb_pt[1], cb_pt[2]);
      t->Scale( scale );
      t->Translate(-cb_pt[0], -cb_pt[1], -cb_pt[2]);
      m_scrollScale[i]->blockSignals(true);
      int nmax = m_scrollScale[i]->maximum();
      if ( dvalue >= 1 )
      {
        m_scrollScale[i]->setValue( nmax/2 + (int)( (dvalue-1.0)*nmax/2 ) );
      }
      else
      {
        m_scrollScale[i]->setValue( nmax/2 - (int)( (1.0-dvalue)*nmax ) );
      }
      m_scrollScale[i]->blockSignals(false);
      excluded << m_scrollScale[i];
      break;
    }
  }
  if (sw != sender())
  {
    ResetUI(excluded);
    sw = sender();
    if (layer)
      layer->AppendIdentityTransform();
  }
  if (layer)
    layer->UpdateLastTransform(t);

  UpdateUI();
}

void DialogTransformSurface::OnButtonRestore()
{
  LayerSurface* layer = ( LayerSurface* )MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if (layer)
    layer->ResetTransform();
  ResetUI();
  UpdateUI();
}

void DialogTransformSurface::OnButtonUndo()
{
  LayerSurface* layer = ( LayerSurface* )MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if (layer)
    layer->UndoLastTransform();
  ResetUI();
  UpdateUI();
}

void DialogTransformSurface::OnButtonSaveTransform()
{
  LayerSurface* layer = ( LayerSurface* )MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if (!layer)
    return;

  QString filename = QFileDialog::getSaveFileName(this, "Save Transform",
                                                  QFileInfo( layer->GetFileName() ).absoluteFilePath(),
                                                  "LTA files (*.lta);;All files (*)");
  if ( !filename.isEmpty() && !layer->SaveTransform(filename))
  {
    QMessageBox::warning(this, "Error", "Failed to save transform");
    qDebug() << "Failed to save transform";
  }
}
