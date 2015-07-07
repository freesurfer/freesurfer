#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QMessageBox>
#include <QSettings>
#include <QKeyEvent>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->mainToolBar->hide();
    QSettings s;
    restoreGeometry(s.value("MainWindow/Geometry").toByteArray());
    connect(ui->widgetRender, SIGNAL(CurrentImageChanged(QImage, int)),
            this, SLOT(OnCurrentImageChanged(QImage, int)));
    connect(ui->sliderFrame, SIGNAL(valueChanged(int)), ui->widgetRender, SLOT(SetCurrentImageIndex(int)));
}

MainWindow::~MainWindow()
{
  QSettings s;
  s.setValue("MainWindow/Geometry", saveGeometry());
  delete ui;
}

void MainWindow::Load(int argc, char *argv[])
{
  QStringList files;
  for (int i = 1; i < argc; i++)
    files << argv[i];

  if (ui->widgetRender->LoadImages(files) == 0)
  {
    QMessageBox::warning(this, "Error Loading", "No image file was loaded successfully.");
    close();
  }
  ui->sliderFrame->setRange(0, ui->widgetRender->GetNumberOfImages()-1);
}

void MainWindow::keyPressEvent(QKeyEvent *e)
{
  if (e->key() == Qt::Key_Left)
    ui->widgetRender->OnBack();
  else if (e->key() == Qt::Key_Right)
    ui->widgetRender->OnForward();
  else if (e->key() == Qt::Key_Escape)
    close();
}

void MainWindow::OnCurrentImageChanged(const QImage &image, int n)
{
  setWindowTitle(image.text("FileName") + " - nmovie");
//  ui->statusBar->showMessage(image.text("FullPath") + QString("  [%1 x %2]").arg(image.width()).arg(image.height()));
  ui->sliderFrame->blockSignals(true);
  ui->sliderFrame->setValue(n);
  ui->sliderFrame->blockSignals(false);
}
