#include "MainWindow.h"
#include <QtWidgets>

MainWindow::MainWindow() : textEdit(new QPlainTextEdit)
{
    setCentralWidget(textEdit);

    createActions();
    createStatusBar();

}

void MainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::createActions()
{

}
