#include "MainWindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  MainWindow w;
  if (argc > 2 && QString(argv[1]) == "-py-cmd")
    w.m_strPythonCmd = argv[2];
  w.ShowDialog();
  return a.exec();
}
