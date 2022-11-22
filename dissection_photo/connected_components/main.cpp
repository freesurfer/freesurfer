#include "MainWindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  MainWindow w;
  for (int i = 1; i < argc; i++)
  {
    if (QString(argv[i]) == "-prof")
      w.m_bProfiling = true;
    else if (QString(argv[i]) == "-py-cmd" && i+1 < argc)
    {
      w.m_strPythonCmd = argv[i+1];
      i++;
    }
  }
  w.ShowDialog();
  return a.exec();
}
