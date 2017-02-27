#include <QtGlobal>
#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
#include <QtWidgets>
#else
#include <QtGui/QApplication>
#endif
#include "MainWindow.h"
#include <QDebug>

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    qWarning() << "nmovie_qt <image file> <image file> ...\n";
    return -1;
  }

  QApplication a(argc, argv);
  MainWindow w;
  w.show();
  w.Load(argc, argv);

  return a.exec();
}
