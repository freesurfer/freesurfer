#include <QtGui/QApplication>
#include "MainWindow.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    a.setApplicationName("dummy_qt");
    MainWindow w;
    w.show();

    return a.exec();
}
