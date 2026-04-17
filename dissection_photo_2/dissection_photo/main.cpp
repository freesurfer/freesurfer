#include "MainWindow.h"
#include "DialogWelcome.h"

#include <QApplication>
#include <QSettings>

int main(int argc, char *argv[])
{
  qputenv("QT_AUTO_SCREEN_SCALE_FACTOR", "1");
  QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
  QApplication a(argc, argv);

  QSettings s;
  if (!s.value("WelcomeWindow/DoNotShow").toBool())
  {
    DialogWelcome dlg;
    dlg.exec();
    if (dlg.GetDoNotShow())
      s.setValue("WelcomeWindow/DoNotShow", true);
  }

  MainWindow w;
  w.show();
  return a.exec();
}
