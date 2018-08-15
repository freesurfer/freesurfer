#include <QApplication>
#include "configurationfileform.h"


int main(int argc, char **argv)
{
    QApplication app (argc, argv);

    ConfigurationFileForm mainWin;
    mainWin.show();

    return app.exec();
}
