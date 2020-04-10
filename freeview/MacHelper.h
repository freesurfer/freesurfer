#ifndef MACHELPER_H
#define MACHELPER_H

#include <QObject>
#include <QPointer>

class QWidget;

class MacHelper : public QObject
{
    Q_OBJECT
public:
    MacHelper(QObject* parent = 0);

    static bool IsDarkMode();

};

#endif // MACHELPER_H
