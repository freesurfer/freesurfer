#ifndef MACHELPER_H
#define MACHELPER_H

#include <QObject>
#include <QPointer>
#include <QPixmap>
#include <QIcon>

class QWidget;

class MacHelper : public QObject
{
    Q_OBJECT
public:
    MacHelper(QObject* parent = 0);

    static bool IsDarkMode();

    static QPixmap InvertPixmap(const QPixmap& pix);

    static QIcon InvertIcon(const QIcon& icn_in, const QSize& size = QSize(), bool bTwoStates = false);

};

#endif // MACHELPER_H
