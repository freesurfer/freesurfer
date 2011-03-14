#ifndef DIALOGPREFERENCES_H
#define DIALOGPREFERENCES_H

#include <QDialog>
#include <QVariantMap>
#include "UIUpdateHelper.h"

namespace Ui {
    class DialogPreferences;
}

class QAbstractButton;

class DialogPreferences : public QDialog, public UIUpdateHelper
{
    Q_OBJECT

public:
    explicit DialogPreferences(QWidget *parent = 0);
    ~DialogPreferences();

    void SetSettings(const QVariantMap& map);

    QVariantMap GetSettings();

protected slots:
    void OnClicked(QAbstractButton* btn);

private:
    Ui::DialogPreferences *ui;
};

#endif // DIALOGPREFERENCES_H
