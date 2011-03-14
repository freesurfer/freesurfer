#ifndef DIALOGSAVEVOLUME_H
#define DIALOGSAVEVOLUME_H

#include <QDialog>

namespace Ui {
    class DialogSaveVolume;
}

class DialogSaveVolume : public QDialog
{
    Q_OBJECT

public:
    explicit DialogSaveVolume(QWidget *parent = 0, const QString& filepath = "");
    ~DialogSaveVolume();

    QString GetFileName();

    bool GetResample();

protected slots:
    void OnOK();
    void OnOpen();

private:
    Ui::DialogSaveVolume *ui;
};

#endif // DIALOGSAVEVOLUME_H
