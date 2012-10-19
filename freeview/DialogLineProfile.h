#ifndef DIALOGLINEPROFILE_H
#define DIALOGLINEPROFILE_H

#include <QDialog>

namespace Ui {
    class DialogLineProfile;
}

class LayerLineProfile;

class DialogLineProfile : public QDialog
{
    Q_OBJECT

public:
    explicit DialogLineProfile(QWidget *parent = 0);
    ~DialogLineProfile();

    double GetResolution();
    int    GetNumberOfSamples();

public slots:
    void UpdatePointSetList();

protected slots:
    void OnCompute();
    void OnExport();
    void OnComboIsoLine(int sel);

private:
    Ui::DialogLineProfile *ui;

    LayerLineProfile*   m_lineProfile;
};

#endif // DIALOGLINEPROFILE_H
