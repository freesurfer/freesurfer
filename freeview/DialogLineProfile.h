#ifndef DIALOGLINEPROFILE_H
#define DIALOGLINEPROFILE_H

#include <QDialog>

namespace Ui {
    class DialogLineProfile;
}

class LayerLineProfile;
class LayerPointSet;

class DialogLineProfile : public QDialog
{
    Q_OBJECT

public:
    explicit DialogLineProfile(QWidget *parent = 0);
    ~DialogLineProfile();

    double GetResolution();
    double GetSpacing();
    double GetOffset();
    int    GetNumberOfSamples();

public slots:
    void UpdatePointSetList();

protected slots:
    void OnCompute();
    void OnExport();
    void OnSave();
    void OnLoad();
    void OnComboIsoLine(int sel);
    void OnSliderOpacity(int);
    void OnEditRadius(const QString& strg);
    void OnColorPicker(const QColor& color);
    void OnLineProfileIdPicked(LayerLineProfile* lp, int nId);

private:
    bool Validate(LayerPointSet*& spline0, LayerPointSet* &spline1);

    Ui::DialogLineProfile *ui;

    LayerLineProfile*   m_lineProfile;
};

#endif // DIALOGLINEPROFILE_H
