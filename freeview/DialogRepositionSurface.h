#ifndef DIALOGREPOSITIONSURFACE_H
#define DIALOGREPOSITIONSURFACE_H

#include <QDialog>

namespace Ui {
    class DialogRepositionSurface;
}

class DialogRepositionSurface : public QDialog
{
    Q_OBJECT

  public:
    explicit DialogRepositionSurface(QWidget *parent = 0);
    ~DialogRepositionSurface();

    int GetVertex();
    int GetNeighborSize();

    double GetIntensity();
    double GetSigma();

    void GetCoordinate( double* pos );

  public slots:
    void OnApply();
    void OnUndo();
    void OnSave();
    void OnSaveAs();
    void OnComboTarget(int n);
    void UpdateUI();

    void OnSurfaceVertexClicked();

  private:
    bool ValidateAll();

    Ui::DialogRepositionSurface *ui;
};

#endif // DIALOGREPOSITIONSURFACE_H
