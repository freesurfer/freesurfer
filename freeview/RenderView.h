#ifndef RENDERVIEW_H
#define RENDERVIEW_H

#include "GenericRenderView.h"
#include "vtkSmartPointer.h"

class Interactor;
class QEvent;
class QMouseEvent;
class QKeyEvent;
class QWheelEvent;
class QFocusEvent;
class vtkActor2D;
class vtkScalarBarActor;

class RenderView : public GenericRenderView
{
    Q_OBJECT
public:
    RenderView( QWidget* parent = NULL );

    enum InteractionMode { IM_Navigate = 0, IM_Measure, IM_VoxelEdit, IM_ROIEdit, IM_PointSetEdit, IM_VolumeCrop };

    void SetWorldCoordinateInfo( const double* origin, const double* size );
    virtual void UpdateViewByWorldCoordinate() {}

    int GetInteractionMode();
    virtual void SetInteractionMode( int nMode );

    int GetAction();

    void SetFocusFrameColor( double r, double g, double b );

    void ViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z );
    void NormalizedViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z );
    void WorldToViewport( double world_x, double world_y, double world_z, double& x, double& y, double& z );
    void WorldToScreen( double world_x, double world_y, double world_z, int& x, int& y );
    void ScreenToWorld( int x, int y, int z, double& world_x, double& world_y, double& world_z );

    virtual void focusInEvent     ( QFocusEvent* event);
    virtual void focusOutEvent    ( QFocusEvent* event );
    virtual void mousePressEvent  ( QMouseEvent* event );
    virtual void mouseReleaseEvent  ( QMouseEvent* event );
    virtual void mouseMoveEvent   ( QMouseEvent* event );
    virtual void wheelEvent       ( QWheelEvent* event );
    virtual void enterEvent       ( QEvent* event );
    virtual void leaveEvent       ( QEvent* event );
    virtual void keyPressEvent    ( QKeyEvent* event );
    virtual void keyReleaseEvent  ( QKeyEvent* event );

    bool GetShowScalarBar();
    virtual void UpdateScalarBar();

    bool SaveScreenShot(const QString& filename, bool bAntiAliasing, int nMag = 1);

    virtual void TriggerContextMenu( QMouseEvent* event ) {}

public slots:
    void RequestRedraw( bool bForce = false );
    void MoveUp();
    void MoveDown();
    void MoveLeft();
    void MoveRight();
    void Zoom( double factor );
    void Reset();
    void SetAction( int nAction );
    void ShowScalarBar( bool bShow );

protected:
    virtual void paintEvent(QPaintEvent *event);

protected slots:
    virtual void OnIdle();
    virtual void OnSlicePositionChanged()
    {
      RequestRedraw();
    }

protected:
    bool    m_bNeedRedraw;
    double  m_dWorldOrigin[3];
    double  m_dWorldSize[3];

    Interactor*   m_interactor;
    int           m_nInteractionMode;

    vtkSmartPointer<vtkActor2D>   m_actorFocusFrame;
    vtkSmartPointer<vtkScalarBarActor>  m_actorScalarBar;
};

#endif // RENDERVIEW_H
