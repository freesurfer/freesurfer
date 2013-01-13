#ifndef LAYERPROPERTYLINEPROFILE_H
#define LAYERPROPERTYLINEPROFILE_H

#include "LayerProperty.h"
#include <QColor>

class LayerPropertyLineProfile : public LayerProperty
{
    Q_OBJECT
public:
    explicit LayerPropertyLineProfile(QObject *parent = 0);

    double GetOpacity()
    {
      return m_dOpacity;
    }

    QColor GetColor()
    {
      return m_color;
    }

    double GetRadius()
    {
      return m_dRadius;
    }

signals:
    void OpacityChanged( double );
    void ColorChanged( const QColor& color );
    void RadiusChanged(double);

public slots:
    void SetColor   (const QColor& color);
    void SetOpacity (double val);
    void SetRadius  (double val);

private:
    QColor  m_color;
    double  m_dOpacity;
    double  m_dRadius;
};

#endif // LAYERPROPERTYLINEPROFILE_H
