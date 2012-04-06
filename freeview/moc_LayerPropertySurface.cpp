/****************************************************************************
** Meta object code from reading C++ file 'LayerPropertySurface.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerPropertySurface.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerPropertySurface.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerPropertySurface[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      29,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       9,       // signalCount

 // signals: signature, parameters, type, tag, flags
      30,   22,   21,   21, 0x05,
      64,   53,   21,   21, 0x05,
      96,   90,   21,   21, 0x05,
     130,  124,   21,   21, 0x05,
     153,   21,   21,   21, 0x05,
     175,   21,   21,   21, 0x05,
     195,   21,   21,   21, 0x05,
     213,   21,   21,   21, 0x05,
     240,  231,   21,   21, 0x05,

 // slots: signature, parameters, type, tag, flags
     278,   22,   21,   21, 0x0a,
     302,  297,   21,   21, 0x0a,
     323,  124,   21,   21, 0x0a,
     355,  349,   21,   21, 0x0a,
     394,  392,   21,   21, 0x0a,
     424,  417,   21,   21, 0x0a,
     453,  417,   21,   21, 0x0a,
     479,  349,   21,   21, 0x0a,
     514,  392,   21,   21, 0x0a,
     535,  349,   21,   21, 0x0a,
     572,  392,   21,   21, 0x0a,
     595,  349,   21,   21, 0x0a,
     632,  392,   21,   21, 0x0a,
     655,  349,   21,   21, 0x0a,
     690,  392,   21,   21, 0x0a,
     711,   90,   21,   21, 0x0a,
     735,  297,   21,   21, 0x0a,
     762,  756,   21,   21, 0x0a,
     781,   53,   21,   21, 0x0a,
     803,   90,   21,   21, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerPropertySurface[] = {
    "LayerPropertySurface\0\0opacity\0"
    "OpacityChanged(double)\0nThickness\0"
    "EdgeThicknessChanged(int)\0nSize\0"
    "VectorPointSizeChanged(int)\0nMode\0"
    "RenderModeChanged(int)\0VertexRenderChanged()\0"
    "MeshRenderChanged()\0ColorMapChanged()\0"
    "PositionChanged()\0dx,dy,dz\0"
    "PositionChanged(double,double,double)\0"
    "SetOpacity(double)\0nMap\0SetCurvatureMap(int)\0"
    "SetSurfaceRenderMode(int)\0r,g,b\0"
    "SetBinaryColor(double,double,double)\0"
    "c\0SetBinaryColor(QColor)\0dvalue\0"
    "SetThresholdMidPoint(double)\0"
    "SetThresholdSlope(double)\0"
    "SetEdgeColor(double,double,double)\0"
    "SetEdgeColor(QColor)\0"
    "SetVectorColor(double,double,double)\0"
    "SetVectorColor(QColor)\0"
    "SetVertexColor(double,double,double)\0"
    "SetVertexColor(QColor)\0"
    "SetMeshColor(double,double,double)\0"
    "SetMeshColor(QColor)\0SetVertexPointSize(int)\0"
    "SetMeshColorMap(int)\0bShow\0"
    "ShowVertices(bool)\0SetEdgeThickness(int)\0"
    "SetVectorPointSize(int)\0"
};

const QMetaObject LayerPropertySurface::staticMetaObject = {
    { &LayerProperty::staticMetaObject, qt_meta_stringdata_LayerPropertySurface,
      qt_meta_data_LayerPropertySurface, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerPropertySurface::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerPropertySurface::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerPropertySurface::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerPropertySurface))
        return static_cast<void*>(const_cast< LayerPropertySurface*>(this));
    return LayerProperty::qt_metacast(_clname);
}

int LayerPropertySurface::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerProperty::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OpacityChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 1: EdgeThicknessChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: VectorPointSizeChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: RenderModeChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: VertexRenderChanged(); break;
        case 5: MeshRenderChanged(); break;
        case 6: ColorMapChanged(); break;
        case 7: PositionChanged(); break;
        case 8: PositionChanged((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 9: SetOpacity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 10: SetCurvatureMap((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: SetSurfaceRenderMode((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: SetBinaryColor((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 13: SetBinaryColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 14: SetThresholdMidPoint((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 15: SetThresholdSlope((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 16: SetEdgeColor((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 17: SetEdgeColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 18: SetVectorColor((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 19: SetVectorColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 20: SetVertexColor((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 21: SetVertexColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 22: SetMeshColor((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 23: SetMeshColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 24: SetVertexPointSize((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 25: SetMeshColorMap((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 26: ShowVertices((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 27: SetEdgeThickness((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 28: SetVectorPointSize((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 29;
    }
    return _id;
}

// SIGNAL 0
void LayerPropertySurface::OpacityChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void LayerPropertySurface::EdgeThicknessChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void LayerPropertySurface::VectorPointSizeChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void LayerPropertySurface::RenderModeChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void LayerPropertySurface::VertexRenderChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void LayerPropertySurface::MeshRenderChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 5, 0);
}

// SIGNAL 6
void LayerPropertySurface::ColorMapChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 6, 0);
}

// SIGNAL 7
void LayerPropertySurface::PositionChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 7, 0);
}

// SIGNAL 8
void LayerPropertySurface::PositionChanged(double _t1, double _t2, double _t3)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)) };
    QMetaObject::activate(this, &staticMetaObject, 8, _a);
}
QT_END_MOC_NAMESPACE
