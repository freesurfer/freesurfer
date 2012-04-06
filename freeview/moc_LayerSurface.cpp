/****************************************************************************
** Meta object code from reading C++ file 'LayerSurface.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerSurface.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerSurface.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerSurface[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      27,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      10,       // signalCount

 // signals: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x05,
      57,   13,   13,   13, 0x05,
      90,   13,   13,   13, 0x05,
     127,   13,   13,   13, 0x05,
     154,   13,   13,   13, 0x05,
     179,   13,   13,   13, 0x05,
     203,  201,   13,   13, 0x05,
     229,  201,   13,   13, 0x05,
     255,  201,   13,   13, 0x05,
     284,  201,   13,   13, 0x05,

 // slots: signature, parameters, type, tag, flags
     321,  308,   13,   13, 0x0a,
     354,  343,   13,   13, 0x0a,
     374,   13,   13,   13, 0x2a,
     399,  390,   13,   13, 0x0a,
     426,  424,   13,   13, 0x0a,
     463,  454,   13,   13, 0x0a,
     491,  454,   13,   13, 0x0a,
     524,   13,   13,   13, 0x09,
     540,   13,   13,   13, 0x09,
     557,   13,   13,   13, 0x09,
     579,   13,   13,   13, 0x09,
     603,   13,   13,   13, 0x09,
     622,   13,   13,   13, 0x09,
     643,   13,   13,   13, 0x09,
     662,   13,   13,   13, 0x09,
     694,  685,   13,   13, 0x09,
     734,   13,   13,   13, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_LayerSurface[] = {
    "LayerSurface\0\0SurfaceAnnotationAdded(SurfaceAnnotation*)\0"
    "SurfaceLabelAdded(SurfaceLabel*)\0"
    "SurfaceOverlayAdded(SurfaceOverlay*)\0"
    "SurfaceOverlyDataUpdated()\0"
    "SurfaceCurvatureLoaded()\0SurfaceVectorLoaded()\0"
    "n\0ActiveSurfaceChanged(int)\0"
    "ActiveOverlayChanged(int)\0"
    "ActiveAnnotationChanged(int)\0"
    "ActiveLabelChanged(int)\0nSurfaceType\0"
    "SetActiveSurface(int)\0bAskRedraw\0"
    "UpdateOverlay(bool)\0UpdateOverlay()\0"
    "bLoadAll\0SetLoadAllSurfaces(bool)\0c\0"
    "SetActiveLabelColor(QColor)\0bOutline\0"
    "SetActiveLabelOutline(bool)\0"
    "SetActiveAnnotationOutline(bool)\0"
    "UpdateOpacity()\0UpdateColorMap()\0"
    "UpdateEdgeThickness()\0UpdateVectorPointSize()\0"
    "UpdateRenderMode()\0UpdateVertexRender()\0"
    "UpdateMeshRender()\0UpdateActorPositions()\0"
    "dx,dy,dz\0UpdateROIPosition(double,double,double)\0"
    "UpdateVectorActor2D()\0"
};

const QMetaObject LayerSurface::staticMetaObject = {
    { &LayerEditable::staticMetaObject, qt_meta_stringdata_LayerSurface,
      qt_meta_data_LayerSurface, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerSurface::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerSurface::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerSurface::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerSurface))
        return static_cast<void*>(const_cast< LayerSurface*>(this));
    return LayerEditable::qt_metacast(_clname);
}

int LayerSurface::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerEditable::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SurfaceAnnotationAdded((*reinterpret_cast< SurfaceAnnotation*(*)>(_a[1]))); break;
        case 1: SurfaceLabelAdded((*reinterpret_cast< SurfaceLabel*(*)>(_a[1]))); break;
        case 2: SurfaceOverlayAdded((*reinterpret_cast< SurfaceOverlay*(*)>(_a[1]))); break;
        case 3: SurfaceOverlyDataUpdated(); break;
        case 4: SurfaceCurvatureLoaded(); break;
        case 5: SurfaceVectorLoaded(); break;
        case 6: ActiveSurfaceChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: ActiveOverlayChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: ActiveAnnotationChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: ActiveLabelChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: SetActiveSurface((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: UpdateOverlay((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: UpdateOverlay(); break;
        case 13: SetLoadAllSurfaces((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: SetActiveLabelColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 15: SetActiveLabelOutline((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: SetActiveAnnotationOutline((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: UpdateOpacity(); break;
        case 18: UpdateColorMap(); break;
        case 19: UpdateEdgeThickness(); break;
        case 20: UpdateVectorPointSize(); break;
        case 21: UpdateRenderMode(); break;
        case 22: UpdateVertexRender(); break;
        case 23: UpdateMeshRender(); break;
        case 24: UpdateActorPositions(); break;
        case 25: UpdateROIPosition((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 26: UpdateVectorActor2D(); break;
        default: ;
        }
        _id -= 27;
    }
    return _id;
}

// SIGNAL 0
void LayerSurface::SurfaceAnnotationAdded(SurfaceAnnotation * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void LayerSurface::SurfaceLabelAdded(SurfaceLabel * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void LayerSurface::SurfaceOverlayAdded(SurfaceOverlay * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void LayerSurface::SurfaceOverlyDataUpdated()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}

// SIGNAL 4
void LayerSurface::SurfaceCurvatureLoaded()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void LayerSurface::SurfaceVectorLoaded()
{
    QMetaObject::activate(this, &staticMetaObject, 5, 0);
}

// SIGNAL 6
void LayerSurface::ActiveSurfaceChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void LayerSurface::ActiveOverlayChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}

// SIGNAL 8
void LayerSurface::ActiveAnnotationChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 8, _a);
}

// SIGNAL 9
void LayerSurface::ActiveLabelChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 9, _a);
}
QT_END_MOC_NAMESPACE
