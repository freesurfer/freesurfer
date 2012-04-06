/****************************************************************************
** Meta object code from reading C++ file 'LayerMRI.h'
**
** Created: Fri Apr 6 15:12:28 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerMRI.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerMRI.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerMRI[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      30,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       7,       // signalCount

 // signals: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x05,
      41,   34,    9,    9, 0x05,
      65,    9,    9,    9, 0x05,
      86,    9,    9,    9, 0x05,
     109,    9,    9,    9, 0x05,
     132,    9,    9,    9, 0x05,
     152,    9,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
     170,   34,    9,    9, 0x0a,
     190,   34,    9,    9, 0x0a,
     217,    9,    9,    9, 0x09,
     237,    9,    9,    9, 0x09,
     253,    9,    9,    9, 0x09,
     282,    9,    9,    9, 0x09,
     317,  307,    9,    9, 0x09,
     336,    9,    9,    9, 0x29,
     352,  307,    9,    9, 0x09,
     376,    9,    9,    9, 0x09,
     397,    9,    9,    9, 0x09,
     411,    9,    9,    9, 0x09,
     435,    9,    9,    9, 0x09,
     472,  455,    9,    9, 0x09,
     516,  509,    9,    9, 0x09,
     539,    9,    9,    9, 0x09,
     563,    9,    9,    9, 0x09,
     584,    9,    9,    9, 0x09,
     607,    9,    9,    9, 0x09,
     629,    9,    9,    9, 0x09,
     649,    9,    9,    9, 0x09,
     676,  666,    9,    9, 0x09,
     710,  705,    9,    9, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_LayerMRI[] = {
    "LayerMRI\0\0ResampleFactorChanged()\0"
    "nFrame\0ActiveFrameChanged(int)\0"
    "SurfaceRegionAdded()\0SurfaceRegionUpdated()\0"
    "SurfaceRegionRemoved()\0IsoSurfaceUpdated()\0"
    "LabelStatsReady()\0SetActiveFrame(int)\0"
    "SetActiveFrameOneBase(int)\0"
    "UpdateDisplayMode()\0UpdateOpacity()\0"
    "UpdateResliceInterpolation()\0"
    "UpdateTextureSmoothing()\0nSegIndex\0"
    "UpdateContour(int)\0UpdateContour()\0"
    "UpdateContourActor(int)\0UpdateContourColor()\0"
    "ShowContour()\0UpdateVolumeRendering()\0"
    "UpdateVectorActor()\0nPlane,imagedata\0"
    "UpdateVectorActor(int,vtkImageData*)\0"
    "nPlane\0UpdateVectorActor(int)\0"
    "ResetSurfaceRegionIds()\0UpdateLabelOutline()\0"
    "UpdateUpSampleMethod()\0UpdateProjectionMap()\0"
    "UpdateTensorActor()\0UpdateColorMap()\0"
    "thread_id\0OnContourThreadFinished(int)\0"
    "vals\0OnAvailableLabels(IntList)\0"
};

const QMetaObject LayerMRI::staticMetaObject = {
    { &LayerVolumeBase::staticMetaObject, qt_meta_stringdata_LayerMRI,
      qt_meta_data_LayerMRI, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerMRI::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerMRI::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerMRI::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerMRI))
        return static_cast<void*>(const_cast< LayerMRI*>(this));
    return LayerVolumeBase::qt_metacast(_clname);
}

int LayerMRI::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerVolumeBase::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ResampleFactorChanged(); break;
        case 1: ActiveFrameChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: SurfaceRegionAdded(); break;
        case 3: SurfaceRegionUpdated(); break;
        case 4: SurfaceRegionRemoved(); break;
        case 5: IsoSurfaceUpdated(); break;
        case 6: LabelStatsReady(); break;
        case 7: SetActiveFrame((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: SetActiveFrameOneBase((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: UpdateDisplayMode(); break;
        case 10: UpdateOpacity(); break;
        case 11: UpdateResliceInterpolation(); break;
        case 12: UpdateTextureSmoothing(); break;
        case 13: UpdateContour((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 14: UpdateContour(); break;
        case 15: UpdateContourActor((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 16: UpdateContourColor(); break;
        case 17: ShowContour(); break;
        case 18: UpdateVolumeRendering(); break;
        case 19: UpdateVectorActor(); break;
        case 20: UpdateVectorActor((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< vtkImageData*(*)>(_a[2]))); break;
        case 21: UpdateVectorActor((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 22: ResetSurfaceRegionIds(); break;
        case 23: UpdateLabelOutline(); break;
        case 24: UpdateUpSampleMethod(); break;
        case 25: UpdateProjectionMap(); break;
        case 26: UpdateTensorActor(); break;
        case 27: UpdateColorMap(); break;
        case 28: OnContourThreadFinished((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 29: OnAvailableLabels((*reinterpret_cast< const IntList(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 30;
    }
    return _id;
}

// SIGNAL 0
void LayerMRI::ResampleFactorChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void LayerMRI::ActiveFrameChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void LayerMRI::SurfaceRegionAdded()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}

// SIGNAL 3
void LayerMRI::SurfaceRegionUpdated()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}

// SIGNAL 4
void LayerMRI::SurfaceRegionRemoved()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void LayerMRI::IsoSurfaceUpdated()
{
    QMetaObject::activate(this, &staticMetaObject, 5, 0);
}

// SIGNAL 6
void LayerMRI::LabelStatsReady()
{
    QMetaObject::activate(this, &staticMetaObject, 6, 0);
}
QT_END_MOC_NAMESPACE
