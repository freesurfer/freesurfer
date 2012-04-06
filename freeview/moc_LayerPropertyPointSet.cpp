/****************************************************************************
** Meta object code from reading C++ file 'LayerPropertyPointSet.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerPropertyPointSet.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerPropertyPointSet.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerPropertyPointSet[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      22,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       8,       // signalCount

 // signals: signature, parameters, type, tag, flags
      29,   23,   22,   22, 0x05,
      68,   60,   22,   22, 0x05,
      98,   22,   22,   22, 0x05,
     128,   22,   22,   22, 0x05,
     147,   22,   22,   22, 0x05,
     165,   22,   22,   22, 0x05,
     188,   22,   22,   22, 0x05,
     210,   22,   22,   22, 0x05,

 // slots: signature, parameters, type, tag, flags
     246,  238,   22,   22, 0x0a,
     267,  265,   22,   22, 0x0a,
     285,  265,   22,   22, 0x0a,
     315,  309,   22,   22, 0x0a,
     342,   60,   22,   22, 0x0a,
     368,  362,   22,   22, 0x0a,
     383,  381,   22,   22, 0x0a,
     408,  401,   22,   22, 0x0a,
     435,  401,   22,   22, 0x0a,
     459,  401,   22,   22, 0x0a,
     483,  401,   22,   22, 0x0a,
     517,  507,   22,   22, 0x0a,
     536,  534,   22,   22, 0x0a,
     553,  534,   22,   22, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerPropertyPointSet[] = {
    "LayerPropertyPointSet\0\0bSnag\0"
    "SnapToVoxelCenterChanged(bool)\0bSpline\0"
    "SplineVisibilityChanged(bool)\0"
    "ScalarLayerChanged(LayerMRI*)\0"
    "ScalarSetChanged()\0ColorMapChanged()\0"
    "OpacityChanged(double)\0RadiusChanged(double)\0"
    "SplineRadiusChanged(double)\0opacity\0"
    "SetOpacity(double)\0r\0SetRadius(double)\0"
    "SetSplineRadius(double)\0bSnap\0"
    "SetSnapToVoxelCenter(bool)\0"
    "SetShowSpline(bool)\0nType\0SetType(int)\0"
    "n\0SetScalarSet(int)\0iValue\0"
    "SetHeatScaleOffset(double)\0"
    "SetHeatScaleMid(double)\0SetHeatScaleMin(double)\0"
    "SetHeatScaleMax(double)\0nColorMap\0"
    "SetColorMap(int)\0c\0SetColor(QColor)\0"
    "SetSplineColor(QColor)\0"
};

const QMetaObject LayerPropertyPointSet::staticMetaObject = {
    { &LayerProperty::staticMetaObject, qt_meta_stringdata_LayerPropertyPointSet,
      qt_meta_data_LayerPropertyPointSet, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerPropertyPointSet::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerPropertyPointSet::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerPropertyPointSet::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerPropertyPointSet))
        return static_cast<void*>(const_cast< LayerPropertyPointSet*>(this));
    return LayerProperty::qt_metacast(_clname);
}

int LayerPropertyPointSet::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerProperty::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SnapToVoxelCenterChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: SplineVisibilityChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: ScalarLayerChanged((*reinterpret_cast< LayerMRI*(*)>(_a[1]))); break;
        case 3: ScalarSetChanged(); break;
        case 4: ColorMapChanged(); break;
        case 5: OpacityChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: RadiusChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 7: SplineRadiusChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 8: SetOpacity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 9: SetRadius((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 10: SetSplineRadius((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 11: SetSnapToVoxelCenter((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: SetShowSpline((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: SetType((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 14: SetScalarSet((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: SetHeatScaleOffset((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 16: SetHeatScaleMid((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 17: SetHeatScaleMin((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 18: SetHeatScaleMax((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 19: SetColorMap((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 20: SetColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 21: SetSplineColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 22;
    }
    return _id;
}

// SIGNAL 0
void LayerPropertyPointSet::SnapToVoxelCenterChanged(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void LayerPropertyPointSet::SplineVisibilityChanged(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void LayerPropertyPointSet::ScalarLayerChanged(LayerMRI * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void LayerPropertyPointSet::ScalarSetChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}

// SIGNAL 4
void LayerPropertyPointSet::ColorMapChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void LayerPropertyPointSet::OpacityChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void LayerPropertyPointSet::RadiusChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void LayerPropertyPointSet::SplineRadiusChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}
QT_END_MOC_NAMESPACE
