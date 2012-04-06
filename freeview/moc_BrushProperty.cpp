/****************************************************************************
** Meta object code from reading C++ file 'BrushProperty.h'
**
** Created: Fri Apr 6 15:12:23 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "BrushProperty.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'BrushProperty.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_BrushProperty[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      21,   15,   14,   14, 0x0a,
      50,   39,   14,   14, 0x0a,
      79,   73,   14,   14, 0x0a,
     123,  115,   14,   14, 0x0a,
     149,  115,   14,   14, 0x0a,
     178,  115,   14,   14, 0x0a,
     205,   73,   14,   14, 0x0a,
     233,  228,   14,   14, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_BrushProperty[] = {
    "BrushProperty\0\0nSize\0SetBrushSize(int)\0"
    "nTolerance\0SetBrushTolerance(int)\0"
    "layer\0SetReferenceLayer(LayerVolumeBase*)\0"
    "bEnable\0SetDrawRangeEnabled(bool)\0"
    "SetExcludeRangeEnabled(bool)\0"
    "SetDrawConnectedOnly(bool)\0"
    "OnLayerRemoved(Layer*)\0bVal\0SetFill3D(bool)\0"
};

const QMetaObject BrushProperty::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_BrushProperty,
      qt_meta_data_BrushProperty, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &BrushProperty::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *BrushProperty::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *BrushProperty::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_BrushProperty))
        return static_cast<void*>(const_cast< BrushProperty*>(this));
    return QObject::qt_metacast(_clname);
}

int BrushProperty::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SetBrushSize((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: SetBrushTolerance((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: SetReferenceLayer((*reinterpret_cast< LayerVolumeBase*(*)>(_a[1]))); break;
        case 3: SetDrawRangeEnabled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 4: SetExcludeRangeEnabled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: SetDrawConnectedOnly((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 6: OnLayerRemoved((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 7: SetFill3D((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
