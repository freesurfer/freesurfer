/****************************************************************************
** Meta object code from reading C++ file 'LayerPropertyTrack.h'
**
** Created: Fri Apr 6 15:12:31 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerPropertyTrack.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerPropertyTrack.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerPropertyTrack[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       5,       // signalCount

 // signals: signature, parameters, type, tag, flags
      20,   19,   19,   19, 0x05,
      42,   19,   19,   19, 0x05,
      70,   19,   19,   19, 0x05,
     101,   99,   19,   19, 0x05,
     127,   19,   19,   19, 0x05,

 // slots: signature, parameters, type, tag, flags
     152,  146,   19,   19, 0x0a,
     175,  170,   19,   19, 0x0a,
     199,  170,   19,   19, 0x0a,
     224,   99,   19,   19, 0x0a,
     246,  170,   19,   19, 0x0a,
     269,  264,   19,   19, 0x0a,
     291,  170,   19,   19, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerPropertyTrack[] = {
    "LayerPropertyTrack\0\0ColorCodeChanged(int)\0"
    "DirectionSchemeChanged(int)\0"
    "DirectionMappingChanged(int)\0c\0"
    "SolidColorChanged(QColor)\0RenderRepChanged()\0"
    "nCode\0SetColorCode(int)\0nVal\0"
    "SetDirectionScheme(int)\0"
    "SetDirectionMapping(int)\0SetSolidColor(QColor)\0"
    "SetRenderRep(int)\0dVal\0SetTubeRadius(double)\0"
    "SetNumberOfSides(int)\0"
};

const QMetaObject LayerPropertyTrack::staticMetaObject = {
    { &LayerProperty::staticMetaObject, qt_meta_stringdata_LayerPropertyTrack,
      qt_meta_data_LayerPropertyTrack, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerPropertyTrack::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerPropertyTrack::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerPropertyTrack::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerPropertyTrack))
        return static_cast<void*>(const_cast< LayerPropertyTrack*>(this));
    return LayerProperty::qt_metacast(_clname);
}

int LayerPropertyTrack::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerProperty::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ColorCodeChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: DirectionSchemeChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: DirectionMappingChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: SolidColorChanged((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 4: RenderRepChanged(); break;
        case 5: SetColorCode((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: SetDirectionScheme((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: SetDirectionMapping((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: SetSolidColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 9: SetRenderRep((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: SetTubeRadius((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 11: SetNumberOfSides((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void LayerPropertyTrack::ColorCodeChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void LayerPropertyTrack::DirectionSchemeChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void LayerPropertyTrack::DirectionMappingChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void LayerPropertyTrack::SolidColorChanged(const QColor & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void LayerPropertyTrack::RenderRepChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}
QT_END_MOC_NAMESPACE
