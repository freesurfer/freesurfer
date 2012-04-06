/****************************************************************************
** Meta object code from reading C++ file 'LayerProperty.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerProperty.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerProperty.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerProperty[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: signature, parameters, type, tag, flags
      21,   15,   14,   14, 0x05,
      43,   14,   14,   14, 0x05,
      61,   14,   14,   14, 0x05,

 // slots: signature, parameters, type, tag, flags
      82,   15,   14,   14, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerProperty[] = {
    "LayerProperty\0\0bShow\0ShowInfoChanged(bool)\0"
    "PropertyChanged()\0DisplayModeChanged()\0"
    "SetShowInfo(bool)\0"
};

const QMetaObject LayerProperty::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_LayerProperty,
      qt_meta_data_LayerProperty, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerProperty::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerProperty::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerProperty::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerProperty))
        return static_cast<void*>(const_cast< LayerProperty*>(this));
    return QObject::qt_metacast(_clname);
}

int LayerProperty::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ShowInfoChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: PropertyChanged(); break;
        case 2: DisplayModeChanged(); break;
        case 3: SetShowInfo((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void LayerProperty::ShowInfoChanged(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void LayerProperty::PropertyChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void LayerProperty::DisplayModeChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}
QT_END_MOC_NAMESPACE
