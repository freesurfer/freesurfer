/****************************************************************************
** Meta object code from reading C++ file 'SurfaceOverlayProperty.h'
**
** Created: Fri Apr 6 15:12:30 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "SurfaceOverlayProperty.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'SurfaceOverlayProperty.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_SurfaceOverlayProperty[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      24,   23,   23,   23, 0x05,
      42,   23,   23,   23, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_SurfaceOverlayProperty[] = {
    "SurfaceOverlayProperty\0\0ColorMapChanged()\0"
    "SmoothChanged()\0"
};

const QMetaObject SurfaceOverlayProperty::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_SurfaceOverlayProperty,
      qt_meta_data_SurfaceOverlayProperty, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &SurfaceOverlayProperty::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *SurfaceOverlayProperty::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *SurfaceOverlayProperty::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_SurfaceOverlayProperty))
        return static_cast<void*>(const_cast< SurfaceOverlayProperty*>(this));
    return QObject::qt_metacast(_clname);
}

int SurfaceOverlayProperty::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ColorMapChanged(); break;
        case 1: SmoothChanged(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void SurfaceOverlayProperty::ColorMapChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void SurfaceOverlayProperty::SmoothChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
