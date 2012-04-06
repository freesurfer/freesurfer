/****************************************************************************
** Meta object code from reading C++ file 'VolumeCropper.h'
**
** Created: Fri Apr 6 15:12:31 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "VolumeCropper.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'VolumeCropper.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_VolumeCropper[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x05,

 // slots: signature, parameters, type, tag, flags
      43,   14,   14,   14, 0x0a,
      51,   14,   14,   14, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_VolumeCropper[] = {
    "VolumeCropper\0\0CropBoundChanged(LayerMRI*)\0"
    "Reset()\0Apply()\0"
};

const QMetaObject VolumeCropper::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_VolumeCropper,
      qt_meta_data_VolumeCropper, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &VolumeCropper::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *VolumeCropper::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *VolumeCropper::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_VolumeCropper))
        return static_cast<void*>(const_cast< VolumeCropper*>(this));
    return QObject::qt_metacast(_clname);
}

int VolumeCropper::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: CropBoundChanged((*reinterpret_cast< LayerMRI*(*)>(_a[1]))); break;
        case 1: Reset(); break;
        case 2: Apply(); break;
        default: ;
        }
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void VolumeCropper::CropBoundChanged(LayerMRI * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
