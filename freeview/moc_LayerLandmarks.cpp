/****************************************************************************
** Meta object code from reading C++ file 'LayerLandmarks.h'
**
** Created: Fri Apr 6 15:12:32 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerLandmarks.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerLandmarks.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerLandmarks[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      21,   16,   15,   15, 0x05,
      51,   15,   15,   15, 0x05,

 // slots: signature, parameters, type, tag, flags
      71,   67,   15,   15, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerLandmarks[] = {
    "LayerLandmarks\0\0n,lm\0LandmarkChanged(int,Landmark)\0"
    "LandmarkAdded()\0mri\0SetMRIRef(LayerMRI*)\0"
};

const QMetaObject LayerLandmarks::staticMetaObject = {
    { &LayerEditable::staticMetaObject, qt_meta_stringdata_LayerLandmarks,
      qt_meta_data_LayerLandmarks, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerLandmarks::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerLandmarks::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerLandmarks::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerLandmarks))
        return static_cast<void*>(const_cast< LayerLandmarks*>(this));
    return LayerEditable::qt_metacast(_clname);
}

int LayerLandmarks::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerEditable::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: LandmarkChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< const Landmark(*)>(_a[2]))); break;
        case 1: LandmarkAdded(); break;
        case 2: SetMRIRef((*reinterpret_cast< LayerMRI*(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void LayerLandmarks::LandmarkChanged(int _t1, const Landmark & _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void LayerLandmarks::LandmarkAdded()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
