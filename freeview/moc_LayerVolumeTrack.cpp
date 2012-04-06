/****************************************************************************
** Meta object code from reading C++ file 'LayerVolumeTrack.h'
**
** Created: Fri Apr 6 15:12:32 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerVolumeTrack.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerVolumeTrack.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerVolumeTrack[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      25,   18,   17,   17, 0x0a,
      40,   17,   17,   17, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerVolumeTrack[] = {
    "LayerVolumeTrack\0\0nLabel\0Highlight(int)\0"
    "RestoreColors()\0"
};

const QMetaObject LayerVolumeTrack::staticMetaObject = {
    { &LayerMRI::staticMetaObject, qt_meta_stringdata_LayerVolumeTrack,
      qt_meta_data_LayerVolumeTrack, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerVolumeTrack::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerVolumeTrack::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerVolumeTrack::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerVolumeTrack))
        return static_cast<void*>(const_cast< LayerVolumeTrack*>(this));
    return LayerMRI::qt_metacast(_clname);
}

int LayerVolumeTrack::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerMRI::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: Highlight((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: RestoreColors(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
