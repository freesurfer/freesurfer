/****************************************************************************
** Meta object code from reading C++ file 'VolumeFilterWorkerThread.h'
**
** Created: Fri Apr 6 15:12:32 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "VolumeFilterWorkerThread.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'VolumeFilterWorkerThread.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_VolumeFilterWorkerThread[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      28,   26,   25,   25, 0x05,
      49,   42,   25,   25, 0x05,

 // slots: signature, parameters, type, tag, flags
      73,   42,   25,   25, 0x0a,
     102,   25,   25,   25, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_VolumeFilterWorkerThread[] = {
    "VolumeFilterWorkerThread\0\0n\0Progress(int)\0"
    "filter\0Finished(VolumeFilter*)\0"
    "ExecuteFilter(VolumeFilter*)\0OnFinished()\0"
};

const QMetaObject VolumeFilterWorkerThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_VolumeFilterWorkerThread,
      qt_meta_data_VolumeFilterWorkerThread, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &VolumeFilterWorkerThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *VolumeFilterWorkerThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *VolumeFilterWorkerThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_VolumeFilterWorkerThread))
        return static_cast<void*>(const_cast< VolumeFilterWorkerThread*>(this));
    return QThread::qt_metacast(_clname);
}

int VolumeFilterWorkerThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: Progress((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: Finished((*reinterpret_cast< VolumeFilter*(*)>(_a[1]))); break;
        case 2: ExecuteFilter((*reinterpret_cast< VolumeFilter*(*)>(_a[1]))); break;
        case 3: OnFinished(); break;
        default: ;
        }
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void VolumeFilterWorkerThread::Progress(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void VolumeFilterWorkerThread::Finished(VolumeFilter * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_END_MOC_NAMESPACE
