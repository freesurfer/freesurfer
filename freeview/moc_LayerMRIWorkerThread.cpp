/****************************************************************************
** Meta object code from reading C++ file 'LayerMRIWorkerThread.h'
**
** Created: Fri Apr 6 15:12:32 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerMRIWorkerThread.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerMRIWorkerThread.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerMRIWorkerThread[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      27,   22,   21,   21, 0x05,

 // slots: signature, parameters, type, tag, flags
      52,   21,   21,   21, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerMRIWorkerThread[] = {
    "LayerMRIWorkerThread\0\0list\0"
    "AvailableLabels(IntList)\0Abort()\0"
};

const QMetaObject LayerMRIWorkerThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_LayerMRIWorkerThread,
      qt_meta_data_LayerMRIWorkerThread, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerMRIWorkerThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerMRIWorkerThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerMRIWorkerThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerMRIWorkerThread))
        return static_cast<void*>(const_cast< LayerMRIWorkerThread*>(this));
    return QThread::qt_metacast(_clname);
}

int LayerMRIWorkerThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: AvailableLabels((*reinterpret_cast< const IntList(*)>(_a[1]))); break;
        case 1: Abort(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void LayerMRIWorkerThread::AvailableLabels(const IntList & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
