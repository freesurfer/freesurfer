/****************************************************************************
** Meta object code from reading C++ file 'Interactor2D.h'
**
** Created: Fri Apr 6 15:12:26 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "Interactor2D.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'Interactor2D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Interactor2D[] = {

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
      28,   14,   13,   13, 0x05,
      58,   50,   13,   13, 0x25,

       0        // eod
};

static const char qt_meta_stringdata_Interactor2D[] = {
    "Interactor2D\0\0message,layer\0"
    "Error(QString,Layer*)\0message\0"
    "Error(QString)\0"
};

const QMetaObject Interactor2D::staticMetaObject = {
    { &Interactor::staticMetaObject, qt_meta_stringdata_Interactor2D,
      qt_meta_data_Interactor2D, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Interactor2D::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Interactor2D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Interactor2D::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Interactor2D))
        return static_cast<void*>(const_cast< Interactor2D*>(this));
    return Interactor::qt_metacast(_clname);
}

int Interactor2D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Interactor::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: Error((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< Layer*(*)>(_a[2]))); break;
        case 1: Error((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void Interactor2D::Error(const QString & _t1, Layer * _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
