/****************************************************************************
** Meta object code from reading C++ file 'DialogLoadVolume.h'
**
** Created: Fri Apr 6 15:12:24 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "DialogLoadVolume.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'DialogLoadVolume.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_DialogLoadVolume[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      18,   17,   17,   17, 0x09,
      27,   17,   17,   17, 0x09,
      53,   48,   17,   17, 0x09,
      69,   48,   17,   17, 0x09,
      80,   17,   17,   17, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_DialogLoadVolume[] = {
    "DialogLoadVolume\0\0OnOpen()\0"
    "OnOpenRegistration()\0nSel\0OnColorMap(int)\0"
    "OnLUT(int)\0OnOK()\0"
};

const QMetaObject DialogLoadVolume::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_DialogLoadVolume,
      qt_meta_data_DialogLoadVolume, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DialogLoadVolume::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DialogLoadVolume::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DialogLoadVolume::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DialogLoadVolume))
        return static_cast<void*>(const_cast< DialogLoadVolume*>(this));
    return QDialog::qt_metacast(_clname);
}

int DialogLoadVolume::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnOpen(); break;
        case 1: OnOpenRegistration(); break;
        case 2: OnColorMap((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: OnLUT((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: OnOK(); break;
        default: ;
        }
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
