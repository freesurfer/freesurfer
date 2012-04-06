/****************************************************************************
** Meta object code from reading C++ file 'DialogLoadDTI.h'
**
** Created: Fri Apr 6 15:12:24 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "DialogLoadDTI.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'DialogLoadDTI.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_DialogLoadDTI[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x09,
      22,   14,   14,   14, 0x09,
      39,   14,   14,   14, 0x09,
      52,   14,   14,   14, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_DialogLoadDTI[] = {
    "DialogLoadDTI\0\0OnOK()\0OnButtonVector()\0"
    "OnButtonFA()\0OnButtonRegistration()\0"
};

const QMetaObject DialogLoadDTI::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_DialogLoadDTI,
      qt_meta_data_DialogLoadDTI, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DialogLoadDTI::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DialogLoadDTI::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DialogLoadDTI::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DialogLoadDTI))
        return static_cast<void*>(const_cast< DialogLoadDTI*>(this));
    return QDialog::qt_metacast(_clname);
}

int DialogLoadDTI::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnOK(); break;
        case 1: OnButtonVector(); break;
        case 2: OnButtonFA(); break;
        case 3: OnButtonRegistration(); break;
        default: ;
        }
        _id -= 4;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
