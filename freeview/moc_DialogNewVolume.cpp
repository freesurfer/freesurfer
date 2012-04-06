/****************************************************************************
** Meta object code from reading C++ file 'DialogNewVolume.h'
**
** Created: Fri Apr 6 15:12:25 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "DialogNewVolume.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'DialogNewVolume.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_DialogNewVolume[] = {

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
      17,   16,   16,   16, 0x09,
      30,   24,   16,   16, 0x09,
      67,   58,   16,   16, 0x09,
     103,   97,   16,   16, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_DialogNewVolume[] = {
    "DialogNewVolume\0\0OnOK()\0bCopy\0"
    "OnToggleCopyVoxelData(bool)\0bChecked\0"
    "OnToggleVoxelDataOption(bool)\0bMask\0"
    "OnToggleMask(bool)\0"
};

const QMetaObject DialogNewVolume::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_DialogNewVolume,
      qt_meta_data_DialogNewVolume, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DialogNewVolume::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DialogNewVolume::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DialogNewVolume::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DialogNewVolume))
        return static_cast<void*>(const_cast< DialogNewVolume*>(this));
    return QDialog::qt_metacast(_clname);
}

int DialogNewVolume::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnOK(); break;
        case 1: OnToggleCopyVoxelData((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: OnToggleVoxelDataOption((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: OnToggleMask((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 4;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
