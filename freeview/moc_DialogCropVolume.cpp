/****************************************************************************
** Meta object code from reading C++ file 'DialogCropVolume.h'
**
** Created: Fri Apr 6 15:12:24 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "DialogCropVolume.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'DialogCropVolume.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_DialogCropVolume[] = {

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
      22,   18,   17,   17, 0x09,
      58,   52,   17,   17, 0x09,
      86,   81,   17,   17, 0x09,
     103,   17,   17,   17, 0x09,
     126,   17,   17,   17, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_DialogCropVolume[] = {
    "DialogCropVolume\0\0mri\0"
    "OnCropBoundChanged(LayerMRI*)\0layer\0"
    "OnLayerRemoved(Layer*)\0nVal\0"
    "OnSpinRange(int)\0showEvent(QShowEvent*)\0"
    "hideEvent(QHideEvent*)\0"
};

const QMetaObject DialogCropVolume::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_DialogCropVolume,
      qt_meta_data_DialogCropVolume, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DialogCropVolume::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DialogCropVolume::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DialogCropVolume::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DialogCropVolume))
        return static_cast<void*>(const_cast< DialogCropVolume*>(this));
    return QDialog::qt_metacast(_clname);
}

int DialogCropVolume::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnCropBoundChanged((*reinterpret_cast< LayerMRI*(*)>(_a[1]))); break;
        case 1: OnLayerRemoved((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 2: OnSpinRange((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: showEvent((*reinterpret_cast< QShowEvent*(*)>(_a[1]))); break;
        case 4: hideEvent((*reinterpret_cast< QHideEvent*(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
