/****************************************************************************
** Meta object code from reading C++ file 'DialogTransformVolume.h'
**
** Created: Fri Apr 6 15:12:25 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "DialogTransformVolume.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'DialogTransformVolume.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_DialogTransformVolume[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      21,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      25,   23,   22,   22, 0x05,

 // slots: signature, parameters, type, tag, flags
      53,   22,   22,   22, 0x09,
      63,   22,   22,   22, 0x09,
      75,   22,   22,   22, 0x09,
      92,   87,   22,   22, 0x09,
     119,   87,   22,   22, 0x09,
     146,   87,   22,   22, 0x09,
     178,  173,   22,   22, 0x09,
     208,  173,   22,   22, 0x09,
     238,  173,   22,   22, 0x09,
     268,   87,   22,   22, 0x09,
     291,   87,   22,   22, 0x09,
     314,   87,   22,   22, 0x09,
     337,  173,   22,   22, 0x09,
     363,  173,   22,   22, 0x09,
     389,  173,   22,   22, 0x09,
     415,   22,   22,   22, 0x09,
     439,   22,   22,   22, 0x09,
     471,  462,   22,   22, 0x09,
     499,   22,   22,   22, 0x09,
     522,   22,   22,   22, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_DialogTransformVolume[] = {
    "DialogTransformVolume\0\0n\0"
    "CurrentLandmarkChanged(int)\0OnApply()\0"
    "OnRestore()\0OnSaveReg()\0nVal\0"
    "OnScrollBarTranslateX(int)\0"
    "OnScrollBarTranslateY(int)\0"
    "OnScrollBarTranslateZ(int)\0text\0"
    "OnLineEditTranslateX(QString)\0"
    "OnLineEditTranslateY(QString)\0"
    "OnLineEditTranslateZ(QString)\0"
    "OnScrollBarScaleX(int)\0OnScrollBarScaleY(int)\0"
    "OnScrollBarScaleZ(int)\0OnLineEditScaleX(QString)\0"
    "OnLineEditScaleY(QString)\0"
    "OnLineEditScaleZ(QString)\0"
    "OnSampleMethodChanged()\0OnActiveLayerChanged()\0"
    "bChecked\0OnRadioButtonLandmark(bool)\0"
    "OnButtonLandmarkPick()\0UpdateLandmarkColors()\0"
};

const QMetaObject DialogTransformVolume::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_DialogTransformVolume,
      qt_meta_data_DialogTransformVolume, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DialogTransformVolume::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DialogTransformVolume::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DialogTransformVolume::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DialogTransformVolume))
        return static_cast<void*>(const_cast< DialogTransformVolume*>(this));
    if (!strcmp(_clname, "UIUpdateHelper"))
        return static_cast< UIUpdateHelper*>(const_cast< DialogTransformVolume*>(this));
    return QDialog::qt_metacast(_clname);
}

int DialogTransformVolume::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: CurrentLandmarkChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: OnApply(); break;
        case 2: OnRestore(); break;
        case 3: OnSaveReg(); break;
        case 4: OnScrollBarTranslateX((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: OnScrollBarTranslateY((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: OnScrollBarTranslateZ((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: OnLineEditTranslateX((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 8: OnLineEditTranslateY((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 9: OnLineEditTranslateZ((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 10: OnScrollBarScaleX((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: OnScrollBarScaleY((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: OnScrollBarScaleZ((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: OnLineEditScaleX((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 14: OnLineEditScaleY((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 15: OnLineEditScaleZ((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 16: OnSampleMethodChanged(); break;
        case 17: OnActiveLayerChanged(); break;
        case 18: OnRadioButtonLandmark((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 19: OnButtonLandmarkPick(); break;
        case 20: UpdateLandmarkColors(); break;
        default: ;
        }
        _id -= 21;
    }
    return _id;
}

// SIGNAL 0
void DialogTransformVolume::CurrentLandmarkChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
