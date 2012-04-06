/****************************************************************************
** Meta object code from reading C++ file 'LayerPLabel.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerPLabel.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerPLabel.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerPLabel[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

static const char qt_meta_stringdata_LayerPLabel[] = {
    "LayerPLabel\0"
};

const QMetaObject LayerPLabel::staticMetaObject = {
    { &LayerMRI::staticMetaObject, qt_meta_stringdata_LayerPLabel,
      qt_meta_data_LayerPLabel, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerPLabel::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerPLabel::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerPLabel::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerPLabel))
        return static_cast<void*>(const_cast< LayerPLabel*>(this));
    return LayerMRI::qt_metacast(_clname);
}

int LayerPLabel::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerMRI::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
QT_END_MOC_NAMESPACE
