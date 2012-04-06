/****************************************************************************
** Meta object code from reading C++ file 'PanelPointSet.h'
**
** Created: Fri Apr 6 15:12:28 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "PanelPointSet.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'PanelPointSet.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PanelPointSet[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      20,   15,   14,   14, 0x09,
      41,   15,   14,   14, 0x09,
      58,   15,   14,   14, 0x09,
      75,   15,   14,   14, 0x09,
      92,   15,   14,   14, 0x09,
     117,  112,   14,   14, 0x09,
     140,  112,   14,   14, 0x09,
     163,  112,   14,   14, 0x09,
     186,  112,   14,   14, 0x09,
     212,  112,   14,   14, 0x09,
     238,  112,   14,   14, 0x09,
     275,  270,   14,   14, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_PanelPointSet[] = {
    "PanelPointSet\0\0nVal\0OnSliderOpacity(int)\0"
    "OnSliderMin(int)\0OnSliderMid(int)\0"
    "OnSliderMax(int)\0OnSliderOffset(int)\0"
    "text\0OnLineEditMin(QString)\0"
    "OnLineEditMid(QString)\0OnLineEditMax(QString)\0"
    "OnLineEditOffset(QString)\0"
    "OnLineEditRadius(QString)\0"
    "OnLineEditSplineRadius(QString)\0nSel\0"
    "OnComboScalarMap(int)\0"
};

const QMetaObject PanelPointSet::staticMetaObject = {
    { &PanelLayer::staticMetaObject, qt_meta_stringdata_PanelPointSet,
      qt_meta_data_PanelPointSet, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PanelPointSet::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PanelPointSet::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PanelPointSet::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PanelPointSet))
        return static_cast<void*>(const_cast< PanelPointSet*>(this));
    return PanelLayer::qt_metacast(_clname);
}

int PanelPointSet::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = PanelLayer::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnSliderOpacity((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: OnSliderMin((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: OnSliderMid((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: OnSliderMax((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: OnSliderOffset((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: OnLineEditMin((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: OnLineEditMid((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: OnLineEditMax((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 8: OnLineEditOffset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 9: OnLineEditRadius((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 10: OnLineEditSplineRadius((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 11: OnComboScalarMap((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 12;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
