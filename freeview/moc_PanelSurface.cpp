/****************************************************************************
** Meta object code from reading C++ file 'PanelSurface.h'
**
** Created: Fri Apr 6 15:12:29 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "PanelSurface.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'PanelSurface.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PanelSurface[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      18,   14,   13,   13, 0x09,
      53,   48,   13,   13, 0x09,
      74,   48,   13,   13, 0x09,
      96,   48,   13,   13, 0x09,
     120,  115,   13,   13, 0x09,
     147,  142,   13,   13, 0x09,
     175,  142,   13,   13, 0x09,
     200,  115,   13,   13, 0x09,
     220,  115,   13,   13, 0x09,
     243,  115,   13,   13, 0x09,
     261,  115,   13,   13, 0x09,
     280,   13,   13,   13, 0x09,
     307,   13,   13,   13, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_PanelSurface[] = {
    "PanelSurface\0\0act\0OnChangeSurfaceType(QAction*)\0"
    "nVal\0OnSliderOpacity(int)\0"
    "OnSliderMidPoint(int)\0OnSliderSlope(int)\0"
    "nSel\0OnComboCurvature(int)\0text\0"
    "OnLineEditMidPoint(QString)\0"
    "OnLineEditSlope(QString)\0OnComboOverlay(int)\0"
    "OnComboAnnotation(int)\0OnComboLabel(int)\0"
    "OnComboVector(int)\0OnButtonConfigureOverlay()\0"
    "OnEditPositionOffset()\0"
};

const QMetaObject PanelSurface::staticMetaObject = {
    { &PanelLayer::staticMetaObject, qt_meta_stringdata_PanelSurface,
      qt_meta_data_PanelSurface, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PanelSurface::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PanelSurface::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PanelSurface::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PanelSurface))
        return static_cast<void*>(const_cast< PanelSurface*>(this));
    return PanelLayer::qt_metacast(_clname);
}

int PanelSurface::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = PanelLayer::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnChangeSurfaceType((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 1: OnSliderOpacity((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: OnSliderMidPoint((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: OnSliderSlope((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: OnComboCurvature((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: OnLineEditMidPoint((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: OnLineEditSlope((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: OnComboOverlay((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: OnComboAnnotation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: OnComboLabel((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: OnComboVector((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: OnButtonConfigureOverlay(); break;
        case 12: OnEditPositionOffset(); break;
        default: ;
        }
        _id -= 13;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
