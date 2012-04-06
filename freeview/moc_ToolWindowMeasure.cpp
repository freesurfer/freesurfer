/****************************************************************************
** Meta object code from reading C++ file 'ToolWindowMeasure.h'
**
** Created: Fri Apr 6 15:12:30 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ToolWindowMeasure.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ToolWindowMeasure.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ToolWindowMeasure[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      23,   19,   18,   18, 0x0a,
      44,   18,   18,   18, 0x2a,
      56,   19,   18,   18, 0x0a,
      89,   18,   18,   18, 0x2a,
     108,   18,   18,   18, 0x09,
     121,  117,   18,   18, 0x09,
     140,   18,   18,   18, 0x09,
     156,   18,   18,   18, 0x09,
     165,   18,   18,   18, 0x09,
     174,   18,   18,   18, 0x09,
     186,   18,   18,   18, 0x09,
     197,   18,   18,   18, 0x09,
     206,   18,   18,   18, 0x09,
     221,  217,   18,   18, 0x09,
     238,  217,   18,   18, 0x09,
     264,  258,   18,   18, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_ToolWindowMeasure[] = {
    "ToolWindowMeasure\0\0reg\0SetRegion(Region2D*)\0"
    "SetRegion()\0SetSurfaceRegion(SurfaceRegion*)\0"
    "SetSurfaceRegion()\0OnIdle()\0act\0"
    "OnAction(QAction*)\0UpdateWidgets()\0"
    "OnLoad()\0OnSave()\0OnSaveAll()\0OnUpdate()\0"
    "OnCopy()\0OnExport()\0val\0OnSpinBoxId(int)\0"
    "OnSpinBoxGroup(int)\0color\0"
    "OnColorGroup(QColor)\0"
};

const QMetaObject ToolWindowMeasure::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ToolWindowMeasure,
      qt_meta_data_ToolWindowMeasure, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ToolWindowMeasure::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ToolWindowMeasure::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ToolWindowMeasure::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ToolWindowMeasure))
        return static_cast<void*>(const_cast< ToolWindowMeasure*>(this));
    if (!strcmp(_clname, "UIUpdateHelper"))
        return static_cast< UIUpdateHelper*>(const_cast< ToolWindowMeasure*>(this));
    return QWidget::qt_metacast(_clname);
}

int ToolWindowMeasure::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SetRegion((*reinterpret_cast< Region2D*(*)>(_a[1]))); break;
        case 1: SetRegion(); break;
        case 2: SetSurfaceRegion((*reinterpret_cast< SurfaceRegion*(*)>(_a[1]))); break;
        case 3: SetSurfaceRegion(); break;
        case 4: OnIdle(); break;
        case 5: OnAction((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 6: UpdateWidgets(); break;
        case 7: OnLoad(); break;
        case 8: OnSave(); break;
        case 9: OnSaveAll(); break;
        case 10: OnUpdate(); break;
        case 11: OnCopy(); break;
        case 12: OnExport(); break;
        case 13: OnSpinBoxId((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 14: OnSpinBoxGroup((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: OnColorGroup((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 16;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
