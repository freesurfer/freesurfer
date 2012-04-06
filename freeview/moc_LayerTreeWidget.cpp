/****************************************************************************
** Meta object code from reading C++ file 'LayerTreeWidget.h'
**
** Created: Fri Apr 6 15:12:28 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerTreeWidget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerTreeWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerTreeWidget[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      17,   16,   16,   16, 0x0a,
      31,   16,   16,   16, 0x0a,
      43,   16,   16,   16, 0x0a,
      55,   16,   16,   16, 0x0a,
      67,   16,   16,   16, 0x0a,
      81,   16,   16,   16, 0x0a,
      97,   16,   16,   16, 0x0a,
     113,   16,   16,   16, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerTreeWidget[] = {
    "LayerTreeWidget\0\0ForceUpdate()\0"
    "OnShowAll()\0OnHideAll()\0OnLockAll()\0"
    "OnUnlockAll()\0OnShowAllInfo()\0"
    "OnHideAllInfo()\0OnSetColorMap()\0"
};

const QMetaObject LayerTreeWidget::staticMetaObject = {
    { &QTreeWidget::staticMetaObject, qt_meta_stringdata_LayerTreeWidget,
      qt_meta_data_LayerTreeWidget, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerTreeWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerTreeWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerTreeWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerTreeWidget))
        return static_cast<void*>(const_cast< LayerTreeWidget*>(this));
    return QTreeWidget::qt_metacast(_clname);
}

int LayerTreeWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QTreeWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ForceUpdate(); break;
        case 1: OnShowAll(); break;
        case 2: OnHideAll(); break;
        case 3: OnLockAll(); break;
        case 4: OnUnlockAll(); break;
        case 5: OnShowAllInfo(); break;
        case 6: OnHideAllInfo(); break;
        case 7: OnSetColorMap(); break;
        default: ;
        }
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
