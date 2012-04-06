/****************************************************************************
** Meta object code from reading C++ file 'PanelLayer.h'
**
** Created: Fri Apr 6 15:12:28 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "PanelLayer.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'PanelLayer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PanelLayer[] = {

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
      12,   11,   11,   11, 0x0a,
      28,   11,   11,   11, 0x09,
      37,   11,   11,   11, 0x09,
      54,   48,   11,   11, 0x09,
      75,   48,   11,   11, 0x09,
      98,   48,   11,   11, 0x09,
     119,   48,   11,   11, 0x09,
     153,  148,   11,   11, 0x09,
     185,  148,   11,   11, 0x09,
     224,   11,   11,   11, 0x09,
     255,  249,   11,   11, 0x09,
     288,   11,   11,   11, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_PanelLayer[] = {
    "PanelLayer\0\0UpdateWidgets()\0OnIdle()\0"
    "OnUpdate()\0layer\0OnLayerAdded(Layer*)\0"
    "OnLayerRemoved(Layer*)\0OnLayerMoved(Layer*)\0"
    "OnActiveLayerChanged(Layer*)\0item\0"
    "OnItemChanged(QTreeWidgetItem*)\0"
    "OnCurrentItemChanged(QTreeWidgetItem*)\0"
    "OnItemSelectionChanged()\0index\0"
    "OnItemDoubleClicked(QModelIndex)\0"
    "OnLayerNameChanged()\0"
};

const QMetaObject PanelLayer::staticMetaObject = {
    { &QScrollArea::staticMetaObject, qt_meta_stringdata_PanelLayer,
      qt_meta_data_PanelLayer, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PanelLayer::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PanelLayer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PanelLayer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PanelLayer))
        return static_cast<void*>(const_cast< PanelLayer*>(this));
    if (!strcmp(_clname, "UIUpdateHelper"))
        return static_cast< UIUpdateHelper*>(const_cast< PanelLayer*>(this));
    return QScrollArea::qt_metacast(_clname);
}

int PanelLayer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QScrollArea::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: UpdateWidgets(); break;
        case 1: OnIdle(); break;
        case 2: OnUpdate(); break;
        case 3: OnLayerAdded((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 4: OnLayerRemoved((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 5: OnLayerMoved((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 6: OnActiveLayerChanged((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 7: OnItemChanged((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1]))); break;
        case 8: OnCurrentItemChanged((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1]))); break;
        case 9: OnItemSelectionChanged(); break;
        case 10: OnItemDoubleClicked((*reinterpret_cast< const QModelIndex(*)>(_a[1]))); break;
        case 11: OnLayerNameChanged(); break;
        default: ;
        }
        _id -= 12;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
