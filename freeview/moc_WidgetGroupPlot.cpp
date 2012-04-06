/****************************************************************************
** Meta object code from reading C++ file 'WidgetGroupPlot.h'
**
** Created: Fri Apr 6 15:12:33 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "WidgetGroupPlot.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'WidgetGroupPlot.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_WidgetGroupPlot[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      23,   17,   16,   16, 0x05,

 // slots: signature, parameters, type, tag, flags
      47,   41,   16,   16, 0x0a,
      68,   66,   16,   16, 0x0a,
      97,   66,   16,   16, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_WidgetGroupPlot[] = {
    "WidgetGroupPlot\0\0frame\0FrameChanged(int)\0"
    "bAuto\0SetAutoScale(bool)\0n\0"
    "SetCurrentVariableIndex(int)\0"
    "SetPlotType(int)\0"
};

const QMetaObject WidgetGroupPlot::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_WidgetGroupPlot,
      qt_meta_data_WidgetGroupPlot, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &WidgetGroupPlot::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *WidgetGroupPlot::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *WidgetGroupPlot::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_WidgetGroupPlot))
        return static_cast<void*>(const_cast< WidgetGroupPlot*>(this));
    return QWidget::qt_metacast(_clname);
}

int WidgetGroupPlot::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: FrameChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: SetAutoScale((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: SetCurrentVariableIndex((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: SetPlotType((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void WidgetGroupPlot::FrameChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
