/****************************************************************************
** Meta object code from reading C++ file 'WidgetHistogram.h'
**
** Created: Fri Apr 6 15:12:31 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "WidgetHistogram.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'WidgetHistogram.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_WidgetHistogram[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      30,   17,   16,   16, 0x05,
      61,   16,   16,   16, 0x05,

 // slots: signature, parameters, type, tag, flags
      84,   77,   16,   16, 0x0a,
     109,  103,   16,   16, 0x0a,
     136,   16,   16,   16, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_WidgetHistogram[] = {
    "WidgetHistogram\0\0button,value\0"
    "MouseButtonPressed(int,double)\0"
    "MarkerChanged()\0bRange\0SetAutoRange(bool)\0"
    "color\0SetForegroundColor(QColor)\0"
    "FlipMarkers()\0"
};

const QMetaObject WidgetHistogram::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_WidgetHistogram,
      qt_meta_data_WidgetHistogram, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &WidgetHistogram::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *WidgetHistogram::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *WidgetHistogram::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_WidgetHistogram))
        return static_cast<void*>(const_cast< WidgetHistogram*>(this));
    return QWidget::qt_metacast(_clname);
}

int WidgetHistogram::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: MouseButtonPressed((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 1: MarkerChanged(); break;
        case 2: SetAutoRange((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: SetForegroundColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 4: FlipMarkers(); break;
        default: ;
        }
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void WidgetHistogram::MouseButtonPressed(int _t1, double _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void WidgetHistogram::MarkerChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
