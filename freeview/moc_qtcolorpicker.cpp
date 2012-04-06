/****************************************************************************
** Meta object code from reading C++ file 'qtcolorpicker.h'
**
** Created: Fri Apr 6 15:12:29 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qtcolorpicker.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qtcolorpicker.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_QtColorPicker[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       1,   34, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x05,

 // slots: signature, parameters, type, tag, flags
      40,   36,   14,   14, 0x0a,
      72,   64,   14,   14, 0x08,
      92,   14,   14,   14, 0x08,

 // properties: name, type, flags
     111,  106, 0x01095003,

       0        // eod
};

static const char qt_meta_stringdata_QtColorPicker[] = {
    "QtColorPicker\0\0colorChanged(QColor)\0"
    "col\0setCurrentColor(QColor)\0toggled\0"
    "buttonPressed(bool)\0popupClosed()\0"
    "bool\0colorDialog\0"
};

const QMetaObject QtColorPicker::staticMetaObject = {
    { &QPushButton::staticMetaObject, qt_meta_stringdata_QtColorPicker,
      qt_meta_data_QtColorPicker, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &QtColorPicker::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *QtColorPicker::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *QtColorPicker::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_QtColorPicker))
        return static_cast<void*>(const_cast< QtColorPicker*>(this));
    return QPushButton::qt_metacast(_clname);
}

int QtColorPicker::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QPushButton::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: colorChanged((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 1: setCurrentColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 2: buttonPressed((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: popupClosed(); break;
        default: ;
        }
        _id -= 4;
    }
#ifndef QT_NO_PROPERTIES
      else if (_c == QMetaObject::ReadProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< bool*>(_v) = colorDialogEnabled(); break;
        }
        _id -= 1;
    } else if (_c == QMetaObject::WriteProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: setColorDialogEnabled(*reinterpret_cast< bool*>(_v)); break;
        }
        _id -= 1;
    } else if (_c == QMetaObject::ResetProperty) {
        _id -= 1;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 1;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 1;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 1;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 1;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 1;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}

// SIGNAL 0
void QtColorPicker::colorChanged(const QColor & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
