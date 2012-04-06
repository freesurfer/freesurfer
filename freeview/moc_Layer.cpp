/****************************************************************************
** Meta object code from reading C++ file 'Layer.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "Layer.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'Layer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Layer[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       6,       // signalCount

 // signals: signature, parameters, type, tag, flags
      12,    7,    6,    6, 0x05,
      33,    6,    6,    6, 0x05,
      53,   47,    6,    6, 0x05,
      66,    6,    6,    6, 0x05,
      81,    6,    6,    6, 0x05,
     105,   96,    6,    6, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_Layer[] = {
    "Layer\0\0name\0NameChanged(QString)\0"
    "Transformed()\0bLock\0Locked(bool)\0"
    "ActorUpdated()\0ActorChanged()\0bVisible\0"
    "VisibilityChanged(bool)\0"
};

const QMetaObject Layer::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Layer,
      qt_meta_data_Layer, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Layer::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Layer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Layer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Layer))
        return static_cast<void*>(const_cast< Layer*>(this));
    return QObject::qt_metacast(_clname);
}

int Layer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: NameChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: Transformed(); break;
        case 2: Locked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: ActorUpdated(); break;
        case 4: ActorChanged(); break;
        case 5: VisibilityChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 6;
    }
    return _id;
}

// SIGNAL 0
void Layer::NameChanged(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void Layer::Transformed()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void Layer::Locked(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void Layer::ActorUpdated()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}

// SIGNAL 4
void Layer::ActorChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void Layer::VisibilityChanged(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}
QT_END_MOC_NAMESPACE
