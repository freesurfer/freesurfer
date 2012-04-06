/****************************************************************************
** Meta object code from reading C++ file 'LayerCollection.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerCollection.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerCollection.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerCollection[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      18,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      14,       // signalCount

 // signals: signature, parameters, type, tag, flags
      17,   16,   16,   16, 0x05,
      44,   16,   16,   16, 0x05,
      63,   16,   16,   16, 0x05,
      84,   16,   16,   16, 0x05,
     104,   16,   16,   16, 0x05,
     123,   16,   16,   16, 0x05,
     143,   16,   16,   16, 0x05,
     163,   16,   16,   16, 0x05,
     186,   16,   16,   16, 0x05,
     211,   16,   16,   16, 0x05,
     234,   16,   16,   16, 0x05,
     250,   16,   16,   16, 0x05,
     269,   16,   16,   16, 0x05,
     295,   16,   16,   16, 0x05,

 // slots: signature, parameters, type, tag, flags
     328,  322,   16,   16, 0x0a,
     346,   16,   16,   16, 0x0a,
     360,   16,   16,   16, 0x0a,
     382,  376,   16,   16, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerCollection[] = {
    "LayerCollection\0\0ActiveLayerChanged(Layer*)\0"
    "LayerAdded(Layer*)\0LayerRemoved(Layer*)\0"
    "LayerCycled(Layer*)\0LayerMoved(Layer*)\0"
    "LayerActorUpdated()\0LayerActorChanged()\0"
    "LayerPropertyChanged()\0LayerVisibilityChanged()\0"
    "LayerShowInfoChanged()\0LayerModified()\0"
    "LayerNameChanged()\0MouseRASPositionChanged()\0"
    "CursorRASPositionChanged()\0bLock\0"
    "LockCurrent(bool)\0MoveLayerUp()\0"
    "MoveLayerDown()\0x,y,z\0"
    "SetMouseRASPosition(double,double,double)\0"
};

const QMetaObject LayerCollection::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_LayerCollection,
      qt_meta_data_LayerCollection, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerCollection::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerCollection::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerCollection::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerCollection))
        return static_cast<void*>(const_cast< LayerCollection*>(this));
    return QObject::qt_metacast(_clname);
}

int LayerCollection::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ActiveLayerChanged((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 1: LayerAdded((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 2: LayerRemoved((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 3: LayerCycled((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 4: LayerMoved((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 5: LayerActorUpdated(); break;
        case 6: LayerActorChanged(); break;
        case 7: LayerPropertyChanged(); break;
        case 8: LayerVisibilityChanged(); break;
        case 9: LayerShowInfoChanged(); break;
        case 10: LayerModified(); break;
        case 11: LayerNameChanged(); break;
        case 12: MouseRASPositionChanged(); break;
        case 13: CursorRASPositionChanged(); break;
        case 14: LockCurrent((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 15: MoveLayerUp(); break;
        case 16: MoveLayerDown(); break;
        case 17: SetMouseRASPosition((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        default: ;
        }
        _id -= 18;
    }
    return _id;
}

// SIGNAL 0
void LayerCollection::ActiveLayerChanged(Layer * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void LayerCollection::LayerAdded(Layer * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void LayerCollection::LayerRemoved(Layer * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void LayerCollection::LayerCycled(Layer * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void LayerCollection::LayerMoved(Layer * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void LayerCollection::LayerActorUpdated()
{
    QMetaObject::activate(this, &staticMetaObject, 5, 0);
}

// SIGNAL 6
void LayerCollection::LayerActorChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 6, 0);
}

// SIGNAL 7
void LayerCollection::LayerPropertyChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 7, 0);
}

// SIGNAL 8
void LayerCollection::LayerVisibilityChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 8, 0);
}

// SIGNAL 9
void LayerCollection::LayerShowInfoChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 9, 0);
}

// SIGNAL 10
void LayerCollection::LayerModified()
{
    QMetaObject::activate(this, &staticMetaObject, 10, 0);
}

// SIGNAL 11
void LayerCollection::LayerNameChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 11, 0);
}

// SIGNAL 12
void LayerCollection::MouseRASPositionChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 12, 0);
}

// SIGNAL 13
void LayerCollection::CursorRASPositionChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 13, 0);
}
QT_END_MOC_NAMESPACE
