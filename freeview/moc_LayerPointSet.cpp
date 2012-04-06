/****************************************************************************
** Meta object code from reading C++ file 'LayerPointSet.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerPointSet.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerPointSet.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerPointSet[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x09,
      32,   14,   14,   14, 0x09,
      48,   14,   14,   14, 0x09,
      64,   14,   14,   14, 0x09,
     101,   90,   14,   14, 0x09,
     121,   14,   14,   14, 0x29,

       0        // eod
};

static const char qt_meta_stringdata_LayerPointSet[] = {
    "LayerPointSet\0\0UpdateColorMap()\0"
    "UpdateOpacity()\0UpdateScalars()\0"
    "UpdateSnapToVoxelCenter()\0bRebuild3D\0"
    "RebuildActors(bool)\0RebuildActors()\0"
};

const QMetaObject LayerPointSet::staticMetaObject = {
    { &LayerEditable::staticMetaObject, qt_meta_stringdata_LayerPointSet,
      qt_meta_data_LayerPointSet, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerPointSet::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerPointSet::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerPointSet::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerPointSet))
        return static_cast<void*>(const_cast< LayerPointSet*>(this));
    return LayerEditable::qt_metacast(_clname);
}

int LayerPointSet::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerEditable::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: UpdateColorMap(); break;
        case 1: UpdateOpacity(); break;
        case 2: UpdateScalars(); break;
        case 3: UpdateSnapToVoxelCenter(); break;
        case 4: RebuildActors((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: RebuildActors(); break;
        default: ;
        }
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
