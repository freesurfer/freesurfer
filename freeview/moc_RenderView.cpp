/****************************************************************************
** Meta object code from reading C++ file 'RenderView.h'
**
** Created: Fri Apr 6 15:12:29 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "RenderView.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'RenderView.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_RenderView[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      14,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      19,   12,   11,   11, 0x0a,
      39,   11,   11,   11, 0x2a,
      55,   11,   11,   11, 0x0a,
      64,   11,   11,   11, 0x0a,
      75,   11,   11,   11, 0x0a,
      86,   11,   11,   11, 0x0a,
     105,   98,   11,   11, 0x0a,
     118,   11,   11,   11, 0x0a,
     134,  126,   11,   11, 0x0a,
     155,  149,   11,   11, 0x0a,
     181,  175,   11,   11, 0x0a,
     211,  207,   11,   11, 0x0a,
     239,   11,   11,   11, 0x09,
     248,   11,   11,   11, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_RenderView[] = {
    "RenderView\0\0bForce\0RequestRedraw(bool)\0"
    "RequestRedraw()\0MoveUp()\0MoveDown()\0"
    "MoveLeft()\0MoveRight()\0factor\0"
    "Zoom(double)\0Reset()\0nAction\0"
    "SetAction(int)\0bShow\0ShowScalarBar(bool)\0"
    "layer\0SetScalarBarLayer(Layer*)\0act\0"
    "SetScalarBarLayer(QAction*)\0OnIdle()\0"
    "OnSlicePositionChanged()\0"
};

const QMetaObject RenderView::staticMetaObject = {
    { &GenericRenderView::staticMetaObject, qt_meta_stringdata_RenderView,
      qt_meta_data_RenderView, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &RenderView::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *RenderView::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *RenderView::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_RenderView))
        return static_cast<void*>(const_cast< RenderView*>(this));
    return GenericRenderView::qt_metacast(_clname);
}

int RenderView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GenericRenderView::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: RequestRedraw((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: RequestRedraw(); break;
        case 2: MoveUp(); break;
        case 3: MoveDown(); break;
        case 4: MoveLeft(); break;
        case 5: MoveRight(); break;
        case 6: Zoom((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 7: Reset(); break;
        case 8: SetAction((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: ShowScalarBar((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: SetScalarBarLayer((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 11: SetScalarBarLayer((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 12: OnIdle(); break;
        case 13: OnSlicePositionChanged(); break;
        default: ;
        }
        _id -= 14;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
