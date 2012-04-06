/****************************************************************************
** Meta object code from reading C++ file 'RenderView2D.h'
**
** Created: Fri Apr 6 15:12:29 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "RenderView2D.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'RenderView2D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_RenderView2D[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x05,
      40,   13,   13,   13, 0x05,
      70,   65,   13,   13, 0x05,

 // slots: signature, parameters, type, tag, flags
     108,   93,   13,   13, 0x0a,
     131,   13,   13,   13, 0x2a,
     150,   13,   13,   13, 0x0a,
     166,   13,   13,   13, 0x0a,
     185,   13,   13,   13, 0x0a,
     209,  203,   13,   13, 0x0a,
     240,   13,   13,   13, 0x09,
     265,   65,   13,   13, 0x09,
     291,   13,   13,   13, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_RenderView2D[] = {
    "RenderView2D\0\0RegionSelected(Region2D*)\0"
    "RegionRemoved(Region2D*)\0view\0"
    "Zooming(RenderView2D*)\0bForScreenShot\0"
    "RefreshAllActors(bool)\0RefreshAllActors()\0"
    "StopSelection()\0UpdateAnnotation()\0"
    "Update2DOverlay()\0bShow\0"
    "ShowCoordinateAnnotation(bool)\0"
    "OnSlicePositionChanged()\0"
    "SyncZoomTo(RenderView2D*)\0OnDuplicateRegion()\0"
};

const QMetaObject RenderView2D::staticMetaObject = {
    { &RenderView::staticMetaObject, qt_meta_stringdata_RenderView2D,
      qt_meta_data_RenderView2D, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &RenderView2D::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *RenderView2D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *RenderView2D::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_RenderView2D))
        return static_cast<void*>(const_cast< RenderView2D*>(this));
    return RenderView::qt_metacast(_clname);
}

int RenderView2D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = RenderView::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: RegionSelected((*reinterpret_cast< Region2D*(*)>(_a[1]))); break;
        case 1: RegionRemoved((*reinterpret_cast< Region2D*(*)>(_a[1]))); break;
        case 2: Zooming((*reinterpret_cast< RenderView2D*(*)>(_a[1]))); break;
        case 3: RefreshAllActors((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 4: RefreshAllActors(); break;
        case 5: StopSelection(); break;
        case 6: UpdateAnnotation(); break;
        case 7: Update2DOverlay(); break;
        case 8: ShowCoordinateAnnotation((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: OnSlicePositionChanged(); break;
        case 10: SyncZoomTo((*reinterpret_cast< RenderView2D*(*)>(_a[1]))); break;
        case 11: OnDuplicateRegion(); break;
        default: ;
        }
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void RenderView2D::RegionSelected(Region2D * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void RenderView2D::RegionRemoved(Region2D * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void RenderView2D::Zooming(RenderView2D * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
QT_END_MOC_NAMESPACE
