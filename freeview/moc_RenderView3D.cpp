/****************************************************************************
** Meta object code from reading C++ file 'RenderView3D.h'
**
** Created: Fri Apr 6 15:12:30 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "RenderView3D.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'RenderView3D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_RenderView3D[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x05,
      37,   13,   13,   13, 0x05,
      75,   13,   13,   13, 0x05,
     123,  112,   13,   13, 0x05,

 // slots: signature, parameters, type, tag, flags
     179,  164,   13,   13, 0x0a,
     202,   13,   13,   13, 0x2a,
     227,  221,   13,   13, 0x0a,
     252,   13,   13,   13, 0x0a,
     277,   13,  272,   13, 0x0a,
     292,   13,   13,   13, 0x0a,
     312,   13,   13,   13, 0x0a,
     343,  221,   13,   13, 0x0a,
     363,   13,   13,   13, 0x2a,

       0        // eod
};

static const char qt_meta_stringdata_RenderView3D[] = {
    "RenderView3D\0\0SurfaceVertexClicked()\0"
    "SurfaceRegionSelected(SurfaceRegion*)\0"
    "SurfaceRegionRemoved(SurfaceRegion*)\0"
    "layer,info\0VolumeTrackMouseOver(Layer*,QVariantMap)\0"
    "bForScreenShot\0RefreshAllActors(bool)\0"
    "RefreshAllActors()\0bShow\0"
    "SetShowSliceFrames(bool)\0UpdateSliceFrames()\0"
    "bool\0UpdateBounds()\0SnapToNearestAxis()\0"
    "UpdateSurfaceCorrelationData()\0"
    "SetShowSlices(bool)\0SetShowSlices()\0"
};

const QMetaObject RenderView3D::staticMetaObject = {
    { &RenderView::staticMetaObject, qt_meta_stringdata_RenderView3D,
      qt_meta_data_RenderView3D, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &RenderView3D::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *RenderView3D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *RenderView3D::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_RenderView3D))
        return static_cast<void*>(const_cast< RenderView3D*>(this));
    return RenderView::qt_metacast(_clname);
}

int RenderView3D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = RenderView::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SurfaceVertexClicked(); break;
        case 1: SurfaceRegionSelected((*reinterpret_cast< SurfaceRegion*(*)>(_a[1]))); break;
        case 2: SurfaceRegionRemoved((*reinterpret_cast< SurfaceRegion*(*)>(_a[1]))); break;
        case 3: VolumeTrackMouseOver((*reinterpret_cast< Layer*(*)>(_a[1])),(*reinterpret_cast< const QVariantMap(*)>(_a[2]))); break;
        case 4: RefreshAllActors((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: RefreshAllActors(); break;
        case 6: SetShowSliceFrames((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 7: UpdateSliceFrames(); break;
        case 8: { bool _r = UpdateBounds();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 9: SnapToNearestAxis(); break;
        case 10: UpdateSurfaceCorrelationData(); break;
        case 11: SetShowSlices((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: SetShowSlices(); break;
        default: ;
        }
        _id -= 13;
    }
    return _id;
}

// SIGNAL 0
void RenderView3D::SurfaceVertexClicked()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void RenderView3D::SurfaceRegionSelected(SurfaceRegion * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void RenderView3D::SurfaceRegionRemoved(SurfaceRegion * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void RenderView3D::VolumeTrackMouseOver(Layer * _t1, const QVariantMap & _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}
QT_END_MOC_NAMESPACE
