/****************************************************************************
** Meta object code from reading C++ file 'GenericRenderView.h'
**
** Created: Fri Apr 6 15:12:26 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "GenericRenderView.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GenericRenderView.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GenericRenderView[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      42,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       8,       // signalCount

 // signals: signature, parameters, type, tag, flags
      19,   18,   18,   18, 0x05,
      44,   18,   18,   18, 0x05,
      83,   18,   18,   18, 0x05,
     114,   18,   18,   18, 0x05,
     130,   18,   18,   18, 0x05,
     163,   18,   18,   18, 0x05,
     197,   18,   18,   18, 0x05,
     231,   18,   18,   18, 0x05,

 // slots: signature, parameters, type, tag, flags
     265,   18,   18,   18, 0x0a,
     295,  292,   18,   18, 0x0a,
     331,  322,   18,   18, 0x0a,
     359,  357,   18,   18, 0x2a,
     389,  380,   18,   18, 0x0a,
     418,  416,   18,   18, 0x2a,
     440,   18,   18,   18, 0x0a,
     466,  458,   18,   18, 0x0a,
     493,  489,   18,   18, 0x0a,
     515,   18,   18,   18, 0x0a,
     533,   18,   18,   18, 0x0a,
     559,   18,   18,   18, 0x0a,
     584,   18,   18,   18, 0x0a,
     612,   18,   18,   18, 0x0a,
     637,   18,   18,   18, 0x0a,
     666,  416,   18,   18, 0x0a,
     697,   18,   18,   18, 0x2a,
     731,  724,   18,   18, 0x0a,
     764,  755,   18,   18, 0x0a,
     800,  798,   18,   18, 0x2a,
     829,  755,   18,   18, 0x0a,
     864,  798,   18,   18, 0x2a,
     894,  755,   18,   18, 0x0a,
     929,  798,   18,   18, 0x2a,
     959,  755,   18,   18, 0x0a,
     994,  798,   18,   18, 0x2a,
    1043, 1024,   18,   18, 0x0a,
    1090,   18,   18,   18, 0x0a,
    1110,   18,   18,   18, 0x0a,
    1144,   18,   18,   18, 0x0a,
    1169, 1154,   18,   18, 0x0a,
    1192,   18,   18,   18, 0x2a,
    1211,   18,   18,   18, 0x08,
    1229,   18,   18,   18, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_GenericRenderView[] = {
    "GenericRenderView\0\0RenderTriggeredByWheel()\0"
    "MouseReleasedWithoutMove(QMouseEvent*)\0"
    "BackgroundColorChanged(QColor)\0"
    "ActorsUpdated()\0KeyLightIntensityChanged(double)\0"
    "FillLightIntensityChanged(double)\0"
    "HeadLightIntensityChanged(double)\0"
    "BackLightIntensityChanged(double)\0"
    "ResetCameraClippingRange()\0qc\0"
    "SetBackgroundColor(QColor)\0n,redraw\0"
    "SetAntialiasing(int,bool)\0n\0"
    "SetAntialiasing(int)\0b,redraw\0"
    "SetAntialiasing(bool,bool)\0b\0"
    "SetAntialiasing(bool)\0CopyToClipboard()\0"
    "bEnable\0EnableInteractor(bool)\0bOn\0"
    "SetStereoRender(bool)\0StereoRenderOff()\0"
    "SetStereoTypeToAnaglyph()\0"
    "SetStereoTypeToRedBlue()\0"
    "SetStereoTypeToInterlaced()\0"
    "SetStereoTypeToDresden()\0"
    "SetStereoTypeToCrystalEyes()\0"
    "SetStereoTypeToLeftRight(bool)\0"
    "SetStereoTypeToLeftRight()\0nAngle\0"
    "SetStereoPairAngle(int)\0d,redraw\0"
    "SetKeyLightIntensity(double,bool)\0d\0"
    "SetKeyLightIntensity(double)\0"
    "SetFillLightIntensity(double,bool)\0"
    "SetFillLightIntensity(double)\0"
    "SetHeadLightIntensity(double,bool)\0"
    "SetHeadLightIntensity(double)\0"
    "SetBackLightIntensity(double,bool)\0"
    "SetBackLightIntensity(double)\0"
    "key,head,fill,back\0"
    "SetLightIntensity(double,double,double,double)\0"
    "SetLightToDefault()\0"
    "EmitRenderTriggeredByInteractor()\0"
    "Refresh()\0bForScreenShot\0"
    "RefreshAllActors(bool)\0RefreshAllActors()\0"
    "UpdateRenderer2()\0UpdateCamera2()\0"
};

const QMetaObject GenericRenderView::staticMetaObject = {
    { &QVTKWidget::staticMetaObject, qt_meta_stringdata_GenericRenderView,
      qt_meta_data_GenericRenderView, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &GenericRenderView::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *GenericRenderView::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *GenericRenderView::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GenericRenderView))
        return static_cast<void*>(const_cast< GenericRenderView*>(this));
    return QVTKWidget::qt_metacast(_clname);
}

int GenericRenderView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QVTKWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: RenderTriggeredByWheel(); break;
        case 1: MouseReleasedWithoutMove((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 2: BackgroundColorChanged((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 3: ActorsUpdated(); break;
        case 4: KeyLightIntensityChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 5: FillLightIntensityChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: HeadLightIntensityChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 7: BackLightIntensityChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 8: ResetCameraClippingRange(); break;
        case 9: SetBackgroundColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 10: SetAntialiasing((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 11: SetAntialiasing((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: SetAntialiasing((*reinterpret_cast< bool(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 13: SetAntialiasing((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: CopyToClipboard(); break;
        case 15: EnableInteractor((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: SetStereoRender((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: StereoRenderOff(); break;
        case 18: SetStereoTypeToAnaglyph(); break;
        case 19: SetStereoTypeToRedBlue(); break;
        case 20: SetStereoTypeToInterlaced(); break;
        case 21: SetStereoTypeToDresden(); break;
        case 22: SetStereoTypeToCrystalEyes(); break;
        case 23: SetStereoTypeToLeftRight((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 24: SetStereoTypeToLeftRight(); break;
        case 25: SetStereoPairAngle((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 26: SetKeyLightIntensity((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 27: SetKeyLightIntensity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 28: SetFillLightIntensity((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 29: SetFillLightIntensity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 30: SetHeadLightIntensity((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 31: SetHeadLightIntensity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 32: SetBackLightIntensity((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 33: SetBackLightIntensity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 34: SetLightIntensity((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< double(*)>(_a[4]))); break;
        case 35: SetLightToDefault(); break;
        case 36: EmitRenderTriggeredByInteractor(); break;
        case 37: Refresh(); break;
        case 38: RefreshAllActors((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 39: RefreshAllActors(); break;
        case 40: UpdateRenderer2(); break;
        case 41: UpdateCamera2(); break;
        default: ;
        }
        _id -= 42;
    }
    return _id;
}

// SIGNAL 0
void GenericRenderView::RenderTriggeredByWheel()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void GenericRenderView::MouseReleasedWithoutMove(QMouseEvent * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void GenericRenderView::BackgroundColorChanged(const QColor & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void GenericRenderView::ActorsUpdated()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}

// SIGNAL 4
void GenericRenderView::KeyLightIntensityChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void GenericRenderView::FillLightIntensityChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void GenericRenderView::HeadLightIntensityChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void GenericRenderView::BackLightIntensityChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}
QT_END_MOC_NAMESPACE
