/****************************************************************************
** Meta object code from reading C++ file 'LayerPropertyMRI.h'
**
** Created: Fri Apr 6 15:12:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LayerPropertyMRI.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LayerPropertyMRI.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_LayerPropertyMRI[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      37,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      13,       // signalCount

 // signals: signature, parameters, type, tag, flags
      18,   17,   17,   17, 0x05,
      36,   17,   17,   17, 0x05,
      66,   17,   17,   17, 0x05,
      92,   17,   17,   17, 0x05,
     115,   17,   17,   17, 0x05,
     136,   17,   17,   17, 0x05,
     171,  164,   17,   17, 0x05,
     190,   17,   17,   17, 0x05,
     207,   17,   17,   17, 0x05,
     229,   17,   17,   17, 0x05,
     273,  264,   17,   17, 0x05,
     307,  299,   17,   17, 0x05,
     334,  164,   17,   17, 0x05,

 // slots: signature, parameters, type, tag, flags
     367,  359,   17,   17, 0x0a,
     402,  386,   17,   17, 0x0a,
     425,  264,   17,   17, 0x0a,
     463,  451,   17,   17, 0x0a,
     503,  495,   17,   17, 0x0a,
     537,  528,   17,   17, 0x0a,
     567,  560,   17,   17, 0x0a,
     592,  586,   17,   17, 0x0a,
     629,  621,   17,   17, 0x0a,
     654,  647,   17,   17, 0x0a,
     673,  671,   17,   17, 0x0a,
     696,  671,   17,   17, 0x0a,
     721,  719,   17,   17, 0x0a,
     745,  719,   17,   17, 0x0a,
     781,  774,   17,   17, 0x0a,
     812,  774,   17,   17, 0x0a,
     843,  560,   17,   17, 0x0a,
     881,  871,   17,   17, 0x0a,
     916,  908,   17,   17, 0x0a,
     947,  941,   17,   17, 0x0a,
     992,  980,   17,   17, 0x0a,
    1028, 1026,   17,   17, 0x0a,
    1058, 1052,   17,   17, 0x0a,
    1085,  941,   17,   17, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LayerPropertyMRI[] = {
    "LayerPropertyMRI\0\0ColorMapChanged()\0"
    "ResliceInterpolationChanged()\0"
    "TextureSmoothingChanged()\0"
    "OpacityChanged(double)\0WindowLevelChanged()\0"
    "GenericWindowLevelChanged()\0bShown\0"
    "ContourShown(bool)\0ContourChanged()\0"
    "ContourColorChanged()\0"
    "ContourSmoothIterationChanged(int)\0"
    "bOutline\0LabelOutlineChanged(bool)\0"
    "nMethod\0UpSampleMethodChanged(int)\0"
    "ProjectionMapShown(bool)\0opacity\0"
    "SetOpacity(double)\0nUpSampleMethod\0"
    "SetUpSampleMethod(int)\0SetShowLabelOutline(bool)\0"
    "nIterations\0SetContourSmoothIterations(int)\0"
    "iSmooth\0SetTextureSmoothing(int)\0"
    "bContour\0SetShowAsContour(bool)\0bClear\0"
    "SetClearZero(bool)\0iMode\0"
    "SetResliceInterpolation(int)\0iWindow\0"
    "SetWindow(double)\0iLevel\0SetLevel(double)\0"
    "b\0SetDisplayVector(bool)\0"
    "SetDisplayTensor(bool)\0n\0"
    "SetVectorInversion(int)\0"
    "SetVectorRepresentation(int)\0dValue\0"
    "SetContourMinThreshold(double)\0"
    "SetContourMaxThreshold(double)\0"
    "SetHeatScaleClearHigh(bool)\0bTruncate\0"
    "SetHeatScaleTruncate(bool)\0bInvert\0"
    "SetHeatScaleInvert(bool)\0bFlag\0"
    "SetContourUseImageColorMap(bool)\0"
    "bExtractAll\0SetContourExtractAllRegions(bool)\0"
    "c\0SetContourColor(QColor)\0bShow\0"
    "SetShowProjectionMap(bool)\0"
    "SetRememberFrameSettings(bool)\0"
};

const QMetaObject LayerPropertyMRI::staticMetaObject = {
    { &LayerProperty::staticMetaObject, qt_meta_stringdata_LayerPropertyMRI,
      qt_meta_data_LayerPropertyMRI, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LayerPropertyMRI::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LayerPropertyMRI::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LayerPropertyMRI::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LayerPropertyMRI))
        return static_cast<void*>(const_cast< LayerPropertyMRI*>(this));
    return LayerProperty::qt_metacast(_clname);
}

int LayerPropertyMRI::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = LayerProperty::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ColorMapChanged(); break;
        case 1: ResliceInterpolationChanged(); break;
        case 2: TextureSmoothingChanged(); break;
        case 3: OpacityChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 4: WindowLevelChanged(); break;
        case 5: GenericWindowLevelChanged(); break;
        case 6: ContourShown((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 7: ContourChanged(); break;
        case 8: ContourColorChanged(); break;
        case 9: ContourSmoothIterationChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: LabelOutlineChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: UpSampleMethodChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: ProjectionMapShown((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: SetOpacity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 14: SetUpSampleMethod((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: SetShowLabelOutline((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: SetContourSmoothIterations((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 17: SetTextureSmoothing((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: SetShowAsContour((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 19: SetClearZero((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 20: SetResliceInterpolation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 21: SetWindow((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 22: SetLevel((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 23: SetDisplayVector((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 24: SetDisplayTensor((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 25: SetVectorInversion((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 26: SetVectorRepresentation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 27: SetContourMinThreshold((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 28: SetContourMaxThreshold((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 29: SetHeatScaleClearHigh((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 30: SetHeatScaleTruncate((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 31: SetHeatScaleInvert((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 32: SetContourUseImageColorMap((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 33: SetContourExtractAllRegions((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 34: SetContourColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 35: SetShowProjectionMap((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 36: SetRememberFrameSettings((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 37;
    }
    return _id;
}

// SIGNAL 0
void LayerPropertyMRI::ColorMapChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void LayerPropertyMRI::ResliceInterpolationChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void LayerPropertyMRI::TextureSmoothingChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}

// SIGNAL 3
void LayerPropertyMRI::OpacityChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void LayerPropertyMRI::WindowLevelChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void LayerPropertyMRI::GenericWindowLevelChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 5, 0);
}

// SIGNAL 6
void LayerPropertyMRI::ContourShown(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void LayerPropertyMRI::ContourChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 7, 0);
}

// SIGNAL 8
void LayerPropertyMRI::ContourColorChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 8, 0);
}

// SIGNAL 9
void LayerPropertyMRI::ContourSmoothIterationChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 9, _a);
}

// SIGNAL 10
void LayerPropertyMRI::LabelOutlineChanged(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 10, _a);
}

// SIGNAL 11
void LayerPropertyMRI::UpSampleMethodChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 11, _a);
}

// SIGNAL 12
void LayerPropertyMRI::ProjectionMapShown(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 12, _a);
}
QT_END_MOC_NAMESPACE
