/****************************************************************************
** Meta object code from reading C++ file 'WindowConfigureOverlay.h'
**
** Created: Fri Apr 6 15:12:31 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "WindowConfigureOverlay.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'WindowConfigureOverlay.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_WindowConfigureOverlay[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      24,   23,   23,   23, 0x0a,
      38,   23,   23,   23, 0x0a,
      55,   49,   23,   23, 0x09,
      90,   86,   23,   23, 0x09,
     123,  118,   23,   23, 0x09,
     149,  144,   23,   23, 0x09,
     174,   23,   23,   23, 0x09,
     195,  193,  188,   23, 0x09,
     242,   23,   23,   23, 0x09,
     280,  267,   23,   23, 0x09,
     322,   23,   23,   23, 0x09,
     349,   23,   23,   23, 0x09,
     372,  367,   23,   23, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_WindowConfigureOverlay[] = {
    "WindowConfigureOverlay\0\0UpdateGraph()\0"
    "UpdateUI()\0layer\0OnActiveSurfaceChanged(Layer*)\0"
    "btn\0OnClicked(QAbstractButton*)\0nVal\0"
    "OnSliderOpacity(int)\0dVal\0"
    "OnSpinBoxOpacity(double)\0OnButtonAdd()\0"
    "bool\0p\0UpdateOverlayProperty(SurfaceOverlayProperty*)\0"
    "UpdateThresholdChanges()\0button,value\0"
    "OnHistogramMouseButtonPressed(int,double)\0"
    "OnHistogramMarkerChanged()\0OnSmoothChanged()\0"
    "strg\0OnTextThresholdChanged(QString)\0"
};

const QMetaObject WindowConfigureOverlay::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_WindowConfigureOverlay,
      qt_meta_data_WindowConfigureOverlay, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &WindowConfigureOverlay::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *WindowConfigureOverlay::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *WindowConfigureOverlay::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_WindowConfigureOverlay))
        return static_cast<void*>(const_cast< WindowConfigureOverlay*>(this));
    if (!strcmp(_clname, "UIUpdateHelper"))
        return static_cast< UIUpdateHelper*>(const_cast< WindowConfigureOverlay*>(this));
    return QWidget::qt_metacast(_clname);
}

int WindowConfigureOverlay::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: UpdateGraph(); break;
        case 1: UpdateUI(); break;
        case 2: OnActiveSurfaceChanged((*reinterpret_cast< Layer*(*)>(_a[1]))); break;
        case 3: OnClicked((*reinterpret_cast< QAbstractButton*(*)>(_a[1]))); break;
        case 4: OnSliderOpacity((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: OnSpinBoxOpacity((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: OnButtonAdd(); break;
        case 7: { bool _r = UpdateOverlayProperty((*reinterpret_cast< SurfaceOverlayProperty*(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 8: UpdateThresholdChanges(); break;
        case 9: OnHistogramMouseButtonPressed((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 10: OnHistogramMarkerChanged(); break;
        case 11: OnSmoothChanged(); break;
        case 12: OnTextThresholdChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 13;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
