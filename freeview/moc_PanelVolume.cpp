/****************************************************************************
** Meta object code from reading C++ file 'PanelVolume.h'
**
** Created: Fri Apr 6 15:12:29 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "PanelVolume.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'PanelVolume.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PanelVolume[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      32,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      19,   13,   12,   12, 0x09,
      49,   44,   12,   12, 0x09,
      75,   70,   12,   12, 0x09,
      96,   70,   12,   12, 0x09,
     125,  120,   12,   12, 0x09,
     179,  174,   12,   12, 0x09,
     209,   12,   12,   12, 0x09,
     229,   12,   12,   12, 0x09,
     248,   12,   12,   12, 0x09,
     265,   12,   12,   12, 0x09,
     282,   12,   12,   12, 0x09,
     299,   12,   12,   12, 0x09,
     324,  319,   12,   12, 0x09,
     350,  319,   12,   12, 0x09,
     375,  319,   12,   12, 0x09,
     398,  319,   12,   12, 0x09,
     421,  319,   12,   12, 0x09,
     444,  319,   12,   12, 0x09,
     470,   12,   12,   12, 0x09,
     494,   12,   12,   12, 0x09,
     518,   12,   12,   12, 0x09,
     545,   12,   12,   12, 0x09,
     569,   12,   12,   12, 0x09,
     585,   12,   12,   12, 0x09,
     602,   12,   12,   12, 0x09,
     620,   12,   12,   12, 0x09,
     643,   12,   12,   12, 0x09,
     671,   12,   12,   12, 0x09,
     703,   12,   12,   12, 0x09,
     722,   12,   12,   12, 0x09,
     758,  751,   12,   12, 0x09,
     786,  784,   12,   12, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_PanelVolume[] = {
    "PanelVolume\0\0bShow\0OnCheckShowContour(bool)\0"
    "nVal\0OnSliderOpacity(int)\0nSel\0"
    "OnComboColorMap(int)\0OnComboLookupTable(int)\0"
    "item\0OnColorTableCurrentItemChanged(QTreeWidgetItem*)\0"
    "strg\0OnLineEditBrushValue(QString)\0"
    "OnSliderWindow(int)\0OnSliderLevel(int)\0"
    "OnSliderMin(int)\0OnSliderMid(int)\0"
    "OnSliderMax(int)\0OnSliderOffset(int)\0"
    "text\0OnLineEditWindow(QString)\0"
    "OnLineEditLevel(QString)\0"
    "OnLineEditMin(QString)\0OnLineEditMid(QString)\0"
    "OnLineEditMax(QString)\0OnLineEditOffset(QString)\0"
    "OnSliderContourMin(int)\0OnSliderContourMax(int)\0"
    "OnSliderContourSmooth(int)\0"
    "OnContourValueChanged()\0OnContourSave()\0"
    "OnCopySettings()\0OnPasteSettings()\0"
    "OnPasteSettingsToAll()\0"
    "OnSliderTrackVolumeMin(int)\0"
    "OnTrackVolumeThresholdChanged()\0"
    "UpdateColorLabel()\0UpdateTrackVolumeThreshold()\0"
    "nFrame\0OnActiveFrameChanged(int)\0b\0"
    "OnShowExistingLabelsOnly(bool)\0"
};

const QMetaObject PanelVolume::staticMetaObject = {
    { &PanelLayer::staticMetaObject, qt_meta_stringdata_PanelVolume,
      qt_meta_data_PanelVolume, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PanelVolume::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PanelVolume::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PanelVolume::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PanelVolume))
        return static_cast<void*>(const_cast< PanelVolume*>(this));
    return PanelLayer::qt_metacast(_clname);
}

int PanelVolume::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = PanelLayer::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnCheckShowContour((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: OnSliderOpacity((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: OnComboColorMap((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: OnComboLookupTable((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: OnColorTableCurrentItemChanged((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1]))); break;
        case 5: OnLineEditBrushValue((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: OnSliderWindow((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: OnSliderLevel((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: OnSliderMin((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: OnSliderMid((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: OnSliderMax((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: OnSliderOffset((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: OnLineEditWindow((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 13: OnLineEditLevel((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 14: OnLineEditMin((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 15: OnLineEditMid((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 16: OnLineEditMax((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 17: OnLineEditOffset((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 18: OnSliderContourMin((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 19: OnSliderContourMax((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 20: OnSliderContourSmooth((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 21: OnContourValueChanged(); break;
        case 22: OnContourSave(); break;
        case 23: OnCopySettings(); break;
        case 24: OnPasteSettings(); break;
        case 25: OnPasteSettingsToAll(); break;
        case 26: OnSliderTrackVolumeMin((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 27: OnTrackVolumeThresholdChanged(); break;
        case 28: UpdateColorLabel(); break;
        case 29: UpdateTrackVolumeThreshold(); break;
        case 30: OnActiveFrameChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 31: OnShowExistingLabelsOnly((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 32;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
