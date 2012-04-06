/****************************************************************************
** Meta object code from reading C++ file 'DialogRepositionSurface.h'
**
** Created: Fri Apr 6 15:12:32 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "DialogRepositionSurface.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'DialogRepositionSurface.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_DialogRepositionSurface[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      25,   24,   24,   24, 0x0a,
      35,   24,   24,   24, 0x0a,
      44,   24,   24,   24, 0x0a,
      53,   24,   24,   24, 0x0a,
      66,   64,   24,   24, 0x0a,
      85,   24,   24,   24, 0x0a,
      96,   24,   24,   24, 0x0a,
     121,   24,   24,   24, 0x0a,
     136,   24,   24,   24, 0x0a,
     154,   24,   24,   24, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_DialogRepositionSurface[] = {
    "DialogRepositionSurface\0\0OnApply()\0"
    "OnUndo()\0OnSave()\0OnSaveAs()\0n\0"
    "OnComboTarget(int)\0UpdateUI()\0"
    "OnSurfaceVertexClicked()\0UpdateVertex()\0"
    "UpdateIntensity()\0OnCoordinateTypeChanged()\0"
};

const QMetaObject DialogRepositionSurface::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_DialogRepositionSurface,
      qt_meta_data_DialogRepositionSurface, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DialogRepositionSurface::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DialogRepositionSurface::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DialogRepositionSurface::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DialogRepositionSurface))
        return static_cast<void*>(const_cast< DialogRepositionSurface*>(this));
    return QDialog::qt_metacast(_clname);
}

int DialogRepositionSurface::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: OnApply(); break;
        case 1: OnUndo(); break;
        case 2: OnSave(); break;
        case 3: OnSaveAs(); break;
        case 4: OnComboTarget((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: UpdateUI(); break;
        case 6: OnSurfaceVertexClicked(); break;
        case 7: UpdateVertex(); break;
        case 8: UpdateIntensity(); break;
        case 9: OnCoordinateTypeChanged(); break;
        default: ;
        }
        _id -= 10;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
