/****************************************************************************
** Meta object code from reading C++ file 'Contour2D.h'
**
** Created: Fri Apr 6 15:12:23 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "Contour2D.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'Contour2D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Contour2D[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      11,   10,   10,   10, 0x05,
      26,   10,   10,   10, 0x05,

 // slots: signature, parameters, type, tag, flags
      48,   41,   10,   10, 0x0a,
      80,   72,   10,   10, 0x0a,
      99,   96,   10,   10, 0x0a,
     125,  119,   10,   10, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Contour2D[] = {
    "Contour2D\0\0ValueChanged()\0ColorChanged()\0"
    "dvalue\0SetContourValue(double)\0bSmooth\0"
    "SetSmooth(bool)\0sd\0SetSmoothSD(double)\0"
    "color\0SetColor(QColor)\0"
};

const QMetaObject Contour2D::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Contour2D,
      qt_meta_data_Contour2D, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Contour2D::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Contour2D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Contour2D::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Contour2D))
        return static_cast<void*>(const_cast< Contour2D*>(this));
    return QObject::qt_metacast(_clname);
}

int Contour2D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ValueChanged(); break;
        case 1: ColorChanged(); break;
        case 2: SetContourValue((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 3: SetSmooth((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 4: SetSmoothSD((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 5: SetColor((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 6;
    }
    return _id;
}

// SIGNAL 0
void Contour2D::ValueChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void Contour2D::ColorChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
