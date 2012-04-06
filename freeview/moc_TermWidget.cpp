/****************************************************************************
** Meta object code from reading C++ file 'TermWidget.h'
**
** Created: Fri Apr 6 15:12:30 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "TermWidget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TermWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_TermWidget[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      18,   12,   11,   11, 0x0a,
      37,   11,   11,   11, 0x2a,
      59,   52,   11,   11, 0x0a,
      86,   11,   11,   11, 0x2a,
     114,  109,   11,   11, 0x0a,
     142,   11,   11,   11, 0x2a,
     170,  166,   11,   11, 0x09,
     198,   11,   11,   11, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_TermWidget[] = {
    "TermWidget\0\0bDark\0SetDarkTheme(bool)\0"
    "SetDarkTheme()\0bRedir\0SetRedirectStdOutput(bool)\0"
    "SetRedirectStdOutput()\0bDup\0"
    "SetDuplicateStdOutput(bool)\0"
    "SetDuplicateStdOutput()\0cmd\0"
    "OnCommandTriggered(QString)\0OnTimeOut()\0"
};

const QMetaObject TermWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_TermWidget,
      qt_meta_data_TermWidget, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &TermWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *TermWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *TermWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_TermWidget))
        return static_cast<void*>(const_cast< TermWidget*>(this));
    return QWidget::qt_metacast(_clname);
}

int TermWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SetDarkTheme((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: SetDarkTheme(); break;
        case 2: SetRedirectStdOutput((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: SetRedirectStdOutput(); break;
        case 4: SetDuplicateStdOutput((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: SetDuplicateStdOutput(); break;
        case 6: OnCommandTriggered((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: OnTimeOut(); break;
        default: ;
        }
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
