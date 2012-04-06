/****************************************************************************
** Meta object code from reading C++ file 'FloatingStatusBar.h'
**
** Created: Fri Apr 6 15:12:25 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "FloatingStatusBar.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'FloatingStatusBar.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_FloatingStatusBar[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      29,   19,   18,   18, 0x0a,
      46,   18,   18,   18, 0x0a,
      61,   18,   18,   18, 0x0a,
      76,   18,   18,   18, 0x0a,
      89,   18,   18,   18, 0x0a,
     102,   18,   18,   18, 0x0a,
     114,   18,   18,   18, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_FloatingStatusBar[] = {
    "FloatingStatusBar\0\0nProgress\0"
    "SetProgress(int)\0ShowProgress()\0"
    "HideProgress()\0Reposition()\0StartTimer()\0"
    "StopTimer()\0OnProgressTimer()\0"
};

const QMetaObject FloatingStatusBar::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_FloatingStatusBar,
      qt_meta_data_FloatingStatusBar, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &FloatingStatusBar::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *FloatingStatusBar::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *FloatingStatusBar::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_FloatingStatusBar))
        return static_cast<void*>(const_cast< FloatingStatusBar*>(this));
    return QWidget::qt_metacast(_clname);
}

int FloatingStatusBar::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SetProgress((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: ShowProgress(); break;
        case 2: HideProgress(); break;
        case 3: Reposition(); break;
        case 4: StartTimer(); break;
        case 5: StopTimer(); break;
        case 6: OnProgressTimer(); break;
        default: ;
        }
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
