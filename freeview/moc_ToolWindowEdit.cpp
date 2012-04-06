/****************************************************************************
** Meta object code from reading C++ file 'ToolWindowEdit.h'
**
** Created: Fri Apr 6 15:12:30 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ToolWindowEdit.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ToolWindowEdit.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ToolWindowEdit[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      16,   15,   15,   15, 0x0a,
      32,   15,   15,   15, 0x09,
      45,   41,   15,   15, 0x09,
      70,   66,   15,   15, 0x09,
      97,   92,   15,   15, 0x09,
     129,   92,   15,   15, 0x09,
     157,   92,   15,   15, 0x09,
     185,   92,   15,   15, 0x09,
     216,   15,   15,   15, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_ToolWindowEdit[] = {
    "ToolWindowEdit\0\0UpdateWidgets()\0"
    "OnIdle()\0act\0OnEditMode(QAction*)\0sel\0"
    "OnComboReference(int)\0strg\0"
    "OnLineEditContourValue(QString)\0"
    "OnLineEditSmoothSD(QString)\0"
    "OnDrawRangeChanged(QString)\0"
    "OnExcludeRangeChanged(QString)\0"
    "OnReplaceLabel()\0"
};

const QMetaObject ToolWindowEdit::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ToolWindowEdit,
      qt_meta_data_ToolWindowEdit, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ToolWindowEdit::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ToolWindowEdit::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ToolWindowEdit::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ToolWindowEdit))
        return static_cast<void*>(const_cast< ToolWindowEdit*>(this));
    if (!strcmp(_clname, "UIUpdateHelper"))
        return static_cast< UIUpdateHelper*>(const_cast< ToolWindowEdit*>(this));
    return QWidget::qt_metacast(_clname);
}

int ToolWindowEdit::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: UpdateWidgets(); break;
        case 1: OnIdle(); break;
        case 2: OnEditMode((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 3: OnComboReference((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: OnLineEditContourValue((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: OnLineEditSmoothSD((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: OnDrawRangeChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: OnExcludeRangeChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 8: OnReplaceLabel(); break;
        default: ;
        }
        _id -= 9;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
