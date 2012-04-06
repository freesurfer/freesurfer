/****************************************************************************
** Meta object code from reading C++ file 'CommandEdit.h'
**
** Created: Fri Apr 6 15:12:23 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "CommandEdit.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'CommandEdit.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_CommandEdit[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      17,   13,   12,   12, 0x05,

 // slots: signature, parameters, type, tag, flags
      43,   12,   12,   12, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_CommandEdit[] = {
    "CommandEdit\0\0cmd\0CommandTriggered(QString)\0"
    "OnSelectionChanged()\0"
};

const QMetaObject CommandEdit::staticMetaObject = {
    { &QTextEdit::staticMetaObject, qt_meta_stringdata_CommandEdit,
      qt_meta_data_CommandEdit, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &CommandEdit::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *CommandEdit::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *CommandEdit::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_CommandEdit))
        return static_cast<void*>(const_cast< CommandEdit*>(this));
    return QTextEdit::qt_metacast(_clname);
}

int CommandEdit::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QTextEdit::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: CommandTriggered((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: OnSelectionChanged(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void CommandEdit::CommandTriggered(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
