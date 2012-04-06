/****************************************************************************
** Meta object code from reading C++ file 'InfoTreeWidget.h'
**
** Created: Fri Apr 6 15:12:26 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "InfoTreeWidget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'InfoTreeWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_InfoTreeWidget[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      22,   16,   15,   15, 0x05,

 // slots: signature, parameters, type, tag, flags
      74,   63,   15,   15, 0x0a,
     122,   15,   15,   15, 0x0a,
     134,   15,   15,   15, 0x09,
     159,   15,   15,   15, 0x09,
     197,  185,   15,   15, 0x09,
     233,   15,   15,   15, 0x09,
     256,  250,   15,   15, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_InfoTreeWidget[] = {
    "InfoTreeWidget\0\0x,y,z\0"
    "RASChangeTriggered(double,double,double)\0"
    "layer,info\0UpdateTrackVolumeAnnotation(Layer*,QVariantMap)\0"
    "UpdateAll()\0OnMousePositionChanged()\0"
    "OnCursorPositionChanged()\0item,column\0"
    "OnItemClicked(QTreeWidgetItem*,int)\0"
    "OnEditFinished()\0bShow\0OnToggleShowInfo(bool)\0"
};

const QMetaObject InfoTreeWidget::staticMetaObject = {
    { &QTreeWidget::staticMetaObject, qt_meta_stringdata_InfoTreeWidget,
      qt_meta_data_InfoTreeWidget, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &InfoTreeWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *InfoTreeWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *InfoTreeWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_InfoTreeWidget))
        return static_cast<void*>(const_cast< InfoTreeWidget*>(this));
    return QTreeWidget::qt_metacast(_clname);
}

int InfoTreeWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QTreeWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: RASChangeTriggered((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 1: UpdateTrackVolumeAnnotation((*reinterpret_cast< Layer*(*)>(_a[1])),(*reinterpret_cast< const QVariantMap(*)>(_a[2]))); break;
        case 2: UpdateAll(); break;
        case 3: OnMousePositionChanged(); break;
        case 4: OnCursorPositionChanged(); break;
        case 5: OnItemClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 6: OnEditFinished(); break;
        case 7: OnToggleShowInfo((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void InfoTreeWidget::RASChangeTriggered(double _t1, double _t2, double _t3)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
