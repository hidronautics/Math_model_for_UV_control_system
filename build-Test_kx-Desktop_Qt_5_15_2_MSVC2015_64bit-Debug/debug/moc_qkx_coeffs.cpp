/****************************************************************************
** Meta object code from reading C++ file 'qkx_coeffs.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../Test_kx/qkx_coeffs.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qkx_coeffs.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Qkx_coeffs_t {
    QByteArrayData data[10];
    char stringdata0[89];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Qkx_coeffs_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Qkx_coeffs_t qt_meta_stringdata_Qkx_coeffs = {
    {
QT_MOC_LITERAL(0, 0, 10), // "Qkx_coeffs"
QT_MOC_LITERAL(1, 11, 11), // "sendSucceed"
QT_MOC_LITERAL(2, 23, 0), // ""
QT_MOC_LITERAL(3, 24, 10), // "sendFailed"
QT_MOC_LITERAL(4, 35, 14), // "receiveSucceed"
QT_MOC_LITERAL(5, 50, 13), // "receiveFailed"
QT_MOC_LITERAL(6, 64, 8), // "setKFile"
QT_MOC_LITERAL(7, 73, 1), // "f"
QT_MOC_LITERAL(8, 75, 8), // "received"
QT_MOC_LITERAL(9, 84, 4) // "tick"

    },
    "Qkx_coeffs\0sendSucceed\0\0sendFailed\0"
    "receiveSucceed\0receiveFailed\0setKFile\0"
    "f\0received\0tick"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Qkx_coeffs[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   49,    2, 0x06 /* Public */,
       3,    0,   50,    2, 0x06 /* Public */,
       4,    0,   51,    2, 0x06 /* Public */,
       5,    0,   52,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       6,    1,   53,    2, 0x0a /* Public */,
       8,    0,   56,    2, 0x08 /* Private */,
       9,    0,   57,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void, QMetaType::QString,    7,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void Qkx_coeffs::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Qkx_coeffs *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->sendSucceed(); break;
        case 1: _t->sendFailed(); break;
        case 2: _t->receiveSucceed(); break;
        case 3: _t->receiveFailed(); break;
        case 4: _t->setKFile((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: _t->received(); break;
        case 6: _t->tick(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (Qkx_coeffs::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Qkx_coeffs::sendSucceed)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (Qkx_coeffs::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Qkx_coeffs::sendFailed)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (Qkx_coeffs::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Qkx_coeffs::receiveSucceed)) {
                *result = 2;
                return;
            }
        }
        {
            using _t = void (Qkx_coeffs::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Qkx_coeffs::receiveFailed)) {
                *result = 3;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject Qkx_coeffs::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_meta_stringdata_Qkx_coeffs.data,
    qt_meta_data_Qkx_coeffs,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *Qkx_coeffs::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Qkx_coeffs::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Qkx_coeffs.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int Qkx_coeffs::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void Qkx_coeffs::sendSucceed()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void Qkx_coeffs::sendFailed()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void Qkx_coeffs::receiveSucceed()
{
    QMetaObject::activate(this, &staticMetaObject, 2, nullptr);
}

// SIGNAL 3
void Qkx_coeffs::receiveFailed()
{
    QMetaObject::activate(this, &staticMetaObject, 3, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
