!greaterThan(QT_MAJOR_VERSION, 4): {

    INCLUDEPATH += \
        $$PWD

    SOURCES += \
        $$PWD/qjson.cpp \
        $$PWD/qjsonarray.cpp \
        $$PWD/qjsondocument.cpp \
        $$PWD/qjsonobject.cpp \
        $$PWD/qjsonparser.cpp \
        $$PWD/qjsonvalue.cpp \
        $$PWD/qjsonwriter.cpp

    HEADERS += \
        $$PWD/qjson_p.h \
        $$PWD/qjsonarray.h \
        $$PWD/qjsondocument.h \
        $$PWD/qjsonobject.h \
        $$PWD/qjsonparser_p.h \
        $$PWD/qjsonvalue.h \
        $$PWD/qjsonwriter_p.h

}
