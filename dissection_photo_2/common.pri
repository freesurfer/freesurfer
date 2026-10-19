QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

INCLUDEPATH += $$PWD
SOURCES += $$PWD/WidgetImageView.cpp $$PWD/cnpy.cpp $$PWD/exif.cpp $$PWD/MyUtils.cpp
HEADERS += $$PWD/WidgetImageView.h $$PWD/cnpy.h \
  $$PWD/CommonDef.h $$PWD/exif.h $$PWD/MyUtils.h
LIBS += -lz

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
