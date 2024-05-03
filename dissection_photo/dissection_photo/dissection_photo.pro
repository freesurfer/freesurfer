include(../common.pri)

SOURCES += \
    ../MyUtils.cpp \
    ../connected_components/MaskProcessor.cpp \
    main.cpp \
    MainWindow.cpp

HEADERS += \
    ../MyUtils.h \
    ../connected_components/MaskProcessor.h \
    MainWindow.h

FORMS += \
    MainWindow.ui

RESOURCES += \
  dissection_photo.qrc
