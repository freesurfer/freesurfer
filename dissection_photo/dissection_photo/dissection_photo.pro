include(../common.pri)

SOURCES += \
    ../connected_components/MaskProcessor.cpp \
    main.cpp \
    MainWindow.cpp

HEADERS += \
    ../connected_components/MaskProcessor.h \
    MainWindow.h

FORMS += \
    MainWindow.ui

RESOURCES += \
  dissection_photo.qrc
