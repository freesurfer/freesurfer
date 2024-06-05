include(../common.pri)

SOURCES += \
    ../MyUtils.cpp \
    ../connected_components/MaskProcessor.cpp \
    ProgressWindow.cpp \
    main.cpp \
    MainWindow.cpp

HEADERS += \
    ../MyUtils.h \
    ../connected_components/MaskProcessor.h \
    MainWindow.h \
    ProgressWindow.h

FORMS += \
    MainWindow.ui \
    ProgressWindow.ui

RESOURCES += \
  dissection_photo.qrc
