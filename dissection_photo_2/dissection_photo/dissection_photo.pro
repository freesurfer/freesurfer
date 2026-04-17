include(../common.pri)

SOURCES += \
    ../MyUtils.cpp \
    MaskProcessor.cpp \
    DialogWelcome.cpp \
    MakeThumbnailThread.cpp \
    ProgressWindow.cpp \
    QThumbnailWidget.cpp \
    main.cpp \
    MainWindow.cpp

HEADERS += \
    ../MyUtils.h \
    MaskProcessor.h \
    DialogWelcome.h \
    MainWindow.h \
    MakeThumbnailThread.h \
    ProgressWindow.h \
    QThumbnailWidget.h

FORMS += \
    DialogWelcome.ui \
    MainWindow.ui \
    ProgressWindow.ui

RESOURCES += \
  dissection_photo.qrc
