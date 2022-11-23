include(../common.pri)

SOURCES += \
    DialogWelcome.cpp \
    main.cpp \
    MainWindow.cpp

HEADERS += \
    DialogWelcome.h \
    MainWindow.h

FORMS += \
    DialogWelcome.ui \
    MainWindow.ui

RESOURCES += \
  calibration.qrc

TARGET = calibration_with_fiducials
