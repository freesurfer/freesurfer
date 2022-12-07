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
  fiducials_calibration.qrc

TARGET = fiducials_calibration
