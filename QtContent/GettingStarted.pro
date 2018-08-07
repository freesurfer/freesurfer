TEMPLATE = app
TARGET = TRACULA_AnatomiCuts_Tool

QT = core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    configurationfileform.cpp

HEADERS += \
    mainwindow.h \
    configurationfileform.h

FORMS += \
    configurationfileform.ui
