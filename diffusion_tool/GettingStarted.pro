TEMPLATE = app
TARGET = TRACULA_AnatomiCuts_Tool

QT = core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

SOURCES += \
    main.cpp \
    configurationfileform.cpp

HEADERS += \
    configurationfileform.h

FORMS += \
    configurationfileform.ui
