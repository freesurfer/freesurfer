#-------------------------------------------------
#
# Project created by QtCreator 2014-04-14T13:14:04
#
#-------------------------------------------------

QT       += core gui opengl

TARGET = nmovie_qt
TEMPLATE = app


SOURCES += main.cpp\
        MainWindow.cpp \
    RenderWidget.cpp

HEADERS  += MainWindow.h \
    RenderWidget.h

FORMS    += MainWindow.ui

OTHER_FILES += \
    Makefile.am
