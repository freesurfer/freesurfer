include(../common.pri)

SOURCES += \
    ReconPage.cpp \
    StepperBar.cpp \
    AppMotion.cpp \
    MaskProcessor.cpp \
    SegmentedControl.cpp \
    CompactSlider.cpp \
    ImageGridView.cpp \
    AppTheme.cpp \
    DialogWelcome.cpp \
    DialogProjectsHub.cpp \
    DialogSegmentationSettings.cpp \
    ProjectCardWidget.cpp \
    ProjectManager.cpp \
    MakeThumbnailThread.cpp \
    ProgressWindow.cpp \
    QThumbnailWidget.cpp \
    main.cpp \
    MainWindow.cpp

HEADERS += \
    ReconPage.h \
    StepperBar.h \
    AppMotion.h \
    MaskProcessor.h \
    SegmentedControl.h \
    CompactSlider.h \
    ImageGridView.h \
    AppTheme.h \
    DialogWelcome.h \
    DialogProjectsHub.h \
    DialogSegmentationSettings.h \
    ProjectCardWidget.h \
    ProjectManager.h \
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

# App icon in Finder and the Dock (main.cpp also sets it at runtime).
macx: ICON = resource/dissection_photo.icns
