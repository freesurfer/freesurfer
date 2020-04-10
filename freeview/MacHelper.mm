#include "MacHelper.h"
#include <AppKit/NSWindow.h>
#include <AppKit/NSApplication.h>
#include <AppKit/NSAppearance.h>
#include <StoreKit/StoreKit.h>
#include <QSysInfo>
#include <QDebug>


MacHelper::MacHelper(QObject* parent) : QObject(parent)
{
}

bool MacHelper::IsDarkMode()
{
    @autoreleasepool {
        if (QSysInfo::productVersion() < "10.14")
            return false;
        NSAppearance *appearance = NSAppearance.currentAppearance;
            return appearance.name == NSAppearanceNameDarkAqua;
    }
}
