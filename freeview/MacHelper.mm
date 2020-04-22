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
        if (@available(macOS 10.14, *)) {
            if (QSysInfo::productVersion() < "10.14")
                return false;
            NSAppearance *appearance = NSAppearance.currentAppearance;
            return appearance.name == NSAppearanceNameDarkAqua;
        } else {
            return false;
        }
    }
}

QPixmap MacHelper::InvertPixmap(const QPixmap &pix)
{
    QImage img = pix.toImage();
    img.invertPixels(QImage::InvertRgb);
    return QPixmap::fromImage(img);
}

QIcon MacHelper::InvertIcon(const QIcon &icn_in, const QSize &size_in, bool bTwoStates)
{
    QIcon icn;
    QSize size = size_in;
    if (!size.isValid())
        size = icn_in.availableSizes().first();
    QPixmap pix = icn_in.pixmap(size, QIcon::Normal, QIcon::Off);
    icn.addPixmap(InvertPixmap(pix));
    QPixmap pix2 = icn_in.pixmap(size, QIcon::Normal, QIcon::On);
    if (!pix2.isNull())
        icn.addPixmap(bTwoStates?InvertPixmap(pix2):pix2, QIcon::Normal, QIcon::On);

    return icn;
}
