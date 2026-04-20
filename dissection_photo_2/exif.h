#ifndef EXIF_H
#define EXIF_H

// This implementation is based on http://www.sentex.net/~mwandel/jhead/
// Rewritten and published in public domain like
// the original code by http://imonad.com

//--------------------------------------------------------------------------
// JPEG markers consist of one or more 0xFF bytes, followed by a marker
// code byte (which is not an FF).  Here are the marker codes of interest
// in this program.  (See jdmarker.c for a more complete list.)
//--------------------------------------------------------------------------

#define M_SOF0  0xC0          // Start Of Frame N
#define M_SOF1  0xC1          // N indicates which compression process
#define M_SOF2  0xC2          // Only SOF0-SOF2 are now in common use
#define M_SOF3  0xC3
#define M_SOF5  0xC5          // NB: codes C4 and CC are NOT SOF markers
#define M_SOF6  0xC6
#define M_SOF7  0xC7
#define M_SOF9  0xC9
#define M_SOF10 0xCA
#define M_SOF11 0xCB
#define M_SOF13 0xCD
#define M_SOF14 0xCE
#define M_SOF15 0xCF
#define M_SOI   0xD8          // Start Of Image (beginning of datastream)
#define M_EOI   0xD9          // End Of Image (end of datastream)
#define M_SOS   0xDA          // Start Of Scan (begins compressed data)
#define M_JFIF  0xE0          // Jfif marker
#define M_EXIF  0xE1          // Exif marker.  Also used for XMP data!
#define M_XMP   0x10E1        // Not a real tag (same value in file as Exif!)
#define M_COM   0xFE          // COMment
#define M_DQT   0xDB
#define M_DHT   0xC4
#define M_DRI   0xDD
#define M_IPTC  0xED          // IPTC marker


typedef unsigned char rint8u;
typedef          char rint8;
typedef unsigned int  rint32u;


#include <QList>
#include <QByteArray>


class Exif{
public:
    Exif();
    ~Exif();
    int getExifOrientation(QFile &file, int *Orientation);
    int readJpegFile(QFile &file, int *Orientation);
    int readJpegSections(QFile &file, int *Orientation);
    int processEXIF(QByteArray *barr, int itemlen, int *Orientation);
    int processEXIFDir(const char *dirStart, const char *offsetBase, rint32u size, rint32u nesting, int MotorolaOrder, int *numOrientations, int *Orientation);
};



#endif // EXIF_H
