#include <stdlib.h>

#include <QtGlobal>
#include <QFile>
#include <QList>
#include <QDebug>
#include <QByteArray>

#include "exif.h"




Exif::Exif(){
}

Exif::~Exif(){
}

rint8u readByte(QFile &file){
    char a;
    file.getChar(&a);
    return (rint8u)a;
}


//--------------------------------------------------------------------------
// Parse the marker stream until SOS or EOI is seen;
//--------------------------------------------------------------------------
int Exif::readJpegSections(QFile &file, int *Orientation){
    QByteArray *data;
    rint8u a;

    a = readByte(file);

    if (a != 0xff){
            return -1;
    }else{
        rint8u b;
        b = readByte(file);
        if(b != M_SOI) return -1;
    }

    //int SectionsRead=0;
    for(;;){
        int itemlen;
        int prev;
        rint8u  marker = 0;
        rint8u ll,lh;

        prev = 0;
        for(int i=0;;i++){
            marker = readByte(file);
            if (marker != 0xff && prev == 0xff) break;
            prev = marker;
        }


        // Read the length of the section.
        lh = readByte(file);
        ll = readByte(file);
        itemlen = (lh << 8) | ll;

        if (itemlen < 2){ // Invalid marker
            return -1;
        }

        data = new QByteArray(file.read(itemlen-2)); // Read the whole section.
        if(data->isEmpty()) return -1; // Could not allocate memory

        if(data->size() != itemlen-2) return -1; // Premature end of file?
        //sections.append(section);


        switch(marker){
            case M_SOS:   // stop before hitting compressed data
                return 0;

            case M_EOI:   // in case it's a tables-only JPEG stream
                return -1;

            case M_COM: // Comment section
                delete(data);
                break;

           case M_JFIF:
                // Regular jpegs always have this tag, exif images have the exif
                // marker instead, althogh ACDsee will write images with both markers.
                // this program will re-create this marker on absence of exif marker.
                // hence no need to keep the copy from the file.
                if (itemlen >= 16){ // if Jfif header not too short
                    // skipped
                }

                delete(data);
                break;

            case M_EXIF:
                // There can be different section using the same marker.
                if(data->left(4) == "Exif"){
                    processEXIF(data, itemlen, Orientation);
                    break;
                }
                // Oterwise, discard this section.
                delete(data);
                break;

            case M_IPTC:
                delete(data);
                break;

            case M_SOF0:
                    case M_SOF1:
                    case M_SOF2:
                    case M_SOF3:
                    case M_SOF5:
                    case M_SOF6:
                    case M_SOF7:
                    case M_SOF9:
                    case M_SOF10:
                    case M_SOF11:
                    case M_SOF13:
                    case M_SOF14:
                    case M_SOF15:
                        //process_SOFn(Data, marker);
                        break;
            default:
                        // Skip any other sections.
                        break;

        } // switch


    } // for(;;)

    return 0;
}



// Convert a 16 bit unsigned value from file's native byte order
int Get16u(const void * Short, int MotorolaOrder){
    if (MotorolaOrder){
        return (((uchar *)Short)[0] << 8) | ((uchar *)Short)[1];
    }else{
        return (((uchar *)Short)[1] << 8) | ((uchar *)Short)[0];
    }
}


// Convert a 32 bit signed value from file's native byte order
int Get32s(const void * Long, int MotorolaOrder){
    if (MotorolaOrder){
        return  ((( char *)Long)[0] << 24) | (((uchar *)Long)[1] << 16)
              | (((uchar *)Long)[2] << 8 ) | (((uchar *)Long)[3] << 0 );
    }else{
        return  ((( char *)Long)[3] << 24) | (((uchar *)Long)[2] << 16)
              | (((uchar *)Long)[1] << 8 ) | (((uchar *)Long)[0] << 0 );
    }
}

// Convert a 32 bit unsigned value from file's native byte order
unsigned Get32u(const void * Long, int MotorolaOrder){
    return (unsigned)Get32s(Long, MotorolaOrder) & 0xffffffff;
}

#define NUM_FORMATS 12
#define FMT_BYTE       1
#define FMT_STRING     2
#define FMT_USHORT     3
#define FMT_ULONG      4
#define FMT_URATIONAL  5
#define FMT_SBYTE      6
#define FMT_UNDEFINED  7
#define FMT_SSHORT     8
#define FMT_SLONG      9
#define FMT_SRATIONAL 10
#define FMT_SINGLE    11
#define FMT_DOUBLE    12


// Evaluate number, be it int, rational, or float from directory.
double ConvertAnyFormat(const void * ValuePtr, int Format, int MotorolaOrder){
    double Value;
     Value = 0;

     switch(Format){
         case FMT_SBYTE:     Value = *(signed char *)ValuePtr;  break;
         case FMT_BYTE:      Value = *(uchar *)ValuePtr;        break;

         case FMT_USHORT:    Value = Get16u(ValuePtr, MotorolaOrder);          break;
         case FMT_ULONG:     Value = Get32u(ValuePtr, MotorolaOrder);          break;

         case FMT_URATIONAL:
         case FMT_SRATIONAL:
             {
                 int Num, Den;
                 Num = Get32s(ValuePtr, MotorolaOrder);
                 Den = Get32s(4+(char *)ValuePtr, MotorolaOrder);
                 if (Den == 0){
                     Value = 0;
                 }else{
                     Value = (double)Num/Den;
                 }
                 break;
             }

         case FMT_SSHORT:    Value = (signed short)Get16u(ValuePtr, MotorolaOrder);  break;
         case FMT_SLONG:     Value = Get32s(ValuePtr, MotorolaOrder);                break;

         // Not sure if this is correct (never seen float used in Exif format)
         case FMT_SINGLE:    Value = (double)*(float *)ValuePtr;      break;
         case FMT_DOUBLE:    Value = *(double *)ValuePtr;             break;

         default: Value = 100;// Illegal format code


     }
     return Value;
}



#define TAG_ORIENTATION        0x0112
#define TAG_INTEROP_OFFSET     0xA005
#define TAG_EXIF_OFFSET        0x8769
const int BytesPerFormat[] = {0,1,1,2,4,8,1,1,2,4,8,4,8};


// Process one of the nested EXIF directories.
int Exif::processEXIFDir(const char *DirStart, const char *OffsetBase, rint32u exifSize, rint32u nesting, int MotorolaOrder, int *NumOrientations, int *Orientation){
    int numDirEntries;

    if(nesting>4) return -1; // Maximum Exif directory nesting exceeded (corrupt Exif header)

    numDirEntries = Get16u(DirStart, MotorolaOrder);
    //qDebug() << "num entries: " << numDirEntries;

    #define DIR_ENTRY_ADDR(Start, Entry) (Start+2+12*(Entry))

    for (int de=0; de<numDirEntries; de++){
        int Tag, Format, Components;
        const char * DirEntry;
        const char * ValuePtr;
        int ByteCount;

        DirEntry = DIR_ENTRY_ADDR(DirStart, de);
        Tag        = Get16u(DirEntry,   MotorolaOrder);
        Format     = Get16u(DirEntry+2, MotorolaOrder);
        Components = Get32u(DirEntry+4, MotorolaOrder);
        //qDebug() << "tag:" << Tag << "format:" << Format << "components:" << Components;

        if(Format-1 >= NUM_FORMATS) continue; // (-1) catches illegal zero case as unsigned underflows to positive large.
        if((unsigned)Components > 0x10000) continue; // Too many components

        ByteCount = Components * BytesPerFormat[Format];
        //qDebug() << "byte count" << ByteCount;

        if (ByteCount > 4){ // If its bigger than 4 bytes, the dir entry contains an offset.
            unsigned OffsetVal = Get32u(DirEntry+8, MotorolaOrder);
            if (OffsetVal+ByteCount > exifSize) continue; // Bogus pointer offset and / or bytecount value
            ValuePtr = OffsetBase+OffsetVal;
        }else{ // 4 bytes or less and value is in the dir entry itself
            ValuePtr = DirEntry+8;
        }


        // Extract useful components of tag
        switch(Tag){
            case TAG_ORIENTATION:
                if (*NumOrientations >= 2){
                            // Can have another orientation tag for the thumbnail, but if there's
                            // a third one, things are stringae.
                            break;
                }

                if (*NumOrientations == 0){
                    *Orientation = (int)ConvertAnyFormat(ValuePtr, Format, MotorolaOrder);
                    //qDebug() << "orientation:" << *Orientation;
                }
                if (*Orientation < 0 || *Orientation > 8){
                    // Undefined rotation value
                    *Orientation = 0;
                }
                *NumOrientations += 1;
                break;

            case TAG_EXIF_OFFSET:

            case TAG_INTEROP_OFFSET:
                const char * SubdirStart;
                SubdirStart = OffsetBase + Get32u(ValuePtr, MotorolaOrder);
                if (SubdirStart < OffsetBase || SubdirStart > OffsetBase+ exifSize){
                    // Illegal Exif or interop ofset directory link
                }else{
                    processEXIFDir(SubdirStart, OffsetBase, exifSize, nesting+1, MotorolaOrder, NumOrientations, Orientation);
                }
                continue;
                break;

            default:
                    // Skip any other sections.
                    break;

        }

    }

    return 0;
}


// Process a EXIF marker
// Describes all the drivel that most digital cameras include...
int Exif::processEXIF(QByteArray *data, int itemlen, int *Orientation){
    int MotorolaOrder = 0;

    // Check the EXIF header component
    if(data->left(6) == "Exif\0\0"){
        qDebug() << data->left(4);
    }

    if(data->mid(6,2) == "II"){ // Exif section in Intel order
        //qDebug() << data->mid(6,2);
        MotorolaOrder = 0;
    }else{
        if(data->mid(6,2) == "MM"){ // Exif section in Motorola order
            //qDebug() << data->mid(6,2);
            MotorolaOrder = 1;
        }else{
            return -1; // Invalid Exif alignment marker.
        }
    }

    // get first offset
    QByteArray ttt(data->mid(10,4));
    const char *ttt2 = ttt.constData();
    rint32u FirstOffset;
    FirstOffset = Get32u(ttt2, MotorolaOrder);
    //qDebug() << "fist offset: " << FirstOffset;
    if (FirstOffset < 8 || FirstOffset > 16){
            if (FirstOffset < 16 || int(FirstOffset) > itemlen-16)  return -1;  // invalid offset for first Exif IFD value ;
    }

    const char *dirStart = data->constData();
    const char *offsetBase = data->constData();

    dirStart   += 6 + FirstOffset;
    offsetBase += 6;
    int numOrientations = 0;
    // First directory starts 16 bytes in.  All offset are relative to 8 bytes in.
    processEXIFDir(dirStart, offsetBase, itemlen-8, 0, MotorolaOrder, &numOrientations, Orientation);
    //qDebug() << "num orientations:" << numOrientations;

    return 0;
}

int Exif::readJpegFile(QFile &file, int *Orientation){
    readJpegSections(file, Orientation);

    return 0;
}



int Exif::getExifOrientation(QFile &file, int *Orientation){
    readJpegFile(file, Orientation);
    return 0;
}
