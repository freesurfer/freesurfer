// Licensed under MIT license; see license.txt.

#include <sbl/image/MotionFieldSeq.h>
#include <sbl/core/StringUtil.h>
#include <sbl/image/MotionFieldUtil.h>
#include <sbl/system/FileSystem.h>
namespace sbl {


//-------------------------------------------
// MOTION FIELD SEQ CLASS
//-------------------------------------------


/// specify the path containing (or that will contain) the motion field files
MotionFieldSeq::MotionFieldSeq( const String &path ) {
    m_path = path;
    m_count = 0;
    Array<String> dirList = dirFileList( path, "motion.", ".mf" );
    for (int i = 0; i < dirList.count(); i++) {
        Array<String> split = dirList[ i ].split( "." );
        if (split.count() == 4) {
            int srcIndex = split[ 2 ].toInt();
            int destIndex = split[ 3 ].toInt();
            if (srcIndex > m_count)
                m_count = srcIndex + 1;
            if (destIndex > m_count)
                m_count = destIndex + 1;
        }
    }
}


/// read a motion field from the sequence
aptr<MotionField> MotionFieldSeq::read( int index ) {
    String fileName = sprintF( "motion.%05d.%05d.mf", index, index + 1 );
    return loadMotionField( fileName );
}


/// write a motion field into the sequence
void MotionFieldSeq::write( int index, MotionField &motionField ) {
    String fileName = sprintF( "motion.%05d.%05d.mf", index, index + 1 );
    saveMotionField( motionField, fileName );
}


} // end namespace sbl

