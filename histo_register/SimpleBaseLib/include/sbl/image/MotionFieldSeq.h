#ifndef _SBL_MOTION_FIELD_SEQ_H_
#define _SBL_MOTION_FIELD_SEQ_H_
#include <sbl/core/String.h>
#include <sbl/core/Pointer.h>
#include <sbl/image/MotionField.h>
namespace sbl {


//-------------------------------------------
// MOTION FIELD SEQ CLASS
//-------------------------------------------


/// The MotionFieldSeq class provides access to motion fields estimated for an image sequence;
/// the motion fields are stored on disk.
class MotionFieldSeq {
public:

    /// specify the path containing (or that will contain) the motion field files
    MotionFieldSeq( const String &path );
    
    /// read a motion field from the sequence
    aptr<MotionField> read( int index );

    /// write a motion field into the sequence
    void write( int index, MotionField &motionField );

    /// the effective sequence length (may include gaps)
    inline int count() const { return m_count; }

private:

    // the path containing (or that will contain) the motion field files
    String m_path;

    // the effective sequence length (may include gaps)
    int m_count;

    // disable copy constructor and assignment operator
    MotionFieldSeq( const MotionFieldSeq &x );
    MotionFieldSeq &operator=( const MotionFieldSeq &x );
};


} // end namespace sbl
#endif // _SBL_MOTION_FIELD_SEQ_H_

