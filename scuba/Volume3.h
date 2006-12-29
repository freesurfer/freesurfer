/**
 * @file  Volume3.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef Volume3_h
#define Volume3_h


template <typename T>
class Volume3 {
public:
  Volume3 ( int izX, int izY, int izZ, T iInitialValue );
  ~Volume3 ();

  void SetAll ( T iValue );

  T Get ( int inX, int inY, int inZ );
  void Set ( int inX, int inY, int inZ, T iValue );

  T Get_Unsafe ( int inX, int inY, int inZ );
  void Set_Unsafe ( int inX, int inY, int inZ, T iValue );

protected:
  T*** mData;
  int mzX;
  int mzY;
  int mzZ;
};


#endif
