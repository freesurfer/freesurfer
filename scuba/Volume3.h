/**
 * @file  Volume3.h
 * @brief Simple class for a 3D volume of elements
 *
 * Template class for creating 3D volumes of a data type with settors
 * and accessors.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/17 23:59:49 $
 *    $Revision: 1.6 $
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
  Volume3 ( int izX, int izY, int izZ, T const iInitialValue );
  ~Volume3 ();

  // Accessors. Unsafe version doesn't check bounds.
  T Get ( int inX, int inY, int inZ ) const;
  T Get_Unsafe ( int inX, int inY, int inZ ) const;
  void GetBounds( int& ozX, int& ozY, int& ozZ ) const;

  // Settors. Unsafe version doesn't check bounds.
  void SetAll ( T const iValue );
  void Set ( int inX, int inY, int inZ, T const iValue );
  void Set_Unsafe ( int inX, int inY, int inZ, T constiValue );

protected:
  T*** mData;
  int mzX;
  int mzY;
  int mzZ;
};


#endif
