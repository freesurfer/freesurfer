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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
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
