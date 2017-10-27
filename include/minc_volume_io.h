/**
 * @file  minc_volume_io.h
 * @brief Wrapper for MNI's volume_io.h, to decouple from MNI lib
 *
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.5 $
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


#ifndef MINC_VOLUME_IO_H
#define MINC_VOLUME_IO_H

#include <stdbool.h>
#include <stdlib.h>

#if defined(BEVIN_EXCLUDE_MINC) || defined(BEVIN_REPLACE_MINC)

// The following is a replacement for some portions of
// mni/1.5/include/minc.h
//
// As such, it needs the following Copyright notice
//
/*
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
*/

#define MIxspace "xspace"
#define MIyspace "yspace"
#define MIzspace "zspace"


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/geom_structs.h
//
// As such, it needs the following Copyright notice
//
/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
@VERSION    : $Header: /private-cvsroot/minc/volume_io/Include/volume_io/geom_structs.h,v 1.20.2.3 2006/11/30 09:15:13 rotor Exp $
---------------------------------------------------------------------------- */

typedef struct
{
    double m[4][4];
} Transform;


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/transforms.h
//
// As such, it needs the following Copyright notice
//
/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
@VERSION    : $Header: /private-cvsroot/minc/volume_io/Include/volume_io/transforms.h,v 1.12.2.2 2005/03/31 17:39:49 bert Exp $
---------------------------------------------------------------------------- */

typedef  void   (*User_transform_function)( void  *user_data,
                                            double  x,
                                            double  y,
                                            double  z,
                                            double  *x_trans,
                                            double  *y_trans,
                                            double  *z_trans );

typedef enum { 
	LINEAR
	//, 
	//THIN_PLATE_SPLINE, 
	//USER_TRANSFORM,
        //CONCATENATED_TRANSFORM, 
	//GRID_TRANSFORM 
} Transform_types;


typedef struct General_transform
{
    Transform_types             type;
    bool                    	inverse_flag;

    /* --- linear transform */

    Transform               *linear_transform;
    Transform               *inverse_linear_transform;

    /* --- non-linear transform */

    int                         n_points;
    int                         n_dimensions;
    double                    **points;
    double                    **displacements;   /* n_points + n_dim + 1 by */
                                                   /* n_dim */

    /* --- grid transform */

    void                        *displacement_volume;

    /* --- user_defined */

    void                        *user_data;
    size_t                      size_user_data;
    User_transform_function     user_transform_function;
    User_transform_function     user_inverse_transform_function;

    /* --- concatenated transform */

    int                         n_transforms;
    struct General_transform    *transforms;

} General_transform;


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/vol_io_prototypes.h
// which did not have its own Copyright notice
//
Transform* get_linear_transform_ptr(
    General_transform   *transform );

Transform* get_inverse_linear_transform_ptr(
    General_transform   *transform );

void copy_general_transform(
    General_transform   *transform,
    General_transform   *copy );

void  delete_general_transform(
    General_transform   *transform );

static const int OK = 0;

int input_transform_file(			// returns OK or <what?>
    const char* filename,
    General_transform   *transform );

void transform_point(
    Transform  *transform,
    double 	x,
    double	y,
    double 	z,
    double	*x_trans,
    double	*y_trans,
    double	*z_trans );


#else


/*
 * Wrapper for MNI's volume_io.h, which has some annoyances which
 * must be circumvented.
 */

/* remove unwanted warnings between hips_basic.h vs. volume_io/basic.h */
#undef ABS
#undef SIGN
#ifdef Darwin
// The result of not defining __MACTYPES__ is a scuba2 build error
// complaining about QT conflicting with MNI over the Point typedef.
// Notes from the file <mni installation>/include/volume_io/geom_structs.h:
/* Th 'Point' typedef is annoying to Mac OS users, since Point has been
 * a basic type on Macs since the beginning.  Testing __MACTYPES__ should
 * work at least with the OS X codebase, I don't know if it existed in
 * earlier versions of the MacTypes.h header.
 */
#define __MACTYPES__
#endif
#ifdef Status
// avoid conflicts with usage of 'Status' in volume_io/basic.h
#undef Status
#endif
#ifdef Windows_NT
#undef ERROR
#endif // Windows_NT

#include <volume_io.h> //from MNI

/* remove unwanted warnings between hips_basic.h vs. volume_io/basic.h */
#undef ABS
#undef SIGN

#endif

#endif // MINC_VOLUME_IO_H
