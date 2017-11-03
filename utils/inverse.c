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
---------------------------------------------------------------------------- */

#if defined(BEVIN_EXCLUDE_MINC)

#include  "minc_volume_io.h"
#include <stdio.h>

bool compute_transform_inverse(
    Transform  *transform,
    Transform  *inverse )
{
    Double4x4 t, inv;

    {
        int i; 
	for( i=0; i < 4; i++ ) {
           int j; 
	   for( j = 0; j < 4; j++ ) {
              t[Index4x4(i,j)] = Transform_elem(*transform,i,j);
	   }
	}
    }

    bool success = invert_4x4_matrix( &t, &inv );

    if( success )
    {
        int i; 
	for( i=0; i < 4; i++ ) {
           int j; 
	   for( j = 0; j < 4; j++ ) {
                Transform_elem(*inverse,i,j) = inv[Index4x4(i,j)];
           }
        }

        static int count = 0;
	if (count++ == 0) {
            Transform  ident;

            concat_transforms( &ident, transform, inverse );

            if( !close_to_identity(&ident) )
            {
                fprintf(stderr, "Error in compute_transform_inverse\n" );
		exit(1);
            }
        }
    }
    else
        make_identity_transform( inverse );

    return( success );
}


#else
#include  <internal_volume_io.h>


#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/minc/volume_io/Geometry/inverse.c,v 1.8.2.1 2004/10/04 20:18:40 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_transform_inverse
@INPUT      : transform
@OUTPUT     : inverse
@RETURNS    : TRUE if successful
@DESCRIPTION: Computes the inverse of the given transformation matrix.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  BOOLEAN   compute_transform_inverse(
    Transform  *transform,
    Transform  *inverse )
{
    int        i, j;
    Real       **t, **inv;
    BOOLEAN    success;

    /* --- copy the transform to a numerical recipes type matrix */

    ALLOC2D( t, 4, 4 );
    ALLOC2D( inv, 4, 4 );

    for_less( i, 0, 4 )
    {
        for_less( j, 0, 4 )
            t[i][j] = Transform_elem(*transform,i,j);
    }

    success = invert_square_matrix( 4, t, inv );

    if( success )
    {
        /* --- copy the resulting numerical recipes matrix to the
               output argument */

        for_less( i, 0, 4 )
        {
            for_less( j, 0, 4 )
            {
                Transform_elem(*inverse,i,j) = inv[i][j];
            }
        }

#ifdef  DEBUG
        /* --- check if this really is an inverse, by multiplying */

        {
            Transform  ident;

            concat_transforms( &ident, transform, inverse );

            if( !close_to_identity(&ident) )
            {
                print_error( "Error in compute_transform_inverse\n" );
            }
        }
#endif
    }
    else
        make_identity_transform( inverse );

    FREE2D( t );
    FREE2D( inv );

    return( success );
}

#endif
