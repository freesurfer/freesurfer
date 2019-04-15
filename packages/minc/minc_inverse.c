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

#include <stdio.h>
#include "minc_internals.h"

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

    bool success = invert_4x4_matrix( (Double4x4 const *)&t, &inv );

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
