/*
 * Overhaul Author: Bevin Brett
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

// BASED ON minc-1.5.1/volume_io/Geometry/gaussian.c

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

#include  "minc_internals.h"

#if !defined(TESTING_GAUSSIAN_C)

static const bool doPrint = false;
    
#else
    
#include <stdbool.h>
typedef double Double4x4[4*4];
#define Index4x4(I,J) (4*(I)+(J))
static const bool doPrint = true;

#endif


#include <math.h>
#include <stdio.h>


static void print(const char* heading, double* mat, int* mapI) {
    printf("%s\n",heading);
    int i,j;
    for (i = 0; i < 4; i++) {
        int mi = mapI ? mapI[i] : i;
        for (j = 0; j < 4; j++) 
            printf("%lf ",mat[Index4x4(mi,j)]);
        printf("\n");
    }
}

#if defined(TESTING_GAUSSIAN_C)
static const char* invert_4x4_matrix_wkr_simple( double* mat, double* inv, int* mapI )
{
    printf("doing inversion simply\n");
    
    int i, i2, j;

    print("mat being inverted", mat, mapI);

    double scaleRow[4];
    for( i = 0; i < 4; i++ ) {
        double max   = fabs(mat[Index4x4(i,0)]);
        for( j = 1; j < 4; j++ ) {
	    double d = fabs(mat[Index4x4(i,j)]);
	    if (max < d) max = d;
        }
	if (max == 0.0) return "Some row is all zeros";
	scaleRow[i] = 1.0 / max;
    }

    for( i = 0; i < 4; i++ ) {
        for( j = 0; j < 4; j++ ) inv[Index4x4(i,j)] = 0.0;
        inv[Index4x4(i,i)] = 1.0;
    }
    print("inv set to id", inv, mapI);

    for( i = 0; i < 4; i++ ) {
        printf("\n\ndoing col %d\n",i);
    
        // choose the best of the remaining rows to use
	//
    	int    mi    = mapI[i];
        int    maxI  = i;
        double max   = fabs(mat[Index4x4(mi, i)] * scaleRow[mi ]);
        for( i2 = i+1; i2 < 4; i2++) {
	    const int mi2 = mapI[i2];
            double d = fabs(mat[Index4x4(mi2,i)] * scaleRow[mi2]);
            if (max < d) {
	        printf("max:%f maxI:%d replaced by %f %d\n", max,maxI, d, i2); 
	        max = d; 
		maxI = i2; 
	    }
        }
        if (max == 0.0) return "Some col is all zeros";

        // swap the best row to be the i'th entry in the map
	//
        {   // it is better to swap than mispredict a test-and-branch
            int tmp    = mapI[maxI];
            mapI[maxI] = mapI[i];
            mapI[i]    = tmp;
        }
		
	// apply the elimination step to all the later rows
	// and to all the rows of the inv
	//
	mi = mapI[i];
        printf("row %d chosen\n",mi);
	
	double scale = 1.0 / mat[Index4x4(mi,i)];
	    // max is the fabs of mat[...] and is known > 0
	
        {
	    for( j = 0; 
	        j < 4; j++ )
                mat[Index4x4(mi,j)] *= scale;
		
            for( j =   0; j < 4; j++ )
                inv[Index4x4(mi,j)] *= scale;
        }

        for( i2 = 0; i2 < 4; i2++ ) {
	    if (i2 == i) continue;
	    const int mi2 = mapI[i2];
	    
            double m = mat[Index4x4(mi2,i)];
            for( j = 0; j < 4; j++ )
                mat[Index4x4(mi2,j)] -= m * mat[Index4x4(mi,j)];
		
            for( j = 0; j < 4; j++ )
                inv[Index4x4(mi2,j)] -= m * inv[Index4x4(mi,j)];
        }
	printf("After elimination:\n");
        print("    mat", mat, NULL);
        print("    inv", inv, NULL);
        print("    mat reordered", mat, mapI);
        print("    inv reordered", inv, mapI);
    }

    return NULL;
}
#endif


static const char* invert_4x4_matrix_wkr_fast( double* mat, double* inv, int* mapI )
{
    if (doPrint) printf("doing inversion optimized\n");
    
    int i, i2, j;

    if (doPrint) print("mat being inverted", mat, mapI);

    double scaleRow[4];
    for( i = 0; i < 4; i++ ) {
        double max   = fabs(mat[Index4x4(i,0)]);
        for( j = 1; j < 4; j++ ) {
	    double d = fabs(mat[Index4x4(i,j)]);
	    if (max < d) max = d;
        }
	if (max == 0.0) return "Some row is all zeros";
	scaleRow[i] = 1.0 / max;
    }

    for( i = 0; i < 4; i++ ) {
        for( j = 0; j < 4; j++ ) inv[Index4x4(i,j)] = 0.0;
        inv[Index4x4(i,i)] = 1.0;
    }
    if (doPrint) print("inv set to id", inv, mapI);

    for( i = 0; i < 3; i++ ) {
        if (doPrint) printf("\n\ndoing col %d\n", i);
    
        // choose the best of the remaining rows to use
	//
    	int    mi    = mapI[i];
        int    maxI  = i;
        double max   = fabs(mat[Index4x4(mi, i)] * scaleRow[mi ]);
        for( i2 = i+1; i2 < 4; i2++) {
	    const int mi2 = mapI[i2];
            double d = fabs(mat[Index4x4(mi2,i)] * scaleRow[mi2]);
            if (max < d) {
	        if (doPrint) printf("max:%f maxI:%d replaced by %f %d\n", max,maxI, d, i2); 
	        max = d; 
		maxI = i2; 
	    }
        }
        if (max == 0.0) return "Some col is all zeros";

        // swap the best row to be the i'th entry in the map
	//
        {   // it is better to swap than mispredict a test-and-branch
            int tmp    = mapI[maxI];
            mapI[maxI] = mapI[i];
            mapI[i]    = tmp;
        }
		
	// apply the elimination step to all the later rows
	// and to all the rows of the inv
	//
	mi = mapI[i];
        if (doPrint) printf("row %d chosen\n",mi);
	
	double scale = 1.0 / mat[Index4x4(mi,i)];
	    // max is the fabs of mat[...] and is known > 0
	
        for( i2 = i+1; i2 < 4; i2++ ) {
	    const int mi2 = mapI[i2];
	    
            double m = mat[Index4x4(mi2,i)] * scale;
            for( j = i+1; j < 4; j++ )
                mat[Index4x4(mi2,j)] -= m * mat[Index4x4(mi,j)];
		
            for( j =   0; j < 4; j++ )
                inv[Index4x4(mi2,j)] -= m * inv[Index4x4(mi,j)];
        }
	if (doPrint) {
	    printf("After elimination:\n");
            print("    mat", mat, NULL);
            print("    inv", inv, NULL);
            print("    mat reordered", mat, mapI);
            print("    inv reordered", inv, mapI);
	}
    }

    if( mat[Index4x4(mapI[3],3)] == 0.0 )
        return "last col all zeros";

    for( i = 3; i >= 0; --i ) {
        int mi = mapI[i];
        for( i2 = i+1; i2 < 4; i2++ ) {
            double scale = mat[Index4x4(mi,i2)];
            int mi2 = mapI[i2];
            for( j = 0; j < 4; j++ ) {
               inv[Index4x4(mi,j)] -= scale * inv[Index4x4(mi2,j)];
            }
        }
        double scale = 1.0 / mat[Index4x4(mi,i)];
        for( j = 0; j < 4; j++)
            inv[Index4x4(mi,j)] *= scale;
    }

    return NULL;
}

bool invert_4x4_matrix( 
	const Double4x4* mat, 
	Double4x4* inv ) 
{
    Double4x4 copyToBeDestroyed;
    Double4x4 permutedInv;

    int i,j;
    for (i = 0; i < 4; i++) 
        for (j = 0; j < 4; j++) 
            copyToBeDestroyed[Index4x4(i,j)] = (*mat)[Index4x4(i,j)];

    int mapI[4];
    for( i = 0; i < 4; i++ ) {
        mapI[i] = i;
    }

    const char* errorMessage = 
        invert_4x4_matrix_wkr_fast(&copyToBeDestroyed[0], &permutedInv[0], &mapI[0] );
	
    if (errorMessage) {
        printf("invert_4x4_matrix failed, %s\n",errorMessage);
    } else {
        for (i = 0; i < 4; i++) {
            int mi = mapI[i];
            for (j = 0; j < 4; j++) 
                (*inv)[Index4x4(i,j)] = permutedInv[Index4x4(mi,j)];
	}
    }
    return errorMessage == NULL;
}


#if defined(TESTING_GAUSSIAN_C)
#include <stdio.h>
int main() {
    Double4x4 mat, inv, mul;
    printf("testing %s\n", __FILE__);
    int i,j,k;

    for (i = 0; i < 4; i++) 
        for (j = 0; j < 4; j++) 
            mat[Index4x4(i,j)] = (double)(i != j) + (i < j);
	    
    if (!invert_4x4_matrix((const Double4x4*)&mat, &inv)) {
        printf("inv failed\n");
    } else {
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
	        double sum = 0.0;
                for (k = 0; k < 4; k++) 
	            sum += mat[Index4x4(i,k)]*inv[Index4x4(k,j)];
	        mul[Index4x4(i,j)] = sum;
	    }
        }
    }
    print("Final inv",&inv[0],NULL);
    print("Final mul",&mul[0],NULL);
    return 0;
}

#endif
