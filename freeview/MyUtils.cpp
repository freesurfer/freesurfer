/**
 * @file  MyUtils.h
 * @brief Misc utility class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/17 00:31:25 $
 *    $Revision: 1.9 $
 *
 * Copyright (C) 2002-2009,
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

#include "MyUtils.h"
#include "vtkFSVolumeSource.h"
#include <math.h>
#include <wx/filename.h>
#include <stddef.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageWriter.h>
#include <vtkBMPWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkRenderLargeImage.h>
#include <vtkVRMLExporter.h>

extern "C" {
#include "matrix.h"
}

typedef struct { int xl, xr, y, dy; } LINESEGMENT;

#define MAXDEPTH 10000

#define PUSH(XL, XR, Y, DY) \
    if( sp < stack+MAXDEPTH && Y+(DY) >= min_y && Y+(DY) <= max_y ) \
{ sp->xl = XL; sp->xr = XR; sp->y = Y; sp->dy = DY; sp++; }

#define POP(XL, XR, Y, DY) \
{ sp--; XL = sp->xl; XR = sp->xr; Y = sp->y + (DY = sp->dy); }

void MyUtils::FloodFill(char** data, int x, int y, int min_x, int min_y, int max_x, int max_y, int fill_value, int border_value)
{
	int left, x1, x2, dy;
	LINESEGMENT stack[MAXDEPTH], *sp = stack;

	if (data[y][x] == border_value || data[y][x] == fill_value)
		return;

	if (x < min_x || x > max_x || y < min_y || y > max_y)
		return;

	PUSH(x, x, y, 1);        /* needed in some cases */
	PUSH(x, x, y+1, -1);    /* seed segment (popped 1st) */

	while (sp > stack ) 
	{
		POP(x1, x2, y, dy);

		for (x = x1; x >= min_x && data[y][x] != border_value && data[y][x] != fill_value; x--)
			data[y][x] = fill_value;

		if( x >= x1 )
			goto SKIP;

		left = x+1;
		if( left < x1 )
			PUSH(left, x1-1, y, -dy);    /* leak on left? */

		x = x1+1;

		do 
		{
			for (; x<=max_x && data[y][x] != border_value && data[y][x] != fill_value; x++)
				data[y][x] = fill_value;

			PUSH(left, x-1, y, dy);

			if (x > x2+1)
				PUSH(x2+1, x-1, y, -dy);    /* leak on right? */

SKIP:		for (x++; x <= x2 && (data[y][x] == border_value || data[y][x] == fill_value); x++) 
			{;}

			left = x;
		} 
		while (x <= x2);
	}
}

bool MyUtils::HasExtension( const wxString& filename, const wxString& ext )
{
	return ( filename.Lower().Right( ext.Len() ) == ext.Lower() );
}

wxString MyUtils::GetNormalizedPath( const wxString& filename )
{
	wxFileName fn( filename );
	fn.Normalize( wxPATH_NORM_ENV_VARS | wxPATH_NORM_DOTS | wxPATH_NORM_ABSOLUTE );
	return fn.GetPath();
}


wxString MyUtils::GetNormalizedFullPath( const wxString& filename )
{
	wxFileName fn( filename );
	fn.Normalize( wxPATH_NORM_ENV_VARS | wxPATH_NORM_DOTS | wxPATH_NORM_ABSOLUTE );
	return fn.GetFullPath();
}

wxArrayString MyUtils::SplitString( const wxString& strg_to_split, const wxString& divider )
{
	wxArrayString sa;
	wxString strg = strg_to_split;
	strg.Trim( true );
	strg.Trim( false );
	int n = strg.Find( divider );
	while ( n != wxNOT_FOUND )
	{
		wxString substr = strg.Left( n );
		substr.Trim( true );
		substr.Trim( false );
		if ( substr.Length() > 0 )
			sa.Add( substr );
		strg = strg.Mid( n + divider.Length() );
		strg.Trim( true );
		strg.Trim( false );
		n = strg.Find( divider );
	}
	if ( strg.Length() > 0 )
		sa.Add( strg );
	
	return sa;
}

wxString MyUtils::GetDateAndTime()
{
	wxString strg = wxString( __DATE__) + " " + __TIME__;
	
	return strg;
}

bool MyUtils::VTKScreenCapture( vtkRenderWindow* renderWnd, vtkRenderer* renderer, 
								const char* filename, bool bAntiAliasing, int nMag )
{
	wxString fn = filename;
	vtkImageWriter* writer = 0;
	if ( HasExtension( filename, "wrl" ) )
	{
		vtkVRMLExporter* exporter = vtkVRMLExporter::New();
		exporter->SetFileName( filename );
		exporter->SetRenderWindow( renderWnd );
		exporter->Write();
		exporter->Delete();
	}
	else if ( HasExtension( filename, "jpg" ) || HasExtension( filename, "jpeg" ) )
		writer = vtkJPEGWriter::New();
	else if ( HasExtension( filename, "bmp" ) )
		writer = vtkBMPWriter::New();
	else if ( HasExtension( filename, "ps" ) )
		writer = vtkPostScriptWriter::New();
	else if ( HasExtension( filename, "tif" ) || HasExtension( filename, "tiff" ) )
		writer = vtkTIFFWriter::New();
	else 
	{
		writer = vtkPNGWriter::New();
		if ( !HasExtension( filename, "png") )
			fn += ".png";
	}
	if (writer)
	{
	//	bool bCurrentAA = GetAntialiasing() > 0;
	//	SetAntialiasing(bAntiAliasing, false);
		vtkRenderLargeImage* image = vtkRenderLargeImage::New();
		image->SetInput( renderer );
		image->SetMagnification( nMag );
		writer->SetInput( image->GetOutput() );
		writer->SetFileName( fn.c_str() );
		writer->Write();
		image->Delete();
		writer->Delete();
	//	SetAntialiasing(bCurrentAA, false);
	}
	return true;	
}


void MyUtils::ViewportToWorld( vtkRenderer* renderer, double x, double y, double& world_x, double& world_y, double& world_z )
{
	world_x = x;
	world_y = y;
	renderer->ViewportToNormalizedViewport( world_x, world_y );
	NormalizedViewportToWorld( renderer, world_x, world_y, world_x, world_y, world_z );
}

void MyUtils::NormalizedViewportToWorld( vtkRenderer* renderer, double x, double y, double& world_x, double& world_y, double& world_z )
{
	world_x = x;
	world_y = y;
	world_z = 0;
	renderer->NormalizedViewportToView( world_x, world_y, world_z );
	renderer->ViewToWorld( world_x, world_y, world_z );
}

template <class T> bool CalculateOptimalVolume_t( int* vox1, int nsize1, int* vox2, int nsize2,
			std::vector<void*> input_volumes, T* output_volume, int vol_size )
{
	int nvars = input_volumes.size();
	MATRIX* m1 = MatrixAlloc( nsize1, nvars, MATRIX_REAL );
	MATRIX* m2 = MatrixAlloc( nsize2, nvars, MATRIX_REAL );
	for ( int i = 0; i < nvars; i++ )
	{
		T* input_vol = (T*)input_volumes[i];
		for ( int j = 0; j < nsize1; j++ )
		{
			*MATRIX_RELT( m1, j+1, i+1 ) = input_vol[vox1[j]];
		}
		for ( int j = 0; j < nsize2; j++ )
		{
			*MATRIX_RELT( m2, j+1, i+1 ) = input_vol[vox2[j]];
		}
	}
	VECTOR* mean1 = VectorAlloc( nvars, m1->type ); 
	VECTOR* mean2 = VectorAlloc( nvars, m2->type ); 
	MATRIX* cov1 = MatrixCovariance( m1, NULL, mean1 );
	if ( cov1 == NULL )
		return false;
	MATRIX* cov2 = MatrixCovariance( m2, NULL, mean2 );
	if ( cov2 == NULL )
		return false;
	MATRIX* cov = MatrixAdd( cov1, cov2, NULL );
	MATRIX* cov_inv = MatrixInverse( cov, NULL );
	if ( cov_inv == NULL )
		return false;
	MATRIX* mean_sub = MatrixSubtract( mean1, mean2, NULL );
	MATRIX* weight = MatrixMultiply( cov_inv, mean_sub, NULL );

	double* w = new double[nvars];
	double sum = 0;
	for ( int i = 0; i < nvars; i++ )
	{
		w[i] = fabs( *MATRIX_RELT( weight, i+1, 1 ) );
		sum += w[i];
		cout << w[i] << endl;
	}
	for ( int i = 0; i < nvars; i++ )
		w[i] /= sum;
	
	double tmp = 0;
	for ( int i = 0; i < vol_size; i++ )
	{
		tmp = 0;
		for ( int j = 0; j < nvars; j++ )
			tmp +=  ((T*)input_volumes[j])[i] * w[j];
		
		output_volume[i] = (T)tmp;
	}

	MatrixFree( &m1 );
	MatrixFree( &m2 );
	MatrixFree( &mean1 );
	MatrixFree( &mean2 );
	MatrixFree( &cov1 );
	MatrixFree( &cov2 );
	MatrixFree( &cov );
	MatrixFree( &cov_inv );
	MatrixFree( &mean_sub );
	MatrixFree( &weight );
	delete[] w;
	
	return true;
}


bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1, int* vox2, int nsize2,
									std::vector<void*> input_volumes, float* output_volume, int vol_size )
{
	return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2, input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1, int* vox2, int nsize2,
									  std::vector<void*> input_volumes, double* output_volume, int vol_size )
{
	return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2, input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1, int* vox2, int nsize2,
									  std::vector<void*> input_volumes, int* output_volume, int vol_size )
{
	return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2, input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1, int* vox2, int nsize2,
									  std::vector<void*> input_volumes, short* output_volume, int vol_size )
{
	return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2, input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1, int* vox2, int nsize2,
									  std::vector<void*> input_volumes, unsigned char* output_volume, int vol_size )
{
	return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2, input_volumes, output_volume, vol_size );
}

bool MyUtils::CalculateOptimalVolume( int* vox1, int nsize1, int* vox2, int nsize2,
									  std::vector<void*> input_volumes, long* output_volume, int vol_size )
{
	return CalculateOptimalVolume_t( vox1, nsize1, vox2, nsize2, input_volumes, output_volume, vol_size );
}

