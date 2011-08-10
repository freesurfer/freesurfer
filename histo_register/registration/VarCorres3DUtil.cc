#include "registration/VarCorres3DUtil.h"
#include <sbl/core/PathConfig.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
#include <sbl/image/MotionFieldUtil.h>
namespace hb {


//-------------------------------------------
// MISC UTILITY FUNCTIONS
//-------------------------------------------


/// returns vector of scales for multi-res optimization
VectorF varScaleSequence( float scaleFactor, float minScale ) {
	VectorF scales( 1 );
	scales[ 0 ] = 1;
	do {
		scales.append( scales.endValue() * scaleFactor );
	} while (scales.endValue() * scaleFactor > minScale);
	return scales;
}


/// build mapping of flow components to system variables;
/// each pixel in mask corresponds to 3 system variables (du, dv, and dw);
/// returns count of pixels in mask
int buildIndexMaps( int startIndex, const ImageGrayU &mask, ImageGrayI &uIndex, ImageGrayI &vIndex, ImageGrayI &wIndex ) {
	int width = mask.width(), height = mask.height();

	// fill in variable indices
	int index = startIndex;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask.data( x, y )) {
				uIndex.data( x, y ) = index++;
				vIndex.data( x, y ) = index++;
				wIndex.data( x, y ) = index++;
			}
		}
	}

	// return count of pixels in index maps
	return (index - startIndex) / 3;
}


//-------------------------------------------
// DISCRETE DERIVATIVES
//-------------------------------------------


/// discrete derivative in x direction
float dx( const ImageGrayF &img, int x, int y ) {
	return img.data( x, y ) - img.data( x - 1, y );
}


/// discrete derivative in y direction
float dy( const ImageGrayF &img, int x, int y ) {
	return img.data( x, y ) - img.data( x, y - 1 );
}


/// discrete derivative in z direction
float dz( const Array<ImageGrayF> &imgSeq, int x, int y, int z ) {
	return imgSeq[ z ].data( x, y ) - imgSeq[ z - 1 ].data( x, y );
}


/// apply discrete derivative in x direction to image
aptr<ImageGrayF> dx( const ImageGrayF &img ) {
	int width = img.width(), height = img.height();
	aptr<ImageGrayF> result( new ImageGrayF( width, height ) );
	result->clear( 0 );
	for (int y = 0; y < height; y++) {
		for (int x = 1; x < width; x++) {
			result->data( x, y ) = img.data( x, y ) - img.data( x - 1, y );
		}
	}
	return result;
}


/// apply discrete derivative in y direction to image
aptr<ImageGrayF> dy( const ImageGrayF &img ) {
	int width = img.width(), height = img.height();
	aptr<ImageGrayF> result( new ImageGrayF( width, height ) );
	result->clear( 0 );
	for (int y = 1; y < height; y++) {
		for (int x = 0; x < width; x++) {
			result->data( x, y ) = img.data( x, y ) - img.data( x, y - 1 );
		}
	}
	return result;
}


/// compute z-axis discrete derivative of a volume
aptr<ImageGrayF> dz( const Array<ImageGrayF> &seq, int i ) {
	int width = seq[ 0 ].width(), height = seq[ 0 ].height();
	aptr<ImageGrayF> result( new ImageGrayF( width, height ) );
	if (i == 0) {
		result->clear( 0 );
	} else {
		const ImageGrayF &img1 = seq[ i - 1 ];
		const ImageGrayF &img2 = seq[ i ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				result->data( x, y ) = img2.data( x, y ) - img1.data( x, y );
			}
		}
	}
	return result;
}


/// apply x discrete derivative to each channel
Array<ImageGrayF> dx( const Array<ImageGrayF> &chan ) {
	Array<ImageGrayF> chanDx;
	for (int i = 0; i < chan.count(); i++) 
		chanDx.append( dx( chan[ i ] ).release() );
	return chanDx;
}


/// apply y discrete derivative to each channel
Array<ImageGrayF> dy( const Array<ImageGrayF> &chan ) {
	Array<ImageGrayF> chanDy;
	for (int i = 0; i < chan.count(); i++) 
		chanDy.append( dy( chan[ i ] ).release() );
	return chanDy;
}


/// apply z discrete derivative to each channel
Array<ImageGrayF> dz( const ImageSetSeq &seq, int i ) {
	Array<ImageGrayF> chanDz;
	int width = seq[ 0 ][ 0 ].width(), height = seq[ 0 ][ 0 ].height();
	for (int k = 0; k < seq[ 0 ].count(); k++) {
		ImageGrayF *img = new ImageGrayF( width, height );
		if (i == 0) {
			img->clear( 0 );
		} else {
			const ImageGrayF &img1 = seq[ i - 1 ][ k ];
			const ImageGrayF &img2 = seq[ i ][ k ];
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					img->data( x, y ) = img2.data( x, y ) - img1.data( x, y );
				}
			}
		}
		chanDz.append( img );
	}
	return chanDz;
}


/// compute discrete derivative of a single slice of a correspondence volume
aptr<ImageGrayF> dz( const Array<CorresField3D> &cfSeq, int i, int axis ) {
	int width = cfSeq[ 0 ].width(), height = cfSeq[ 0 ].height();
	aptr<ImageGrayF> cfDz( new ImageGrayF( width, height ) );
	if (i == 0) {
		cfDz->clear( 0 );
	} else {
		const ImageGrayF *img1 = NULL; // seq[ i - 1 ][ k ];
		const ImageGrayF *img2 = NULL; // seq[ i ][ k ];
		if (axis == 0) {
			img1 = &cfSeq[ i - 1 ].u();
			img2 = &cfSeq[ i ].u();
		} else if (axis == 1) {
			img1 = &cfSeq[ i - 1 ].v();
			img2 = &cfSeq[ i ].v();
		} else if (axis == 2) {
			img1 = &cfSeq[ i - 1 ].w();
			img2 = &cfSeq[ i ].w();
		} else {
			fatalError( "invalid axis" );
		}
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				cfDz->data( x, y ) = img2->data( x, y ) - img1->data( x, y );
			}
		}
	}
	return cfDz;
}


//-------------------------------------------
// VARIATIONAL OBJECTIVE FUNCTION
//-------------------------------------------


/// this is the robust norm used in Brox et al. (denoted as Psi)
float psi( float diffSqd, bool robust ) {
	return robust ? sqrtf( diffSqd + 1e-6f ) : diffSqd;
}


/// derivative of robust norm w.r.t. squared diff
float psiDeriv( float diffSqd, bool robust ) {
	return robust ? 0.5f / sqrtf( diffSqd + 1e-6f ) : 1.0f;
}


//-------------------------------------------
// CONFIG WITH PARAMETERS FOR VAR MOTION
//-------------------------------------------


/// returns config object with parameters for variational correspondence algorithm
aptr<Config> _varCorres3DConfig() {
	String fileName = dataPath() + "varCorres3D.conf";
	aptr<Config> conf( new Config );
	conf->load( fileName );
	if (conf->entryCount() == 0) 
		warning( "failed to load varCorres3D config: %s", fileName.c_str() );
	return conf;
}


} // end namespace hb
