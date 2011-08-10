#include "Normalization.h"
#include <sbl/math/MathUtil.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/image/ImageUtil.h>
namespace hb {


/// compute a entropy image in which each pixel is set according to the amount of local entropy;
/// see: C. Wachinger and N. Navab, "Structural image representation for image registration," CVPRW 2010
aptr<ImageGrayU> entropyImage( const ImageGrayU &image, const ImageGrayU &mask, int patchRadius, float patchSigma ) {
	int width = image.width(), height = image.height();
	ImageGrayF entImage( width, height );
	entImage.clear( -1 );
	float patchGaussFactor = gaussFactor( patchSigma );
	float minEntropy = 1e8; 

	// loop over image computing entropy for each pixel
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask.data( x, y )) {

				// compute neighborhood bounds
				int xMin = x - patchRadius; 
				int xMax = x + patchRadius;
				int yMin = y - patchRadius;
				int yMax = y + patchRadius;
				if (xMin < 0) xMin = 0;
				if (xMax > width - 1) xMax = width - 1;
				if (yMin < 0) yMin = 0;
				if (yMax > height - 1) yMax = height - 1;

				// compute histogram of pixel values in neighborhood
				// fix(later): normalize so that pixel values in neighborhood span histogram?
				VectorF histogram( 64 );
				histogram.clear( 0 );
				int count = 0, outsideCount = 0;
				for (int yp = yMin; yp <= yMax; yp++) {
					for (int xp = xMin; xp <= xMax; xp++) {
						if (mask.data( xp, yp )) {
							int index = image.data( xp, yp ) / 4;
							int dx = xp - x;
							int dy = yp - y;
							float distSqd = (float) (dx * dx + dy * dy);
							float weight = gauss( distSqd, patchGaussFactor ); // note: might be faster to store in a table
							histogram[ index ] += weight;
							if (index - 1 >= 0)
								histogram[ index - 1 ] += weight * 0.5f;
							if (index + 1 < 64)
								histogram[ index + 1 ] += weight * 0.5f;
							count++;
						} else {
							outsideCount++;
						}
					}
				}

				// if have enough values
				if (count) {
	
					// normalize the histogram
					float factor = 1.0f / histogram.sum();
					multiply( histogram, factor, histogram );

					// compute entropy from histogram
					float e = (float) entropy( toDouble( histogram ) );
					entImage.data( x, y ) = e;
					if (e < minEntropy)
						minEntropy = e;
				}
			}
		}
	}

	// clear unititialzed values to the minimum value
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (entImage.data( x, y ) == -1)
				entImage.data( x, y ) = minEntropy;
		}
	}

	// display stats
	float min = 0, mean = 0, max = 0;
	imageStats( entImage, min, mean, max );
	disp( 1, "entropy: %f / %f / %f", min, mean, max );

	// transform to [0, 255] range
	aptr<ImageGrayU> entImageGray = toUChar( entImage );
	return entImageGray;
}


/// compute a normalized image by combining a blurred mask with an entropy image
aptr<ImageGrayU> normalize( const ImageGrayU &image, const ImageGrayU &mask, const ImageGrayU &maskBlurred ) {

	// parameters
	float beta = 0.5f;
	int patchRadius = 12;
	float patchSigma = 6.0f;
	float blendSigma = 10.0f;
	int blendPreBlurSize = 45;

	// compute blurred mask and entropy image
	aptr<ImageGrayU> entropy = entropyImage( image, mask, patchRadius, patchSigma );

	// compute blending mask
	aptr<ImageGrayU> maskBlend = blurBoxAndThreshold( mask, blendPreBlurSize, 250 );
	maskBlend = blurGauss( *maskBlend, blendSigma );

	// allocate normalized image
	int width = image.width(), height = image.height();
	aptr<ImageGrayU> norm( new ImageGrayU( width, height ) );

	// combine blurred mask and entropy image
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			float e = 1.0f - (float) entropy->data( x, y ) / 255.0f;
			float m = (float) maskBlurred.data( x, y ) / 255.0f;
			float b = beta * maskBlend->data( x, y ) / 255.0f;
			float v = b * e + (1.0f - b) * m;
			norm->data( x, y ) = bound( round( v * 255.0f ), 0, 255 );
		}
	}
	return norm;
}


} // end namespace hb
