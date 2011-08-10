#include "prep/MPrep.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
using namespace sbl;
namespace hb {


// a multi-channel image
typedef Array<ImageGrayU> MultiImage;


/// compute the sum-abs-difference between a pair of image patches
int patchDiff( MultiImage &image1, MultiImage &image2, int x1, int y1, int x2, int y2, int radius ) {
	int sumDiff = 0;
	for (int dy = -radius; dy <= radius; dy++) {
		for (int dx = -radius; dx <= radius; dx++) {
			for (int j = 0; j < image1.count(); j++) {
				int diff = image1[ j ].data( x1 + dx, y1 + dy ) - image2[ j ].data( x2 + dx, y2 + dy );
				if (diff < 0)
					diff = -diff;
				sumDiff += diff;
			}
		}
	}
	return sumDiff;
}


/// given a set of source and destination training pairs, predict a new destination image for the given source image
aptr<ImageGrayU> predict( Array<MultiImage> &srcTrainSet, Array<ImageGrayU> &destTrainSet, MultiImage &srcPredImage ) {
	int width = srcPredImage[ 0 ].width();
	int height = srcPredImage[ 0 ].height();
	aptr<ImageGrayU> destPredImage( new ImageGrayU( width, height ) );
	int radius = 2;
	
	// loop over pred pixels
	for (int yp = radius; yp < height - radius; yp++) {
		for (int xp = radius; xp < width - radius; xp++) {

			// will hold best match
			int bestDiff = -1;
			int bestVal = 0;

			// find best match for this pixel among all training images
			for (int i = 0; i < srcTrainSet.count(); i++) {
				MultiImage srcTrainImage = srcTrainSet[ i ];
				assertAlways( srcTrainImage.count() == srcPredImage.count() );
				assertAlways( srcTrainImage[ 0 ].width() == width && srcTrainImage[ 0 ].height() == height );
				for (int yt = radius; yt < height - radius; yt++) {
					for (int xt = radius; xt < width - radius; xt++) {
						int diff = patchDiff( srcTrainImage, srcPredImage, xt, yt, xp, yp, radius );
						if (diff < bestDiff || bestDiff == -1) {
							bestDiff = diff;
							bestVal = destTrainSet[ i ].data( xt, yt );
						}
					}
				}
			}

			// store best match's pixel value
			destPredImage->data( xp, yp ) = bestVal;
		}
	}
	return destPredImage;
}


/// given a set of source and destination training pairs, predict a new destination image for each test source image
void predictImages( Config &conf ) {

	// load training images

	// load test images

	// perform prediction

	// if have destination images for test images, compute prediction error

}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initPrediction() {
	registerCommand( "pred", predictImages );
}


} // end namespace hb
