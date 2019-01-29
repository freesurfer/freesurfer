#include "prep/BlockPrep.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/TensorUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/image/ImageRegister.h>
#include <sbl/image/Video.h>
#include <sbl/other/Plot.h>
#include <vector>
#include <string> 
#include <iostream>
#ifdef HAVE_OPENMP
  #include <omp.h>
#endif

using namespace sbl;
namespace hb {


/// takes a subset of the blockface images (to handle the case that multiple images were captured for each slice)
void selectBlockFaceImageSubset( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/extra" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/raw" ) );
	int imagesPerSlice = conf.readInt( "imagesPerSlice", 6 );
	int offset = conf.readInt( "offset", 5 );

	// create output dir
	createDir( outputPath );

	// loop over images
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );
	if (fileList.count() == 0)
    fileList = dirFileList( inputPath, "", ".jpg" ); 
	if (fileList.count() == 0)
    fileList = dirFileList( inputPath, "", ".PNG" ); 
	if (fileList.count() == 0)
    fileList = dirFileList( inputPath, "", ".png" ); 

	for (int inputIndex = offset; inputIndex < fileList.count(); inputIndex += imagesPerSlice) {
		disp( 1, "index: %d, file: %s", inputIndex, fileList[ inputIndex ].c_str() );

		// move file to output path
		moveFile( inputPath + "/" + fileList[ inputIndex ], outputPath + "/" + fileList[ inputIndex ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}


/// compute a rough mask of a block-face image
aptr<ImageGrayU> roughMask( const ImageColorU &image, int pixelstep = 1) {

	// blur to remove outliers
	aptr<ImageColorU> blurImage = blurBox( image, 3 );
	int width = image.width(), height = image.height();

/*	// compute mask of candidate pixels
	int candCount = 0;
	aptr<ImageGrayU> mask( new ImageGrayU( width, height ) );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int r = blurImage->r( x, y );
			int b = blurImage->b( x, y );
			if (r - b > 30) {
				mask->data( x, y ) = 255;
				candCount++;
			} else {
				mask->data( x, y ) = 0;
			}
		}
	}*/

  int histogramSize = 50;
	std::vector < int >  HistogramRaw(histogramSize,0);
	int histoCount = 0;

	// first run to compute min max of  red-blue
  int minrmb = 257;
  int maxrmb = -257;
	for (int y = 0; y < height; y+=pixelstep) {
		for (int x = 0; x < width; x+=pixelstep) {
			int rmb = blurImage->r( x, y ) - blurImage->b( x, y );
      if (rmb < minrmb) minrmb = rmb;
      if (rmb > maxrmb) maxrmb = rmb;
		}
	}
  //disp( 1, "rhisto: minrmb: %d, maxrmb: %d", minrmb, maxrmb );
  //disp( 1, "histosize: %d", histogramSize );
  
	// second run to compute histogram
	for (int y = 0; y < height; y+=pixelstep) {
		for (int x = 0; x < width; x+=pixelstep) {
			double rmb = blurImage->r( x, y ) - blurImage->b( x, y );
      rmb = double(rmb-minrmb)/(maxrmb-minrmb);  // scale to  0..1
      int rmbh = (int) (rmb * (histogramSize - 1));
      
      if (rmbh < 0 || rmbh >= histogramSize)
      {
        disp( 1, "Error out of bounds histo: %f %d",rmb, rmbh );
        exit(1);
      } 
      
      HistogramRaw[ rmbh ]++;
      histoCount++;
		}
	}

  // compute split in histogram
  int s = HistogramRaw.size();
  int max1 = 0; double max1val = HistogramRaw.front();
  int max2 = s -1 ; double max2val = HistogramRaw.back();
  for (int i = 1 ; i< histogramSize/2 ; i++)
  {
    if (HistogramRaw[i] > max1val) 
    {
      max1val = HistogramRaw[i];
      max1 = i;
    }
    if (HistogramRaw[s-1-i] > max2val)
    {
      max2val = HistogramRaw[s-1-i];
      max2 = s-1-i;
    }
  }
  int min = max1; double minval = HistogramRaw[max1];
  for (int i = max1+1; i<max2;i++)
    if (HistogramRaw[i] < minval)
    {
      minval = HistogramRaw[i];
      min = i;
    }
    

	// compute mask of candidate pixels based on histogram split at min
	int candCount = 0;
	aptr<ImageGrayU> mask( new ImageGrayU( width, height ) );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double rmb = blurImage->r( x, y ) - blurImage->b( x, y );
      rmb = double(rmb-minrmb)/(maxrmb-minrmb);  // scale to  0..1
      int rmbh = (int) (rmb * (histogramSize - 1));
			if (rmbh > min) {
				mask->data( x, y ) = 255;
				candCount++;
			} else {
				mask->data( x, y ) = 0;
			}
		}
	}


	// clean up the mask
	mask = blurBoxAndThreshold( *mask, 5, 128 );

	// mark any component near the center
	int xMin = width / 2 - width / 8;
	int xMax = width / 2 + width / 8;
	int yMin = height / 2 - height / 8;
	int yMax = height / 2 + height / 8;
	for (int y = yMin; y <= yMax; y++) {
		for (int x = xMin; x <= xMax; x++) {
			if (mask->data( x, y ) == 255) {
				floodFill( *mask, 255, 255, 100, x, y );
			}
		}
	}

	// only keep marked components
	int finalCount = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask->data( x, y ) == 100) {
				mask->data( x, y ) = 255;
				finalCount++;
			} else if (mask->data( x, y )) {
				mask->data( x, y ) = 0;
			}
		}
	}
	disp( 1, "Rough mask: candCount: %d, finalCount: %d", candCount, finalCount );
      
	return mask;
}

double normalizeBackground( ImageColorU& image, const ImageGrayU& mask)
{
	int width = image.width(), height = image.height();
  long int sumr = 0;
  long int sumg = 0;
  long int sumb = 0;
  long int bcount = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask.data( x, y ) == 0) {
				sumr += image.r( x, y );
        sumg += image.g( x, y );
        sumb += image.b( x, y );
        bcount++;
			}
    }
  }

  double scaler = (200.0 * bcount) / sumr;
  double scaleg = (200.0 * bcount) / sumg;
  double scaleb = (200.0 * bcount) / sumb;
  disp(1,"mean background rgb %f %f %f",((double)sumr)/bcount,((double)sumg)/bcount,((double)sumb)/bcount);
  //disp(1,"scale rgb %f %f %f",scaler,scaleg,scaleb);
  double rval,gval,bval;
  double mval = 0.0;
  long int fcount = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
		  rval= image.r(x,y) * scaler;
      if (rval > 255) rval = 255;
      image.setR(x,y,(unsigned char)rval);
		  gval= image.g(x,y) * scaleg;
      if (gval > 255) gval = 255;
      image.setG(x,y,(unsigned char)gval);
		  bval= image.b(x,y) * scaleb;
      if (bval > 255) bval = 255;
      image.setB(x,y,(unsigned char)bval);
      if (mask.data( x, y) ==255){
        mval += (rval+gval+bval)/3.0;
        fcount++;
      }
    }
  }
  mval /= fcount;
  disp(1,"mean foreground intensity %f",mval);
  return mval;

}


/// estimate a transformation that registers the images before/after the microtome re-adjustment
aptr<ImageTransform> computeSplitTransform( const String &fileName1, const String &fileName2, const String &visPath ) {
	bool useInternal = true;

	// load boundary images
	aptr<ImageColorU> image1 = load<ImageColorU>( fileName1 );
	aptr<ImageColorU> image2 = load<ImageColorU>( fileName2 );

	// compute masks
	aptr<ImageGrayU> mask1 = roughMask( *image1 );
	aptr<ImageGrayU> mask2 = roughMask( *image2 );
	mask1 = blurBox( *mask1, 15 );
	mask2 = blurBox( *mask2, 15 );

	// perform initial registration using masks
	int paramCount = 6;
	int step = 2;
	int xBorder = 10;
	int yBorder = 10;
	float offsetBound = 200;
	aptr<ImageTransform> transform = registerUsingImageTransform( *mask1, *mask2, paramCount, step, xBorder, yBorder, offsetBound );
	disp( 1, "init transform:" );
	transform->display( 2 );

	// perform final registration using images inside masks
	if (useInternal) {
		aptr<ImageGrayU> regMask1 = threshold( *mask1, 250, false );
		aptr<ImageGrayU> regMask2 = threshold( *mask2, 250, false );
		aptr<ImageGrayU> regImage1 = blurBox( *toGray( *image1 ), 3 );
		aptr<ImageGrayU> regImage2 = blurBox( *toGray( *image2 ), 3 );
		transform = registerUsingImageTransform( *regImage1, *regImage2, paramCount, step, xBorder, yBorder, offsetBound, transform.get(), regMask1.get(), regMask2.get() );
		disp( 1, "final transform:" );
		transform->display( 2 );
	}

	// save diagnostic images
	int width = image1->width(), height = image1->height();
	aptr<ImageColorU> mapped = transform->mapForward( *image1, width, height, 200 );
	saveImage( *mapped, visPath + "/mapped.png" );
	saveImage( *image2, visPath + "/image.png" );
	saveImage( *mask1, visPath + "/mask1.png" );
	saveImage( *mask2, visPath + "/mask2.png" );
	return transform;
}


/// compute the approximate x-y-plane diameter of the MRI volume in pixels (assumes images are already segmented)
int computeMriSize( const String &mriPath ) {
	Array<String> fileList = dirFileList( mriPath, "", ".png" );
	int xMin = 100000, xMax = 0, yMin = 100000, yMax = 0;
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
		aptr<ImageGrayU> image = load<ImageGrayU>( mriPath + "/" + fileList[ inputIndex ] );
		int width = image->width(), height = image->height();
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (image->data( x, y )) {
					if (x < xMin) 
						xMin = x;
					if (x > xMax)
						xMax = x;
					if (y < yMin)
						yMin = y;
					if (y > yMax)
						yMax = y;
				}
			}
		}
	}
	int xSize = xMax - xMin + 1;
	int ySize = yMax - yMin + 1;
	return (xSize + ySize) / 2;
}


/// crops, aligns, and resizes a set of block-face imges
void cropBlockFaceImages( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/raw" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/crop" ) );
	String outSegPath = addDataPath( conf.readString( "outputPath", "blockface/roughseg" ) );
	String mriPath = addDataPath( conf.readString( "mriPath", "mri/seg" ) );
  int splitIndex = conf.readInt( "splitIndex", -2 );
  bool findsplit = (splitIndex == -2 ); //automatically estimate split
	int visScale = 8;
	int pad = 100; // note: the pad will end up being larger than this because of the split transform

	// create output dir
	createDir( outputPath );
	createDir( outSegPath );

	// get input file list
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );
	if (fileList.count() == 0)
    fileList = dirFileList( inputPath, "", ".jpg" ); 
	if (fileList.count() == 0)
    fileList = dirFileList( inputPath, "", ".PNG" ); 
	if (fileList.count() == 0)
    fileList = dirFileList( inputPath, "", ".png" ); 
	if (fileList.count() == 0) { 
		warning( "no input files at %s", inputPath.c_str() );
		return;
	}

	// prepare visualization path
	String visPath = outputPath + "/vis";
	createDir( visPath );
  //createDir( outSegPath + "/vis");

	// load first image to get dimensions
	aptr<ImageColorU> firstImage = load<ImageColorU>( inputPath + "/" + fileList[ 0 ] );
	int width = firstImage->width(), height = firstImage->height();
	firstImage.release();

	// open visualization video
	int visWidth = ((width / visScale) / 8) * 8; // make sure video size divisible by 8
	int visHeight = ((height / visScale) / 8) * 8;
	OutputVideo outputVideo( visPath + "/vis.avi", visWidth, visHeight );

	// compute MR bounds
	int mriSize = computeMriSize( mriPath );

	// first pass: compute bounds
	int xMin = 100000, xMax = 0, yMin = 100000, yMax = 0;
	VectorD countVect(fileList.count());
  bool keepgoing = true;
  bool doreturn = false;
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
	
    if (! keepgoing) continue; // can't break in parallel loop
		// load the input image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
		assertAlways( image->width() == width && image->height() == height );

		// compute rough mask of tissue area
		aptr<ImageGrayU> mask = roughMask( *image, 1 );
    saveImage( *mask,  outSegPath + "/" +  fileList[ inputIndex ].leftOfLast( '.' ) + ".png" );
    
		// find bounds of mask area
		int insideCount = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (mask->data( x, y ) == 255) {
					insideCount++;
					if (x < xMin) 
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
{	  				if (x < xMin) xMin = x; }
					if (x > xMax)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
{						if (x > xMax) xMax = x; }
					if (y < yMin)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
{						if (y < yMin) yMin = y; }
					if (y > yMax)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
{						if (y > yMax) yMax = y; }
					image->setRGB( x, y, 255, 0, 0 );
				}
			}
		}
		    
		// store for counts
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
		disp( 1, "inputIndex: %d, file: %s, insideCount: %d", inputIndex, fileList[ inputIndex ].c_str(), insideCount );
		//disp( 1, "inputIndex: %d, file: %s, insideCount: %d, diff: %d", inputIndex, fileList[ inputIndex ].c_str(), insideCount, diff );
		countVect[inputIndex] =  (double) insideCount ;

		// check for bad slice
		if (insideCount < 100000 ) {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
{
			warning( "bad slice (low inside count %i): %s; stopping command",insideCount, fileList[ inputIndex ].c_str() );
      keepgoing = false;
      doreturn = true;
//			return;
}
		}
    if (! keepgoing) continue; // can't break in parallel loop

		// create visualization image
//		aptr<ImageColorU> visImage = resize( *image, visWidth, visHeight, true );
//		outputVideo.append( *visImage );
////		saveImage( *visImage, visPath + "/" + fileList[ inputIndex ] );

		// check for user cancel
#ifdef HAVE_OPENMP
    if(omp_get_thread_num() == 0)
#endif
    {
		  if (checkCommandEvents())
      {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
          keepgoing = false;
//			  break;
      }
    }
	}
  if (doreturn) return;
	if (checkCommandEvents())
		return;
	disp( 1, "orig xMin: %d, xMax: %d, yMin: %d, yMax: %d", xMin, xMax, yMin, yMax );
	xMin = bound( xMin - pad, 0, width - 1 );
	xMax = bound( xMax + pad, 0, width - 1 );
	yMin = bound( yMin - pad, 0, height - 1 );
	yMax = bound( yMax + pad, 0, height - 1 );
	disp( 1, "crop xMin: %d, xMax: %d, yMin: %d, yMax: %d", xMin, xMax, yMin, yMax );

	// save plot of counts
	Plot plot;
	plot.add( countVect );
	plot.save( visPath + "/counts.svg" );

  // find split and transform
	int maxDiff = 0;
  aptr<ImageTransform> transform;
	for (int inputIndex = 1; inputIndex < fileList.count(); inputIndex++)
  {
    // check for split
    int diff = iAbs( (int) countVect[inputIndex] - (int) countVect[inputIndex-1] );
    if ( diff > 1000000) {
			  warning( "bad slice: %s; stopping command", fileList[ inputIndex ].c_str() );
			  return;
    }
    if (diff > maxDiff)
    {
  		maxDiff = diff;
      if (findsplit) splitIndex = inputIndex - 1;
    }
  }

  
  if (splitIndex > 0)
  {
	  // compute transformation
	  String fileName1 = fileList[ splitIndex ];
	  String fileName2 = fileList[ splitIndex + 1 ];
    int diff = iAbs( (int) countVect[splitIndex] - (int) countVect[splitIndex+1]);
	  disp( 1, "split index: %d, diff: %d, max diff: %d, file 1: %s, file 2: %s", splitIndex, diff, maxDiff, fileName1.c_str(), fileName2.c_str() );
	  disp( 1, "Computing Transform (microtome adjustment) ...");
	  transform = computeSplitTransform( inputPath + "/" + fileName1, inputPath + "/" + fileName2, visPath );
  }

	// compute final scale factor
	int xSize = xMax - xMin + 1;
	int ySize = yMax - yMin + 1;
	int blockSize = (xSize + ySize) / 2;
	int outputWidth = xSize * mriSize / blockSize;
	int outputHeight = ySize * mriSize / blockSize;
	outputWidth -= outputWidth & 7;
	outputHeight -= outputHeight & 7;
	disp( 1, "block size: %d, mri size: %d, output width: %d, output height: %d", blockSize, mriSize, outputWidth, outputHeight );

	// second pass: crop images
  keepgoing = true;
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
	
    if (! keepgoing) continue;
  
		// load the image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
		disp( 1, "file: %s", fileList[ inputIndex ].c_str() );
		aptr<ImageGrayU> mask = load<ImageGrayU>( outSegPath + "/" + fileList[ inputIndex ].leftOfLast( '.' ) + ".png" );
    normalizeBackground(*image,*mask);

		// if before split, apply transform
		if (inputIndex <= splitIndex) 
			image = transform->mapForward( *image, width, height, 255 );

		// crop and shrink the image
		image = crop( *image, xMin, xMax, yMin, yMax );
		image = resize( *image, outputWidth, outputHeight, true );

		// save to output path
		saveImage( *image, outputPath + "/" + fileList[ inputIndex ].leftOfLast( '.' ) + ".png" );

		// check for user cancel
#ifdef HAVE_OPENMP
    if(omp_get_thread_num() == 0)
#endif
    {
		  if (checkCommandEvents())
      {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
          keepgoing = false;
//			  break;
      }
    }
	}
/*	// second pass: crop masks
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
	
		// load the image
		aptr<ImageColorU> image = load<ImageColorU>( outSegPath + "/" + fileList[ inputIndex ].leftOfLast( '.' ) + ".png" );
		disp( 1, "Seg file: %s", fileList[ inputIndex ].c_str() );

		// if before split, apply transform
		if (inputIndex <= splitIndex) 
			image = transform->mapForward( *image, width, height, 255 );

		// crop and shrink the image
		image = crop( *image, xMin, xMax, yMin, yMax );
		image = resize( *image, outputWidth, outputHeight, true );

		// save to output path
		saveImage( *image, outSegPath + "/" + fileList[ inputIndex ].leftOfLast( '.' ) + ".png" );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}*/
}


// data type used for block-face histogram
typedef Tensor3F Histogram;


/// initialize a histogram with the specified number of bins in each dimension
void initHistogram( Histogram &histogram, int histogramSize ) {
	for (int i = 0; i < histogram.dimCount(); i++)
		histogram.setSize( i, histogramSize );
	histogram = 0;
}


/// compute a feature vector and transform to histogram indices
VectorI featureToIndices( int r, int g, int b, int histogramSize ) {
	VectorI featInd( 3 );
	featInd[ 0 ] = r * (histogramSize - 1) / 255;
	featInd[ 1 ] = g * (histogramSize - 1) / 255;
	featInd[ 2 ] = b * (histogramSize - 1) / 255;
	return featInd;
}


/// segment a set of block-face images
void segmentBlockFaceImages( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/crop" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/seg" ) );
	int histogramSize = conf.readInt( "histogramSize", 30 );
	int histogramBlurSize = conf.readInt( "histogramBlurSize", 5 );

	// create output dir
	createDir( outputPath );

	// if already run, stop here
//	if (dirFileList( outputPath, "", ".png" ).count()) {
//		disp( 1, "output path already has files; stopping" );
//		return;
//	}

	// get input file list
	Array<String> fileList = dirFileList( inputPath, "", ".png" );
	if (fileList.count() == 0) { 
		warning( "no input files at %s", inputPath.c_str() );
		return;
	}

	// prepare visualization path
	String visPath = outputPath + "/vis";
	createDir( visPath );

	// initialize histograms for foreground and background distributions
	Histogram fgHistogramRaw, bgHistogramRaw;
	initHistogram( fgHistogramRaw, histogramSize );
	initHistogram( bgHistogramRaw, histogramSize );
	int fgCount = 0, bgCount = 0;

	// we'll skip over some pixels and some images for the sake of efficiency
	int pixelStep = 3;
	int imageStep = 3;

	// first pass: build models
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex += imageStep) {

		// load image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
		int width = image->width(), height = image->height();

		// compute rough mask of tissue area
		aptr<ImageGrayU> mask = roughMask( *image );

		// blur the mask to obtain inside, outside, and uncertain areas
		mask = blurBox( *mask, 35 );

		// update histograms
		for (int y = 0; y < height; y += pixelStep) {
			for (int x = 0; x < width; x += pixelStep) {
				int r = image->r( x, y );
				int g = image->g( x, y );
				int b = image->b( x, y );
				VectorI featInd = featureToIndices( r, g, b, histogramSize );
				int m = mask->data( x, y );
				if (m == 0) {
					bgHistogramRaw.elem( featInd.dataPtr() )++;
					bgCount++;
				} else if (m == 255) {
					fgHistogramRaw.elem( featInd.dataPtr() )++;
					fgCount++;
				}
			}
		}
		status( "inputIndex: %d, fgCount: %d, bgCount: %d", inputIndex, fgCount, bgCount );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	status( "\n" );
	disp( 1, "fgCount: %d, bgCount: %d", fgCount, bgCount );

	// create blurred histograms
	Histogram fgHistogram, bgHistogram;
	initHistogram( fgHistogram, histogramSize );
	initHistogram( bgHistogram, histogramSize );
	blurBox( fgHistogramRaw, fgHistogram, histogramBlurSize );
	blurBox( bgHistogramRaw, bgHistogram, histogramBlurSize );

	// normalize the histograms
	float fgNorm = 1.0f / fgCount;
	float bgNorm = 1.0f / bgCount;
	fgHistogram *= fgNorm;
	bgHistogram *= bgNorm;

	// second pass: perform (soft) segmentation
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {

		// load image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
		int width = image->width(), height = image->height();

		// init mask
		aptr<ImageGrayU> mask( new ImageGrayU( width, height ) );

		// compute soft mask value using histograms
		int fgOutCount = 0, bgOutCount = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int r = image->r( x, y );
				int g = image->g( x, y );
				int b = image->b( x, y );
				VectorI featInd = featureToIndices( r, g, b, histogramSize );
				float fgVal = fgHistogram.elem( featInd.dataPtr() );
				float bgVal = bgHistogram.elem( featInd.dataPtr() );
				if (fgVal > 10 * bgVal) {
					mask->data( x, y ) = 255;
					fgOutCount++;
				} else {
					mask->data( x, y ) = 0;
					bgOutCount++;
				}
			}
		}
		status( "inputIndex: %d, fgCount: %d, bgCount: %d", inputIndex, fgOutCount, bgOutCount );

		// save the mask
		saveImage( *mask, outputPath + "/" + fileList[ inputIndex ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	status( "\n" );
}

/// segment a set of cropped block-face images (using only rough mask)
void segmentBlockFaceImagesSimple( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/crop" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/seg" ) );

	// create output dir
	createDir( outputPath );

	// get input file list
	Array<String> fileList = dirFileList( inputPath, "", ".png" );
	if (fileList.count() == 0) { 
		warning( "no input files at %s", inputPath.c_str() );
		return;
	}

	// prepare visualization path
	String visPath = outputPath + "/vis";
	createDir( visPath );

	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
    disp(1,"file name %s ",fileList[ inputIndex ].c_str());
  
		// load image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );

		// compute rough mask of tissue area
		aptr<ImageGrayU> mask = roughMask( *image);

		// save the mask
		saveImage( *mask, outputPath + "/" + fileList[ inputIndex ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
  
	status( "\n" );
}

/// perform all steps needed to prepare a set of block-face images for registration
void prepareBlockFaceImages( Config &conf ) {
	cropBlockFaceImages( conf );
	execCommand( "vcross blockface/crop blockface/crop/vis/cross", false );
	segmentBlockFaceImagesSimple( conf );
	execCommand( "vcross blockface/seg blockface/seg/vis/cross", false );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initBlockPrep() {
	registerCommand( "bsubset", selectBlockFaceImageSubset );
	registerCommand( "bcrop", cropBlockFaceImages ); 
	registerCommand( "bseg", segmentBlockFaceImages ); 
	registerCommand( "bsegsimple", segmentBlockFaceImagesSimple ); 
	registerCommand( "bprep", prepareBlockFaceImages ); 
}


} // end namespace hb
