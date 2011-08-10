#include "registration/VarCorres3D.h"
#include <sbl/core/Command.h> // for checkCommandEvents
#include <sbl/core/PathConfig.h> 
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageSeqUtil.h>
#include <sbl/image/ImageTransform.h>
#include <sbl/image/MotionFieldUtil.h>
#include <pvl/SparseSystem.h>
#include "registration/VarCorres3DUtil.h"
#include "registration/ImageSetSeq.h"
using namespace pvl;
namespace hb {


// dealing with edge cases:
// - the mask specifies which pixels are variables in the sparse system (to have flow estimated)
// - a flow vector that projects outside the image has no data penalty
// - the lower (y=0) and left (x=0) pixels have no data term and only a basic smoothness term


//-------------------------------------------
// DIAGNOSTIC UTILS
//-------------------------------------------


/// display statistics about the given image
void dispStats( const String &name, const ImageGrayF &img ) {
	int width = img.width(), height = img.height();
	double min = 1e10, max = -1e10, sum = 0;
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++) {
			float val = img.data( x, y );
			if (val < min)
				min = val;
			if (val > max)
				max = val;
			sum += val;
		}
	double mean = sum / (double) (width * height);
	disp( 5, "[%s]: min: %f, mean: %f, max: %f", name.c_str(), min, mean, max );
}


/// display statistics about the given image
void dispStats( const String &name, const ImageGrayFSeq &imgSeq ) {
	int length = imgSeq.count();
	int width = imgSeq[ 0 ].width(), height = imgSeq[ 0 ].height();
	double min = 1e10, max = -1e10, sum = 0;
	for (int z = 0; z < length; z++) {
		const ImageGrayF &img = imgSeq[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				float val = img.data( x, y );
				if (val < min)
					min = val;
				if (val > max)
					max = val;
				sum += val;
			}
		}
	}
	double mean = sum / (double) (width * height * length);
	disp( 5, "[%s]: min: %f, mean: %f, max: %f", name.c_str(), min, mean, max );
}


/// The VarCorresInfo class holds a collection of information used to update the correspondence vectors. 
class VarCorresInfo {
public:

	// basic constructor
	VarCorresInfo() { pixelCount = 0; }

	// number of pixels being optimized
	int pixelCount;

	// mask of pixels being optimized
	Array<ImageGrayU> mask;

	// index (if any) of each pixel within sparse system
	Array<ImageGrayI> uIndex;
	Array<ImageGrayI> vIndex;
	Array<ImageGrayI> wIndex;

	// derivatives of correspondence components
	Array<ImageGrayF> uDx; 
	Array<ImageGrayF> uDy;
	Array<ImageGrayF> uDz;
	Array<ImageGrayF> vDx;
	Array<ImageGrayF> vDy; 
	Array<ImageGrayF> vDz; 
	Array<ImageGrayF> wDx;
	Array<ImageGrayF> wDy; 
	Array<ImageGrayF> wDz; 

	// the current update to the correspondence vectors
	Array<ImageGrayF> du;
	Array<ImageGrayF> dv;
	Array<ImageGrayF> dw;

	// other optimizations terms
	Array<ImageGrayF> ix;
	Array<ImageGrayF> iy;
	Array<ImageGrayF> iz;
	Array<ImageGrayF> ic;
	Array<ImageGrayF> dPsiData;
	Array<ImageGrayF> dPsiSmoothness;

	/// check that all the data sequences have the correct dimensions
	void check( int width, int height, int depth );

	/// display stats about the data sequences
	void dispStats( int indent, int z );
};


// check the image sequence has the dimensions
template<typename ImageType> void checkImageSeqDimensions( const Array<ImageType> &seq, int width, int height, int depth ) {
	assertAlways( seq.count() == depth && seq[ 0 ].width() == width && seq[ 0 ].height() == height );
}


/// check that all the data sequences have the correct dimensions
void VarCorresInfo::check( int width, int height, int depth ) {
	checkImageSeqDimensions( mask, width, height, depth );
	checkImageSeqDimensions( uIndex, width, height, depth );
	checkImageSeqDimensions( vIndex, width, height, depth );
	checkImageSeqDimensions( wIndex, width, height, depth );
	checkImageSeqDimensions( uDx, width, height, depth );
	checkImageSeqDimensions( uDy, width, height, depth );
	checkImageSeqDimensions( uDz, width, height, depth );
	checkImageSeqDimensions( vDx, width, height, depth );
	checkImageSeqDimensions( vDy, width, height, depth );
	checkImageSeqDimensions( vDz, width, height, depth );
	checkImageSeqDimensions( wDx, width, height, depth );
	checkImageSeqDimensions( wDy, width, height, depth );
	checkImageSeqDimensions( wDz, width, height, depth );
	checkImageSeqDimensions( ix, width, height, depth );
	checkImageSeqDimensions( iy, width, height, depth );
	checkImageSeqDimensions( iz, width, height, depth );
	checkImageSeqDimensions( ic, width, height, depth );
	checkImageSeqDimensions( dPsiData, width, height, depth );
	checkImageSeqDimensions( dPsiSmoothness, width, height, depth );
	checkImageSeqDimensions( du, width, height, depth );
	checkImageSeqDimensions( dv, width, height, depth );
	checkImageSeqDimensions( dw, width, height, depth );
}


/// display stats about the data sequences 
void VarCorresInfo::dispStats( int indent, int z ) {
	hb::dispStats( "uDx", uDx[ z ] );
	hb::dispStats( "uDy", uDy[ z ] );
	hb::dispStats( "uDz", uDz[ z ] );
	hb::dispStats( "vDx", vDx[ z ] );
	hb::dispStats( "vDy", vDy[ z ] );
	hb::dispStats( "vDz", vDz[ z ] );
	hb::dispStats( "wDx", wDx[ z ] );
	hb::dispStats( "wDy", wDy[ z ] );
	hb::dispStats( "wDz", wDz[ z ] );
	hb::dispStats( "ix", ix[ z ] );
	hb::dispStats( "iy", iy[ z ] );
	hb::dispStats( "iz", iz[ z ] );
	hb::dispStats( "ic", ic[ z ] );
	hb::dispStats( "du", du[ z ] );
	hb::dispStats( "dv", dv[ z ] );
	hb::dispStats( "dw", dw[ z ] );
}


//-------------------------------------------
// VARIATIONAL OPTIMIZATION
//-------------------------------------------


/// solves for du, dv, dw (for all pixels in the mask)
void varSolveSystem( VarCorresInfo &vci, Config &varConf, bool fast ) {

	// check parameters
	assertAlways( vci.pixelCount );
	int depth = vci.mask.count();
	assertAlways( depth );
	int width = vci.mask[ 0 ].width(), height = vci.mask[ 0 ].height();
	vci.check( width, height, depth );

	// get parameters
	float smoothness = varConf.readFloat( "smoothness" ); // aka alpha (in Brox papers)
	int solverItersPerMegapixel = varConf.readInt( "solverItersPerMegapixel" );
	float sorFactor = varConf.readFloat( "sorFactor" );
	float maxMove = varConf.readFloat( "maxMove" );
	int adaptIterMask = varConf.readInt( "adaptIterMask" );
	float adaptUpdateThresh = varConf.readFloat( "adaptUpdateThresh" );
	float wSmooth = 1.0f;
//	bool zTerms = false;

	// add more iters if small
	int solverIters = sbl::round( (float) solverItersPerMegapixel * (float) vci.pixelCount / 1000000.0f);
	if (fast) {
		adaptUpdateThresh *= 10.0f;
		solverIters /= 2;
	}

	// allocate sparse linear system
	SparseSystem sparseSystem( vci.pixelCount * 3 );
	VectorF init( vci.pixelCount * 3 );
	init.clear( -7777.7777f );

	// loop over pixels, adding equation for each variable (three variables per pixel)
	for (int z = 0; z < depth; z++) {
		const ImageGrayU &mask = vci.mask[ z ];
		const ImageGrayF &du = vci.du[ z ];
		const ImageGrayF &dv = vci.dv[ z ];
		const ImageGrayF &dw = vci.dw[ z ];
		const ImageGrayF &dPsiSmoothness = vci.dPsiSmoothness[ z ];
		const ImageGrayF &dPsiData = vci.dPsiData[ z ];
		const ImageGrayI &uIndex = vci.uIndex[ z ];
		const ImageGrayI &vIndex = vci.vIndex[ z ];
		const ImageGrayI &wIndex = vci.wIndex[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (mask.data( x, y )) {

					// get indices of variables used in this pair of equation
					int uInd = uIndex.data( x, y );
					int vInd = vIndex.data( x, y );
					int wInd = wIndex.data( x, y );
					int uMxInd = x > 0 ? uIndex.data( x - 1, y ) : -1;
					int uPxInd = x < width - 1 ? uIndex.data( x + 1, y ) : -1;
					int uMyInd = y > 0 ? uIndex.data( x, y - 1 ) : -1;
					int uPyInd = y < height - 1 ? uIndex.data( x, y + 1 ) : -1;
					int uMzInd = z > 0 ? vci.uIndex[ z - 1 ].data( x, y ) : -1;
					int uPzInd = z < depth - 1 ? vci.uIndex[ z + 1 ].data( x, y ) : -1;

					int vMxInd = x > 0 ? vIndex.data( x - 1, y ) : -1;
					int vPxInd = x < width - 1 ? vIndex.data( x + 1, y ) : -1;
					int vMyInd = y > 0 ? vIndex.data( x, y - 1 ) : -1;
					int vPyInd = y < height - 1 ? vIndex.data( x, y + 1 ) : -1;
					int vMzInd = z > 0 ? vci.vIndex[ z - 1 ].data( x, y ) : -1;
					int vPzInd = z < depth - 1 ? vci.vIndex[ z + 1 ].data( x, y ) : -1;

					int wMxInd = x > 0 ? wIndex.data( x - 1, y ) : -1;
					int wPxInd = x < width - 1 ? wIndex.data( x + 1, y ) : -1;
					int wMyInd = y > 0 ? wIndex.data( x, y - 1 ) : -1;
					int wPyInd = y < height - 1 ? wIndex.data( x, y + 1 ) : -1;
					int wMzInd = z > 0 ? vci.wIndex[ z - 1 ].data( x, y ) : -1;
					int wPzInd = z < depth - 1 ? vci.wIndex[ z + 1 ].data( x, y ) : -1;

					// check indices
					if (uInd < 0 || uInd >= init.length() || vInd < 0 || vInd >= init.length() || wInd < 0 || wInd >= init.length()) 
						fatalError( "SolveInnerVarSystem: sanity check failed" );
		
					// store current value of du, dv, dw as initial values for iterative solver
					init[ uInd ] = du.data( x, y );
					init[ vInd ] = dv.data( x, y );
					init[ wInd ] = dw.data( x, y );

					// these coefs are for the first equation (derivative w.r.t. du)
					float duCoef1 = 0, dvCoef1 = 0, dwCoef1 = 0, b1 = 0;

					// these coefs are for the second equation (derivative w.r.t. dv)
					float duCoef2 = 0, dvCoef2 = 0, dwCoef2 = 0, b2 = 0;

					// these coefs are for the third equation (derivative w.r.t. dw)
					float duCoef3 = 0, dvCoef3 = 0, dwCoef3 = 0, b3 = 0;

					// prep for data coefs
					float ix = vci.ix[ z ].data( x, y );
					float iy = vci.iy[ z ].data( x, y );
					float iz = vci.iz[ z ].data( x, y );
					float ic = vci.ic[ z ].data( x, y );
					float dataFactor = dPsiData.data( x, y );

					// compute data coefs
					duCoef1 += dataFactor * ix * ix;
					dvCoef1 += dataFactor * ix * iy;
					dwCoef1 += dataFactor * ix * iz;
					b1 += -dataFactor * ix * ic; // negate because subtract from both sides
					duCoef2 += dataFactor * iy * ix;
					dvCoef2 += dataFactor * iy * iy;
					dwCoef2 += dataFactor * iy * iz;
					b2 += -dataFactor * iy * ic; // negate because subtract from both sides
					duCoef3 += dataFactor * iz * ix;
					dvCoef3 += dataFactor * iz * iy;
					dwCoef3 += dataFactor * iz * iz;
					b3 += -dataFactor * iz * ic; // negate because subtract from both sides

					// prep for smoothness term
					float alpha = smoothness; // + localSmoothness->data( x, y );
					float smoothFactor = alpha * dPsiSmoothness.data( x, y );
					float smoothFactorPx = x < width - 1 ? alpha * dPsiSmoothness.data( x + 1, y ) : 0;
					float smoothFactorPy = y < height - 1 ? alpha * dPsiSmoothness.data( x, y + 1 ) : 0;
					float smoothFactorPz = z < depth - 1 ? alpha * vci.dPsiSmoothness[ z + 1 ].data( x, y ) : 0;

					// smoothness: u eqn coefs
					duCoef1 += 3.0f * smoothFactor;
					duCoef1 += smoothFactorPx;
					duCoef1 += smoothFactorPy;
					duCoef1 += smoothFactorPz;
					float duMxCoef = -smoothFactor;
					float duMyCoef = -smoothFactor;
					float duMzCoef = -smoothFactor;
					float duPxCoef = -smoothFactorPx;
					float duPyCoef = -smoothFactorPy;
					float duPzCoef = -smoothFactorPz;
					b1 -= (vci.uDx[ z ].data( x, y ) + vci.uDy[ z ].data( x, y ) + vci.uDz[ z ].data( x, y )) * smoothFactor;
					if (x < width - 1) b1 += vci.uDx[ z ].data( x + 1, y ) * smoothFactorPx;
					if (y < height - 1) b1 += vci.uDy[ z ].data( x, y + 1 ) * smoothFactorPy;
					if (z < depth - 1) b1 += vci.uDz[ z + 1 ].data( x, y ) * smoothFactorPz;

					// smoothness: v eqn coefs
					dvCoef2 += 3.0f * smoothFactor;
					dvCoef2 += smoothFactorPx;
					dvCoef2 += smoothFactorPy;
					dvCoef2 += smoothFactorPz;
					float dvMxCoef = -smoothFactor;
					float dvMyCoef = -smoothFactor;
					float dvMzCoef = -smoothFactor;
					float dvPxCoef = -smoothFactorPx;
					float dvPyCoef = -smoothFactorPy;
					float dvPzCoef = -smoothFactorPz;
					b2 -= (vci.vDx[ z ].data( x, y ) + vci.vDy[ z ].data( x, y ) + vci.vDz[ z ].data( x, y )) * smoothFactor;
					if (x < width - 1) b2 += vci.vDx[ z ].data( x + 1, y ) * smoothFactorPx;
					if (y < height - 1) b2 += vci.vDy[ z ].data( x, y + 1 ) * smoothFactorPy;
					if (z < depth - 1) b2 += vci.vDz[ z + 1 ].data( x, y ) * smoothFactorPz;
					
					// smoothness: w eqn coefs
					dwCoef3 += 3.0f * smoothFactor * wSmooth;
					dwCoef3 += smoothFactorPx * wSmooth;
					dwCoef3 += smoothFactorPy * wSmooth;
					dwCoef3 += smoothFactorPz * wSmooth;
					float dwMxCoef = -smoothFactor * wSmooth;
					float dwMyCoef = -smoothFactor * wSmooth;
					float dwMzCoef = -smoothFactor * wSmooth;
					float dwPxCoef = -smoothFactorPx * wSmooth;
					float dwPyCoef = -smoothFactorPy * wSmooth;
					float dwPzCoef = -smoothFactorPz * wSmooth;
					b3 -= (vci.wDx[ z ].data( x, y ) + vci.wDy[ z ].data( x, y ) + vci.wDz[ z ].data( x, y )) * smoothFactor * wSmooth;
					if (x < width - 1) b3 += vci.wDx[ z ].data( x + 1, y ) * smoothFactorPx * wSmooth;
					if (y < height - 1) b3 += vci.wDy[ z ].data( x, y + 1 ) * smoothFactorPy * wSmooth;
					if (z < depth - 1) b3 += vci.wDz[ z + 1 ].data( x, y ) * smoothFactorPz * wSmooth;

					// smoothness: handle edge cases 
					if (uMxInd == -1) {
						uMxInd = vMxInd = wMxInd = 0;
						duMxCoef = dvMxCoef = dwMxCoef = 0;
					}
					if (uPxInd == -1) {
						uPxInd = vPxInd = wPxInd = 0;
						duPxCoef = dvPxCoef = dwPxCoef = 0;
					}
					if (uMyInd == -1) {
						uMyInd = vMyInd = wMyInd = 0;
						duMyCoef = dvMyCoef = dwMyCoef = 0;
					}
					if (uPyInd == -1) {
						uPyInd = vPyInd = wPyInd = 0;
						duPyCoef = dvPyCoef = dwPyCoef = 0;
					}
					if (uMzInd == -1) {
						uMzInd = vMzInd = wMzInd = 0;
						duMzCoef = dvMzCoef = dwMzCoef = 0;
					}
					if (uPzInd == -1) {
						uPzInd = vPzInd = wPzInd = 0;
						duPzCoef = dvPzCoef = dwPzCoef = 0;
					}

					// add equations to sparse system
					sparseSystem.addEquation(  uInd,    vInd,    wInd,    uMxInd,   uPxInd,   uMyInd,   uPyInd,   uMzInd,   uPzInd, 
											  duCoef1, dvCoef1, dwCoef1, duMxCoef, duPxCoef, duMyCoef, duPyCoef, duMzCoef, duPzCoef, b1 );
					sparseSystem.addEquation(  vInd,    uInd,    wInd,    vMxInd,   vPxInd,   vMyInd,   vPyInd,   vMzInd,   vPzInd, 
											  dvCoef2, duCoef2, dwCoef2, dvMxCoef, dvPxCoef, dvMyCoef, dvPyCoef, dvMzCoef, dvPzCoef, b2 );
					sparseSystem.addEquation(  wInd,    uInd,    vInd,    wMxInd,   wPxInd,   wMyInd,   wPyInd,   wMzInd,   wPzInd, 
											  dwCoef3, duCoef3, dvCoef3, dwMxCoef, dwPxCoef, dwMyCoef, dwPyCoef, dwMzCoef, dwPzCoef, b3 );
				}
			}
		}
	}

	// set initial value for iterative solver
	sparseSystem.setInit( init );

	// solve the system
	sparseSystem.setMaxIter( solverIters );
	sparseSystem.setSORFactor( sorFactor );
	sparseSystem.enableAdaptive( adaptIterMask, adaptUpdateThresh );
	VectorF result = sparseSystem.solve();

	// copy solution from result vector into du and dv (returned results)
	for (int z = 0; z < depth; z++) {
		const ImageGrayU &mask = vci.mask[ z ];
		ImageGrayF &du = vci.du[ z ];
		ImageGrayF &dv = vci.dv[ z ];
		ImageGrayF &dw = vci.dw[ z ];
		const ImageGrayI &uIndex = vci.uIndex[ z ];
		const ImageGrayI &vIndex = vci.vIndex[ z ];
		const ImageGrayI &wIndex = vci.wIndex[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (mask.data( x, y )) {
					float duBounded = bound( result.data( uIndex.data( x, y ) ), -maxMove, maxMove );
					float dvBounded = bound( result.data( vIndex.data( x, y ) ), -maxMove, maxMove );
					float dwBounded = bound( result.data( wIndex.data( x, y ) ), -maxMove, maxMove );
					du.data( x, y ) = duBounded;
					dv.data( x, y ) = dvBounded;
					dw.data( x, y ) = dwBounded;
				}
			}
		}
	}
}


/// incrementally compute du, dv
void varInnerLoop( VarCorresInfo &vci, Config &varConf, bool fast ) {

	// get parameters
	int innerIters = varConf.readInt( "innerIters" );
	int verbosity = varConf.readInt( "verbosity" );
	bool robustDataTerm = true;
	bool robustSmoothnessTerm = false;

	// allocate and initialize workspace
	int depth = vci.mask.count();
	int width = vci.mask[ 0 ].width(), height = vci.mask[ 0 ].height();	
	initImageSeq( vci.dPsiData, width, height, depth, true, 1 );
	initImageSeq( vci.dPsiSmoothness, width, height, depth, true, 1 );

	// run inner loop
	for (int innerIter = 0; innerIter < innerIters; innerIter++) {

		// check for cancel
		if (checkCommandEvents())
			break;
	
		// compute dPsiData, dPsiSmoothness using current values of du, dv (skip x = 0 and y = 0 pixels)
		if (robustDataTerm || robustSmoothnessTerm) {
			for (int z = 1; z < depth; z++) {
				const ImageGrayU &mask = vci.mask[ z ];
				const ImageGrayF &du = vci.du[ z ];
				const ImageGrayF &dv = vci.dv[ z ];
				const ImageGrayF &dw = vci.dw[ z ];
				ImageGrayF &dPsiSmoothness = vci.dPsiSmoothness[ z ];
				ImageGrayF &dPsiData = vci.dPsiData[ z ];
				for (int y = 1; y < height; y++) {
					for (int x = 1; x < width; x++) {

						// compute dPsi terms for each pixel
						if (mask.data( x, y )) {

							// data term: compute dPsi for each channel
							if (robustDataTerm) {
								float diff = vci.ix[ z ].data( x, y ) * du.data( x, y ) 
										   + vci.iy[ z ].data( x, y ) * dv.data( x, y ) 
										   + vci.iz[ z ].data( x, y ) * dw.data( x, y ) 
										   + vci.ic[ z ].data( x, y );
								dPsiData.data( x, y ) = psiDeriv( diff * diff, true );
							}

							// smoothness term (using gradient of (u + du, v + dv, w + dw))
							if (robustSmoothnessTerm) {
								float uNewDx = vci.uDx[ z ].data( x, y ) + dx( du, x, y );
								float uNewDy = vci.uDy[ z ].data( x, y ) + dy( du, x, y );
								float uNewDz = vci.uDz[ z ].data( x, y ) + dz( vci.du, x, y, z );
								float vNewDx = vci.vDx[ z ].data( x, y ) + dx( dv, x, y );
								float vNewDy = vci.vDy[ z ].data( x, y ) + dy( dv, x, y );
								float vNewDz = vci.vDz[ z ].data( x, y ) + dz( vci.dv, x, y, z );
								float wNewDx = vci.wDx[ z ].data( x, y ) + dx( dw, x, y );
								float wNewDy = vci.wDy[ z ].data( x, y ) + dy( dw, x, y );
								float wNewDz = vci.wDz[ z ].data( x, y ) + dz( vci.dw, x, y, z );
								float sum = uNewDx * uNewDx + uNewDy * uNewDy + uNewDz * uNewDz 
										  + vNewDx * vNewDx + vNewDy * vNewDy + vNewDz * vNewDz 
										  + wNewDx * wNewDx + wNewDy * wNewDy + wNewDz * wNewDz;
								dPsiSmoothness.data( x, y ) = psiDeriv( sum, true );
							}
						}
					}
				}
			}
		}

		// diagnostics
		if (verbosity > 12) {
//			for (int k = 0; k < chanCount; k++)
//				DispStats( "dPsiData", &dPsiData->ref( k ));
			dispStats( "dPsiData", vci.dPsiData[ depth / 2 ] );
			dispStats( "dPsiSmoothness", vci.dPsiSmoothness[ depth / 2 ] );
		}

		// solve the system
		varSolveSystem( vci, varConf, fast );
	}
}


/// computes an update (du, dv) to the flow field
void varIteration( Array<ImageGrayF> &uSeq, Array<ImageGrayF> &vSeq, Array<ImageGrayF> &wSeq, 
				   const Array<ImageGrayF> &srcSeqScaled, 
				   const Array<ImageGrayF> &destSeqScaled, 
				   Config &varConf, bool fast ) {	

	// get config parameters
	int verbosity = varConf.readInt( "verbosity" );

	// info for this iteration	
	int depth = srcSeqScaled.count();
	int width = srcSeqScaled[ 0 ].width(), height = srcSeqScaled[ 0 ].height();
	VarCorresInfo vci;
	initImageSeq( vci.mask, width, height, depth, true, 255 );
	initImageSeq( vci.ix, width, height, depth, true, 0 );
	initImageSeq( vci.iy, width, height, depth, true, 0 );
	initImageSeq( vci.iz, width, height, depth, true, 0 );
	initImageSeq( vci.ic, width, height, depth, true, 0 );
	initImageSeq( vci.du, width, height, depth, true, 0 );
	initImageSeq( vci.dv, width, height, depth, true, 0 );
	initImageSeq( vci.dw, width, height, depth, true, 0 );
	initImageSeq( vci.uIndex, width, height, depth, true, 0 );
	initImageSeq( vci.vIndex, width, height, depth, true, 0 );
	initImageSeq( vci.wIndex, width, height, depth, true, 0 );

	// get derivatives of destination channels
	Array<ImageGrayF> destSeqDx, destSeqDy, destSeqDz;
	for (int z = 0; z < depth; z++) {
		destSeqDx.append( dx( destSeqScaled[ z ] ).release() );
		destSeqDy.append( dy( destSeqScaled[ z ] ).release() );
		destSeqDz.append( dz( destSeqScaled, z ).release() );
	}

	// loop over corres fields
	for (int z = 0; z < depth; z++) {

		// get flow components for quick reference
		ImageGrayF &u = uSeq[ z ];
		ImageGrayF &v = vSeq[ z ];
		ImageGrayF &w = wSeq[ z ];
		ImageGrayU &mask = vci.mask[ z ];

		// use u and v to compute values needed by inner loop
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {

				// if inside optimization region
				if (mask.data( x, y )) {
					float xProj = x + u.data( x, y );
					float yProj = y + v.data( x, y );
					float zProj = z + w.data( x, y );
					int xProjInt = (int) xProj;
					int yProjInt = (int) yProj;
					int zProjInt = (int) zProj;

					// if projects inside image (no penalty if project outside image)
					if (       ((xProj >= 0 && xProj <= width - 1)  || (xProjInt == 0 && width == 1))
							&& ((yProj >= 0 && yProj <= height - 1) || (yProjInt == 0 && height == 1))
							&& ((zProj >= 0 && zProj <= depth - 1)  || (zProjInt == 0 && depth == 1))) {
						vci.ix[ z ].data( x, y ) = interp( destSeqDx, xProj, yProj, zProj );
						vci.iy[ z ].data( x, y ) = interp( destSeqDy, xProj, yProj, zProj );
						vci.iz[ z ].data( x, y ) = interp( destSeqDz, xProj, yProj, zProj );
						float ic = interp( destSeqScaled, xProj, yProj, zProj ) 
								 - srcSeqScaled[ z ].data( x, y );
						vci.ic[ z ].data( x, y ) = ic;
					}
				}
			}
		}

		// compute flow derivatives
		vci.uDx.append( dx( u ).release() );
		vci.uDy.append( dy( u ).release() );
		vci.uDz.append( dz( uSeq, z ).release() );
		vci.vDx.append( dx( v ).release() );
		vci.vDy.append( dy( v ).release() );
		vci.vDz.append( dz( vSeq, z ).release() );
		vci.wDx.append( dx( w ).release() );
		vci.wDy.append( dy( w ).release() );
		vci.wDz.append( dz( wSeq, z ).release() );

		// compute index maps (indices of solver variables)
		ImageGrayI &uIndex = vci.uIndex[ z ];
		ImageGrayI &vIndex = vci.vIndex[ z ];
		ImageGrayI &wIndex = vci.wIndex[ z ];
		int startIndex = vci.pixelCount * 3;
		vci.pixelCount += buildIndexMaps( startIndex, mask, uIndex, vIndex, wIndex );
	}

	// run inner loop to compute du, dv
	varInnerLoop( vci, varConf, fast );

	// increment u and v according to du and dv
	for (int z = 0; z < depth; z++) {
		const ImageGrayU &mask = vci.mask[ z ];
		ImageGrayF &u = uSeq[ z ];
		ImageGrayF &v = vSeq[ z ];
		ImageGrayF &w = wSeq[ z ];
		const ImageGrayF &du = vci.du[ z ];
		const ImageGrayF &dv = vci.dv[ z ];
		const ImageGrayF &dw = vci.dw[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (mask.data( x, y )) {
					u.data( x, y ) = u.data( x, y ) + du.data( x, y );
					v.data( x, y ) = v.data( x, y ) + dv.data( x, y );
					w.data( x, y ) = w.data( x, y ) + dw.data( x, y );
				}
			}
		}
	}

	// display diagnostics
	if (verbosity >= 8) {
		if (depth <= 2) {
			for (int z = 0; z < depth; z++) 
				vci.dispStats( 3, z );
		} else {
			vci.dispStats( 3, depth / 2 );
		}
	}
	if (verbosity >= 6) {
		dispStats( "u", uSeq );
		dispStats( "v", vSeq );
		dispStats( "w", wSeq );
	}
}


/// resize and blur an image volume (image sequence)
void resizeAndBlur( const ImageGrayFSeq &inSeq, ImageGrayFSeq &outSeq, int newWidth, int newHeight, int newLength, float sigma ) {
	ImageGrayFSeq tempSeq1, tempSeq2; // note: could use outSeq as tempSeq1 (unless inSeq is same as outSeq)
	for (int z = 0; z < inSeq.count(); z++) 
		tempSeq1.append( resize( inSeq[ z ], newWidth, newHeight, true ).release() );
	outSeq.reset(); // needed for corres field components where outSeq is same as inSeq
	if (sigma) {
		resizeSeqZ( tempSeq1, tempSeq2, newLength );
		blurGaussSeqZ( tempSeq2, outSeq, sigma );
		blurGaussSeqXY( outSeq, sigma );
		assertAlways( outSeq.count() == newLength );
		assertAlways( outSeq[ 0 ].width() == newWidth && outSeq[ 0 ].height() == newHeight );
	} else {
		resizeSeqZ( tempSeq1, outSeq, newLength );
	}
}


/// loops over scales; for each scale; executes outer loop
void varScaleLoop( const Array<ImageGrayF> &srcSeq, const Array<ImageGrayF> &destSeq, 
				   Array<CorresField3D> &cfSeq, Config &varConf ) {

	// get config parameters
	float sigma = varConf.readFloat( "sigma" );
	float scaleFactor = varConf.readFloat( "scaleFactor" );
	float minScale = varConf.readFloat( "minScale" );
	int verbosity = varConf.readInt( "verbosity" );
	bool fast = true; // check this

	// compute scales 
	VectorF scales = varScaleSequence( scaleFactor, minScale );
	int scaleCount = scales.length();
	if (verbosity >= 4)
		disp( 1, "iterating over %d scales (scale-factor = %f)", scaleCount, scaleFactor );

	// will hold correspondence field components
	ImageGrayFSeq uSeq, vSeq, wSeq;

	// loop over scales (from smallest to largest)
	for (int scaleIndex = scaleCount - 1; scaleIndex >= 0; scaleIndex--) {
		float scale = scales[ scaleIndex ];

		// compute motion and image dimensions
		int smallWidth = sbl::round( srcSeq[ 0 ].width() * scale );
		int smallHeight = sbl::round( srcSeq[ 0 ].height() * scale );
		int smallLength = sbl::round( srcSeq.count() * scale );
		if (verbosity >= 4)
			disp( 2, "scale %f: %d x %d x %d", scale, smallWidth, smallHeight, smallLength );

		// shrink and blur the channels to current size
		Array<ImageGrayF> srcSeqScaled, destSeqScaled;
		if (scales.length() > 1) {
			resizeAndBlur( srcSeq, srcSeqScaled, smallWidth, smallHeight, smallLength, sigma );
			resizeAndBlur( destSeq, destSeqScaled, smallWidth, smallHeight, smallLength, sigma );
		}
	
		// create or resize corres fields
		if (uSeq.count() == 0) {
			initImageSeq( uSeq, smallWidth, smallHeight, smallLength, true, 0 );
			initImageSeq( vSeq, smallWidth, smallHeight, smallLength, true, 0 );
			initImageSeq( wSeq, smallWidth, smallHeight, smallLength, true, 0 );
		} else {
			resizeAndBlur( uSeq, uSeq, smallWidth, smallHeight, smallLength, sigma );
			resizeAndBlur( vSeq, vSeq, smallWidth, smallHeight, smallLength, sigma ); 
			resizeAndBlur( wSeq, wSeq, smallWidth, smallHeight, smallLength, sigma ); 
		}

		// run outer iteration (inner loop)
		if (scales.length() > 1) 
			varIteration( uSeq, vSeq, wSeq, srcSeqScaled, destSeqScaled, varConf, fast );
		else
			varIteration( uSeq, vSeq, wSeq, srcSeq, destSeq, varConf, fast );
	}

	// populate cf seq from components
	int length = srcSeq.count();
	int width = srcSeq[ 0 ].width(), height = srcSeq[ 0 ].height();
	for (int z = 0; z < length; z++) {
		CorresField3D *cf = new CorresField3D( width, height );
		const ImageGrayF &u = uSeq[ z ];
		const ImageGrayF &v = vSeq[ z ];
		const ImageGrayF &w = wSeq[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				cf->u( x, y ) = u.data( x, y );
				cf->v( x, y ) = v.data( x, y );
				cf->w( x, y ) = w.data( x, y );
			}
		}
		cfSeq.append( cf );
	}
}


//-------------------------------------------
// TOP-LEVEL VARIATIONAL ALGORITHM
//-------------------------------------------


/// estimate a 3D mapping (analogous to optical flow) from one image volume to another image volume
void varCorres3D( const Array<ImageGrayU> &src, const Array<ImageGrayU> &dest, 
				  Array<CorresField3D> &cfSeq, int scaleFactor ) {

	// load config with parameters for this algorithm
	String configFileName = "varCorres3D.conf";
	Config varConf;
	varConf.load( configFileName );
	if (varConf.entryCount() == 0) {
		warning( "failed to load varCorres3D config: %s", configFileName.c_str() );
		return;
	}

	// check inputs
	assertAlways( src.count() && src.count() == dest.count() );
	assertAlways( src[ 0 ].width() == dest[ 0 ].width() && src[ 0 ].height() == dest[ 0 ].height() );
	assertAlways( cfSeq.count() == 0 );
	assertAlways( scaleFactor >= 1 );

	// determine working image size
	int newWidth = src[ 0 ].width() / scaleFactor;
	int newHeight = src[ 0 ].height() / scaleFactor;

	// prepare source images
	Array<ImageGrayF> srcSeq;
	for (int i = 0; i < src.count(); i++) {
		aptr<ImageGrayU> srcSmall = resize( src[ i ], newWidth, newHeight, true );
		srcSeq.append( toFloat( *srcSmall, 1.0f / 255.0f ).release() );
	}

	// prepare dest images
	Array<ImageGrayF> destSeq;
	for (int i = 0; i < dest.count(); i++) {
		aptr<ImageGrayU> destSmall = resize( dest[ i ], newWidth, newHeight, true );
		destSeq.append( toFloat( *destSmall, 1.0f / 255.0f ).release() );
	}

	// run motion estimation
	varScaleLoop( srcSeq, destSeq, cfSeq, varConf );

	// resize flow to full image size
	for (int i = 0; i < cfSeq.count(); i++) 
		cfSeq[ i ].resize( src[ 0 ].width(), src[ 0 ].height(), true );
}


} // end namespace hb
