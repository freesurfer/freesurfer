#include "CorresField3D.h"
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageSeqUtil.h>
#include <sbl/image/ImageTransform.h>
namespace hb {


//-------------------------------------------
// CORRES FIELD 3D CLASS
//-------------------------------------------


/// create an uninitialized correspondence field
CorresField3D::CorresField3D( int width, int height ) {
	m_u = new ImageGrayF( width, height );
	m_v = new ImageGrayF( width, height );
	m_w = new ImageGrayF( width, height );
}


/// load correspondence field from file
CorresField3D::CorresField3D( File &file ) {
	int width = file.readInt();
	int height = file.readInt();
	assertAlways( width > 0 && width < 100000 );
	assertAlways( height > 0 && height < 100000 );
	m_u = new ImageGrayF( width, height );
	m_v = new ImageGrayF( width, height );
	m_w = new ImageGrayF( width, height );
	for (int y = 0; y < height; y++) 
		for (int x = 0; x < width; x++)
			m_u->data( x, y ) = file.readFloat();
	for (int y = 0; y < height; y++) 
		for (int x = 0; x < width; x++)
			m_v->data( x, y ) = file.readFloat();
	for (int y = 0; y < height; y++) 
		for (int x = 0; x < width; x++)
			m_w->data( x, y ) = file.readFloat();
}


// basic destructor
CorresField3D::~CorresField3D() {
	delete m_u;
	delete m_v;
	delete m_w;
}


/// set all correspondence vectors to zero
void CorresField3D::clear() {
	m_u->clear( 0 );
	m_v->clear( 0 );
	m_w->clear( 0 );
}


/// shrink/zoom correspondence field
void CorresField3D::resize( int width, int height, bool rescale ) {

	// if no scale change, don't do anything
	if (width == m_u->width() && height == m_u->height())
		return;

	// rescale if requested
	if (rescale) {
		float xScaleFactor = (float) width / (float) m_u->width();
		float yScaleFactor = (float) height / (float) m_v->height();
		multiply( *m_u, xScaleFactor, *m_u );
		multiply( *m_v, yScaleFactor, *m_v ); // note: we don't rescale w(x,y)
	}

	// resize
	m_u = sbl::resize( *m_u, width, height, true ).release();
	m_v = sbl::resize( *m_v, width, height, true ).release();
	m_w = sbl::resize( *m_w, width, height, true ).release();
}


/// save correspondence field to file
void CorresField3D::save( File &file ) const {
	int width = m_u->width();
	int height = m_u->height();
	file.writeInt( width );
	file.writeInt( height );
	for (int y = 0; y < height; y++) 
		for (int x = 0; x < width; x++)
			file.writeFloat( m_u->data( x, y ) );
	for (int y = 0; y < height; y++) 
		for (int x = 0; x < width; x++)
			file.writeFloat( m_v->data( x, y ) );
	for (int y = 0; y < height; y++) 
		for (int x = 0; x < width; x++)
			file.writeFloat( m_w->data( x, y ) );
}


/// create a color visualization of the correspondence vectors
aptr<ImageColorU> CorresField3D::colorize() {
	int width = m_u->width(), height = m_u->height();
	aptr<ImageColorU> vis( new ImageColorU( width, height ) );

	// compute max motion
	float max = 0;
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++) {
			float val = m_u->data( x, y );
			if (val < 0) val = -val;
			if (val > max) max = val;
			val = m_v->data( x, y );
			if (val < 0) val = -val;
			if (val > max) max = val;
			val = m_w->data( x, y );
			if (val < 0) val = -val;
			if (val > max) max = val;
		}

	// create image
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int u = bound( 128 + sbl::round( m_u->data( x, y ) * 128.0f / max ), 0, 255 );
			int v = bound( 128 + sbl::round( m_v->data( x, y ) * 128.0f / max ), 0, 255 );
			int w = bound( 128 + sbl::round( m_w->data( x, y ) * 128.0f / max ), 0, 255 );
			vis->setRGB( x, y, u, v, w );
		}
	}

	// return motion plot
	return vis; 
}


//-------------------------------------------
// CORRES FIELD UTILS
//-------------------------------------------


/// map source images forward according to corres seq using a splatting method
void mapForward( const Array<CorresField3D> &corresSeq, float zOffset,
				 const Array<ImageGrayU> &source, Array<ImageGrayU> &dest ) {

	// get info for quick reference
	int sourceCount = source.count();
	int destCount = dest.count();
	assertAlways( sourceCount && sourceCount == corresSeq.count() );
	assertAlways( destCount );
	int sourceWidth = source[ 0 ].width();
	int sourceHeight = source[ 0 ].height();
	int destWidth = dest[ 0 ].width();
	int destHeight = dest[ 0 ].height();

	// loop from large splats to small splats
	for (int radius = 2; radius >= 0; radius--) {

		// loop over source slices
		for (int z = 0; z < sourceCount; z++) {
			const ImageGrayU &sourceImage = source[ z ];
			const ImageGrayF &u = corresSeq[ z ].u();
			const ImageGrayF &v = corresSeq[ z ].v();
			const ImageGrayF &w = corresSeq[ z ].w();

			// loop over this source slice, splatting into dest
			for (int y = 0; y < sourceHeight; y++) {
				for (int x = 0; x < sourceWidth; x++) {
					int xProj = sbl::round( x + u.data( x, y ) );
					int yProj = sbl::round( y + v.data( x, y ) );
					int zProj = sbl::round( z + w.data( x, y ) + zOffset );

					// if projects in bounds, copy small patch to dest
					if (xProj > 2 && xProj < destWidth - 2 && yProj > 2 && yProj < destHeight - 2 && zProj >= 0 && zProj < destCount) {
						ImageGrayU &destImage = dest[ zProj ];
						for (int yOffset = -radius; yOffset <= radius; yOffset++) {
							for (int xOffset = -radius; xOffset <= radius; xOffset++) {
								destImage.data( xProj + xOffset, yProj + yOffset ) = sourceImage.data( x, y );
							}
						}
					}
				}
			}
		}
	}
}


/// compute mean difference between the correspondence vectors
float meanAbsDiff( const CorresField3D &cf1, const CorresField3D &cf2 ) {
	int width = cf1.width(), height = cf2.height();
	assertAlways( width == cf2.width() && height == cf2.height() );
	double sum = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			sum += fAbs( cf1.u( x, y ) - cf2.u( x, y ) );
			sum += fAbs( cf1.v( x, y ) - cf2.v( x, y ) );
			sum += fAbs( cf1.w( x, y ) - cf2.w( x, y ) );
		}
	}
	return (float) (sum / (double) (width * height * 3));
}


/// display statistics about the correspondence vectors
void dispStats( int indent, const Array<CorresField3D> &cfSeq ) {
	int depth = cfSeq.count();
	int width = cfSeq[ 0 ].width(), height = cfSeq[ 0 ].height();
	double uSum = 0, uMin = 1e10, uMax = -1e10;
	double vSum = 0, vMin = 1e10, vMax = -1e10;
	double wSum = 0, wMin = 1e10, wMax = -1e10;
	for (int z = 0; z < depth; z++) {
		const CorresField3D &cf = cfSeq[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				double u = cf.u( x, y );
				double v = cf.v( x, y );
				double w = cf.w( x, y );
				uSum += u;
				vSum += v;
				wSum += w;
				if (u < uMin) uMin = u;
				if (u > uMax) uMax = u;
				if (v < vMin) vMin = v;
				if (v > vMax) vMax = v;
				if (w < wMin) wMin = w;
				if (w > wMax) wMax = w;
			}
		}
	}
	double count = (double) (width * height * depth);
	disp( indent, "u: %f/%f/%f", uMin, uSum / count, uMax );
	disp( indent, "v: %f/%f/%f", vMin, vSum / count, vMax );
	disp( indent, "w: %f/%f/%f", wMin, wSum / count, wMax );
}


/// maps the image sequence according to the correspondences 
void mapBack( const Array<CorresField3D> &cfSeq, 
			  const Array<ImageGrayU> &seq, Array<ImageGrayU> &seqMapped ) {

	// get dimensions
	int length = cfSeq.count();
	assertAlways( length );
	assertAlways( length == seq.count() && length == seqMapped.count() );
	int width = cfSeq[ 0 ].width(), height = cfSeq[ 0 ].height();
	assertAlways( width == seq[ 0 ].width() && width == seqMapped[ 0 ].width() );
	assertAlways( height == seq[ 0 ].height() && height == seqMapped[ 0 ].height() );

	// loop over sequence, mapping each pixel
	for (int z = 0; z < length; z++) {
		const CorresField3D &cf = cfSeq[ z ];
		ImageGrayU &imageMapped = seqMapped[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				float xDest = x + cf.u( x, y );
				float yDest = y + cf.v( x, y );
				float zDest = z + cf.w( x, y );
				int v = 0;
				if (zDest >= 0 && zDest <= length - 1 && imageMapped.inBounds( xDest, yDest )) {
					v = round( interp( seq, xDest, yDest, zDest ) );
				}
				imageMapped.data( x, y ) = v;
			}
		}
	}
}


/// compute the inverse mapping of the given correspondence volume; assumes the correspodnences are smooth
Array<CorresField3D> invert( const Array<CorresField3D> &cfSeq ) {

	// initialize a sequence with the same dimensions as input and filled with a special marker value
	float marker = -7777.77f;
	int depth = cfSeq.count();
	assertAlways( depth );
	int width = cfSeq[ 0 ].width();
	int height = cfSeq[ 0 ].height();
	Array<CorresField3D> cfSeqInv;
	for (int z = 0; z < depth; z++) {
		CorresField3D *cf = new CorresField3D( width, height );
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				cf->u( x, y ) = marker;
				cf->v( x, y ) = marker;
				cf->w( x, y ) = marker;
			}
		}
		cfSeqInv.append( cf );
	}

	// map forward
	for (int z = 0; z < depth; z++) {
		const CorresField3D &cf = cfSeq[ z ];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				float u = cf.u( x, y );
				float v = cf.v( x, y );
				float w = cf.w( x, y );
				float xDest = x + u;
				float yDest = y + v;
				float zDest = z + w;
				int xDestInt = (int) xDest;
				int yDestInt = (int) yDest;
				int zDestInt = (int) zDest;
				if (zDestInt >= 0 && zDestInt < depth) {
					if (xDestInt >= 0 && xDestInt < width && yDestInt >= 0 && yDestInt < height) {
						cfSeqInv[ zDestInt ].set( xDestInt, yDestInt, -u, -v, -w );
/*
						if (xDestInt + 1 < width) 
							cfSeqInv[ zDestInt ].set( xDestInt + 1, yDestInt, -u, -v, -w );
						if (yDestInt + 1 < height) {
							cfSeqInv[ zDestInt ].set( xDestInt, yDestInt + 1, -u, -v, -w );
							if (xDestInt + 1 < width) 
								cfSeqInv[ zDestInt ].set( xDestInt + 1, yDestInt + 1, -u, -v, -w );
						}
						if (zDestInt + 1 < depth) {
							cfSeqInv[ zDestInt + 1 ].set( xDestInt, yDestInt, -u, -v, -w );
							if (xDestInt + 1 < width) 
								cfSeqInv[ zDestInt + 1 ].set( xDestInt + 1, yDestInt, -u, -v, -w );
							if (yDestInt + 1 < height) {
								cfSeqInv[ zDestInt + 1 ].set( xDestInt, yDestInt + 1, -u, -v, -w );
								if (xDestInt + 1 < width) 
									cfSeqInv[ zDestInt + 1 ].set( xDestInt + 1, yDestInt + 1, -u, -v, -w );
							}
						}
*/
					}
				}
			}
		}
	}

	// fill holes using neighbors
	int emptyCount = 0;
	int lastEmptyCount = -1;
	do {
		emptyCount = 0;
		for (int z = 0; z < depth; z++) {
			CorresField3D &cf = cfSeqInv[ z ];
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					if (cf.u( x, y ) == marker) {

						// neighbor stats
						int foundCount = 0;
						double uSum = 0, vSum = 0, wSum = 0;

						// compute bounds of neighbor region
						int dxMin = -1, dxMax = 1, dyMin = -1, dyMax = 1, dzMin = -1, dzMax = 1;
						if (x == 0) dxMin = 0;
						if (x == width - 1) dxMax = 0;
						if (y == 0) dyMin = 0;
						if (y == height - 1) dyMax = 0;
						if (z == 0) dzMin = 0;
						if (z == depth - 1) dzMax = 0;

						// loop over neighbors
						for (int dz = dzMin; dz <= dzMax; dz++) {
							CorresField3D &cfOther = cfSeqInv[ z + dz ];
							for (int dy = dyMin; dy <= dyMax; dy++) {
								for (int dx = dxMin; dx <= dxMax; dx++) {
									if (cfOther.u( x + dx, y + dy ) != marker) {
										uSum += cfOther.u( x + dx, y + dy );
										vSum += cfOther.v( x + dx, y + dy );
										wSum += cfOther.w( x + dx, y + dy );
										foundCount++;
									}
								}
							}
						}
						if (foundCount) {
							cf.u( x, y ) = (float) (uSum / (double) foundCount);
							cf.v( x, y ) = (float) (vSum / (double) foundCount);
							cf.w( x, y ) = (float) (wSum / (double) foundCount);
						} else {
							emptyCount++;
						}
					}
				}
			}
		}
		if (emptyCount == lastEmptyCount) {
			fatalError( "unable to fill holes when inverting correspondence volume" );
		}
		lastEmptyCount = emptyCount;
	} while (emptyCount);
	return cfSeqInv;
}


/// load a correspondence volume from file
Array<CorresField3D> loadCorresSeq( const String &fileName ) {
	Array<CorresField3D> cfSeq;
	File file( fileName, FILE_READ, FILE_BINARY );
	if (file.openSuccess())
		file.readArray( cfSeq );
	return cfSeq;
}


/// save a correspondence volume to file
void saveCorresSeq( const Array<CorresField3D> &cfSeq, const String &fileName ) {
	File file( fileName, FILE_WRITE, FILE_BINARY );
	if (file.openSuccess())
		file.writeArray( cfSeq );
}


} // end namespace hb
