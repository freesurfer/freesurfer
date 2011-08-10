#ifndef _CORRES_FIELD_3D_H_
#define _CORRES_FIELD_3D_H_
#include <sbl/core/Pointer.h>
#include <sbl/core/File.h>
#include <sbl/image/Image.h>
using namespace sbl;
namespace hb {


//-------------------------------------------
// CORRES FIELD 3D CLASS
//-------------------------------------------

	
/// The CorresField3D class represents a 2D field of 3D vectors, representing offsets from a volume slice into another volume.
class CorresField3D {
public:

	/// create an uninitialized correspondence field
	CorresField3D( int width, int height );

	/// load correspondence field from file
	CorresField3D( File &file );

	// basic destructor
	~CorresField3D();

	/// access correspondence field images
	inline ImageGrayF &u() { return *m_u; }
	inline ImageGrayF &v() { return *m_v; }
	inline ImageGrayF &w() { return *m_w; }
	inline const ImageGrayF &u() const { return *m_u; }
	inline const ImageGrayF &v() const { return *m_v; }
	inline const ImageGrayF &w() const { return *m_w; }

	/// get correspondence field elements
	inline float u( int x, int y ) const { return m_u->data( x, y ); }
	inline float v( int x, int y ) const { return m_v->data( x, y ); }
	inline float w( int x, int y ) const { return m_w->data( x, y ); }

	/// access correspondence field elements
	inline float &u( int x, int y ) { return m_u->data( x, y ); }
	inline float &v( int x, int y ) { return m_v->data( x, y ); }
	inline float &w( int x, int y ) { return m_w->data( x, y ); }

	/// set all components at a single location
	inline void set( int x, int y, float u, float v, float w ) { m_u->data( x, y ) = u; m_v->data( x, y ) = v; m_w->data( x, y ) = w; }

	/// dimensions of correspondence field
	inline int width() const { return m_u->width(); }
	inline int height() const { return m_u->height(); }

	/// set all correspondence vectors to zero
	void clear();

	/// shrink/zoom correspondence field
	void resize( int width, int height, bool rescale );

	/// save correspondence field to file
	void save( File &file ) const;

	/// create a color visualization of the correspondence vectors
	aptr<ImageColorU> colorize();

private:

	// an image for each 
	ImageGrayF *m_u;
	ImageGrayF *m_v;
	ImageGrayF *m_w;

	// disable copy constructor and assignment operator
	CorresField3D( const CorresField3D &x );
	CorresField3D &operator=( const CorresField3D &x );
};


//-------------------------------------------
// CORRES FIELD UTILS
//-------------------------------------------


/// compute mean difference between the correspondence vectors
float meanAbsDiff( const CorresField3D &cf1, const CorresField3D &cf2 );


/// display statistics about the correspondence vectors
void dispStats( int indent, const Array<CorresField3D> &cfSeq );


/// maps the image sequence according to the correspondences 
void mapBack( const Array<CorresField3D> &cfSeq, 
			  const Array<ImageGrayU> &seq, Array<ImageGrayU> &seqMapped );


/// compute the inverse mapping of the given correspondence volume; assumes the correspodnences are smooth
Array<CorresField3D> invert( const Array<CorresField3D> &cfSeq );


/// load a correspondence volume from file
Array<CorresField3D> loadCorresSeq( const String &fileName );


/// save a correspondence volume to file
void saveCorresSeq( const Array<CorresField3D> &cfSeq, const String &fileName );


} // end namespace hb
#endif // _CORRES_FIELD_3D_H_
