#ifndef TAGBUFFERIO_H
#define TAGBUFFERIO_H

// forward declarations
struct VOL_GEOM;
class  MRI;
struct MATRIX;
struct COLOR_TABLE;
class  FSbufferIO;

class TAGbufferIO
{
public:
  TAGbufferIO(int niftiheaderextension=false);
  ~TAGbufferIO();

  int getbufferlen() { return tagbuflen; }
  unsigned char *getbuffer();

  // write data to buffer allocated for the object
  int writedata(void *data, long long dlen);

  // methods to write various TAGs
  int writetag(int tag, void *data, long long dlen);
  int writematrix(MATRIX *M, int tag);
  int write_old_colortable(COLOR_TABLE *ctab);
  int write_mri_frames(MRI *mri);
  
  int write_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target);
  int write_gcamorph_meta(int warpFieldFormat, int gcamorphSpacing, double gcamorphExp_k);
  int write_gcamorph_labels(int x, int y, int z, int   ***gcamorphLabel);
  
private:
  //int tagbufsize;  // 16k
  int tagbuflen;
  
  //unsigned char *tagbuf;

  FSbufferIO *fsbufio;
};

#endif // TAGBUFFERIO_H
