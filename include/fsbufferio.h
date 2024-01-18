#ifndef FSBUFFERIO_H
#define FSBUFFERIO_H

struct MATRIX;
struct VOL_GEOM;

class FSbufferIO
{
public:
  FSbufferIO(int size);
  ~FSbufferIO();

  int getbufferlen() { return fsbuflen; }
  unsigned char *getbuffer() { return fsbuf; }
  
  int write_int(int v);
  int write_long(long long v);
  int write_float(float f);
  int write_double(double d);
  int write_buf(void *buf, int buflen);

  int write_matrix(MATRIX *M);
  int write_geom(VOL_GEOM *volgeom);

private:
  void _checkbufsize(int buflen);
  
private:
  int fsbufsize;
  int fsbuflen;
  
  unsigned char *fsbuf;

};

#endif // FSBUFFERIO_H
