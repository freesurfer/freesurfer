#ifndef DTK_FS_INC
#define DTK_FS_INC

#ifdef X
#undef X
#endif

#define DTK_N_SCALARS_MAX 20

/*---------------------------------------------------------*/
// http://trackvis.org/docs/?subsect=fileformat
/*---------------------------------------------------------*/
typedef struct
{
  int npoints; // number of points in the track
  float *c, *r, *s; // Col, Row, Slice coordinates
  int n_scalars; // May be inherited from HDR
  float *scalars[DTK_N_SCALARS_MAX]; // n_s by npoints
  int *label;
  int n_properties;
  float *properties; // May be inherited from HDR
} DTK_TRACK;

// DO NOT ADD, REMOVE, OR REARRANGE ANYTHING FROM THIS STRUCT!!!!
typedef struct
{
  char id_string[6];
  short int dim[3];
  float voxel_size[3];
  float origin[3];
  short int n_scalars;
  char scalar_name[10][20];
  short int n_properties;
  char property_name[10][20];
  char reserved[508];
  char voxel_order[4];
  char pad2[4];
  float image_orientation_patient[6];
  char pad1[2];
  unsigned char invert_x;
  unsigned char invert_y;
  unsigned char swap_xy;
  unsigned char swap_yz;
  unsigned char swap_zx;
  int n_count;
  int version;
  int hdr_size;
}  DTK_HDR;

typedef struct
{
  DTK_HDR *hdr; // header
  DTK_TRACK     **trk; // array of tracks (n_count)
  MRI           *mri_template; // only need geometry
} DTK_TRACK_SET;

const char *DTKFSSrcVersion(void);
DTK_TRACK_SET *DTKloadTrackSet(char *trkfile, char *mrifile);
DTK_TRACK *DTKreadTrack(FILE *fp, int n_scalars, int n_properties, 
			float dc, float dr, float ds);
DTK_TRACK *DTKallocTrack(int npoints, int n_scalars, int n_properties);
int DTKprintTrack(FILE *fp, DTK_TRACK *trk);
int DTKprintHeader(FILE *fp, DTK_HDR *dtkhdr);
MRI *DTKmapEndPoints(DTK_TRACK_SET *dtkset);
int DTKwriteTrack(FILE *fp, DTK_TRACK *trk,float dc, float dr, float ds);
int DTKwriteTrackSet(char *trkfile, DTK_TRACK_SET *dtkset);
int DTKlabelTracks(DTK_TRACK_SET *trkset, MRI *seg);
DTK_TRACK_SET *DTKextractCC(DTK_TRACK_SET *trkset);
MRI *DTKmapTrackNos(DTK_TRACK_SET *trkset);
MRI *DTKsegROI(DTK_TRACK_SET *trkset, MRI *seg, int segid);
DTK_TRACK_SET *DTKextractSeg(DTK_TRACK_SET *trkset, int segid);
DTK_TRACK_SET *DTKextractSegEndPoints(DTK_TRACK_SET *trkset, int segid);

//DTKwriteTrackSet(char *fname, DTK_TRACK_SET *dtkset)
//DTKfree(DTK_TRACK_SET **pdtkset)
//DTK_HDR *DTKloadHeader(char *trkfile);
//DTKallocTracks(DTK_TRACK_SET *dtkset)
//DTKloadTracks(FILE *fp, DTK_TRACK_SET *dtkset)
//DTKtoCRS(DTK_TRACK_SET *dtkset)
//DTKtoXYZ(DTK_TRACK_SET *dtkset)


#endif
