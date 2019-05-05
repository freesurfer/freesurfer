// Included by mrisurf.h
//
#if defined(COMPILING_MRISURF_TOPOLOGY) || defined(COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED)
#define CONST_EXCEPT_MRISURF_TOPOLOGY 
#else
#define CONST_EXCEPT_MRISURF_TOPOLOGY const
#endif

#if defined(COMPILING_MRISURF_METRIC_PROPERTIES) || defined(COMPILING_MRISURF_METRIC_PROPERTIES_FRIEND)
#define CONST_EXCEPT_MRISURF_METRIC_PROPERTIES 
#else
#define CONST_EXCEPT_MRISURF_METRIC_PROPERTIES const
#endif
    //
    // Used to find and control where various fields are written
    
struct face_topology_type_ {    // not used much yet
  vertices_per_face_t v;
};

struct face_type_ 
{
#define LIST_OF_FACE_ELTS_1    \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY vertices_per_face_t,v) SEP               /* vertex numbers of this face */    \
  ELTT(float,area) SEP    \
  ELTT(angles_per_triangle_t,angle) SEP    \
  ELTT(angles_per_triangle_t,orig_angle) SEP    \
  ELTT(char,ripflag) SEP                        /* ripped face */    \
  ELTT(char,oripflag) SEP                       /* stored version */    \
  ELTT(int,marked) SEP                          /* marked face */    \
  ELTP(DMATRIX,norm) SEP  /* 3x1 normal vector */ \
  ELTP(DMATRIX,gradNorm[3]) SEP  /* 3x3 Gradient of the normal wrt each of the 3 vertices*/ 
    // end of macro

#if 0
  float logshear,shearx,sheary;  /* compute_shear */
#endif

// Why does mrishash need these?  Where else are they used?
#if 0
#define LIST_OF_FACE_ELTS_2    \
  ELTT(float,cx) SEP    \
  ELTT(float,cy) SEP    \
  ELTT(float,cz) SEP         /* coordinates of centroid */   \
    // end of macro
#define LIST_OF_FACE_ELTS \
    LIST_OF_FACE_ELTS_1 SEP \
    LIST_OF_FACE_ELTS_2 \
    // end of macro
#else
#define LIST_OF_FACE_ELTS \
    LIST_OF_FACE_ELTS_1
#endif

#define ELTT(T,N) T N;
#define ELTP(TARGET,NAME) TARGET *NAME ;
#define SEP
LIST_OF_FACE_ELTS
#undef SEP
#undef ELTT
#undef ELTP

};



#define LIST_OF_VERTEX_TOPOLOGY_ELTS \
  /* put the pointers before the ints, before the shorts, before uchars, to reduce size  */ \
  /* the whole fits in much less than one cache line, so further ordering is no use      */ \
  ELTP(int,f) SEP                                               /* array[v->num] the fno's of the neighboring faces         */ \
  ELTP(uchar,n) SEP           	                                /* array[v->num] the face.v[*] index for this vertex        */ \
  ELTP(int,e) SEP                                               /* edge state for neighboring vertices                      */ \
  ELTP(CONST_EXCEPT_MRISURF_TOPOLOGY int,v) SEP                 /* array[v->vtotal or more] of vno, head sorted by hops     */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,vnum)                /* number of 1-hop neighbots    should use [p]VERTEXvnum(i) */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,v2num) SEP           /* number of 1, or 2-hop neighbors                          */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,v3num) SEP           /* number of 1,2,or 3-hop neighbors                         */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,vtotal) SEP          /* total # of neighbors. copy of vnum.nsizeCur              */ \
  ELTX(CONST_EXCEPT_MRISURF_TOPOLOGY short,nsizeMaxClock) SEP   /* copy of mris->nsizeMaxClock when v#num                   */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY uchar,nsizeMax) SEP        /* the max nsize that was used to fill in vnum etc          */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY uchar,nsizeCur) SEP        /* index of the current v#num in vtotal                     */ \
  ELTT(uchar,num) SEP                                           /* number of neighboring faces                              */ \
  // end of macro

//static short* pVERTEXvnum(VERTEX_TOPOLOGY * v, int i);
//static short VERTEXvnum(VERTEX_TOPOLOGY const * v, int i);

// The above elements historically were in the VERTEX
// and can still be there by
//  having VERTEX_TOPOLOGY be a typedef of VERTEX
//  having the mris->vertices and the mris->vertices_topology be the same pointer
// and this is what the code was doing until the separation was completed.
//
#define SEPARATE_VERTEX_TOPOLOGY
#ifndef SEPARATE_VERTEX_TOPOLOGY

#define LIST_OF_VERTEX_TOPOLOGY_ELTS_IN_VERTEX LIST_OF_VERTEX_TOPOLOGY_ELTS SEP

#else

struct VERTEX_TOPOLOGY {
    // The topology of the vertex describes its neighbors
    // but not its position nor any properties derived from its position.
    //
    // This data is not changed as the vertices are moved during distortion of the polyhedra.
    //

#define SEP
#define ELTX(TYPE,NAME) TYPE NAME ;
#define ELTT(TYPE,NAME) TYPE NAME ;
#define ELTP(TARGET,NAME) TARGET *NAME ;
  LIST_OF_VERTEX_TOPOLOGY_ELTS
#undef ELTP
#undef ELTT
#undef ELTX
#undef SEP
};

#define LIST_OF_VERTEX_TOPOLOGY_ELTS_IN_VERTEX

#endif


struct vertex_type_
{
// The LIST_OF_VERTEX_ELTS macro used here enables the the mris_hash
// and other algorithms to process all the elements without having to explicitly name them there and here
//
// The order is important because
//      . it affects the hash (until a better hash algorithm is implemented)
//      . it affects the size (each item must be aligned on its appropriate boundary for its size)
//      . it affects the number of cache lines that must be read and written
//
// By putting the ripflag at the end, reading it will cause its whole cache line to be read, perhaps the 
// first few elements of the next vertex earlier, and the ripflag test will probably be correctly predicted
// and so the cpu won't wait.
//
#define LIST_OF_VERTEX_ELTS_1    \
  LIST_OF_VERTEX_TOPOLOGY_ELTS_IN_VERTEX \
  \
  /* managed by MRISfreeDists[_orig] and MRISmakeDists[_orig] */ \
  ELTX(float* const,dist)      SEP                                              /* distance to neighboring vertices based on  xyz    */ \
  ELTX(float* const,dist_orig) SEP                                              /* distance to neighboring vertices based on origxyz */ \
  ELTX(int,dist_capacity)      SEP \
  ELTX(int,dist_orig_capacity) SEP \
  \
  ELTT(/*CONST_EXCEPT_MRISURF_METRIC_PROPERTIES*/ float,x)          SEP             /* current coordinates */                       \
  ELTT(/*CONST_EXCEPT_MRISURF_METRIC_PROPERTIES*/ float,y)          SEP             /* use MRISsetXYZ() to set */                   \
  ELTT(/*CONST_EXCEPT_MRISURF_METRIC_PROPERTIES*/ float,z)          SEP                                                             \
  \
  ELTT(const float,origx)                                       SEP             /* original coordinates, see also MRIS::origxyz_status */   \
  ELTT(const float,origy)                                       SEP             /* use MRISsetOriginalXYZ() */                      \
  ELTT(const float,origz)                                       SEP             /* or MRISsetOriginalXYZfromXYZ to set */           \
  \
  ELTT(float,nx) SEP    \
  ELTT(float,ny) SEP    \
  ELTT(float,nz) SEP        /* curr normal */    \
  ELTT(float,pnx) SEP    \
  ELTT(float,pny) SEP    \
  ELTT(float,pnz) SEP     /* pial normal */    \
  /* the above is the first cache line */ \
  ELTT(float,wnx) SEP    \
  ELTT(float,wny) SEP    \
  ELTT(float,wnz) SEP     /* white normal */    \
  ELTT(float,onx) SEP    \
  ELTT(float,ony) SEP    \
  ELTT(float,onz) SEP     /* original normal */    \
  ELTT(float,dx) SEP    \
  ELTT(float,dy) SEP    \
  ELTT(float,dz) SEP     /* current change in position */    \
  ELTT(float,odx) SEP    \
  ELTT(float,ody) SEP    \
  ELTT(float,odz) SEP  /* last change of position (for momentum) */    \
  ELTT(float,tdx) SEP    \
  ELTT(float,tdy) SEP    \
  ELTT(float,tdz) SEP  /* temporary storage for averaging gradient */    \
  ELTT(float,curv) SEP            /* curr curvature */    \
  ELTT(float,curvbak) SEP    \
  ELTT(float,val) SEP             /* scalar data value (file: rh.val, sig2-rh.w) */    \
  ELTT(float,imag_val) SEP       /* imaginary part of complex data value */    \
  ELTT(float,cx) SEP    \
  ELTT(float,cy) SEP    \
  ELTT(float,cz) SEP     /* coordinates in canonical coordinate system */    \
  ELTT(float,tx) SEP    \
  ELTT(float,ty) SEP    \
  ELTT(float,tz) SEP     /* tmp coordinate storage */    \
  ELTT(float,tx2) SEP    \
  ELTT(float,ty2) SEP    \
  ELTT(float,tz2) SEP  /* tmp coordinate storage */    \
  ELTT(float,targx) SEP    \
  ELTT(float,targy) SEP    \
  ELTT(float,targz) SEP   /* target coordinates */    \
  ELTT(float,pialx) SEP    \
  ELTT(float,pialy) SEP    \
  ELTT(float,pialz) SEP   /* pial surface coordinates */    \
  ELTT(float,whitex) SEP    \
  ELTT(float,whitey) SEP    \
  ELTT(float,whitez) SEP/* white surface coordinates */    \
  ELTT(float,l4x) SEP    \
  ELTT(float,l4y) SEP    \
  ELTT(float,l4z) SEP   /* layerIV surface coordinates */    \
  ELTT(float,infx) SEP    \
  ELTT(float,infy) SEP    \
  ELTT(float,infz) SEP /* inflated coordinates */    \
  ELTT(float,fx) SEP    \
  ELTT(float,fy) SEP    \
  ELTT(float,fz) SEP      /* flattened coordinates */    \
  ELTT(int,px) SEP    \
  ELTT(int,qx) SEP    \
  ELTT(int,py) SEP    \
  ELTT(int,qy) SEP    \
  ELTT(int,pz) SEP    \
  ELTT(int,qz) SEP /* rational coordinates for exact calculations */    \
  ELTT(float,e1x) SEP    \
  ELTT(float,e1y) SEP    \
  ELTT(float,e1z) SEP  /* 1st basis vector for the local tangent plane */    \
  ELTT(float,e2x) SEP    \
  ELTT(float,e2y) SEP    \
  ELTT(float,e2z) SEP  /* 2nd basis vector for the local tangent plane */    \
  ELTT(float,pe1x) SEP    \
  ELTT(float,pe1y) SEP    \
  ELTT(float,pe1z) SEP  /* 1st basis vector for the local tangent plane */    \
  ELTT(float,pe2x) SEP    \
  ELTT(float,pe2y) SEP    \
  ELTT(float,pe2z) SEP  /* 2nd basis vector for the local tangent plane */    \
  // end of macro

#if 0
#define LIST_OF_VERTEX_ELTS_2 \
  float dipx;    \
  float ipy;    \
  float ipz;  /* dipole position */    \
  float dipnx;    \
  float ipny;    \
  float ipnz; /* dipole orientation */    \
  // end of macro
#else
#define LIST_OF_VERTEX_ELTS_2
#endif

#define LIST_OF_VERTEX_ELTS_3   \
  ELTT(float,nc) SEP              /* curr length normal comp */    \
  ELTT(float,val2) SEP            /* complex comp data value (file: sig3-rh.w) */    \
  ELTT(float,valbak) SEP          /* scalar data stack */    \
  ELTT(float,val2bak) SEP         /* complex comp data stack */    \
  ELTT(float,stat) SEP            /* statistic */    \
    \
  ELTT(int,undefval) SEP          /* [previously dist=0] */    \
  ELTT(int,old_undefval) SEP      /* for smooth_val_sparse */    \
  ELTT(int,fixedval) SEP          /* [previously val=0] */    \
    \
  ELTT(float,fieldsign) SEP       /* fieldsign--final: -1,0,1 (file: rh.fs) */    \
  ELTT(float,fsmask) SEP          /* significance mask (file: rh.fm) */    \
  ELTT(float,d) SEP              /* for distance calculations */    \
  // end of macro
  
#if 0
#define LIST_OF_VERTEX_ELTS_4 \
  ELTP(float,tri_area) SEP      /* array of triangle areas - num long */    \
  ELTP(float,orig_tri_area) SEP /* array of original triangle areas - num long */    \
  ELTT(float,dist) SEP            /* dist from sampled point [defunct: or 1-cos(a)] */    \
  ELTT(float,ox) SEP    \
  ELTT(float,y) SEP    \
  ELTT(float,z) SEP        /* last position (for undoing time steps) */    \
  ELTT(float,mx) SEP    \
  ELTT(float,y) SEP    \
  ELTT(float,z) SEP        /* last movement */    \
  ELTT(float,onc) SEP             /* last length normal comp */    \
  ELTT(float,oval) SEP            /* last scalar data (for smooth_val) */    \
  ELTP(float,fnx) SEP           /* face normal - x component */    \
  ELTP(float,fny) SEP           /* face normal - y component */    \
  ELTP(float,fnz) SEP           /* face normal - z component */    \
  ELTT(float,bnx) SEP    \
  ELTT(float,ny) SEP    \
  ELTT(float,bnx) SEP    \
  ELTT(float,bny) SEP /* boundary normal */    \
  ELTP(float,tri_angle) SEP     /* angles of each triangle this vertex belongs to */    \
  ELTP(float,orig_tri_angle) SEP/* original values of above */    \
  ELTT(float,stress) SEP          /* explosion */    \
  ELTT(float,logshear) SEP    \
  ELTT(float,hearx) SEP    \
  ELTT(float,heary) SEP    \
  ELTT(float,shearx) SEP    \
  ELTT(float,sheary) SEP  /* for shear term */    \
  ELTT(float,ftmp) SEP           /* temporary floating pt. storage */    \
  ELTT(float,logarat) SEP    \
  ELTT(float,logarat) SEP    \
  ELTT(float,qrtarat) SEP /* for area term */    \
  ELTT(float,smx) SEP    \
  ELTT(float,my) SEP    \
  ELTT(float,mz) SEP    \
  ELTT(float,smx) SEP    \
  ELTT(float,smy) SEP    \
  ELTT(float,smz) SEP/* smoothed curr,last move */    \
  // end of macro
#else
#define LIST_OF_VERTEX_ELTS_4
#endif

#define LIST_OF_VERTEX_ELTS_5   \
  ELTT(int,annotation) SEP      /* area label (defunct--now from label file name!) */    \
  ELTT(char,oripflag) SEP    \
  ELTT(char,origripflag) SEP     /* cuts flags */    \
  // end of macro

#if 0
#define LIST_OF_VERTEX_ELTS_6
  float coords[3];
#else
#define LIST_OF_VERTEX_ELTS_6

#endif

#define LIST_OF_VERTEX_ELTS_7    \
  ELTP(void,vp) SEP                     /* to store user's information */    \
  ELTT(float,theta) SEP    \
  ELTT(float,phi) SEP               /* parameterization */    \
  ELTT(float,area) SEP    \
  ELTT(float,origarea) SEP    \
  ELTT(float,group_avg_area) SEP    \
  ELTT(float,K) SEP                 /* Gaussian curvature */    \
  ELTT(float,H) SEP                 /* mean curvature */    \
  ELTT(float,k1) SEP    \
  ELTT(float,k2) SEP                    /* the principal curvatures */    \
  ELTT(float,mean) SEP    \
  ELTT(float,mean_imag) SEP         /* imaginary part of complex statistic */    \
  ELTT(float,std_error) SEP    \
  ELTT(unsigned int,flags) SEP    \
  ELTT(int,fno) SEP                 /* face that this vertex is in */    \
  ELTT(int,cropped) SEP \
  ELTT(short,marked) SEP            /* for a variety of uses */    \
  ELTT(short,marked2) SEP    \
  ELTT(short,marked3) SEP    \
  ELTT(char,neg) SEP                /* 1 if the normal vector is inverted */    \
  ELTT(char,border) SEP             /* flag */    \
  ELTT(char,ripflag)                /* vertex no longer exists - placed last to load the next vertex into cache */ \
  // end of macro
  
#define LIST_OF_VERTEX_ELTS \
  LIST_OF_VERTEX_ELTS_1   \
  LIST_OF_VERTEX_ELTS_2   \
  LIST_OF_VERTEX_ELTS_3   \
  LIST_OF_VERTEX_ELTS_4   \
  LIST_OF_VERTEX_ELTS_5   \
  LIST_OF_VERTEX_ELTS_6   \
  LIST_OF_VERTEX_ELTS_7   \
  // end of macro

#define SEP
#define ELTX(TYPE,NAME) TYPE NAME ;
#define ELTT(TYPE,NAME) TYPE NAME ;
#define ELTP(TARGET,NAME) TARGET *NAME ;
  LIST_OF_VERTEX_ELTS
#undef ELTP
#undef ELTT
#undef ELTX
#undef SEP

  vertex_type_() : dist(nullptr), dist_orig(nullptr), x(0), y(0), z(0), origx(0), origy(0), origz(0) {}

};


#ifndef SEPARATE_VERTEX_TOPOLOGY
typedef vertex_type VERTEX_TOPOLOGY; 
#endif


struct MRIS
{
// The LIST_OF_MRIS_ELTS macro used here enables the the mris_hash
// and other algorithms to process all the elements without having to explicitly name them there and here
//
#define LIST_OF_MRIS_ELTS_1     \
    \
  ELTT(const int,nverticesFrozen) SEP           /* # of vertices on surface is frozen */                                                    \
  ELTT(const int,nvertices) SEP                 /* # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al */         \
  ELTT(const int,nfaces) SEP                    /* # of faces on surface,    change by calling MRISreallocVerticesAndFaces et al */         \
  ELTT(const bool,faceAttachmentDeferred) SEP   /* defer connecting faces to vertices, for performance reasons                   */         \
  ELTT(int,nedges) SEP                          /* # of edges on surface*/    \
  ELTT(int,nstrips) SEP    \
  ELTP(VERTEX_TOPOLOGY,vertices_topology) SEP    \
  ELTP(VERTEX,vertices) SEP    \
  ELTP(void*,dist_storage) SEP                  /* the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored */ \
  ELTP(void*,dist_orig_storage) SEP             /* the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored */ \
  ELTT(const int,tempsAssigned) SEP             /* State of various temp fields that can be borrowed if not already in use   */    \
  ELTP(FACE,faces) SEP    \
  ELTP(MRI_EDGE,edges) SEP    \
  ELTP(FaceNormCacheEntry,faceNormCacheEntries) SEP \
  ELTP(FaceNormDeferredEntry,faceNormDeferredEntries) SEP \
  ELTP(STRIP,strips) SEP    \
  ELTT(float,xctr) SEP    \
  ELTT(float,yctr) SEP    \
  ELTT(float,zctr) SEP    \
  ELTT(float,xlo) SEP    \
  ELTT(float,ylo) SEP    \
  ELTT(float,zlo) SEP    \
  ELTT(float,xhi) SEP    \
  ELTT(float,yhi) SEP    \
  ELTT(float,zhi) SEP    \
  ELTT(float,x0) SEP             /* center of spherical expansion */    \
  ELTT(float,y0) SEP    \
  ELTT(float,z0) SEP    \
  ELTP(VERTEX,v_temporal_pole) SEP    \
  ELTP(VERTEX,v_frontal_pole) SEP    \
  ELTP(VERTEX,v_occipital_pole) SEP    \
  ELTT(float,max_curv) SEP    \
  ELTT(float,min_curv) SEP    \
  ELTT(float,total_area) SEP    \
  ELTT(double,avg_vertex_area) SEP    \
  ELTT(const double,avg_vertex_dist) SEP  /* set by MRIScomputeAvgInterVertexDist */ \
  ELTT(double,std_vertex_dist) SEP    \
  ELTT(float,orig_area) SEP    \
  ELTT(float,neg_area) SEP    \
  ELTT(float,neg_orig_area) SEP   /* amount of original surface in folds */    \
  ELTT(int,zeros) SEP    \
  ELTT(int,hemisphere) SEP      /* which hemisphere */    \
  ELTT(int,initialized) SEP \
  // end of macro

#if 0

#define LIST_OF_MRIS_ELTS_2 \
  ELTT(General_transform,transform) SEP   /* the next two are from this struct (MNI transform) */    \
  ELTP(Transform,linear_transform) SEP    \
  ELTP(Transform,inverse_linear_transform) SEP \
  // end of macro

#else

#define LIST_OF_MRIS_ELTS_2 \
  // end of macro

#endif

#define LIST_OF_MRIS_ELTS_3     \
  ELTP(LTA,lta) SEP    \
  ELTP(MATRIX,SRASToTalSRAS_) SEP    \
  ELTP(MATRIX,TalSRASToSRAS_) SEP    \
  ELTT(int,free_transform) SEP    \
  ELTT(double,radius) SEP           /* radius (if status==MRIS_SPHERE) */    \
  ELTT(float,a) SEP    \
  ELTT(float,b) SEP    \
  ELTT(float,c) SEP                 /* ellipsoid parameters */    \
  ELTT(MRIS_fname_t,fname) SEP      /* file it was originally loaded from */    \
  ELTT(float,Hmin) SEP              /* min mean curvature */    \
  ELTT(float,Hmax) SEP              /* max mean curvature */    \
  ELTT(float,Kmin) SEP              /* min Gaussian curvature */    \
  ELTT(float,Kmax) SEP              /* max Gaussian curvature */    \
  ELTT(double,Ktotal) SEP           /* total Gaussian curvature */    \
  ELTT(MRIS_Status,status) SEP          /* type of surface (e.g. sphere, plane) */    \
  ELTT(MRIS_Status,origxyz_status) SEP  /* type of surface (e.g. sphere, plane) that this origxyz were obtained from */    \
  ELTT(int,patch) SEP               /* if a patch of the surface */    \
  ELTT(int,nlabels) SEP    \
  ELTP(MRIS_AREA_LABEL,labels) SEP  /* nlabels of these (may be null) */    \
  \
  ELTT(char,nsize) SEP              /* size of neighborhoods, or -1 */    \
  ELTT(uchar,vtotalsMightBeTooBig) SEP /* MRISsampleDistances sets this */ \
  ELTX(short,nsizeMaxClock) SEP     /* changed whenever an edge is added or removed, which invalidates the vertex v#num values */ \
  ELTT(char,max_nsize) SEP          /* max the neighborhood size has been set to (typically 3) */    \
  ELTT(char,dist_nsize) SEP         /* max mrisComputeVertexDistances has computed distances out to */ \
  ELTT(char,dist_orig_nsize) SEP    /* max mrisComputeOriginalVertexDistances has computed distances out to */ \
  ELTT(char,dist_alloced_flags) SEP /* two flags, set when any dist(1) or dist_orig(2) allocated */ \
  \
  ELTT(float,avg_nbrs) SEP          /* mean # of vertex neighbors */    \
  ELTP(void,vp) SEP                 /* for misc. use */    \
  ELTT(float,alpha) SEP             /* rotation around z-axis */    \
  ELTT(float,beta) SEP             /* rotation around y-axis */    \
  ELTT(float,gamma) SEP            /* rotation around x-axis */    \
  ELTT(float,da) SEP    \
  ELTT(float,db) SEP    \
  ELTT(float,dg) SEP                /* old deltas */    \
  ELTT(int,type) SEP                /* what type of surface was this initially*/    \
  ELTT(const int,max_vertices) SEP  /* may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces */    \
  ELTT(const int,max_faces) SEP     /* may be bigger than nfaces, set by calling MRISreallocVerticesAndFaces */    \
  ELTT(MRIS_subject_name_t,subject_name) SEP /* name of the subject */    \
  ELTT(float,canon_area) SEP    \
  ELTT(int,noscale) SEP          /* don't scale by surface area if true */    \
  ELTP(float,dx2) SEP             /* an extra set of gradient (not always alloced) */    \
  ELTP(float,dy2) SEP    \
  ELTP(float,dz2) SEP    \
  ELTP(COLOR_TABLE,ct) SEP    \
  ELTT(int,useRealRAS) SEP        /* if 0, vertex position is a conformed volume RAS with c_(r,a,s)=0       */    \
                                  /* if 1, vertex position is a real RAS (volume stored RAS)                */    \
                                  /* The default is 0.                                                      */    \
  ELTT(VOL_GEOM,vg) SEP           /* volume info from which this surface is created. valid iff vg.valid = 1 */    \
  ELTX(MRIS_cmdlines_t, cmdlines) SEP    \
  ELTT(int,ncmds) SEP    \
  ELTT(float,group_avg_surface_area) SEP    /* average of total surface area for group */       \
  ELTT(int,group_avg_vtxarea_loaded) SEP    /* average vertex area for group at each vertex */  \
  ELTT(int,triangle_links_removed) SEP      /* for quad surfaces                         */     \
  ELTP(void,user_parms) SEP                 /* for whatever the user wants to hang here  */     \
  ELTP(MATRIX,m_sras2vox) SEP               /* for converting surface ras to voxel       */     \
  ELTP(MRI,mri_sras2vox) SEP                /* volume that the above matrix is for       */     \
  ELTP(void,mht) SEP \
  ELTP(void,temps)  \
  // end of macro
  
#define LIST_OF_MRIS_ELTS       \
    LIST_OF_MRIS_ELTS_1         \
    LIST_OF_MRIS_ELTS_2         \
    LIST_OF_MRIS_ELTS_3         \
    // end of macro

#define SEP ;
#define ELTP(TARGET, MBR)   TARGET *MBR     // pointers 
#define ELTT(TYPE,   MBR)   TYPE    MBR     // other members that should     be included in the hash
#define ELTX(TYPE,   MBR)   TYPE    MBR     // other members that should NOT be included in the hash
LIST_OF_MRIS_ELTS ;
#undef ELTX
#undef ELTT
#undef ELTP
#undef SEP

};

