
#pragma once
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#include "mrisurf_aaa.h"
#define SEPARATE_VERTEX_TOPOLOGY
struct face_type_ {
      vertices_per_face_t v          ;
    float                 area       ;
    angles_per_triangle_t angle      ;
    angles_per_triangle_t orig_angle ;
    char                  ripflag    ;
    char                  oripflag   ;
    int                   marked     ;
    PDMATRIX              norm       ;
    A3PDMATRIX            gradNorm   ;
};		// face_type_

struct VERTEX_TOPOLOGY {
    //  put the pointers before the ints, before the shorts, before uchars, to reduce size
    //  the whole fits in much less than one cache line, so further ordering is no use
    pSeveralInt   f             ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
    pSeveralUchar n             ;  // size() is num.    array[v->num] the face.v[*] index for this vertex        
    pSeveralInt   e             ;  //  edge state for neighboring vertices                      
    pSeveralInt   v             ;  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
    short         vnum          ;  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
    short         v2num         ;  //  number of 1, or 2-hop neighbors                          
    short         v3num         ;  //  number of 1,2,or 3-hop neighbors                         
    short         vtotal        ;  //  total # of neighbors. copy of vnum.nsizeCur              
    short         nsizeMaxClock ;  //  copy of mris->nsizeMaxClock when v#num                   
    uchar         nsizeMax      ;  //  the max nsize that was used to fill in vnum etc          
    uchar         nsizeCur      ;  //  index of the current v#num in vtotal                     
    uchar         num           ;  //  number of neighboring faces                              
};		// VERTEX_TOPOLOGY

struct vertex_type_ {
    //  managed by MRISfreeDists[_orig] and MRISmakeDists[_orig]
    pSeveralFloat dist               ;  // size() is vtotal.    distance to neighboring vertices based on  xyz   
    pSeveralFloat dist_orig          ;  // size() is vtotal.    distance to neighboring vertices based on origxyz
    int           dist_capacity      ;  //  -- should contain at least vtx_vtotal elements   
    int           dist_orig_capacity ;  //  -- should contain at least vtx_vtotal elements   
    float         x                  ;  //  current coordinates	
    float         y                  ;  //  use MRISsetXYZ() to set
    float         z                  ;
    float         origx              ;  //  original coordinates, see also MRIS::origxyz_status
    float         origy              ;  //  use MRISsetOriginalXYZ(, 
    float         origz              ;  //  or MRISsetOriginalXYZfromXYZ to set
    float         nx                 ;
    float         ny                 ;
    float         nz                 ;  //  curr normal
    float         pnx                ;
    float         pny                ;
    float         pnz                ;  //  pial normal
    float         wnx                ;
    float         wny                ;
    float         wnz                ;  //  white normal
    float         onx                ;
    float         ony                ;
    float         onz                ;  //  original normal
    float         dx                 ;
    float         dy                 ;
    float         dz                 ;  //  current change in position
    float         odx                ;
    float         ody                ;
    float         odz                ;  //  last change of position (for momentum, 
    float         tdx                ;
    float         tdy                ;
    float         tdz                ;  //  temporary storage for averaging gradient
    float         curv               ;  //  curr curvature
    float         curvbak            ;
    float         val                ;  //  scalar data value (file: rh.val, sig2-rh.w)
    float         imag_val           ;  //  imaginary part of complex data value
    float         cx                 ;
    float         cy                 ;
    float         cz                 ;  //  coordinates in canonical coordinate system
    float         tx                 ;
    float         ty                 ;
    float         tz                 ;  //  tmp coordinate storage
    float         t2x                ;
    float         t2y                ;
    float         t2z                ;  //  another tmp coordinate storage
    float         targx              ;
    float         targy              ;
    float         targz              ;  //  target coordinates
    float         pialx              ;
    float         pialy              ;
    float         pialz              ;  //  pial surface coordinates
    float         whitex             ;
    float         whitey             ;
    float         whitez             ;  //  white surface coordinates
    float         l4x                ;
    float         l4y                ;
    float         l4z                ;  //  layerIV surface coordinates
    float         infx               ;
    float         infy               ;
    float         infz               ;  //  inflated coordinates
    float         fx                 ;
    float         fy                 ;
    float         fz                 ;  //  flattened coordinates
    int           px                 ;
    int           qx                 ;
    int           py                 ;
    int           qy                 ;
    int           pz                 ;
    int           qz                 ;  //  rational coordinates for exact calculations
    float         e1x                ;
    float         e1y                ;
    float         e1z                ;  //  1st basis vector for the local tangent plane
    float         e2x                ;
    float         e2y                ;
    float         e2z                ;  //  2nd basis vector for the local tangent plane
    float         pe1x               ;
    float         pe1y               ;
    float         pe1z               ;  //  1st basis vector for the local tangent plane
    float         pe2x               ;
    float         pe2y               ;
    float         pe2z               ;  //  2nd basis vector for the local tangent plane
    float         nc                 ;  //  curr length normal comp 
    float         val2               ;  //  complex comp data value (file: sig3-rh.w) 
    float         valbak             ;  //  scalar data stack 
    float         val2bak            ;  //  complex comp data stack 
    float         stat               ;  //  statistic 
    int           undefval           ;  //  [previously dist=0] 
    int           old_undefval       ;  //  for smooth_val_sparse 
    int           fixedval           ;  //  [previously val=0] 
    float         fieldsign          ;  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
    float         fsmask             ;  //  significance mask (file: rh.fm) 
    float         d                  ;  //  for distance calculations 
    int           annotation         ;  //  area label (defunct--now from label file name!) 
    char          oripflag           ;
    char          origripflag        ;  //  cuts flags 
    p_void        vp                 ;  //  to store user's information 
    float         theta              ;
    float         phi                ;  //  parameterization 
    float         area               ;
    float         origarea           ;
    float         group_avg_area     ;
    float         K                  ;  //  Gaussian curvature 
    float         H                  ;  //  mean curvature 
    float         k1                 ;
    float         k2                 ;  //  the principal curvatures 
    float         mean               ;
    float         mean_imag          ;  //  imaginary part of complex statistic 
    float         std_error          ;
    uint          flags              ;
    int           fno                ;  //  face that this vertex is in 
    int           cropped            ;
    short         marked             ;  //  for a variety of uses 
    short         marked2            ;
    short         marked3            ;
    char          neg                ;  //  1 if the normal vector is inverted 
    char          border             ;  //  flag 
    char          ripflag            ;  //  vertex no longer exists - placed last to load the next vertex into cache
};		// vertex_type_

struct MRIS {
    //  Fields being maintained by specialist functions
    int                           nverticesFrozen          ;  //  # of vertices on surface is frozen
    int                           nvertices                ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
    int                           nfaces                   ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
    bool                          faceAttachmentDeferred   ;  //  defer connecting faces to vertices for performance reasons
    int                           nedges                   ;  //  # of edges on surface
    int                           ncorners                 ;  //  # of triangle corners
    int                           nstrips                  ;
    pSeveralVERTEX_TOPOLOGY       vertices_topology        ;
    pSeveralVERTEX                vertices                 ;
    p_p_void                      dist_storage             ;  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
    p_p_void                      dist_orig_storage        ;  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
    int                           tempsAssigned            ;  //  State of various temp fields that can be borrowed if not already in use
    pSeveralFACE                  faces                    ;
    pSeveralMRI_EDGE              edges                    ;
    pSeveralMRI_CORNER            corners                  ;
    pSeveralFaceNormCacheEntry    faceNormCacheEntries     ;
    pSeveralFaceNormDeferredEntry faceNormDeferredEntries  ;
    pSeveralSTRIP                 strips                   ;
    float                         xctr                     ;
    float                         yctr                     ;
    float                         zctr                     ;
    float                         xlo                      ;
    float                         ylo                      ;
    float                         zlo                      ;
    float                         xhi                      ;
    float                         yhi                      ;
    float                         zhi                      ;
    float                         x0                       ;  //  center of spherical expansion
    float                         y0                       ;
    float                         z0                       ;
    //  v_temporal_pole, v_frontal_pole, and v_occipital_pole don't appear to be used, and are unusual being pointers to vertices
    PVERTEX                       v_temporal_pole          ;
    PVERTEX                       v_frontal_pole           ;
    PVERTEX                       v_occipital_pole         ;
    float                         max_curv                 ;
    float                         min_curv                 ;
    float                         total_area               ;
    double                        avg_vertex_area          ;
    double                        avg_vertex_dist          ;  //  set by MRIScomputeAvgInterVertexDist
    double                        std_vertex_dist          ;
    float                         orig_area                ;
    float                         neg_area                 ;
    float                         neg_orig_area            ;  //  amount of original surface in folds
    int                           zeros                    ;
    int                           hemisphere               ;  //  which hemisphere
    int                           initialized              ;
    PLTA                          lta                      ;
    PMATRIX                       SRASToTalSRAS_           ;
    PMATRIX                       TalSRASToSRAS_           ;
    int                           free_transform           ;
    double                        radius                   ;  //  radius (if status==MRIS_SPHERE)
    float                         a                        ;
    float                         b                        ;
    float                         c                        ;  //  ellipsoid parameters
    MRIS_fname_t                  fname                    ;  //  file it was originally loaded from
    float                         Hmin                     ;  //  min mean curvature
    float                         Hmax                     ;  //  max mean curvature
    float                         Kmin                     ;  //  min Gaussian curvature
    float                         Kmax                     ;  //  max Gaussian curvature
    double                        Ktotal                   ;  //  total Gaussian curvature
    MRIS_Status                   status                   ;  //  type of surface (e.g. sphere, plane)
    MRIS_Status                   origxyz_status           ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
    int                           patch                    ;  //  if a patch of the surface
    int                           nlabels                  ;
    PMRIS_AREA_LABEL              labels                   ;  //  nlabels of these (may be null)
    char                          nsize                    ;  //  size of neighborhoods or -1
    uchar                         vtotalsMightBeTooBig     ;  //  MRISsampleDistances sets this
    short                         nsizeMaxClock            ;  //  changed whenever an edge is added or removed, which invalidates the vertex v#num values
    char                          max_nsize                ;  //  max the neighborhood size has been set to (typically 3)
    char                          dist_nsize               ;  //  max mrisComputeVertexDistances has computed distances out to
    char                          dist_orig_nsize          ;  //  max mrisComputeOriginalVertexDistances has computed distances out to
    char                          dist_alloced_flags       ;  //  two flags, set when any dist(1) or dist_orig(2) allocated
    float                         avg_nbrs                 ;  //  mean # of vertex neighbors
    p_void                        vp                       ;  //  for misc. use
    float                         alpha                    ;  //  rotation around z-axis
    float                         beta                     ;  //  rotation around y-axis
    float                         gamma                    ;  //  rotation around x-axis
    float                         da                       ;
    float                         db                       ;
    float                         dg                       ;  //  old deltas
    int                           type                     ;  //  what type of surface was this initially
    int                           max_vertices             ;  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
    int                           max_faces                ;  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
    MRIS_subject_name_t           subject_name             ;  //  name of the subject
    float                         canon_area               ;
    int                           noscale                  ;  //  don't scale by surface area if true
    pSeveralFloat                 dx2                      ;  //  an extra set of gradient (not always alloced)
    pSeveralFloat                 dy2                      ;
    pSeveralFloat                 dz2                      ;
    PCOLOR_TABLE                  ct                       ;
  int                             orig_xyzspace            ;  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
    int                           useRealRAS               ;  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
    VOL_GEOM                      vg                       ;  //  volume info from which this surface is created. valid iff vg.valid = 1
    MRIS_cmdlines_t               cmdlines                 ;
    int                           ncmds                    ;
    float                         group_avg_surface_area   ;  //  average of total surface area for group
    int                           group_avg_vtxarea_loaded ;  //  average vertex area for group at each vertex
    int                           triangle_links_removed   ;  //  for quad surfaces
    p_void                        user_parms               ;  //  for whatever the user wants to hang here
    PMATRIX                       m_sras2vox               ;  //  for converting surface ras to voxel
    PMRI                          mri_sras2vox             ;  //  volume that the above matrix is for
    p_void                        mht                      ;
    p_void                        temps                    ;
};		// MRIS

#define LIST_OF_FACE_ELTS \
    ELTT(vertices_per_face_t,v)  SEP \
    ELTT(float,area)  SEP \
    ELTT(angles_per_triangle_t,angle)  SEP \
    ELTT(angles_per_triangle_t,orig_angle)  SEP \
    ELTT(char,ripflag)  SEP \
    ELTT(char,oripflag)  SEP \
    ELTT(int,marked)  SEP \
    ELTP(DMATRIX,norm)  SEP \
    ELTX(A3PDMATRIX,gradNorm)  \
// end of macro

#define LIST_OF_VERTEX_TOPOLOGY_ELTS \
    ELTP(int,f)  SEP \
    ELTP(uchar,n)  SEP \
    ELTP(int,e)  SEP \
    ELTP(int,v)  SEP \
    ELTT(short,vnum)  SEP \
    ELTT(short,v2num)  SEP \
    ELTT(short,v3num)  SEP \
    ELTT(short,vtotal)  SEP \
    ELTX(short,nsizeMaxClock)  SEP \
    ELTT(uchar,nsizeMax)  SEP \
    ELTT(uchar,nsizeCur)  SEP \
    ELTT(uchar,num)  \
// end of macro

#define LIST_OF_VERTEX_ELTS_1 \
    ELTP(float,dist)  SEP \
    ELTP(float,dist_orig)  SEP \
    ELTT(int,dist_capacity)  SEP \
    ELTT(int,dist_orig_capacity)  SEP \
    ELTT(float,x)  SEP \
    ELTT(float,y)  SEP \
    ELTT(float,z)  SEP \
    ELTT(float,origx)  SEP \
    ELTT(float,origy)  SEP \
    ELTT(float,origz)  SEP \
    ELTT(float,nx)  SEP \
    ELTT(float,ny)  SEP \
    ELTT(float,nz)  SEP \
    ELTT(float,pnx)  SEP \
    ELTT(float,pny)  SEP \
    ELTT(float,pnz)  SEP \
    ELTT(float,wnx)  SEP \
    ELTT(float,wny)  SEP \
    ELTT(float,wnz)  SEP \
    ELTT(float,onx)  SEP \
    ELTT(float,ony)  SEP \
    ELTT(float,onz)  SEP \
    ELTT(float,dx)  SEP \
    ELTT(float,dy)  SEP \
    ELTT(float,dz)  SEP \
    ELTT(float,odx)  SEP \
    ELTT(float,ody)  SEP \
    ELTT(float,odz)  SEP \
    ELTT(float,tdx)  SEP \
    ELTT(float,tdy)  SEP \
    ELTT(float,tdz)  SEP \
    ELTT(float,curv)  SEP \
    ELTT(float,curvbak)  SEP \
    ELTT(float,val)  SEP \
    ELTT(float,imag_val)  SEP \
    ELTT(float,cx)  SEP \
    ELTT(float,cy)  SEP \
    ELTT(float,cz)  SEP \
    ELTT(float,tx)  SEP \
    ELTT(float,ty)  SEP \
    ELTT(float,tz)  SEP \
    ELTT(float,t2x)  SEP \
    ELTT(float,t2y)  SEP \
    ELTT(float,t2z)  SEP \
    ELTT(float,targx)  SEP \
    ELTT(float,targy)  SEP \
    ELTT(float,targz)  SEP \
    ELTT(float,pialx)  SEP \
    ELTT(float,pialy)  SEP \
    ELTT(float,pialz)  SEP \
    ELTT(float,whitex)  SEP \
    ELTT(float,whitey)  SEP \
    ELTT(float,whitez)  SEP \
    ELTT(float,l4x)  SEP \
    ELTT(float,l4y)  SEP \
    ELTT(float,l4z)  SEP \
    ELTT(float,infx)  SEP \
    ELTT(float,infy)  SEP \
    ELTT(float,infz)  SEP \
    ELTT(float,fx)  SEP \
    ELTT(float,fy)  SEP \
    ELTT(float,fz)  SEP \
    ELTT(int,px)  SEP \
    ELTT(int,qx)  SEP \
    ELTT(int,py)  SEP \
    ELTT(int,qy)  SEP \
    ELTT(int,pz)  SEP \
    ELTT(int,qz)  SEP \
    ELTT(float,e1x)  SEP \
    ELTT(float,e1y)  SEP \
    ELTT(float,e1z)  SEP \
    ELTT(float,e2x)  SEP \
    ELTT(float,e2y)  SEP \
    ELTT(float,e2z)  SEP \
    ELTT(float,pe1x)  SEP \
    ELTT(float,pe1y)  SEP \
    ELTT(float,pe1z)  SEP \
    ELTT(float,pe2x)  SEP \
    ELTT(float,pe2y)  SEP \
    ELTT(float,pe2z)  SEP \
// end of macro

#define LIST_OF_VERTEX_ELTS_3 \
    ELTT(float,nc)  SEP \
    ELTT(float,val2)  SEP \
    ELTT(float,valbak)  SEP \
    ELTT(float,val2bak)  SEP \
    ELTT(float,stat)  SEP \
    ELTT(int,undefval)  SEP \
    ELTT(int,old_undefval)  SEP \
    ELTT(int,fixedval)  SEP \
    ELTT(float,fieldsign)  SEP \
    ELTT(float,fsmask)  SEP \
    ELTT(float,d)  SEP \
// end of macro

#define LIST_OF_VERTEX_ELTS_5 \
    ELTT(int,annotation)  SEP \
    ELTT(char,oripflag)  SEP \
    ELTT(char,origripflag)  SEP \
// end of macro

#define LIST_OF_VERTEX_ELTS_7 \
    ELTX(p_void,vp)  SEP \
    ELTT(float,theta)  SEP \
    ELTT(float,phi)  SEP \
    ELTT(float,area)  SEP \
    ELTT(float,origarea)  SEP \
    ELTT(float,group_avg_area)  SEP \
    ELTT(float,K)  SEP \
    ELTT(float,H)  SEP \
    ELTT(float,k1)  SEP \
    ELTT(float,k2)  SEP \
    ELTT(float,mean)  SEP \
    ELTT(float,mean_imag)  SEP \
    ELTT(float,std_error)  SEP \
    ELTT(uint,flags)  SEP \
    ELTT(int,fno)  SEP \
    ELTT(int,cropped)  SEP \
    ELTT(short,marked)  SEP \
    ELTT(short,marked2)  SEP \
    ELTT(short,marked3)  SEP \
    ELTT(char,neg)  SEP \
    ELTT(char,border)  SEP \
    ELTT(char,ripflag)  SEP \
// end of macro

#define LIST_OF_VERTEX_ELTS \
    LIST_OF_VERTEX_ELTS_1 SEP \
    LIST_OF_VERTEX_ELTS_3 SEP \
    LIST_OF_VERTEX_ELTS_5 SEP \
    LIST_OF_VERTEX_ELTS_7 \
// end of macro

#define LIST_OF_MRIS_ELTS_1 \
    ELTT(int,nverticesFrozen)  SEP \
    ELTT(int,nvertices)  SEP \
    ELTT(int,nfaces)  SEP \
    ELTT(bool,faceAttachmentDeferred)  SEP \
    ELTT(int,nedges)  SEP \
    ELTT(int,ncorners)  SEP \
    ELTT(int,nstrips)  SEP \
    ELTP(VERTEX_TOPOLOGY,vertices_topology)  SEP \
    ELTP(VERTEX,vertices)  SEP \
    ELTX(p_p_void,dist_storage)  SEP \
    ELTX(p_p_void,dist_orig_storage)  SEP \
    ELTT(int,tempsAssigned)  SEP \
    ELTP(FACE,faces)  SEP \
    ELTP(MRI_EDGE,edges)  SEP \
    ELTP(MRI_CORNER,corners)  SEP \
    ELTP(FaceNormCacheEntry,faceNormCacheEntries)  SEP \
    ELTP(FaceNormDeferredEntry,faceNormDeferredEntries)  SEP \
    ELTP(STRIP,strips)  SEP \
    ELTT(float,xctr)  SEP \
    ELTT(float,yctr)  SEP \
    ELTT(float,zctr)  SEP \
    ELTT(float,xlo)  SEP \
    ELTT(float,ylo)  SEP \
    ELTT(float,zlo)  SEP \
    ELTT(float,xhi)  SEP \
    ELTT(float,yhi)  SEP \
    ELTT(float,zhi)  SEP \
    ELTT(float,x0)  SEP \
    ELTT(float,y0)  SEP \
    ELTT(float,z0)  SEP \
    ELTP(VERTEX,v_temporal_pole)  SEP \
    ELTP(VERTEX,v_frontal_pole)  SEP \
    ELTP(VERTEX,v_occipital_pole)  SEP \
    ELTT(float,max_curv)  SEP \
    ELTT(float,min_curv)  SEP \
    ELTT(float,total_area)  SEP \
    ELTT(double,avg_vertex_area)  SEP \
    ELTT(double,avg_vertex_dist)  SEP \
    ELTT(double,std_vertex_dist)  SEP \
    ELTT(float,orig_area)  SEP \
    ELTT(float,neg_area)  SEP \
    ELTT(float,neg_orig_area)  SEP \
    ELTT(int,zeros)  SEP \
    ELTT(int,hemisphere)  SEP \
    ELTT(int,initialized)  SEP \
// end of macro

#define LIST_OF_MRIS_ELTS_3 \
    ELTP(LTA,lta)  SEP \
    ELTP(MATRIX,SRASToTalSRAS_)  SEP \
    ELTP(MATRIX,TalSRASToSRAS_)  SEP \
    ELTT(int,free_transform)  SEP \
    ELTT(double,radius)  SEP \
    ELTT(float,a)  SEP \
    ELTT(float,b)  SEP \
    ELTT(float,c)  SEP \
    ELTT(MRIS_fname_t,fname)  SEP \
    ELTT(float,Hmin)  SEP \
    ELTT(float,Hmax)  SEP \
    ELTT(float,Kmin)  SEP \
    ELTT(float,Kmax)  SEP \
    ELTT(double,Ktotal)  SEP \
    ELTT(MRIS_Status,status)  SEP \
    ELTT(MRIS_Status,origxyz_status)  SEP \
    ELTT(int,patch)  SEP \
    ELTT(int,nlabels)  SEP \
    ELTP(MRIS_AREA_LABEL,labels)  SEP \
    ELTT(char,nsize)  SEP \
    ELTT(uchar,vtotalsMightBeTooBig)  SEP \
    ELTX(short,nsizeMaxClock)  SEP \
    ELTT(char,max_nsize)  SEP \
    ELTT(char,dist_nsize)  SEP \
    ELTT(char,dist_orig_nsize)  SEP \
    ELTT(char,dist_alloced_flags)  SEP \
    ELTT(float,avg_nbrs)  SEP \
    ELTX(p_void,vp)  SEP \
    ELTT(float,alpha)  SEP \
    ELTT(float,beta)  SEP \
    ELTT(float,gamma)  SEP \
    ELTT(float,da)  SEP \
    ELTT(float,db)  SEP \
    ELTT(float,dg)  SEP \
    ELTT(int,type)  SEP \
    ELTT(int,max_vertices)  SEP \
    ELTT(int,max_faces)  SEP \
    ELTT(MRIS_subject_name_t,subject_name)  SEP \
    ELTT(float,canon_area)  SEP \
    ELTT(int,noscale)  SEP \
    ELTP(float,dx2)  SEP \
    ELTP(float,dy2)  SEP \
    ELTP(float,dz2)  SEP \
    ELTP(COLOR_TABLE,ct)  SEP \
    ELTT(int,orig_xyzspace)  SEP \
    ELTT(int,useRealRAS)  SEP \
    ELTT(VOL_GEOM,vg)  SEP \
    ELTX(MRIS_cmdlines_t,cmdlines)  SEP \
    ELTT(int,ncmds)  SEP \
    ELTT(float,group_avg_surface_area)  SEP \
    ELTT(int,group_avg_vtxarea_loaded)  SEP \
    ELTT(int,triangle_links_removed)  SEP \
    ELTX(p_void,user_parms)  SEP \
    ELTP(MATRIX,m_sras2vox)  SEP \
    ELTP(MRI,mri_sras2vox)  SEP \
    ELTX(p_void,mht)  SEP \
    ELTX(p_void,temps)  SEP \
// end of macro

#define LIST_OF_MRIS_ELTS \
    LIST_OF_MRIS_ELTS_1 SEP \
    LIST_OF_MRIS_ELTS_3 \
// end of macro

