#if 0
/*------------------------ STATIC PROTOTYPES -------------------------*/
static void notifyActiveRealmTreesChangedNFacesNVertices(MRIS const * const mris);
int MRIScomputeAllDistances(MRI_SURFACE *mris);
#if 0
static MRI_SP *MRISPiterative_blur(MRI_SURFACE *mris,
                                   MRI_SP *mrisp_source,
                                   MRI_SP *mrisp_dst,
                                   float sigma, int frame) ;
#endif
static double MRISavgInterVertexDist(MRIS *Surf, double *StdDev);
static int mrisReadAsciiCurvatureFile(MRI_SURFACE *mris, const char *fname);
static double mrisComputeSSE_MEF(
    MRI_SURFACE *mris, INTEGRATION_PARMS *parms, MRI *mri30, MRI *mri5, double weight30, double weight5, MHT *mht);
static int mrisMarkIntersections(MRI_SURFACE *mris);
static int mrisAverageSignedGradients(MRI_SURFACE *mris, int num_avgs);
#if 0
static int mrisAverageWeightedGradients(MRI_SURFACE *mris, int num_avgs) ;
#endif
static int mrisComputeDuraTerm(MRI_SURFACE *mris, double l_dura, MRI *mri_dura, double dura_thresh);
static double mrisComputeHistoNegativeLikelihood(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisComputeNegativeLogPosterior(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int *pnvox);
static double mrisComputeNegativeLogPosterior2D(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int *pnvox);
static int mrisComputeHistoTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputePosteriorTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputePosterior2DTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
int MRISrestoreExtraGradients(MRI_SURFACE *mris);
static int mrisComputePositioningGradients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisFindGrayWhiteBorderMean(MRI_SURFACE *mris, MRI *mri);
static int mrisDumpDefectiveEdge(MRI_SURFACE *mris, int vno1, int vno2);
static int mrisMarkBadEdgeVertices(MRI_SURFACE *mris, int mark);
static int mrisCheckSurface(MRI_SURFACE *mris);
#if 0
static int mrisComputeCanonicalBasis(MRI_SURFACE *mris, int fno,
                                     double origin[3],double e0[3],
                                     double e1[3]);
#endif
static int mrisInitializeNeighborhood(MRI_SURFACE *mris, int vno);
static int mrisSetVertexFaceIndex(MRI_SURFACE *mris, int vno, int fno);
static int isFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2);
static int findFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2);
static int mrisAddFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2);
#if (!SPHERE_INTERSECTION)
static int mrisComputeCanonicalEdgeBasis(
    MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2, double origin[3], double e0[3], double e1[3]);
#endif

static int triangleMarked(MRI_SURFACE *mris, int fno);

static int mrisCalculateOriginalFaceCentroid(MRI_SURFACE *mris, int fno, float *px, float *py, float *pz);
static int mrisCalculateFaceCentroid(MRI_SURFACE *mris, int fno, float *px, float *py, float *pz);
static int mrisCalculateCanonicalFaceCentroid(MRI_SURFACE *mris, int fno, float *px, float *py, float *pz);
static int mrisDirectionTriangleIntersection(
    MRI_SURFACE *mris, float x0, float y0, float z0, float nx, float ny, float nz, MHT *mht, double *pdist, int vno);
static int mrisComputeCurvatureMinMax(MRI_SURFACE *mris);
static int mrisAllNormalDirectionCurrentTriangleIntersections(
    MRI_SURFACE *mris, VERTEX *v, MHT *mht, double *pdist, int *flist);
static int load_triangle_vertices(MRI_SURFACE *mris, int fno, double U0[3], double U1[3], double U2[3], int which);
static int load_orig_triangle_vertices(MRI_SURFACE *mris, int fno, double U0[3], double U1[3], double U2[3]);
static void mrisDumpFace(MRI_SURFACE *mris, int fno, FILE *fp);
static int mrisAddEdge(MRI_SURFACE *mris, int vno1, int vno2);

static int mrisComputeSurfaceDimensions(MRI_SURFACE *mris);
// static int   mrisFindNeighbors(MRI_SURFACE *mris) ;
static float mrisNormalize(float v[3]);
static float mrisTriangleArea(MRIS *mris, int fac, int n);
static void mrisFaceAreaNormal(MRIS *mris, int fac, float norm[]);
static int mrisComputeOrigNormal(MRIS *mris, int vno, float norm[]);
static int mrisComputeWhiteNormal(MRIS *mris, int vno, float norm[]);
static int mrisComputeWhichSurfaceRepulsionTerm(
    MRI_SURFACE *mris, double l_repulse, MHT *mht, int which, float max_dot);
static int mrisComputePialNormal(MRIS *mris, int vno, float norm[]);
static int mrisOrigNormalFace(MRIS *mris, int fac, int n, float norm[]);
static int mrisPialNormalFace(MRIS *mris, int fac, int n, float norm[]);
static int mrisWhiteNormalFace(MRIS *mris, int fac, int n, float norm[]);
static int mrisReadTransform(MRIS *mris, const char *mris_fname);
static MRI_SURFACE *mrisReadAsciiFile(const char *fname);
static MRI_SURFACE *mrisReadGeoFile(const char *fname);
static MRI_SURFACE *MRISreadVTK(MRI_SURFACE *mris, const char *fname);
static MRI_SURFACE *mrisReadSTLfile(const char *fname);
static int mrisReadGeoFilePositions(MRI_SURFACE *mris, const char *fname);
static MRI_SURFACE *mrisReadTriangleFile(const char *fname, double pct_over);
static int mrisReadTriangleFilePositions(MRI_SURFACE *mris, const char *fname);
static SMALL_SURFACE *mrisReadTriangleFileVertexPositionsOnly(const char *fname);
/*static int   mrisReadFieldsign(MRI_SURFACE *mris, const char *fname) ;*/
static double mrisComputeNonlinearAreaSSE(MRI_SURFACE *mris);
static double mrisComputeNonlinearDistanceSSE(MRI_SURFACE *mris);
static double mrisComputeSpringEnergy(MRI_SURFACE *mris);
static double mrisComputeLaplacianEnergy(MRI_SURFACE *mris);
static double mrisComputeThicknessSmoothnessEnergy(MRI_SURFACE *mris, double l_repulse, INTEGRATION_PARMS *parms);
static double mrisComputeThicknessMinimizationEnergy(MRI_SURFACE *mris, double l_thick_min, INTEGRATION_PARMS *parms);
static double mrisComputeThicknessNormalEnergy(MRI_SURFACE *mris, double l_thick_normal, INTEGRATION_PARMS *parms);
static double mrisComputeThicknessSpringEnergy(MRI_SURFACE *mris, double l_thick_spring, INTEGRATION_PARMS *parms);
static double mrisComputeThicknessParallelEnergy(MRI_SURFACE *mris, double l_thick_parallel, INTEGRATION_PARMS *parms);
static double mrisComputeAshburnerTriangleEnergy(MRI_SURFACE *mris,
                                                 double l_ashburner_triangle,
                                                 INTEGRATION_PARMS *parms);
static double mrisComputeRepulsiveEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht_v_current, MHT *mht_f_current);
static int mrisComputeRepulsiveTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht_v, MHT *mht_f);
static double mrisComputeRepulsiveRatioEnergy(MRI_SURFACE *mris, double l_repulse);
static int mrisComputeRepulsiveRatioTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht_v);
static int mrisComputeSurfaceRepulsionTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht);
static double mrisComputeSurfaceRepulsionEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht);
static int mrisComputeThicknessSmoothnessTerm(MRI_SURFACE *mris, double l_tsmooth, INTEGRATION_PARMS *parms);
static int mrisComputeThicknessMinimizationTerm(MRI_SURFACE *mris, double l_thick_min, INTEGRATION_PARMS *parms);
static int mrisComputeThicknessNormalTerm(MRI_SURFACE *mris, double l_thick_normal, INTEGRATION_PARMS *parms);
static int mrisComputeThicknessSpringTerm(MRI_SURFACE *mris, double l_thick_spring, INTEGRATION_PARMS *parms);
static int mrisComputeThicknessParallelTerm(MRI_SURFACE *mris, double l_thick_parallel, INTEGRATION_PARMS *parms);
static int mrisComputeAshburnerTriangleTerm(MRI_SURFACE *mris, double l_ashburner_triangle, INTEGRATION_PARMS *parms);
static double mrisComputeNonlinearSpringEnergy(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisComputeTangentialSpringEnergy(MRI_SURFACE *mris);
static double mrisComputeIntensityError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisComputeTargetLocationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisComputeDuraError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);

static int mrisCheckSurfaceNbrs(MRI_SURFACE *mris);
static double mrisComputeIntensityGradientError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisComputeSphereError(MRI_SURFACE *mris, double l_sphere, double a);
static double mrisComputeDistanceError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);

static double mrisComputeCorrelationErrorTraceable(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds, bool trace);
static double mrisComputeCorrelationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds);

static int mrisComputeVertexDistances(MRI_SURFACE *mris);
static int mrisComputeOriginalVertexDistances(MRI_SURFACE *mris);
static double mrisComputeError(MRI_SURFACE *mris,
                               INTEGRATION_PARMS *parms,
                               float *parea_rms,
                               float *pangle_rms,
                               float *pcurv_rms,
                               float *pdist_rms,
                               float *pcorr_rms);
static int mrisIntegrationEpoch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_avgs);
static int mrisRemoveNegativeArea(
    MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_avgs, float min_area_pct, int max_passes);
static double mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisLineMinimizeSearch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisAsynchronousTimeStep(MRI_SURFACE *mris, float momentum, float dt, MHT *mht, float max_mag);
static double mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, float momentum, float dt, MHT *mht, float max_mag);
static double mrisAdaptiveTimeStep(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisOrientEllipsoid(MRI_SURFACE *mris);
static int mrisOrientPlane(MRI_SURFACE *mris);
#if AVERAGE_AREAS
static int mrisAverageAreas(MRI_SURFACE *mris, int num_avgs, int which);
#endif
static int transform(float *xptr, float *yptr, float *zptr, float nx, float ny, float nz, float d);

static int mrisComputeTangentPlanes(MRI_SURFACE *mris);
static int mrisRemoveLink(MRI_SURFACE *mris, int vno1, int vno2);
static int mrisRemoveEdge(MRI_SURFACE *mris, int vno1, int vno2);
static int mrisRemoveFace(MRI_SURFACE *mris, int fno);
static int mrisCountTotalNeighbors(MRI_SURFACE *mris);
static int mrisCountValidLinks(MRI_SURFACE *mris, int vno1, int vno2);
static int mrisComputeSpringTerm(MRI_SURFACE *mris, double l_spring);
static int mrisComputeBorderTerm(MRI_SURFACE *mris, double l_border);
static int mrisComputeMaxSpringTerm(MRI_SURFACE *mris, double l_spring);
static int mrisComputeLaplacianTerm(MRI_SURFACE *mris, double l_laplacian);
static int mrisComputeLinkTerm(MRI_SURFACE *mris, double l_spring, int pial);
static int mrisComputeNormalizedSpringTerm(MRI_SURFACE *mris, double l_spring);
static int mrisComputeIntensityTerm(
    MRI_SURFACE *mris, double l_intensity, MRI *mri_brain, MRI *mri_smooth, double sigma, INTEGRATION_PARMS *parms);
static int mrisComputeTargetLocationTerm(MRI_SURFACE *mris, double l_location, INTEGRATION_PARMS *parms);
static int mrisComputeIntensityGradientTerm(MRI_SURFACE *mris, double l_grad, MRI *mri_brain, MRI *mri_smooth);
static int mrisComputeSphereTerm(MRI_SURFACE *mris, double l_sphere, float radius, int explode_flag);
static int mrisComputeConvexityTerm(MRI_SURFACE *mris, double l_convex);
static int mrisComputeExpansionTerm(MRI_SURFACE *mris, double l_expand);
static int mrisComputeDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputeNonlinearDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputeCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv);
static int mrisComputeSurfaceNormalIntersectionTerm(MRI_SURFACE *mris, MHT *mht, double l_norm, double max_dist);
static double mrisComputeQuadraticCurvatureSSE(MRI_SURFACE *mris, double l_curv);
static int mrisComputePolarCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisComputeVectorCorrelationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds);
static int mrisComputeVectorCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputePolarVectorCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputeAngleAreaTerms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int mrisComputeNonlinearAreaTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static int MRISclearD(MRI_SURFACE *mris);
static int mrisClearExtraGradient(MRI_SURFACE *mris);
static int mrisClearMomentum(MRI_SURFACE *mris);
static int mrisValidFaces(MRI_SURFACE *mris);
static int mrisLabelVertices(MRI_SURFACE *mris, float cx, float cy, float cz, int label, float radius);
static int mrisComputeShrinkwrapTerm(MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap);
static double mrisComputeShrinkwrapError(MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap);
static int mrisComputeExpandwrapTerm(MRI_SURFACE *mris, MRI *mri_brain, double l_expandwrap);
static double mrisComputeExpandwrapError(MRI_SURFACE *mris, MRI *mri_brain, double l_expandwrap, double target_radius);

static int project_point_onto_sphere(float cx, float cy, float cz, float radius, float *pcx, float *pcy, float *pcz);
static int mrisProjectOntoSurface(MRI_SURFACE *mris, int which_vertices);
static int mrisProjectSurface(MRI_SURFACE *mris);
static int mrisOrientSurface(MRI_SURFACE *mris);
static int mrisComputeBoundaryNormals(MRI_SURFACE *mris);
static int mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter);
static int mrisFlipPatch(MRI_SURFACE *mris);

static int mrisPlaceVertexInOrigFace(MRI_SURFACE *mris, VERTEX *v, int fno);
static int vertexInFace(MRI_SURFACE *mris, int vno, int fno);

static int mrisComputeNonlinearSpringTerm(MRI_SURFACE *mris, double l_nlspring, INTEGRATION_PARMS *parms);
static int mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring);
static int mrisComputeNonlinearTangentialSpringTerm(MRI_SURFACE *mris, double l_spring, double min_dist);
static int mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring);
static float minNeighborDistance(MRI_SURFACE *mris);
static bool mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno, 
                                               MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context* ctx);
static int mrisComputeVariableSmoothnessCoefficients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);

static int mrisLogStatus(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, FILE *fp, float dt, float old_sse);
static int mrisWriteSnapshots(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t);
static int mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t);
static int mrisTrackTotalDistance(MRI_SURFACE *mris);
static int mrisTrackTotalDistanceNew(MRI_SURFACE *mris);
static bool mrisLimitGradientDistance(MRI_SURFACE *mris, MHT const *mht, int vno,
                                     MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context* ctx);
static int mrisFillFace(MRI_SURFACE *mris, MRI *mri, int fno);
static int mrisHatchFace(MRI_SURFACE *mris, MRI *mri, int fno, int on);

static double mrisRmsDistanceError(MRI_SURFACE *mris);
static int mrisRemoveVertexLink(MRI_SURFACE *mris, int vno1, int vno2);
static int mrisStoreVtotalInV3num(MRI_SURFACE *mris);
static int mrisFindAllOverlappingFaces(MRI_SURFACE *mris, MHT *mht, int fno, int *flist);

// The following two functions added for processing two channel MEF
static int mrisComputeIntensityTerm_mef(MRI_SURFACE *mris,
                                        double l_intensity,
                                        MRI *mri_30,
                                        MRI *mri_5,
                                        double sigma_global,
                                        float weight30,
                                        float weight5,
                                        INTEGRATION_PARMS *parms);
static double mrisRmsValError_mef(MRI_SURFACE *mris, MRI *mri_30, MRI *mri_5, float weight30, float weight5);

static int mrisDivideEdge(MRI_SURFACE *mris, int vno1, int vno2);
static int mrisDivideFace(MRI_SURFACE *mris, int fno, int vno1, int vno2, int vnew_no);

static MATRIX *getSRASToTalSRAS(LT *lt);

#endif

#if 0
static int mrisSamplePialCoordsInTangentPlane(MRI_SURFACE *mris, VERTEX *v,
    float x, float y, float z, float *pxp,
    float *pyp, float *pzp);
#endif
// static int mrisSoapBubbleIntersectingDefects(MRI_SURFACE *mris);


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisRemoveNormalGradientComponent(MRI_SURFACE *mris, int vno)
{
  VERTEX   *v ;
  float    dot ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
  {
    return(NO_ERROR) ;
  }

  dot = v->nx*v->odx + v->ny*v->ody + v->nz*v->odz ;
  v->odx -= dot*v->nx ;
  v->ody -= dot*v->ny ;
  v->odz -= dot*v->nz ;

  return(NO_ERROR) ;
}
#endif


#if 0
static float computeAngleSign(MRIS *mris, int v,int v2, int v3)
{
  float v0[3],v1[3],d1,d2,d3,dot;

  v0[0] = mris->vertices[v2].x - mris->vertices[v].x;
  v0[1] = mris->vertices[v2].y - mris->vertices[v].y;
  v0[2] = mris->vertices[v2].z - mris->vertices[v].z;
  v1[0] = mris->vertices[v3].x - mris->vertices[v].x;
  v1[1] = mris->vertices[v3].y - mris->vertices[v].y;
  v1[2] = mris->vertices[v3].z - mris->vertices[v].z;
  d1 = -v1[1]*v0[2] + v0[1]*v1[2];
  d2 = v1[0]*v0[2] - v0[0]*v1[2];
  d3 = -v1[0]*v0[1] + v0[0]*v1[1];

  dot =  mris->vertices[v].x * d1
         + mris->vertices[v].y * d2
         + mris->vertices[v].z * d3 ;

  return dot;
}
#endif

#if ADD_EXTRA_VERTICES
/* verify if some defects are not overlapping : to be implemented */
static DEFECT_LIST *mrisDefectAnalysis(MRIS *mris, DEFECT_LIST *dl)
{
  int i, n;
  DEFECT *defect;

  fprintf(stderr, "analyzing defects...\n");
  /* initiliaze flags */
  for (n = 0; n < mris->nvertices; n++) {
    mris->vertices[n].ripflag = 0;
    mris->vertices[n].marked = 0;
    mris->vertices[n].fixedval = 0;
    mris->vertices[n].undefval = 0;
    mris->vertices[n].old_undefval = 0;
  }

  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    // fprintf(stderr,"\n%d  : ",i);
    for (n = 0; n < defect->nvertices; n++) {
      // fprintf(stderr,"%d ",defect->vertices[n]);
      if (mris->vertices[defect->vertices[n]].marked) {
        ErrorExit(ERROR_BADPARM,
                  "mrisDefectAnalysis : defect %d overlap defect %d",
                  i,
                  mris->vertices[defect->vertices[n]].marked - 1);
      }
      mris->vertices[defect->vertices[n]].marked = i + 1;
    }
  }
  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    for (n = 0; n < defect->nborder; n++) {
      if (mris->vertices[defect->border[n]].marked) {
        ErrorExit(ERROR_BADPARM,
                  "mrisDefectAnalysis : defect border %d overlap inside defect %d",
                  i,
                  mris->vertices[defect->border[n]].marked - 1);
      }
    }
  }
  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    for (n = 0; n < defect->nborder; n++) {
      if (mris->vertices[defect->border[n]].marked) {
        ErrorExit(ERROR_BADPARM,
                  "mrisDefectAnalysis : defect border %d overlap border defect %d",
                  i,
                  mris->vertices[defect->border[n]].marked - 1);
      }
      mris->vertices[defect->border[n]].marked = i + 1;
    }
  }

  // fprintf(stderr,"\nend analysis\n");
  return dl;
}
#endif

#if ADD_EXTRA_VERTICES
static int pushApartDefects(MRIS *mris, DEFECT_LIST *dl)
{
  int i, n, p, m, vnop, total_nadded, nadded;
  DEFECT *defect;
  VERTEX *v, *vp, *vadded;

  fprintf(stderr, "separating defects...\n");
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    v->fixedval = 0;
    v->marked = 0;
  }
  for (total_nadded = i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];

    for (nadded = n = 0; n < defect->nborder; n++) {
      v = &mris->vertices[defect->border[n]];
      for (p = 0; p < v->vnum; p++) {
        vnop = v->v[p];
        vp = &mris->vertices[vnop];
        if (vp->fixedval && vp->fixedval != i + 1) {
          /* border vertex from another defect : add one vertex */
          if (mris->nvertices == mris->max_vertices) {
            ErrorExit(ERROR_BADPARM, "pushApartDefect : could not allocate extra vertex\n");
          }
          /* just making sure */
          for (m = 0; m < mris->nvertices; m++) {
            mris->vertices[m].marked = 0;
          }
          mrisDivideEdge(mris, defect->border[n], vnop);
          vadded = &mris->vertices[mris->nvertices - 1];
          // fprintf(stderr,"adding vertex %d(%d)\n",mris->nvertices-1,vn->fixedval);
          vadded->fixedval = 0;

          /* spherical projection */
          sphericalProjection(vadded->cx, vadded->cy, vadded->cz, &vadded->cx, &vadded->cy, &vadded->cz);
          vadded->x = vadded->cx;
          vadded->y = vadded->cy;
          vadded->z = vadded->cz;
          nadded++;
        }
      }
      mris->vertices[defect->border[n]].fixedval = i + 1;
    }
    //          if(nadded) fprintf(stderr,"defect %d : %d vertices have been added to the
    //          surface\n",defect->defect_number,nadded);
    total_nadded += nadded;
  }
  fprintf(stderr, "   total of %d vertices have been added to the surface\n", total_nadded);
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    v->fixedval = 0;
    v->marked = 0;
  }
  return NO_ERROR;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisCountNegativeVertices(MRI_SURFACE *mris)
{
  int     vno, neg ;
  VERTEX  *v ;

  for (neg = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (v->neg)
    {
      neg++ ;
    }
  }

  return(neg) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI_SURFACE *
MRISremoveNegativeVertices(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                           int min_neg, float min_neg_pct)
{
  int   t, niterations, write_iterations, neg, total_vertices, base_avgs ;
  float pct_neg, delta_t, scale, pct_neg_area, l_dist ;

  if (min_neg < 0)
  {
    min_neg = 0 ;
  }
  if (min_neg_pct < 0.0f)
  {
    min_neg_pct = 0.0f ;
  }

  if (Gdiag & DIAG_WRITE && parms->fp == NULL)
  {
    char fname[STRLEN] ;

    sprintf
    (fname, "%s.%s.out",
     mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    parms->fp = fopen(fname, "w") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
                Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    mrisLogIntegrationParms(stderr, mris, parms) ;
  }

  parms->start_t = 0 ;
  mrisProjectSurface(mris) ;
  if (Gdiag & DIAG_WRITE)
  {
    mrisLogStatus(mris, parms, parms->fp, 0.0f, -1) ;
  }
  mrisClearMomentum(mris) ;
  niterations = parms->niterations ;
  write_iterations = parms->write_iterations ;
  if (Gdiag & DIAG_WRITE && write_iterations > 0)
  {
    mrisWriteSnapshot(mris, parms, 0) ;
  }
  total_vertices = MRISvalidVertices(mris) ;
  neg = mrisCountNegativeVertices(mris) ;
  pct_neg = (float)neg / (float)total_vertices ;
  l_dist = parms->l_dist ;
  pct_neg_area =
    (float)mris->neg_area / (float)(mris->total_area+mris->neg_area) ;
  base_avgs = parms->n_averages ;
  for (t = 0 ;
       (t < niterations) && (neg > min_neg) && (pct_neg_area > min_neg_pct) ;
       t++)
  {
    if (pct_neg_area < 0.001)  /* hack!!, but it speeds things up */
    {
      parms->l_dist *= 1.1 ;
      if (parms->l_dist > 10*l_dist)
      {
        parms->l_dist = 10*l_dist ;
      }
      if (parms->l_dist > 1.0)
      {
        parms->l_dist = 1.0 ;
      }
    }
    if (pct_neg_area < 0.001)  /* another hack!!, but it speeds things up */
    {
      static int first = 1 ;
      /* don't want big steps or momentum for fine-scale stuff */
      parms->momentum = 0.0f ;
      parms->dt = 0.1 ;
      if (Gdiag & DIAG_SHOW && first)
        fprintf(stdout, "setting momentum=%2.1f, dt=%2.1f, l_dist=%2.2f\n",
                parms->momentum, parms->dt, parms->l_dist) ;
      first = 0 ;
    }

    if (mris->patch)  /* area is constant so spring force doesn't decrease */
    {
      scale = sqrt(mris->orig_area / (mris->total_area+mris->neg_area)) ;
      MRISscaleBrain(mris, mris, scale) ;
      MRIScomputeMetricProperties(mris) ;
    }
    else
    {
      mrisComputeVertexDistances(mris) ;
    }
    MRISclearGradient(mris) ;
    mrisComputeSpringTerm(mris, parms->l_spring) ;
    mrisComputeLaplacianTerm(mris, parms->l_lap) ;
    mrisComputeDistanceTerm(mris, parms) ;
    mrisComputeAngleAreaTerms(mris, parms) ;
    /*    mrisAverageGradient(mris, parms->n_averages) ;*/

    switch (parms->integration_type)
    {
    case INTEGRATE_LM_SEARCH:
      delta_t = mrisLineMinimizeSearch(mris, parms) ;
      break ;
    default:
    case INTEGRATE_LINE_MINIMIZE:
      delta_t = mrisLineMinimize(mris, parms) ;
      break ;
    case INTEGRATE_MOMENTUM:
      delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt,
                                     parms->tol, parms->n_averages) ;
      break ;
    case INTEGRATE_ADAPTIVE:
      delta_t = mrisAdaptiveTimeStep(mris, parms);
      break ;
    }
    mrisProjectSurface(mris) ;
    MRIScomputeMetricProperties(mris) ;
    neg = mrisCountNegativeVertices(mris) ;
    pct_neg = (float)neg / (float)total_vertices ;
    pct_neg_area =
      (float)mris->neg_area / (float)(mris->total_area+mris->neg_area) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout,
              "%3.3d: count: %d (%2.2f%%), area: %2.2f (%2.2f%%)   \n",
              t, neg, 100.0f*pct_neg, mris->neg_area, 100.0f*pct_neg_area) ;
    if ((write_iterations > 0) &&
        !((t+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
    {
      mrisWriteSnapshot(mris, parms, t+1) ;
    }
    if (parms->n_averages == 0)
    {
      parms->n_averages = base_avgs ;
    }
    else
    {
      parms->n_averages /= 2 ;
    }
  }

  parms->n_averages = base_avgs ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "\n") ;
    if (Gdiag & DIAG_WRITE)
    {
      fclose(parms->fp) ;
      parms->fp = NULL ;
    }
  }
  mrisProjectSurface(mris) ;
  return(mris) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI_SURFACE  *
MRISflatten(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     base_averages, n_averages, done, steps, total_steps ;
  float   base_tol ;
  double  starting_sse, ending_sse ;

  base_averages = parms->n_averages ;


  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    parms->fp = fopen(fname, "w") ;
    if (!parms->fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
                Progname, fname) ;
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    mrisLogIntegrationParms(stderr, mris, parms) ;
  }


  parms->start_t = 0 ;
  base_tol = parms->tol ;
  do
  {
    done = 0 ;
    mrisClearMomentum(mris) ;
    starting_sse = MRIScomputeSSE(mris, parms) ;
    for (total_steps = 0, n_averages = base_averages; !done ; n_averages /= 2)
    {
      steps = MRISintegrate(mris, parms, n_averages) ;
      parms->start_t += steps ;
      total_steps += steps ;
      done = n_averages == 0 ;   /* finished integrating
                                                    at smallest scale */
    }
    parms->dt = parms->base_dt ;         /* reset time step */
    ending_sse = MRIScomputeSSE(mris, parms) ;
  }
  while (!FZERO(ending_sse) && ((starting_sse-ending_sse) > parms->tol)) ;


  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  mrisProjectSurface(mris) ;
  return(mris) ;
}
#endif


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
static int
mrisSmoothNormals(MRI_SURFACE *mris, int niterations)
{
  int     vno, i, vn ;
  double  g ;
  VERTEX  *vertex, *vnb ;
  VECTOR  *v_n, *v_n2 ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_n2 = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */

  for (i = 0 ; i < niterations ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
      {
        continue ;
      }
      VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
      V3_SCALAR_MUL(v_n, kernel[0], v_n) ;
      g = kernel[1] ;
      for (vn = 0 ; vn < vertex->vnum ; vn++)
      {
        vnb = &mris->vertices[vertex->v[vn]] ;
        VECTOR_LOAD(v_n2, vnb->nx, vnb->ny, vnb->nz) ;
        V3_SCALAR_MUL(v_n2, g, v_n2) ;
        V3_ADD(v_n, v_n2, v_n) ;
      }
      g = kernel[2] ;
      for ( ; vn < vertex->v2num ; vn++)
      {
        vnb = &mris->vertices[vertex->v[vn]] ;
        VECTOR_LOAD(v_n2, vnb->nx, vnb->ny, vnb->nz) ;
        V3_SCALAR_MUL(v_n2, g, v_n2) ;
        V3_ADD(v_n, v_n2, v_n) ;
      }
      V3_NORMALIZE(v_n, v_n) ;
      vertex->tdx = V3_X(v_n) ;
      vertex->tdy = V3_Y(v_n) ;
      vertex->tdz = V3_Z(v_n) ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      vertex = &mris->vertices[vno] ;
      if (vertex->ripflag)
      {
        continue ;
      }
      vertex->nx = vertex->tdx ;
      vertex->ny = vertex->tdy ;
      vertex->nz = vertex->tdz ;
    }
  }

  VectorFree(&v_n) ;
  VectorFree(&v_n2) ;
  return(NO_ERROR) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisSmoothNormalOutliers(MRI_SURFACE *mris, double ndist)
{
  int     vno, n, m, smooth, nsmoothed ;
  VERTEX  *v, *vn ;
  float   sx, sy, sz, nx, ny, nz, nc, x, y, z, dist, dx, dy, dz ;


  for (nsmoothed = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    x = v->x ;
    y = v->y ;
    z = v->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    /* see if we are futher than ndist from all our neighbors in
       the normal direction. If so, smooth.
    */
    for (smooth = 1, m = 0 ; smooth && m < v->vnum ; m++)
    {
      vn = &mris->vertices[v->v[m]] ;
      if (!vn->ripflag)
      {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        nc = dx*nx+dy*ny+dz*nz;   /* projection onto normal */
        dist = sqrt(nc) ;         /* distance in normal direction */
        if (dist < ndist)
        {
          smooth = 0 ;
        }
        sx += dx ;
        sy += dy ;
        sz += dz ;
        n++;
      }
    }
    if (!smooth)
    {
      continue ;
    }
    nsmoothed++ ;
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
#if 0
    nc = sx*nx+sy*ny+sz*nz;   /* projection onto normal */
    sx = nc*nx ;              /* move in normal direction */
    sy = nc*ny ;
    sz = nc*nz;
#endif

    v->dx += sx ;
    v->dy += sy ;
    v->dz += sz ;
  }

  fprintf(stdout, "%d smoothed (%2.2f%%)\n",
          nsmoothed, 100.0f*(float)nsmoothed / (float)mris->nvertices) ;
  return(NO_ERROR) ;
}


static int
mrisComputeAverageNormalTerm(MRI_SURFACE *mris, int navgs, double l_normal)
{
  VERTEX  *v, *vn ;
  double  nc_avg, nc, vnum, delta ;
  int     n, vno, marked ;
  float   x, y, z, dx, dy, dz, nx, ny, nz, s ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    marked = v->marked ;

    /* compute projection onto normal */
    x = v->x ;
    y = v->y ;
    z = v->z ;
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    for (vnum = nc_avg = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked != marked)
      {
        continue ;
      }
      vnum++ ;
      dx = vn->x - x ;
      dy = vn->y - y ;
      dz = vn->z - z ;
      nc_avg += dx*nx+dy*ny+dz*nz;   /* projection onto normal */
    }
    if (vnum > 0.0)
    {
      nc_avg /= vnum ;
    }
    else
    {
      nc_avg = 0.0 ;
    }
    v->d = nc_avg ;
  }

  mrisAverageDs(mris, navgs) ;

  /* now move each vertex in the direction of the local average */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    marked = v->marked ;

    /* compute projection onto normal */
    nc_avg = v->d ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    for (nc = vnum = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked != marked)
      {
        continue ;
      }
      vnum++ ;
      dx = vn->x - x ;
      dy = vn->y - y ;
      dz = vn->z - z ;
      nc += dx*nx+dy*ny+dz*nz;   /* projection onto normal */
    }
    if (vnum > 0.0)
    {
      nc /= (float)vnum ;
    }
    else
    {
      nc = 0.0 ;
    }

    s = nc_avg < 0.0 ? -1.0 : 1.0 ;
    nc_avg = sqrt(fabs(nc_avg)) * s ;
    s = nc < 0.0 ? -1.0 : 1.0 ;
    nc = sqrt(fabs(nc)) * s ;
    delta = -l_normal * (nc_avg - nc) ;
    dx = nx * delta ;
    dy = ny * delta ;
    dz = nz * delta ;
    v->dx += dx ;
    v->dy += dy ;
    v->dz += dz ;
  }
  return(NO_ERROR) ;
}
static int
mrisComputeCurvatureTerm(MRI_SURFACE *mris, double l_curv)
{
  int     vno, n ;
  VERTEX  *v, *vn ;
  float   nc, nx, ny, nz, x, y, z, ndx, ndy, ndz, tdx, tdy, tdz, dx, dy, dz,
          tdist, ndist, nc_avg ;
  double  curv ;

  if (FZERO(l_curv))
  {
    return(NO_ERROR) ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    x = v->x ;
    y = v->y ;
    z = v->z ;

    /* first mean tangential and normal spacing */
    nc_avg = ndx = ndy = ndz = tdist = ndist = 0.0f ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dx = vn->x - x ;
      dy = vn->y - y ;
      dz = vn->z - z ;
      nc = (dx*nx + dy*ny + dz*nz) ;
      nc_avg += nc ;
#if 0
      ndx += nc*nx ;
      ndy += nc*ny ;
      ndz += nc*nz ;
#endif
      tdx = dx-nc*nx ;
      tdy = dy-nc*ny ;
      tdz = dz-nc*nz ;
      tdist += sqrt(tdx*tdx+tdy*tdy+tdz*tdz) ;
    }
#if 0
    ndx /= (float)v->vnum ;
    ndy /= (float)v->vnum ;
    ndz /= (float)v->vnum ;
    ndist = sqrt(ndx*ndx+ndy*ndy+ndz*ndz) ;
    if (nc_avg < 0.0f)
    {
      ndist *= -1 ;  /* neighbors are predominantly below tangent plane */
    }
#else
    ndist = nc_avg ;
#endif
    tdist /= (float)v->vnum ;
    ndist /= (float)v->vnum ;
    if (FZERO(tdist))
    {
      continue ;
    }
    curv = ndist / tdist ;

    if (fabs(curv) < 0.25f)
    {
      continue ;
    }
    curv *= l_curv ;
    v->dx += curv * nx ;
    v->dy += curv * ny ;
    v->dz += curv * nz ;
    if (vno == Gdiag_no)
      fprintf(stdout,
              "vertex %d normal curvature term: (%2.3f, %2.3f, %2.3f)\n"
              , vno, curv*nx, curv*ny, curv*nz) ;
  }

  return(NO_ERROR) ;
}
#endif


#if 0
*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static double
mrisComputeAverageHeight(MRI_SURFACE *mris)
{
  int    vno, n, nv ;
  VERTEX *vertex, *vn ;
  double height ;
  float  dx, dy, dz, nx, ny, nz, x, y, z ;

  for (height = 0.0, nv = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ;
    y = vertex->y ;
    z = vertex->z ;
    nx = vertex->nx ;
    ny = vertex->ny ;
    nz = vertex->nz ;
    if (vertex->ripflag)
    {
      continue ;
    }
    nv += vertex->vtotal ;
    for (n = 0 ; n < vertex->vtotal ; n++)
    {
      vn = &mris->vertices[vertex->v[n]] ;
      dx = x - vn->x ;
      dy = y - vn->y ;
      dz = z - vn->z ;
      height += fabs(dx * nx + dy * ny + dz * nz) ;
    }
  }

  return(sqrt(height/(double)nv)) ;
}
#endif

#if 0
static int
mrisSamplePialCoordsInTangentPlane(MRI_SURFACE *mris, VERTEX *v, float x, float y, float z, float *pxp, float *pyp, float *pzp)
{
  double  d, xi, yi, zi, u0, v0, u1, v1, wt, wt_total ;
  int     n ;
  VERTEX  *vn ;

  u0 = v->pe1x * x + v->pe1y*y + v->pe1z * z ;  // project point onto tangent plane
  v0 = v->pe2x * x + v->pe2y*y + v->pe2z * z ;
  for (wt_total = 0.0, n = -1 ; n < v->vnum ; n++)
  {
    if (n < 0)
    {
      vn = v ;
    }
    else
    {
      vn = &mris->vertices[v->v[n]] ;
    }
    u1 = v->pe1x * vn->pialx + v->pe1y*vn->pialy + v->pe1z * vn->pialz ;  // project nbr onto tangent plane
    v1 = v->pe2x * vn->pialx + v->pe2y*vn->pialy + v->pe2z * vn->pialz ;
    d = sqrt(SQR(u1-u0) + SQR(v1-v0)) ;
#define WT_VAR (0.1 * 0.1)
    wt = exp(-d*d/(2*WT_VAR)) ;
    wt_total += wt ;
  }
  xi = yi = zi = 0 ;
  for (n = -1 ; n < v->vnum ; n++)
  {
    if (n < 0)
    {
      vn = v ;
    }
    else
    {
      vn = &mris->vertices[v->v[n]] ;
    }
    u1 = v->pe1x * vn->pialx + v->pe1y*vn->pialy + v->pe1z * vn->pialz ;  // project nbr onto tangent plane
    v1 = v->pe2x * vn->pialx + v->pe2y*vn->pialy + v->pe2z * vn->pialz ;
    d = sqrt(SQR(u1-u0) + SQR(v1-v0));
    wt = exp(-d*d/(2*WT_VAR)) ;
    xi += wt * vn->pialx ;
    yi += wt * vn->pialy ;
    zi += wt * vn->pialz ;
  }
  *pxp = xi/wt_total  ;
  *pyp = yi/wt_total ;
  *pzp = zi/wt_total ;
  return(NO_ERROR) ;
}
#endif

#if 0
static double
ashburnerTriangleEnergy(MRI_SURFACE *mris, int fno, double lambda)
{
  static MATRIX *m_x, *m_x_inv, *m_y, *m_J, *m_U, *m_V, *m_m, 
    *m_evectors = NULL ;
  FACE   *f ;
  VERTEX *v1, *v2, *v3 ;
  float  e1x, e1y, e1z, e2x, e2y, e2z , x11, x12, x13, x21, x22, x23, 
    y11, y12, y13, y21, y22, y23, evalues[2], det, s11, s22 ;
  double energy /*, a, b, c, d*/ ;

  if (m_x == NULL)
  {
    m_x = MatrixAlloc(3,3,MATRIX_REAL) ;
    m_m = MatrixAlloc(3,3,MATRIX_REAL) ;
    m_x_inv = MatrixAlloc(3,3,MATRIX_REAL) ;
    m_y = MatrixAlloc(3,3,MATRIX_REAL) ;
    m_J = MatrixAlloc(2,2,MATRIX_REAL) ;
    m_U = MatrixAlloc(2,2,MATRIX_REAL) ;
    m_V = MatrixAlloc(2,2,MATRIX_REAL) ;
  }

  f = &mris->faces[fno] ;
  if (fno == Gdiag_no || f->area < 0)
  {
    DiagBreak();
  }
  get_face_axes(mris, f, &e1x, &e1y, &e1z, &e2x, &e2y, &e2z) ;
  // first index refers to coord # and second to which vertex
  v1 = &mris->vertices[f->v[0]] ;
  v2 = &mris->vertices[f->v[1]] ;
  v3 = &mris->vertices[f->v[2]] ;

  // original coords in tangent plane
  x11 = v1->cx * e1x + v1->cy * e1y + v1->cz * e1z ;
  x21 = v1->cx * e2x + v1->cy * e2y + v1->cz * e2z ;
  x12 = v2->cx * e1x + v2->cy * e1y + v2->cz * e1z ;
  x22 = v2->cx * e2x + v2->cy * e2y + v2->cz * e2z ;
  x13 = v3->cx * e1x + v3->cy * e1y + v3->cz * e1z ;
  x23 = v3->cx * e2x + v3->cy * e2y + v3->cz * e2z ;

  // mapped coords in tangent plane
  y11 = v1->x * e1x + v1->y * e1y + v1->z * e1z ;
  y21 = v1->x * e2x + v1->y * e2y + v1->z * e2z ;
  y12 = v2->x * e1x + v2->y * e1y + v2->z * e1z ;
  y22 = v2->x * e2x + v2->y * e2y + v2->z * e2z ;
  y13 = v3->x * e1x + v3->y * e1y + v3->z * e1z ;
  y23 = v3->x * e2x + v3->y * e2y + v3->z * e2z ;


  *MATRIX_RELT(m_x, 1, 1) = x11 ;
  *MATRIX_RELT(m_x, 1, 2) = x12 ;
  *MATRIX_RELT(m_x, 1, 3) = x13 ;
  *MATRIX_RELT(m_x, 2, 1) = x21 ;
  *MATRIX_RELT(m_x, 2, 2) = x22 ;
  *MATRIX_RELT(m_x, 2, 3) = x23 ;
  *MATRIX_RELT(m_x, 3, 1) = 1.0 ;
  *MATRIX_RELT(m_x, 3, 2) = 1.0 ;
  *MATRIX_RELT(m_x, 3, 3) = 1.0 ;

  *MATRIX_RELT(m_y, 1, 1) = y11 ;
  *MATRIX_RELT(m_y, 1, 2) = y12 ;
  *MATRIX_RELT(m_y, 1, 3) = y13 ;
  *MATRIX_RELT(m_y, 2, 1) = y21 ;
  *MATRIX_RELT(m_y, 2, 2) = y22 ;
  *MATRIX_RELT(m_y, 2, 3) = y23 ;
  *MATRIX_RELT(m_y, 3, 1) = 1.0 ;
  *MATRIX_RELT(m_y, 3, 2) = 1.0 ;
  *MATRIX_RELT(m_y, 3, 3) = 1.0 ;

#define SMALL 1e-8
  if (MatrixInverse(m_x, m_x_inv) == NULL)
  {
    return(log(SMALL)*log(SMALL)) ;
  }
  MatrixMultiply(m_y, m_x_inv, m_m) ;
  MatrixCopyRegion(m_m, m_J, 1, 1, 2, 2, 1,1) ;
  //  MatrixSVDEigenValues(m_J, evalues) ;
  det = MatrixDeterminant(m_J) ;
  if (det < 0)
  {
    return(log(SMALL)*log(SMALL)) ;
  }
  m_evectors = MatrixEigenSystem(m_J, evalues, m_evectors) ;
  s11 = evalues[0] ;
  s22 = evalues[1] ;

#if 0
  a = *MATRIX_RELT(m_m, 1, 1) ;
  b = *MATRIX_RELT(m_m, 1, 2) ;
  c = *MATRIX_RELT(m_m, 2, 1) ;
  d = *MATRIX_RELT(m_m, 2, 2) ;
  det = a*d - b*c ;
  s11 = ((a+d)/2) + sqrt((4*b*c + (a-d)*(a-d)) / 2);
  s22 = ((a+d)/2) - sqrt((4*b*c + (a-d)*(a-d)) / 2);
#endif
  if (!isfinite(s11))
    return(log(SMALL)*log(SMALL)) ;

  if (!isfinite(s22))
    return(log(SMALL)*log(SMALL)) ;

  if (s11 <= 0)
    s11 = SMALL ;

  s11 = log(s11) ;
  s11 *= s11 ;

  if (s22 <= 0)
  {
    s22 = SMALL ;
  }
  s22 = log(s22) ;
  s22 *= s22 ;

  energy = lambda * (1 + det) * (s11 + s22) ;  // log and square of snn already taken

  if (!isfinite(energy))
    DiagBreak() ;

  return(energy/2) ;
}
#else
#endif


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Expand the list of neighbors of each vertex, reallocating
  the v->v array to hold the expanded list.
  ------------------------------------------------------*/
#define MAX_V 1000 /* max for any one node, actually way too big */
int
MRISsampleDistances(MRI_SURFACE *mris, int *nbrs, int max_nbhd)
{
  int          i,n, vno, vnum, old_vnum, total_nbrs, max_possible,max_v,vtotal;
  VERTEX       *v, *vn, *vn2 ;
  int          *vnbrs, *vall, found, n2, vnbrs_num, vall_num, nbhd_size,done ;
  float        xd, yd, zd, min_dist, dist, dist_scale, old_dist[MAX_V],
               old_v[MAX_V], min_angle, angle ;
  VECTOR       *v1, *v2 ;

  printf("Starting MRISsampleDistances %d ------ ()\n",max_nbhd);
  for (i = 0; i < max_nbhd ; i++)
  {
    printf("%d %d\n",i,nbrs[i]);
  }
  printf("randSeed %ld",getRandSeed());

  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;

  dist_scale = (1.0 + sqrt(2.0)) / 2 ;  /* adjust for Manhattan distance */
  vnbrs = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int)) ;
  vall = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int)) ;
  total_nbrs = 0 ;
  for (vtotal = max_possible = 0, n = 1 ; n <= max_nbhd ; n++)
  {
    max_possible += nbrs[n] ;
    if (n > mris->nsize)
    {
      vtotal += nbrs[n] ;
    }
  }

  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stdout,
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal,
            (float)vtotal/((float)max_nbhd-(float)mris->nsize),
            (float)vtotal*mris->nvertices*sizeof(float)*3.0f /
            (1024.0f*1024.0f)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices/25))))
      fprintf(stdout, " %%%1.0f",
              100.0f*(float)vno / (float)mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      DiagBreak()  ;
    }

    if (v->ripflag)
    {
      continue ;
    }

    /* small neighborhood is always fixed, don't overwrite them */
    vtotal = v->vtotal ;
    if (v->nsize == 3)
    {
      v->vtotal = v->v3num ;
    }
    else if (v->nsize == 2)
    {
      v->vtotal = v->v2num ;
    }
    else
    {
      v->vtotal = v->vnum ;
    }

    max_v = v->vtotal+max_possible ;
    if (vtotal < max_v) /* won't fit in current allocation,
                                         reallocate stuff */
    {
      /* save and restore neighbor list */
      memmove(old_v, v->v, v->vtotal*sizeof(v->v[0])) ;
      free(v->v) ;
      v->v = (int *)calloc(max_v, sizeof(int)) ;
      if (!v->v)
        ErrorExit(ERROR_NO_MEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "nbrs at v=%d", max_v, vno) ;
      memmove(v->v, old_v, v->vtotal*sizeof(v->v[0])) ;

      /* save and restore distance vector */
      memmove(old_dist, v->dist, v->vtotal*sizeof(v->dist[0])) ;
      free(v->dist) ;
      v->dist = (float *)calloc(max_v, sizeof(float)) ;
      if (!v->dist)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d", max_v, vno) ;
      memmove(v->dist, old_dist, v->vtotal*sizeof(v->dist[0])) ;

      /* save and restore original distance vector */
      memmove(old_dist, v->dist_orig, v->vtotal*sizeof(v->dist_orig[0])) ;
      free(v->dist_orig) ;
      v->dist_orig = (float *)calloc(max_v, sizeof(float)) ;
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d", max_v, vno) ;
      memmove(v->dist_orig, old_dist, v->vtotal*sizeof(v->dist_orig[0])) ;
    }

    vall[0] = vno ;
    vall_num = 1 ;
    old_vnum = 0 ;
    v->marked = 1 ;  /* a hack - it is a zero neighbor */
    for (nbhd_size = 1 ; vall_num < MAX_NBHD_VERTICES && nbhd_size <= max_nbhd ;
         nbhd_size++)
    {
      /* expand neighborhood outward by a ring of vertices */
      vnbrs_num = 0 ;  /* will count neighbors in this ring */
      vnum = vall_num ;
      for (found = 0, n = old_vnum;
           vall_num<MAX_NBHD_VERTICES && n < vall_num;
           n++)
      {
        vn = &mris->vertices[vall[n]] ;
        if (vn->ripflag)
        {
          continue ;
        }

        /* search through vn's neighbors to find an unmarked vertex */
        for (n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (vn2->ripflag || vn2->marked)
          {
            continue ;
          }

          /* found one, mark it and put it in the vall list */
          found++ ;
          vn2->marked = nbhd_size ;
          vall[vnum++] = vn->v[n2] ;
          if (nbrs[nbhd_size] > 0)  /* want to store this distance */
          {
            vnbrs[vnbrs_num++] = vn->v[n2] ;
          }
        }
      }  /* done with all neighbors at previous distance */

      /* found all neighbors at this extent - calculate distances */
      old_vnum = vall_num ;
      vall_num += found ;
      for (n = old_vnum ; n < vall_num ; n++)
      {
        vn = &mris->vertices[vall[n]] ;
        for (min_dist = 10000.0, n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (!vn2->marked)
          {
            continue ;
          }
          xd = vn2->x - vn->x ;
          yd = vn2->y - vn->y ;
          zd = vn2->z - vn->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
          if (nbhd_size > 1)
          {
            dist /= dist_scale ;
          }
          if (vn2->d + dist < min_dist)
          {
            min_dist = vn2->d + dist ;
          }
        }
        vn->d = min_dist  ;
        if (nbhd_size <= 2)
        {
          xd = vn->x - v->x ;
          yd = vn->y - v->y ;
          zd = vn->z - v->z ;
          dist = sqrt(xd*xd + yd*yd + zd*zd) ;
          vn->d = dist ;
        }
      }

      /* if this set of neighbors are to be stored, sample from them */
      if (nbrs[nbhd_size] <= 0)
      {
        continue ;
      }

      /* make sure the points are not too close together */
      min_angle = 0.9*2.0*M_PI / (float)nbrs[nbhd_size] ;

      /*
        at this point, the vall list contains all
        the neighbors currently found
        at ALL distances, while the vnbrs list contains ONLY the
        nbhd_size-neighbors.
      */
      if (found <= nbrs[nbhd_size])  /* just copy them all in */
      {
        for (n = 0 ; n < found ; n++, v->vtotal++)
        {
          v->v[v->vtotal] = vnbrs[n] ;
          v->dist_orig[v->vtotal] = mris->vertices[vnbrs[n]].d ;
        }
      }
      else                   /* randomly sample from them */
      {
        int vstart = v->vtotal ;
        for (n = 0 ; n < nbrs[nbhd_size] ; n++, v->vtotal++)
        {
          int j, niter = 0 ;
          do
          {
            do
            {
              i = nint(randomNumber(0.0, (double)found-1)) ;
            }
            while (vnbrs[i] < 0) ;
            /*
              now check to make sure that the angle between this
              point and the others already selected is not too
              small to make sure the points are not bunched.
            */
            vn = &mris->vertices[vnbrs[i]] ;
            VECTOR_LOAD(v1, vn->x-v->x, vn->y-v->y, vn->z-v->z) ;
            done = 1 ;
            for (j = vstart ; done && j < v->vtotal ; j++)
            {
              vn2 = &mris->vertices[v->v[j]] ;
              VECTOR_LOAD(v2,
                          vn2->x-v->x, vn2->y-v->y, vn2->z-v->z) ;
              angle = Vector3Angle(v1, v2) ;
              if (angle < min_angle)
              {
                done = 0 ;
              }
            }
            if (++niter > found)  /* couldn't find
                                     enough at this difference */
            {
              min_angle *= 0.75f ;  /* be more liberal */
              niter = 0 ;
            }
          }
          while (!done && !FZERO(min_angle)) ;
          vn = &mris->vertices[vnbrs[i]] ;
          v->v[v->vtotal] = vnbrs[i] ;
          v->dist_orig[v->vtotal] = vn->d ;
          vnbrs[i] = -1 ;
        }
      }
    }


    if ((Gdiag_no == vno) && DIAG_VERBOSE_ON)
    {
      FILE  *fp ;
      char  fname[STRLEN] ;

      sprintf(fname, "v%d", vno) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "%d\n", vall_num) ;
      for (n = 0 ; n < vall_num ; n++)
      {
        fprintf(fp, "%d\n", vall[n]) ;
      }
      fclose(fp) ;

      sprintf(fname, "vn%d", vno) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "%d\n", v->vtotal) ;
      for (n = 0 ; n < v->vtotal ; n++)
      {
        fprintf(fp, "%d\n", v->v[n]) ;
      }
      fclose(fp) ;
      for (n = 0 ; n < mris->nvertices ; n++)
      {
        vn = &mris->vertices[n] ;
#if 0
        if (vn->ripflag)
        {
          continue ;
        }
#endif
        vn->curv = vn->d ;
      }
      sprintf(fname, "%s.dist",
              mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh");
      MRISwriteCurvature(mris, fname) ;
    }

    /*
      done building arrays - allocate distance vectors and
      sample from the found neighbors list.
    */
    /* now unmark them all */
    for (n = 0 ; n < vall_num ; n++)
    {
      mris->vertices[vall[n]].marked = 0 ;
      mris->vertices[vall[n]].d = 0.0 ;
    }

    total_nbrs += v->vtotal ;
  }

  /* now fill in immediate neighborhood(Euclidean) distances */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      DiagBreak()  ;
    }
    if (v->ripflag)
    {
      continue ;
    }
    if (v->nsize == 3)
    {
      vtotal = v->v3num ;
    }
    else if (v->nsize == 2)
    {
      vtotal = v->v2num ;
    }
    else
    {
      vtotal = v->vnum ;
    }
    for (n = 0 ; n < vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
      {
        continue ;
      }
      xd = v->x - vn->x ;
      yd = v->y - vn->y ;
      zd = v->z - vn->z ;
      v->dist_orig[n] = sqrt(xd*xd+yd*yd+zd*zd) ;
    }
  }

  mris->avg_nbrs = (float)total_nbrs / (float)MRISvalidVertices(mris) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    fprintf(stdout, "avg_nbrs = %2.1f\n", mris->avg_nbrs) ;
  }

  free(vnbrs) ;
  free(vall) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  if (Gdiag & DIAG_HEARTBEAT)
  {
    fprintf(stdout, " done.\n") ;
  }
  printf("randSeed %ld",getRandSeed());
  printf("MRISsampleDistances(): done\n");

  return(NO_ERROR) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisClipGradient(MRI_SURFACE *mris, float max_len)
{
  int     vno ;
  VERTEX  *v ;
  float   dx, dy, dz, len ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    dx = v->dx ;
    dy = v->dy ;
    dz = v->dz ;
    len = sqrt(dx*dx+dy*dy+dz*dz) ;
    if (len > max_len)
    {
      len = max_len / len ;
      v->dx *= len ;
      v->dy *= len ;
      v->dz *= len ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisClipMomentumGradient(MRI_SURFACE *mris, float max_len)
{
  int     vno ;
  VERTEX  *v ;
  float   dx, dy, dz, len ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    dx = v->odx ;
    dy = v->ody ;
    dz = v->odz ;
    len = sqrt(dx*dx+dy*dy+dz*dz) ;
    if (len > max_len)
    {
      len = max_len / len ;
      v->odx *= len ;
      v->ody *= len ;
      v->odz *= len ;
    }
  }
  return(NO_ERROR) ;
}
#endif


#if 1
#else
// these versions use a 1d fit y = a*r^2 + b
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Fit a 1-d quadratic to the surface locally and move the
  vertex in the normal direction to improve the fit.
  ------------------------------------------------------*/
static int mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv)
{
  MATRIX *m_R, *m_R_inv;
  VECTOR *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr;
  int vno, n;
  VERTEX *v, *vn;
  float ui, vi, rsq, a, b;
  FILE *fp;  // for debugging

  if (FZERO(l_curv)) {
    return (NO_ERROR);
  }

  mrisComputeTangentPlanes(mris);
  v_n = VectorAlloc(3, MATRIX_REAL);
  v_A = VectorAlloc(2, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v_nbr = VectorAlloc(3, MATRIX_REAL);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    if (vno == Gdiag_no) fp = fopen("qcurv.dat", "w");

    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL);    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 2, MATRIX_REAL); /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz);
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z);
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z);
    for (n = 0; n < v->vtotal; n++) /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]];
      VERTEX_EDGE(v_nbr, v, vn);
      VECTOR_ELT(v_Y, n + 1) = V3_DOT(v_nbr, v_n);
      ui = V3_DOT(v_e1, v_nbr);
      vi = V3_DOT(v_e2, v_nbr);
      rsq = ui * ui + vi * vi;
      *MATRIX_RELT(m_R, n + 1, 1) = rsq;
      *MATRIX_RELT(m_R, n + 1, 2) = 1;

      if (vno == Gdiag_no) fprintf(fp, "%d %f %f %f %f\n", v->v[n], ui, vi, rsq, VECTOR_ELT(v_Y, n + 1));
    }
    if (vno == Gdiag_no) fclose(fp);
    m_R_inv = MatrixPseudoInverse(m_R, NULL);
    if (!m_R_inv) {
      MatrixFree(&m_R);
      VectorFree(&v_Y);
      continue;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A);
    a = VECTOR_ELT(v_A, 1);
    b = VECTOR_ELT(v_A, 2);

    v->dx += b * v->nx;
    v->dy += b * v->ny;
    v->dz += b * v->nz;

    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "a=%2.2f, b=%2.1f\n",
              vno,
              b * v->nx,
              b * v->ny,
              b * v->nz,
              a,
              b);
    VectorFree(&v_Y);
    MatrixFree(&m_R);
    MatrixFree(&m_R_inv);
  }

  VectorFree(&v_n);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  VectorFree(&v_nbr);
  VectorFree(&v_A);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Fit a 1-d quadratic to the surface locally and move the
  vertex in the normal direction to improve the fit.
  ------------------------------------------------------*/
static double mrisComputeQuadraticCurvatureSSE(MRI_SURFACE *mris, double l_curv)
{
  MATRIX *m_R, *m_R_inv;
  VECTOR *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr;
  int vno, n;
  VERTEX *v, *vn;
  float ui, vi, rsq, a, b;
  double sse = 0.0;

  if (FZERO(l_curv)) {
    return (0.0);
  }

  mrisComputeTangentPlanes(mris);
  v_n = VectorAlloc(3, MATRIX_REAL);
  v_A = VectorAlloc(2, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v_nbr = VectorAlloc(3, MATRIX_REAL);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL);    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 2, MATRIX_REAL); /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz);
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z);
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z);
    for (n = 0; n < v->vtotal; n++) /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]];
      VERTEX_EDGE(v_nbr, v, vn);
      VECTOR_ELT(v_Y, n + 1) = V3_DOT(v_nbr, v_n);
      ui = V3_DOT(v_e1, v_nbr);
      vi = V3_DOT(v_e2, v_nbr);
      rsq = ui * ui + vi * vi;
      *MATRIX_RELT(m_R, n + 1, 1) = rsq;
      *MATRIX_RELT(m_R, n + 1, 2) = 1;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL);
    if (!m_R_inv) {
      MatrixFree(&m_R);
      VectorFree(&v_Y);
      continue;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A);
    a = VECTOR_ELT(v_A, 1);
    b = VECTOR_ELT(v_A, 2);
    sse += b * b;
    if (vno == Gdiag_no) {
      printf("v %d: curvature sse %2.2f\n", vno, b * b);
    }
    MatrixFree(&m_R);
    VectorFree(&v_Y);
    MatrixFree(&m_R_inv);
  }

  VectorFree(&v_n);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  VectorFree(&v_nbr);
  VectorFree(&v_A);
  return (sse);
}
#endif

#if 0
static int
mrisAverageDs(MRI_SURFACE *mris, int num_avgs)
{
  VERTEX  *v, *vn ;
  double  d, vnum ;
  int     n, vno, i, marked ;

  /* now average them */
  for (i = 0 ; i < num_avgs ; i++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }
      marked = v->marked ;
      d = v->d ;
      for (vnum = 1.0, n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked != marked)
        {
          continue ;
        }
        vnum++ ;
        d += vn->d ;
      }
      if (vnum > 0.0)
      {
        d /= vnum ;
      }
      else
      {
        d = 0.0 ;
      }
      v->tdx = d ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }
      v->d = v->tdx ;
    }
  }

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
#define USE_TANGENTIAL_TERM 0
#define SCALE_BY_N 1
static int
mrisComputeCurvatureTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_n, *v_delta ;
  int     vno ;
  VERTEX  *vertex ;
  double  l_curv, deltaH, Hdesired ;

  l_curv = parms->l_curv ;
  if (FZERO(l_curv))
  {
    return(NO_ERROR) ;
  }

  Hdesired = parms->Hdesired ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;

    deltaH = (Hdesired - vertex->H) ;  /* sign will be reversed in MUL */
#if 1
#define TAN_SCALE 0.5
    deltaH = tanh(deltaH * TAN_SCALE) ;
#endif
    V3_SCALAR_MUL(v_n, -deltaH, v_delta) ;

    vertex->dx += l_curv * V3_X(v_delta) ;
    vertex->dy += l_curv * V3_Y(v_delta) ;
    vertex->dz += l_curv * V3_Z(v_delta) ;
    if (vno == Gdiag_no && DIAG_VERBOSE_ON)
      fprintf(stdout, "Hdes=%2.3f, dH = %2.3f, tanh= %2.3f, dx=%2.3f\n",
              Hdesired,Hdesired - vertex->H, deltaH, vertex->dx) ;
  }

  VectorFree(&v_delta) ;
  VectorFree(&v_n) ;

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
static int
mrisComputeSethianCurvatureTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_n, *v_delta ;
  int     vno ;
  VERTEX  *vertex ;
  double  l_curv, deltaH, Hdesired ;

  Hdesired = parms->Hdesired ;
  l_curv = parms->l_curv ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }


    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    /* Osher and Sethian */
#if 1
    deltaH = (1 - parms->epsilon*tanh(Hdesired-vertex->H)*.5) ;
#else
    deltaH = (tanh(vertex->H-Hdesired)*.5) ;
#endif
    V3_SCALAR_MUL(v_n, deltaH, v_delta) ;

    vertex->dx += l_curv * V3_X(v_delta) ;
    vertex->dy += l_curv * V3_Y(v_delta) ;
    vertex->dz += l_curv * V3_Z(v_delta) ;
  }

  VectorFree(&v_delta) ;
  VectorFree(&v_n) ;

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
static int
mrisComputeCurvatureGradientTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_n, *v_g1, *v_g2, *v_y, *v_delta, *v_tmp, *v_n_r2, *v_u_g1,*v_v_g2;
  int     vno, i, N ;
  VERTEX  *vertex, *vnb ;
  double  r2, r3, z, u, v, l_curv, deltaH, Hdesired ;
  FILE    *fp = NULL ;

  Hdesired = parms->Hdesired ;
  l_curv = parms->l_curv ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_g1 = VectorAlloc(3, MATRIX_REAL) ;
  v_g2 = VectorAlloc(3, MATRIX_REAL) ;
  v_y = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;
  v_tmp = VectorAlloc(3, MATRIX_REAL) ;
  v_n_r2 = VectorAlloc(3, MATRIX_REAL) ;
  v_u_g1 = VectorAlloc(3, MATRIX_REAL) ;
  v_v_g2 = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    /*
      first compute term which comes from moving the this vertex in its
      own tangent plane.
    */
#if 1
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    VECTOR_LOAD(v_g1, vertex->e1x, vertex->e1y, vertex->e1z) ;
    VECTOR_LOAD(v_g2, vertex->e2x, vertex->e2y, vertex->e2z) ;
    V3_CLEAR(v_delta) ;

    if (Gdiag & DIAG_SHOW && vno == Gdiag_no && DIAG_VERBOSE_ON)
    {
      char fname[STRLEN] ;

      sprintf(fname, "v%d_%d.m", vno, parms->t) ;
      fp = fopen(fname, "w") ;
      fprintf(fp, "U%d_%d = [", vno, parms->t) ;
    }


    deltaH = (Hdesired - vertex->H) ;
    N = 0 ;
    for (i = 0 ; i < vertex->v2num ; i++)  /* for each neighbor */
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag)
      {
        continue ;
      }
      VECTOR_LOAD(v_y,
                  vnb->x-vertex->x, vnb->y-vertex->y, vnb->z-vertex->z) ;

      /* calculate projection onto tangent plane */
      u = V3_DOT(v_y, v_g1) ;
      v = V3_DOT(v_y, v_g2) ;
      r2 = u*u + v*v ;
      r3 = r2 * sqrt(r2) ;
      if (FZERO(r3))
      {
        continue ;
      }
      z = V3_DOT(v_y, v_n) ;     /* height above tangent plane */
      if (Gdiag & DIAG_SHOW && vno == Gdiag_no && DIAG_VERBOSE_ON)
      {
        fprintf(fp, "%2.3f  %2.3f  %2.3f; ", u, v, z) ;
      }
      V3_SCALAR_MUL(v_n, -1.0/r2, v_n_r2) ;  /* -n/r^2 */
      V3_SCALAR_MUL(v_g1, u, v_u_g1) ;       /* ui g1 */
      V3_SCALAR_MUL(v_g2, v, v_v_g2) ;       /* vi g2 */
      V3_ADD(v_u_g1, v_v_g2, v_tmp) ;        /* ui g1 + vi g2 */
      V3_SCALAR_MUL(v_tmp, 4*z/r3, v_tmp) ;/*  4 z / n^3 (ui g1 + vi g2) */
#if USE_TANGENTIAL_TERM
      V3_ADD(v_tmp, v_delta, v_delta) ;      /* add it into total delta */
#endif
      V3_ADD(v_n_r2, v_delta, v_delta) ;     /* add it into total delta */
      N++ ;
    }
    if (N > 0)
#if SCALE_BY_N
      V3_SCALAR_MUL(v_delta, deltaH * 2.0/N, v_delta) ;
#else
      V3_SCALAR_MUL(v_delta, deltaH * 2.0, v_delta) ;
#endif

#endif

    if (Gdiag & DIAG_SHOW && vno == Gdiag_no)
    {
      fprintf(fp, "] ;\n") ;
      fclose(fp) ;
    }

    /*
      now add terms which come from this vertex's appearance in
      neighboring tangent planes.
    */
    for (i = 0 ; i < vertex->v2num ; i++)  /* for each neighbor */
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag || !vnb->v2num)
      {
        continue ;
      }

      /* load coordinate system for neighbor's tangent plane */
      VECTOR_LOAD(v_n, vnb->nx, vnb->ny, vnb->nz) ;
      VECTOR_LOAD(v_g1, vnb->e1x, vnb->e1y, vnb->e1z) ;
      VECTOR_LOAD(v_g2, vnb->e2x, vnb->e2y, vnb->e2z) ;

      deltaH = (Hdesired - vnb->H) ;
      VECTOR_LOAD(v_y,
                  vertex->x-vnb->x, vertex->y-vnb->y, vertex->z-vnb->z) ;

      /* calculate projection onto tangent plane */
      u = V3_DOT(v_y, v_g1) ;
      v = V3_DOT(v_y, v_g2) ;
      r2 = u*u + v*v ;
      r3 = r2 * sqrt(r2) ;
      if (FZERO(r3))
      {
        continue ;
      }
      z = V3_DOT(v_y, v_n) ;        /* height above tangent plane */
      V3_SCALAR_MUL(v_n, 1.0/r2, v_n_r2) ;   /* n/r^2 */
      V3_SCALAR_MUL(v_g1, u, v_u_g1) ;       /* ui g1 */
      V3_SCALAR_MUL(v_g2, v, v_v_g2) ;       /* vi g2 */
      V3_ADD(v_u_g1, v_v_g2, v_tmp) ;        /* ui g1 + vi g2 */

#if USE_TANGENTIAL_TERM
      V3_SCALAR_MUL(v_tmp, -4*z/r3, v_tmp) ;  /*  -4z / n^3 (ui g1 + vi g2) */
      V3_ADD(v_n_r2, v_tmp, v_tmp) ;
#else
      V3_SCALAR_MUL(v_n_r2, 1.0, v_tmp) ;
#endif
#if SCALE_BY_N
      V3_SCALAR_MUL(v_tmp, deltaH*2/(double)vnb->v2num, v_tmp) ;
#else
      V3_SCALAR_MUL(v_tmp, deltaH*2, v_tmp) ;
#endif
      V3_ADD(v_tmp, v_delta, v_delta) ;
    }

    vertex->dx += l_curv * V3_X(v_delta) ;
    vertex->dy += l_curv * V3_Y(v_delta) ;
    vertex->dz += l_curv * V3_Z(v_delta) ;
    if (vno == Gdiag_no)
      fprintf(stdout, "moving v %d by (%2.3f, %2.3f, %2.3f)\n",
              vno, vertex->dx, vertex->dy, vertex->dz) ;
  }

  VectorFree(&v_tmp) ;
  VectorFree(&v_delta) ;
  VectorFree(&v_y) ;
  VectorFree(&v_n) ;
  VectorFree(&v_g1) ;
  VectorFree(&v_g2) ;
  VectorFree(&v_n_r2) ;
  VectorFree(&v_u_g1) ;
  VectorFree(&v_v_g2) ;

  return(NO_ERROR) ;
}
#endif

#if 1
#else
static int mrisComputeIntensityGradientTerm(MRI_SURFACE *mris, double l_grad, MRI *mri_brain, MRI *mri_smooth)
{
  int vno;
  VERTEX *v;
  float x, y, z, nx, ny, nz;
  double xw, yw, zw, dx, dy, dz, val, mag, mag_next, scale;

  if (FZERO(l_grad)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;
    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    MRISsurfaceRASToVoxelCached(mris, mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz);
    mag = sqrt(dx * dx + dy * dy + dz * dz);
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val);

    /* sample outward from surface */
    xw = x + v->dx;
    yw = y + v->dy;
    zw = z + v->dz;
    MRISsurfaceRASToVoxelCached(mris, mri_smooth, xw, yw, zw, &xw, &yw, &zw);
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz);
    mag_next = sqrt(dx * dx + dy * dy + dz * dz);

    /* gradient decreasing in outwards direction */
    scale = 1.0 + l_grad;
    if (fabs(mag) > fabs(mag_next)) {
      scale /= 1.0f;
    }

    v->dx *= scale;
    v->dy *= scale;
    v->dz *= scale;
    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d intensity gradient term: (%2.3f, %2.3f, %2.3f) "
              "(mag = %2.1f, mag_next=%2.1f)\n",
              vno,
              v->dx,
              v->dy,
              v->dz,
              mag,
              mag_next);
  }

  return (NO_ERROR);
}
#endif

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  See if the line in the normal direction passing
  through vertex v intersects any of the triangles at the
  given location.
  ------------------------------------------------------*/
#if 0
static int
mrisNormalDirectionTriangleIntersection(MRI_SURFACE *mris, VERTEX *v,
                                        MHT *mht, double *pdist, int *flist, int which)
{
  double  dist, min_dist, U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3], dot ;
  float   nx, ny, nz, x, y, z, dx, dy, dz ;
  MHBT    *bucket ;
  MHB     *bin ;
  int     i, found, fno, ret ;


  dist = *pdist ;
  nx = v->nx ;
  ny = v->ny ;
  nz = v->nz ;
  dir[0] = v->nx ;
  dir[1] = v->ny ;
  dir[2] = v->nz ;
  pt[0] = v->origx  ;
  pt[1] = v->origy  ;
  pt[2] = v->origz  ;
  x = v->origx + nx * dist ;
  y = v->origy + ny * dist ;
  z = v->origz + nz * dist ;

  bucket = MHTgetBucket(mht, x, y, z) ;
  if (bucket == NULL)
  {
    return(-1) ;
  }

  min_dist = 10000.0f ;
  for (bin = bucket->bins, found = i = 0 ; i < bucket->nused ; i++, bin++)
  {
    fno = bin->fno ;

    load_triangle_vertices(mris, fno, U0, U1, U2, which) ;
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
    if (ret)
    {
      dx = int_pt[0] - v->origx ;
      dy = int_pt[1] - v->origy ;
      dz = int_pt[2] - v->origz ;
      dist = sqrt(dx*dx + dy*dy + dz*dz) ;
      dot = dx*nx + dy*ny + dz*nz ;
      if (dot >= 0 && dist < min_dist)
      {
        if (flist)
        {
          flist[found] = fno ;
        }
        found++ ;
        if (found ==1000)
          ErrorExit
          (ERROR_BADPARM,
           "too many intersecting faces.  "
           "check filled volume for correctness\n");
        *pdist = min_dist = dist ;
      }
    }
  }
  return(found) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  See if the line in the normal direction passing
  through vertex v intersects any of the triangles at the
  given location.
  ------------------------------------------------------*/
static int
mrisAllCurrentTriangleIntersections(MRI_SURFACE *mris, float x, float y,
                                    float z, float nx, float ny, float nz,
                                    MHT *mht, int *flist)
{
  double  U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3] ;
  MHBT    *bucket ;
  MHB     *bin ;
  int     i, found, fno, ret ;


  dir[0] = nx ;
  dir[1] = ny ;
  dir[2] = nz ;
  pt[0] = x  ;
  pt[1] = y  ;
  pt[2] = z  ;

  bucket = MHTgetBucket(mht, x, y, z) ;
  if (bucket == NULL)
  {
    return(-1) ;
  }

  for (bin = bucket->bins, found = i = 0 ; i < bucket->nused ; i++, bin++)
  {
    fno = bin->fno ;

    load_triangle_vertices(mris, fno, U0, U1, U2, CURRENT_VERTICES) ;
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt) ;
    if (ret)
    {
      flist[found] = fno ;
      found++ ;
      if (found ==1000)
        ErrorExit
        (ERROR_BADPARM,
         "too many intersecting faces.  "
         "check filled volume for correctness\n");
    }
  }
  return(found) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static double
mrisFindClosestFilledVoxel(MRI_SURFACE *mris, MRI *mri_filled, int vno,
                           double max_dist)
{
  VERTEX  *v ;
  double  min_dist, dist, xd, yd, zd ;
  int     xoff, yoff, zoff, whalf, x, y, z, xi, yi, zi, window_size ;
  double    xw, yw, zw ;

  v = &mris->vertices[vno] ;
  if (v->ripflag)
  {
    return(0.0) ;
  }
  MRISvertexToVoxel(mris, v, mri_filled, &xw, &yw, &zw) ;
  x = nint(xw) ;
  y = nint(yw) ;
  ;
  z = nint(zw) ;

  window_size = nint(max_dist/mri_filled->xsize + 0.5) ;
  whalf = (window_size-1)/2 ;

  min_dist = 10.0*max_dist*max_dist ;
  for (zoff = -whalf ; zoff <= whalf ; zoff++)
  {
    zi = mri_filled->zi[zoff+z] ;
    zd = zi - z ;
    for (yoff = -whalf ; yoff <= whalf ; yoff++)
    {
      yi = mri_filled->yi[yoff+y] ;
      yd = yi-y ;
      for (xoff = -whalf ; xoff <= whalf ; xoff++)
      {
        xi = mri_filled->xi[xoff+x] ;
        if (MRItest_bit(mri_filled, xi, yi, zi))
        {
          xd = xi-x ;
          dist = xd*xd+yd*yd+zd*zd ;
          if (dist < min_dist)
          {
            min_dist = dist ;
          }
        }
      }
    }
  }

  min_dist = sqrt(min_dist) * mri_filled->xsize ;

  return(min_dist) ;
}
#endif


#if 0
static int
mrisComputeThicknessNormalTerm(MRI_SURFACE *mris, double l_thick_normal, INTEGRATION_PARMS *parms)
{
  int     vno, fno ;
  FACE    *face ;
  double  dxc, dyc, dzc, fdist ;
  float   dx, dy, dz, xw, yw, zw, dxn, dyn, dzn ;
  float   dxp_dxc, dyp_dxc, dzp_dxc,
          dxp_dyc, dyp_dyc, dzp_dyc,
          dxp_dzc, dyp_dzc, dzp_dzc ;
  float   e1x, e1y, e1z, e2x, e2y, e2z, norm, e1dot, e2dot ;
  float   xp, yp, zp, xpP, xpM, ypP, ypM, zpP, zpM ;
  VERTEX  *v ;

  if (FZERO(l_thick_normal))
  {
    return(0.0) ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }

    MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw) ;
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,v->x,v->y,v->z, PIAL_VERTICES, &xp, &yp, &zp);
    dxn = xp-xw ;
    dyn = yp-yw ;
    dzn = zp-zw ;
    norm = sqrt(dxn*dxn + dyn*dyn + dzn*dzn) ;
    if (!FZERO(norm))
    {
      dxn /= norm ;
      dyn /= norm ;
      dzn /= norm ;
    }
    dx = dxn-v->wnx ;
    dy = dyn-v->wny ;
    dz = dzn-v->wnz ;
    dx *= -1 ;
    dy *= -1 ;
    dz *= -1 ; // move in -gradient direction

    e1x = v->e1x ;
    e1y = v->e1y ;
    e1z = v->e1z ;
    e2x = v->e2x ;
    e2y = v->e2y ;
    e2z = v->e2z ;

#define NSAMPLE_DIST .1
    // d / dxc
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                  v->x+NSAMPLE_DIST,v->y,v->z, PIAL_VERTICES, &xpP, &ypP, &zpP);
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                  v->x-NSAMPLE_DIST,v->y,v->z, PIAL_VERTICES, &xpM, &ypM, &zpM);
    MHTfindClosestFaceGeneric((MHT *)(parms->mht), mris, v->x,v->y,v->z, 1000, -1, 1, &face, &fno, &fdist) ;
    MRISsampleFaceCoords(mris, fno,
                         v->x+NSAMPLE_DIST,v->y,v->z, WHITE_VERTICES, &xpP, &ypP, &zpP);
    MRISsampleFaceCoords(mris, fno,
                         v->x-NSAMPLE_DIST,v->y,v->z, WHITE_VERTICES, &xpM, &ypM, &zpM);
    dxp_dxc = (xpP - xpM) / (2*NSAMPLE_DIST) ;
    dyp_dxc = (ypP - ypM) / (2*NSAMPLE_DIST) ;
    dzp_dxc = (zpP - zpM) / (2*NSAMPLE_DIST) ;


    // d / dyc
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                  v->x,v->y+NSAMPLE_DIST,v->z, PIAL_VERTICES, &xpP, &ypP, &zpP);
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                  v->x,v->y-NSAMPLE_DIST,v->z, PIAL_VERTICES, &xpM, &ypM, &zpM);
    MRISsampleFaceCoords(mris, fno,
                         v->x,v->y+NSAMPLE_DIST,v->z, WHITE_VERTICES, &xpP, &ypP, &zpP);
    MRISsampleFaceCoords(mris, fno,
                         v->x,v->y-NSAMPLE_DIST,v->z, WHITE_VERTICES, &xpM, &ypM, &zpM);
    dxp_dyc = (xpP - xpM) / (2*NSAMPLE_DIST) ;
    dyp_dyc = (ypP - ypM) / (2*NSAMPLE_DIST) ;
    dzp_dyc = (zpP - zpM) / (2*NSAMPLE_DIST) ;

    // d / dzc
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                  v->x,v->y,v->z+NSAMPLE_DIST, PIAL_VERTICES, &xpP, &ypP, &zpP);
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                  v->x,v->y,v->z-NSAMPLE_DIST, PIAL_VERTICES, &xpM, &ypM, &zpM);
    MRISsampleFaceCoords(mris, fno,
                         v->x,v->y,v->z+NSAMPLE_DIST, WHITE_VERTICES, &xpP, &ypP, &zpP);
    MRISsampleFaceCoords(mris, fno,
                         v->x,v->y,v->z-NSAMPLE_DIST, WHITE_VERTICES, &xpM, &ypM, &zpM);
    dxp_dzc = (xpP - xpM) / (2*NSAMPLE_DIST) ;
    dyp_dzc = (ypP - ypM) / (2*NSAMPLE_DIST) ;
    dzp_dzc = (zpP - zpM) / (2*NSAMPLE_DIST) ;

    dxc = l_thick_normal * ((dx * dxp_dxc) + (dy * dyp_dxc) + (dz * dzp_dxc)) ;
    dyc = l_thick_normal * ((dx * dxp_dyc) + (dy * dyp_dyc) + (dz * dzp_dyc)) ;
    dzc = l_thick_normal * ((dx * dxp_dzc) + (dy * dyp_dzc) + (dz * dzp_dzc)) ;
    norm = sqrt(dxc*dxc + dyc*dyc + dzc*dzc) ;
    if (norm > MAX_NORM)
    {
      dxc = MAX_NORM * dxc / norm ;
      dyc = MAX_NORM * dyc / norm ;
      dzc = MAX_NORM * dzc / norm ;
    }
    // project onto tangent direction
    e1dot = e1x*dxc + e1y*dyc + e1z*e1z ;
    e2dot = e2x*dxc + e2y*dyc + e2z*e2z ;
    //    dxc = e1x*e1dot + e2x*e2dot ; dyc = e1y*e1dot + e2y*e2dot ; dzc = e1z*e1dot + e2z*e2dot ;
    v->dx += dxc ;
    v->dy += dyc ;
    v->dz += dzc ;
    if (Gdiag_no == vno)
    {
      double d=NSAMPLE_DIST ;

      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                    v->x+d*v->dx,v->y+d*v->dy,v->z+d*v->dz, PIAL_VERTICES,
                                    &xpP,&ypP,&zpP);
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris,
                                    v->x-d*v->dx,v->y-d*v->dy,v->z-d*v->dz, PIAL_VERTICES,
                                    &xpM,&ypM,&zpM);

      if (((xpP-xpM) * dx < 0) ||
          ((ypP-ypM) * dy < 0) ||
          ((zpP-zpM) * dz < 0))
      {
        DiagBreak() ;
      }
      printf("l_thick_normal: v %d: fno = %d, DX = (%2.2f, %2.2f, %2.2f)\n", vno, fno, dxc, dyc, dzc) ;
      printf("    N = (%2.2f, %2.2f, %2.2f), D = (%2.2f, %2.2f, %2.2f), dot= %f\n",
             v->wnx, v->wny, v->wnz, dxn, dyn, dzn,
             v->wnx*dxn + v->wny*dyn + v->wnz*dzn);
      DiagBreak() ;
    }
  }
  return(NO_ERROR) ;
}
#else
#endif


#if 0
static double
mrisComputeNegativeLogPosterior2D(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int *pnvox)
{
  MRI       *mri  = parms->mri_brain ;
  double    sse = 0.0, ll, wm_frac, gm_frac, out_frac, Ig, Ic, pval, dist ;
  float     vmin, vmax, val ;
  HISTO2D   *hin, *hout ;
  int       x, y, z, nvox, label ;
  MRI        *mri_ll = NULL ;
  static double last_sse = 0.0 ;
  HISTOGRAM2D *hs ;

  if (parms->mri_dtrans == NULL)
    MRIScomputeDistanceToSurface(mris, NULL, parms->resolution) ;

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    mri_ll = MRIcloneDifferentType(mri, MRI_FLOAT) ;
    sprintf(fname, "%s.vfrac.%4.4d.mgz", parms->base_name, parms->t) ;
    MRIwrite(parms->mri_volume_fractions, fname) ;
  }

  MRIvalRange(mri, &vmin, &vmax) ;
  if (mri->type == MRI_UCHAR)
  {
    vmin = 0 ; vmax = 255 ;
  }
  if (parms->h2d_gm == NULL)
  {
    printf("creating 2D intensity histograms\n") ;
    parms->h2d_gm = hin = HISTO2Dinit(parms->h2d_gm, 256, 256, (double)vmin, (double)vmax, -10, 10) ;
    parms->h2d_out = hout = HISTO2Dinit(parms->h2d_out, 256, 256, (double)vmin,  (double)vmax, -10, 10) ;
    MRISclearMarks(mris) ;
    
    // build histogram estimates of PDFs of interior and exterior of ribbon
    for (x = 0 ; x < mri->width; x++)
      for (y = 0 ; y < mri->height; y++)
    for (z = 0 ; z < mri->depth; z++)
    {
      if (Gx == x && Gy == y && Gz == z)
        DiagBreak() ;
      if (Gx2 == x && Gy2 == y && Gz2 == z)
        DiagBreak() ;
      val = MRIgetVoxVal(mri, x, y, z, 0) ;
      if (FZERO(val))
        continue ;

      dist = MRIgetVoxVal(parms->mri_dtrans, x, y, z, 0) ;
      wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0) ;
      gm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 1) ;
      out_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 2) ;
      if (parms->mri_aseg)
      {
        label = MRIgetVoxVal(parms->mri_aseg, x, y, z, 0) ;
        if (FZERO(gm_frac) && IS_CORTEX(label))  // aseg thinks it is but outside ribbon - ambiguous
          continue ;
      }
      HISTO2DaddFractionalSample(hout, val, dist, 0, 0, 0, 0, out_frac) ;
      HISTO2DaddFractionalSample(hin, val, dist, 0, 0, 0, 0, gm_frac) ;
    }
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      HISTOGRAM2D *hs ;

      hs = HISTO2DsoapBubbleZeros(hin, NULL, 10) ;
      sprintf(fname, "hin.%s.%3.3d.plt", parms->base_name, parms->t) ;
      printf("writing histogram %s\n", fname) ;
      HISTO2Dwrite(hin, fname) ;
      sprintf(fname, "hin.%s.%3.3d.soap.plt", parms->base_name, parms->t) ;
      printf("writing histogram %s\n", fname) ;
      HISTO2Dwrite(hs, fname) ;
      HISTO2Dfree(&hs) ;

      hs = HISTO2DsoapBubbleZeros(hout, NULL, 10) ;
      HISTO2Dfree(&hs) ;
    }

    HISTO2DmakePDF(hin, hin) ;
    HISTO2DmakePDF(hout, hout) ;
    hs = HISTO2Dsmooth(hin, NULL, 5) ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      sprintf(fname, "hin.%s.%3.3d.plt", parms->base_name, parms->t) ;
      HISTO2Dwrite(hin, fname) ; 
      sprintf(fname, "hin.%s.%3.3d.smooth.plt", parms->base_name, parms->t) ;
      HISTO2Dwrite(hs, fname) ; 
    }
    HISTO2Dfree(&hin) ; hin = parms->h2d_gm = hs ;

    hs = HISTO2Dsmooth(hout, NULL, 5) ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      sprintf(fname, "hout.%s.%3.3d.plt", parms->base_name, parms->t) ;
      HISTO2Dwrite(hout, fname) ; 
      sprintf(fname, "hout.%s.%3.3d.smooth.plt", parms->base_name, parms->t) ;
      HISTO2Dwrite(hs, fname) ; 
    }
    HISTO2Dfree(&hout) ; hout = parms->h2d_out = hs ;
    

    if (Gx > 0)
    {
      int b1, b2 ;
      float wm_frac, gm_frac, out_frac, val, cin, cout ;

      x = Gx ; y = Gy ; z = Gz ;
      dist = MRIgetVoxVal(parms->mri_dtrans, x, y, z, 0) ;
      wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0) ;
      gm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 1) ;
      out_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 2) ;
      val = MRIgetVoxVal(mri, x, y, z, 0) ;
      b1 = HISTO2DfindBin1(hin, val) ;
      b2 = HISTO2DfindBin2(hin, dist) ;
      cin = HISTO2DgetCount(hin, val, dist) ;
      cout = HISTO2DgetCount(hout, val, dist) ;
      printf("voxel (%d, %d, %d) - vfracts (%2.1f, %2.1f, %2.1f) = %2.0f, dist=%2.1f, bin=%d,%d, cin=%f, cout=%f\n",
         Gx, Gy, Gz, wm_frac, gm_frac, out_frac, val, dist, b1, b2, cin, cout) ;
      
      DiagBreak() ;
    }
  }
  else // use previously computed ones
  {
    hin = parms->h2d_gm ;
    hout = parms->h2d_out ;
  }

  for (sse = 0.0, nvox = x = 0 ; x < mri->width; x++)
    for (y = 0 ; y < mri->height; y++)
      for (z = 0 ; z < mri->depth; z++)
      {
    if (Gx == x && Gy == y && Gz == z)
      DiagBreak() ;
    if (Gx2 == x && Gy2 == y && Gz2 == z)
      DiagBreak() ;
    val = MRIgetVoxVal(mri, x, y, z, 0) ;
    if (FZERO(val))
      continue ;
    wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0) ;
    if (wm_frac > 0)
      continue ;

    gm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 1) ;
    out_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 2) ;
    dist = MRIgetVoxVal(parms->mri_dtrans, x, y, z, 0) ;

    if (FZERO(out_frac))   // all gm
      pval = HISTO2DgetCount(hin, val, dist) ;
    else if (FZERO(gm_frac))
      pval = HISTO2DgetCount(hout, val, dist) ;
    else  // partial volume voxel
    {
      for (pval = 0.0, Ig = 0 ; Ig <= 256 ; Ig++)
      {
        Ic = (val - gm_frac*Ig) / out_frac ;
        pval += HISTO2DgetCount(hout, Ic, dist) * HISTO2DgetCount(hin, Ig, dist) ;
      }
    }
    
    if (pval > 1 || pval < 0)
      DiagBreak() ;
    if (DZERO(pval))
      pval = 1e-20 ;
    ll = -log10(pval) ;
    sse += ll ;
    nvox++ ;
    if (mri_ll)
      MRIsetVoxVal(mri_ll, x, y, z, 0, ll) ;
    if (Gx == x && Gy == y && Gz == z)
      printf("voxel(%d, %d, %d) = %d, vfracs = (%2.1f, %2.1f, %2.1f), ll = %2.1f\n",
         x, y, z, nint(val), wm_frac, gm_frac, out_frac, ll) ;
      }

  if (!FZERO(last_sse) && sse > last_sse)
    DiagBreak() ;
  if (mri_ll)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.ll.%4.4d.mgz", parms->base_name, parms->t) ;
    printf("writing log likelihood volume to %s\n", fname) ;
    MRIwrite(mri_ll, fname) ;
    MRIfree(&mri_ll) ;
  }
//  HISTOfree(&hin) ; HISTOfree(&hout) ; 
  if (pnvox)
    *pnvox = nvox ;
  if (Gdiag_no >= 0)
    printf("E_logPosterior, %3.3d: %.1f (nvox=%d)\n", parms->t, sse, nvox) ;

  last_sse = sse ; // diagnostics
  return(sse) ;
}

static int
mrisComputePosterior2DTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  MRI       *mri  = parms->mri_brain ;
  static MRI *mri_white_dist ;
  FILE      *fp ;
  double    dist, xs, ys, zs, dn, best_dn, best_ll, ll ;
  float     vmin, vmax, val, wdist ;
  HISTO2D   *hin, *hout ;
  int       x, y, z, vno, n ;
  VECTOR    *v1, *v2 ;
  MATRIX    *m_vox2vox ;
  MHT       *mht ;
  VERTEX    *v, *vn ;
  VOXEL_LIST **vl, **vl2 ;

  if (FZERO(parms->l_map2d))
    return(NO_ERROR) ;
  vl = vlst_alloc(mris, MAX_VOXELS) ;
  vl2 = vlst_alloc(mris, MAX_VOXELS) ;

  hin = parms->h2d_gm ; hout = parms->h2d_out ;  // created in mrisComputeNegativeLogPosterior2D

  if (mri_white_dist == NULL)
  {
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
    mri_white_dist = MRIScomputeDistanceToSurface(mris, NULL, 0.5) ;
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  }
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 10);

  m_vox2vox = MRIgetVoxelToVoxelXform(mri, mri_white_dist) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ; v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1,4) = 1.0; VECTOR_ELT(v2,4) = 1.0;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_white_dist, "wd.mgz") ;

  MRIvalRange(mri, &vmin, &vmax) ;
  if (mri->type == MRI_UCHAR)
  {
    vmin = 0 ; vmax = 255 ;
  }
  MRISclearMarks(mris) ;

  // build histogram estimates of PDFs of interior and exterior of ribbon
  for (x = 0 ; x < mri->width; x++)
    for (y = 0 ; y < mri->height; y++)
      for (z = 0 ; z < mri->depth; z++)
      {
    val = MRIgetVoxVal(mri, x, y, z, 0) ;
    if (FZERO(val))
      continue ;
    V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
    MatrixMultiply(m_vox2vox, v1, v2) ;

    if (MRIindexNotInVolume(mri_white_dist, V3_X(v2), V3_Y(v2), V3_Z(v2)))
      continue ;
    wdist = MRIgetVoxVal(mri_white_dist, V3_X(v2), V3_Y(v2), V3_Z(v2), 0) ;
    if (wdist < 0)
      continue ;

    // add this voxel to the list of voxels of the vertex it is closest to
    MRIvoxelToSurfaceRAS(mri, x, y, z, &xs, &ys, &zs) ;
    MHTfindClosestVertexGeneric(mht, mris,  xs, ys, zs, 10, 4, &v, &vno, &dist) ;
    if (v == NULL)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->marked++ ;
    VLSTadd(vl[vno], x, y, z, xs, ys, zs) ;
      }

  // find vertices that don't have enough data and pool across nbrs
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vlst_add_to_list(vl[vno], vl2[vno]) ;
    if (v->ripflag || vlst_enough_data(mris, vno, vl[vno], 0.0))
      continue ;
    for (n = 0; n < v->vnum ; n++)
    {
      if (vl2[vno]->nvox + vl[v->v[n]]->nvox >= vl2[vno]->max_vox)
    break ;
      vlst_add_to_list(vl[v->v[n]], vl2[vno]) ;
    }
    v->marked = vl[vno]->nvox ;
  }

  vlst_free(mris, &vl) ; vl = vl2 ;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri, parms->mri_dtrans) ; // dist to pial surface
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag || v->marked == 0)
      continue ;
    if (vlst_enough_data(mris, vno, vl[vno], 0.0) == 0)
    {
      v->marked = 0 ;
      continue ;
    }

    if (vno == Gdiag_no)
    {
      char fname[STRLEN] ;
      sprintf(fname, "vno%d.%d.l.dat", vno, parms->t) ;
      fp = fopen(fname, "w") ;
    }
    else
      fp = NULL ;
    best_ll = -1e10 ; best_dn = 0 ;
    for (dn = -MAX_DISPLACEMENT ; dn <= MAX_DISPLACEMENT ; dn += DISPLACEMENT_DELTA)
    {
      if (fp)
    fprintf(fp, "%f ", dn) ;
      ll = -vlst_loglikelihood2D(mris, mri, vno, dn, vl[vno], hin, hout, fp) ;
      if (fp)
    fprintf(fp, "\n") ;
      if (devIsnan(ll))
    DiagBreak() ;
      if (fp)
    fprintf(fp, "%f %f\n", dn, ll) ;
      if (ll > best_ll)
      {
    best_dn = dn ;
    best_ll = ll ;
      }
    }
    if (fp)
      fclose(fp) ;

//#undef MAX_MOVEMENT
//#define MAX_MOVEMENT 1

#if 1
    if (fabs(best_dn) > MAX_MOVEMENT)
      best_dn = MAX_MOVEMENT*best_dn / fabs(best_dn) ;
#endif
    if (vno == Gdiag_no)
    {
      int i ;
      char fname[STRLEN] ;
      double dx, dy, dz, dist, dot, pin, pout ;

      sprintf(fname, "vno%d.%d.vox.l.dat", vno, parms->t) ;
      fp = fopen(fname, "w") ;

      for (i = 0 ; i < vl[vno]->nvox ; i++)
      {
    val = MRIgetVoxVal(mri, vl[vno]->xi[i], vl[vno]->yi[i], vl[vno]->zi[i], 0) ;
    dx = vl[vno]->xd[i]-v->x ; dy = vl[vno]->yd[i]-v->y ; dz = vl[vno]->zd[i]-v->z ;
    dist = dx*dx + dy*dy + dz*dz ;
    dot = dx*v->nx + dy*v->ny + dz*v->nz ;
    if (dot <0)
      dist *= -1 ;
    pin = HISTO2DgetCount(hin, val, dot) ; pout = HISTO2DgetCount(hout, val, dot) ;
    fprintf(fp, "%d %d %d %d %d %f %f %f %f\n",
        vno, i, vl[vno]->xi[i], vl[vno]->yi[i], vl[vno]->zi[i], val, dist, pin, pout) ;
      }
      fclose(fp) ;

      printf("l_map: vno %d, best displacement %2.3f, ll = %2.3f, D = (%2.3f, %2.3f, %2.3f)\n",  
         vno, best_dn, best_ll, 
         best_dn * v->nx * parms->l_map2d,
         best_dn * v->ny * parms->l_map2d,
         best_dn * v->nz * parms->l_map2d) ;
      DiagBreak() ;
    }
    v->dx += best_dn * v->nx * parms->l_map2d ;
    v->dy += best_dn * v->ny * parms->l_map2d ;
    v->dz += best_dn * v->nz * parms->l_map2d ;
    v->d = best_dn ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    int    num ;
    double dn = 0.0 ;

    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->marked > 0)  // already estimated
      continue ;
    for (n = num = 0; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked == 0)
    continue ;
      num++ ;
      dn += vn->d ;
    }
    if (num > 0)
    {
      dn /= num ;
      v->dx += dn * v->nx * parms->l_map2d ;
      v->dy += dn * v->ny * parms->l_map2d ;
      v->dz += dn * v->nz * parms->l_map2d ;
      v->d = dn ;
      if (vno == Gdiag_no)
    printf("l_map: vno %d, soap bubble displacement %2.3f, D = (%2.3f, %2.3f, %2.3f)\n",  
           vno, dn, dn * v->nx * parms->l_map2d, dn * v->ny * parms->l_map2d,
           dn * v->nz * parms->l_map2d) ;
    }
  }

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN], path[STRLEN] ;

    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s.%d.dist.mgz", path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : 
        mris->hemisphere == BOTH_HEMISPHERES ? "both" : "rh", parms->t) ;
    MRISwriteD(mris, fname) ;
    DiagBreak() ;
  }
  MHTfree(&mht) ;
//  HISTOfree(&hin) ; HISTOfree(&hout) ; 
  VectorFree(&v1) ; VectorFree(&v2) ; MatrixFree(&m_vox2vox) ;
  vlst_free(mris, &vl) ;
  return(NO_ERROR) ;
}
#endif


#if 0
static double
vlst_loglikelihood2D(MRI_SURFACE *mris, MRI *mri, int vno, double displacement, VOXEL_LIST *vl, HISTOGRAM2D *hin, HISTOGRAM2D *hout)
{
  double   ll = 0.0, dot, dx, dy, dz, pval, dist, Ig, Ic, gm_frac, out_frac ;
  int      i ;
  float    val ;
  VERTEX   *v ;
  double     xs, ys, zs ;

  v = &mris->vertices[vno] ;
  xs = v->x + displacement*v->nx ; ys = v->y + displacement*v->ny ; zs = v->z + displacement*v->nz ;
  for (i = 0 ; i < vl->nvox ; i++)
  {
    dx = vl->xd[i]-xs ; dy = vl->yd[i]-ys ; dz = vl->zd[i]-zs ; 
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    dot = dx*v->nx + dy*v->ny + dz*v->nz ;
    val = MRIgetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0) ;
    if (dist < .5) // distance to center<.5 --> distance to edge <1
    {
      if (dot > 0)
      {
    out_frac = dist+.5 ;
    gm_frac = 1-out_frac ;
      }
      else
      {
    gm_frac = dist+.5 ;
    out_frac = 1-gm_frac ;
      }
      for (pval = 0.0, Ig = 0 ; Ig <= 256 ; Ig++)
      {
    Ic = (val - gm_frac*Ig) / out_frac ;
    pval += HISTO2DgetCount(hout, Ic, dot) * HISTO2DgetCount(hin, Ig, dot) ;
      }
    }
    else if (dot > 0)  // outside surface
      pval = HISTO2DgetCount(hout, val, dot) ;
    else         // inside the surface
      pval = HISTO2DgetCount(hin, val, dot) ;
    if (DZERO(pval))
      pval = 1e-10 ;
    ll += -log(pval) ;
  }
  
  return(ll) ;
}
#else
#endif

#if 0
/////////////// not used ////////////////////////////////////////////////
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Find the 3 poles (temporal, occipital, and frontal) of
  the cortical surface.
  ------------------------------------------------------*/
#define MIN_Z_DISTANCE 30.0f
#define MIN_Y_DISTANCE 30.0f
#define MIN_ANGLE_VARIATION RADIANS(30.0f)

static int
mrisFindPoles(MRIS *mris)
{
  int     vno, n, neigh, local_max ;
  VERTEX  *vertex, *v_neigh ;
  double    x, y, z, xt, yt, zt ;
  float   temporal_y_hi = -1000, zfront, yfront, angle ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "(%2.0f, %2.0f, %2.0f) --> (%2.0f, %2.0f, %2.0f), ctr "
            "(%2.0f, %2.0f, %2.0f)\n",
            mris->xlo, mris->ylo, mris->zlo, mris->xhi, mris->yhi, mris->zhi,
            mris->xctr, mris->yctr, mris->zctr);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    fprintf(stdout, "finding cortical poles...") ;
  }

  /* first find frontal and occipital poles */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->y >= mris->yhi)
    {
      mris->v_frontal_pole = vertex ;
    }
    else if (vertex->y <= mris->ylo)
    {
      mris->v_occipital_pole = vertex ;
    }
  }

  /* now find temporal pole */
  zfront = mris->v_frontal_pole->z ;
  yfront = mris->v_frontal_pole->y ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = (double)vertex->x ;
    y = (double)vertex->y ;
    z = (double)vertex->z ;
    // get the talairach point (RAS)
    if (mris->lta)
    {
      TransformWithMatrix(mris->SRASToTalSRAS_, x, y, z, &xt, &yt, &zt);
    }
    else
    {
      // assume that talRAS and RAS are the same
      xt = x;
      yt = y;
      zt = z;
    }
#if 0
    // I don't understand the reasoning of flip x, exchange y and z
    // since surface vertex is already a RAS.
    if (mris->linear_transform)
    {
      transform_point(mris->linear_transform, -x, z, y, &xt, &yt, &zt) ;
    }
    else
    {
      xt = -x ;
      yt = z ;
      zt = y ;
    }
#endif
    /*
      some rules for finding the temporal pole:

      1. must be a local max in y (posterior-anterior) direction.    y <-- A
      2. Must be below some absolute y talairach coordinate.
      3. Must be a minimum distance from the frontal pole
      in y and z directions.
      4. Must have a normal vector within 30 degrees of (0,1,0).
    */
    if ((yt < MAX_TALAIRACH_Y) &&  /* low enough to be temporal pole */
        ((zfront - vertex->z) > MIN_Z_DISTANCE) &&
        ((yfront - vertex->y) > MIN_Y_DISTANCE))
    {
      local_max = 1 ;
      if (vertex->y > temporal_y_hi)  /* check neighbors positions */
      {
        for (n = 0 ; n < vertex->vnum ; n++)
        {
          neigh = vertex->v[n] ;
          v_neigh = &mris->vertices[neigh] ;
          if (v_neigh->y > vertex->y)
          {
            local_max = 0 ;
            break ;
          }
        }

        angle = acos(vertex->ny) ;
        if (local_max && (angle < MIN_ANGLE_VARIATION))
        {
          mris->v_temporal_pole = vertex ;
          temporal_y_hi = vertex->y ;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    if (mris->v_temporal_pole)
      fprintf(stdout, "F: (%2.0f,%2.0f,%2.0f), T: (%2.0f,%2.0f,%2.0f) "
              "O: (%2.0f,%2.0f,%2.0f).\n",
              mris->v_frontal_pole->x, mris->v_frontal_pole->y,
              mris->v_frontal_pole->z,
              mris->v_temporal_pole->x, mris->v_temporal_pole->y,
              mris->v_temporal_pole->z,
              mris->v_occipital_pole->x, mris->v_occipital_pole->y,
              mris->v_occipital_pole->z) ;
    else
      fprintf(stdout, "F: (%2.0f,%2.0f,%2.0f), T: (NOT FOUND), "
              "O: (%2.0f,%2.0f,%2.0f).\n",
              mris->v_frontal_pole->x, mris->v_frontal_pole->y,
              mris->v_frontal_pole->z,
              mris->v_occipital_pole->x, mris->v_occipital_pole->y,
              mris->v_occipital_pole->z) ;
  }
  return(NO_ERROR) ;
}
///////////// not used ////////////////////////////////////////////////////
#endif

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0

#if 0
static float area_scoefs[] =
{
  1.0f,    1.0f,   1.0f,  1.0f, 1.0f, 0.1f, 0.01f, 0.001f
};
static float dist_scoefs[] =
{
  0.0001f, 0.001f, 0.01f, 0.1f, 1.0f, 1.0f, 1.0f,  1.0f
} ;
#else
static float area_scoefs[] =
{
  1.0f,  0.1f, 0.01f
} ;
static float dist_scoefs[] =
{
  1.0f,  1.0f, 1.0f
} ;
#endif

#define NSCOEFS sizeof(area_scoefs) / sizeof(area_scoefs[0])

MRI_SURFACE *
MRISunfoldOnSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int max_passes)
{
  int     base_averages, i, nbrs[MAX_NBHD_SIZE], niter, passno, msec ;
  double  starting_sse, ending_sse ;
  struct  timeb start ;

  starting_sse = ending_sse = 0.0f ;   /* compiler warning */
  memset(nbrs, 0, MAX_NBHD_SIZE*sizeof(nbrs[0])) ;
  for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
  {
    nbrs[i] = parms->max_nbrs ;
  }

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s.out", parms->base_name) ;
    if (!parms->fp)
    {
      if (!parms->start_t)
      {
        parms->fp = fopen(fname, "w") ;
      }
      else
      {
        parms->fp = fopen(fname, "a") ;
      }
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
                  Progname, fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris,parms) ;
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
      {
        fprintf(parms->fp, "%d: %d | ", i, nbrs[i]) ;
      }
    fprintf(parms->fp, "\n") ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      if (nbrs[i])
      {
        fprintf(stdout, "%d: %d | ", i, nbrs[i]) ;
      }
    fprintf(stdout, "\n") ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    mrisLogIntegrationParms(stderr, mris, parms) ;
  }

  /*
    integrate until no improvement can be made at ANY scale, or until
    the error has asymptoted.
  */
  base_averages = parms->n_averages ;
  niter = parms->niterations ;
  passno = 0 ;
  mrisProjectSurface(mris) ;
  MRIScomputeMetricProperties(mris) ;
  do
  {
    if (passno++ >= max_passes)
    {
      break ;
    }

    /* first time through only - use big ratio to remove folds */
    TimerStart(&start) ;
    for (i = 0 ; i < NSCOEFS ; i++)
    {
      if (mris->nsize < parms->nbhd_size)  /* resample distances on surface */
      {
        if (Gdiag & DIAG_SHOW)
        {
          fprintf(stdout, "resampling long-range distances") ;
        }
        MRISsaveVertexPositions(mris, TMP_VERTICES) ;
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
        MRISsampleDistances(mris, nbrs, parms->nbhd_size) ;
        MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
        mrisClearMomentum(mris) ;
      }
      if (!i)  /* record starting error */
      {
        parms->l_nlarea = area_scoefs[NSCOEFS-1] ;
        parms->l_dist = dist_scoefs[NSCOEFS-1] ;
        starting_sse = MRIScomputeSSE(mris, parms) ;
      }

      /* remove any folds in the surface */
      mrisRemoveNegativeArea(mris, parms, base_averages,
                             MAX_NEG_AREA_PCT, 2) ;

      parms->l_dist = dist_scoefs[i] ;
      parms->l_nlarea = area_scoefs[i] ;
      parms->l_angle = ANGLE_AREA_SCALE * area_scoefs[i] ;
      mrisIntegrationEpoch(mris, parms, base_averages) ;
    }

    parms->l_area = area_scoefs[NSCOEFS-1] ;
    parms->l_dist = dist_scoefs[NSCOEFS-1] ;
    ending_sse = MRIScomputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
    {
#if 0
      fprintf(stdout, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
              passno, starting_sse, ending_sse,
              (starting_sse-ending_sse)/ending_sse) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,
                "pass %d: start=%2.4f, end=%2.4f, ratio=%2.4f\n",
                passno, starting_sse, ending_sse,
                (starting_sse-ending_sse)/starting_sse) ;
#endif
    }
    msec = TimerStop(&start) ;
    if (Gdiag & DIAG_SHOW)
    {
      fprintf(stdout,
              "epoch took %2.2f minutes\n",(float)msec/(1000.0f*60.0f));
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,
                "epoch took %2.2f minutes\n",(float)msec/(1000.0f*60.0f));
    }
  }
  while (!FZERO(ending_sse) &&
         (((starting_sse-ending_sse)/starting_sse) > parms->tol)) ;

#if 0
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout,
            "initial unfolding complete - settling to equilibrium...\n") ;
  parms->niterations = 1000 ;  /* let it go all the way to equilibrium */
  mrisIntegrationEpoch(mris, parms, parms->n_averages = 0) ;
#endif

  /* finally, remove all the small holes */
  parms->l_nlarea = 1.0f ;
  parms->l_area = 0.0 ;
  parms->l_dist = 0.001f ;
  parms->l_angle = ANGLE_AREA_SCALE * area_scoefs[0] ;
  parms->niterations = niter ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "removing remaining folds...\n") ;
  }
  mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 3);
  if (Gdiag & DIAG_SHOW)
  {
    mrisLogStatus(mris, parms, stderr, 0,-1) ;
  }
  if (Gdiag & DIAG_WRITE)
  {
    mrisLogStatus(mris, parms, parms->fp, 0,-1) ;
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "unfolding complete.\n") ;
  }
  if (Gdiag & DIAG_WRITE)
  {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  mrisProjectSurface(mris) ;
  return(mris) ;
}
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if AVERAGE_AREAS
#error AVERAGE_AREAS assumed obsolete - move this to mrisurf_deform.c

static int mrisAverageAreas(MRI_SURFACE *mris, int num_avgs, int which)
{
  int i, vno, vnb, *pnb, fno, num, vf, nfno;
  float area;
  VERTEX *v, *vn;
  FACE *f;

  for (i = 0; i < num_avgs; i++) {
    switch (which) {
      case ORIG_AREAS:
        for (vno = 0; vno < mris->nvertices; vno++) {
          v = &mris->vertices[vno];
          if (v->ripflag) {
            continue;
          }

          for (fno = 0; fno < v->num; fno++) /* each face of this vertex */
          {
            f = &mris->faces[v->f[fno]]; /* pointer to the face */
            if (f->ripflag) {
              continue;
            }
            area = 0.0f;

            /* now go through each vertex associated with this face */
            for (vf = 0; vf < VERTICES_PER_FACE; vf++) {
              vn = &mris->vertices[f->v[vf]];
              num += vn->num;
              for (nfno = 0; nfno < vn->num; nfno++) {
                area += vn->orig_tri_area[nfno];
              }
            }
            area /= (float)num;
            v->orig_tri_area[fno] = area;
          }
        }
        break;
      case CURRENT_AREAS:
        for (vno = 0; vno < mris->nvertices; vno++) {
          v = &mris->vertices[vno];
          if (v->ripflag) {
            continue;
          }

          for (fno = 0; fno < v->num; fno++) /* each face of this vertex */
          {
            f = &mris->faces[v->f[fno]]; /* pointer to the face */
            if (f->ripflag) {
              continue;
            }
            area = 0.0f;

            /* now go through each vertex associated with this face */
            for (vf = 0; vf < VERTICES_PER_FACE; vf++) {
              vn = &mris->vertices[f->v[vf]];
              num += vn->num;
              for (nfno = 0; nfno < vn->num; nfno++) {
                area += vn->tri_area[nfno];
              }
            }
            area /= (float)num;
            v->tri_area[fno] = area;
          }
        }
        break;
    }
  }
  return (NO_ERROR);
}
#endif

#if 1
#else
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteD(MRI_SURFACE *mris, const char *sname)
{
  int k, num; /* loop counters */
  float f;
  char fname[STRLEN], *cp;
  FILE *fp;
  double sum = 0, sum2 = 0, max = -1000, min = 1000;
  int ftype, err;
  MRI *TempMRI;

  MRISbuildFileName(mris, sname, fname);

  // Try saving it in a "volume" format -- but not img or nifti
  // as they use shorts for the number of vertices. Should add
  // a reshape.
  ftype = mri_identify(sname);
  if (ftype != MRI_VOLUME_TYPE_UNKNOWN) {
    TempMRI = MRIcopyMRIS(NULL, mris, 0, "d");
    if (TempMRI == NULL) {
      printf("ERROR: MRISwriteD: could not alloc MRI\n");
      return (1);
    }
    printf("Saving surf d vals in 'volume' format to %s\n", sname);
    err = MRIwrite(TempMRI, sname);
    return (err);
  }

  cp = strrchr(fname, '.');
  if (!cp || *(cp + 1) != 'w') {
    strcat(fname, ".w");
  }
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "writing surface d values to %s.\n", fname);
  }

  fp = fopen(fname, "wb");
  if (fp == NULL) {
    ErrorExit(ERROR_NOFILE, "Can't create file %s\n", fname);
  }

  for (k = 0, num = 0; k < mris->nvertices; k++)
    if (mris->vertices[k].d != 0) {
      num++;
    }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("num = %d\n", num);
  }
  fwrite2(0, fp);
  fwrite3(num, fp);

  for (k = 0; k < mris->nvertices; k++) {
    if (mris->vertices[k].d != 0) {
      fwrite3(k, fp);
      f = mris->vertices[k].d;
      if (!isfinite(f)) ErrorPrintf(ERROR_BADPARM, "MRISwriteD(%s): v->d at vertex %d is not finite", fname, k);

      fwriteFloat(f, fp);
      sum += f;
      sum2 += f * f;
      if (f > max) {
        max = f;
      }
      if (f < min) {
        min = f;
      }
    }
  }
  fclose(fp);

  if (num > 0) {
    sum /= num;
    sum2 = (sum2 / num - sum * sum);
    if (sum2 > 0) {
      sum2 = sqrt(sum2);
    }
    else {
      sum2 = 0;
    }
    printf("avg = %2.3f, stdev = %2.3f, min = %2.3f, max = %2.3f\n", sum, sum2, min, max);
  }
  else {
    printf("Warning: all vertex d values are zero\n");
  }

  return (NO_ERROR);
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisRipVertices(MRI_SURFACE *mris)
{
  int     fno, n ;
  VERTEX  *v ;
  FACE    *f ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag == 0)
    {
      continue ;
    }
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;
      v->ripflag = 1 ;
    }
  }
  return(NO_ERROR) ;
}
#endif

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
static int
mrisStoreCurrentGradient(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    v->odx = v->dx ;
    v->ody = v->dy ;
    v->odz = v->dz ;
  }
  return(NO_ERROR) ;
}
#endif

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
static int
mrisComputeBoundaryTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
#if 0
  int    vno ;
  VERTEX *v ;
  double l_boundary ;

  l_boundary = parms->l_boundary ;

  if ((mris->status != MRIS_PLANE) || FZERO(l_boundary))
  {
    return(NO_ERROR) ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || !v->border)
    {
      continue ;
    }
    if (v->neg)
    {
#if 0
      v->dx -= l_boundary * v->bnx ;
      v->dy -= l_boundary * v->bny ;
#endif
    }
    else
    {
#if 1
      v->dx += l_boundary * v->bnx ;
      v->dy += l_boundary * v->bny ;
#endif
    }
  }

#endif
  return(NO_ERROR) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisComputeNegTerm(MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  int    vno, n, neg ;
  VERTEX *v, *vn ;
  double l_neg, dx, dy, dz, len ;

  l_neg = parms->l_neg ;

  if (FZERO(l_neg))
  {
    return(NO_ERROR) ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (v->neg)
    {
      neg = 1 ;
    }
    else
    {
      neg = 0 ;
    }
    dx = dy = dz = 0.0 ;
    for (n = 0 ; n < v->vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->neg)
      {
        neg++ ;
        continue ;
      }
      dx += vn->x - v->x ;
      dy += vn->y - v->y ;
      dz += vn->z - v->z ;
    }
    len = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (!FZERO(len) && neg > 0)
    {
      dx /= len ;
      dy /= len ;
      dz /= len ;
      v->dx += l_neg*dx ;
      v->dy += l_neg*dy ;
      v->dz += l_neg*dz ;
    }
  }

  return(NO_ERROR) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRISsoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs, float pct_fixed)
{
  int    i, vno, vnb, *pnb, vnum, j ;
  float  x, y, z, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++)
  {
    MRISclearMarks(mris) ;
    MRISmarkRandomVertices(mris, pct_fixed) ;
    for (j = 0 ; j < 10 ; j++)
    {
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag || v->marked)
        {
          continue ;
        }
        x = v->x ;
        y = v->y ;
        z = v->z ;
        pnb = v->v ;
        vnum = v->vnum ;
        for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++)
        {
          vn = &mris->vertices[*pnb++] ;/*neighboring vertex pointer */
          if (vn->ripflag)
          {
            continue ;
          }
          num++ ;
          x += vn->x ;
          y += vn->y ;
          z += vn->z ;
        }
        num++ ;   /* account for central vertex */
        v->tdx = x / num ;
        v->tdy = y / num ;
        v->tdz = z / num ;
      }
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag || v->marked)
        {
          continue ;
        }
        v->x = v->tdx ;
        v->y = v->tdy ;
        v->z = v->tdz ;
      }
    }
  }
  return(NO_ERROR) ;
}
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
#define WHALF (5 - 1) / 2
static int
mrisComputeWhiteSurfaceValues(MRI_SURFACE *mris, MRI *mri_brain,
                              MRI *mri_wm, float nsigma)
{
  double    val, x, y, z ;
  int     total_vertices, vno, xv, yv, zv, xo, yo, zo, xi, yi, zi, nvox ;
  float   total, total_sq, sigma, mean_wm, mean_gray, mean ;
  VERTEX  *v ;

  /* first compute intensity of local gray/white boundary */
  mean_wm = mean_gray = 0.0f ;

  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    v->d = 0.0f ;
    MRISvertexToVoxel(mris, v, mri_wm, &x, &y, &z) ;
    xv = nint(x) ;
    yv = nint(y) ;
    zv = nint(z) ;

    /* compute mean and variance of wm in a neighborhood of voxel */
    total = total_sq = 0.0f ;
    nvox = 0 ;
    for (zo = zv-WHALF ; zo <= zv + WHALF ; zo++)
    {
      zi = mri_wm->zi[zo] ;
      for (yo = yv-WHALF ; yo <= yv + WHALF ; yo++)
      {
        yi = mri_wm->yi[yo] ;
        for (xo = xv-WHALF ; xo <= xv + WHALF ; xo++)
        {
          xi = mri_wm->xi[xo] ;
          val = (double)MRIvox(mri_wm, xi, yi, zi) ;
          if (val > WM_MIN_VAL)
          {
#if 0
            if (MRIneighborsOff(mri_wm, xi, yi, zi, WM_MIN_VAL) == 0)
            {
              continue ;  /* not a border voxel */
            }
#endif
            val =
              (double)MRIvox(mri_brain, xi, yi, zi) ;
            /* use smoothed val */
            total += val ;
            total_sq += val * val ;
            nvox++ ;
          }
        }
      }
    }
    if (!nvox)
    {
      v->val = 0.0f ;
    }
    else
    {
      mean = total / (float)nvox ;
      sigma = sqrt(total_sq / (float)nvox - mean*mean) ;
      MRISvertexToVoxel(mris, v, mri_wm, &x, &y, &z) ;
      MRIsampleVolume(mri_brain, x, y, z, &val) ;
      v->val = mean - nsigma * sigma ;
      mean_gray += v->val ;
      mean_wm += mean ;
      total_vertices++ ;
    }

  }
  mean_wm /= (float)total_vertices ;
  mean_gray /= (float)total_vertices ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (FZERO(v->val))   /* no border voxels nearby */
    {
      v->val = mean_gray ;
    }
  }

#if 0
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    v->val = mean_gray ;
  }
#endif

  fprintf(stdout, "mean wm=%2.1f, gray=%2.1f, averaging targets and "
          "smoothing surface...", mean_wm, mean_gray) ;
  return(NO_ERROR) ;
}

#else

#if 1

#else
/*-----------------------------------------------------
Parameters:

Returns value:

Description
------------------------------------------------------*/
int MRIScomputeWhiteSurfaceValues(MRI_SURFACE *mris, MRI *mri_brain, MRI *mri_wm, float nsigma)
{
  double val, x, y, z;
  int total_vertices, vno, xv, yv, zv, xo, yo, zo, xi, yi, zi, nwhite_vox, ngray_vox;
  float total_white, total_gray, total_sq_white, total_sq_gray, std_white, mean_white, mean_gray, total_mean_gray,
      total_mean_white, std_gray;
  VERTEX *v;

  /* first compute intensity of local gray/white boundary */
  total_vertices = 0;
  total_mean_white = total_mean_gray = 0.0f;
  total_sq_white = total_sq_gray = 0.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->d = 0.0f;
    MRISvertexToVoxel(mris, v, mri_wm, &x, &y, &z);
    xv = nint(x);
    yv = nint(y);
    zv = nint(z);

    /* compute mean and variance of white in a neighborhood of voxel */
    total_white = 0.0f;
    nwhite_vox = 0;
    total_gray = 0.0f;
    ngray_vox = 0;
    for (zo = zv - WHALF; zo <= zv + WHALF; zo++) {
      zi = mri_wm->zi[zo];
      for (yo = yv - WHALF; yo <= yv + WHALF; yo++) {
        yi = mri_wm->yi[yo];
        for (xo = xv - WHALF; xo <= xv + WHALF; xo++) {
          xi = mri_wm->xi[xo];
          val = (double)MRIvox(mri_wm, xi, yi, zi);
          if (val > WM_MIN_VAL) /* white matter */
          {
            if (MRIneighborsOff(mri_wm, xi, yi, zi, WM_MIN_VAL) == 0) {
              continue; /* not a border voxel */
            }
            val = (double)MRIvox(mri_brain, xi, yi, zi); /* use smoothed
              val */
            total_white += val;
            total_sq_white += val * val;
            nwhite_vox++;
          }
          else /* gray matter */
          {
            if (MRIneighborsOn(mri_wm, xi, yi, zi, WM_MIN_VAL + 1) == 0) {
              continue; /* not a border voxel */
            }
            val = (double)MRIvox(mri_brain, xi, yi, zi); /* use
          smoothed
          val */
            total_gray += val;
            total_sq_gray += val * val;
            ngray_vox++;
          }
        }
      }
    }
    if (!nwhite_vox || !ngray_vox) {
      v->val = 0.0f;
    }
    else {
      mean_white = total_white / (float)nwhite_vox;
      mean_gray = total_gray / (float)ngray_vox;
      std_white = sqrt(total_sq_white / (float)nwhite_vox - mean_white * mean_white);
      std_gray = sqrt(total_sq_gray / (float)ngray_vox - mean_gray * mean_gray);
      std_white = 0.0f;           /* only use white matter */
      if (mean_gray > mean_white) /* shouldn't happen */
      {
        if (DIAG_VERBOSE_ON)
          fprintf(stdout, "mean gray (%2.1f) > mean white (%2.1f) at v %d!\n", mean_gray, mean_white, vno);
        v->val = mean_white;
      }
      else
        v->val = (std_gray * mean_white + std_white * mean_gray) / (std_gray + std_white);
      total_vertices++;
      total_mean_white += mean_white;
      total_mean_gray += mean_gray;
    }
    v->mean = 20.0f;
    if (vno == Gdiag_no) fprintf(stdout, "v %d, target value = %2.1f, mag = %2.1f\n", Gdiag_no, v->val, v->mean);
  }
  total_mean_white /= (float)total_vertices;
  total_mean_gray /= (float)total_vertices;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (FZERO(v->val)) /* no border voxels nearby */
    {
      v->val = (total_mean_gray + total_mean_white) / 2;
    }
  }

#if 0
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    v->val = mean_gray ;
  }
#endif

  fprintf(stdout, "mean white=%2.1f, gray=%2.1f\n", total_mean_white, total_mean_gray);
  return (NO_ERROR);
}
#endif
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
static int
mrisNeighborAtVoxel(MRI_SURFACE *mris, MRI *mri, int vno, int xv,int yv,int zv)
{
  int      n, xnv, ynv, znv ;
  VERTEX   *v, *vn ;
  double     xn, yn, zn ;

  v = &mris->vertices[vno] ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    MRISvertexToVoxel(mris, v, mri, &xn, &yn, &zn) ;
    xnv = nint(xn) ;
    ynv = nint(yn) ;
    znv = nint(zn) ;
    if (xnv == xv && ynv == yv && znv == zv)
    {
      return(1) ;
    }
  }
  return(0) ;
}
#endif

#if 0
static double
mrisComputeDefectCurvatureEnergy
(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans)
{
  double sse = 0.0, nx, ny, nz, dx, dy, dz, vtotal, dot ;
  int    i, vno, n, nvertices = 0 ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
    {
      vno = vertex_trans[defect->vertices[i]] ;
    }
    else
    {
      vno = vertex_trans[defect->border[i-defect->nvertices]] ;
    }
    if (vno < 0)
    {
      continue ;
    }
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    nvertices++ ;
    mrisComputeDefectVertexNormal(mris, vno, &nx, &ny, &nz) ;
    vtotal = 0.0f ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dx = vn->origx - v->origx ;
      dy = vn->origy - v->origy ;
      dz = vn->origz - v->origz ;
      dot = dx*nx + dy*ny + dz*nz ;
      vtotal += (dot*dot) ;
      if (!isfinite(vtotal))
      {
        DiagBreak() ;
      }
    }
    if (v->vnum == 0)
    {
      continue ;
    }
    sse += (vtotal / (float)v->vnum) ;
  }
  return(sse) ;
}
#endif
#if 0
static int
mrisFindSecondNeighborhood
(MRI_SURFACE *mris, int vno, int *nbrs, int *num_nbrs)
{
  int     n, n2 ;
  VERTEX  *v, *vn, *vn2 ;

  *num_nbrs = 0 ;
  v = &mris->vertices[vno] ;
  v->marked = 1 ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    vn->marked = 1 ;
    nbrs[*num_nbrs] = v->v[n] ;
    *num_nbrs += 1 ;
  }

  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    for (n2 = 0 ; n2 < vn->vnum ; n2++)
    {
      vn2 = &mris->vertices[vn->v[n2]] ;
      if (vn2->marked)
      {
        continue ;
      }
      vn2->marked = 1 ;
      nbrs[*num_nbrs] = vn->v[n2] ;
      *num_nbrs += 1 ;
    }
  }

  v->marked = 0 ;
  for (n = 0 ; n < *num_nbrs ; n++)
  {
    mris->vertices[nbrs[n]].marked = 0 ;
  }
  return(NO_ERROR) ;
}
#endif

#if 0
static double
mrisComputeDefectQuadraticCurvatureEnergy
(MRI_SURFACE *mris, DEFECT *defect, int *vertex_trans)
{
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n, i ;
  VERTEX   *v, *vn ;
  float    ui, vi, rsq, a, b ;
  double   sse = 0.0 ;
  int      nbrs[MAX_NBRS], num_nbrs ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

  mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(2, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
    {
      vno = vertex_trans[defect->vertices[i]] ;
    }
    else
    {
      vno = vertex_trans[defect->border[i-defect->nvertices]] ;
    }
    if (vno < 0)
    {
      continue ;
    }
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->vnum <= 0)
    {
      continue ;
    }
    mrisFindSecondNeighborhood(mris, vno, nbrs, &num_nbrs) ;
    if (num_nbrs < 3)
    {
      continue ;
    }
    MRIScomputeSecondFundamentalFormAtVertex(mris, vno, nbrs, num_nbrs) ;
    v_Y = VectorAlloc(num_nbrs, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(num_nbrs, 2, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < num_nbrs ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[nbrs[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ;
      vi = V3_DOT(v_e2, v_nbr) ;
      rsq = ui*ui + vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsq ;
      *MATRIX_RELT(m_R, n+1, 2) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    if (!m_R_inv)
    {
      MatrixFree(&m_R) ;
      VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
    if (!isfinite(b))
    {
      DiagBreak() ;
    }
    sse += b*b ;
    if (vno == Gdiag_no && Gdiag > 0)
    {
      printf("v %d: curvature sse %2.2f\n", vno, b*b) ;
    }
    MatrixFree(&m_R) ;
    VectorFree(&v_Y) ;
    MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_nbr) ;
  VectorFree(&v_A) ;

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  return(sse) ;
}
#endif

#if 0
static int writeOverlap(MRIS *mris,int nclusters)
{
  int n,m,found;
  float diameter,dist,curv;
  VERTEX *vn,*vm;

  diameter=0.0f;
  for (n = 0 ; n < mris->nvertices ; n++)
  {
    vn=&mris->vertices[n];
    for (m = 0 ; m < mris->nvertices ; m++)
    {
      vm=&mris->vertices[m];
      diameter=MAX(diameter,SQR(vn->x-vm->x)+
                   SQR(vn->y-vm->y)+SQR(vn->z-vm->z));
    }
  }

  fprintf(stderr,"the diameter of the defect is %2.3f ",sqrt(diameter));

  if (nclusters<2)
  {
    nclusters=2;
  }
  dist=MAX(1.0,diameter/(float)SQR(nclusters)); /* work with the
                                                   square of the distance */

  fprintf(stderr,"(clustering distance : %2.3f)\n",sqrt(dist));

  curv=0.0f;
  found=1;
  while (found)
  {
    found=0;
    for (n = 0 ; n < mris->nvertices ; n++)
    {
      vn=&mris->vertices[n];
      if (vn->marked)
      {
        continue;
      }
      curv += 1.0f;
      vn->curv = curv;
      vn->marked = 1;
      found=1;
      for (m = 0 ; m < mris->nvertices ; m++)
      {
        vm=&mris->vertices[m];
        if (vm->marked)
        {
          continue;
        }
        if (SQR(vn->x-vm->x)+SQR(vn->y-vm->y)+SQR(vn->z-vm->z)>dist)
        {
          continue;
        }
        vm->curv = curv;
        vm->marked=1;
      }
    }
  }
  /* unmark vertices */
  for (n = 0 ; n < mris->nvertices ; n++)
  {
    mris->vertices[n].marked=0;
  }

  return NO_ERROR;
}
#endif

#if 0
static int generateMovie(MRIS *mris_src)
{
  int n,counter;
  static int nbr=0;
  float t,x,y,z;
  VERTEX *vd,*vs;
  char fname[100];

  MRIS *mris;
  mris=MRISclone(mris_src);
  nbr++; /* static */

  for (counter = 0 , t=0.0f ; t<=1.0f ; t += 0.01f)
  {
    /* compute new coords */
    for (n = 0 ; n < mris->nvertices; n++)
    {
      vd=&mris->vertices[n];
      vs=&mris_src->vertices[n];
      vd->x=t*vs->cx+(1.0f-t)*vs->origx;
      vd->y=t*vs->cy+(1.0f-t)*vs->origy;
      vd->z=t*vs->cz+(1.0f-t)*vs->origz;
    }
    /* center coords to (0,0,0) */
    x=y=z=0.0f;
    for (n = 0 ; n < mris->nvertices ; n++)
    {
      vd=&mris->vertices[n];
      x += vd->x;
      y += vd->y;
      z += vd->z;
    }
    x /= (float)mris->nvertices;
    y /= (float)mris->nvertices;
    z /= (float)mris->nvertices;
    for (n = 0 ; n < mris->nvertices ; n++)
    {
      vd=&mris->vertices[n];
      vd->x -= x;
      vd->y -= y;
      vd->z -= z;
    }
    /* write out image */
    sprintf(fname,"./mov%d_%03d",nbr,counter);
    MRISwrite(mris,fname);
    counter++;
  }

  MRISfree(&mris);
  return NO_ERROR;
}
#endif

#if 0
static int
mrisAverageWeightedGradients(MRI_SURFACE *mris, int num_avgs)
{
  int    vno, vlist[MAX_NBRS], nbrs, n, n2 ;
  float  nx, ny, nz, dx, dy, dz, sigma, wts[MAX_NBRS], total_wt, wt ;
  VERTEX *v, *vn, *vn2 ;
  MRI_SP *mrisp, *mrisp_blur ;

  if (num_avgs <= 0)
  {
    return(NO_ERROR) ;
  }

  sigma = sqrt((double)num_avgs) * M_PI / 2.0 ;
  if (Gdiag_no >= 0)
  {
    v = &mris->vertices[Gdiag_no] ;
    fprintf(stdout, "before averaging dot = %2.2f ",
            v->dx*v->nx+v->dy*v->ny+v->dz*v->nz) ;
  }
  if (0 && mris->status == MRIS_PARAMETERIZED_SPHERE)  /* use convolution */
  {
    sigma = sqrt((float)num_avgs) / 4.0 ;
    mrisp = MRISgradientToParameterization(mris, NULL, 1.0) ;
    mrisp_blur = MRISPblur(mrisp, NULL, sigma, -1) ;
    MRISgradientFromParameterization(mrisp_blur, mris) ;
    MRISPfree(&mrisp) ;
    MRISPfree(&mrisp_blur) ;
  }
  else
  {
    MRISclearMarks(mris) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }

      if (v->ripflag)
      {
        continue ;
      }

      nx = v->nx ;
      ny = v->ny ;
      nz = v->nz ;
      dx = v->dx ;
      dy = v->dy ;
      dz = v->dz ;

      /* find all 1-neighbors */
      nbrs = 1 ;
      vlist[0] = vno ;
      wts[0] = 1.0 ;
      for (n = 0  ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked)
        {
          continue ;
        }
        vn->marked = 1 ;
        wt = vn->nx*nx + vn->ny*ny + vn->nz*nz ;
        if (wt < 0)
        {
          wt = 0 ;
        }
        vlist[nbrs] = v->v[n] ;
        wts[nbrs] = exp(-1.0/(2.0f*sigma*sigma))*wt ;
        nbrs++ ;
      }
      /* find all 2-neighbors */
      for (n = v->vnum ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked || v->v[n] == vno)
        {
          continue ;
        }
        vn->marked = 2;
        wt = vn->nx*nx + vn->ny*ny + vn->nz*nz ;
        if (wt < 0)
        {
          wt = 0 ;
        }
        vlist[nbrs] = v->v[n] ;
        wts[nbrs] = exp(-4.0/(2.0f*sigma*sigma))*wt ;
        nbrs++ ;
      }
      /* find all 3-neighbors */
      for (n = v->vnum ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked != 2)  /* a two-neighbor */
        {
          continue ;
        }
        for (n2 = 0 ; n2 < vn->vnum ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (vn2->marked || vn->v[n2] == vno)
          {
            continue ;
          }
          vn2->marked = 3 ;
          wt = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (wt < 0)
          {
            wt = 0 ;
          }
          vlist[nbrs] = vn->v[n2] ;
          wts[nbrs] = exp(-9.0/(2.0f*sigma*sigma))*wt ;
          nbrs++ ;
        }
        for (n2 = vn->vnum ; n2 < vn->vtotal ; n2++)
        {
          vn2 = &mris->vertices[vn->v[n2]] ;
          if (vn2->marked || vn->v[n2] == vno)
          {
            continue ;
          }
          vn2->marked = 4 ;
          wt = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (wt < 0)
          {
            wt = 0 ;
          }
          vlist[nbrs] = vn->v[n2] ;
          wts[nbrs] = exp(-16.0/(2.0f*sigma*sigma))*wt ;
          nbrs++ ;
        }
      }

      v->tdx = v->tdy = v->tdz = 0.0 ;
      for (total_wt = 0.0, n = 0 ; n < nbrs ; n++)
      {
        wt = wts[n] ;
        total_wt += wt ;
        vn = &mris->vertices[vlist[n]] ;
        if (vlist[n] == Gdiag_no || vlist[n] == Gx)
        {
          DiagBreak() ;
        }
        vn->marked = 0 ;
        v->tdx += wt * vn->dx ;
        v->tdy += wt * vn->dy ;
        v->tdz += wt * vn->dz ;
      }
      v->tdx /= total_wt ;
      v->tdy /= total_wt ;
      v->tdz /= total_wt ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }
      v->dx = v->tdx ;
      v->dy = v->tdy ;
      v->dz = v->tdz ;
    }
  }
  if (Gdiag_no >= 0)
  {
    float dot ;
    v = &mris->vertices[Gdiag_no] ;
    dot = v->nx*v->dx + v->ny*v->dy + v->nz*v->dz ;
    fprintf(stdout, " after dot = %2.2f (%2.3f, %2.3f, %2.3f)\n",dot, v->dx, v->dy, v->dz) ;
    if (fabs(dot) > 50)
    {
      DiagBreak() ;
    }
  }
  return(NO_ERROR) ;
}
#endif

#if 0
static int
mrisMarkSulcalVertices(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z ;
  double    val0, xw,yw,zw ;
  double  del0, dot ;

  MRISclearFlags(mris, VERTEX_SULCAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    x = v->x ;
    y = v->y ;
    z = v->z ;

    // MRIworldToVoxel(parms->mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRISsurfaceRASToVoxelCached(mris, parms->mri_brain, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolume(parms->mri_brain, xw, yw, zw, &val0) ;
    dot = v->dx * v->nx + v->dy * v->ny + v->dz * v->nz ;

    del0 = v->val - val0 ;
    if (dot < 0 && del0 > 0)  /* too bright and moving inward */
    {
      v->flags |= VERTEX_SULCAL ;
      if (vno == Gdiag_no)
        printf("v %d: intensity %2.1f darker than "
               "target %2.1f - marked as sulcal\n",
               vno, val0, v->val) ;
    }
  }

  return(NO_ERROR) ;
}
static int
mrisUpdateSulcalGradients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int     vno, num ;
  VERTEX  *v ;
  double  dot ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    if (v->flags & VERTEX_SULCAL)
    {
      dot = v->dx * v->nx + v->dy * v->ny + v->dz * v->nz ;
      if (dot > 0)  /* now moving outwards - take out normal component */
      {
        num++ ;
        v->dx -= dot*v->nx ;
        v->dy -= dot*v->ny ;
        v->dz -= dot*v->nz ;
        if (vno == Gdiag_no)
          printf
          ("v %d: removing normal component "
           "%2.3f to prevent sulcal crossing\n",
           vno, dot) ;
      }
    }
  }

  printf("%d vertices detected in sulcal-crossing\n", num) ;
  return(NO_ERROR) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisAddVertices(MRI_SURFACE *mris, double thresh)
{
  double   dist ;
  int      vno, nadded, n,nvertices, nfaces, nedges, eno ;
  VERTEX   *v, *vn ;
  float    x, y, z ;

  /* make it squared so we don't need sqrts later */
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "dividing edges more than %2.2f mm long.\n", thresh) ;
  }
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    x = v->origx ;
    y = v->origy ;
    z = v->origz ;

    /*
      only add vertices if average neighbor vector is in
      normal direction, that is, if the region is concave or sulcal.
    */
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dist = sqrt(SQR(vn->origx-x) + SQR(vn->origy - y) + SQR(vn->origz - z));
      if (dist > thresh)
      {
        if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
        {
          nadded++ ;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "%d vertices added: # of vertices=%d, # of faces=%d.\n",
            nadded, mris->nvertices, mris->nfaces) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
  }
  return(nadded) ;
}
#endif
#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Orient the faces of the tessellation so that the
  point outward (i.e. in the same direction as the original
  surface).
  ------------------------------------------------------*/
static int
mrisOrientFaces(MRI_SURFACE *mris, DEFECT *defect, int *vtrans)
{
  int    vno, vno0, vno1, i, n, fno, oriented = 0 ;
  FACE   *f ;
  VERTEX *v, *vn ;
  float  norm[3], dot ;

  /* first orient vertices */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    vno = vtrans[defect->vertices[i]] ;
    if (vno < 0)   /* not in new tessellation */
    {
      continue ;
    }
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    for (dot = 0.0, n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dot += v->nx*vn->nx + v->ny*vn->ny + v->nz*vn->nz ;
    }
    if (dot < 0)
    {
      oriented++ ;
      if (Gdiag & DIAG_SHOW)
      {
        fprintf(stdout, "reversing normal for vertex %d\n", vno) ;
      }
      v->nx *= -1.0f ;
      v->ny *= -1.0f ;
      v->nz *= -1.0f ;
    }
  }
  /*  mrisRestoreDefectPositions(mris, defect, ORIGINAL_VERTICES) ;*/
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
    {
      vno = vtrans[defect->vertices[i]] ;
    }
    else
    {
      vno = vtrans[defect->border[i-defect->nvertices]] ;
    }
    if (vno < 0)   /* not in new tessellation */
    {
      continue ;
    }
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    for (n = 0 ; n < v->num ; n++)
    {
      fno = v->f[n] ;
      if (mris->faces[fno].ripflag)  /* only consider it once */
      {
        continue ;
      }
      if (fno == Gdiag_no)
      {
        DiagBreak() ;
      }
      mris->faces[fno].ripflag = 1 ;
      mrisNormalFace(mris, fno, v->n[n], norm) ;
      dot = norm[0]*v->nx + norm[1]*v->ny + norm[2]*v->nz ;
      if (dot < 0)   /* change order of vertices in face */
      {
        oriented++ ;
        f = &mris->faces[fno] ;
        vno0 = f->v[0] ;
        vno1 = f->v[1] ;
        f->v[0] = vno1 ;
        f->v[1] = vno0 ;
        mrisSetVertexFaceIndex(mris, vno0, fno) ;
        mrisSetVertexFaceIndex(mris, vno1, fno) ;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        {
          fprintf(stdout, "reversing face %d orientation\n", fno) ;
        }
      }
    }
  }
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    if (i < defect->nvertices)
    {
      vno = vtrans[defect->vertices[i]] ;
    }
    else
    {
      vno = vtrans[defect->border[i-defect->nvertices]] ;
    }
    if (vno < 0)   /* not in new tessellation */
    {
      continue ;
    }
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->num ; n++)
    {
      fno = v->f[n] ;
      mris->faces[fno].ripflag = 0 ;
    }
  }
  /*  mrisRestoreDefectPositions(mris, defect, TMP_VERTICES) ;*/
  return(oriented) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Set the positions of all vertices in a defect to the
  original, canonical, or tmp vertices.
  ------------------------------------------------------*/
static int
mrisRestoreDefectPositions(MRI_SURFACE *mris, DEFECT *defect, int which)
{
  int     vno, i ;
  VERTEX  *v ;

  for (i = 0 ; i < defect->nvertices ; i++)
  {
    vno = defect->vertices[i] ;
    v = &mris->vertices[vno] ;
#if 0
    if (v->ripflag)
    {
      continue ;
    }
#endif
    switch (which)
    {
    case CANONICAL_VERTICES:
      v->x = v->cx ;
      v->y = v->cy ;
      v->z = v->cz ;
      break ;
    case ORIGINAL_VERTICES:
      v->x = v->origx ;
      v->y = v->origy ;
      v->z = v->origz ;
      break ;
    case TMP2_VERTICES:
      v->x = v->tx2 ;
      v->y = v->ty2 ;
      v->z = v->tz2 ;
      break ;
    default:
    case TMP_VERTICES:
      v->x = v->tx ;
      v->y = v->ty ;
      v->z = v->tz ;
      break ;
    }
  }
  return(NO_ERROR) ;
}
#endif
#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Add the appropriate face to the tessellation.
  ------------------------------------------------------*/
static int
mrisAddDefectFaces(MRI_SURFACE *mris, double e0[3], double e1[3],
                   double origin[3], EDGE *et, int nedges)
{
  int     i, vno1, vno2, nfaces, n ;
  VERTEX  *v ;

  for (i = 0 ; i < nedges ; i++)
  {
    if (et[i].used == 0)
    {
      continue ;
    }
    vno1 = et[i].vno1 ;
    vno2 = et[i].vno2 ;
#if (!SPHERE_INTERSECTION)
    mrisComputeCanonicalEdgeBasis(mris, et+i, et+i, origin, e0, e1) ;
#endif
    if (vno1 == 108332 && vno2 == 109240)
    {
      DiagBreak() ;
    }

    /* for every vertex which is a neighbor of
       both of these, add 1 triangle */
    v = &mris->vertices[vno1] ;
    for (nfaces = n = 0 ; n < v->vnum ; n++)
    {
      if (v->v[n] == vno2)
      {
        continue ;
      }
      if (vertexNeighbor(mris, vno2, v->v[n]) &&
          !isFace(mris,vno1, vno2, v->v[n]) &&
#if SPHERE_INTERSECTION
          !containsAnotherVertexOnSphere(mris,vno1,vno2,v->v[n],0))
      {
#else
          !containsAnotherVertex(mris,vno1,vno2,v->v[n],e0,e1,origin))
      {
#endif
        if (nfaces++ > 1)
        {
          DiagBreak() ;
        }
        if (mris->nfaces == Gdiag_no)
        {
          DiagBreak() ;
        }
#if 0
        if (((vno1 != 110718) && (vno1 != 110732) && (vno1 != 110748)) ||
            ((vno2 != 110718) && (vno2 != 110732) && (vno2 != 110748)) ||
            ((v->v[n] != 110718) && (v->v[n] != 110732) &&
             (v->v[n]!=110748)))
#endif
          mrisAddFace(mris, vno1, vno2, v->v[n]) ;
      }
    }
  }
  return(NO_ERROR) ;
}
#endif


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale the surface so that it's max dimension has length maxr
  ------------------------------------------------------*/
static int
mrisScaleMaxDimension(MRI_SURFACE *mris, float maxr)
{
  maxr /= 2.0 ;
  mrisComputeSurfaceDimensions(mris) ;
  if (mris->xhi >= maxr)
  {
    MRISscaleBrain(mris, mris, maxr / mris->xhi) ;
  }
  if (mris->yhi >= maxr)
  {
    MRISscaleBrain(mris, mris, maxr / mris->yhi) ;
  }
  if (mris->yhi >= maxr)
  {
    MRISscaleBrain(mris, mris, maxr / mris->yhi) ;
  }
  if (fabs(mris->xlo) >= maxr)
  {
    MRISscaleBrain(mris, mris, maxr / fabs(mris->xlo)) ;
  }
  if (fabs(mris->ylo) >= maxr)
  {
    MRISscaleBrain(mris, mris, maxr / fabs(mris->ylo)) ;
  }
  if (fabs(mris->zlo) >= maxr)
  {
    MRISscaleBrain(mris, mris, maxr / fabs(mris->zlo)) ;
  }
  return(NO_ERROR) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Return count of # of vertices shared by 2 triangles
  ------------------------------------------------------*/
static int
triangleNeighbors(MRI_SURFACE *mris, int fno1, int fno2)
{
  int  n1, n2, num ;
  FACE *f1, *f2 ;

  f1 = &mris->faces[fno1] ;
  f2 = &mris->faces[fno2] ;
  for (num = n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
  {
    for (n2 = 0 ; n2 < VERTICES_PER_FACE ; n2++)
      if (f1->v[n1] == f2->v[n2])
      {
        num++ ;
      }
  }
  return(num) ;
}
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  See if the edge e intersects any edges in the defect or
  it's border. Sorry this code is such a hatchet job. I'm sure
  there are far more elegant ways of doing intersection
  (e.g. sorting).
  ------------------------------------------------------*/
#if 0
static int
colinearDefectEdges(MRI_SURFACE *mris, DEFECT *defect, EDGE *e,
                    int *vertex_trans)
{
  VERTEX *v ;
  int    i, n, vno ;
  EDGE   edge2 ;

  if ((e->vno1 == 109259 && e->vno2 == 108332) ||
      (e->vno2 == 109259 && e->vno1 == 108332))
  {
    DiagBreak() ;
  }

  for (i = 0 ; i < defect->nvertices+defect->nchull ; i++)
  {
    if (i < defect->nvertices)
    {
      vno = vertex_trans[defect->vertices[i]] ;
    }
    else
    {
      vno = vertex_trans[defect->chull[i-defect->nvertices]] ;
    }

    if (vno < 0)
    {
      continue ;
    }

    edge2.vno1 = vno ;
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }

      if ((vno == 53153 && v->v[n] == 53158) ||
          (v->v[n] == 53153 && vno == 53158))
      {
        DiagBreak() ;
      }
      if ((vno == 52024 && v->v[n] == 53158) ||
          (v->v[n] == 52024 && vno == 53158))
      {
        DiagBreak() ;
      }

      edge2.vno2 = v->v[n] ;

      if ((v->v[n] == e->vno1 || v->v[n] == e->vno2) ||
          (vno == e->vno1 || vno == e->vno2))  /* check for colinearity */
      {
        int vno0, vno1, vno2, ncoords ;
        VERTEX *v0, *v1, *v2 ;
        float  dx01, dy01, dz01, dx02, dy02, dz02, len01, len02 ;

        /*
          set vno0 and vno1 so that they are the existing edge, with
          vno0 being the shared vertex, and vno2 is the vertex for the
          putative edge
        */
        if (vno == e->vno1 || vno == e->vno2)  /* vno is shared vertex */
        {
          vno0 = vno ;
          vno1 = v->v[n] ;
        }
        else   /* v->v[n] is shared vertex */
        {
          vno0 = v->v[n] ;
          vno1 = vno ;
        }
        if (e->vno1 == vno0)
        {
          vno2 = e->vno2 ;
        }
        else
        {
          vno2 = e->vno1 ;
        }
        v0 = &mris->vertices[vno0] ;
        v1 = &mris->vertices[vno1] ;
        v2 = &mris->vertices[vno2] ;
        dx01 = v1->origx - v0->origx ;
        dy01 = v1->origy - v0->origy ;
        dz01 = v1->origz - v0->origz ;
        len01 = sqrt(dx01*dx01+dy01*dy01+dz01*dz01) ;
        if (FZERO(len01))
        {
          len01 = 1 ;
        }
        dx01 /= len01 ;
        dy01 /= len01 ;
        dz01 /= len01 ;
        dx02 = v2->origx - v0->origx ;
        dy02 = v2->origy - v0->origy ;
        dz02 = v2->origz - v0->origz ;
        len02 = sqrt(dx02*dx02+dy02*dy02+dz02*dz02) ;
        if (FZERO(len02))
        {
          len02 = 1 ;
        }
        dx02 /= len02 ;
        dy02 /= len02 ;
        dz02 /= len02 ;


        /* see if vno1 and vno2 are colinear. If not, no problemo */
        ncoords = FZERO(dx01-dx02)+FZERO(dy01-dy02)+FZERO(dz01-dz02);
        if (ncoords == 3)
        {
          if (e->vno1 == 16968 || e->vno2 == 16968)
          {
            DiagBreak() ;
          }
          return(1) ;   /* colinear */
        }
      }
    }
  }

  return(0) ;
}
#endif


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisComputeCanonicalBasis(MRI_SURFACE *mris, int fno, double origin[3],
                          double e0[3], double e1[3])
{
  FACE    *f ;
  double   len, normal[3] ;
  float    fx, fy, fz ;

  f = &mris->faces[fno] ;
  mrisCalculateCanonicalFaceCentroid(mris, fno, &fx, &fy, &fz) ;
  origin[0] = (double)fx ;
  origin[1] = (double)fy ;
  origin[2] = (double)fz ;
  normal[0] = origin[0] ;
  normal[1] = origin[1] ;
  normal[2] = origin[2] ;
  len = 1.0f/VLEN(normal) ;
  SCALAR_MUL(normal, len, normal) ;

  /* pick any non-parallel vector and cross it with normal */
  e1[0] = normal[1] ;
  e1[1] = normal[2] ;
  e1[2] = normal[0] ;
  CROSS(e0, normal, e1) ;
  if ((VZERO(e0)))  /* happened to pick parallel vector */
  {
    e1[0] = normal[1] ;
    e1[1] = -normal[2] ;
    e1[2] = normal[0] ;
    CROSS(e0, normal, e1) ;
  }
  CROSS(e1, e0, normal) ;
  len = 1.0f/VLEN(e0) ;
  SCALAR_MUL(e0, len, e0) ;
  len = 1.0f/VLEN(e1) ;
  SCALAR_MUL(e1, len, e1) ;
  len = DOT(e0, e1) ;
  if ((VZERO(e0)) || (VZERO(e1)))
  {
    fprintf(stdout, "face %d, canonical basis degenerate!\n", fno) ;
  }
  if (fabs(len) > 0.001)
  {
    fprintf(stdout, "face %d, canonical basis not orthogonal!\n", fno) ;
  }
  return(NO_ERROR) ;
}
#endif

#if 0
static int
mrisSoapBubbleIntersectingDefects(MRI_SURFACE *mris)
{
  int      vno, num, vno2, n ;
  VERTEX   *v, *v2 ;

  /* coming in the v->marked field lists the defect #. If not part
     of a defect v->marked == -1.
  */
  MRIScopyMarkedToMarked2(mris) ;
  n = 0 ;
  while ((num = mrisMarkIntersections(mris)) > 0)
  {
    // mark all of each defect that has an intersection
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->marked)   // it is intersecting - mark the whole defect
      {
        if (v->marked2 >= 0)  // part of a defect
        {
          for (vno2 = 0 ; vno2 < mris->nvertices ; vno2++)
          {
            v2 = &mris->vertices[vno2] ;
            if (v2->marked2 == v->marked2)  // part of same defect
            {
              v2->marked = 1 ;
            }
          }
        }
      }
    }

    printf("%03d: %d intersecting\n", n, num) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->marked = !v->marked ;  // soap bubble will fix the marked ones
    }
    MRISsoapBubbleVertexPositions(mris, 5) ;
    if (n++ > 100)
    {
      break ;
    }
  }
  return(NO_ERROR) ;
}
#endif

#if 0
static MRI_SP *
MRISPiterative_blur(MRI_SURFACE *mris,
                    MRI_SP *mrisp_src,
                    MRI_SP *mrisp_dst,
                    float sigma, int frame)
{
  int niter ;
  float *curvs = (float *)calloc(mris->nvertices, sizeof(float)) ;

  if (!mrisp_dst)
  {
    mrisp_dst = MRISPclone(mrisp_src) ;
  }
  mrisp_dst->sigma = sigma ;

  MRISextractCurvatureVector(mris, curvs) ;

  MRISfromParameterization(mrisp_src, mris, frame);

  niter = nint(sigma*sigma*M_PI/2) ;
  MRISaverageCurvatures(mris, niter) ;
  MRIStoParameterization(mris, mrisp_dst, 1, frame) ;

  MRISimportCurvatureVector(mris, curvs) ;
  free(curvs) ;
  return(mrisp_dst) ;
}
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if 0
static int
mrisReadFieldsign(MRI_SURFACE *mris, const char *mris_fname)
{
  int k,i,vnum;
  float f;
  FILE *fp;

  printf("surfer: read_fieldsign(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  fread(&vnum,1,sizeof(int),fp);
  printf("surfer: vertex_index = %d, vnum = %d\n",vertex_index,vnum);
  if (vnum!=vertex_index)
  {
    printf("surfer: Warning: incompatible vertex number in file %s\n",fname);
  }
  for (k=0; k<vnum; k++)
  {
    fread(&f,1,sizeof(float),fp);
    vertex[k].fieldsign = f;
  }
  fclose(fp);
  fieldsignflag = TRUE;
  PR
  return(NO_ERROR) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
normalize_binary_curvature(MRI_SURFACE *mris)
{
  int k;
  float curv,min,max;
  float sum,avg,sq,sum_sq,sd,n;
  FILE *fp;

  if (!curvloaded)
  {
    printf("surfer: ### curv not loaded!\n");
    PR return;
  }

  sum = 0;
  for (k=0; k<vertex_index; k++)
  {
    sum += vertex[k].curv;
  }
  avg = sum/vertex_index;

  n = (float)vertex_index;
  sum = sum_sq = 0.0;
  for (k=0; k<vertex_index; k++)
  {
    vertex[k].curv -= avg;
    curv = vertex[k].curv;
    sum += curv;
    sum_sq += curv*curv;
  }
  sd = sqrt((n*sum_sq - sum*sum)/(n*(n-1.0)));

  for (k=0; k<vertex_index; k++)
  {
    curv = (vertex[k].curv)/sd;
    if (k==0)
    {
      min=max=curv;
    }
    if (curv<min)
    {
      min=curv;
    }
    if (curv>max)
    {
      max=curv;
    }
    if (curv<CURVIM_NORM_MIN)
    {
      curv = CURVIM_NORM_MIN;
    }
    if (curv>CURVIM_NORM_MAX)
    {
      curv = CURVIM_NORM_MAX;
    }
    vertex[k].curv = curv;
  }
  curvmin = CURVIM_NORM_MIN;
  curvmax = CURVIM_NORM_MAX;
  printf("surfer: curvature normalized: avg=%f sd=%f\n",avg,sd);
  printf("surfer: min=%f max=%f trunc to (%f,%f)\n",min,max,curvmin,curvmax);
  PR
}

#endif


#if 0

static double mrisComputeVertexNormalSpacingStats(MRI_SURFACE *mris,
    static int mrisTessellateFace(MRI_SURFACE *mris, int fno) ;
    static int VertexReplaceNeighbor(VERTEX *v, int vno_old, int vno_new) ;
    double *psigma);

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisAddVertices(MRI_SURFACE *mris, double nsigma)
{
  double   mean, sigma, dist, thresh, dot ;
  int      vno, nadded, n,nvertices, nfaces, nedges, eno, added ;
  VERTEX   *v, *vn ;
  float    x, y, z, sx, sy, sz, nx, ny, nz ;

#if 1
  mean = MRIScomputeVertexSpacingStats(mris, &sigma, NULL, NULL, NULL, NULL, CURRENT_VERTICES) ;
#else
  mean = mrisComputeVertexNormalSpacingStats(mris, &sigma) ;
  thresh *= thresh ;   /* make it squared so we don't need sqrts later */
#endif
  thresh = mean + sigma * nsigma ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "dividing edges more than %2.2f mm long.\n", thresh) ;
  }
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }

    /* calculate average normal displacement to neighbors */
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    sx = sy = sz = 0.0 ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (!vn->ripflag)
      {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
      }
    }
    dot = sx*nx+sy*ny+sz*nz;   /* projection onto normal */

    /*
      only add vertices if average neighbor vector is in
      normal direction, that is, if the region is concave or sulcal.
    */
    for (added = n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dist = sqrt(SQR(vn->x - x) + SQR(vn->y - y) + SQR(vn->z - z)) ;
      if (dist > thresh)
      {
        if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
        {
          nadded++ ;
        }
      }
    }

    /* check for sulcal vertices that have asymptoted */
    if (!added && v->marked && dot >= 0.0f)
    {
      dot = v->odx * nx + v->ody * ny + v->odz * nz ;
      if (dot > 0.0f)
      {
        for (n = 0 ; n < v->vnum ; n++)
        {
          vn = &mris->vertices[v->v[n]] ;
          dist = sqrt(SQR(vn->x - x) + SQR(vn->y - y) + SQR(vn->z - z)) ;
          if (dist > mean)
          {
            if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
            {
              added++ ;
            }
          }
        }
      }
    }
    nadded += added ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "%d vertices added: # of vertices=%d, # of faces=%d.\n",
            nadded, mris->nvertices, mris->nfaces) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
  }
  return(mean) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  NOT FINISHED!!!!!!
  ------------------------------------------------------*/
static  int
mrisTessellateFace(MRI_SURFACE *mris, int fno)
{
  VERTEX  *vc, *v1, *v2, *v, *vnew[3] ;
  int     n, vno1, vno2, vnew_no, nv, fnew_no, vc_no, vlist[15] ;
  float   x, y, z, dx, dy, dz, ox, oy, oz, val ;
  FACE    *f, *fnew ;

  if (mris->nvertices + 4 >= mris->max_vertices)
  {
    return(ERROR_NO_MEMORY) ;
  }
  if (mris->nfaces + 5 >= mris->max_faces)
  {
    return(ERROR_NO_MEMORY) ;
  }

  f = &mris->faces[fno] ;

  /* find centroid of current face and put a vertex there */
  ox = oy = oz = x = y = z = dx = dy = dz = val = 0.0 ;
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    v = &mris->vertices[f->v[n]] ;
    x += v->x ;
    y += v->y ;
    z += v->z ;
    dx += v->odx ;
    dy += v->ody ;
    dz += v->odz ;
    ox += v->origx ;
    oy += v->origy ;
    oz += v->origz ;
    val += v->val ;
  }

  vc_no = mris->nvertices++ ;
  vc = &mris->vertices[vc_no] ;
  vc->val = val / (float)VERTICES_PER_FACE ;
  vc->x = x / (float)VERTICES_PER_FACE ;
  vc->y = y / (float)VERTICES_PER_FACE ;
  vc->z = z / (float)VERTICES_PER_FACE ;
  vc->odx = dx / (float)VERTICES_PER_FACE ;
  vc->ody = dy / (float)VERTICES_PER_FACE ;
  vc->odz = dz / (float)VERTICES_PER_FACE ;
  vc->origx = ox / (float)VERTICES_PER_FACE ; CHANGES_ORIG
  vc->origy = oy / (float)VERTICES_PER_FACE ;
  vc->origz = oz / (float)VERTICES_PER_FACE ;
  vc->vnum = 0 ;
  vc->num = 6 ;

  /* now allocate 3 vertices which bisect each edge of the face */
  vlist[0] = f->v[0] ;
  nv = 1 ;
  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    vno1 = f->v[n] ;
    vno2 = n < VERTICES_PER_FACE-1 ? f->v[n+1] : f->v[0] ;
    v1 = &mris->vertices[vno1] ;
    v2 = &mris->vertices[vno2] ;
    vnew_no = mris->nvertices++ ;
    v = vnew[n] = &mris->vertices[vnew_no] ;
    v->val = (v1->val+v2->val) / 2.0f ;

    v->x = (v1->x+v2->x) / 2.0f ;
    v->y = (v1->y+v2->y) / 2.0f ;
    v->z = (v1->z+v2->z) / 2.0f ;

    v->odx = (v1->odx+v2->odx) / 2.0f ;
    v->ody = (v1->ody+v2->ody) / 2.0f ;
    v->odz = (v1->odz+v2->odz) / 2.0f ;

    v->origx = (v1->origx+v2->origx) / 2.0f ; CHANGES_ORIG
    v->origy = (v1->origy+v2->origy) / 2.0f ;
    v->origz = (v1->origz+v2->origz) / 2.0f ;
    v->num = 0 ;
    VertexReplaceNeighbor(v1, vno2, vnew_no) ;
    VertexReplaceNeighbor(v2, vno1, vnew_no) ;
    vlist[nv++] = vnew_no ;
    vlist[nv++] = vno2 ;

    /* now build the new vertex's neighbor list */
    v->vnum = 3 ;
    v->v[0] = vno1 ;
    v->v[1] = vno2 ;
    v->v[2] = vc_no ;
    vc->v[vc->vnum++] = vno1 ;
    vc->v[vc->vnum++] = vnew_no ;
    vc->vtotal = vc->vnum ;
  }

  /*
    at this point all the vertices and edges are in place. Now
    put in new faces, reusing the one we are supertessellating.
  */
  for (n = 0 ; n < nv-1 ; n++)
  {
    if (!n)
    {
      fnew_no = fno ;
    }
    else
    {
      fnew_no = mris->nfaces++ ;
    }
    fnew = &mris->faces[fnew_no] ;
    fnew->v[0] = vlist[n] ;
    fnew->v[1] = vlist[n+1] ;
    fnew->v[2] = vc_no ;
  }

  return(NO_ERROR) ;
}
static int
VertexReplaceNeighbor(VERTEX *v, int vno_old, int vno_new)
{
  int n ;

  for (n = 0 ; n < v->vnum ; n++)
  {
    if (v->v[n] == vno_old)
    {
      v->v[n] = vno_new ;
      break ;
    }
  }
  return(NO_ERROR) ;
}
#endif


#if 0
// defined in tritri.h
#define CROSS(dest, v1, v2)                \
  dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
  dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
  dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])
#endif
#if 0
static float mrisComputeFaceStretch(MRI_SURFACE *mris, int fno) ;
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  mris is the original (hi-res) surface (from smoothwm).
  v is a vertex on the icosahedron.
  ------------------------------------------------------*/
static int
mrisChooseFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int    i, /*fno,*/ flist[10000], nfaces, min_fno = 0, l ;
  MHBT   *bucket, bucket_bak  ;
  MHB    *bin ;
  FACE   *f ;
  VERTEX *v1, *v2, *v3 ;
  float  /*ut, vt, stretch,*/ n[3], leg[3], p[3], d[3], dot, ldot, leg2[3] ;

  bucket = MHTgetBucket(mht, v->x, v->y, v->z) ;
  bucket_bak = *bucket ;
  bin = bucket->bins ;
  for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)  /* check each face */
  {
    n[0] = n[1] = n[2] = 0.0f ;
    f = &mris->faces[bin->fno] ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      n[0] += v1->nx ;
      n[1] += v1->ny ;
      n[2] += v1->nz ;
    }
    n[0] /= VERTICES_PER_FACE ;
    n[1] /= VERTICES_PER_FACE ;
    n[2] /= VERTICES_PER_FACE ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      switch (l)
      {
      default:
      case 0:
        v2 = &mris->vertices[f->v[1]] ;
        v3 = &mris->vertices[f->v[2]] ;
        break ;
      case 1:
        v2 = &mris->vertices[f->v[2]] ;
        v3 = &mris->vertices[f->v[0]] ;
        break ;
      case 2:
        v2 = &mris->vertices[f->v[0]] ;
        v3 = &mris->vertices[f->v[1]] ;
        break ;
      }

      VERTEX_DIF(leg, v1, v2) ;   /* leg of triangle */
      VERTEX_DIF(leg2, v1, v3) ;  /* leg of triangle */

      /* express p as point in triangle plane */
      VERTEX_DIF(p, v1, v) ;     /* vector from vertex to point in question */
      dot = DOT(p,n) ;
      p[0] -= dot*n[0] ;
      p[1] -= dot*n[1] ;
      p[2] -= dot*n[2] ;
#if 0
      p[0] = ut*leg[0] + vt*leg2[0] ;
      p[1] = ut*leg[1] + vt*leg2[1] ;
      p[2] = ut*leg[2] + vt*leg2[2] ;
      ut = DOT(p, leg) ;
      vt = DOT(p, leg2) ;
#endif

      CROSS(d, leg, n) ;
      dot = DOT(d, p) ;
      ldot = DOT(d, leg2) ;

      /* on different side of leg from 3rd vertex */
      if (!FZERO(ldot) && !FZERO(dot) && dot*ldot < 0)
      {
        break ;
      }
    }
    if (l >= VERTICES_PER_FACE)
    {
      if (nfaces == 10000)
      {
        ErroExit(ERROR_BADPARM, "Too many faces");
      }
      flist[nfaces++] = bin->fno ;
    }
  }
  if (!nfaces)  /* something went wrong, but Anders will fix it */
  {
    float dist, min_dist ;

    fprintf(stderr, "no faces found on sphere!\n") ;
    bin = bucket->bins ;
    min_dist = 1000000.0f ;
    min_fno = 0 ;
    for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      f = &mris->faces[bin->fno] ;
      v1 = &mris->vertices[f->v[0]] ;
      v2 = &mris->vertices[f->v[1]] ;
      v3 = &mris->vertices[f->v[2]] ;
#define VDIST(v1, v2) (sqrt(SQR(v1->x - v2->x) + SQR(v1->y - v2->y) + SQR(v1->z - v2->z)))
      dist = VDIST(v1,v) + VDIST(v2, v) + VDIST(v3,v) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = bin->fno ;
      }
    }
  }
  else   /* pix the face that is closest to the soap bubble location */
  {
    float min_dist, x, y, z, dist, curv ;

    if (nfaces > 1)
    {
      DiagBreak() ;
      curv = 1.0f ;
    }
    else
    {
      curv = 0.0f ;
    }
    for ( i = 0 ; i < nfaces ; i++)/* check each face */
    {
      f = &mris->faces[flist[i]] ;
      for (l = 0 ; l < VERTICES_PER_FACE ; l++)
      {
        mris->vertices[f->v[l]].curv = curv ;
      }
    }

    min_dist = 100000.0f ;
    for (i = 0 ; i < nfaces ; i++)
    {
      mrisCalculateCanonicalFaceCentroid(mris, flist[i], &x, &y, &z) ;
      dist = sqrt(SQR(v->x-x)+SQR(v->y-y)+SQR(v->z-z)) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = flist[i] ;
      }
    }
  }
  return(min_fno) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale a surface anisotropically.
  ------------------------------------------------------*/
static int
mrisFindUnambiguousFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int    i, /*fno,*/ flist[10000], nfaces, min_fno = 0, l ;
  MHBT   *bucket ;
  MHB    *bin ;
  FACE   *f ;
  VERTEX *v1, *v2, *v3 ;
  float  /*ut, vt, */stretch, n[3], leg[3], p[3], d[3], dot, ldot, leg2[3] ;

  bucket = MHTgetBucket(mht, v->x, v->y, v->z) ;
  bin = bucket->bins ;
  for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)  /* check each face */
  {
    n[0] = n[1] = n[2] = 0.0f ;
    f = &mris->faces[bin->fno] ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      n[0] += v1->nx ;
      n[1] += v1->ny ;
      n[2] += v1->nz ;
    }
    n[0] /= VERTICES_PER_FACE ;
    n[1] /= VERTICES_PER_FACE ;
    n[2] /= VERTICES_PER_FACE ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      switch (l)
      {
      default:
      case 0:
        v2 = &mris->vertices[f->v[1]] ;
        v3 = &mris->vertices[f->v[2]] ;
        break ;
      case 1:
        v2 = &mris->vertices[f->v[2]] ;
        v3 = &mris->vertices[f->v[0]] ;
        break ;
      case 2:
        v2 = &mris->vertices[f->v[0]] ;
        v3 = &mris->vertices[f->v[1]] ;
        break ;
      }

      VERTEX_DIF(leg, v1, v2) ;   /* leg of triangle */
      VERTEX_DIF(leg2, v1, v3) ;  /* leg of triangle */

      /* express p as point in triangle plane */
      VERTEX_DIF(p, v1, v) ;     /* vector from vertex to point in question */
      dot = DOT(p,n) ;
      p[0] -= dot*n[0] ;
      p[1] -= dot*n[1] ;
      p[2] -= dot*n[2] ;
#if 0
      p[0] = ut*leg[0] + vt*leg2[0] ;
      p[1] = ut*leg[1] + vt*leg2[1] ;
      p[2] = ut*leg[2] + vt*leg2[2] ;
      ut = DOT(p, leg) ;
      vt = DOT(p, leg2) ;
#endif

      CROSS(d, leg, n) ;
      dot = DOT(d, p) ;
      ldot = DOT(d, leg2) ;

      /* on different side of leg from 3rd vertex */
      if (!FZERO(ldot) && !FZERO(dot) && dot*ldot < 0)
      {
        break ;
      }
    }
    if (l >= VERTICES_PER_FACE)
    {
      if (nfaces == 10000)
      {
        ErrorExit(ERROR_BADPARM, "Too many faces");
      }
      flist[nfaces++] = bin->fno ;
    }
  }
  if (!nfaces)
  {
    float dist, min_dist ;

    fprintf(stderr, "no faces found on sphere!\n") ;
    bin = bucket->bins ;
    min_dist = 1000000.0f ;
    min_fno = 0 ;
    for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      f = &mris->faces[bin->fno] ;
      v1 = &mris->vertices[f->v[0]] ;
      v2 = &mris->vertices[f->v[1]] ;
      v3 = &mris->vertices[f->v[2]] ;
#define VDIST(v1, v2) (sqrt(SQR(v1->x - v2->x) + SQR(v1->y - v2->y) + SQR(v1->z - v2->z)))
      dist = VDIST(v1,v) + VDIST(v2, v) + VDIST(v3,v) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = bin->fno ;
      }
    }
    printf("min_dist = %f (min_fno=%d)\n",min_dist,min_fno);

  }
  else
  {
    float min_stretch, curv ;

    if (nfaces > 1)
    {
      DiagBreak() ;
      curv = 1.0f ;
    }
    else
    {
      curv = 0.0f ;
    }
    for ( i = 0 ; i < nfaces ; i++)/* check each face */
    {
      f = &mris->faces[flist[i]] ;
      for (l = 0 ; l < VERTICES_PER_FACE ; l++)
      {
        mris->vertices[f->v[l]].curv = curv ;
      }
    }

    min_stretch = 100000.0f ;
    for (i = 0 ; i < nfaces ; i++)
    {
      stretch = mrisComputeFaceStretch(mris, flist[i]) ;
      if (stretch < min_stretch)
      {
        min_stretch = stretch ;
        min_fno = flist[i] ;
      }
    }
  }
  if (nfaces <= 1)
  {
    return(min_fno) ;
  }
  else
  {
    return(-1) ;
  }
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale a surface anisotropically.
  ------------------------------------------------------*/
static float
mrisComputeFaceStretch(MRI_SURFACE *mris, int fno)
{
  int    fv, n0, n1 ;
  float  stretch, max_stretch, dist, inflated_dist ;
  FACE   *f ;
  VERTEX *v0, *v1 ;

  f = &mris->faces[fno] ;
  max_stretch = -1.0f ;
  for (fv = 0 ; fv < VERTICES_PER_FACE ; fv++)
  {
    n0 = f->v[fv] ;
    n1 = fv < VERTICES_PER_FACE - 1 ? f->v[fv+1] : f->v[0] ;
    v0 = &mris->vertices[n0] ;
    v1 = &mris->vertices[n1] ;
    inflated_dist =
      SQR(v0->cx-v1->cx) + SQR(v0->cy-v1->cy) + SQR(v0->cz-v1->cz);
    dist =
      SQR(v0->origx-v1->origx) +
      SQR(v0->origy-v1->origy) + SQR(v0->origz-v1->origz);
    if (!FZERO(dist))
    {
      stretch = inflated_dist / dist ;
      if (stretch > max_stretch)
      {
        max_stretch = stretch ;
      }
    }
  }
  return(max_stretch) ;
}
#endif


#if 0
/*-----------------------------------------------------
  MRISfastRead() just calls MRISRead()
  Parameters:
  Returns value:
  Description
  ------------------------------------------------------*/
MRI_SURFACE *MRISfastRead(const char *fname)
{
/********* why you keep the rest ? ******************/
#if 1
  return (MRISread(fname));
#else
  MRI_SURFACE *mris;
  int nquads, nvertices, magic, version, ix, iy, iz, vno, fno, n, m;
  int imnr, imnr0, imnr1, type, vertices[4], num;
  float x, y, z, xhi, xlo, yhi, ylo, zhi, zlo;
  FILE *fp;
  VERTEX *vertex;
  FACE *face;

  mris = NULL;
  fp = NULL;
  type = MRISfileNameType(fname);
  if (type == MRIS_ASCII_TRIANGLE_FILE) {
    mris = mrisReadAsciiFile(fname);
    if (!mris) {
      return (NULL);
    }
    version = -3;
  }
  else if (type == MRIS_ICO_FILE) {
    mris = ICOread(fname);
    if (!mris) {
      return (NULL);
    }
    return (mris);
    version = -2;
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE) {
    mris = mrisReadGeoFile(fname);
    if (!mris) {
      return (NULL);
    }
    version = -4;
  }
  else /* custom binary file - find out which type using magic # */
  {
    fp = fopen(fname, "rb");
    if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "MRISread(%s): could not open file", fname));

    fread3(&magic, fp);
    if (magic == TRIANGLE_FILE_MAGIC_NUMBER) {
      fclose(fp);
      mris = mrisReadTriangleFile(fname, 0.0);
      if (!mris) {
        ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n"));
      }
      version = -3;
    }
    else if (magic == QUAD_FILE_MAGIC_NUMBER) {
      version = -1;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        fprintf(stdout, "new surface file format\n");
      }
    }
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) {
      version = -2;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        fprintf(stdout, "new surface file format\n");
      }
    }
    else {
      rewind(fp);
      version = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        printf("surfer: old surface file format\n");
      }
    }
  }
  if (version >= -2) {
    fread3(&nvertices, fp);
    fread3(&nquads, fp);

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "reading %d vertices and %d faces.\n", nvertices, 2 * nquads);

    mris = MRISalloc(nvertices, 2 * nquads);
    mris->type = MRIS_BINARY_QUADRANGLE_FILE;

    imnr0 = 1000;
    imnr1 = 0;
    for (vno = 0; vno < nvertices; vno++) {
      vertex = &mris->vertices[vno];
      if (version == -1) {
        fread2(&ix, fp);
        fread2(&iy, fp);
        fread2(&iz, fp);
        vertex->x = ix / 100.0;
        vertex->y = iy / 100.0;
        vertex->z = iz / 100.0;
      }
      else /* version == -2 */
      {
        vertex->x = freadFloat(fp);
        vertex->y = freadFloat(fp);
        vertex->z = freadFloat(fp);
      }
#if 0
      vertex->label = NO_LABEL ;
#endif
      imnr = (int)((vertex->y - START_Y) / SLICE_THICKNESS + 0.5);
      if (imnr > imnr1) {
        imnr1 = imnr;
      }
      if (imnr < imnr0) {
        imnr0 = imnr;
      }
      if (version == 0) /* old surface format */
      {
        fread1(&num, fp);
        vertex->num = num;
        vertex->f = (int *)calloc(vertex->num, sizeof(int));
        if (!vertex->f) ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces", vertex->num);
        vertex->n = (uchar *)calloc(vertex->num, sizeof(uchar));
        if (!vertex->n) ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d nbrs", vertex->n);
        for (n = 0; n < vertex->num; n++) {
          fread3(&vertex->f[n], fp);
        }
      }
      else {
        vertex->num = 0;
      }
    }

    for (fno = 0; fno < mris->nfaces; fno += 2) {
      for (n = 0; n < 4; n++) /* read quandrangular face */
      {
        fread3(&vertices[n], fp);
      }

      /* 1st triangle */
      mris->faces[fno].v[0] = vertices[0];
      mris->faces[fno].v[1] = vertices[1];
      mris->faces[fno].v[2] = vertices[3];
      if (version < 0)
        for (n = 0; n < VERTICES_PER_FACE; n++) {
          mris->vertices[mris->faces[fno].v[n]].num++;
        }

      /* 2nd triangle */
      mris->faces[fno + 1].v[0] = vertices[2];
      mris->faces[fno + 1].v[1] = vertices[3];
      mris->faces[fno + 1].v[2] = vertices[1];
      if (version < 0)
        for (n = 0; n < VERTICES_PER_FACE; n++) {
          mris->vertices[mris->faces[fno + 1].v[n]].num++;
        }
    }
    fclose(fp);
  }
  strcpy(mris->fname, fname);
  {
    char *surf_name;

    surf_name = strrchr(fname, '/');
    if (surf_name == NULL) {
      surf_name = fname;
    }
    else {
      surf_name++; /* past the last slash */
    }
    if (toupper(*surf_name) == 'R') {
      mris->hemisphere = RIGHT_HEMISPHERE;
    }
    else {
      mris->hemisphere = LEFT_HEMISPHERE;
    }
  }
  if ((version < 0) || type == MRIS_ASCII_TRIANGLE_FILE) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      vertex = &mris->vertices[vno];
      mris->vertices[vno].f = (int *)calloc(mris->vertices[vno].num, sizeof(int));
      if (!mris->vertices[vno].f)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d faces at %dth vertex",
                  fname,
                  vno,
                  mris->vertices[vno].num);

      mris->vertices[vno].n = (uchar *)calloc(mris->vertices[vno].num, sizeof(uchar));
      if (!mris->vertices[vno].n)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d indices at %dth vertex",
                  fname,
                  vno,
                  mris->vertices[vno].num);
      mris->vertices[vno].num = 0;
    }
    for (fno = 0; fno < mris->nfaces; fno++) {
      face = &mris->faces[fno];
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        vertex = &mris->vertices[face->v[n]];
        vertex->f[vertex->num++] = fno;
      }
    }
  }

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].curv = 0;
    mris->vertices[vno].origarea = -1;
    mris->vertices[vno].border = 0;
#if 0
    mris->vertices[vno].origripflag = 0;
    mris->vertices[vno].ripflag = 0;
    mris->vertices[vno].val = 0;
    mris->vertices[vno].dist = 0;
    mris->vertices[vno].mx = 0;
    mris->vertices[vno].my = 0;
    mris->vertices[vno].mz = 0;
    mris->vertices[vno].fieldsign = 0;
    mris->vertices[vno].fsmask = 1;
    mris->vertices[vno].nc = 0;
    mris->vertices[vno].marked = 0;
#endif
    for (n = 0; n < mris->vertices[vno].num; n++) {
      for (m = 0; m < VERTICES_PER_FACE; m++) {
        if (mris->faces[mris->vertices[vno].f[n]].v[m] == vno) {
          mris->vertices[vno].n[n] = m;
        }
      }
    }
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;
    if (x > xhi) {
      xhi = x;
    }
    if (x < xlo) {
      xlo = x;
    }
    if (y > yhi) {
      yhi = y;
    }
    if (y < ylo) {
      ylo = y;
    }
    if (z > zhi) {
      zhi = z;
    }
    if (z < zlo) {
      zlo = z;
    }
  }
  mris->xlo = xlo;
  mris->ylo = ylo;
  mris->zlo = zlo;
  mris->xhi = xhi;
  mris->yhi = yhi;
  mris->zhi = zhi;
  mris->xctr = (xhi + xlo) / 2;
  mris->yctr = (yhi + ylo) / 2;
  mris->zctr = (zhi + zlo) / 2;
  mrisFindNeighbors(mris);
  MRIScomputeNormals(mris);
#if 0
  mrisComputeVertexDistances(mris) ;
  mrisReadTransform(mris, fname) ;
#endif
  if (type == MRIS_ASCII_TRIANGLE_FILE || type == MRIS_GEO_TRIANGLE_FILE) {
    MRISsetNeighborhoodSizeAndDist(mris, 2);
    MRIScomputeSecondFundamentalForm(mris);
    MRISuseMeanCurvature(mris);
  }
  else {
    if (MRISreadBinaryCurvature(mris, fname) != NO_ERROR) {
      fprintf(stdout, "ignoring curvature file...\n"); /*return(NULL) ;*/
    }
#if 0
    if (MRISreadBinaryAreas(mris, fname) != NO_ERROR)
    {
      return(NULL) ;
    }
#endif
  }

#if 0
  if (IS_QUADRANGULAR(mris))
  {
    MRISremoveTriangleLinks(mris) ;
  }
#endif
  MRISstoreCurrentPositions(mris);
  return (mris);
#endif
}

#endif


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisDebugVertex(MRI_SURFACE *mris, int vno)
{
  int     n ;
  VERTEX  *v, *vn ;
  float   d, dx, dy, dz ;

  v = &mris->vertices[vno] ;
  fprintf(stdout,
          "vertex #%d @ (%2.2f, %2.2f, %2.2f), n = (%2.2f, %2.2f, %2.2f) "
          "(%2.2f, %2.2f, %2.2f), val=%2.1f\n",
          vno, v->x, v->y, v->z, v->nx, v->ny, v->nz, v->dx, v->dy, v->dz,
          v->val) ;

  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    dx = vn->x - v->x ;
    dy = vn->y - v->y ;
    dz = vn->z - v->z ;
    d = sqrt(dx*dx + dy*dy + dz*dz) ;
    fprintf(stdout,
            "\tn %d: %6.6d, delta = (%2.3f, %2.3f, %2.3f), dist = %2.3f "
            "(%2.2f, %2.2f, %2.2f), val=%2.1f\n",
            n, v->v[n], dx, dy, dz, d, vn->dx, vn->dy, vn->dz, vn->val) ;
  }
  return(NO_ERROR) ;
}
#endif

#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
mrisAddVertices(MRIS *mris, double thresh)
{
  double   dist ;
  int      vno, nadded, n,nvertices, nfaces, nedges, eno ;
  VERTEX   *v, *vn ;
  float    x, y, z ;

  /* make it squared so we don't need sqrts later */
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "dividing edges more than %2.2f mm long.\n", thresh) ;
  }
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    x = v->origx ;
    y = v->origy ;
    z = v->origz ;

    /*
      only add vertices if average neighbor vector is in
      normal direction, that is, if the region is concave or sulcal.
    */
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      dist = sqrt(SQR(vn->origx-x) + SQR(vn->origy - y) + SQR(vn->origz - z));
      if (dist > thresh)
      {
        if (mrisDivideEdge(mris, vno, v->v[n]) == NO_ERROR)
        {
          nadded++ ;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "%d vertices added: # of vertices=%d, # of faces=%d.\n",
            nadded, mris->nvertices, mris->nfaces) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
  }
  return(nadded) ;
}
#endif
