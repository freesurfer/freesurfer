static bool FUNCTION_NAME(
#ifdef COMPILING_MRIS_MP
    MRIS_MP* mris, 
#else
    MRIS*    mris, 
#endif
    int new_dist_nsize, 
    bool check)
{
  cheapAssert(0 < new_dist_nsize);
  cheapAssert(new_dist_nsize <= 3);
  
#ifndef COMPILING_MRIS_MP
  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d %s ",__FILE__,__LINE__, __MYFUNCTION__);
    mris_print_hash(stdout, mris, "mris ", "\n");
  }
#endif

  if (false && mris->status != mris->INPUT_STATUS) {
    auto oldStatusFormula = MRIS_Status_distanceFormula(mris->status);
    auto newStatusFormula = MRIS_Status_distanceFormula(mris->INPUT_STATUS);
    fprintf(stdout, "%s:%d %s mris->status:%s != mris->INPUT_STATUS:%s\n",__FILE__,__LINE__, __MYFUNCTION__,
      MRIS_Status_text(mris->status), 
      MRIS_Status_text(mris->INPUT_STATUS));
    if (oldStatusFormula != newStatusFormula) 
      fprintf(stdout, "  executing formula_%d instead of _%d\n", newStatusFormula, oldStatusFormula);
    else
      fprintf(stdout, "  but using the same formula_%d so not important\n", oldStatusFormula);
  }

  // if v_dist_buffer isn't allocated here, OUTPUT_MAKER produces memory leaks when multithreaded
#ifdef COMPILING_MRIS_MP
  if (!mris->v_dist_buffer) mris->v_dist_buffer = (float**)calloc(mris->nvertices, sizeof(float*));
#endif

  int nonZeroInputXCount = 0;
  int errors = 0;
  int vno;

  switch (MRIS_Status_distanceFormula(mris->INPUT_STATUS)) {
    default: 
      cheapAssert(false);

    case MRIS_Status_DistanceFormula_0:
    {
      ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(shown_reproducible) reduction(+:errors)  reduction(+:nonZeroInputXCount)
#endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin
    
        if (vno == Gdiag_no) DiagBreak();

        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
#ifndef COMPILING_MRIS_MP
        VERTEX                * const v  = &mris->vertices         [vno];
        if (v->ripflag) continue;
        float const vx = v->INPUT_X, vy = v->INPUT_Y, vz = v->INPUT_Z;
#else
        if (mris->v_ripflag[vno]) continue;
        float const vx = mris->INPUT_X[vno], vy = mris->INPUT_Y[vno], vz = mris->INPUT_Z[vno];
#endif

        if (!check) OUTPUT_MAKER(mris,vno);

#ifndef COMPILING_MRIS_MP
        float* const output_dist = v->OUTPUT_DIST;
#else
        float* const output_dist = mris->v_dist[vno];
#endif
        if (check && !output_dist) continue;

	if (vx != 0.0f) nonZeroInputXCount++;

        int const vtotal = check ? VERTEXvnum(vt,new_dist_nsize) : vt->vtotal;

        int n;
        for (n = 0; n < vtotal; n++) {
          int const vnon = vt->v[n];
#ifndef COMPILING_MRIS_MP
          VERTEX const * const vn = &mris->vertices[vnon];
          float const vnx = vn->INPUT_X, vny = vn->INPUT_Y, vnz = vn->INPUT_Z;
#else
          float const vnx = mris->INPUT_X[vnon], vny = mris->INPUT_Y[vnon], vnz = mris->INPUT_Z[vnon];
#endif
          float xd = vx - vnx;
          float yd = vy - vny;
          float zd = vz - vnz;
          float d = xd * xd + yd * yd + zd * zd;

          if (!check) 
            output_dist[n] = sqrt(d);
          else if (output_dist[n] != (float)(sqrt(d))) 
            errors++;

        }

        ROMP_PFLB_end
      }
      ROMP_PF_end
    } break;

    case MRIS_Status_DistanceFormula_1:
    {
      // Sped up by computing the normalized vectors for all the vertices once rather than repeatedly
      // and storing them along with their ripflag in a structure!
      //
      struct PerVertexInfo {
        char  ripped;
        int   rippedOrZeroLength;
        XYZ   xyz_unnormalized;
        float xyz_length;
        XYZ   xyz_normalized;
      }* vertexInfo = (struct PerVertexInfo*)malloc(mris->nvertices * sizeof(*vertexInfo));

      // In parallel, load the normalized xyz so only computed once per vertex
      // 
      ROMP_PF_begin
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(shown_reproducible) reduction(+:errors) reduction(+:nonZeroInputXCount)
#endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin

#ifndef COMPILING_MRIS_MP
        VERTEX const * const v = &mris->vertices[vno];
        char  const v_ripflag = (char)v->ripflag;
        float const vx = v->INPUT_X, vy = v->INPUT_Y, vz = v->INPUT_Z;
#else
        char  const v_ripflag = (char)mris->v_ripflag[vno];
        float const vx = mris->INPUT_X[vno], vy = mris->INPUT_Y[vno], vz = mris->INPUT_Z[vno];
#endif
        vertexInfo[vno].ripped             = v_ripflag;
        vertexInfo[vno].rippedOrZeroLength = v_ripflag;

        if (v_ripflag) continue;

	if (vx != 0.0f) nonZeroInputXCount++;

        vertexInfo[vno].xyz_unnormalized.x = vx;
        vertexInfo[vno].xyz_unnormalized.y = vy;
        vertexInfo[vno].xyz_unnormalized.z = vz;
          // length xyz_length
          
        XYZ_NORMALIZED_LOAD(&vertexInfo[vno].xyz_normalized, &vertexInfo[vno].xyz_length, vx, vy, vz);  
          // length 1 along radius vector

        if (FZERO(vertexInfo[vno].xyz_length)) vertexInfo[vno].rippedOrZeroLength = (char)true;

        ROMP_PFLB_end
      }
      ROMP_PF_end

      // Do the edges
      //
      ROMP_PF_begin
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(shown_reproducible) reduction(+:errors)
#endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin

        struct PerVertexInfo* vi = &vertexInfo[vno];

        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
#ifndef COMPILING_MRIS_MP
        VERTEX const * const v = &mris->vertices[vno];
#else
#endif

        if (vi->ripped) continue;
          // Non-ripflag with a radius of zero still need to get their distances set to zero

        if (!check) OUTPUT_MAKER(mris,vno);

#ifndef COMPILING_MRIS_MP
        float* const output_dist = v->OUTPUT_DIST;
#else
        float* const output_dist = mris->v_dist[vno];
#endif
        if (check && !output_dist) continue;


        float const radius = vi->xyz_length;

        int const vtotal = check ? VERTEXvnum(vt,new_dist_nsize) : vt->vtotal;
        int n;
        for (n = 0; n < vtotal; n++) {
          int const vno2 = vt->v[n];
          struct PerVertexInfo* vi2 = &vertexInfo[vno2];

          float d = 0.0;

          if (!vi2->rippedOrZeroLength) {

            float angle = 
              fabs(XYZApproxAngle_knownLength(
                  &vi->xyz_normalized, 
                  vi2->xyz_unnormalized.x, vi2->xyz_unnormalized.y, vi2->xyz_unnormalized.z, vi2->xyz_length));
              // radians, so 2pi around the circumference

            d = angle * radius;
              // the length of the arc, rather than the straight line distance
              // but this definition is still weird because the two points are not equal radius from center
              // so the distance from a to b might not be the same as b to a!
              // the vertices are usually very close and far from the origin, making this almost irrelevant.
          }

          if (!check) 
            output_dist[n] = d;
          else if (output_dist[n] != (float)(d)) 
            errors++;
        }
        ROMP_PFLB_end
      }
      ROMP_PF_end

      free(vertexInfo);

      break;
    }
  }

  if (nonZeroInputXCount == 0) {
    copeWithLogicProblem(NULL, "all xyz inputs zero - cant fix");
  }
  
  return (errors == 0);
}

#undef FUNCTION_NAME 
#undef INPUT_STATUS
#undef INPUT_X 
#undef INPUT_Y 
#undef INPUT_Z
#undef OUTPUT_DIST
#undef OUTPUT_MAKER
