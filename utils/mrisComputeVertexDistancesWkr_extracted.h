static bool FUNCTION_NAME(MRIS *mris, int new_dist_nsize, bool check)
{
  cheapAssert(0 < new_dist_nsize);
  cheapAssert(new_dist_nsize <= 3);
  
  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d %s stdout ",__FILE__,__LINE__, 
        FUNCTION_NAME == mrisComputeVertexDistancesWkr 
        ? "mrisComputeVertexDistancesWkr"
        : "mrisComputeOriginalVertexDistancesWkr");
    mris_print_hash(stdout, mris, "mris ", "\n");
  }

  int nonZeroInputXCount = 0;
  
  int errors = 0;
  
  int vno;

  switch (mris->status) {
    default: /* don't really know what to do in other cases */

    case MRIS_PLANE:

      ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(shown_reproducible) reduction(+:errors)  reduction(+:nonZeroInputXCount)
#endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin
    
        if (vno == Gdiag_no) DiagBreak();

        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
        VERTEX                * const v  = &mris->vertices         [vno];
        if (v->ripflag) continue;

        if (!check) OUTPUT_MAKER(mris,vno);
        else if (!v->OUTPUT_DIST) continue;

	if (v->INPUT_X != 0.0f) nonZeroInputXCount++;

        int *pv;
        int const vtotal = check ? VERTEXvnum(vt,new_dist_nsize) : vt->vtotal;
        int n;
        for (pv = vt->v, n = 0; n < vtotal; n++) {
          VERTEX const * const vn = &mris->vertices[*pv++];
          // if (vn->ripflag) continue;
          float xd = v->INPUT_X - vn->INPUT_X;
          float yd = v->INPUT_Y - vn->INPUT_Y;
          float zd = v->INPUT_Z - vn->INPUT_Z;
          float d = xd * xd + yd * yd + zd * zd;

          if (!check) 
            v->OUTPUT_DIST[n] = sqrt(d);
          else if (v->OUTPUT_DIST[n] != (float)(sqrt(d))) 
            errors++;
        }

        ROMP_PFLB_end
      }
      ROMP_PF_end
  
      break;

    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE: {

      ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(shown_reproducible) reduction(+:errors) reduction(+:nonZeroInputXCount)
#endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin
    
        if (vno == Gdiag_no) DiagBreak();

        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
        VERTEX          const * const v  = &mris->vertices         [vno];
        if (v->ripflag) continue;

        if (!check) OUTPUT_MAKER(mris,vno);
        else if (!v->OUTPUT_DIST) continue;

	if (v->INPUT_X != 0.0f) nonZeroInputXCount++;

        XYZ xyz1_normalized;
        float xyz1_length;
        XYZ_NORMALIZED_LOAD(&xyz1_normalized, &xyz1_length, v->INPUT_X, v->INPUT_Y, v->INPUT_Z);  // length 1 along radius vector

        float const radius = xyz1_length;

        int *pv;
        int const vtotal = check ? VERTEXvnum(vt,new_dist_nsize) : vt->vtotal;
        int n;
        for (pv = vt->v, n = 0; n < vtotal; n++) {
          VERTEX const * const vn = &mris->vertices[*pv++];
          if (vn->ripflag) continue;
          
          float angle = fabs(XYZApproxAngle(&xyz1_normalized, vn->INPUT_X, vn->INPUT_Y, vn->INPUT_Z));
            // radians, so 2pi around the circumference

          float d = angle * radius;
            // the length of the arc, rather than the straight line distance

          if (!check) 
            v->OUTPUT_DIST[n] = d;
          else if (v->OUTPUT_DIST[n] != (float)(d)) 
            errors++;
        }
        ROMP_PFLB_end
      }
      ROMP_PF_end

      break;
    }
  }

  if (nonZeroInputXCount == 0) {
    copeWithLogicProblem(NULL, "all xyz inputs zero - cant fix");
  }
  
  return (errors == 0);
}

#undef FUNCTION_NAME 
#undef INPUT_X 
#undef INPUT_Y 
#undef INPUT_Z
#undef OUTPUT_DIST
#undef OUTPUT_MAKER
