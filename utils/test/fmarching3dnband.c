/* fmarching3dnband.c : functions to implement the Fast Marching Method in 3D */
/* In this version, only compute the signed distance function till a certain narrowband */

/* Fast Marching Method to compute signed distance function in 3D Cube. 
 * Assume the initial contour is given by an array Ori, whose values 
 * are zero at the contour,	
 * +1 outside, and -1 inside. Fast marching is then used to build 
 * the signed distance function		
    --------------------> j (XN)
    |
    |
    |
    V 
    i (YN)
 *  k (ZN) indicates the depth
 */

#include "fmarching3dnband.h"
#define IFij  1  /* 1/F[k][i][j] */
#define ISFij 1  /* 1/(F[k][i][j]*F[k][i][j]) */
#define MYMIN(a,b) (a <b? a: b)

float ReCompute(float Nv,float Sv,float Ev,float Wv,float Fv, float Bv, 
		unsigned char Nl, unsigned char Sl, unsigned char El, unsigned char Wl, unsigned char Fl, unsigned char Bl){
/* Compute the new value at a NarrowBand point using values of its 
 * ALIVE(Accepted) 6-connected neighbours
 * Note: at least one of its 6 neighbours must be ALIVE 
 * If ALIVE neighbours exist only in the north-south (or east-west, or
 * front-back) neighbours, the updated value is equal
 * to the value of the smaller one plus 1 (1/F(k,i,j) if F(k,i,j) ne 1). 
 * Otherwise, a quadratic equation
 * need to be solved, and the bigger root is taken as the updated value
 */
/* The following program assumes F(k,i,j) = 1; if not, replace IFij and
 * ISFij with the true values. Note, F(k,i,j) should be positive!
 */ 

  /* Suppose a, b, and c are the minimum T value in three directions */	
  float s, s2; /* s = a + b +c; s2 = a*a + b*b +c*c */
  float tmp;
  int count;
  
  s = 0; s2 =0; count = 0;
  
  if (Nl == (unsigned char)ALIVE && Sl == (unsigned char)ALIVE){ 
    tmp = MYMIN(Nv, Sv); /* Take the smaller one if both ALIVE */
    s += tmp;
    s2 += tmp*tmp;
    count++;
  }
  else if (Nl == (unsigned char)ALIVE){
    s += Nv;	 /* Else, take the ALIVE one */
    s2 += Nv*Nv;
    count++;
  } 
  else if (Sl == (unsigned char)ALIVE){
    s += Sv;
    s2 += Sv*Sv;
    count++;
  }
  
  /* Similarly in the east-west direction to get correct 
     approximation to the derivative in the x-direction */
  if (El == (unsigned char)ALIVE && Wl == (unsigned char)ALIVE){ 
    tmp = MYMIN(Ev, Wv); /* Take the smaller one if both ALIVE */
    s += tmp;
    s2 += tmp*tmp;
    count++;
  }
  else if (El == (unsigned char)ALIVE){
    s += Ev;	 /* Else, take the ALIVE one */
    s2 += Ev*Ev;
    count++;
  } 
  else if (Wl == (unsigned char)ALIVE){
    s += Wv;
    s2 += Wv*Wv;
    count++;
  }

  /* Similarly in the front-back direction to get correct 
     approximation to the derivative in the z-direction */
  if (Fl == (unsigned char)ALIVE && Bl == (unsigned char)ALIVE){ 
    tmp = MYMIN(Fv, Bv); /* Take the smaller one if both ALIVE */
    s += tmp;
    s2 += tmp*tmp;
    count++;
  }
  else if (Fl == (unsigned char)ALIVE){
    s += Fv;	 /* Else, take the ALIVE one */
    s2 += Fv*Fv;
    count++;
  } 
  else if (Bl == (unsigned char)ALIVE){
    s += Bv;
    s2 += Bv*Bv;
    count++;
  }
  
  /* count must be greater than zero since there must be one ALIVE pt 
     in the neighbors */
  
  tmp = (s + (float)sqrt((double)(s*s - count*(s2 - ISFij))))/count; 
  /* The larger root */
  return tmp;
}

void fmarching3d(MRI *Ori, MRI *T, float Thred){
  /* Ori :Original level set function embedding the original contour 
   *      location 1: outside; -1: inside 
   * T	 :Signed distance function computed using the Fast Maching 
   *      Levelset Method
   * XN	 :Image width in 3D
   * YN	 :Image height in 3D
   * ZN	 :Image depth in 3D
   * Thred: distance threshold to stop marching
   */
  
  /* The unsigned distance function is computed first, then it is 
   * multiplied by the sign of Ori to get the signed distance function 
   */

  /* Offset indexing for 6-connected neighbours */
  /* East, South, West, North, Front, Back */
  static int xoff[6] = {1, 0, -1, 0, 0, 0}; 
  static int yoff[6] = {0, 1, 0, -1, 0, 0};
  static int zoff[6] = {0, 0, 0, 0, 1, -1};
  int XN, YN, ZN;

  MRI *label; /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */
  MRI *BackPointer;	
  /* 3D cube storing the backpointer to the narrowband heap */
  
  Xheap H;				/* Narrrow band heap */
  XheapElement he;  
  
  int heapindex;
  
  int i,j,d;			/* Iteration index, d for depth */
  int koff;				/* Neighbourhood index */
  int newi,newj, newd;
  float Nv,Sv,Wv,Ev, Fv, Bv;    /* Value at six neighours of a pixel */
  float Cv;									 /* Value at local point */	
  unsigned char Nl, Sl, Wl, El, Fl, Bl; /* Label at six neighours of a pixel */
  
  int NSFlag, WEFlag, FBFlag; 
       /* Logical variables aid in locating sign change */
  float s,t,w;		 /* Used in initialization */
  float result;          /* 1/s*s + 1/t*t + 1/w*w */		

  float value, newvalue, oldvalue;
  
  int LX,HX,LY,HY,LZ,HZ; /* Logical variables aid in dealing with boundary */ 

  XN = Ori->width; YN = Ori->height; ZN = Ori->depth;

  /* Memeory allocation */
  label = MRIallocSequence(Ori->width, Ori->height,Ori->depth, 
                       MRI_UCHAR, Ori->nframes);

  MRIcopyHeader(Ori, label);

  BackPointer = MRIallocSequence(Ori->width, Ori->height,Ori->depth, 
                       MRI_INT, Ori->nframes);
  MRIcopyHeader(Ori, BackPointer);

  /* Initialize heap structure */
  H = xhInitEmpty();
  
  /* Initialization for marching inwards and outwards simultaneously */
  /* It is assumed the original contour lies exactly at grid points. 
   * If not, more complicated initialization method need to be performed
   */
  for(d = 0; d < ZN; d++){
    for(i=0; i<YN; i++){
      for(j=0; j<XN;j++){
	MRIseq_vox(label, j, i, d, 0) = (unsigned char) FAWAY;/*All points are labelled as FAR AWAY */
      }
    }
  }
  
  /* First, find grid points that lie no more than 1 pixel away from the 
   * contour, label them as ALIVE
   */
  for(d=0; d<ZN;d++){
    for(i=0; i<YN; i++){
      for (j=0; j<XN;j++){
	value = MRIgetVoxVal(Ori, j, i, d, 0); 
	if( value < 0.0000001 && value > -0.0000001)
	  { /* This grid point lies exactly on the contour */
	    MRIsetVoxVal(T, j,  i, d, 0, 0); /* Its distance value is exactly zero */
	    MRIseq_vox(label, j, i, d, 0) =  (unsigned char)ALIVE;
	  }
	else{ /* Judge whether it's within one pixel away from the contour */
	  LY = (i==0); HY = (i==(YN-1));
	  LX = (j==0); HX = (j==(XN-1));
	  LZ = (d==0); HZ = (d==(ZN-1));
	  
	  NSFlag = 0; WEFlag = 0; FBFlag = 0;
	  
	  Nv = MRIgetVoxVal(Ori, j, i-1+LY, d, 0); // Ori[d][i-1+LY][j];
	  Sv = MRIgetVoxVal(Ori, j, i+1-HY, d, 0); //Ori[d][i+1-HY][j];
	  Wv = MRIgetVoxVal(Ori, j-1+LX, i, d, 0); // Ori[d][i][j-1+LX];
	  Ev = MRIgetVoxVal(Ori, j+1-HX, i, d, 0); //Ori[d][i][j+1-HX];
	  Fv = MRIgetVoxVal(Ori, j, i, d+1-HZ, 0); //Ori[d+1-HZ][i][j];
	  Bv = MRIgetVoxVal(Ori, j, i, d-1+LZ, 0); //Ori[d-1+LZ][i][j];
	  Cv = MRIgetVoxVal(Ori, j, i, d, 0); // Ori[d][i][j];
	  
	  
	  if(Nv*Cv < 0){
	    NSFlag = 1;
	    s = Nv;
	  }
	  if(Sv*Cv < 0){
	    if(NSFlag == 0){
	      NSFlag = 1;
	      s = Sv;
	    }else{
	      s = (fabs(Nv) > fabs(Sv))? Nv : Sv;
	    }  
	  }	
	  
	  if(Wv*Cv < 0){
	    WEFlag = 1;
	    t = Wv;
	  }
	  if(Ev*Cv < 0){
	    if(WEFlag == 0){
	      WEFlag = 1;
	      t = Ev;
	    }else{
	      t = (fabs(Ev) > fabs(Wv))? Ev : Wv;
	    }  
	  }
	  
	  if(Fv*Cv < 0){
	    FBFlag = 1;
	    w = Fv;
	  }
	  if(Bv*Cv < 0){
	    if(FBFlag == 0){
	      FBFlag = 1;
	      w = Bv;
	    }else{
	      w = (fabs(Fv) > fabs(Bv))? Fv : Bv;
	    }  
	  }
	  
	  result = 0;
	  if(NSFlag){
	    s = Cv/(Cv-s);
	    result += 1/(s*s);
	  }
	  if(WEFlag){
	    t = Cv/(Cv-t);
	    result += 1/(t*t);
	  }
	  if(FBFlag){
	    w = Cv/(Cv-w);
	    result += 1/(w*w);
	  }
	  
	  if(result == 0) continue; /* Not within 1 pixel away */

	  MRIseq_vox(label, j, i, d, 0) = (unsigned char)ALIVE; /*Within 1 pixel away, put into ALIVE */
	  
	  result = (float)sqrt((double)result);
	  
	  MRIsetVoxVal(T, j, i, d, 0, IFij/result);
	}
      }
    }
  }
  
  /* Initialize NarrowBand Heap */
  for(d=0; d<ZN;d++){
    for(i=0; i<YN; i++){
      for(j=0; j<XN;j++){
	if(MRIseq_vox(label,j,i,d,0) != (unsigned char)ALIVE) continue;
	
	/* Put its 6 neighbors into NarrowBand */
	
	for(koff=0; koff <= 5; koff++){/*Find six neighbouring points */ 
	  newi = i + yoff[koff];
	  newj = j + xoff[koff];
	  newd = d + zoff[koff];
	  
	  if(newi <0 || newi >= YN || newj < 0 || 
	     newj >= XN || newd < 0 || newd >= ZN) 
	    continue; /* Out of computational Boundary*/
	  
	  if(MRIseq_vox(label, newj, newi, newd, 0) != (unsigned char)FAWAY) continue;
	  
	  /* Compute the distance at this point and add to the heap */
	  MRIseq_vox(label, newj, newi, newd, 0) = (unsigned char)NBAND;
	  /* Note: Only ALIVE points contribute to the distance computation */
	  /* Neighbour to the north*/
	  if(newi > 0) { 
	    Nv = MRIgetVoxVal(T, newj, newi-1, newd,0); // T[newd][newi-1][newj];
	    Nl = MRIseq_vox(label, newj, newi-1, newd,0); //label[newd][newi-1][newj];
	  }else Nl = 0;
	  
	  /* Neighbour to the south*/
	  if(newi < YN-1) {
	    Sv = MRIgetVoxVal(T, newj, newi+1, newd,0); // T[newd][newi+1][newj];
	    Sl = MRIseq_vox(label, newj, newi+1, newd,0); // label[newd][newi+1][newj];
	  } else Sl = 0;
	  
	  /* Neighbour to the east*/
	  if(newj < XN-1) {
	    Ev = MRIgetVoxVal(T, newj+1, newi, newd,0); // T[newd][newi][newj+1];
	    El = MRIseq_vox(label, newj+1, newi, newd,0); //label[newd][newi][newj+1];
	  }else El = 0;
	  
	  /*Neighbour to the west*/
	  if(newj > 0){
	    Wv = MRIgetVoxVal(T, newj-1, newi, newd,0);// T[newd][newi][newj-1];
	    Wl = MRIseq_vox(label, newj-1, newi, newd,0); // label[newd][newi][newj-1];
	  }else Wl = 0;
	  
	  /* Neighbour to the front */
	  if(newd < ZN-1) {
	    Fv = MRIgetVoxVal(T, newj, newi, newd+1,0); // T[newd+1][newi][newj];
	    Fl = MRIseq_vox(label, newj, newi, newd+1,0);// label[newd+1][newi][newj];
	  }else Fl = 0;
	  
	  /*Neighbour to the back */
	  if(newd > 0){
	    Bv = MRIgetVoxVal(T, newj, newi, newd-1,0); // T[newd-1][newi][newj];
	    Bl = MRIseq_vox(label, newj, newi, newd-1,0); //label[newd-1][newi][newj];
	  }else Bl = 0;
	  
	  /*Update the value of this to-be-updated NarrowBand point*/
	  newvalue = ReCompute(Nv,Sv,Ev,Wv,Fv,Bv,Nl,Sl,El,Wl,Fl,Bl);	
	  MRIsetVoxVal(T, newj, newi, newd, 0, newvalue);
	  xhInsert(newvalue, newj, newi, newd, &MRIIseq_vox(BackPointer, newj, newi, newd, 0), H);
	  
	}
	
      }
    }
  }
  /* End of Initialization */
  
  /* Begin Fast Marching to get the unsigned distance function inwords 
   * and outwards simultaneously
   * since points inside and outside the contour won't interfere with 
   * each other
   */
  while(!xhIsEmpty(H)){ /* There are still points not yet accepted */
    he = xhRemove(H); /* Label the point with smallest value among all 
			 NarrowBand points as ALIVE */
    
    /* Put the smallest heap element to ALIVE */
    d = he.z;
    i = he.y;
    j = he.x;

    if(he.value > Thred) break;

    MRIsetVoxVal(T, j, i, d, 0, he.value);
    MRIseq_vox(label, j, i,d , 0) = (unsigned char)ALIVE; 
    
    /* Update its neighbor */
    /* Put FARAWAY neighbour into NarrowBand, Recompute values at 
     * NarrowBand neighbours,
     * Keep ALIVE (Accepted) neighbour unchanged
     */
    for(koff=0; koff <= 5; koff++){
      newi = i + yoff[koff];
      newj = j + xoff[koff];
      newd = d + zoff[koff];

      if(newi <0 || newi >= YN || newj < 0 || newj >= XN || newd <0 || newd >=ZN) 
	continue; /* Out of boundary */
      // if(label[newd][newi][newj] == (unsigned char)ALIVE) 
      if(MRIseq_vox(label, newj, newi, newd, 0) == (unsigned char)ALIVE) 
	continue; /* Don't change ALIVE neighbour */
      
      /* ReCompute the value at (newj, newi, newd) */
      /* Get the values and labels of six neighbours of the to-be-updated 
       * point. The labels are needed since only values at ALIVE 
       * neighbours will be used to update
       * the value of the to-be-updated point
       */
      
      /* Note: Only ALIVE points contribute to the distance computation */
      /* Neighbour to the north */
      if(newi > 0) { 
	Nv = MRIgetVoxVal(T, newj, newi-1, newd,0); // T[newd][newi-1][newj];
	Nl = MRIseq_vox(label, newj, newi-1, newd,0); //label[newd][newi-1][newj];
      }else Nl = 0;
      
      /* Neighbour to the south*/
      if(newi < YN-1) {
	Sv = MRIgetVoxVal(T, newj, newi+1, newd,0); // T[newd][newi+1][newj];
	Sl = MRIseq_vox(label, newj, newi+1, newd,0); // label[newd][newi+1][newj];
      } else Sl = 0;
      
      /* Neighbour to the east*/
      if(newj < XN-1) {
	Ev = MRIgetVoxVal(T, newj+1, newi, newd,0); // T[newd][newi][newj+1];
	El = MRIseq_vox(label, newj+1, newi, newd,0); //label[newd][newi][newj+1];
      }else El = 0;
      
      /*Neighbour to the west*/
      if(newj > 0){
	Wv = MRIgetVoxVal(T, newj-1, newi, newd,0);// T[newd][newi][newj-1];
	Wl = MRIseq_vox(label, newj-1, newi, newd,0); // label[newd][newi][newj-1];
      }else Wl = 0;
      
      /* Neighbour to the front */
      if(newd < ZN-1) {
	Fv = MRIgetVoxVal(T, newj, newi, newd+1,0); // T[newd+1][newi][newj];
	Fl = MRIseq_vox(label, newj, newi, newd+1,0);// label[newd+1][newi][newj];
      }else Fl = 0;
      
      /*Neighbour to the back */
      if(newd > 0){
	Bv = MRIgetVoxVal(T, newj, newi, newd-1,0); // T[newd-1][newi][newj];
	Bl = MRIseq_vox(label, newj, newi, newd-1,0); //label[newd-1][newi][newj];
      }else Bl = 0;
      
      
      /* Update the value of this to-be-updated NarrowBand point */
      newvalue = ReCompute(Nv,Sv,Ev,Wv,Fv,Bv,Nl,Sl,El,Wl,Fl,Bl);	
      
      /* If it was a FARAWAY point, add it to the NarrowBand Heap; 
       * otherwise, just update its value
       * using the backpointer
       */
      if(MRIseq_vox(label, newj, newi, newd, 0) == (unsigned char)NBAND) 
	xhChangeValue(MRIIseq_vox(BackPointer, newj, newi, newd, 0),newvalue, H);
      else{
	xhInsert(newvalue, newj, newi, newd, &(MRIIseq_vox(BackPointer, newj, newi, newd,0)), H);
	MRIseq_vox(label, newj, newi, newd, 0) = (unsigned char)NBAND;
      }
    } /* End of updating 6 neighbours */
    
  }/* End of marching loop */
  
  /* Add signs to the unsigned distance function */
  for(d=0; d<ZN; d++){
    for(i=0; i<YN; i++){
      for(j=0; j<XN;j++){
	if(MRIseq_vox(label, j, i, d , 0) != (unsigned char)ALIVE) 
	  MRIsetVoxVal(T, j, i, d, 0, Thred);
	
	if(MRIgetVoxVal(Ori, j, i,d ,0) < 0 ){ 
	  newvalue = -1*MRIgetVoxVal(T, j, i, d, 0);
	  MRIsetVoxVal(T, j, i, d, 0, newvalue);
	}
      }
    }
  }

  /* Free memory */
  MRIfree(&label);
  MRIfree(&BackPointer);
  
  xhDestroy(H);
  
  return;
}










