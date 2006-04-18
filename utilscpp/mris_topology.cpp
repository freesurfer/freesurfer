#include "mris_topology.h"

extern "C" bool MRIScorrectDefect(MRIS *mris, int defect_number){
	
	return true; 
}


extern "C"  int MRISgetEulerNumber(const MRIS *mris, const int *list_of_faces, int nfs){

	int nv,nf,ne;

	//we need to allocate the edge structure for the faces in the list
	//mark and count the vertices
	for(int n = 0 ; n < mris->nvertices ; n++){
		mris->vertices[n].marked = 0;
		mris->vertices[n].e=0;
	}
	nv = 0 ;
	for(int n = 0 ; n < nfs ; n++)
		for(int i = 0 ; i < 3 ; i++){
			VERTEX *v=&mris->vertices[mris->faces[list_of_faces[n]].v[i]];
			if(v->marked==0){
				nv++;
				v->marked=1;
				if(v->e) delete [] v->e; 
				v->e = (int*)calloc(v->vnum,sizeof(int));
				//for(int p = 0 ; p < v->vnum ; p++)
				//v->e[p]=0;
			};
		}

	//mark and count the faces
	for(int n = 0 ; n < mris->nfaces ; n++)
		mris->faces[n].marked=0;
	nf=ne=0;
	for(int n = 0 ; n < nfs ; n++){
		FACE *face= &mris->faces[list_of_faces[n]];
		if(face->marked==0){
			nf++;
			face->marked=1;
			int vn0,vn1;
			VERTEX *v0,*v1;
			for(int i = 0 ; i < 3 ; i++){
				vn0 = face->v[i];
				if(i==2) vn1 = face->v[0];
				else vn1 = face->v[i+1];
				v0 = &mris->vertices[vn0];
				v1 = &mris->vertices[vn1];
				//edge vn0 <--> vn1 ?
				for(int p = 0 ; p < v0->vnum ; p++){
					if(v0->v[p] == vn1 && v0->e[p]==0){
						ne++;
						v0->e[p]=1;
						//mark the other edge
						for(int m = 0 ; m < v1->vnum ; m++){
							if(v1->v[m]==vn0){
								v1->e[m]=1;
								break;
							}
						}
						break;
					}
				}
			}
		}
	}
	
	//free the edge structure - reset marks to 0
	for(int n = 0 ; n < nfs ; n++){
		FACE *face= &mris->faces[list_of_faces[n]];
		face->marked = 0;
		for(int i = 0 ; i < 3 ; i++){
			VERTEX *v=&mris->vertices[face->v[i]];
			if(v->marked){
				v->marked=0;
				if(v->e) delete [] v->e;
				v->e=0;
			}
		}
	}

	return (nv-ne+nf);
}

extern "C" MRIP* MRIPextractFromMRIS(MRIS *mris, int defect_number){
	int *list_of_faces, nfaces;

	/* count the number of faces and vertices */
	//counting the number of vertices
	int nvertices = 0 ; 
	for(int n = 0 ; n < mris->nvertices ; n++)
		if(mris->vertices[n].marked2 == defect_number)
			nvertices++;
	if(nvertices==0) return NULL;
	//counting the number of faces
	nfaces = 0 ;
	for(int n = 0 ; n < mris->nfaces ; n++){
		bool is_face = true;
		FACE *face=&mris->faces[n];
		face->marked=0;
		for(int i = 0 ; i < 3 ; i++)
			if(mris->vertices[face->v[i]].marked2 != defect_number){
				is_face = false;
				break;
			};
		if(is_face) {
			nfaces++;
			face->marked=1;
		}
	}
	list_of_faces = new int[nfaces];
	nfaces=0;
	for(int n = 0 ; n < mris->nfaces ; n++)
		if(mris->faces[n].marked) list_of_faces[nfaces++]=n;

	/* first compute euler number of defect : if one return  0 */
	int euler = MRISgetEulerNumber(mris, list_of_faces, nfaces);
	
	//checking the value of the euler number
	if(euler==1) {
		delete [] list_of_faces;
		return NULL;
	}
	if(euler%2) {
		//ici	
	};
	int ncorrections = (1-euler)/2;

	/* allocate the patch with extra space */
	int max_vertices,max_faces;
//	_OverAlloc(2*loop.npoints+2*pdisk->disk.nvertices,2*(2*loop.npoints+pdisk->disk.nfaces+pdisk->ring.npoints)); 
	max_vertices = nvertices + ncorrections*(nvertices + 2*33);
	max_faces = nfaces + ncorrections*(2*nvertices+ 52);

	MRIP *mrip = MRIPalloc(max_vertices, max_faces);
	mrip->n_vertices = nvertices;
	mrip->n_faces = nfaces;

	MRIS *mris_dst=mrip->mris;
	/* copy the necessary information */
	mris_dst->nvertices=nvertices;
	mris_dst->nfaces=nfaces;

	//the corresponding tables
	int *vt_to,*vt_from;
	vt_to = new int[max_vertices];
	for(int n = 0 ; n < max_vertices ; n++)
		vt_to[n]=-1;
	vt_from = new int[nvertices];

	nvertices=0;
	for(int n = 0 ; n < mris->nvertices ; n++)
		if(mris->vertices[n].marked2 == defect_number){
			VERTEX *vsrc = &mris->vertices[n];
			VERTEX *v = &mris_dst->vertices[nvertices];
			vt_from[n]=nvertices;
			vt_to[nvertices++]=n;
			//copy the strict necessary
			v->x = vsrc->x;
			v->y = vsrc->y;
			v->z = vsrc->z;
		};
	mrip->vtrans_to = vt_to;
	mrip->vtrans_from = vt_from;

	//the corresponding tables
	int *ft_to,*ft_from;
	ft_to = new int[max_faces];
	for(int n = 0 ; n < max_faces ; n++)
		ft_to[n]=-1;
	ft_from = new int[mris->nfaces];

	nfaces=0;
	for(int n = 0 ; n < mris->nfaces ; n++){
		bool is_face = true;
		for(int i = 0 ; i < 3 ; i++)
			if(mris->vertices[mris->faces[n].v[i]].marked2 != defect_number){
				is_face = false;
				break;
			};
		if(is_face){
			FACE *fsrc = &mris->faces[n];
			FACE *f = &mris_dst->faces[nfaces];
			ft_from[n]=nfaces;
			ft_to[nfaces++]=n;
			for(int i = 0 ; i < 3 ; i++)
				f->v[i] = vt_from[fsrc->v[i]];
		}
	}
	mrip->ftrans_to = ft_to;
	mrip->ftrans_from = ft_from;

	//init the rest
	MRISinitSurface(mris_dst);

	return mrip;
}

void MRISinitSurface(MRIS *mris){
	VERTEX *v;
	
	for(int n = 0 ; n < mris->nvertices ; n++){
		v=&mris->vertices[n];
		v->num=0;
		v->marked=0;
	}

	// counting the number of faces per vertex
	for(int n = 0 ; n < mris->nfaces ; n++)
		for(int i = 0 ; i < 3 ; i++)
			mris->vertices[mris->faces[n].v[i]].num++;
	// allocate the list of faces
	for(int n = 0 ; n < mris->nvertices ; n++){
		mris->vertices[n].f=(int *)calloc(mris->vertices[n].num,sizeof(int));
		mris->vertices[n].n=(uchar *)calloc(mris->vertices[n].num,sizeof(uchar));
		mris->vertices[n].num=0;
	}
	// initialize the list of faces
	for(int n = 0 ; n < mris->nfaces ; n++){
		for(int i = 0 ; i < 3 ; i++){
			int vno = mris->faces[n].v[i];
			v=&mris->vertices[vno];
			v->f[v->num]=n;
			v->n[v->num++]=i;
		}
	}

	// counting the list of vertices
	for(int n = 0 ; n < mris->nvertices ; n++){
		v=&mris->vertices[n];
		v->vnum=0;
		for(int p = 0 ; p < v->num ; p++){
			FACE *face=&mris->faces[v->f[p]];
			for( int i = 0 ; i < 3 ; i++){
				int vn=face->v[i];
				if(vn==n) continue;
				if(mris->vertices[vn].marked) continue;
				mris->vertices[vn].marked=1;
				v->vnum++;
			}
		}
		// allocate the list of vertices
		v->v=(int*)calloc(mris->vertices[n].vnum,sizeof(int));
		v->vnum=0;
		for(int p = 0 ; p < v->num ; p++){
			FACE *face=&mris->faces[v->f[p]];
			for( int i = 0 ; i < 3 ; i++){
				int vn=face->v[i];
				if(vn==n) continue;
				if(mris->vertices[vn].marked == 0 ) continue;
				mris->vertices[vn].marked=0;
				v->v[v->vnum++]=vn;
			}
		}
	}

}


MRIP *MRIPalloc(int nvertices, int nfaces){
	MRIP *patch;

	/* allocate the patch */
	patch = (MRIP*)calloc(1,sizeof(MRIP));
	/* allocate the surface */
	MRIS *mris = MRISalloc(nvertices,nfaces);
	for(int n = 0 ; n < mris->nvertices ; n++){ //making sure...
		VERTEX *v=&mris->vertices[n];
		v->f=NULL;
		v->n=NULL;
		v->v=NULL;
	}
	patch->mris=mris;

	return patch;
}

void MRIPfree(MRIP **mrip){
	MRIP *patch;

	patch = *mrip;
	*mrip = NULL;

	MRISfree(&patch->mris);
	if(patch->vtrans_to) free(patch->vtrans_to);
	if(patch->ftrans_to) free(patch->ftrans_to);
	if(patch->vtrans_from) free(patch->vtrans_from);
	if(patch->ftrans_from) free(patch->ftrans_from);
	free(patch);
}

MRIP *MRIPclone(MRIP *src){
	MRIP *dst;

	/* allocate the patch */
	dst = (MRIP*)calloc(1,sizeof(MRIP));
	/* clone the surface */
	dst->mris = MRISclone(src->mris);
	
	// do not clone the rest
	return dst;
}

Surface *MRIStoSurface(MRIS *mris){

	Surface *surface=new Surface(mris->max_vertices,mris->max_faces);
	surface->nvertices = mris->nvertices;
	surface->nfaces = mris->nfaces;

	for(int n = 0 ; n < mris->nvertices;n++){
		Vertex *vdst = &surface->vertices[n];
		VERTEX *vsrc = &mris->vertices[n];
		vdst->x=vsrc->x;
		vdst->y=vsrc->y;
		vdst->z=vsrc->z;
	}
	for(int n = 0 ; n < mris->nfaces ;n++){
		Face *fdst = &surface->faces[n];
		FACE *fsrc = &mris->faces[n];
		fdst->v[0]=fsrc->v[0];
		fdst->v[1]=fsrc->v[1];
		fdst->v[2]=fsrc->v[2];
	}
	
	surface->InitSurface();

	return surface;
}

void SurfaceToMRIP(Surface &surface, MRIS &mris){
	mris.nvertices=surface.nvertices;
	mris.nfaces=surface.nfaces;
	//Vertices
	for(int n = 0 ; n < surface.nvertices ; n++){
		VERTEX *vdst=&mris.vertices[n];
		Vertex *vsrc = &surface.vertices[n];
		//vertices
		if(vdst->v) delete [] vdst->v;
		vdst->v = (int*)calloc(vsrc->vnum , sizeof(int));
		for(int p = 0 ; p < vsrc->vnum ; p++)
			vdst->v[p]=vsrc->v[p];
		vdst->vnum=vsrc->vnum;
		//faces
		if(vdst->f) delete [] vdst->f;
		if(vdst->n) delete [] vdst->n;
		vdst->f = (int*)calloc(vsrc->fnum , sizeof(int));
		vdst->n = (uchar*)calloc(vsrc->fnum , sizeof(uchar));
		for(int p = 0 ; p < vsrc->fnum ; p++){
			vdst->f[p]=vsrc->f[p];
			vdst->n[p]=vsrc->n[p];
		}
		vdst->num=vsrc->fnum;
	}
	//Faces
	for(int n = 0 ; n < surface.nfaces ; n++){
		FACE *fdst=&mris.faces[n];
		Face *fsrc = &surface.faces[n];
		fdst->v[0]=fsrc->v[0];
		fdst->v[1]=fsrc->v[1];
		fdst->v[2]=fsrc->v[2];
	}
}


bool MRISaddMRIP(MRIS *mris_dst, MRIP *mrip){
	int nvertices=mrip->n_vertices,nfaces=mrip->n_faces;
	int *vto = mrip->vtrans_to,*fto=mrip->ftrans_to;
	
	MRIS *mris = mrip->mris;

	//vertices
	for(int n = 0 ; n < mris->nvertices ; n++){
		VERTEX *vsrc = &mris->vertices[n];
		if(n >= nvertices) vto[n] = mris_dst->nvertices++;
		VERTEX *vdst = &mris_dst->vertices[vto[n]];
		//coordonnees
		vdst->x=vsrc->x;
		vdst->y=vsrc->y;
		vdst->z=vsrc->z;
	}
	//faces
	for(int n = 0 ; n < mris->nfaces ; n++){
		FACE *fsrc = &mris->faces[n];
		if(n >= nfaces) fto[n] = mris_dst->nfaces++;
		FACE *fdst = &mris_dst->faces[vto[n]];
		//vertex indices
		fdst->v[0]=vto[fsrc->v[0]];
		fdst->v[1]=vto[fsrc->v[1]];
		fdst->v[2]=vto[fsrc->v[2]];
	}
	//stuff in vertices and faces
	for(int n = 0 ; n < mris->nvertices ; n++){
		VERTEX *vsrc = &mris->vertices[n];
		VERTEX *vdst = &mris_dst->vertices[vto[n]];
		vsrc=vdst=NULL;
	//...

	}


	return true;
}
