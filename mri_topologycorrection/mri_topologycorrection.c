#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "mrisurf.h"
#include "matrix.h"
#include "proto.h" 
#include "stats.h"
#include "timer.h"
#include "const.h"
#include "mrishash.h"
#include "icosahedron.h"
#include "tritri.h"
#include "timer.h"
#include "chklc.h"
#include "diag.h"
#include "mri_topology.h"
#include "mri_tess.h"
#include "mrisutils.h"
#include "gca.h"
#include "transform.h"

#define SQR(x) ((x)*(x))

char *Progname;

typedef struct Cell
{
  struct Cell *next,*previous;
  int x,y,z;
  unsigned char type;
  float prior;
  float cost;
}Cell;


//for list(cell)
#define BIN_NBR 10000

#define NONSIMPLE 1
#define SIMPLE 2
#define UNKNOWN 3
#define KNOWN 4

#define MAXCELLS 10000
#define MAXCOMPONENTS 1000
#define MAXCOST 20.0f

#define NBR_OF_COMPTS 500

//for initSegmentation (single_mode) and CorrectSegmentation (multiple_mode)
#define SINGLE_MODE 0
#define MULTIPLE_MODE 1

//flags for the segmentation
#define BODY 1
#define VISITED 3 
#define RESIDUE 2
#define ISOLATED 4

#define F_B 11
#define F_R 12
#define F_I 14

#define B_B 21
#define B_R 22
#define B_I 24


#define DINTMAX 10
#define DINTMIN 2
#define DOUTMAX 10
#define DOUTMIN 2
#define PDMAX 1.0f

//threshold for correcting the volume
#define INITIAL_THRESHOLD 0.0
#define THRESHOLD_INCREASE 0.05

//different mode of corrections
#define NORMAL_MODE 0
#define VOXEL_MODE 1 
#define PROB_MODE 2
#define PROB_MAP_MODE 3
#define MAP_MODE 4

//define the minimum probability (zero doesn't exist!)
#define MINPROB 0.0001f
//define the maximum probability (one doesn't exist!)
#define MAXPROB 0.999f

//define the minmax value
#define MINMAX(a,b,c) MIN(a,MAX(b,c))
#define PROB(x) MINMAX(MAXPROB,MINPROB,x)

typedef struct
{
  int number_of_bins;
  int reference_max_cell;
  int ncells;
  Cell** sorting;

  int width,height,depth;
  Cell*** table;
}List;

static List *initList(int nbr_of_bins,int width,int height,int depth);
static int listFree(List **list);
static int addCell(List *list,Cell *new_cell);
static int removeCell(List *list, Cell *old_cell);
static float f(float x);

typedef struct ConnectedComponent
{
  struct ConnectedComponent *next,*previous;
  int ncells;
  int maxcells;
  Cell **cells;
  float cost;
  float map;
  unsigned char found;
}ConnectedComponent;

typedef struct SEGMENTATION
{
  int ncomponents;
  int maxcomponents;
  ConnectedComponent *components;
  void* ccs;
}SEGMENTATION;

typedef struct CCSorting
{
  int number_of_bins;
  int reference_max_component;
  int ncomponents;
  float maxcost;
  ConnectedComponent** sorting;

  SEGMENTATION *segmentation;
}CCS;
 
static CCS *initCCS(int nbr_of_components,SEGMENTATION *segmentation);
static int CCSfree(CCS **ccs);
static int ccsSortComponents(CCS* ccs);
static int removeComponent(CCS *ccs,ConnectedComponent *old_cpmt);
static int addComponent(CCS *ccs,ConnectedComponent *new_component);
static int reallocateComponent(ConnectedComponent *cc,int maxcells);
static int reallocateSegmentation(SEGMENTATION *segmentation, int maxcomponents);
static int componentNew(SEGMENTATION *segmentation);
static SEGMENTATION* segmentationAlloc(int maxcomponents,int maxcells);
static int segmentationFree(SEGMENTATION **segmentation);
static int addCellToComponent(ConnectedComponent *cc,Cell *cell);


/*parameters used by the process*/
typedef struct TC_PARMS    
{ 
  /*volume with orig values*/
  MRI *mri_orig;
  /*volume with labels*/
  MRI *mri_seg;
  /*output volume*/
  MRI *mri_output;
  /*binary volume*/
  MRI *mri_bin;
  /*distance map*/
  MRI *mri_dist;
  /*probability map*/
  MRI *mri_prob;
  MRI *mri_cprob;
  /*cost map for foreground and background*/
  MRI *mri_fcost;
  MRI *mri_bcost;
  /*prior map for foreground and background*/
  MRI *mri_fprior;
  MRI *mri_bprior;
  /*label map*/
  MRI *mri_labeled;
  //max & min distance inside and outside
  float maxindist,maxoutdist;
  /*weighting parameters between dist map and prior map for prioritization and prior*/
  float alpha,beta;
  
  /*label that we want to segment*/
  int nblabels;
  int labels[50];
  /*connectivity*/
  int c_c; //current
  int f_c; //foreground
  int b_c; //background

  /*final surface output name*/
  char *surfname;
  /*writing out the different maps*/
  char *mapsfname;
  /*current surface*/
  MRIS *mris;
  /*compute the initial param of the surface or not*/
  int  initsurface;
  
  /* to fit the surface onto the volume*/
  int fit;
  /*to use a specific tesselation*/
  int tess;
  /*compute the cost from the number of voxels*/
  int costmode;

  int guess;

  char *priormap;

  /*voxel coord for mri_orig,mri_seg*/   
  MRI_REGION  region;
  /*dimensions for mri_bin, mri_dist, ...*/
  int width,height,depth;

  List *list;
  int current_label;
  SEGMENTATION *segmentation;
  CCS *ccs;

  //mode used (single region growing or multiple merging)
  int multiplemode;
  
  //to specify foreground or background
  int mode;
  //specify only mode
  int only;
  //float threshold
  float threshold;

  int border_labels[27];
  int nlabels;

  //the four segmentation and sorting (F_B,F_R,B_B,B_R)
  SEGMENTATION *F_Bseg,*F_Rseg,*B_Bseg,*B_Rseg;
  CCS *F_Bccs,*F_Rccs,*B_Bccs,*B_Rccs;

  //using priors
  int priors;
  char* transform_fname;
  char *gca_fname;

} TC_PARMS ;

static TC_PARMS *initTC_PARMS(void);
static int TC_PARMSfree(TC_PARMS **parms);

static void Error(char *string);
static int get_option(int argc, char *argv[],TC_PARMS *parms);
static void initImages(TC_PARMS* parms);
static int initSegmentation(TC_PARMS *parms);
static int mriChangeLabel(MRI *mri,int src,int dst);
static int segmentBody(TC_PARMS *parms);
static void CTExpansion(TC_PARMS *parms);
static int initCCSSEG(TC_PARMS *parms);
static int initCellsFromMap(TC_PARMS *parms);
static int computeResidualSegmentation(TC_PARMS *parms);
static int componentMerge(TC_PARMS *parms,int s0 ,int s1);
static void SaveOrigMaps(TC_PARMS *parms);

static List *initList(int nbr_of_bins,int width,int height,int depth)
{
  int k,i,j;
  List *list;
  if(nbr_of_bins<1)
    nbr_of_bins=BIN_NBR;

  list=(List*)calloc(1,sizeof(List));
  
  list->width=width;
  list->height=height;
  list->depth=depth;
  list->table=(Cell***)malloc(depth*sizeof(Cell**));
  for(k=0;k<depth;k++)
    {
      list->table[k]=(Cell**)malloc(height*sizeof(Cell*));
      for(j=0;j<height;j++)
	{
	  list->table[k][j]=(Cell*)calloc(width,sizeof(Cell));
	  for(i=0;i<width;i++)
	    {
	      list->table[k][j][i].prior=-1;
	      list->table[k][j][i].cost=0;
	      list->table[k][j][i].type=UNKNOWN;
	      list->table[k][j][i].x=i;
	      list->table[k][j][i].y=j;
	      list->table[k][j][i].z=k;
	    }
	}
    }
  
  list->number_of_bins=nbr_of_bins;
  list->sorting=(Cell**)calloc(nbr_of_bins,sizeof(Cell*));
  
  return list;
}

static int listFree(List **list)
{
  int k,j,height,depth;
  List *l=*list;
  *list=NULL;
  
  height=l->height;
  depth=l->depth;

  free(l->sorting);
  for(k=0;k<depth;k++)
    {
      for(j=0;j<height;j++)
	free(l->table[k][j]);
      free(l->table[k]);
    }
  free(l->table);
  free(l);

  return NO_ERROR;
}

static int removeCell(List *list, Cell *old_cell)
{
  int ref;
  Cell *ncell;
  float old_prior=old_cell->prior;

  ref=(int)(old_prior*list->number_of_bins);
  if(ref==list->number_of_bins)  //if prior==1!
    ref--;

  ncell=old_cell->next;
  if(old_cell->previous==NULL)  //first element of the list
    {
      list->sorting[ref]=ncell;
      if(ncell)
	ncell->previous=NULL;
    }
  else               //non first element of the list
    {
      old_cell->previous->next=ncell;
      if(ncell)
	ncell->previous=old_cell->previous;
    }
  
  
  if(ref==list->reference_max_cell && (!list->sorting[ref]))
    {
      int k;
      for(k=ref-1;k>=0;k--)
	if(list->sorting[k])
	  break;
      list->reference_max_cell=k; //eventually -1 if last point!
    }

  old_cell->previous=NULL;
  old_cell->next=NULL;

  list->ncells--;

  return NO_ERROR;
}

static int addCell(List *list,Cell *new_cell)
{
  Cell *cell,*pcell;
  int ref;
  float new_prior=new_cell->prior;

  ref=(int)(new_prior*list->number_of_bins);
  if(ref==list->number_of_bins)  //if prior==1!
    ref--;

  cell=list->sorting[ref];
  if(!cell)  //new_cell will be the first and unique element of the list
    {
      new_cell->next=NULL;
      new_cell->previous=NULL;
      list->sorting[ref]=new_cell;
    }
  else  //find the correct location
    {
      if(new_prior>cell->prior)  //must be the first element
	{
	  new_cell->next=cell;
	  new_cell->previous=NULL;
	  list->sorting[ref]=new_cell;
	  cell->previous=new_cell;
	}
      else    //must find the correct location 
	{     //in this case, we put the cell at the last correct location
	  pcell=cell;cell=cell->next;
	  while((cell!=NULL) && (cell->prior>=new_prior))
	    {
	      pcell=cell;
	      cell=cell->next;
	}
	  new_cell->next=cell;
	  new_cell->previous=pcell;
	  pcell->next=new_cell;
	  if(cell)
	    cell->previous=new_cell;
	}
    }
  
  
  if(ref>list->reference_max_cell)
    list->reference_max_cell=ref;

  list->ncells++;

  return NO_ERROR;
}

#if 0
static int sortingList(List *list)
{
  int k,j,i,width,height,depth;
  Cell *cell;
  Cell*** table=list->table;
 
  width=list->width;
  height=list->height;
  depth=list->depth;
  
 //reinit
  memset(list->sorting,0,list->number_of_bins*sizeof(Cell*));
  list->reference_max_cell=0;
  list->ncells=0;
  //sorting
  for(k=0;k<depth;k++)
    for(j=0;j<height;j++)
      for(i=0;i<width;i++)
	{
	  cell=&table[k][j][i];
	  if(cell->prior>=0) //the prior is never negative
	    addCell(list,cell);
	}

  return NO_ERROR;
}
#endif




static CCS *initCCS(int nbr_of_bins,SEGMENTATION *segmentation)
{
  CCS *ccs;
  if(nbr_of_bins<1)
    nbr_of_bins=NBR_OF_COMPTS;

  ccs=(CCS*)calloc(1,sizeof(CCS));
  
  ccs->number_of_bins=nbr_of_bins;
  ccs->segmentation=segmentation;
  ccs->sorting=(ConnectedComponent**)calloc(nbr_of_bins,sizeof(ConnectedComponent*));
  ccs->maxcost=MAXCOST;
  
  return ccs;
}

static int CCSfree(CCS **ccs)
{
  CCS *c=*ccs;
  *ccs=NULL;
  
  free(c->sorting);
  
  free(c);  
  return NO_ERROR;
}

static int ccsSortComponents(CCS* ccs)
{
  int k;
  ConnectedComponent *cmpt;
  
  //reinit
  memset(ccs->sorting,0,ccs->number_of_bins*sizeof(ConnectedComponent*));
  ccs->ncomponents=0;
  ccs->reference_max_component=0;
  //sort the different segments
  for(k=0;k<ccs->segmentation->maxcomponents;k++)
    {
      cmpt=&ccs->segmentation->components[k];
      if(cmpt->cost>=0)
	addComponent(ccs,cmpt);
    }

  return NO_ERROR;
}

static int removeComponent(CCS *ccs,ConnectedComponent *old_cmpt)
{
  int ref;
  ConnectedComponent *ncmpt;
  float old_cost=old_cmpt->cost;

  ref=(int)(old_cost*ccs->number_of_bins/ccs->maxcost);
  if(ref>=ccs->number_of_bins)  //if cost==1!
    ref=ccs->number_of_bins-1;

  ncmpt=old_cmpt->next;
  if(old_cmpt->previous==NULL)  //first element of the list
    {
      ccs->sorting[ref]=ncmpt;
      if(ncmpt)
	ncmpt->previous=NULL;
    }
  else               //non first element of the list
    {
      old_cmpt->previous->next=ncmpt;
      if(ncmpt)
	ncmpt->previous=old_cmpt->previous;
    }
  
  
  if(ref==ccs->reference_max_component && (!ccs->sorting[ref]))
    {
      int k;
      for(k=ref-1;k>=0;k--)
	if(ccs->sorting[k])
	  break;
      ccs->reference_max_component=k; //eventually -1 if last point!
    }

  old_cmpt->previous=NULL;
  old_cmpt->next=NULL;
  
  ccs->ncomponents--;

  return NO_ERROR;
}

static int addComponent(CCS *ccs,ConnectedComponent *new_component)
{
  ConnectedComponent *component,*pcomponent;
  int ref;
  float new_cost=new_component->cost;

  ref=(int)(new_cost*ccs->number_of_bins/ccs->maxcost);
  if(ref>=ccs->number_of_bins)  //if cost==1!
    ref=ccs->number_of_bins-1;

  component=ccs->sorting[ref];
  if(!component)  //new_component will be the first and unique element of the ccs
    {
      new_component->next=NULL;
      new_component->previous=NULL;
      ccs->sorting[ref]=new_component;
    }
  else  //find the correct location
    {
      if(new_cost>component->cost)  //must be the first element
	{
	  new_component->next=component;
	  new_component->previous=NULL;
	  ccs->sorting[ref]=new_component;
	  component->previous=new_component;
	}
      else    //must find the correct location 
	{     //in this case, we put the component at the last correct location
	  pcomponent=component;component=component->next;
	  while((component!=NULL) && (component->cost>=new_cost))
	    {
	      pcomponent=component;
	      component=component->next;
	}
	  new_component->next=component;
	  new_component->previous=pcomponent;
	  pcomponent->next=new_component;
	  if(component)
	    component->previous=new_component;
	}
    }
  
  
  if(ref>ccs->reference_max_component)
    ccs->reference_max_component=ref;

  ccs->ncomponents++;

  return NO_ERROR;
}


static int reallocateComponent(ConnectedComponent *cc,int maxcells)
{
  int k;
  Cell **oldcells;

  if(maxcells<=0)
    maxcells=MAXCELLS;

  oldcells=cc->cells;

  cc->cells=(Cell**)calloc(maxcells,sizeof(Cell*));
  cc->maxcells=maxcells;
  for(k=0;k<cc->ncells;k++)
    cc->cells[k]=oldcells[k];
    
  if(oldcells)
    free(oldcells);
  return NO_ERROR;
}

static int reallocateSegmentation(SEGMENTATION *segmentation, int maxcomponents)
{
  int k,n;
  CCS* ccs;
  ConnectedComponent *oldcc;

  oldcc=segmentation->components;
  
  segmentation->components=(ConnectedComponent*)calloc(maxcomponents,sizeof(ConnectedComponent));
  for(k=0;k<segmentation->maxcomponents;k++)
    {
      segmentation->components[k].ncells=oldcc[k].ncells;
      segmentation->components[k].maxcells=oldcc[k].maxcells;
      for(n=0;n<oldcc[k].ncells;n++)
	segmentation->components[k].cells[n]=oldcc[k].cells[n];
      segmentation->components[k].found=oldcc[k].found;
    }
  segmentation->maxcomponents=maxcomponents;
  
  
  ccs=(CCS*)segmentation->ccs;
  if(ccs)
    ccsSortComponents(ccs);

  if(oldcc)
    free(oldcc);

  return NO_ERROR;
}

static int addCellToComponent(ConnectedComponent *cc,Cell *cell)
{
  if(cc->ncells>=cc->maxcells)
    reallocateComponent(cc,nint(cc->maxcells*1.5));
  
  cc->cells[cc->ncells]=cell;
  cc->ncells++;
  if(cc->ncells==1)
    cc->map=0;
  cc->map+=f(cell->cost);
  cc->cost=MAXCOST*(1.0f+tanh(cc->map))/2.;
  return NO_ERROR;
}

static int componentNew(SEGMENTATION *segmentation)
{
  int s;
  if(segmentation->ncomponents>=segmentation->maxcomponents)
    reallocateSegmentation(segmentation,nint(segmentation->maxcomponents*1.5));
  segmentation->ncomponents++;
   
  for(s=0;s<segmentation->maxcomponents;s++)
    if(segmentation->components[s].ncells == 0)
      return (s);
  return (s);
}

static SEGMENTATION* segmentationAlloc(int maxcomponents,int maxcells)
{
  SEGMENTATION *segmentation;
  ConnectedComponent *cc;
  int k;

  segmentation=(SEGMENTATION*)calloc(1,sizeof(SEGMENTATION));
  segmentation->maxcomponents=maxcomponents;
  segmentation->ncomponents=0;
  segmentation->components=(ConnectedComponent*)calloc(maxcomponents,sizeof(ConnectedComponent));
  for(k=0;k<maxcomponents;k++)
    {
      cc=&segmentation->components[k];
      cc->maxcells=maxcells;
      cc->cells=(Cell**)calloc(maxcells,sizeof(Cell*));
      cc->cost=-1.0;
    }

  return segmentation;
}

static int segmentationFree(SEGMENTATION **segmentation)
{
  int k;
  SEGMENTATION *seg=*segmentation;

  if(!seg)
    return NO_ERROR;
  
  *segmentation=NULL;

  for(k=0;k<seg->maxcomponents;k++)
    if(seg->components[k].maxcells)
      free(seg->components[k].cells);
  free(seg->components);
  free(seg);

  return NO_ERROR;
}

static TC_PARMS *initTC_PARMS(void)
{
  TC_PARMS *parms;

  parms=(TC_PARMS*)calloc(1,sizeof(TC_PARMS));
  parms->f_c=4;
  parms->b_c=associatedConnectivity(parms->f_c);
  parms->c_c=parms->f_c;
  parms->alpha=1.0f;
  parms->beta=1.0f;
  parms->costmode=VOXEL_MODE;
  parms->tess=-1;

  parms->threshold=INITIAL_THRESHOLD;
  
  return parms;
}

static int initCCSSEG(TC_PARMS *parms)
{
  parms->segmentation=segmentationAlloc(MAXCOMPONENTS,MAXCELLS);
  parms->ccs=initCCS(NBR_OF_COMPTS,parms->segmentation);
  parms->segmentation->ccs=parms->ccs;
  parms->ccs->segmentation=parms->segmentation;
  return NO_ERROR;
}

static int TC_PARMSfree(TC_PARMS **parms)
{
  TC_PARMS *p=*parms;
  *parms=NULL;

  if(p->list)
    listFree(&p->list);

  if(p->segmentation)
    segmentationFree(&p->segmentation);

  if(p->ccs)
    CCSfree(&p->ccs);

  free(p);

  return NO_ERROR;
}


//--------------------------------------------
/*Error routine*/
static void Error(char *string)
{
  fprintf(stderr, "\nError %s\n",string) ;
  exit(0) ;
}

/*get_option routine*/
static int get_option(int argc, char *argv[],TC_PARMS *parms)
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  
  if (!strcmp(option, "nothing"))
  {
    fprintf(stderr,"Mode:          NOTHING\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "surf"))
  {
    parms->surfname=argv[2];
    fprintf(stderr,"Mode:          Writing out final surface into %s\n",parms->surfname) ;
    nargs = 1 ;
  }else if (!strcmp(option, "foreground"))
  {
    parms->mode=1;
    fprintf(stderr,"Mode:          Making corrections on the background first\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "fit"))
  {
    parms->fit=1;
    fprintf(stderr,"Mode:          Fitting surface on volume\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "tess"))
  {
    parms->tess=atoi(argv[2]);
    fprintf(stderr,"Mode:          Tesselating with mode %d\n",parms->tess) ;
    nargs = 1 ;
  }else if (!strcmp(option, "priormap"))
  {
    parms->priormap=argv[2];
    fprintf(stderr,"Mode:          Loading prior maps from %s\n",parms->priormap) ;
    nargs = 1 ;
  }else if (!strcmp(option, "VOXEL"))
  {
    parms->costmode=VOXEL_MODE;
    fprintf(stderr,"Mode:          Cost computed from the number of voxels\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "MAP"))
  {
    parms->costmode=MAP_MODE;
    fprintf(stderr,"Mode:          Cost computed from map estimate\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "PROB"))
  {
    parms->costmode=PROB_MODE;
    fprintf(stderr,"Mode:          Cost computed from the sum of the voxel prob\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "PROB_MAP"))
  {
    parms->costmode=PROB_MAP_MODE;
    fprintf(stderr,"Mode:          Cost computed from map estimate\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "only"))
  {
    parms->only=1;
    fprintf(stderr,"Mode:          Making corrections on the foreground/background only\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "guess"))
  {
    parms->guess=1;
    fprintf(stderr,"Mode:          Guess the initial Segmentation\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "background"))
  {
    parms->mode=2;
    fprintf(stderr,"Mode:          Making corrections on the foreground first\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "priors"))
  {
    parms->priors=1;
    parms->transform_fname=argv[2];
    parms->gca_fname=argv[3];
    fprintf(stderr,"Mode:          Using prior information: "
	    "\n               transform %s"
	    "\n               gca %s\n"
	    ,parms->transform_fname,parms->gca_fname) ;
    nargs = 2 ;
  }else if (!strcmp(option, "maps"))
  {
    parms->mapsfname=argv[2];
    fprintf(stderr,"Mode:          Writing out the maps into folder %s\n",parms->mapsfname) ;
    nargs = 1 ;
  }else if (!strcmp(option, "init"))
  {
    parms->initsurface=1;
    fprintf(stderr,"Mode:          Computing the initial surface parameters\n") ;
    nargs = 0 ;
  }else if (!strcmp(option, "connectivity"))
  {
    parms->f_c=atoi(argv[2]);
    parms->b_c=associatedConnectivity(parms->f_c);
    fprintf(stderr,"Mode:          Connectivity %d\n",parms->f_c) ;
    nargs = 1 ;
  }else if (!strcmp(option, "label"))
  {
    parms->labels[parms->nblabels]=atoi(argv[2]);
    parms->nblabels++;
    fprintf(stderr,"Mode:          Correction topology for label %d\n",atoi(argv[2])) ;
    nargs = 1 ;
  }else if (!strcmp(option, "beta"))
  {
    parms->beta=MAX(0.0,MIN(1.0,atof(argv[2])));
    fprintf(stderr,"Mode:          Mixing parameters set to %3.3f\n",parms->beta) ;
    nargs = 1 ;
  }else if (!strcmp(option, "alpha"))
  {
    parms->alpha=MAX(0.0,MIN(1.0,atof(argv[2])));
    fprintf(stderr,"Mode:          Mixing parameters set to %3.3f\n",parms->alpha) ;
    nargs = 1 ;
  }else switch (toupper(*option))
  {
  case 'L':
    parms->labels[parms->nblabels]=atoi(argv[2]);
    parms->nblabels++;
    fprintf(stderr,"Mode:          Correction topology for label %d\n",atoi(argv[2])) ;
    nargs = 1 ;
    break;
  default:
    printf("Mode:          unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
} 



static void setBorderValueToSegmentation(TC_PARMS *parms,int s,int dst)
{
  SEGMENTATION *segmentation=parms->segmentation;
  Cell ***table=parms->list->table;
  int i,j,k;
  int width,height,depth;
  
  width=parms->mri_bin->width;
  height=parms->mri_bin->height;
  depth=parms->mri_bin->depth;

  k=dst;
  for(i=dst;i<width-dst;i++)
    for(j=dst;j<height-dst;j++)
      addCellToComponent(&segmentation->components[s],&table[k][j][i]);
      
  k=depth-1-dst;
  for(i=dst;i<width-dst;i++)
    for(j=dst;j<height-dst;j++)
      addCellToComponent(&segmentation->components[s],&table[k][j][i]);
      
  j=dst;
  for(i=dst;i<width-dst;i++)
    for(k=dst;k<depth-dst;k++)
      addCellToComponent(&segmentation->components[s],&table[k][j][i]);
      
  j=height-1-dst;
  for(i=dst;i<width-dst;i++)
    for(k=dst;k<depth-dst;k++)
      addCellToComponent(&segmentation->components[s],&table[k][j][i]);
      
  i=dst;
  for(k=dst;k<depth-dst;k++)
    for(j=dst;j<height-dst;j++)
      addCellToComponent(&segmentation->components[s],&table[k][j][i]);
      
  i=width-1-dst;
  for(k=dst;k<depth-dst;k++)
    for(j=dst;j<height-dst;j++)
      addCellToComponent(&segmentation->components[s],&table[k][j][i]);
      
}

static void setBorderValue(MRI *mri,int val,int dst)
{
  int i,j,k;
  int width,height,depth;
  
  width=mri->width;
  height=mri->height;
  depth=mri->depth;

  switch(mri->type)
    {
    case MRI_UCHAR:
      k=dst;
      for(i=dst;i<width-dst;i++)
	for(j=dst;j<height-dst;j++)
	  MRIvox(mri,i,j,k)=val;
      k=depth-1-dst;
      for(i=dst;i<width-dst;i++)
	for(j=dst;j<height-dst;j++)
	  MRIvox(mri,i,j,k)=val;
      j=dst;
      for(i=dst;i<width-dst;i++)
	for(k=dst;k<depth-dst;k++)
	  MRIvox(mri,i,j,k)=val;
      j=height-1-dst;
      for(i=dst;i<width-dst;i++)
	for(k=dst;k<depth-dst;k++)
	  MRIvox(mri,i,j,k)=val;
      i=dst;
      for(k=dst;k<depth-dst;k++)
	for(j=dst;j<height-dst;j++)
	  MRIvox(mri,i,j,k)=val;
      i=width-1-dst;
      for(k=dst;k<depth-dst;k++)
	for(j=dst;j<height-dst;j++)
	  MRIvox(mri,i,j,k)=val;
      break;
    case MRI_SHORT:
      k=dst;
      for(i=dst;i<mri->width-dst;i++)
	for(j=dst;j<height-dst;j++)
	  MRISvox(mri,i,j,k)=val;
      k=depth-1-dst;
      for(i=dst;i<width-dst;i++)
	for(j=dst;j<height-dst;j++)
	  MRISvox(mri,i,j,k)=val;
      j=dst;
      for(i=dst;i<width-dst;i++)
	for(k=dst;k<depth-dst;k++)
	  MRISvox(mri,i,j,k)=val;
      j=height-1-dst;
      for(i=dst;i<width-dst;i++)
	for(k=dst;k<depth-dst;k++)
	  MRISvox(mri,i,j,k)=val;
      i=dst;
      for(k=dst;k<depth-dst;k++)
	for(j=dst;j<height-dst;j++)
	  MRISvox(mri,i,j,k)=val;
      i=width-1-dst;
      for(k=dst;k<depth-dst;k++)
	for(j=dst;j<height-dst;j++)
	  MRISvox(mri,i,j,k)=val;
      break;
    }
}

static float f(float x)
{
  return x;
} 

static int initCellsFromMap(TC_PARMS *parms)
{
  int i,j,k;
  int width,height,depth;
  float cost,dval,val,pval,A,prior;
  char fname[100];
  MRI *mri_fprior=parms->mri_fprior
    ,*mri_bprior=parms->mri_bprior,*mri_dist=parms->mri_dist
    ,*mri_fcost=parms->mri_fcost,*mri_bcost=parms->mri_bcost;
  Cell ***table=parms->list->table; 
  MRI *mri_c,*mri_prior;

  width=mri_dist->width;
  height=mri_dist->height;
  depth=mri_dist->depth;
  

  A=parms->alpha;

  switch(parms->costmode)
    {
    case VOXEL_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      dval=fabs(MRIFvox(mri_dist,i,j,k));
	      prior=dval;
	      cost=1.0f;
	      table[k][j][i].cost=cost;
	      table[k][j][i].prior=prior;
	    }
      break;
    case PROB_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      val=MRIFvox(mri_dist,i,j,k);
	      dval=fabs(val);
	      if(val>0) //inside
		{  
		  pval=MRIFvox(mri_fprior,i,j,k);
		  cost=MRIFvox(mri_fcost,i,j,k);
		}
	      else
		{
		  pval=MRIFvox(mri_bprior,i,j,k);
		  cost=MRIFvox(mri_bcost,i,j,k);
		}
	      if(A!=1)
		dval=0.5*(1.+fabs(val));
	      
	      prior=A*dval+(1.0f-A)*pval;

	      if(A==0)
		{
		  if(prior==MAXPROB)
		    prior=prior*dval;
		  else
		    prior=prior/2.;
		}

	      table[k][j][i].prior=prior;   
	      table[k][j][i].cost=cost;
	    }
      break;
    case MAP_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      val=MRIFvox(mri_dist,i,j,k);
	      dval=0.5*(1.+fabs(val));
	      if(val>0) //inside
		{  
		  prior=MRIFvox(mri_fprior,i,j,k);
		  cost=MRIFvox(mri_fcost,i,j,k);
		}
	      else
		{
		  prior=MRIFvox(mri_bprior,i,j,k);
		  cost=MRIFvox(mri_bcost,i,j,k);
		}
	      if(prior>=(0.5+0.5*tanh(0.1*log(MAXPROB/MINPROB))))
		  prior=prior*dval;
	      else
		prior=prior/2.;

	      table[k][j][i].prior=prior;
	      table[k][j][i].cost=cost;
	    }
      break;
   case PROB_MAP_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      val=MRIFvox(mri_dist,i,j,k);
	      dval=fabs(val);
	      if(val>0) //inside
		{  
		  pval=MRIFvox(mri_fprior,i,j,k);
		  cost=MRIFvox(mri_fcost,i,j,k);
		}
	      else
		{
		  pval=MRIFvox(mri_bprior,i,j,k);
		  cost=MRIFvox(mri_bcost,i,j,k);
		}
	      if(A!=1)
		dval=0.5*(1.+fabs(val));      
	      prior=A*dval+(1.0f-A)*pval;
	      if(A==0)
		{
		  if(prior==MAXPROB)
		    prior=prior*dval;
		  else
		    prior=prior/2.;
		}
	      table[k][j][i].prior=prior;   
	      table[k][j][i].cost=cost;
	    }
      break;
    default:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      val=MRIFvox(mri_dist,i,j,k);
	      dval=fabs(val);
	      cost=0.5*(1.+dval);
	      prior=dval;
	      table[k][j][i].prior=prior;   
	      table[k][j][i].cost=cost;
	    }
      break;
    }

  if(parms->mapsfname)
    {
      mri_c=MRIalloc(width,height,depth,MRI_FLOAT);
      mri_prior=MRIalloc(width,height,depth,MRI_FLOAT);

      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_c,i,j,k)=table[k][j][i].cost;
	      MRIFvox(mri_prior,i,j,k)=table[k][j][i].prior;
	    }
      

      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_d.mgh");
      MRIwrite(mri_dist,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_c.mgh");
      MRIwrite(mri_c,fname);
      MRIfree(&mri_c);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_pr.mgh");
      MRIwrite(mri_prior,fname);
      MRIfree(&mri_prior);
    }

  return NO_ERROR;
}

static void addBorderVoxels(TC_PARMS *parms,int dst)
{
  int i,j,k;
  MRI *mri=parms->mri_bin;
  int width,height,depth;
  List *list=parms->list;
  Cell *cell; 
  Cell *** table=parms->list->table;
 
  width=mri->width;
  height=mri->height;
  depth=mri->depth;
  

  k=dst;
  for(i=dst;i<(width-dst);i++)
    for(j=dst;j<(height-dst);j++)
      if(MRIvox(mri,i,j,k)==RESIDUE)
      {
	cell=&table[k][j][i];
	cell->type=UNKNOWN;
	addCell(list,cell);
	MRIvox(mri,i,j,k)=VISITED;    
      }
  k=depth-1-dst;
  for(i=dst;i<(width-dst);i++)
    for(j=dst;j<(height-dst);j++)
      if(MRIvox(mri,i,j,k)==RESIDUE)
      {
	cell=&table[k][j][i];
	cell->type=UNKNOWN;
	addCell(parms->list,cell);
	MRIvox(mri,i,j,k)=VISITED;  
      }
  j=dst;
  for(i=dst;i<(width-dst);i++)
    for(k=dst+1;k<(depth-dst-1);k++)
      if(MRIvox(mri,i,j,k)==RESIDUE)
      {
	cell=&table[k][j][i];
	cell->type=UNKNOWN;
	addCell(parms->list,cell);
	MRIvox(mri,i,j,k)=VISITED;  
      }
  j=height-1-dst;
  for(i=dst;i<(width-dst);i++)
    for(k=dst+1;k<(depth-dst-1);k++)
      if(MRIvox(mri,i,j,k)==RESIDUE)
      {
	cell=&table[k][j][i];
	cell->type=UNKNOWN;
	addCell(parms->list,cell);
	MRIvox(mri,i,j,k)=VISITED;  
      }
  i=dst;
  for(k=dst+1;k<(depth-dst-1);k++)
    for(j=dst+1;j<(height-dst-1);j++)
       if(MRIvox(mri,i,j,k)==RESIDUE)
       {
	cell=&table[k][j][i];
	cell->type=UNKNOWN;
	addCell(parms->list,cell);
	MRIvox(mri,i,j,k)=VISITED;  
      }
  i=width-1-dst;
  for(k=dst+1;k<(depth-dst-1);k++)
    for(j=dst+1;j<(height-dst-1);j++)
      if(MRIvox(mri,i,j,k)==RESIDUE)
      {
	cell=&table[k][j][i];
	cell->type=UNKNOWN;
	addCell(parms->list,cell);
	MRIvox(mri,i,j,k)=VISITED;  
      }

}

static float computeMap(TC_PARMS *parms)
{
  int width,height,depth,i,j,k,x,y,z,xinit,yinit,zinit,val,n,nlabels,*tab;
  MRI *mri_seg=parms->mri_seg,*mri=parms->mri_bin,*mri_psi=parms->mri_prob,*mri_pcsi=parms->mri_cprob;
  double map;

  width=mri->width;
  height=mri->height;
  depth=mri->depth;

  xinit=parms->region.x;
  yinit=parms->region.y;
  zinit=parms->region.z;

  nlabels=parms->nblabels;
  tab=parms->labels;

  map=0;
  for(k=2;k<depth-2;k++)
    for(j=2;j<height-2;j++)
      for(i=2;i<width-2;i++)
	{
	  x=xinit+i-2;
	  y=yinit+j-2;
	  z=zinit+k-2;
	  val=1;
	  for(n=0;n<nlabels;n++)
	    if(MRIvox(mri_seg,x,y,z)==tab[n])
	      val=0;
	  if(val && (MRIvox(mri,i,j,k)==F_B||MRIvox(mri,i,j,k)==F_R))
	     map+=log(MRIFvox(mri_psi,i,j,k)/MRIFvox(mri_pcsi,i,j,k));
	  if(!val && (MRIvox(mri,i,j,k)==B_B||MRIvox(mri,i,j,k)==B_R))
	     map-=log(MRIFvox(mri_psi,i,j,k)/MRIFvox(mri_pcsi,i,j,k));
	}
  return map;
}

static int computeNLabels(TC_PARMS *parms)
{
  int width,height,depth,i,j,k,x,y,z,xinit,yinit,zinit,val,n,nlabels,*tab;
  MRI *mri_seg=parms->mri_seg,*mri=parms->mri_bin;
  int nb,count;

  width=mri->width;
  height=mri->height;
  depth=mri->depth;

  xinit=parms->region.x;
  yinit=parms->region.y;
  zinit=parms->region.z;

  nlabels=parms->nblabels;
  tab=parms->labels;

  nb=0;count=0;
  for(k=2;k<depth-2;k++)
    for(j=2;j<height-2;j++)
      for(i=2;i<width-2;i++)
	{
	  x=xinit+i-2;
	  y=yinit+j-2;
	  z=zinit+k-2;
	  val=1;
	  for(n=0;n<nlabels;n++)
	    if(MRIvox(mri_seg,x,y,z)==tab[n])
	      val=0;
	  if(val && (MRIvox(mri,i,j,k)==F_B||MRIvox(mri,i,j,k)==F_R))
	    nb+=1;
	  if(!val && (MRIvox(mri,i,j,k)==B_B||MRIvox(mri,i,j,k)==B_R))
	     nb+=1;
	  if(MRIvox(mri,i,j,k)==F_B||MRIvox(mri,i,j,k)==F_R)
	    count++;
	}

  fprintf(stderr,"\n**********************************************************"
	         "\n**********************************************************"
	  "\n%d voxels have been changed: %3.3f for the label %d\n"
	  ,nb,100.*nb/count,parms->labels[0]);

  return nb;
}

static void SmoothImage(MRI *mri,TC_PARMS *parms)
{
  MRI *mri_bin=parms->mri_bin,*mri_tmp;
  int width,height,depth,i,j,k,a,b,c;
  float average,val;

  mri_tmp=MRIclone(mri,NULL);
  MRIcopy(mri,mri_tmp);

  width=mri->width;
  height=mri->height;
  depth=mri->depth;

  for(k=1;k<depth-1;k++)
    for(j=1;j<height-1;j++)
      for(i=1;i<width-1;i++)
	if(MRIvox(mri_bin,i,j,k)==F_R) //inside the volume
	  {
	    average=0;val=0.0f;
	    for(a=-1;a<2;a++)
	      for(b=-1;b<2;b++)
		for(c=-1;c<2;c++)
		  if(MRIvox(mri_bin,i+a,j+b,k+c)==F_R)
		    {
		      average+=MRIFvox(mri_tmp,i+a,j+b,k+c);
		      val+=1.0f;
		    }
	    MRIFvox(mri,i,j,k)=average/val;
	  }	
	else if (MRIvox(mri_bin,i,j,k)==B_R)
	  {
	    average=0;val=0.0f;
	    for(a=-1;a<2;a++)
	      for(b=-1;b<2;b++)
		for(c=-1;c<2;c++)
		  if(MRIvox(mri_bin,i+a,j+b,k+c)==B_R)
		    {
		      average+=MRIFvox(mri_tmp,i+a,j+b,k+c);
		      val+=1.0f;
		    }
	    MRIFvox(mri,i,j,k)=average/val;
	  }	
  
  MRIfree(&mri_tmp);
}

#define PI 3.14159
#if 0
static double gauss(double x,double m,double v)
{
  return (1./(sqrt((double)2.*PI*v))*exp(-(double)0.5*SQR((x-m))/v));
}
#endif
#define NLABELS 50
static void initProb(TC_PARMS *parms)
{
  GCA *gca;
  GCA_PRIOR *gcap;
  GCA_NODE *gcan;
  TRANSFORM *transform;
  MRI *mri_orig=parms->mri_orig;
  MRI *mri_prob=parms->mri_prob;
  int label,width,height,depth,i,j,k,l,m,a,b,c,x,y,z;
  int xinit,yinit,zinit,n,nlabels=parms->nblabels,*tab=parms->labels;
  float ps,pcs,pis,pi,psi,pcsi;
  char fname[100];
  MRI *mri_ps,*mri_pcs,*mri_pis,*mri_pi,*mri_psi,*mri_pcsi;
  float intensity;
  
  double Ps1[NLABELS],Pis1[NLABELS],Psi1[NLABELS];
  double Ps2[NLABELS],Pis2[NLABELS],Psi2[NLABELS];
  int nlabels1,nlabels2,Labels1[NLABELS],Labels2[NLABELS],test;
  float maxp1,maxp2;

  if(!parms->priors)
    return;

  gca=GCAread(parms->gca_fname);
  if(!gca)
    {
      fprintf(stderr,"\nIMPOSSIBLE TO READ GCA from %s",parms->gca_fname);
      parms->alpha=1.0f;
      parms->beta=1.0f;
      parms->priors=0;
      return;
    }
  transform=TransformRead(parms->transform_fname);
  if(!transform)
    {
      GCAfree(&gca);
      fprintf(stderr,"\nIMPOSSIBLE TO READ TRANSFORM from %s",parms->transform_fname);
      parms->beta=1.0f;
      parms->alpha=1.0f;
      parms->priors=0;
      return;
    }
  
  width=parms->width;
  height=parms->height;
  depth=parms->depth;

  xinit=parms->region.x;
  yinit=parms->region.y;
  zinit=parms->region.z;

  if(parms->mapsfname)
    {
      mri_ps=MRIalloc(width,height,depth,MRI_FLOAT);
      mri_pcs=MRIalloc(width,height,depth,MRI_FLOAT);
      mri_pis=MRIalloc(width,height,depth,MRI_FLOAT);
      mri_pi=MRIalloc(width,height,depth,MRI_FLOAT);
      mri_psi=MRIalloc(width,height,depth,MRI_FLOAT);
      mri_pcsi=MRIalloc(width,height,depth,MRI_FLOAT);
    }
  

  TransformInvert(transform, mri_orig) ;
  for(k=0;k<=depth-2;k++)
    for(j=0;j<=height-2;j++)
      for(i=0;i<=width-2;i++)
	{
	  x=xinit+i-1;y=yinit+j-1;z=zinit+k-1;

	  ps=pis=0;
	  GCAsourceVoxelToPrior(gca,mri_orig,transform, x, y, z,&a,&b,&c) ;
	  gcap=&gca->priors[a][b][c];
	  GCAsourceVoxelToNode(gca, mri_orig, transform, x, y, z, &a, &b, &c) ;
	  gcan=&gca->nodes[a][b][c];
	  intensity=(float)MRIvox(mri_orig,x,y,z);

	  nlabels1=nlabels2=0;
	  //compute prior prob
	  maxp1=maxp2=0;
	  for(l=0;l<gcap->nlabels;l++)
	    {
	      test=1;
	      label=gcap->labels[l];
	      for(n=0;n<nlabels;n++)
		if(tab[n]==label)
		  {
		    Labels1[nlabels1]=label;
		    Ps1[nlabels1]=gcap->priors[l];
		    nlabels1++;
		    test=0;
		    if(gcap->priors[l]>maxp1)
		      maxp1=gcap->priors[l];
		  }
	      if(test)
		{
		  Labels2[nlabels2]=label;
		  Ps2[nlabels2]=gcap->priors[l];
		  nlabels2++;
		  if(gcap->priors[l]>maxp2)
		    maxp2=gcap->priors[l];
		} 
	    }
	  ps=maxp1;
	  pcs=maxp2;
	  //compute p(i|s)
	  maxp1=0;
	  for(n=0;n<nlabels1;n++) 
	    for(m=0;m<gcan->nlabels;m++)
	      if(Labels1[n]==gcan->labels[m])
		{
			// Pis1[n]=PROB(gauss(intensity,gcan->gcs[m].mean,gcan->gcs[m].var)); 
		  if(Pis1[n]>maxp1)
		    maxp1=Pis1[n];
		}
	  maxp2=0;
	  for(n=0;n<nlabels2;n++) 
	    for(m=0;m<gcan->nlabels;m++)
	      if(Labels2[n]==gcan->labels[m])
		{
		  //Pis2[n]=PROB(gauss(intensity,gcan->gcs[m].mean,gcan->gcs[m].var)); 
		  if(Pis2[n]>maxp2)
		    maxp2=Pis2[n];
		}
	  pis=maxp1;
	  //compute p(i)
	  pi=0.0f;
	  for(l=0;l<nlabels1;l++)
	      pi+=Pis1[l]*Ps1[l];
	  for(l=0;l<nlabels2;l++)
	      pi+=Pis2[l]*Ps2[l];
	  pi=PROB(pi);

	  //compute p(s|i)
	  maxp1=maxp2=0;
	  for(l=0;l<nlabels1;l++)
	    {
	      Psi1[l]=Pis1[l]*Ps1[l]/pi;
	      if(Psi1[l]>maxp1)
		maxp1=Psi1[l];
	    }
	  for(l=0;l<nlabels2;l++)
	    {
	      Psi2[l]=Pis2[l]*Ps2[l]/pi;
	      if(Psi2[l]>maxp2)
		maxp2=Psi2[l];
	    }
	  psi=maxp1;
	  pcsi=maxp2;

	  if(!nlabels)
	    {
	      psi=0;
	      pcsi=1;
	    }
    	  
	  psi=PROB(psi);
	  pcsi=PROB(pcsi);

	  MRIFvox(mri_prob,i+1,j+1,k+1)=psi;
	  MRIFvox(parms->mri_cprob,i+1,j+1,k+1)=pcsi;

	  if(parms->mapsfname)
	    {
	      MRIFvox(mri_ps,i+1,j+1,k+1)=ps;
	      MRIFvox(mri_pcs,i+1,j+1,k+1)=pcs;
	      MRIFvox(mri_pis,i+1,j+1,k+1)=pis;
	      MRIFvox(mri_pi,i+1,j+1,k+1)=pi;
	      MRIFvox(mri_psi,i+1,j+1,k+1)=psi;
	      MRIFvox(mri_pcsi,i+1,j+1,k+1)=pcsi;
	    }
	}


  if(parms->mapsfname)
    {
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_ps.mgh");
      MRIwrite(mri_ps,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_pcs.mgh");
      MRIwrite(mri_pcs,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_pis.mgh");
      MRIwrite(mri_pis,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_pi.mgh");
      MRIwrite(mri_pi,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_psi.mgh");
      MRIwrite(mri_psi,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_pcsi.mgh");
      MRIwrite(mri_pcsi,fname);
      
      MRIfree(&mri_ps);
      MRIfree(&mri_pcs);
      MRIfree(&mri_pis);
      MRIfree(&mri_pi);
      MRIfree(&mri_psi);
      MRIfree(&mri_pcsi);
    }

  GCAfree(&gca);
  free(transform);
}

static void initCostMaps(TC_PARMS *parms)
{
  int i,j,k,width,depth,height;
  MRI *mri=parms->mri_bin,*mri_fcost,*mri_bcost,*mri_psi,*mri_pcsi,*mri_dist;
  char fname[100];

  width=mri->width;
  height=mri->height;
  depth=mri->depth;
  
  mri_fcost=MRIalloc(width, height, depth, MRI_FLOAT);
  parms->mri_fcost=mri_fcost;
  mri_bcost=MRIalloc(width, height, depth, MRI_FLOAT);
  parms->mri_bcost=mri_bcost;
  
  mri_psi=parms->mri_prob;
  mri_pcsi=parms->mri_cprob;
  mri_dist=parms->mri_dist;

  switch(parms->costmode)
    {
    case VOXEL_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fcost,i,j,k)=1;
	      MRIFvox(mri_bcost,i,j,k)=1;
	    }
      break;
    case PROB_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fcost,i,j,k)=MRIFvox(mri_psi,i,j,k);
	      MRIFvox(mri_bcost,i,j,k)=MRIFvox(mri_pcsi,i,j,k);
	    }
      break;
    case MAP_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fcost,i,j,k)=log(MRIFvox(mri_psi,i,j,k)/MRIFvox(mri_pcsi,i,j,k));
	      MRIFvox(mri_bcost,i,j,k)=-MRIFvox(mri_fcost,i,j,k);
	    }
      break;
    case PROB_MAP_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fcost,i,j,k)=log(MRIFvox(mri_psi,i,j,k)/MRIFvox(mri_pcsi,i,j,k));
	      MRIFvox(mri_bcost,i,j,k)=-MRIFvox(mri_fcost,i,j,k);
	    }
      break;
    case NORMAL_MODE:   //only based on the distance function
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fcost,i,j,k)=0.5*(1.+MRIFvox(mri_dist,i,j,k));
	      MRIFvox(mri_bcost,i,j,k)=1.0f-MRIFvox(mri_fcost,i,j,k);
	    }
      break;
    }

  if(parms->mapsfname)
    {
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_fc.mgh");
      MRIwrite(mri_fcost,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_bc.mgh");
      MRIwrite(mri_bcost,fname);
    }

}

static void initPriorMaps(TC_PARMS *parms)
{
  int i,j,k,width,depth,height;
  MRI *mri=parms->mri_bin,*mri_dist=parms->mri_dist,*mri_fprior,*mri_bprior,*mri_psi,*mri_pcsi;
  char fname[100];

  width=mri->width;
  height=mri->height;
  depth=mri->depth;
  
  mri_fprior=MRIalloc(width, height, depth, MRI_FLOAT);
  parms->mri_fprior=mri_fprior;
  mri_bprior=MRIalloc(width, height, depth, MRI_FLOAT);
  parms->mri_bprior=mri_bprior;
  
  mri_psi=parms->mri_prob;
  mri_pcsi=parms->mri_cprob;

  switch(parms->costmode)
    {
    case MAP_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fprior,i,j,k)=0.5
		+0.5*tanh(0.1*log(MRIFvox(mri_psi,i,j,k)/MRIFvox(mri_pcsi,i,j,k)));
	      MRIFvox(mri_bprior,i,j,k)=0.5
		+0.5*tanh(0.1*log(MRIFvox(mri_pcsi,i,j,k)/MRIFvox(mri_psi,i,j,k)));
	    }
      break;
    case NORMAL_MODE:
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fprior,i,j,k)=0.5*(1.+MRIFvox(mri_dist,i,j,k));
	      MRIFvox(mri_bprior,i,j,k)=1.0f-MRIFvox(mri_fprior,i,j,k);
	    }
      break;
    default:          
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      MRIFvox(mri_fprior,i,j,k)=MRIFvox(mri_psi,i,j,k);
	      MRIFvox(mri_bprior,i,j,k)=MRIFvox(mri_pcsi,i,j,k);
	    }
      break;
    }
  if(parms->mapsfname)
    {
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_fpr.mgh");
      MRIwrite(mri_fprior,fname);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_bpr.mgh");
      MRIwrite(mri_bprior,fname);
    }
}

static void guessSegmentation(TC_PARMS *parms)
{
  GCA *gca;
  GCA_PRIOR *gcap;
  GCA_NODE *gcan;
  TRANSFORM *transform;
  MRI *mri_orig=parms->mri_orig;
  MRI *mri_seg;
  int width,height,depth,i,j,k,l,m,a,b,c,x,y,z,label;
  int xmin,ymin,zmin,xmax,ymax,zmax;
  float ps,pcs,pis,pi,psi,pcsi;
  float intensity;
  
  double Ps1[NLABELS],Pis1[NLABELS],Psi1[NLABELS];
  double Ps2[NLABELS],Pis2[NLABELS],Psi2[NLABELS];
  int nlabels=parms->nblabels,*tab=parms->labels,Labels1[NLABELS],Labels2[NLABELS];
  int nlabels1,nlabels2,test,n,val;
  float maxp1,maxp2;

  if(!parms->priors)
    return;

  if(parms->mri_seg)
    MRIfree(&parms->mri_seg);

  gca=GCAread(parms->gca_fname);
  if(!gca)
    {
      Error("IMPOSSIBLE TO READ GCA");
    }
  transform=TransformRead(parms->transform_fname);
  if(!transform)
    {
      GCAfree(&gca);
      Error("IMPOSSIBLE TO READ TRANSFORM");
    }
  
  width=mri_orig->width;
  height=mri_orig->height;
  depth=mri_orig->depth;

  mri_seg=MRIalloc(width,height,depth,MRI_UCHAR);

  TransformInvert(transform, mri_orig) ;
	  
  xmin=ymin=zmin=256;
  xmax=ymax=zmax=0;
  //find the bounding box for which prior > 0
  for(k=0;k<=depth;k++)
    for(j=0;j<=height;j++)
      for(i=0;i<=width;i++)
	{
	  x=i;y=j;z=k;
	  GCAsourceVoxelToPrior(gca,mri_orig,transform, x, y, z,&a,&b,&c) ;
	  gcap=&gca->priors[a][b][c];
	  //compute prior prob
	  ps=0;
	  for(l=0;l<gcap->nlabels;l++)
	    for(n=0;n<nlabels;n++)
	      if(gcap->labels[l]==tab[n] && ps<gcap->priors[l])
		ps=gcap->priors[l];	      
	  if(ps>0)
	    {
	      if(xmin>i) xmin=i;
	      if(ymin>j) ymin=j;
	      if(zmin>k) zmin=k;
	      if(xmax<i) xmax=i;
	      if(ymax<j) ymax=j;
	      if(zmax<k) zmax=k;
	    }
	}

  for(k=zmin;k<=zmax;k++)
    for(j=ymin;j<=ymax;j++)
      for(i=xmin;i<=xmax;i++)
	{
	  x=i;y=j;z=k;
	  
	  GCAsourceVoxelToPrior(gca,mri_orig,transform, x, y, z,&a,&b,&c) ;
	  gcap=&gca->priors[a][b][c];
	  GCAsourceVoxelToNode(gca, mri_orig, transform, x, y, z, &a, &b, &c) ;
	  gcan=&gca->nodes[a][b][c];
	  intensity=(float)MRIvox(mri_orig,x,y,z);
	  
	  nlabels1=nlabels2=0;
	  //compute prior prob
	  maxp1=maxp2=0;
	  for(l=0;l<gcap->nlabels;l++)
	    {
	      test=1;
	      label=gcap->labels[l];
	      for(n=0;n<nlabels;n++)
		if(tab[n]==label)
		  {
		    Labels1[nlabels1]=label;
		    Ps1[nlabels1]=gcap->priors[l];
		    nlabels1++;
		    test=0;
		    if(gcap->priors[l]>maxp1)
		      maxp1=gcap->priors[l];
		  }
	      if(test)
		{
		  Labels2[nlabels2]=label;
		  Ps2[nlabels2]=gcap->priors[l];
		  nlabels2++;
		  if(gcap->priors[l]>maxp2)
		    maxp2=gcap->priors[l];
		} 
	    }
	  ps=maxp1;
	  pcs=maxp2;
	  //compute p(i|s)
	  maxp1=0;
	  for(n=0;n<nlabels1;n++) 
	    for(m=0;m<gcan->nlabels;m++)
	      if(Labels1[n]==gcan->labels[m])
		{
		  //		  Pis1[n]=PROB(gauss(intensity,gcan->gcs[m].mean,gcan->gcs[m].var)); 
		  if(Pis1[n]>maxp1)
		    maxp1=Pis1[n];
		}
	  maxp2=0;
	  for(n=0;n<nlabels2;n++) 
	    for(m=0;m<gcan->nlabels;m++)
	      if(Labels2[n]==gcan->labels[m])
		{
		  // Pis2[n]=PROB(gauss(intensity,gcan->gcs[m].mean,gcan->gcs[m].var)); 
		  if(Pis2[n]>maxp2)
		    maxp2=Pis2[n];
		}
	  pis=maxp1;
	  //compute p(i)
	  pi=0.0f;
	  for(l=0;l<nlabels1;l++)
	      pi+=Pis1[l]*Ps1[l];
	  for(l=0;l<nlabels2;l++)
	      pi+=Pis2[l]*Ps2[l];
	  pi=PROB(pi);

	  //compute p(s|i)
	  maxp1=maxp2=0;
	  val=0;
	  for(l=0;l<nlabels1;l++)
	    {
	      Psi1[l]=Pis1[l]*Ps1[l]/pi;
	      if(Psi1[l]>maxp1)
		{
		  val=Labels1[l];
		  maxp1=Psi1[l];
		}
	    }
	  for(l=0;l<nlabels2;l++)
	    {
	      Psi2[l]=Pis2[l]*Ps2[l]/pi;
	      if(Psi2[l]>maxp2)
		maxp2=Psi2[l];
	    }
	  psi=maxp1;
	  pcsi=maxp2;

	  if(!nlabels)
	    {
	      psi=0;
	      pcsi=1;
	    }
    	  
	  psi=PROB(psi);
	  pcsi=PROB(pcsi);

	  if(psi>pcsi)
	    MRIvox(mri_seg,i,j,k)=val;
	}
  parms->mri_seg=mri_seg;
  GCAfree(&gca);
  free(transform);

  MRIwrite(mri_seg,"/tmp/tmp");
}

static void initImages(TC_PARMS* parms)
{
  int i,j,k,width,depth,height,xinit,yinit,zinit,a,b,c;
  MRI *mri,*mri_bin,*mri_dist,*mri_prob,*mri_cprob,*mri_labeled;
  int changed,con,sum,absval,val;
  int max_dist_int,max_dist_out;
  float fval;
  int count,nlabels,n;
  
  nlabels=parms->nblabels;
  
  fprintf(stderr,"\n****************************************************");
  fprintf(stderr,"\nINITIALIZATION OF THE IMAGES");


  if(parms->guess)
    guessSegmentation(parms);

  mri=parms->mri_seg;

  width=mri->width;
  height=mri->height;
  depth=mri->depth;

#if 0
  MRIvox(mri,140,85,109)=parms->labels[0];
  MRIvox(mri,140,84,109)=parms->labels[0];
  MRIvox(mri,140,84,110)=parms->labels[0];
  MRIvox(mri,140,84,111)=parms->labels[0];
  MRIvox(mri,140,85,112)=parms->labels[0];
  MRIvox(mri,140,85,113)=parms->labels[0];
  MRIvox(mri,140,85,114)=parms->labels[0];
  MRIvox(mri,140,85,115)=parms->labels[0];
  MRIvox(mri,140,85,114)=parms->labels[0];
  MRIvox(mri,140,86,115)=parms->labels[0];
#if 0  
  MRIvox(mri,136,89,99)=parms->labels[0];
  MRIvox(mri,136,89,98)=parms->labels[0];
  MRIvox(mri,146,88,98)=parms->labels[0];

  MRIvox(mri,129,96,124)=parms->labels[0];
  MRIvox(mri,129,97,124)=parms->labels[0];
  MRIvox(mri,152,96,100)=parms->labels[0];
  MRIvox(mri,152,95,100)=parms->labels[0];
  MRIvox(mri,140,99,86)=parms->labels[0];
  MRIvox(mri,140,99,85)=parms->labels[0];
  MRIvox(mri,140,98,85)=parms->labels[0];
  MRIvox(mri,141,98,86)=parms->labels[0];
#endif
#endif

  count=0;
  //find region and extract MRI
  parms->region.x=parms->region.y=parms->region.z=1000;
  parms->region.dx=parms->region.dy=parms->region.dz=0;
  for(k=0;k<depth;k++)
    for(j=0;j<height;j++)
      for(i=0;i<width;i++)
	{
	  val=MRIvox(mri,i,j,k);
	  for(n=0;n<nlabels;n++)
	    if(val==parms->labels[n])
	      {
		count++;
		if(parms->region.x>i)
		  parms->region.x=i;
		if(parms->region.dx<i)
		  parms->region.dx=i;
		if(parms->region.y>j)
		  parms->region.y=j;
		if(parms->region.dy<j)
		  parms->region.dy=j;
		if(parms->region.z>k)
		  parms->region.z=k;
		if(parms->region.dz<k)
		  parms->region.dz=k;
	      }
	}
  if(count==0)
    Error("NO LABEL FOUND!\n");

  //  parms->region.z=80;
  // parms->region.dz=200;

  parms->region.dz-=parms->region.z;
  parms->region.dy-=parms->region.y;
  parms->region.dx-=parms->region.x;
  
  //allocate and init binary volume
  width=parms->region.dx+5;
  height=parms->region.dy+5;
  depth=parms->region.dz+5;
  
  parms->width=width;
  parms->height=height;
  parms->depth=depth;

  parms->mri_bin=MRIalloc(width, height, depth, MRI_UCHAR);
  mri_bin=parms->mri_bin;
  
  xinit=parms->region.x;
  yinit=parms->region.y;
  zinit=parms->region.z;
  
  for(k=0;k<depth-4;k++)
    for(j=0;j<height-4;j++)
      for(i=0;i<width-4;i++)
	{	  
	  val=MRIvox(mri,xinit+i,yinit+j,zinit+k);
	  MRIvox(mri_bin,i+2,j+2,k+2)=B_R;
	  for(n=0;n<nlabels;n++)
	    if(val==parms->labels[n])
	      MRIvox(mri_bin,i+2,j+2,k+2)=F_R;
	}	  
  
  //take care of the exterior voxels
  setBorderValue(mri_bin,0,0); 

  //allocate and init dist map volume
  parms->mri_dist=MRIalloc(width, height, depth, MRI_FLOAT);
  mri_dist=parms->mri_dist;
  for(k=1;k<depth-1;k++)
    for(j=1;j<height-1;j++)
      for(i=1;i<width-1;i++)
	{
	  val=MRIvox(mri_bin,i,j,k);
	  if(val==F_R && checkNbh(mri_bin,i,j,k,B_R,parms->f_c))
	    MRIFvox(mri_dist,i,j,k)=1.0f;
	  else if(val==B_R && checkNbh(mri_bin,i,j,k,F_R,parms->b_c))
	    MRIFvox(mri_dist,i,j,k)=-1.0f;
	  else
	    MRIFvox(mri_dist,i,j,k)=0.0f;
	}

  changed=1;
  max_dist_int=0;
  max_dist_out=0;
  while(changed)
    {
      changed=0;
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    {
	      if(MRIFvox(mri_dist,i,j,k))
		continue;
	      if(MRIvox(mri_bin,i,j,k)==F_R)  //inside the volume
		{
		  con=associatedConnectivity(parms->f_c);
		  absval=1000;
		  for(a=-1;a<2;a++)
		    for(b=-1;b<2;b++)
		      for(c=-1;c<2;c++)
			{
			  sum=abs(a)+abs(b)+abs(c);
			  if(sum>con || (!sum))
			    continue;
			  val=MRIFvox(mri_dist,i+a,j+b,k+c);
			  if(!val)
			    continue;
			  if(absval>abs(val))
			    absval=abs(val);
			}
		  if(absval!=1000)
		    {
		      MRIFvox(mri_dist,i,j,k)=absval+1;
		      changed=1;
		      if(absval+1>max_dist_int)
			max_dist_int=absval+1;
		    }
		}
	      else if(MRIvox(mri_bin,i,j,k)==B_R) //outside the volume
		{
		  con=associatedConnectivity(parms->b_c);
		  absval=1000;
		  for(a=-1;a<2;a++)
		    for(b=-1;b<2;b++)
		      for(c=-1;c<2;c++)
			{
			  sum=abs(a)+abs(b)+abs(c);
			  if(sum>con || (!sum))
			    continue;
			  val=MRIFvox(mri_dist,i+a,j+b,k+c);
			  if(!val)
			    continue;
			  if(absval>abs(val))
			    absval=abs(val);
			}
		  if(absval!=1000)
		    {
		      MRIFvox(mri_dist,i,j,k)=-absval-1;
		      changed=1;
		      if(absval+1>max_dist_out)
			max_dist_out=absval+1;
		    }
		}
	    }
    }
+
  parms->maxindist=max_dist_int;
  parms->maxoutdist=max_dist_out;

  SmoothImage(mri_dist,parms);

  //finish updating distance image
  for(k=1;k<depth-1;k++)
    for(j=1;j<height-1;j++)
      for(i=1;i<width-1;i++)
	{
	  fval=MRIFvox(mri_dist,i,j,k);
	  if(fval>=0)
	    MRIFvox(mri_dist,i,j,k)=fval/parms->maxindist;
	  else
	    MRIFvox(mri_dist,i,j,k)=fval/parms->maxoutdist;
	}
	  

  //allocate and init mri_prob volume (linear scaling so far)
  mri_prob=MRIalloc(width, height, depth, MRI_FLOAT);
  parms->mri_prob=mri_prob;
  mri_cprob=MRIalloc(width, height, depth, MRI_FLOAT);
  parms->mri_cprob=mri_cprob;
  
  if(parms->priormap)
    {
      MRI *mri_tmp;
      mri_tmp=MRIread(parms->priormap);
      for(k=0;k<depth-4;k++)
	for(j=0;j<height-4;j++)
	  for(i=0;i<width-4;i++)
	    {
	      MRIFvox(mri_prob,i+2,j+2,k+2)=MRIFvox(mri_tmp,xinit+i,yinit+j,zinit+k);
	      MRIFvox(mri_cprob,i+2,j+2,k+2)=1.0f-MRIFvox(mri_prob,i+2,j+2,k+2);
	    }
      MRIfree(&mri_tmp);
    }
  else
    initProb(parms);
 
  //allocate and init mri_labeled volume
  mri_labeled=MRIalloc(width, height, depth, MRI_SHORT);
  for(k=0;k<depth;k++)
    for(j=0;j<height;j++)
      for(i=0;i<width;i++)
	MRISvox(mri_labeled,i,j,k)=-1;
  parms->mri_labeled=mri_labeled;

  //finish allocating
  parms->list=initList(BIN_NBR,parms->width,parms->height,parms->depth);
  initCostMaps(parms);
  initPriorMaps(parms);
  initCellsFromMap(parms);
  if(parms->mapsfname)
    SaveOrigMaps(parms);
}

static int mriChangeLabel(MRI *mri,int src,int dst)
{
  int i,j,k,width,height,depth,count;
  
  width=mri->width;
  height=mri->height;
  depth=mri->depth;
  
  count=0;
  for(k=0;k<depth;k++)
    for(j=0;j<height;j++)
      for(i=0;i<width;i++)
	if(MRIvox(mri,i,j,k)==src)
	  {
	    MRIvox(mri,i,j,k)=dst;
	    count++;
	  }
  //    fprintf(stderr,"\n%d voxels changed of type %d to type %d \n",count,src,dst);
  return NO_ERROR;
}

static void modifyImage(TC_PARMS *parms)
{
  int i,j,k;
  int width,height,depth;
  int xinit,yinit,zinit;
  char fname[100];
  MRI *mri=parms->mri_bin,*mri_output=parms->mri_output;
  MRI *mri_final;

  if(!mri_output)
    mri_output=MRIclone(parms->mri_seg,NULL);
  parms->mri_output=mri_output;

  width=parms->region.dx+5;
  height=parms->region.dy+5;
  depth=parms->region.dz+5;
  
  xinit=parms->region.x;
  yinit=parms->region.y;
  zinit=parms->region.z;
  
  for(k=0;k<depth-4;k++)
    for(j=0;j<height-4;j++)
      for(i=0;i<width-4;i++)
	if(MRIvox(mri,i+2,j+2,k+2)==F_B)
	  MRIvox(mri_output,xinit+i,yinit+j,zinit+k)=1;
  
  if(parms->mapsfname)
    {
      mri_final=MRIalloc(width,height,depth,MRI_FLOAT);
      for(k=1;k<depth-1;k++)
	for(j=1;j<height-1;j++)
	  for(i=1;i<width-1;i++)
	    if(MRIvox(mri,i,j,k)==F_B)
	      MRIFvox(mri_final,i,j,k)=1;
      
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_output.mgh");
      MRIwrite(mri_final,fname);
      MRIfree(&mri_final);
    }

	
}

static MSV* findPoint(TC_PARMS *parms, MSV *point)
{
  MSV* pt;
  MRI *mri=parms->mri_bin,*mri_labeled=parms->mri_labeled;
  int border_labels[27],nlabels;
  int x,y,z,width,height,depth,s,a,b,c,val,n,i,j,k;
  ConnectedComponent *cc;
  int connectivity=parms->c_c,associatedconnectivity;
  Nbh fnbh,bnbh;
  int fgtp,bgtp,con,sum;
  float maxprior,topo;
  Cell*** table=parms->list->table;

  con=connectivityNumber(connectivity);
  parms->nlabels=0;
  associatedconnectivity=associatedConnectivity(connectivity);
  if(!point)
    pt=(MSV*)malloc(sizeof(MSV));
  else
    pt=point;
  
  width=mri->width;
  height=mri->height;
  depth=mri->depth;

  if(!parms->multiplemode)
    {
      maxprior=-1;
      for (z = 1 ; (z<depth-1) ; z++)
	for (y = 1 ; (y < height-1) ; y++)
	  for (x = 1 ; (x < width-1) ; x++)
	    if(MRIvox(mri,x,y,z)==RESIDUE && table[z][y][x].prior>maxprior)
	      {	 
		if(checkNbh(mri,x,y,z,BODY,connectivity))
		  {
		    loadNbh(mri,&fnbh,x,y,z,BODY);
		    reverseNbh(&fnbh,&bnbh);
		    
		    fgtp=checkTn(&fnbh,&fnbh,connectivity);
		    bgtp=checkTn(&bnbh,&bnbh,associatedconnectivity);
		    
		    if(fgtp==1 && bgtp==1)
		      {
			maxprior=table[z][y][x].prior;
			pt->x=x;
			pt->y=y;
			pt->z=z;
		      }
		  }
		else
		  {
		    loadNbh(mri,&fnbh,x,y,z,RESIDUE);
		    if(checkTn(&fnbh,&fnbh,associatedconnectivity)==0)
		      {
			MRIvox(mri,x,y,z)=BODY;
			s=componentNew(parms->segmentation);
			addCellToComponent(&parms->segmentation->components[s],&parms->list->table[z][y][x]);
			MRISvox(parms->mri_labeled,x,y,z)=s;
			addComponent(parms->ccs,&parms->segmentation->components[s]);
		      }
		    else
		      {
			maxprior=table[z][y][x].prior;
			pt->x=x;
			pt->y=y;
			pt->z=z;
		      } 
		  }
	      };
    }
  else
    {
      maxprior=-1;
      for (z = 1 ; (z<depth-1) ; z++)
	for (y = 1 ; (y < height-1) ; y++)
	  for (x = 1 ; (x < width-1) ; x++)
	    if(MRIvox(mri,x,y,z)==RESIDUE && table[z][y][x].prior>maxprior)
	      {	 
		if(checkNbh(mri,x,y,z,BODY,connectivity))
		  {
		    //first we have to find all the neighbors
		    nlabels=0;
		    for(a=-1;a<2;a++)
		      for(b=-1;b<2;b++)
			for(c=-1;c<2;c++)
			  {
			    sum=abs(a)+abs(b)+abs(c);
			    if(!sum && sum>con)
			      continue;
			    i=x+a;j=y+b;k=z+c;
			    if(MRIvox(mri,i,j,k)!=BODY)
			      continue;
			    val=MRISvox(mri_labeled,i,j,k);
			    cc=&parms->segmentation->components[val];
			    if(cc->found==0)
			      {
				border_labels[nlabels]=val;
				nlabels++;
				cc->found=1;
			      }
			  }
		    //then check if topology is correct with each of them
		    topo=1;
		    for(n=0;topo&&n<nlabels;n++)
		      {
			val=border_labels[n];
			//load Nbh
			for(a=-1;a<2;a++)
			  for(b=-1;b<2;b++)
			    for(c=-1;c<2;c++)
			      {
				i=x+a;j=y+b;k=z+c;
				if(MRIvox(mri,i,j,k)==BODY &&
				   MRISvox(mri_labeled,i,j,k)==val)
				  fnbh[1+a][1+b][1+c]=1;
				else
				  fnbh[1+a][1+b][1+c]=0;
			      }
			reverseNbh(&fnbh,&bnbh);
			if((checkTn(&fnbh,&fnbh,connectivity)!=1) ||
			   (checkTn(&bnbh,&bnbh,associatedconnectivity)!=1))
			  topo=0;
		      }
		    //then reinit the components to found=0;
		    for(n=0;n<nlabels;n++)
		      parms->segmentation->components[border_labels[n]].found=0;
		    
		    if(topo)   //in this case every comp doesn't change the topology
		      {	
			maxprior=table[z][y][x].prior;
			for(n=0;n<nlabels;n++)
			  parms->border_labels[n]=border_labels[n];
			parms->nlabels=nlabels;
			pt->x=x;
			pt->y=y;
			pt->z=z;
		      }
		  }
		else
		  {
		    loadNbh(mri,&fnbh,x,y,z,RESIDUE);
		    if(checkTn(&fnbh,&fnbh,associatedconnectivity)==0)
		      {
			MRIvox(mri,x,y,z)=BODY;
			s=componentNew(parms->segmentation);
			addCellToComponent(&parms->segmentation->components[s],&parms->list->table[z][y][x]);
			addComponent(parms->ccs,&parms->segmentation->components[s]);
			MRISvox(parms->mri_labeled,x,y,z)=s;
		      }
		    else
		      {
			parms->nlabels=0;
			maxprior=table[z][y][x].prior;
			pt->x=x;
			pt->y=y;
			pt->z=z;
		      } 
		  }
	      };
    }    
  if(maxprior==-1)
    {
      //      fprintf(stderr,"\nno point found!\n");
      return NULL;
    }
  //  fprintf(stderr,"\npoint (%d , %d , %d) found %d!\n",pt->x,pt->y,pt->z,parms->nlabels);
  return pt; 
}


static MSV* findMLVoxel(TC_PARMS *parms,MSV *point)
{
  MSV* pt;
  List *list=parms->list;
  Cell *cell;int k;
  Nbh fnbh,bnbh;
  int connectivity=parms->c_c;
  int associatedconnectivity,loop;
  MRI *mri=parms->mri_bin,*mri_labeled;
  int *border_labels;
  int nlabels,a,b,c,sum,con,x,y,z,val,n,topo;
  ConnectedComponent *cc;
  
  if(!list->ncells)
    return NULL;

  associatedconnectivity=associatedConnectivity(connectivity);

  if(!point)
    pt=(MSV*)malloc(sizeof(MSV));
  else
    pt=point;

  parms->nlabels=0;
  if(!parms->multiplemode) 
    {
      loop=1;
      for(k=list->reference_max_cell;loop&&k>=0;k--)
	{
	  cell=list->sorting[k];
	  while(loop && cell!=NULL)  //check if this cell is appropriate
	    {
	      switch(cell->type)
		{
		case SIMPLE:
		  pt->x=cell->x;
		  pt->y=cell->y;
		  pt->z=cell->z;
		  loop=0;
		  break;
		case NONSIMPLE:
		  cell=cell->next;
		  break;
		case UNKNOWN:
		  loadNbh(mri,&fnbh,cell->x,cell->y,cell->z,BODY);
		  reverseNbh(&fnbh,&bnbh);
		  if(checkTn(&fnbh,&fnbh,connectivity)==1 &&
		     checkTn(&bnbh,&bnbh,associatedconnectivity)==1)
		    {
		      cell->type=SIMPLE;
		      pt->x=cell->x;
		      pt->y=cell->y;
		      pt->z=cell->z;
		      loop=0;
		    }
		  else
		    {
		      cell->type=NONSIMPLE;
		      cell=cell->next;
		    }
		  break;
		}
	    }
	}
    }
  else
    {
      border_labels=parms->border_labels;
      con=connectivityNumber(connectivity);
      mri_labeled=parms->mri_labeled;
      loop=1;
      for(k=list->reference_max_cell;loop&&k>=0;k--)
	{
	  cell=list->sorting[k];
	  while(loop && cell!=NULL)  //check if this cell is appropriate
	    {
	      switch(cell->type)
		{
		case SIMPLE:
		  pt->x=cell->x;
		  pt->y=cell->y;
		  pt->z=cell->z;
		  loop=0;
		  break;
		case NONSIMPLE:
		  cell=cell->next;
		  break;
		case UNKNOWN:
		  //first we have to find all the neighbors
		  nlabels=0;
		  for(a=-1;a<2;a++)
		    for(b=-1;b<2;b++)
		      for(c=-1;c<2;c++)
			{
			  sum=abs(a)+abs(b)+abs(c);
			  if(!sum && sum>con)
			    continue;
			  x=cell->x+a;y=cell->y+b;z=cell->z+c;
			  if(MRIvox(mri,x,y,z)!=BODY)
			    continue;
			  val=MRISvox(mri_labeled,x,y,z);
			  cc=&parms->segmentation->components[val];
			  if(cc->found==0)
			    {
			      border_labels[nlabels]=val;
			      nlabels++;
			      cc->found=1;
			    }
			}
		  //then check if topology is correct with each of them
		  topo=1;
		  parms->nlabels=nlabels;
		  for(n=0;topo&&n<nlabels;n++)
		    {
		      val=border_labels[n];
		      //load Nbh
		      for(a=-1;a<2;a++)
			for(b=-1;b<2;b++)
			  for(c=-1;c<2;c++)
			    {
			      x=cell->x+a;y=cell->y+b;z=cell->z+c;
			      if(MRIvox(mri,x,y,z)==BODY &&
				 MRISvox(mri_labeled,x,y,z)==val)
				fnbh[1+a][1+b][1+c]=1;
			      else
				fnbh[1+a][1+b][1+c]=0;
			    }
		      reverseNbh(&fnbh,&bnbh);
		      if((checkTn(&fnbh,&fnbh,connectivity)!=1) ||
			 (checkTn(&bnbh,&bnbh,associatedconnectivity)!=1))
			topo=0;
		    }
		  //then reinit the components to found=0;
		  for(n=0;n<nlabels;n++)
		    parms->segmentation->components[border_labels[n]].found=0;
		  if(topo)   //in this case every comp doesn't change the topology
		    {
		      cell->type=SIMPLE;
		      pt->x=cell->x;
		      pt->y=cell->y;
		      pt->z=cell->z;
		      loop=0;
		    }
		  else
		    {
		      cell->type=NONSIMPLE;
		      cell=cell->next;
		    }
		  break;
		}
	    }
	}
    }

  if(loop)
    {
      if(!point)
	free(pt);
      return NULL;
    }
  return pt;
}

static int findNeighbors(TC_PARMS *parms,MSV *pt)
{
  int a,b,c,sum,con,val,x,y,z;
  List *list=parms->list;
  MRI *mri=parms->mri_bin;

  con=connectivityNumber(parms->c_c);
  for(a=-1;a<2;a++)
    for(b=-1;b<2;b++)
      for(c=-1;c<2;c++)
	{
	  x=pt->x+a;
	  y=pt->y+b;
	  z=pt->z+c;
	  val=MRIvox(mri,x,y,z);
	  if(val==RESIDUE || val==VISITED)
	    {
	      list->table[z][y][x].type=UNKNOWN;
	      //check
	    }
	  sum=abs(a)+abs(b)+abs(c);
	  if(!sum || sum>con)
	      continue;
 
	  if(val==RESIDUE)
	    {
	      addCell(list,&list->table[z][y][x]);
	      MRIvox(mri,x,y,z)=VISITED;
	    } 
	}

  return NO_ERROR;
}

static void changeLabelfromList(TC_PARMS *parms)
{
  List *list=parms->list;
  int k,x,y,z;
  Cell *cell,*ncell;
  MRI *mri=parms->mri_bin;

  for(k=0;k<list->number_of_bins;k++)
    {
      cell=list->sorting[k];
      while(cell)
	{
	  x=cell->x;
	  y=cell->y;
	  z=cell->z;
	  MRIvox(mri,x,y,z)=RESIDUE;
	  ncell=cell->next;
	  removeCell(list,cell);
	  cell=ncell;
	}
    }
}

static int segmentVoxel(TC_PARMS *parms,MSV *pt)
{
  int nlabels=parms->nlabels,*border_labels=parms->border_labels;
  SEGMENTATION *segmentation=parms->segmentation,*seg;
  ConnectedComponent *cc;
  int max,width,height,depth,nvox,x,y,z,label,a,b,c,val,sval,i,j,k,n;
  MRI*mri_labeled=parms->mri_labeled,*mri_bin=parms->mri_bin;
  CCS *ccs=parms->ccs;

  width=mri_labeled->width;
  height=mri_labeled->height;
  depth=mri_labeled->depth;

  x=pt->x;y=pt->y;z=pt->z;
  //first take care of the body component  
  switch(nlabels)
    {
    case 0:
      label=componentNew(segmentation);
      cc=&segmentation->components[label];
      break;    
    case 1:          /* assign this voxel to the one that it borders */
      label = border_labels[0] ;
      cc=&segmentation->components[label];
      removeComponent(ccs,cc);
      break ;
    default:         /* merge segments and assign to largest number of cells */
      max=-1;
      label=-1;
      for(nvox = 0 ; nvox < nlabels ; nvox++) //find largest number
	if (segmentation->components[border_labels[nvox]].ncells>max)
	  {
	    label=border_labels[nvox];
	    max=segmentation->components[border_labels[nvox]].ncells;
	  }
      cc=&segmentation->components[label];
      for(nvox = 0 ; nvox < nlabels ; nvox++) //merge into lowest number
	{
	  if(label==border_labels[nvox])
	    continue;
	  componentMerge(parms, label, border_labels[nvox]);
	}
      removeComponent(ccs,cc);
      break ;
    }
  /* add it to the existing list */
  parms->current_label=label;
  addCellToComponent(cc,&parms->list->table[z][y][x]);
  addComponent(ccs,cc);
  MRISvox(mri_labeled, x, y, z) = label ;
  MRIvox(mri_bin,x,y,z)=BODY;
  
  //then take care of the residual neighbors...
  if(parms->c_c==parms->f_c)
    {
      seg=parms->F_Rseg;
      ccs=parms->F_Rccs;
    }
  else
    {
      seg=parms->B_Rseg;
      ccs=parms->B_Rccs;
    }  

  for(a=-1;a<2;a++)
    for(b=-1;b<2;b++)
      for(c=-1;c<2;c++)
	{
	  val=MRIvox(mri_bin,x+a,y+b,z+c);
	  if(val==RESIDUE || val==VISITED)
	    {
	      sval=MRISvox(mri_labeled,x+a,y+b,z+c);
	      if(sval>=0)
		{
		  cc=&seg->components[sval];
		  for(n=0;n<cc->ncells;n++)
		    {
		      i=cc->cells[n]->x;
		      j=cc->cells[n]->y;
		      k=cc->cells[n]->z;
		      cc->cells[n]->type=UNKNOWN;
		      MRISvox(mri_labeled,i,j,k)=-1;
		    }
		  cc->ncells=0;
		  removeComponent(ccs,cc);
		  seg->ncomponents--;
		  cc->cost=-1;  //CHECK
		  cc->ncells=0;
		}
	    }
	}
  
  return NO_ERROR;
}

static void CTExpansion(TC_PARMS *parms)
{
  MSV pt;
  int x,y,z,count=0;
  List*list=parms->list;

  if(!parms->multiplemode)
    while(findMLVoxel(parms,&pt)) //find the max of the list
    {
      count++;
      fprintf(stderr,"\r   iteration n=%5d: ",count);
      //merge this voxel
      x=pt.x;y=pt.y;z=pt.z;
      addCellToComponent(&parms->segmentation->components[parms->current_label]
			 ,&parms->list->table[z][y][x]);
      MRISvox(parms->mri_labeled,x,y,z)=parms->current_label;
      MRIvox(parms->mri_bin,x,y,z)=BODY;
      removeCell(list,&list->table[z][y][x]);
      //find the neighbors of this point and update the list
      findNeighbors(parms,&pt);      
    }
  else
    while(findMLVoxel(parms,&pt)) //find the max of the list
    {
      count++;
      //      fprintf(stderr,"\r   iteration n=%5d: ",count);
      //merge this voxel
      x=pt.x;y=pt.y;z=pt.z;
      segmentVoxel(parms,&pt);
      removeCell(list,&list->table[z][y][x]);
      //find the neighbors of this point and update the list
      findNeighbors(parms,&pt);      
    }
  changeLabelfromList(parms);
  
}

static int componentMerge(TC_PARMS *parms,int s0 ,int s1)
{
  SEGMENTATION *segmentation=parms->segmentation;
  MRI *mri_labeled=parms->mri_labeled;
  ConnectedComponent *cc0,*cc1;
  int  v, total_cells, x, y, z,ncells ;
  CCS*ccs=parms->ccs;


  cc0=&segmentation->components[s0];
  cc1=&segmentation->components[s1];

  removeComponent(ccs,cc0);

  total_cells = cc0->ncells+cc1->ncells ;

  if(total_cells>=cc0->maxcells)
    reallocateComponent(cc0,total_cells+10);
  
  ncells=cc0->ncells;
  for(v = ncells ; v < total_cells ; v++)
  {
    x = cc1->cells[v-ncells]->x ;
    y = cc1->cells[v-ncells]->y ;
    z = cc1->cells[v-ncells]->z ;
    MRISvox(mri_labeled, x, y, z) = s0 ;
    addCellToComponent(cc0,cc1->cells[v-ncells]);
  }
  
  removeComponent(ccs,cc1);
  addComponent(ccs,cc0);
  segmentation->ncomponents--;
  
  cc1->ncells = 0 ;
  cc1->cost=-1.0;
  
  return NO_ERROR;
}

//assumes that the important labels are BODY,RESIDUE
//eventually assignes the ISOLATED label to some voxels
static int segmentBody(TC_PARMS *parms)
{
  MSV msv;
  int s,x,y,z;

  if(!parms->multiplemode)
    while(findPoint(parms,&msv))
    {
      x=msv.x;
      y=msv.y;
      z=msv.z;
      s=componentNew(parms->segmentation);
      addCellToComponent(&parms->segmentation->components[s],&parms->list->table[z][y][x]);    
      MRISvox(parms->mri_labeled,x,y,z)=s;
      parms->current_label=s;
      MRIvox(parms->mri_bin,x,y,z)=BODY;//modify in the CTExpansion function
      //find the neighbors of this point and update the list
      findNeighbors(parms,&msv);
      CTExpansion(parms);
    }
  else
    while(findPoint(parms,&msv))
    {
      x=msv.x;
      y=msv.y;
      z=msv.z;
      segmentVoxel(parms,&msv);
      findNeighbors(parms,&msv);
      CTExpansion(parms);
    }
  return NO_ERROR;
}

static int segmentConnectedComponents(TC_PARMS *parms,Cell **list,int ncells)
{
  SEGMENTATION *segmentation=parms->segmentation;
  MRI *mri=parms->mri_bin,*mri_labeled=parms->mri_labeled;
  ConnectedComponent *cc;
  int sum,x, y, z, width, height, depth, xi, yi, zi, xk, yk, zk
    , border_labels[27], nlabels, label, nvox,m,connectivity,con,max;

  width=mri->width;
  height=mri->height;
  depth=mri->depth;
  
  connectivity=parms->c_c;
  con=3;
  //modification florent 4th March!
  //con=connectivityNumber(connectivity);
  
  for(m=0;m<ncells;m++)
    {
      x=list[m]->x;
      y=list[m]->y;
      z=list[m]->z;
      memset(border_labels, -1, 27*sizeof(int)) ;
      for (nvox = 0, zk = -1 ; zk <= 1 ; zk++)
	{
	  zi = z+zk ; 
	  for (yk = -1 ; yk <= 1 ; yk++)
	    {
	      yi = y+yk ; 
	      for (xk = -1 ; xk <= 1 ; xk++)
		{	
		  sum=abs(xk) + abs(yk) + abs(zk);
		  if (!sum || sum>con)
		    continue ;
		  xi = x+xk ; 
		  if(MRIvox(mri, xi, yi, zi)!=RESIDUE)
		    continue;
		  label = MRISvox(mri_labeled, xi, yi, zi) ;
		  if ((label >= 0) && (!segmentation->components[label].found))
		    {
		      segmentation->components[label].found = 1 ;
		      border_labels[nvox] = label ;
		      nvox++;
		    }
		}
	    }
	}  
      nlabels=nvox;
      for (nvox = 0 ; nvox < nlabels ; nvox++)
	segmentation->components[border_labels[nvox]].found= 0;/* for next time */
      label = 0 ;
      switch (nlabels)
	{
	case 0:          /* allocate a new segment */
	  label=componentNew(segmentation);
	  cc=&segmentation->components[label];
	  break ;
	case 1:          /* assign this voxel to the one that it borders */
	  label = border_labels[0] ;
	  cc=&segmentation->components[label];
	  break ;
	default:         /* merge segments and assign to largest number of cells */
	  max=-1;
	  for(nvox = 0 ; nvox < nlabels ; nvox++) //find largest number
	    if (segmentation->components[border_labels[nvox]].ncells>max)
	      {
		label=border_labels[nvox];
		max=segmentation->components[border_labels[nvox]].ncells;
	      }
	  cc=&segmentation->components[label];
	  for(nvox = 0 ; nvox < nlabels ; nvox++) //merge into lowest number
	    {
	      if(label==border_labels[nvox])
		continue;
	      componentMerge(parms, label, border_labels[nvox]);
	    }
	  break ;
	}
     /* add it to the existing list */
      addCellToComponent(cc,&parms->list->table[z][y][x]);
      MRISvox(mri_labeled, x, y, z) = label ;
    } 

return NO_ERROR;
}

//analyze two RCCs to see if they need to be merged together!
static int analyzeRCCs(TC_PARMS *parms,SEGMENTATION *bccseg,int rcc0,int rcc1)
{
  int x0,y0,z0,x1,y1,z1,ncells0,ncells1,a,b,c;
  Cell **list0,**list1;
  int loop,m,n,l,p,sum;
  
  int con,connectivity=parms->c_c,val;
  MRI *mri_bin=parms->mri_bin,*mri_labeled=parms->mri_labeled;
  int nbrnbh;
  int pos[26][4]; //for x,y,z and label
  int border_labels[10],nlabels;

  //compilator warnings
  x0=y0=z0=x1=y1=z1=0;

  //size of the list 
  if(parms->segmentation->components[rcc0].ncells<parms->segmentation->components[rcc1].ncells)
    {
      list0=parms->segmentation->components[rcc0].cells;
      list1=parms->segmentation->components[rcc1].cells;
      ncells0=parms->segmentation->components[rcc0].ncells;
      ncells1=parms->segmentation->components[rcc1].ncells;
    }
  else
    {
      list0=parms->segmentation->components[rcc1].cells;
      list1=parms->segmentation->components[rcc0].cells;
      ncells0=parms->segmentation->components[rcc1].ncells;
      ncells1=parms->segmentation->components[rcc0].ncells;
    }

  con=connectivityNumber(connectivity);

  for(loop=1,m=0;loop&&m<ncells0;m++)
    {
      x0=list0[m]->x;y0=list0[m]->y;z0=list0[m]->z;
      
      //find the neighbors in BCCs
      //check if the first voxel has at least two BCC neighbors!
      nlabels=0;nbrnbh=0;
      for(a=-1;a<2;a++)
	for(b=-1;b<2;b++)
	  for(c=-1;c<2;c++)
	    {
	      sum=abs(a)+abs(b)+abs(c);
	      if(!sum || sum>con)
		continue;
	      if(MRIvox(mri_bin,x0+a,y0+b,z0+c)!=BODY)
		continue;
	      val=MRISvox(mri_labeled,x0+a,y0+b,z0+c);
	      if(val>=0)
		{
		  pos[nbrnbh][0]=x0+a;
		  pos[nbrnbh][1]=y0+b;
		  pos[nbrnbh][2]=z0+c;
		  pos[nbrnbh][3]=val;
		  nbrnbh++;
		  if(!bccseg->components[val].found)
		    {
		      bccseg->components[val].found=1;
		      border_labels[nlabels]=val;
		      nlabels++;
		    }
		}
	    }
      for(n=0;n<nlabels;n++)
	bccseg->components[border_labels[n]].found=0;
      if(nlabels<2)  //not enough neighbors!
	  continue;
      //at least 2 neighbors -> check in the second RCC rcc1
      for(n=0;loop&&n<ncells1;n++)
	{
	  //first check if the two voxels are N26 neighbors
	  x1=list1[n]->x;
	  if(abs(x1-x0)>1)
	    continue;
	  y1=list1[n]->y;
	  if(abs(y1-y0)>1)
	    continue;
	  z1=list1[n]->z;
	  if(abs(z1-z0)>1)
	    continue;
	  //then check if they have two commons neighbors
	  for(l=0;loop&&l<nbrnbh;l++)
	    {
	      //check if this point is neighbor of x1
	      if((abs(pos[l][0]-x1)+abs(pos[l][1]-y1)+abs(pos[l][2]-z1))>con)
		continue;
	      val=pos[l][3];
	      for(p=l+1;loop&&p<nbrnbh;p++)
		{
		  if((pos[p][3]==val)||
		     ((abs(pos[p][0]-x1)+abs(pos[p][1]-y1)+abs(pos[p][2]-z1))>con))
		    continue;
		  loop=0;
		}
	    }
	} 
    }
  
  if(!loop) //merge rcc0 into rcc1
    {
      componentMerge(parms, rcc1, rcc0);
      return -1;
    }
     
  return 1;
}




//return !=-1 only if an RCC is found with a label smaller than rcc
static int findRCC(TC_PARMS *parms,int rcc)
{
  int label;
  int x,y,z,a,b,c,k,loop,val1,val2;
  MRI *mri_bin=parms->mri_bin,*mri_labeled=parms->mri_labeled;
  ConnectedComponent* component=&parms->segmentation->components[rcc];
  SEGMENTATION *segmentation=parms->segmentation;

  //check for N26
  for(label=-1,loop=1,k=0;loop&&k<component->ncells;k++)
    {
      x=component->cells[k]->x;
      y=component->cells[k]->y;
      z=component->cells[k]->z;

      for(a=-1;loop&&a<2;a++)
	for(b=-1;loop&&b<2;b++)
	  for(c=-1;loop&&c<2;c++)
	    {
	      val1=MRIvox(mri_bin,a+x,b+y,c+z);
	      val2=MRISvox(mri_labeled,a+x,b+y,c+z);
	      if(val1==RESIDUE && val2>=0 
		 && val2<rcc && (!segmentation->components[val2].found))
		{
		  segmentation->components[val2].found=1;
		  loop=0;
		  label=val2;
		}
	    }
    }
  return label;
}

#define MAX_BORDER_LABELS 50
static int mergeRCC(TC_PARMS *parms,SEGMENTATION *segmentation)
{
  int k,nrcccomponents,rcclabel,nlabels,maxlabels,*border_labels,*tmp,m;
  ConnectedComponent *component;
  int count;
 
  nrcccomponents=parms->segmentation->ncomponents;
 
  maxlabels=MAX_BORDER_LABELS;
  border_labels=(int*)malloc(maxlabels*sizeof(int));

  for(count=0,k=nrcccomponents-1;k>=0;k--)
    {
      if(!parms->multiplemode)
	fprintf(stderr,"\r   iteration %5d:",nrcccomponents-k);
      component=&parms->segmentation->components[k];
      //analyze the current rcc with the ones, which indice is smaller
      rcclabel=0;
      nlabels=0;
      while(rcclabel>=0)
	{
	  rcclabel=findRCC(parms,k);
	  if(rcclabel>=0)
	    {
	      border_labels[nlabels]=rcclabel;
	      nlabels++;
	      rcclabel=analyzeRCCs(parms,segmentation,k,rcclabel); //if return -1, means k -> rcclabel
	      if(nlabels==maxlabels)
		{
		  maxlabels=(int)(maxlabels*1.5);
		  tmp=(int*)malloc(maxlabels*sizeof(int));
		  for(m=0;m<nlabels;m++)
		    tmp[m]=border_labels[m];
		  free(border_labels);
		  border_labels=tmp;
		}
	      if(rcclabel==-1)
		count++;
	    }
	}
      for(m=0;m<nlabels;m++)
	parms->segmentation->components[border_labels[m]].found=0;
    }

  free(border_labels);
  if(nrcccomponents && !parms->multiplemode)
    fprintf(stderr," %d components merged together: %d RCC components\n"
	    ,count,parms->segmentation->ncomponents);

  return NO_ERROR;
}



#define CELLNBR 10000
static int computeResidualSegmentation(TC_PARMS *parms)
{
  int i,j,k,n;
  int width,height,depth;
  MRI *mri=parms->mri_bin;
  Cell **table,**tabtmp;
  int cellnbr,count;

  width=mri->width;
  height=mri->height;
  depth=mri->depth;


  //first list all the RESIDUE voxels
  cellnbr=CELLNBR;
  table=(Cell**)calloc(cellnbr,sizeof(Cell*));
  count=0;
  for(k=1;k<depth-1;k++)
    for(j=1;j<height-1;j++)
      for(i=1;i<width-1;i++)
	if(MRIvox(mri,i,j,k)==RESIDUE)
	  {
	    if(count==cellnbr)
	      {
		cellnbr=cellnbr*2;
		tabtmp=(Cell**)calloc(cellnbr,sizeof(Cell*));
		for(n=0;n<count;n++)
		  tabtmp[n]=table[n];
		free(table);
		table=tabtmp;
	      }
	    MRISvox(parms->mri_labeled,i,j,k)=-1;
	    table[count]=&parms->list->table[k][j][i];
	    count++;
	  }
  cellnbr=count;

  //then segment from this list
  segmentConnectedComponents(parms,table,cellnbr);


  //finally check if I should merge some RCCs together!!!
  if(parms->c_c==parms->f_c)
    mergeRCC(parms,parms->F_Bseg);
  else
    mergeRCC(parms,parms->B_Bseg);
  

  return NO_ERROR;
}

#if 0
static void PrintStatistics(SEGMENTATION *segmentation)
{
  int k;
  fprintf(stderr,"\nthe segmentation has %d components",segmentation->ncomponents);
  if((CCS*)segmentation->ccs)
    fprintf(stderr," (et %d components)",((CCS*)segmentation->ccs)->ncomponents);
  for(k=0;k<segmentation->maxcomponents;k++)
    if(segmentation->components[k].cost>=0 && segmentation->components[k].ncells)
      ;//fprintf(stderr,"\n  the cost of the component nbr %d =is %f with %d cells",k,segmentation->components[k].cost,segmentation->components[k].ncells);
}


#endif

static void PrintSurfaceStatistics(MRI *mri,int label,int connectivity)
{
    int euler,pnvertices,  pnfaces, pnedges;
    MRIS *mris;
    mris=MRIScreateSurfaceFromVolume(mri,label,connectivity);
    euler=MRIScomputeEulerNumber(mris,&pnvertices,&pnfaces,&pnedges);
    fprintf(stderr," \neuler number = %d v=%d f=%d e=%d\n",
	  euler,pnvertices,pnfaces,pnedges);
    MRISwrite(mris,"./surface");
    MRISfree(&mris);
}

static int backgroundSegmentation(TC_PARMS *parms)
{
  int s;
  
  //first segmentation of the main connecting component
  s=componentNew(parms->segmentation);
  parms->current_label=s;
  setBorderValue(parms->mri_bin,BODY,1);
  setBorderValue(parms->mri_labeled,s,1);
  setBorderValueToSegmentation(parms,s,1);
  addBorderVoxels(parms,2);
  parms->segmentation->components[s].map=10000000000;//parms->segmentation->components[s].ncells;

  CTExpansion(parms);
  return NO_ERROR;
}

static void SaveOrigMaps(TC_PARMS *parms)
{
  MRI *mri;
  char fname[100];
  int width,height,depth,x,y,z;
  int xinit,yinit,zinit,n,val;

  width=parms->width;
  height=parms->height;
  depth=parms->depth;

  xinit=parms->region.x-2;
  yinit=parms->region.y-2;
  zinit=parms->region.z-2;

  mri=MRIalloc(width,height,depth,MRI_FLOAT);
  for(z=1;z<depth-1;z++)
    for(y=1;y<height-1;y++)
      for(x=1;x<width-1;x++)
	{
	  val=MRIvox(parms->mri_seg,xinit+x,yinit+y,zinit+z);
	  for(n=0;n<parms->nblabels;n++)
	    if(val==parms->labels[n])
	      MRIFvox(mri,x,y,z)=1;
	}

  sprintf(fname,parms->mapsfname);
  strcat(fname,"/out_input.mgh");
  MRIwrite(mri,fname);

  if(parms->mri_orig)
    {
      for(z=1;z<depth-1;z++)
	for(y=1;y<height-1;y++)
	  for(x=1;x<width-1;x++)
	    MRIFvox(mri,x,y,z)=(float)MRIvox(parms->mri_orig,xinit+x,yinit+y,zinit+z);
      sprintf(fname,parms->mapsfname);
      strcat(fname,"/out_orig.mgh");
      MRIwrite(mri,fname);
    }
  MRIfree(&mri);

}

static void SaveInitMaps(TC_PARMS *parms)
{
  MRI *mri,*mri_cost,*mri_seg;
  char fname[100];
  int width,height,depth,k,l,x,y,z;
  SEGMENTATION *seg;

  width=parms->width;
  height=parms->height;
  depth=parms->depth;

  mri_cost=MRIalloc(width,height,depth,MRI_FLOAT);
  mri=MRIalloc(width,height,depth,MRI_FLOAT);
  mri_seg=MRIalloc(width,height,depth,MRI_FLOAT);
  seg=parms->F_Bseg;
  for(k=0;k<seg->maxcomponents;k++)
    if(seg->components[k].cost>=0)
      for(l=0;l<seg->components[k].ncells;l++)
	{
	  x=seg->components[k].cells[l]->x;
	  y=seg->components[k].cells[l]->y;
	  z=seg->components[k].cells[l]->z;
	  MRIFvox(mri,x,y,z)=k+1;
	  MRIFvox(mri_cost,x,y,z)=seg->components[k].map;
	  MRIFvox(mri_seg,x,y,z)=1;
	}
  seg=parms->F_Rseg;
  for(k=0;k<seg->maxcomponents;k++)
    if(seg->components[k].cost>=0)
      for(l=0;l<seg->components[k].ncells;l++)
	{
	  x=seg->components[k].cells[l]->x;
	  y=seg->components[k].cells[l]->y;
	  z=seg->components[k].cells[l]->z;
	  MRIFvox(mri,x,y,z)=-k-1;
	  MRIFvox(mri_cost,x,y,z)=seg->components[k].map;
	  MRIFvox(mri_seg,x,y,z)=3;
	}
  sprintf(fname,parms->mapsfname);
  strcat(fname,"/out_fs.mgh");
  MRIwrite(mri,fname);
  MRIfree(&mri);
  mri=MRIalloc(width,height,depth,MRI_FLOAT);
  seg=parms->B_Bseg;
  for(k=0;k<seg->maxcomponents;k++)
    if(seg->components[k].cost>=0)
      for(l=0;l<seg->components[k].ncells;l++)
	{
	  x=seg->components[k].cells[l]->x;
	  y=seg->components[k].cells[l]->y;
	  z=seg->components[k].cells[l]->z;
	  MRIFvox(mri,x,y,z)=k+1;
	  MRIFvox(mri_cost,x,y,z)=seg->components[k].map;
	}
  seg=parms->B_Rseg;
  for(k=0;k<seg->maxcomponents;k++)
    if(seg->components[k].cost>=0)
      for(l=0;l<seg->components[k].ncells;l++)
	{
	  x=seg->components[k].cells[l]->x;
	  y=seg->components[k].cells[l]->y;
	  z=seg->components[k].cells[l]->z;
	  MRIFvox(mri,x,y,z)=-k-1;
	  MRIFvox(mri_cost,x,y,z)=seg->components[k].map;
	  MRIFvox(mri_seg,x,y,z)=2;
	}
  sprintf(fname,parms->mapsfname);
  strcat(fname,"/out_bs.mgh");
  MRIwrite(mri,fname);
  MRIfree(&mri);

  sprintf(fname,parms->mapsfname);
  strcat(fname,"/out_seg.mgh");
  MRIwrite(mri_seg,fname);
  MRIfree(&mri_seg);

  

  sprintf(fname,parms->mapsfname);
  strcat(fname,"/out_ccost.mgh");
  MRIwrite(mri_cost,fname);
  MRIfree(&mri_cost);

}

static int initSegmentation(TC_PARMS *parms)
{ 
  fprintf(stderr,"\n***********************************************");
  fprintf(stderr,"\nINITIALIZATION OF THE SEGMENTATION");
  //FOREGROUND SEGMENTATION
  //BODY
  fprintf(stderr,"\n   FOREGROUND SEGMENTATION:      BODY\n");
  parms->c_c=parms->f_c;
  mriChangeLabel(parms->mri_bin,F_R,RESIDUE);  
  initCCSSEG(parms);
  segmentBody(parms);
  ccsSortComponents(parms->ccs);
  parms->F_Bseg=parms->segmentation;
  parms->F_Bccs=parms->ccs;

  fprintf(stderr,"\r   %d component(s)                       "
	  ,parms->segmentation->ncomponents);
  //PrintSurfaceStatistics(parms->mri_bin,parms->c_c);
  //PrintStatistics(parms->segmentation);
      
  //RESIDUE
  fprintf(stderr,"\n   FOREGROUND SEGMENTATION:      RESIDUE\n");
  initCCSSEG(parms);
  computeResidualSegmentation(parms);
  ccsSortComponents(parms->ccs);
  parms->F_Rseg=parms->segmentation;
  parms->F_Rccs=parms->ccs;

  fprintf(stderr,"   %d component(s)     "
	  ,parms->segmentation->ncomponents);
  //  PrintStatistics(parms->segmentation);

  mriChangeLabel(parms->mri_bin,BODY,F_B);
  mriChangeLabel(parms->mri_bin,RESIDUE,F_R);
  
  //BACKGROUND SEGMENTATION
  //BODY
  fprintf(stderr,"\n   BACKGROUND SEGMENTATION:      BODY\n");
  parms->c_c=parms->b_c;
  mriChangeLabel(parms->mri_bin,B_R,RESIDUE);
  initCCSSEG(parms);
  backgroundSegmentation(parms);
  segmentBody(parms);
  ccsSortComponents(parms->ccs);
  parms->B_Bseg=parms->segmentation;
  parms->B_Bccs=parms->ccs;
  
  fprintf(stderr,"\r   %d component(s)                      "
	  ,parms->segmentation->ncomponents);  
  //PrintSurfaceStatistics(parms->mri_bin,parms->c_c)
  //PrintStatistics(parms->segmentation);
  
 //RESIDUE
  fprintf(stderr,"\n   FOREGROUND SEGMENTATION:      RESIDUE\n");
  initCCSSEG(parms);
  computeResidualSegmentation(parms);
  ccsSortComponents(parms->ccs);
  parms->B_Rseg=parms->segmentation;
  parms->B_Rccs=parms->ccs;
  
  fprintf(stderr,"   %d component(s)     "
	  ,parms->segmentation->ncomponents);
  //PrintStatistics(parms->segmentation);

  mriChangeLabel(parms->mri_bin,BODY,B_B);
  mriChangeLabel(parms->mri_bin,RESIDUE,B_R);

  if(parms->mapsfname)
    SaveInitMaps(parms);
  
  return NO_ERROR;
}

static int findandRemoveComponents(TC_PARMS *parms)
{
  CCS*ccs;
  ConnectedComponent *cc,*ncc;
  MRI *mri=parms->mri_bin,*mri_labeled=parms->mri_labeled,*mri_prior,*mri_cost,*mri_dist;
  int ref,k,n,x,y,z;
  int label;
  float threshold=parms->threshold;
  SEGMENTATION *seg;
  int check;

  mri_dist=parms->mri_dist;
  
  check=1;
  if(parms->c_c==parms->b_c)  //work on the foreground
  {
    ccs=parms->F_Bccs;  //body
    seg=parms->F_Bseg;
    label=B_R;
    mri_prior=parms->mri_bprior;
    mri_cost=parms->mri_bcost;
  }
  else
  {
    ccs=parms->B_Bccs;  //body
    seg=parms->B_Bseg;
    label=F_R;
    mri_prior=parms->mri_fprior;
    mri_cost=parms->mri_fcost;
  }  
  
  ref=(int)(threshold*ccs->number_of_bins/ccs->maxcost);
  if(ref>=ccs->number_of_bins)  
    ref=ccs->number_of_bins-1;
  for(k=0;k<=ref;k++)
    {
      cc=ccs->sorting[k];
      if(parms->c_c==parms->b_c && k==ccs->reference_max_component)
	cc=cc->next;
      while(cc)
	{
	  ncc=cc->next;
	  if(parms->c_c==parms->f_c && 
	     MRISvox(mri_labeled,cc->cells[0]->x,cc->cells[0]->y,cc->cells[0]->z)==0)
	      cc=ncc;
	  else
	    { 
	      if(cc->cost<=threshold)
		{
		  for(n=0;n<cc->ncells;n++)
		    {
		      x=cc->cells[n]->x;
		      y=cc->cells[n]->y;
		      z=cc->cells[n]->z;
		     
		      cc->cells[n]->cost=MRIFvox(mri_cost,x,y,z);
		      cc->cells[n]->prior=MRIFvox(mri_prior,x,y,z);

		      MRIvox(mri,x,y,z)=label;
		      MRISvox(mri_labeled,x,y,z)=-1;
		      cc->cells[n]->type=UNKNOWN;
		    }
		  removeComponent(ccs,cc);
		  seg->ncomponents--;
		  cc->ncells=0;
		  cc->cost=-1;
		}
	      cc=ncc;
	    }
	}
    }

  if(ccs->ncomponents>1)
    check=0;

  if(parms->c_c==parms->b_c)  //work on the foreground
  {
    ccs=parms->F_Rccs;  //residue
    seg=parms->F_Rseg;
    label=B_R;
    mri_prior=parms->mri_bprior;
    mri_cost=parms->mri_bcost;
  }
  else
  {
    ccs=parms->B_Rccs;  //residue
    seg=parms->B_Rseg;
    label=F_R;
    mri_prior=parms->mri_fprior;
    mri_cost=parms->mri_fcost;
  }  

  ref=(int)(threshold*ccs->number_of_bins/ccs->maxcost);
  if(ref>=ccs->number_of_bins)  
    ref=ccs->number_of_bins-1;
  for(k=0;k<=ref;k++)
    {
      cc=ccs->sorting[k];
      while(cc)
	{
	  ncc=cc->next;
	  if(cc->cost<=threshold)
	    {
	      for(n=0;n<cc->ncells;n++)
		{
		  x=cc->cells[n]->x;
		  y=cc->cells[n]->y;
		  z=cc->cells[n]->z;
		  
		  cc->cells[n]->cost=MRIFvox(mri_cost,x,y,z);
		  cc->cells[n]->prior=MRIFvox(mri_prior,x,y,z);

		  MRIvox(mri,x,y,z)=label;
		  MRISvox(mri_labeled,x,y,z)=-1;
		  cc->cells[n]->type=UNKNOWN;
		}
	      removeComponent(ccs,cc);
	      seg->ncomponents--;
	      cc->cost=-1;
	      cc->ncells=0;
	    }
	  cc=ncc;
	}
    }
 if(ccs->ncomponents)
    check=0;

 seg->ncomponents=ccs->ncomponents;


   if(parms->c_c==parms->b_c)  //work on the foreground
    {
      ccs=parms->B_Rccs;  //residue
      seg=parms->B_Rseg;
      label=B_R;
      mri_prior=parms->mri_bprior;
      mri_cost=parms->mri_bcost;
    }
  else
    {
      ccs=parms->F_Rccs;  //residue
      seg=parms->F_Rseg;
      label=F_R;
      mri_prior=parms->mri_fprior;
      mri_cost=parms->mri_fcost;
    }

  for(k=0;k<ccs->number_of_bins;k++)
    {
      cc=ccs->sorting[k];
      while(cc)
	{
	  ncc=cc->next;
	  for(n=0;n<cc->ncells;n++)
	    {
	      x=cc->cells[n]->x;
	      y=cc->cells[n]->y;
	      z=cc->cells[n]->z;
	      
	      cc->cells[n]->cost=MRIFvox(mri_cost,x,y,z);
	      cc->cells[n]->prior=MRIFvox(mri_prior,x,y,z);
		  
	      MRIvox(mri,x,y,z)=label;
	      MRISvox(mri_labeled,x,y,z)=-1;
	      cc->cells[n]->type=UNKNOWN;
	    }
	  removeComponent(ccs,cc);
	  seg->ncomponents--;
	  cc->cost=-1;
	  cc->ncells=0;
	  cc=ncc;
	}
    }
 
  seg->ncomponents=0;
  for(k=0;k<seg->maxcomponents;k++)
    {
      seg->components[k].cost=0;
      seg->components[k].ncells=0;
    }
  
  return check;
}

static int updateVolume(TC_PARMS *parms)
{
  if(parms->c_c==parms->f_c)
    {
      mriChangeLabel(parms->mri_bin,F_B,BODY);
      mriChangeLabel(parms->mri_bin,F_R,RESIDUE);
      parms->segmentation=parms->F_Bseg;
      parms->ccs=parms->F_Bccs;
    }
  else
    {
      mriChangeLabel(parms->mri_bin,B_B,BODY);
      mriChangeLabel(parms->mri_bin,B_R,RESIDUE);
      parms->segmentation=parms->B_Bseg;
      parms->ccs=parms->B_Bccs;
    }

  return NO_ERROR;
}

static int resetVolume(TC_PARMS *parms)
{

  //  PrintSurfaceStatistics(parms->mri_bin,BODY,parms->c_c);

  if(parms->c_c==parms->f_c)
    {
      mriChangeLabel(parms->mri_bin,BODY,F_B);
      mriChangeLabel(parms->mri_bin,RESIDUE,F_R);
    }
  else
    {
      mriChangeLabel(parms->mri_bin,BODY,B_B);
      mriChangeLabel(parms->mri_bin,RESIDUE,B_R);
    }

  return NO_ERROR;
}

static int updateResidualSegmentation(TC_PARMS *parms)
{
  int i,j,k;
  int width,height,depth;
  MRI *mri=parms->mri_bin,*mri_labeled=parms->mri_labeled;

  width=mri->width;
  height=mri->height;
  depth=mri->depth;

  if(parms->c_c==parms->f_c)
    {
      parms->segmentation=parms->F_Rseg;
      parms->ccs=parms->F_Rccs;
    }
  else
    {
      parms->segmentation=parms->B_Rseg;
      parms->ccs=parms->B_Rccs;
    }

  for(k=0;k<parms->segmentation->maxcomponents;k++)
    {
      parms->segmentation->components[k].ncells=0;
      parms->segmentation->components[k].cost=-1;
    }
  parms->segmentation->ncomponents=0;
  
  
  for(k=1;k<depth-1;k++)
    for(j=1;j<height-1;j++)
      for(i=1;i<width-1;i++)
	if(MRIvox(mri,i,j,k)==RESIDUE)
	  MRISvox(mri_labeled,i,j,k)=-1;

  computeResidualSegmentation(parms);
  ccsSortComponents(parms->ccs);

  return NO_ERROR;
}

#define INITCOST 100000000.0
static float updateThreshold(TC_PARMS *parms,int type)
{
  int ref,k;
  CCS *bccs,*rccs;  
  ConnectedComponent *cc;
  float threshold=parms->threshold,cost,mincost,th1,th2;

  if(type==parms->b_c)   //work on the foreground ->priority background
    {
      bccs=parms->F_Bccs;
      rccs=parms->F_Rccs;
    }
  else
    {
      bccs=parms->B_Bccs;
      rccs=parms->B_Rccs;
    }

  //body
  ref=(int)(threshold*bccs->number_of_bins/bccs->maxcost);
  if(ref>=bccs->number_of_bins)  
    ref=bccs->number_of_bins-1;
  mincost=INITCOST;
  for(k=0;k<=ref;k++)
    {
      cc=bccs->sorting[k];
      while(cc)
	{
	  cost=cc->cost;
	  if(cost>threshold && cost<mincost)
	    {
	      mincost=cost;
	    }
	  if(cost<=threshold)
	    break;
	  cc=cc->next;
	}
    }
  if(mincost==INITCOST)
    mincost=threshold+THRESHOLD_INCREASE;
  th1=mincost;
  //residue
  ref=(int)(threshold*rccs->number_of_bins/rccs->maxcost);
  if(ref>=rccs->number_of_bins)  
    ref=rccs->number_of_bins-1;
  mincost=INITCOST;
  for(k=0;k<=ref;k++)
    {
      cc=rccs->sorting[k];
      while(cc)
	{
	  cost=cc->cost;
	  if(cost>threshold && cost<mincost)
	    {
	      mincost=cost;
	    }
	  if(cost<=threshold)
	    break;
	  cc=cc->next;
	}
    }
  if(mincost==INITCOST)
    mincost=threshold+THRESHOLD_INCREASE;
  th2=mincost;
  
  return MIN(th1,th2);
}

static int correctSegmentation(TC_PARMS *parms)
{
  float threshold,th1,th2;
  fprintf(stderr,"\n****************************************************");
  fprintf(stderr,"\nCORRECTION OF THE SEGMENTATION");

  computeMap(parms);
  //initialization step
  parms->multiplemode=1;
  if(parms->mode==1)   //priority to the foreground !
    {
      parms->c_c=parms->f_c;
      parms->threshold=updateThreshold(parms,parms->c_c);
    }
  else if(parms->mode==2)               //priority to the background !
    {    
      parms->c_c=parms->b_c;
      parms->threshold=updateThreshold(parms,parms->c_c);
    }
  else
    {
      th1=updateThreshold(parms,parms->f_c); //foreground correction -> priority background
      th2=updateThreshold(parms,parms->b_c); //background
      threshold=MIN(th1,th2);
      if(threshold==th1)
	parms->c_c=parms->f_c;
      else
	parms->c_c=parms->b_c;
      parms->threshold=threshold;
    }

  fprintf(stderr,"\nThreshold Mode  F/B   F/R   B/B   B/R  MAP\n"); 
  while(parms->F_Rseg->ncomponents || parms->B_Rseg->ncomponents)
  {
    fprintf(stderr,"\r");
    fprintf(stderr,"%5.5f   ",parms->threshold);
    if(parms->c_c==parms->b_c)
      fprintf(stderr,"bgd ");
    else
      fprintf(stderr,"fgd ");

    //find the components smaller than threshold and remove them
    if(findandRemoveComponents(parms))   //we stop the segmentation process!
      {
	fprintf(stderr,"    1     0     1     0  ");
	if(parms->c_c==parms->b_c)
	  mriChangeLabel(parms->mri_bin,B_R,B_B);
	else
	  mriChangeLabel(parms->mri_bin,F_R,F_B);

	if(parms->priors)
	  fprintf(stderr," %5.5f    (forced exit)          ",computeMap(parms));
	else
	  fprintf(stderr," (forced exit)       ");

	  break;
      }
     
    //update the volume
    updateVolume(parms);

    //expansion step for the complementary volume
    segmentBody(parms);

    //PrintSurfaceStatistics(parms->mri_bin,F_R,parms->c_c);     

    //update the residual segmentation
    updateResidualSegmentation(parms);

    //back to orig volume
    resetVolume(parms);
    
    fprintf(stderr," %4d  %4d  %4d  %4d  "
	    ,parms->F_Bseg->ncomponents,parms->F_Rseg->ncomponents
	    ,parms->B_Bseg->ncomponents,parms->B_Rseg->ncomponents);

    if(parms->priors)
      fprintf(stderr," %5.5f              ",computeMap(parms));
    else
      fprintf(stderr,"               ");

    if(parms->only)
      parms->threshold=updateThreshold(parms,parms->c_c);
    else 
      {
	th1=updateThreshold(parms,parms->c_c); 
	th2=updateThreshold(parms,associatedConnectivity(parms->c_c)); //background
	threshold=MIN(th1,th2);
	if(threshold==th2)  //change of connectivity
	  parms->c_c=associatedConnectivity(parms->c_c);
	parms->threshold=threshold;
	
	//parms->c_c=associatedConnectivity(parms->c_c);
	//if((parms->mode==1 && parms->c_c==parms->b_c)
	   // || (parms->mode==0 && parms->c_c==parms->f_c))
	  //parms->threshold+=THRESHOLD_INCREASE;//updateThreshold(parms)
      }
  }
  fprintf(stderr,"\n");
  return NO_ERROR;
}

//at this point the volume F_R is topologically correct
//Just need to add point without changing the topology
static int finalConditionalExpansion(TC_PARMS *parms)
{  
  int width,height,depth,n,label,ref,x,y,z,xinit,yinit,zinit,i,j,k,val;
  int nlabels=parms->nblabels,*tab=parms->labels,test;
  ConnectedComponent *cc;
  MRI* mri_seg=parms->mri_seg,*mri=parms->mri_bin;
  MSV pt;
  List *list=parms->list;
  Cell ***table=parms->list->table;

  fprintf(stderr,"\n****************************************************");
  fprintf(stderr,"\nFINAL TOPOLOGICAL EXPANSION");

  parms->multiplemode=0;
  label=parms->labels[0];
  xinit=parms->region.x;
  yinit=parms->region.y;
  zinit=parms->region.z;
  
  width=mri->width;
  height=mri->height;
  depth=mri->depth;

  //foreground
  fprintf(stderr,"\n   FOREGROUND\n");
  mriChangeLabel(mri,F_B,BODY);
  parms->c_c=parms->f_c;
  ref=parms->F_Bccs->reference_max_component;
  cc=parms->F_Bccs->sorting[ref];
  x=cc->cells[0]->x;y=cc->cells[0]->y;z=cc->cells[0]->z;
  parms->current_label=MRISvox(parms->mri_labeled,x,y,z);
 
  for(k=0;k<depth-4;k++)
    for(j=0;j<height-4;j++)
      for(i=0;i<width-4;i++)
	{
	  val=MRIvox(mri_seg,xinit+i,yinit+j,zinit+k);
	  for(n=0;n<nlabels;n++)
	    if(val==tab[n] && 
	       MRIvox(mri,i+2,j+2,k+2)!=BODY)
	      {
		MRIvox(mri,i+2,j+2,k+2)=RESIDUE;
		list->table[k+2][j+2][i+2].type=UNKNOWN;
		list->table[k+2][j+2][i+2].prior=table[k+2][j+2][i+2].prior;
	      }
	}
  for(k=1;k<depth;k++)
    for(j=1;j<height;j++)
      for(i=1;i<width;i++)
	if(MRIvox(mri,i,j,k)==BODY && checkNbh(mri,i,j,k,RESIDUE,parms->c_c))
	  {
	    pt.x=i;pt.y=j;pt.z=k;
	    findNeighbors(parms,&pt);
	  }
  CTExpansion(parms);
 
  mriChangeLabel(mri,BODY,F_B);
  mriChangeLabel(mri,RESIDUE,B_B);
  //background
  fprintf(stderr,"\r   BACKGROUND                   \n");
  mriChangeLabel(mri,B_B,BODY);
  parms->c_c=parms->b_c;
  parms->current_label=0;
  
  for(k=0;k<depth-4;k++)
    for(j=0;j<height-4;j++)
      for(i=0;i<width-4;i++)
	{
	  test=1;
	  val=MRIvox(mri_seg,xinit+i,yinit+j,zinit+k);
	  for(n=0;n<nlabels;n++)
	    if(val==tab[n])
	      test=0;
	  if(test && MRIvox(mri,i+2,j+2,k+2)!=BODY)
	    {
	      MRIvox(mri,i+2,j+2,k+2)=RESIDUE;
	      list->table[k+2][j+2][i+2].type=UNKNOWN;
	      list->table[k+2][j+2][i+2].prior=table[k+2][j+2][i+2].prior;
	    }
	}
 
  for(k=1;k<depth;k++)
    for(j=1;j<height;j++)
      for(i=1;i<width;i++)
	if(MRIvox(mri,i,j,k)==BODY && checkNbh(mri,i,j,k,RESIDUE,parms->c_c))
	  {
	    pt.x=i;pt.y=j;pt.z=k;
	    findNeighbors(parms,&pt);
	  }
  CTExpansion(parms);

  mriChangeLabel(mri,BODY,B_B);
  mriChangeLabel(mri,RESIDUE,F_B);
  
  return NO_ERROR;
}

int MRIcorrectTopology(MRI *mri_orig,MRI *mri_seg,MRI **mri_output
			,MRIS **mris,int *labeltable,int nlabels,int connectivity,TC_PARMS *parms)
{
  int freeparms=0,allocmrioutput=0,n=0;
  
  if(!parms)
    {
      parms=initTC_PARMS();
      freeparms=1;
    }

  parms->f_c=connectivity;
  if(parms->tess==-1)
    parms->tess=parms->f_c;
  parms->b_c=associatedConnectivity(parms->f_c);
  parms->mri_orig=mri_orig;
  parms->mri_seg=mri_seg;

  if(mri_output)
    parms->mri_output=*mri_output;
  else
    {
      parms->mri_output=NULL;
      allocmrioutput=1;
    }

  if(mris)
    parms->mris=*mris;

  for(n=0;n<nlabels;n++)
    parms->labels[n]=labeltable[n];
  parms->nblabels=nlabels;

  initImages(parms);
  
#if 0
  (*mris)=MRIScreateSurfaceFromVolume(parms->mri_seg,parms->labels[0],parms->tess);
  
  MRISsmoothSurface2(*mris,2,0.5,0);
  MRISsmoothSurface2(*mris,5,0.15,1);
  MRISsmoothSurface2(*mris,10,0.05,5);
  MRISwrite(*mris,"./testestests");
#endif

  if(parms->priors)
    computeMap(parms);

  initSegmentation(parms);

  correctSegmentation(parms);

  if(parms->priors)
    computeMap(parms);

  finalConditionalExpansion(parms);

  computeNLabels(parms);

  if(parms->priors)
    computeMap(parms);

  if(mri_output)
    modifyImage(parms);

  fprintf(stderr,"\r****************************************************\n");
  if(mris)
    {
      int euler,pnvertices,  pnfaces, pnedges;
      fprintf(stderr,"\nComputing the topologically correct surface");
      if(allocmrioutput)
	modifyImage(parms);
      
      (*mris)=MRIScreateSurfaceFromVolume(parms->mri_output,1,parms->tess);
      
#if 0  
      MRISsmoothSurface2(*mris,2,0.5,0);
      MRISsmoothSurface2(*mris,5,0.15,1);
      MRISsmoothSurface2(*mris,10,0.05,5);
      MRISwrite(*mris,"./testestests1");
#endif

      euler=MRIScomputeEulerNumber(*mris,&pnvertices,&pnfaces,&pnedges);
      fprintf(stderr,"\neuler number = %d v=%d f=%d e=%d\n",
	      euler,pnvertices,pnfaces,pnedges);
      if(allocmrioutput)
	MRIfree(&parms->mri_output);
    }

  if(freeparms)
    TC_PARMSfree(&parms);

  return NO_ERROR;
}

int main(int argc, char *argv[])
{
  TC_PARMS* parms;
  MRIS **mris;
  char  *in_orig_fname, *in_seg_fname,*out_fname;
  int nargs,n;
  char fname[100];
  
  Progname=argv[0];
  parms=initTC_PARMS();
  fprintf(stderr,"\n");
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv,parms) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if(argc<4)
    {
      fprintf(stderr, "\nUsage: %s options input_orig_file input_segmented_file output_folder\n", Progname);   
      exit(1);
    };

  in_orig_fname=argv[argc-3];
  in_seg_fname = argv[argc-2];
  out_fname = argv[argc-1];

  fprintf(stderr,"************************************************************"
          "\nThe input orig volume is %s"
	  "\nThe input segmented volume is %s"
          "\nThe output folder is %s"
          "\nIf this is incorrect, please exit quickly the program (Ctl-C)\n",in_orig_fname,in_seg_fname,out_fname);
  for(n=0;n<parms->nlabels;n++)
    fprintf(stderr,"label = %d: %s \n",parms->labels[n],cma_label_to_name(parms->labels[n]));
  if(parms->priors)
    fprintf(stderr,"mixing parameters: alpha=%1.3f , beta=%1.3f \n",parms->alpha,parms->beta);
  else
    {
      parms->beta=1.0f;
      parms->alpha=1.0f;
    }
  fprintf(stderr,"connectivity = %d\n",parms->f_c);

  if(parms->tess==-1)
    parms->tess=parms->f_c;

  parms->mri_orig=MRIread(in_orig_fname);
  if (!parms->mri_orig && parms->priors)
    Error("orig volume: read failed\n");
  parms->mri_seg=MRIread(in_seg_fname);
  if (!parms->mri_seg)
    Error("segmented volume: read failed\n");

  
  //  if(parms->initsurface)
  fprintf(stderr,"\n doing it...");
  PrintSurfaceStatistics(parms->mri_seg,parms->labels[0],parms->tess);
  fprintf(stderr,"done\n");

  if(parms->surfname)
    mris=&parms->mris;
  else
    mris=NULL;

  MRIcorrectTopology(parms->mri_orig,parms->mri_seg,&parms->mri_output,mris
		     ,parms->labels,parms->nblabels,parms->f_c,parms);

  MRIwrite(parms->mri_output,out_fname);

#if 0
  //validation of the algo
  {
    FILE *f;
    MRIS *mristb[20],*mrisr;
    int n,i,j,k,depth,height,width,count,count2;
    int tab[20]={4,43,51,12,52,13,54,18,53,17,49,10,50,11};//,6,7,10,11,12,13,17,18,43,44,45,46,49,50,51,52,53,54};
    MRI *mri_val=MRIclone(parms->mri_seg,NULL);
    parms->nlabels=1;

    depth=parms->mri_seg->depth;
    height=parms->mri_seg->height;
    width=parms->mri_seg->width;
    for(n=0;n<14;n++)
      {
	MRIfree(&parms->mri_output);
	MRIfree(&parms->mri_bin);
	MRIfree(&parms->mri_dist);
	MRIfree(&parms->mri_fcost);
	MRIfree(&parms->mri_bcost);
	MRIfree(&parms->mri_fprior);
	MRIfree(&parms->mri_bprior);
	MRIfree(&parms->mri_labeled);
	segmentationFree(&parms->F_Bseg);
	segmentationFree(&parms->F_Rseg);
	segmentationFree(&parms->B_Bseg);
	segmentationFree(&parms->B_Rseg);
	CCSfree(&parms->F_Bccs);
	CCSfree(&parms->F_Rccs);
	CCSfree(&parms->B_Bccs);
	CCSfree(&parms->B_Rccs);

	parms->labels[0]=tab[n];
	MRIcorrectTopology(parms->mri_orig,parms->mri_seg,&parms->mri_output,mris
		     ,parms->labels,parms->nblabels,parms->f_c,parms);
	
	

	MRISwrite(*mris,"./tmp");
	mristb[n]=MRISread("./tmp");
#if 0
	count=0;count2=0;
	for(k=0;k<depth;k++)
	  for(j=0;j<height;j++)
	    for(i=0;i<width;i++)
	      {
		if(MRIvox(parms->mri_seg,i,j,k)==parms->labels[0])
		  count2++;
		if(MRIvox(parms->mri_output,i,j,k)==1)
		  {
		    MRIvox(mri_val,i,j,k)++;
		    if(MRIvox(parms->mri_seg,i,j,k)!=parms->labels[0])
		      count++;
		  }
		else if(MRIvox(parms->mri_seg,i,j,k)==parms->labels[0])
		      count++;
	      }
		fprintf(stderr,"\n yeh %d %d %f \n",count,count2,100.*count/count2);
	sprintf(fname,"./label%d",tab[n]);
	f=fopen(fname,"a+");
	fprintf(f,"\n %d %d %f ",count,count2,(float)100.*count/count2);
	fclose(f);
#endif

#if 0
	sprintf(fname,"./surf%d",n);
	MRISwrite(mristb[n],fname);
	MRISsmoothSurface2(mristb[n],5,0.5,0);
	MRISsmoothSurface2(mristb[n],5,0.25,2);
	MRISsmoothSurface2(mristb[n],10,0.05,5);
	sprintf(fname,"./surfsmooth%d",n);
	mristb[n]->type=MRIS_TRIANGULAR_SURFACE;//MRIS_BINARY_QUADRANGLE_FILE;
	MRISwrite(mristb[n],fname);

	MRISsetNeighborhoodSize(mristb[n],3) ;
	MRIScomputeMetricProperties(mristb[n]) ;
	MRIScomputeSecondFundamentalForm(mristb[n]) ;
	MRISuseMeanCurvature(mristb[n]);
	MRISaverageCurvatures(mristb[n],2) ;
	MRISnormalizeCurvature(mristb[n]) ;
	sprintf(fname,"./curv%d",n);
	MRISwriteCurvature(mristb[n],fname);
#endif
      }

#if 0
    mrisr=MRISconcatenateQuadSurfaces(n,mristb);
    mrisr->type=MRIS_TRIANGULAR_SURFACE;
    MRISwrite(mrisr,"./lh.ZURFACE");
    

    //    for(k=0;k<mrisr->nvertices;k++)
    // mrisr->vertices[k].curv=0.3;

    //MRISnormalizeCurvature(mrisr) ;
    MRISwriteCurvature(mrisr,"./ZURFACE_CURVATURE");
    for(k=0;k<mrisr->nvertices;k++)
      mrisr->vertices[k].curv=mrisr->vertices[k].val;
    MRISwriteCurvature(mrisr,"./ZURFACE_VAL");
#endif

    n=0;count=0;
    for(k=0;k<depth;k++)
      for(j=0;j<height;j++)
	for(i=0;i<width;i++)
	  {
	    if(MRIvox(mri_val,i,j,k)>=1)
	      {
		n++;
		if(MRIvox(mri_val,i,j,k)>1)
		  count++;
	      }
	  }
    //    sprintf(fname,"./labeltotal");
    /// f=fopen(fname,"a+");
    //fprintf(f,"\n %s %d %d %f ",in_seg_fname,count,n,(float)100.*count/n);
    //fclose(f);




#if 0
    MRIwrite(mri_val,"/tmp/tmp");
#endif

    fprintf(stderr,"\n WE HAVE %d %d %f   \n",count,n,100.*count/n);
    
  }
#endif

  if(parms->surfname)
    {    
      sprintf(fname,parms->surfname);
      MRISwrite(parms->mris,fname);

      MRISsmoothSurface(parms->mris,7,0.2);
      strcat(fname,"_smooth");
      MRISwrite(parms->mris,fname);

      if(parms->fit)
	{
	  sprintf(fname,parms->surfname);
	  strcat(fname,"_fit");
	  MRISmatchSurfaceToLabel(parms->mris,parms->mri_output,1,NULL,NULL,parms->f_c);
	  MRISwrite(parms->mris,fname);
	}
    }

  TC_PARMSfree(&parms);

  fprintf(stderr,"\n");
  return NO_ERROR;
}










