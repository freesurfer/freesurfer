#include "mri_topology.h"

#define TMP    100
#define MAX_COMP 10 //maximum number of components in a Nbh

//////////////////////////////////////////////////////////////////////
//                  CONNECTIVITY / TOPOLOGICAL NUMBERS 
//
//     TOPOLOGICAL CONVENTION
//     1:            (6+,18)
//     2:            (18,6+)
//     3:            (6,26)
//     4:            (26,6)   
//
////////////////////////////////////////////////////////////////////

int connectivityNumber(int connectivity)
{
  int con;
  switch(connectivity)
    {
    case 1:               //T6+(x,X)=#C6(N_6_3(x,X))
      con=1;

      break;
    case 2:              //T6+(x,X)=#C6(N_18_2)
      con=2;

      break;
    case 3:
      con=1;
      break;
    case 4:
      con=3;
      break;
    default:
      con=1;
      break;
    }
  return con;
}

int associatedConnectivity(int connectivity)
{
  switch(connectivity)
    {
    case 1:               
      return 2;
      break;
    case 2:    
      return 1;
      break;
    case 3:
      return 4;
      break;
    case 4:
      return 3;
      break;
    default:
      return 1;
      break;
    }
}

Nbh* loadNbh(MRI* mri,Nbh* nbh_dst,int i,int j, int k,int label)
{
  int a,b,c;
  Nbh *nbh;

  if(!nbh_dst)
    nbh=(Nbh*)malloc(sizeof(Nbh));
  else
    nbh=nbh_dst;

  for(a=-1;a<2;a++)
    for(b=-1;b<2;b++)
      for(c=-1;c<2;c++)
	if(MRIvox(mri,a+i,b+j,c+k)==label)
	  (*nbh)[a+1][b+1][c+1]=1;
	else
	  (*nbh)[a+1][b+1][c+1]=0;

  return nbh;
}

Nbh* loadSNbh(MRI* mri,Nbh* nbh_dst,int i,int j, int k,int label)
{
  int a,b,c;
  Nbh *nbh;

  if(!nbh_dst)
    nbh=(Nbh*)malloc(sizeof(Nbh));
  else
    nbh=nbh_dst;

  for(a=-1;a<2;a++)
    for(b=-1;b<2;b++)
      for(c=-1;c<2;c++)
	if(MRISvox(mri,a+i,b+j,c+k)==label)
	  (*nbh)[a+1][b+1][c+1]=1;
	else
	  (*nbh)[a+1][b+1][c+1]=0;

  return nbh;
}

Nbh *N_6_1(Nbh* nbh_src,Nbh* nbh_dst)
{
  int i,j,k;
  Nbh* nbh;
  if(!nbh_dst)
    {
      nbh=(Nbh*)calloc(1,sizeof(Nbh));
      if((*nbh_src)[0][1][1])
	(*nbh)[0][1][1]=1;
      if((*nbh_src)[2][1][1])
	(*nbh)[2][1][1]=1;
      if((*nbh_src)[1][0][1])
	(*nbh)[1][0][1]=1;
      if((*nbh_src)[1][2][1])
	(*nbh)[1][2][1]=1;
      if((*nbh_src)[1][1][0])
	(*nbh)[1][1][0]=1;
      if((*nbh_src)[1][1][2])
	(*nbh)[1][1][2]=1;
      return nbh;
    }
  else
    nbh=nbh_dst;

  if((*nbh_src)[0][1][1])
    (*nbh)[0][1][1]=TMP;
  if((*nbh_src)[2][1][1])
    (*nbh)[2][1][1]=TMP;
  if((*nbh_src)[1][0][1])
    (*nbh)[1][0][1]=TMP;
  if((*nbh_src)[1][2][1])
    (*nbh)[1][2][1]=TMP;
  if((*nbh_src)[1][1][0])
    (*nbh)[1][1][0]=TMP;
  if((*nbh_src)[1][1][2])
    (*nbh)[1][1][2]=TMP;
 
  for(k=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
	if((*nbh)[i][j][k]==TMP)
	  (*nbh)[i][j][k]=1;
	else
	  (*nbh)[i][j][k]=0;

  return nbh;
}

Nbh* N_6_2(Nbh* nbh_src,Nbh* nbh_dst)
{
  int i,j,k;
  Nbh* nbh;
  
  nbh=N_6_1(nbh_src,NULL);
  
  for(i=0;i<3;i=i+2)
    if((*nbh)[i][1][1])
      {
	if((*nbh_src)[i][0][1])
	  (*nbh)[i][0][1]=1;
	if((*nbh_src)[i][2][1])
	  (*nbh)[i][2][1]=1;
	if((*nbh_src)[i][1][0])
	  (*nbh)[i][1][0]=1;
	if((*nbh_src)[i][1][2])
	  (*nbh)[i][1][2]=1;
      }

  for(j=0;j<3;j=j+2)
    if((*nbh)[1][j][1])
      {
	if((*nbh_src)[0][j][1])
	  (*nbh)[0][j][1]=1;
	if((*nbh_src)[2][j][1])
	  (*nbh)[2][j][1]=1;
	if((*nbh_src)[1][j][0])
	  (*nbh)[1][j][0]=1;
	if((*nbh_src)[1][j][2])
	  (*nbh)[1][j][2]=1;
      }

  for(k=0;k<3;k=k+2)
    if((*nbh)[1][1][k])
      {
	if((*nbh_src)[0][1][k])
	  (*nbh)[0][1][k]=1;
	if((*nbh_src)[2][1][k])
	  (*nbh)[2][1][k]=1;
	if((*nbh_src)[1][0][k])
	  (*nbh)[1][0][k]=1;
	if((*nbh_src)[1][2][k])
	  (*nbh)[1][2][k]=1;
      }

  if(nbh_dst)
  {
    for(k=0;k<3;k++)
      for(j=0;j<3;j++)
	for(i=0;i<3;i++)
	  (*nbh_dst)[i][j][k]=(*nbh)[i][j][k];
    free(nbh);
    nbh=nbh_dst;
  }

  return nbh;
}

Nbh* N_6_3(Nbh* nbh_src,Nbh* nbh_dst)
{
  int i,j,k;
  Nbh* nbh;

  nbh=N_6_2(nbh_src,NULL);

  i=0;j=0;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;
  i=0;j=0;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;    
  i=0;j=2;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;    
  i=0;j=2;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;  
  i=2;j=0;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;  
  i=2;j=0;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;  
  i=2;j=2;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;  
  i=2;j=2;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1])
      (*nbh)[i][j][k]=1;   

  if(nbh_dst)
  {
    for(k=0;k<3;k++)
      for(j=0;j<3;j++)
	for(i=0;i<3;i++)
	  (*nbh_dst)[i][j][k]=(*nbh)[i][j][k];
    free(nbh);
    nbh=nbh_dst;
  }

  return nbh;
}

Nbh *N_18_1(Nbh* nbh_src,Nbh* nbh_dst)
{
  int i,j,k;
  Nbh* nbh;

  if(!nbh_dst)
    nbh=(Nbh*)calloc(1,sizeof(Nbh));
  else
    nbh=nbh_dst;

  (*nbh)[0][0][0]=TMP;(*nbh)[2][0][0]=TMP;
  (*nbh)[0][0][2]=TMP;(*nbh)[2][0][2]=TMP;
  (*nbh)[0][2][0]=TMP;(*nbh)[2][2][0]=TMP;
  (*nbh)[0][2][2]=TMP;(*nbh)[2][2][2]=TMP;
  (*nbh)[1][1][1]=TMP;

  for(k=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
	if((*nbh)[i][j][k]!=TMP)
	  {
	    if((*nbh_src)[i][j][k])
	      (*nbh)[i][j][k]=1;
	    else
	      (*nbh)[i][j][k]=0;
	  }
	else
	  (*nbh)[i][j][k]=0;

  return nbh;
}

Nbh* N_18_2(Nbh* nbh_src,Nbh* nbh_dst)
{
  int i,j,k;
  Nbh* nbh;
  
  nbh=N_18_1(nbh_src,NULL);

  i=0;j=0;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;
  i=0;j=0;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;  
  i=0;j=2;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;
  i=0;j=2;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;
  i=2;j=0;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;
  i=2;j=0;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;
  i=2;j=2;k=0;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;
  i=2;j=2;k=2;
  if((*nbh_src)[i][j][k])
    if((*nbh)[1][j][k]||(*nbh)[i][1][k]||(*nbh)[i][j][1]
       ||(*nbh)[1][1][k]||(*nbh)[1][j][1]||(*nbh)[i][1][1])
      (*nbh)[i][j][k]=1;

  if(nbh_dst)
  {
    for(k=0;k<3;k++)
      for(j=0;j<3;j++)
	for(i=0;i<3;i++)
	  (*nbh_dst)[i][j][k]=(*nbh)[i][j][k];
    free(nbh);
    nbh=nbh_dst;
  }

  return nbh;
}

Nbh *N_26_1(Nbh* nbh_src,Nbh* nbh_dst)
{
  int i,j,k;
  Nbh* nbh;

  if(!nbh_dst)
    nbh=(Nbh*)calloc(1,sizeof(Nbh));
  else
    nbh=nbh_dst;

  for(k=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
	if((*nbh_src)[i][j][k])
	  (*nbh)[i][j][k]=1;
	else
	  (*nbh)[i][j][k]=0;
  
  (*nbh)[1][1][1]=0;

  return nbh;
}

Nbh* Nnk(Nbh* nbh_src,Nbh *nbh_dst,int connectivity)
{
  Nbh* nbh;

  if(!nbh_dst)
    nbh=(Nbh*)calloc(1,sizeof(Nbh));
  else
    nbh=nbh_dst;

  switch(connectivity)
  {
    case 1:               //T6+(x,X)=#C6(N_6_3(x,X))
      N_6_3(nbh_src,nbh);
      break;
    case 2:              //T6+(x,X)=#C18(N_18_2)
      N_18_2(nbh_src,nbh);
      break;
    case 3:             //T6(x,X)=#C6(N_6_2) 
      N_6_2(nbh_src,nbh);
      break;
    case 4:             //T26(x,X)=#C6(N_26_1) 
      N_26_1(nbh_src,nbh);
      break;
    default:                 //T6+(x,X)=#C6(N_6_3(x,X))
      N_6_3(nbh_src,nbh);
      break;
    }

  return nbh;
}

int checkNbh(MRI *mri,int i,int j,int k,int label,int connectivity)
{
  int a,b,c;
  int con,sum;
  con=connectivityNumber(connectivity);
 
  for(a=-1;a<2;a++)
    for(b=-1;b<2;b++)
      for(c=-1;c<2;c++)
	{
	  sum=abs(a)+abs(b)+abs(c);
	  if(sum>con || (!sum))
	    continue;
	  if(MRIvox(mri,i+a,j+b,k+c)==label)
	    return 1;
	}

  return 0;
}

int checkSNbh(MRI *mri,int i,int j,int k,int label,int connectivity)
{
  int a,b,c;
  int con,sum;
  con=connectivityNumber(connectivity);
 
  for(a=-1;a<2;a++)
    for(b=-1;b<2;b++)
      for(c=-1;c<2;c++)
	{
	  sum=abs(a)+abs(b)+abs(c);
	  if(sum>con || (!sum))
	    continue;
	  if(MRISvox(mri,i+a,j+b,k+c)==label)
	    return 1;
	}

  return 0;
}

//compute the topological number associated with a Nbh and a certain connectivity
int checkTn(Nbh *nbh_src,Nbh *nbh_dst,int connectivity)
{
  int i,j,k,a,b,c,ik,jk,kk,ct;
  int con,nvox,label,sum;
  int comp_table[MAX_COMP];
  int min_val;
  int x,y,z;

  Nbh* nbh;

  if(!nbh_dst)
    nbh=(Nbh*)calloc(1,sizeof(Nbh));
  else
    nbh=nbh_dst;

  con=connectivityNumber(connectivity);
  Nnk(nbh_src,nbh,connectivity);
  
  memset(comp_table,0,MAX_COMP*sizeof(int));

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	if((*nbh)[i][j][k])
	  {
	    for(nvox=0,ik=-1;ik<=1;ik++)
	      {
		a=i+ik;
		if(a<0 || a>=3)
		  continue;
		for(jk=-1;jk<=1;jk++)
		  {
		    b=j+jk;
		    if(b<0 || b>=3)
		      continue;
		    for(kk=-1;kk<=1;kk++)
		      {
			sum=abs(ik)+abs(jk)+abs(kk);
			if(sum>con || (!sum))
			  continue;
			c=k+kk;
			if(c<0 || c>=3)
			  continue;
			label=(*nbh)[a][b][c];
			if(label>1)
			  {
			    comp_table[label-1]=2;
			    nvox++;
			  }
		      }
		  }
	      }
	    if(!nvox)  //find new basin!
	      {
		for(ct=1;comp_table[ct] && ct<MAX_COMP;ct++)
		  ;
		(*nbh)[i][j][k]=ct+1;  //label the new basin
		comp_table[ct]=1;   //note that this number is taken
	      }
	    else
	      {
		min_val=MAX_COMP+1;

		//merging into the smallest value
		for(ct=1;ct<MAX_COMP;ct++)
		  if(comp_table[ct]==2)
		    {
		      min_val=ct;
		      break;
		    }

		(*nbh)[i][j][k]=min_val+1;
		comp_table[min_val]=1;
 
		//merging of the other neighboring values into the smallest one
		for(ct=min_val+1;ct<MAX_COMP;ct++)
		  if(comp_table[ct]==2)
		    {
		      for(x=0;x<3;x++)
			for(y=0;y<3;y++)
			  for(z=0;z<3;z++)
			    if((*nbh)[x][y][z]==ct+1)
			      (*nbh)[x][y][z]=min_val+1;		      
		      //specify that this basin nbr 
		      //doesn't exist anymore
		      comp_table[ct]=0;
		    }	
	      }
	  }

  for(nvox=0,ct=1;ct<MAX_COMP;ct++)
    if(comp_table[ct])
      nvox++;

  if(!nbh_dst)
    free(nbh);

  return nvox;
}

int checkSP(Nbh *fgnbh_src,Nbh *fgnbh_dst,int *fgtp,Nbh *bgnbh_src,Nbh *bgnbh_dst,int *bgtp,int connectivity)
{
  (*fgtp)=checkTn(fgnbh_src,fgnbh_dst,connectivity);
  (*bgtp)=checkTn(bgnbh_src,bgnbh_dst,associatedConnectivity(connectivity));

  if(((*fgtp)==1) && (*bgtp)==1)
    return 1;
  else
    return 0;
}


//Check if two points are strongly X-connected with respect to Y
//for the Strong Connectivity, the two points are assumed to be inside X, Y-adjacent
// and x1 is insice N26(x1)
int checkSC(MRI* mri,int i0, int j0, int k0,int i1,int j1,int k1,int inside_label,int outside_label,int connectivity)
{
  Nbh nbh1,nbh2;
  int u,v,w,a,b,c,loop;

  u=i1-i0;
  v=j1-j0;
  w=k1-k0;

  //check if x1 is inside Nnk(x0,X)
  loadNbh(mri,&nbh1,i0,j0,k0,inside_label);
  Nnk(&nbh1,&nbh1,connectivity);
  if(!nbh1[1+u][1+v][1+w])
    return 0;
  
  //check if Nnk(x0,Y) intersects Nnk(x1,Y)
  loadNbh(mri,&nbh1,i0,j0,k0,outside_label);
  Nnk(&nbh1,&nbh1,connectivity);
  loadNbh(mri,&nbh2,i1,j1,k1,outside_label);
  Nnk(&nbh2,&nbh2,connectivity);
  
  loop=1;
  for(a=0;loop&&a<2;a++)
    for(b=0;loop&&b<2;b++)
      for(c=0;loop&&c<2;c++)
	if((nbh1[1+a*u][1+b*v][1+c*w])&&(nbh2[1+(a-1)*u][1+(b-1)*v][1+(c-1)*w]))
	  loop=0;

  if(loop)
    return 0;
  return 1;
}

//Check if two points are weakly X-connected with respect to Y
//for the Weak Connectivity, the two points are assumed to be inside X, Y-adjacent
int checkWC(MRI* mri,int i0, int j0, int k0,int i1,int j1,int k1,int inside_label,int outside_label,int connectivity)
{
  int a,b,c,loop;
  Nbh nbh;
  
  loadNbh(mri,&nbh,i0,j0,k0,inside_label);
  Nnk(&nbh,&nbh,connectivity);
  loop=1;
  for(a=-1;loop&&a<2;a++)
    for(b=-1;loop&&b<2;b++)
      for(c=-1;loop&&c<2;c++)
	if(nbh[1+a][1+b][1+c] && checkNbh(mri,i0+a,j0+b,k0+c,outside_label,connectivity) 
	   && checkSC(mri,i0,j0,k0,i0+a,j0+b,k0+c,inside_label,outside_label,connectivity)
	   && checkSC(mri,i1,j1,k1,i0+a,j0+b,k0+c,inside_label,outside_label,connectivity))
	  loop=0;
     
  if(loop)
    return 0;

  return 1;
}

Nbh* loadNbh2(unsigned char*** im,int x,int y, int z,int label)
{
  Nbh *nbh=(Nbh*)malloc(sizeof(Nbh));
  int i,j,k;

  if(label==-1)
    for(k=-1;k<2;k++)
      for(j=-1;j<2;j++)
	for(i=-1;i<2;i++)
	  if(im[z+k][y+j][x+i])
	    (*nbh)[i+1][j+1][k+1]=1;
	  else
	    (*nbh)[i+1][j+1][k+1]=0;  
  else
    for(k=-1;k<2;k++)
      for(j=-1;j<2;j++)
	for(i=-1;i<2;i++)
	  if(im[z+k][y+j][x+i]==label)
	    (*nbh)[i+1][j+1][k+1]=1;
	  else
	    (*nbh)[i+1][j+1][k+1]=0;
  
  return nbh;
}



Nbh* reverseNbh(Nbh* nbh_src,Nbh *nbh_dst)
{
  Nbh *nbh;
  int i,j,k;
  if(nbh_dst)
    nbh=nbh_dst;
  else
    nbh=(Nbh*)malloc(sizeof(Nbh));
  
  for(k=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
	if((*nbh_src)[i][j][k])
	  (*nbh)[i][j][k]=0;
	else
	  (*nbh)[i][j][k]=1;

  return nbh;
}
