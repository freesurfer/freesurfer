#include <map>
#include <algorithm>
#include "mrisurf_vals.h"
#include <random>
#include <math.h>
#include "cma.h"
#include <vnl/vnl_matops.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

class MRIS_MultimodalRefinement {
	public:
		void refine(MRIS* surface, float lambda, int iterations);
		void  getTarget(MRIS* surface);
		void  getTargetO(MRIS* surface);
		void addImage(MRI* image)
		{
			images.push_back(image);
		}
		void SetVertexDebug(int i){ this->vertexDebug =i;}
		void SetStep(float i){ this->step =i;}
		void SetNumberOfSteps(int i){ this->numberOfSteps =i;}
		void SetGradientSigma(float i){ this->gradientSigma =i;}
		void SetSegmentation(MRI* s) { this->segmentation=s;}
		void SetWhiteMR(MRI* s) { this->whiteMR=s;}
		void SetVesselMR(MRI* s) { this->vesselMR=s;}
		void SetPosteriors(MRI* s) { this->posteriorMR= s;}
		void FindMaximumGradient(bool b){this->maximumGradient=b;}
		void SetWhite(MRIS* white){ this->white=white;}
		void SetSphere(MRIS* sphere){ this->sphere=sphere;}
		void getNormal(MRIS* surf,MRI* image, int vtxno, double* x, double* y, double* z);
		std::vector<std::pair<float, float>> GetMeanAndVariance(MRIS* surf, MRI* image);
		void  SegmentWM(MRI* t1, MRI* t2, MRI* output, int contrast_type);
		void  SegmentVessel(MRI* t1, MRI* t2, MRI* output, int contrast_type);
		void GeneratePlanes();
		float DistanceToMidline(double x, double y, double z);
		void SetMinPGrey(float p){ this->minPGrey = p;}

	private:
		MRIS* white;
		MRIS* sphere;
		int vertexDebug= -1;
		float step=0.4;
		int numberOfSteps =20;
		double T2_max=200;
		// double T2_min_gray=120;  // ATH commenting out because this is an unused field
		double T2_min=120;
		double MAX_CSF_=80;
		double AVG_CSF=80;
		float gradientSigma=0.25;
		bool maximumGradient= True;
		std::vector<MRI*> images;
		MRI* segmentation;
		MRI* whiteMR;
		MRI* vesselMR;
		MRI* posteriorMR= nullptr;
		double center[3];	
		double normal[3];
		float GetVentricleIntensity();	
		float minPGrey=20;
		
};
void MRIS_MultimodalRefinement::getTarget(MRIS* surf )
{
	std::cout << " getTarget "<< std::endl;
	std::cout << " number of steps "<< this->numberOfSteps <<std::endl;
	std::cout << " step "<< step << std::endl;
	std::cout <<" gradient sigma " << this->gradientSigma << std::endl;
	std::cout <<" Min P Grey  " << this->minPGrey << std::endl;
	std::cout << " min p " << 1.0 / pow(10,this->minPGrey) << " " << 1.0 / pow(10,this->minPGrey/3)<< std::endl;
	//std::map<std::tuple<int,int,int>, int> priors =getPs(surf);
			
  	MRIS_HASH_TABLE *mht ;
	mht = MHTcreateVertexTable_Resolution(surf, CURRENT_VERTICES, 10);

	//orientationFilter->SetBabyMode(false);
	//this->GeneratePlanes();
	std::vector<float> intensities(this->numberOfSteps*5+1,0);
	std::vector<float> magnitudes(this->numberOfSteps*5+1,0);
	std::vector<float> ps(this->numberOfSteps*5+1,1);
	std::vector<float> masks(this->numberOfSteps*5+1,0);
	std::vector<float> masksW(this->numberOfSteps*5+1,0);
	std::vector<float> segs(this->numberOfSteps*5+1,0);

	std::cout << " ps size " << ps.size() << " images size " << images.size() << std::endl;	

	std::vector<std::vector<std::pair<float, float>>> stats_vertex_images;
	
	std::cout << " holi " << std::endl;
	
	if( this->posteriorMR == nullptr)
	{
		for(int i=0;i<images.size();i++)
		{
			stats_vertex_images.push_back(  GetMeanAndVariance(surf, images[i]));
			std::cout << " ah " <<  stats_vertex_images[i][vertexDebug].first <<  " " << stats_vertex_images[i][vertexDebug].second << std::endl;
		}
	}
	for (unsigned int j=0;j<surf->nvertices;j++)
	{
		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
		surf->vertices[j].ripflag =1 ;
	
		for(int i=0;i<ps.size(); i++)
		{
			ps[i] =intensities[i] = magnitudes[i] =0;
		}
		double xv, yv, zv; //xvp, yvp, zvp, xvn, yvn, zvn;//, nxv, nyv,nzv;
		float x,y,z;
		double nx=0,ny=0,nz=0;
		getNormal(surf, images[0],j, &nx, &ny, &nz);
		//std::cout << nx << " "<<ny<< " " << nz << std::endl;
		

		if(j==vertexDebug)
			std::cout << nx << " " << ny << " " << nz << std::endl;

		float dx=pow( surf->vertices[j].x - surf->vertices[j].whitex,2);
		float dy= pow(surf->vertices[j].y - surf->vertices[j].whitey,2);
		float dz= pow(surf->vertices[j].z - surf->vertices[j].whitez,2);

		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &xv,&yv,&zv);
	 	double dmin=0;
		int numberOfHops = 2;
		SURFHOPLIST* hops = SetSurfHopList(j, surf, numberOfHops);
		for(int t=1; t<5 ; t++)
		{
			int vno1, nthface, faceno2,vno2,k; 
			float dminv;
			double d,dL=0.2;
			double pv1[3], pf1[3], pf2[3], pf3[3], pmin[3], pmin0[3];
			VERTEX v1, *v2;


			StuffVertexCoords(surf, j, pv1);
			v1 = surf->vertices[j];    
			
			double nnx = surf->vertices[j].nx;
			double nny = surf->vertices[j].ny;
			double nnz = surf->vertices[j].nz;
			double sum = sqrt(pow( nnx,2) +pow( nny,2) +  pow(nnz,2));

			v1.x = surf->vertices[j].x + (t*.5)*surf->vertices[j].nx/sum;
			v1.y = surf->vertices[j].y + (t*.5) *surf->vertices[j].ny/sum;
			v1.z = surf->vertices[j].z + (t*.5) *surf->vertices[j].nz/sum;


			// Get the closest vertex in surf2
    			vno2 = mht->findClosestVertexNoXYZ(v1.x,v1.y,v1.z, &dminv);
			//vno2 = MHTfindClosestVertexNo2(mht, surf, surf, &v1, &dminv);
	//		vno2 = MHTfindClosestVertexNo2(mht, surf, surf, &v1, &dminv);
			if( j == vertexDebug)
				std::cout << "v1  " <<  v1.x <<" " << surf->vertices[j].x <<  " "  << surf->vertices[j].nx << " " << sum<< " " << dminv <<  std::endl;
			if( vno2 > 0 && vno2 <surf->nvertices )
			{	

				int norm=0;
				bool good= true;
				for(int h=0; h<numberOfHops; h++)
				{
					for(int n=0; n< hops->nperhop[h];n++)
					{
						int vtxno = hops->vtxlist[h][n];
						if( vtxno == vno2 )
						{
							good=false;			
							break;
						}
					}
				}
				if(good ) 
					dmin+=1;
			}
		
		}
		SurfHopListFree(&hops);
		if( j == vertexDebug)
			std::cout << "distance " <<  dmin <<  " distance to white "<< dx + dy + dz <<  std::endl;
		surf->vertices[j].curv = dmin;

		double mag, val;
		bool count =false;
		int counting = 0;
		for(int t=-this->numberOfSteps; (t<=this->numberOfSteps || (ps[t+numberOfSteps-1] > (1.0/pow(10,this->minPGrey))  && t +numberOfSteps< ps.size()));t++)
		{
			//ps[t+numberOfSteps] =1;

			x=xv +nx*(t)*step;
			y=yv +ny*(t)*step;
			z=zv +nz*(t)*step;

			double whiteIntensity, vesselIntensity;
			MRIsampleVolume(this->whiteMR, x, y, z, &whiteIntensity);
			MRIsampleVolume(this->vesselMR, x, y, z, &vesselIntensity);
			masks[t+numberOfSteps]= vesselIntensity;
			masksW[t+numberOfSteps]= whiteIntensity;
			double label=0;	
			MRIsampleVolumeFrameType(this->segmentation, x, y, z, 0, SAMPLE_NEAREST, &label);
			segs[t+numberOfSteps]= label;
			for(int k=0;k<images.size();k++)
			{
				MRIsampleVolume(images[k], x, y, z, &val);

				MRIsampleVolumeDerivativeScale(images[k], x, y, z, nx, ny, nz, &mag, this->gradientSigma); //, step, vertexDebug==j);				
				//probability of grey
				//float p = (1.0/sqrt(2*varGrey)) *exp( - pow(val-meanGrey,2) /(2*varGrey));
				float p = exp( - pow(val-stats_vertex_images[k][j].first,2) /(2*stats_vertex_images[k][j].second));

				magnitudes[t+numberOfSteps]+= fabs(mag) /images.size();
				intensities[t+numberOfSteps]+=val/images.size();
				ps[t+numberOfSteps]+=(1-whiteIntensity) * p; ///images.size();
			}
			if(vertexDebug==j)
			{
				std::cout << "distance from ori " << t*step <<  " " << magnitudes[t+numberOfSteps] <<   " " << val << " " << ps[t+numberOfSteps] << "  white intensity "<< whiteIntensity << " vesselIntensity " << vesselIntensity  << std::endl;
			}
			if(ps[t+numberOfSteps]>.5)
				count = true;
			if(count && ps[t+numberOfSteps] ==0)
				counting++;
			else
				counting =0;

			if(counting ==2)	
				break;
		}
		if( vertexDebug == j)
			std::cout << "dist to white" <<  (dx + dy + dz   ) << " dist to next sulcus " <<   dmin  <<  std::endl;
		if(((dx + dy + dz >1  ) &&  (dmin <=3))) 
		{
			

			surf->vertices[j].ripflag =false;
			float opt_mag=0,opt_val;
			int opt_t;
			/*float leftW= this->DistanceToMidline(surf->vertices[j].whitex, surf->vertices[j].whitey, surf->vertices[j].whitez);
			leftW /= fabs(leftW);*/
			int  label=0, prevLabel=0, lastCortexLabel=0;
			float touchedStructure =0;
			float changeHemis=0;
			float changeHemisAseg=0;
			bool good=true; 
			int zeroLabel=0;
			for(int i=1;(i< numberOfSteps+1 || ( ps[i-1] > (1.0/pow(10,this->minPGrey)) && i < ps.size() -1));i++)
			{
				if (vertexDebug ==j)
					std::cout <<" opt mag " << opt_mag << " mag " << magnitudes[i] << " " << ps[i] << " " << ps[i-1] << " " << ps[i-2] << std::endl;
				x=xv  +nx*(i -numberOfSteps)*step;
				y= yv +ny*(i-numberOfSteps)*step;
				z=zv +nz*(i-numberOfSteps)*step;

				double point[3] = {x,y,z};

				//avoid white matter, brainstem, cerebellum and hyppocamps
				double xt,yt,zt;
				MRIvoxelToSurfaceRAS(images[0], x,y,z, &xt,&yt,&zt);

				/*float leftP = this->DistanceToMidline(xt,yt,zt);				
				leftP /= fabs(leftP);*/
				label=segs[i];
				if(  label ==53 || label == 17 || label == 16 || label == 47 || label == 8 || ((label != 2 && label !=  41 )&& masksW[i] <.1 && masks[i+1] > .1 && masksW[i+1]<.01) ) 
				{
					touchedStructure ++;
				}
				else if (( segs[i+1]>1000 && lastCortexLabel > 1000  && fabs(segs[i+1] -lastCortexLabel) > 900 ) || ( segs[i+1] ==42  && lastCortexLabel ==3) || (segs[i+1] ==3 && lastCortexLabel ==42)  )
				{
					changeHemisAseg ++;
					changeHemis ++;
				}
				if (label== 0)
				{ 	/*if(fabs(leftW + leftP)  ==0)  
					{
						changeHemis ++;
					}*/
					zeroLabel++;		
				}
				
				if (vertexDebug ==j)
					std::cout << lastCortexLabel << " " << segs[i+1] << " " << touchedStructure <<   std::endl;
			//label = (label>.5)?1:0;
				if((fabs(magnitudes[i])>opt_mag && ps[i] > (1.0/pow(10,this->minPGrey)))  || touchedStructure>0 || (changeHemis > 0  && fabs(magnitudes[i]) *1.5 > opt_mag) ||( changeHemisAseg >0&& zeroLabel <3 &&opt_mag < 5))
				{

					int w=i;
					good = true;
					for( ; w>0 && good && (ps[w-1] < (1.0/pow(10, this->minPGrey /3.0)) || fabs(i-w) <1);w--)
					{
						good = ps[w-1] >= ps [w] && masks[w] <.1 ;
						if(vertexDebug==j)
							std::cout << ps[w-1 ] << " " << ps[w] <<  " " << w << " " << masks[w]<< std::endl;
					}
					if((((ps[w-1] > (1.0/pow(10, this->minPGrey/3.0)) )|| good ) || (touchedStructure >0  && ps[i-1] > 0.0001 )) && label != 41 && label !=2 )  
					{
						opt_mag = fabs( magnitudes[i]);
						opt_val= intensities[i];
						opt_t =  i;
						surf->vertices[j].targx = xt ;
						surf->vertices[j].targy = yt;
						surf->vertices[j].targz = zt;

						if(vertexDebug==j)
						{
							std::cout << "TAKING: distance from ori " << (i-numberOfSteps)*step <<  " " << opt_mag <<   " " << opt_val << std::endl;

						}
					}
				}		
				if (  touchedStructure>0  ||changeHemis >1  || changeHemisAseg >0) //|| (!IS_BRAIN(label) && masks[i] >.5 ))
					break;
	
				if (label!= 0)
				{
					lastCortexLabel=label;
					zeroLabel=0;
				}
				prevLabel= label;
			}
		}
	}
	MHTfree(&mht);
}
void MRIS_MultimodalRefinement::getTargetO(MRIS* surf )
{
	std::cout << " getTarget "<< std::endl;
	std::cout << " number of steps "<< this->numberOfSteps <<std::endl;
	std::cout << " step "<< step << std::endl;
	std::cout <<" gradient sigma " << this->gradientSigma << std::endl;
	//std::map<std::tuple<int,int,int>, int> priors =getPs(surf);
			
  	MRIS_HASH_TABLE *mht ;
	mht = MHTcreateVertexTable_Resolution(surf, CURRENT_VERTICES, 10);

	//orientationFilter->SetBabyMode(false);
	this->GeneratePlanes();
	std::vector<float> intensities(this->numberOfSteps*5+1,0);
	std::vector<float> magnitudes(this->numberOfSteps*5+1,0);
	std::vector<float> ps(this->numberOfSteps*5+1,0);
	std::vector<float> masks(this->numberOfSteps*5+1,0);
	std::vector<float> masksW(this->numberOfSteps*5+1,0);
	std::vector<float> segs(this->numberOfSteps*5+1,0);

	std::cout << " ps size " << ps.size() << " images size " << images.size() << std::endl;	

	std::vector<std::vector<std::pair<float, float>>> stats_vertex_images;
	
	std::cout << " holi " << std::endl;

	for(int i=0;i<images.size();i++)
	{
		stats_vertex_images.push_back(  GetMeanAndVariance(surf, images[i]));
		std::cout << " ah " <<  stats_vertex_images[i][vertexDebug].first <<  " " << stats_vertex_images[i][vertexDebug].second << std::endl;
	}
	for (unsigned int j=0;j<surf->nvertices;j++)
	{
		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
		surf->vertices[j].ripflag =1 ;
	
		for(int i=0;i<ps.size(); i++)
		{
			ps[i] =intensities[i] = magnitudes[i] =0;
		}
		double xv, yv, zv; //xvp, yvp, zvp, xvn, yvn, zvn;//, nxv, nyv,nzv;
		float x,y,z;
		double nx=0,ny=0,nz=0;
		getNormal(surf, images[0],j, &nx, &ny, &nz);
		//std::cout << nx << " "<<ny<< " " << nz << std::endl;
		

		if(j==vertexDebug)
			std::cout << nx << " " << ny << " " << nz << std::endl;

		float dx=pow( surf->vertices[j].x - surf->vertices[j].whitex,2);
		float dy= pow(surf->vertices[j].y - surf->vertices[j].whitey,2);
		float dz= pow(surf->vertices[j].z - surf->vertices[j].whitez,2);

		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &xv,&yv,&zv);
	 	double dmin=0;
		int numberOfHops = 2;
		SURFHOPLIST* hops = SetSurfHopList(j, surf, numberOfHops);
		for(int t=1; t<5 ; t++)
		{
			int vno1, nthface, faceno2,vno2,k; 
			float dminv;
			double d,dL=0.2;
			double pv1[3], pf1[3], pf2[3], pf3[3], pmin[3], pmin0[3];
			VERTEX v1, *v2;


			StuffVertexCoords(surf, j, pv1);
			v1 = surf->vertices[j];    
			
			double nnx = surf->vertices[j].nx;
			double nny = surf->vertices[j].ny;
			double nnz = surf->vertices[j].nz;
			double sum = sqrt(pow( nnx,2) +pow( nny,2) +  pow(nnz,2));

			v1.x = surf->vertices[j].x + (t*.5)*surf->vertices[j].nx/sum;
			v1.y = surf->vertices[j].y + (t*.5) *surf->vertices[j].ny/sum;
			v1.z = surf->vertices[j].z + (t*.5) *surf->vertices[j].nz/sum;


			// Get the closest vertex in surf2
    			vno2 = mht->findClosestVertexNoXYZ(v1.x,v1.y,v1.z, &dminv);
			//vno2 = MHTfindClosestVertexNo2(mht, surf, surf, &v1, &dminv);
	//		vno2 = MHTfindClosestVertexNo2(mht, surf, surf, &v1, &dminv);
			if( j == vertexDebug)
				std::cout << "v1  " <<  v1.x <<" " << surf->vertices[j].x <<  " "  << surf->vertices[j].nx << " " << sum<< " " << dminv <<  std::endl;
			if( vno2 > 0 && vno2 <surf->nvertices )
			{	

				int norm=0;
				bool good= true;
				for(int h=0; h<numberOfHops; h++)
				{
					for(int n=0; n< hops->nperhop[h];n++)
					{
						int vtxno = hops->vtxlist[h][n];
						if( vtxno == vno2 )
						{
							good=false;			
							break;
						}
					}
				}
				if(good ) 
					dmin+=1;
			}
		
		}
		SurfHopListFree(&hops);
		if( j == vertexDebug)
			std::cout << "distance " <<  dmin <<  " distance to white "<< dx + dy + dz <<  std::endl;
		surf->vertices[j].curv = dmin;

		double mag, val;
		bool count =false;
		int counting = 0;
		double vertexMask =0;
		for(int t=-this->numberOfSteps; (t<=this->numberOfSteps || (ps[t+numberOfSteps-1] >1.0e-20  && t +numberOfSteps< ps.size()));t++)
		{
			//ps[t+numberOfSteps] =1;

			x=xv +nx*(t)*step;
			y=yv +ny*(t)*step;
			z=zv +nz*(t)*step;

			double whiteIntensity, vesselIntensity;
			MRIsampleVolume(this->whiteMR, x, y, z, &whiteIntensity);
			MRIsampleVolume(this->vesselMR, x, y, z, &vesselIntensity);
			masks[t+numberOfSteps]= vesselIntensity;
			masksW[t+numberOfSteps]= whiteIntensity;
			if(t <this->numberOfSteps && vesselIntensity > .1)
				vertexMask += vesselIntensity;
			double label=0;	
			MRIsampleVolumeFrameType(this->segmentation, x, y, z, 0, SAMPLE_NEAREST, &label);
			segs[t+numberOfSteps]= label;
			for(int k=0;k<images.size();k++)
			{
				MRIsampleVolume(images[k], x, y, z, &val);

				MRIsampleVolumeDerivativeScale(images[k], x, y, z, nx, ny, nz, &mag, this->gradientSigma); //, step, vertexDebug==j);				
				//probability of grey
				//float p = (1.0/sqrt(2*varGrey)) *exp( - pow(val-meanGrey,2) /(2*varGrey));
				float p = exp( - pow(val-stats_vertex_images[k][j].first,2) /(2*stats_vertex_images[k][j].second));

				magnitudes[t+numberOfSteps]+= fabs(mag) /images.size();
				intensities[t+numberOfSteps]+=val/images.size();
				ps[t+numberOfSteps]+=(1-whiteIntensity) * p/images.size();
			}
			if(vertexDebug==j)
			{
				std::cout << "distance from ori " << t*step <<  " " << magnitudes[t+numberOfSteps] <<   " " << val << " " << ps[t+numberOfSteps] << "  white intensity "<< whiteIntensity << " vesselIntensity " << vesselIntensity  << std::endl;
			}
			if(ps[t+numberOfSteps]>.5)
				count = true;
			if(count && ps[t+numberOfSteps] ==0)
				counting++;
			else
				counting =0;

			if(counting ==2)	
				break;
		}
		if( vertexDebug == j)
			std::cout << "dist to white" <<  (dx + dy + dz   ) << " dist to next sulcus " <<   dmin  << " vessel ?  " << vertexMask  << std::endl;
		if(((dx + dy + dz >1  ) &&  (dmin <=3))) // || vertexMask >.5)
		{
			

			surf->vertices[j].ripflag =false;
			float opt_mag=0,opt_val;
			int opt_t;
			float leftW= this->DistanceToMidline(surf->vertices[j].whitex, surf->vertices[j].whitey, surf->vertices[j].whitez);
			leftW /= fabs(leftW);
			int  label=0, prevLabel=0, lastCortexLabel=0;
			float touchedStructure =0;
			float changeHemis=0;
			float changeHemisAseg=0;
			bool good=true; 
			int zeroLabel=0;
			for(int i=1;(i< numberOfSteps+1 || ( ps[i-1] > 1e-15 && i < ps.size() -1));i++)
			{
				if (vertexDebug ==j)
					std::cout <<" opt mag " << opt_mag << " mag " << magnitudes[i] << " " << ps[i] << " " << ps[i-1] << " " << ps[i-2] << std::endl;
				x=xv  +nx*(i -numberOfSteps)*step;
				y= yv +ny*(i-numberOfSteps)*step;
				z=zv +nz*(i-numberOfSteps)*step;

				double point[3] = {x,y,z};

				//avoid white matter, brainstem, cerebellum and hyppocamps
				double xt,yt,zt;
				MRIvoxelToSurfaceRAS(images[0], x,y,z, &xt,&yt,&zt);

				float leftP = this->DistanceToMidline(xt,yt,zt);				
				leftP /= fabs(leftP);
				label=segs[i];
				if(  label ==53 || label == 17 || label == 16 || label == 47 || label == 8 || ((label != 2 && label !=  41 )&& masksW[i] <.1 && masks[i+1] > .1 && masksW[i+1]<.01) ) 
				{
					touchedStructure ++;
				}
				else if (( segs[i+1]>1000 && lastCortexLabel > 1000  && fabs(segs[i+1] -lastCortexLabel) > 900 ) || ( segs[i+1] ==42  && lastCortexLabel ==3) || (segs[i+1] ==3 && lastCortexLabel ==42)  )
				{
					changeHemisAseg ++;
					changeHemis ++;
				}
				if (label== 0)
				{ 	if(fabs(leftW + leftP)  ==0)  
					{
						changeHemis ++;
					}
					zeroLabel++;		
				}
				
				if (vertexDebug ==j)
					std::cout << lastCortexLabel << " " << segs[i+1] << " " << touchedStructure << " distance midplane " << leftP << " " <<leftW << std::endl;
			//label = (label>.5)?1:0;
				if((fabs(magnitudes[i])>opt_mag && ps[i] > 1e-15)  || touchedStructure>0 || (changeHemis > 0  && fabs(magnitudes[i]) *1.5 > opt_mag) ||( changeHemisAseg >0&& zeroLabel <3 &&opt_mag < 5))
				{

					int w=i;
					good = true;
					for( ; w>0 && good && (ps[w-1] <1e-05  || fabs(i-w) <1);w--)
					{
						good = ps[w-1] >= ps [w] && masks[w] <.1 ;
						if(vertexDebug==j)
							std::cout << ps[w-1 ] << " " << ps[w] <<  " " << w << " " << masks[w]<< std::endl;
					}
					if((((ps[w-1] > 1e-05 )|| good ) || (touchedStructure >0  && ps[i-1] > 1e-1 )) && label != 41 && label !=2 )  
					{
						opt_mag = fabs( magnitudes[i]);
						opt_val= intensities[i];
						opt_t =  i;
						surf->vertices[j].targx = xt ;
						surf->vertices[j].targy = yt;
						surf->vertices[j].targz = zt;

						if(vertexDebug==j)
						{
							std::cout << "TAKING: distance from ori " << (i-numberOfSteps)*step <<  " " << opt_mag <<   " " << opt_val << std::endl;

						}
					}
				}		
				if (  touchedStructure>0  ||changeHemis >1  || changeHemisAseg >0) //|| (!IS_BRAIN(label) && masks[i] >.5 ))
					break;
	
				if (label!= 0)
				{
					lastCortexLabel=label;
					zeroLabel=0;
				}
				prevLabel= label;
			}
		}
	}
	MHTfree(&mht);
}
std::vector<std::pair<float, float>> MRIS_MultimodalRefinement::GetMeanAndVariance(MRIS* surf, MRI* image)
{

	int numberOfHops =10;
	std::vector<std::pair<float,float>> mean_variance;
	std::map<int,std::vector<double>> values;
	for(int j = 0;j<surf->nvertices;j++)
	{
		SURFHOPLIST* hops = SetSurfHopList(j, surf, numberOfHops);
		double meanGrey = 0;
		double varGrey = 0;
		
		std::vector<float> greyValues ;
		for(int h=0; h<numberOfHops; h++)
		{
			for(int n=0; n< hops->nperhop[h];n++)
			{
				int vtxno = hops->vtxlist[h][n];
				if( values.find(vtxno) == values.end() )
				{
					std::vector<double> hola;

					float dx= surf->vertices[vtxno].x - surf->vertices[vtxno].whitex;
					float dy= surf->vertices[vtxno].y - surf->vertices[vtxno].whitey;
					float dz= surf->vertices[vtxno].z - surf->vertices[vtxno].whitez;
					for(int i =2;i<7; i+=2)
					{
						float x  = surf->vertices[vtxno].x -dx/i;
						float y  = surf->vertices[vtxno].y -dy/i;
						float z  = surf->vertices[vtxno].z -dz/i;
						double xv, yv,zv;	
						double val, valmask;
						MRISsurfaceRASToVoxel(surf, image, x,y,z, &xv,&yv,&zv);
						MRIsampleVolume(image, xv, yv, zv, &val);
						MRIsampleVolume(whiteMR, xv, yv, zv, &valmask);
						if( val >30 && val <300 && valmask <.3)
						{
							hola.push_back(val);
							greyValues.push_back(val);
							meanGrey+=val;
						}	
					}
					values[vtxno]=hola;
				}
				else	
				{
					for(int i=0;i<values[vtxno].size();i++)
					{
						//std::cout << (values[vtxno])[i] << std::endl;
						greyValues.push_back((values[vtxno])[i]);
						meanGrey += (values[vtxno])[i];
						
					}
				}
	
			} 
		}
		meanGrey/= greyValues.size();
		for(int n=0; n<greyValues.size();n++)
		{
			varGrey  += pow(greyValues[n] - meanGrey,2);
		}	

		varGrey/= (greyValues.size());
		varGrey =  sqrt(varGrey);
		mean_variance.push_back(std::pair<float,float> (meanGrey,varGrey));

		SurfHopListFree(&hops);
	}
	std::cout << " hi " << std::endl;
	return mean_variance;
}
void  MRIS_MultimodalRefinement::SegmentWM(MRI* t1, MRI* t2, MRI* output, int contrast_type)
{
	for (int x = 0 ; x < t1->width ; x++)
	{
		for (int y = 0 ; y < t1->height ; y++)
		{
			for (int z = 0 ; z < t1->depth ; z++)
			{
				//int label =(imageAseg)? MRIgetVoxVal(imageAseg, x, y, z, 0) :0;
				float T1val = MRIgetVoxVal(t1, x, y, z, 0) ;
				float T2val = MRIgetVoxVal(t2, x, y, z, 0) ;

				int val = 0;
				if(contrast_type <0)
				{
					double label;
					MRIsampleVolumeFrameType(t2, x, y, z, 0, SAMPLE_NEAREST, &label);
					if(label==2 || label==41)
					{
						val=1;
					}
				} else 	if(1.3*T1val > T2val &&  (T1val +  T2val) >= 180 && T1val > 80 && T2val > 80)
				{
					val =1;
				}

				MRIsetVoxVal(output, x, y, z, 0, val) ;
			}
		}
	}


}
void  MRIS_MultimodalRefinement::SegmentVessel(MRI* t1, MRI* t2, MRI* output, int contrast_type)
{

	int sum=180;
	int max=80;
	int min=110;
	float mult=1.4;
	std::cout << contrast_type << " " << CONTRAST_FLAIR<<" flair" <<std::endl;
	if (contrast_type == CONTRAST_FLAIR)
	{
		mult=1;
		max=50;
	}
	for (int x = 0 ; x < t1->width ; x++)
	{
		for (int y = 0 ; y < t1->height ; y++)
		{
			for (int z = 0 ; z < t1->depth ; z++)
			{
				//int label =(imageAseg)? MRIgetVoxVal(imageAseg, x, y, z, 0) :0;
				float T1val = MRIgetVoxVal(t1, x, y, z, 0) ;
				float T2val = MRIgetVoxVal(t2, x, y, z, 0) ;
				int val = 0;
				if( contrast_type <0)
				{  
					double label;
					MRIsampleVolumeFrameType(t2, x, y, z, 0, SAMPLE_NEAREST, &label);
					/*if(label==0)
					{
						val=1;
					}*/
				
				} else if (  (mult*T1val > T2val && (T1val +T2val <sum)) || (T1val <max && T2val <max) || (T1val > min && T2val>min))
				{
					val=1;
				}

				MRIsetVoxVal(output, x, y, z, 0, val) ;
			}
		}
	}


}
// Calculate the unit-length normal to the vertex in VOXEL space
void MRIS_MultimodalRefinement::getNormal(MRIS* surf, MRI* image, int vtxno,  double* nx, double* ny, double* nz)
{
	double x,y,z;
	double xw, yw, zw;
	double xw1, yw1, zw1;
	x = surf->vertices[vtxno].x;
	y = surf->vertices[vtxno].y;
	z = surf->vertices[vtxno].z;
	MRISsurfaceRASToVoxel(surf, image, x , y,z, &xw,&yw,&zw);

	x = surf->vertices[vtxno].x + surf->vertices[vtxno].nx;
	y = surf->vertices[vtxno].y + surf->vertices[vtxno].ny;
	z = surf->vertices[vtxno].z + surf->vertices[vtxno].nz;

	MRISsurfaceRASToVoxel(surf, image, x , y,z, &xw1,&yw1,&zw1);
	*nx = xw1 - xw;
	*ny = yw1 - yw;
	*nz = zw1 - zw;
	float dist = sqrt(SQR(*nx) + SQR(*ny) + SQR(*nz));
	*nx /= dist;
	*ny /= dist;
	*nz /= dist;
}
void MRIS_MultimodalRefinement::GeneratePlanes()
{
	std::vector<double*> values;
	this->center[0] =0;
	this->center[1] =0;
	this->center[2]=0;
	bool m_baby =false;
 	for(int i=0;i<segmentation->width;i++)
	{
		for(int j=0;j<segmentation->height;j++)
		{
			for(int k=0;k<segmentation->depth;k++)
			{	
				int pixel = MRIgetVoxVal(segmentation,i,j,k,0);
				double v[3];
				MRIvoxelToSurfaceRAS(segmentation, i,j,k,&v[0] , &v[1], &v[2]);
				if(pixel == 14 || pixel == 254 || pixel ==253 || pixel ==252 ||   (pixel ==3026 && m_baby)|| (pixel ==4026 && m_baby) || (pixel== 3010 && m_baby)|| (pixel== 4010 && m_baby) || (pixel ==175 && m_baby)) // || pixel == 255 || pixel == 251)
				{
					values.push_back(v);
					for(int w =0;w<3;w++)
						this->center[w] += v[w];
				}
			}
		}
	}
	//This is my center of mass
	for(int w =0;w<3;w++)
		this->center[w] = this->center[w]/values.size();
	

	vnl_matrix<double> points(values.size(), 4);     

	double weightSum = 0;
	for(int i = 0; i<values.size();i++ )
	{
		for(int w =0;w<3;w++)
			points(i,w) = values[i][w];
		points(i,3)=1;
	}

	vnl_symmetric_eigensystem<double> eig(points.transpose()*points);

	this->normal[0] = eig.get_eigenvector(0)[2];
	this->normal[1] = eig.get_eigenvector(1)[2];
	this->normal[2] = eig.get_eigenvector(2)[2];
	std::cout << normal[0] << normal[1] << normal[2] <<std::endl;
}

float MRIS_MultimodalRefinement::DistanceToMidline(double x, double y, double z)
{
	float d = 0;
	float norm =0;
	for (int w=0; w<3;w++)
	{	
		d  -= normal[w]*center[w];
		norm += normal[w]* normal[w];
	}
	norm = std::sqrt(norm);
	float dist =  (normal[0]*x + normal[1]*y + normal[2] *z + d ) /norm;
	
	return dist;
}


