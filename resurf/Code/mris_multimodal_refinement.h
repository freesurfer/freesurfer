#include <map>

class MRIS_MultimodalRefinement {
	public:
		void refine(MRIS* surface, float lambda, int iterations);
		void  getTarget(MRIS* surface);
		std::map<std::tuple<int, int, int>, int>  getPs(MRIS* surface);
		void  getTarget2(MRIS* surface);
		void  getTargetAdHoc(MRIS* surface);
		void addImage(MRI* image)
		{
			images.push_back(image);
		}
		void SetVertexDebug(int i){ this->vertexDebug =i;}
		void SetStep(float i){ this->step =i;}
		void SetNumberOfSteps(int i){ this->numberOfSteps =i;}
		void SetGradientSigma(float i){ this->gradientSigma =i;}
		void SetSegmentation(MRI* s) { this->segmentation=s;}
		void FindMaximumGradient(bool b){this->maximumGradient=b;}
		void SetWhite(MRIS* white){ this->white=white;}
	private:
		MRIS* white;
		int vertexDebug= -1;
		float step=0.4;
		int numberOfSteps =20;
		double T2_max=200;
		// double T2_min_gray=120;  unused
		double T2_min=120;
		double MAX_CSF=80;
		double AVG_CSF=80;
		float gradientSigma=0.25;
		bool maximumGradient= True;
		std::vector<MRI*> images;
		MRI* segmentation;

		float GetVentricleIntensity();	
		bool isCSF(double val)
		{
			return (val<MAX_CSF);
		}
		bool continuousDecrease(double a, double b, double c)
		{
			//return (a<b && b<c);
			return (a - b + b -c )>0;
		}
		bool continuousIncrease(double a, double b, double c)
		{
			//return (a>b && b>c);
			return (a-b+ b -c)<0;
		}
		bool isGrey(double val)
		{
			return val > T2_min && val < T2_max;
		}
		bool isPial(double val)
		{
			return val >T2_max;
		}
		float score(double val, double max_val, double mag, double max_mag, int t, int max_t)
		{

			//		std::cout << val/(max_val+1.1)+ fabs(mag)/(max_mag+1.1)  << " " ; //+ t/float(max_t);
			return val/(max_val+101)+ fabs(mag)/(max_mag+1.1) ; //+ t/float(max_t);

			//	((val >= .9*max_val && mag > max_mag ) || ( val > max_val && mag >=.9 *max_mag)|| t<max_t) ;

		}
		
};
void MRIS_MultimodalRefinement::getTarget(MRIS* surf )
{
	std::cout << " getTarget "<< std::endl;
	std::cout << " number of steps "<< this->numberOfSteps <<std::endl;
	std::cout << " step "<< step << std::endl;
	std::cout <<" gradient sigma " << this->gradientSigma << std::endl;
	AVG_CSF = GetVentricleIntensity();
	std::cout << "Average intensity " << AVG_CSF << std::endl;
	std::map<std::tuple<int,int,int>, int> priors =getPs(surf);
	//P(CSF)
	float pCSF=0, countCSF=0, pCount = 0;
	for ( std::map<std::tuple<int,int,int>, int>::iterator it= priors.begin(); it!= priors.end(); ++it)
	{	
		if (std::get<2>(it->first) ==1)
			countCSF += it->second;
		pCount += it->second;
	}	
	pCSF= (float)countCSF / (float)pCount;
	std::cout << pCSF << " " << countCSF << " " << pCount  << std::endl;
	for (unsigned j=0;j<surf->nvertices;j++)
	{
		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
	

		double xv, yv, zv; //xvp, yvp, zvp, xvn, yvn, zvn;//, nxv, nyv,nzv;
		float x,y,z;
		double nx,ny,nz;
		std::vector<float> intensities(numberOfSteps*2+1,0);
		std::vector<float> magnitudes(numberOfSteps*2+1,0);
		// Calculate the unit-length normal to the vertex in VOXEL space
		{
			double x,y,z;
			double xw, yw, zw;
			double xw1, yw1, zw1;
			x = surf->vertices[j].x;
			y = surf->vertices[j].y;
			z = surf->vertices[j].z;
			MRISsurfaceRASToVoxel(surf, images[0], x , y,z, &xw,&yw,&zw);

			x = surf->vertices[j].x + surf->vertices[j].nx;
			y = surf->vertices[j].y + surf->vertices[j].ny;
			z = surf->vertices[j].z + surf->vertices[j].nz;

			MRISsurfaceRASToVoxel(surf, images[0], x , y,z, &xw1,&yw1,&zw1);
			nx = xw1 - xw;
			ny = yw1 - yw;
			nz = zw1 - zw;
			float dist = sqrt(SQR(nx) + SQR(ny) + SQR(nz));
			if (FZERO(dist)) break;  
			nx /= dist;
			ny /= dist;
			nz /= dist;
		}

		if(j==vertexDebug)
			std::cout << nx << " " << ny << " " << nz << std::endl;

		double mag, val;

		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &xv,&yv,&zv);
		for(int t=-this->numberOfSteps; t<=this->numberOfSteps;t++)
		{
			//for(int k=0;k<images.size();k++)
			//int k=0;
			{
				x=xv +nx*(t)*step;
				y=yv +ny*(t)*step;
				z=zv +nz*(t)*step;

				MRIsampleVolume(images[0], x, y, z, &val);

				//MRIsampleVolume(images[k], xv, yv, zv, &val);
//				MRIsampleVolumeDerivativeScaleV(images[0], x, y, z, nx, ny, nz, &mag, this->gradientSigma, step, vertexDebug==j);				
				MRIsampleVolumeDerivativeScale(images[0], x, y, z, nx, ny, nz, &mag, this->gradientSigma); //, step, vertexDebug==j);				
				//MRIsampleVolumeDerivative(images[k], x, y, z, nx, ny, nz, &mag);  

				
				/*double next_val, prev_val;
				float norm=0, valn=0, valp=0;
				for (int i = 1; i<2 ;i++) 
				{
					float w = exp(-pow((i-1)*step,2) / (2 * pow(this->gradientSigma ,2)));
					norm += w;
					x= xv +nx*(i+t)*step;
					y=yv +ny*(i+t)*step;
					z=zv +nz*(i+t)*step;

					MRIsampleVolume(images[k], x, y, z, &next_val);
					valn += w * next_val;

					x=xv +nx*(t-i)*step;
					y=yv +ny*(t-i)*step;
					z=zv +nz*(t-i)*step;

					MRIsampleVolume(images[k], x, y, z, &prev_val);
					valp += w * prev_val;
				}
				if(vertexDebug==j)
					std::cout << valn << " "<<val <<" "<< norm <<std::endl;
				valn /= norm;
				valp /= norm;
				mag = (valn - valp) ; */
				magnitudes[t+numberOfSteps]+= mag;
				intensities[t+numberOfSteps]+=val;
			}
					if(vertexDebug==j)
					{
						std::cout << "distance from ori " << t*step <<  " " << mag <<   " " << val << std::endl;

					}
							
		}
		float in=0, out=0; 
		for(int i=0;i<intensities.size()/2;i++)
		{
			in+= intensities[i];		
			out+= intensities[i+this->numberOfSteps];		
		}
		in= (int) (in)/ (this->numberOfSteps*100); //put it between 0 and 3
		out = (int) (out)/ (this->numberOfSteps*100);
		//std::cout << in << " " <<out << std::endl;
		if(in >3)
			in=3;
		if(out >3)
			out=3;
	
		
		std::tuple<int,int,int> tuple = std::make_tuple((int)in, (int)out, 1);
		std::tuple<int,int,int> tuple2 = std::make_tuple((int)in, (int)out, 0);
		float p=0;
		int forCSF = (priors.count(tuple) >0) ? priors[tuple]:0;
		int forSkull = (priors.count(tuple2) >0) ? priors[tuple2]:0;
		if(forCSF +forSkull >0)
			p = forCSF/ (float)(forCSF +forSkull);
		
		surf->vertices[j].curv = p;
		float opt_mag=0,opt_val;
		int opt_t;
		for(int i=0;i<intensities.size();i++)
		{
			if((magnitudes[i]>opt_mag && this->maximumGradient && p>=.6) ||
				(magnitudes[i]<opt_mag && (!this->maximumGradient || p<.4)))
			{
				x=xv  +nx*(i -numberOfSteps)*step;
				y= yv +ny*(i-numberOfSteps)*step;
				z=zv +nz*(i-numberOfSteps)*step;

				//MRISsurfaceRASToVoxel(surf, images[0],x,y,z, &xv,&yv,&zv);
				//MRIsampleVolumeType(segmentation, xv, yv, zv, &val, SAMPLE_NEAREST);
		//		MRIsampleVolume(segmentation, x, y, z, &val);
	
		//		if(val <.5)
				{
					opt_mag = magnitudes[i];
					opt_val= intensities[i];
					opt_t =  i;
					double xt,yt,zt;
					MRIvoxelToSurfaceRAS(images[0], x,y,z, &xt,&yt,&zt);
					surf->vertices[j].targx = xt ;
					surf->vertices[j].targy = yt;
					surf->vertices[j].targz = zt;

					if(vertexDebug==j)
					{
						std::cout << "TAKING: distance from ori " << (i-numberOfSteps)*step <<  " " << opt_mag <<   " " << opt_val << std::endl;

					}
				}
			}		
		}

		
		/*
		float opt_mag=0,opt_val;
		int opt_t;
		bool foundCSF=false;
		for(int i=0;i<intensities.size();i++)
		{
			//if(fabs(magnitudes[i])>=fabs(max_mag))
			if(intensities[i] >250)
			{
				foundCSF=true;
			}	
			if((magnitudes[i]>opt_mag && this->maximumGradient) ||
				(magnitudes[i]<opt_mag && !this->maximumGradient ))
			{
				x=xv  +nx*(i -numberOfSteps)*step;
				y= yv +ny*(i-numberOfSteps)*step;
				z=zv +nz*(i-numberOfSteps)*step;

				//MRISsurfaceRASToVoxel(surf, images[0],x,y,z, &xv,&yv,&zv);
				//MRIsampleVolumeType(segmentation, xv, yv, zv, &val, SAMPLE_NEAREST);
		//		MRIsampleVolume(segmentation, x, y, z, &val);
	
		//		if(val <.5)
				{
					opt_mag = magnitudes[i];
					opt_val= intensities[i];
					opt_t =  i;
					double xt,yt,zt;
					MRIvoxelToSurfaceRAS(images[0], x,y,z, &xt,&yt,&zt);
					surf->vertices[j].targx = xt ;
					surf->vertices[j].targy = yt;
					surf->vertices[j].targz = zt;

					if(vertexDebug==j)
					{
						std::cout << "TAKING: distance from ori " << (i-numberOfSteps)*step <<  " " << opt_mag <<   " " << opt_val << std::endl;

					}
				}
			}		
		}
 		opt_mag=0;
		if(!foundCSF && this->maximumGradient)
		{
			for(int i=0;i<intensities.size();i++)
			{
				if((magnitudes[i])<opt_mag )
				{
					x=xv  +nx*(i -numberOfSteps)*step;
					y= yv +ny*(i-numberOfSteps)*step;
					z=zv +nz*(i-numberOfSteps)*step;

					opt_mag = magnitudes[i];
					opt_val= intensities[i];
					opt_t =  i;
					double xt,yt,zt;
					MRIvoxelToSurfaceRAS(images[0], x,y,z, &xt,&yt,&zt);
					surf->vertices[j].targx = xt ;
					surf->vertices[j].targy = yt;
					surf->vertices[j].targz = zt;

					if(vertexDebug==j)
					{
						std::cout << "TAKING: distance from ori " << (i-numberOfSteps)*step <<  " " << opt_mag <<   " " << opt_val << std::endl;

					}
				}		

			}
		}*/

	}
}
std::map<std::tuple<int,int,int>, int> MRIS_MultimodalRefinement::getPs(MRIS* surf )
{
	std::map<std::tuple<int, int, int>, int>  priors ; 
	for (unsigned j=0;j<surf->nvertices;j++)
	{
		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
	

		double xv, yv, zv; //xvp, yvp, zvp, xvn, yvn, zvn;//, nxv, nyv,nzv;
		float x,y,z;
		double nx,ny,nz;
		std::vector<float> intensities(numberOfSteps*2+1,0);
		std::vector<float> magnitudes(numberOfSteps*2+1,0);
		{
			double x,y,z;
			double xw, yw, zw;
			double xw1, yw1, zw1;
			x = surf->vertices[j].x;
			y = surf->vertices[j].y;
			z = surf->vertices[j].z;
			MRISsurfaceRASToVoxel(surf, images[0], x , y,z, &xw,&yw,&zw);

			x = surf->vertices[j].x + surf->vertices[j].nx;
			y = surf->vertices[j].y + surf->vertices[j].ny;
			z = surf->vertices[j].z + surf->vertices[j].nz;

			MRISsurfaceRASToVoxel(surf, images[0], x , y,z, &xw1,&yw1,&zw1);
			nx = xw1 - xw;
			ny = yw1 - yw;
			nz = zw1 - zw;
			float dist = sqrt(SQR(nx) + SQR(ny) + SQR(nz));
			if (FZERO(dist)) break;  
			nx /= dist;
			ny /= dist;
			nz /= dist;
		}


		double mag, val;

		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &xv,&yv,&zv);
		for(int t=-this->numberOfSteps; t<=this->numberOfSteps;t++)
		{
				x=xv +nx*(t)*step;
				y=yv +ny*(t)*step;
				z=zv +nz*(t)*step;

				MRIsampleVolume(images[0], x, y, z, &val);

				MRIsampleVolumeDerivativeScale(images[0], x, y, z, nx, ny, nz, &mag, this->gradientSigma); 

				magnitudes[t+numberOfSteps]+= mag;
				intensities[t+numberOfSteps]+=val;
							
		}
		int foundCSF=0;
		float in=0, out=0;
	
		for(int i=0;i<intensities.size()/2;i++)
		{
			//if(fabs(magnitudes[i])>=fabs(max_mag))
			if(intensities[i] > AVG_CSF|| intensities[i+this->numberOfSteps] >AVG_CSF)
			{
				foundCSF=1;
			}	
			in+= intensities[i];		
			out+= intensities[i+this->numberOfSteps];		
		}
		in=(int) in/ (this->numberOfSteps*100); //put it between 0 and 3
		out = (int)out / (this->numberOfSteps*100);
		if( in >3)
			in =3;
		if (out >3)
			out=3;
		
		std::tuple<int,int,int> tuple = std::make_tuple((int)in, (int)out, foundCSF);
		if(priors.count(tuple) >0)
			priors[tuple] +=1;
		else
			priors[tuple] =1;
	}
	return priors;
}

float MRIS_MultimodalRefinement::GetVentricleIntensity()
{
	float intensity=0;
	int count=0;
	for (int x = 0 ; x < this->images[0]->width ; x++)
	{
		for (int y = 0 ; y < this->images[0]->height ; y++)
		{
			for (int z = 0 ; z < this->images[0]->depth ; z++)
			{

				int label = MRIgetVoxVal(this->segmentation, x, y, z, 0) ;
				if(label==4 ||label==43)
				{
					intensity+=MRIgetVoxVal(this->images[0], x,y,z,0);
					count ++;
				}
			}			
		}
	}
	return .7*intensity/count;
}
void MRIS_MultimodalRefinement::getTarget2(MRIS* surf )
{
	std::cout << " getTarget "<< std::endl;
	std::cout << " number of steps "<< this->numberOfSteps <<std::endl;
	std::cout << " step "<< step << std::endl;
	std::cout <<" gradient sigma " << this->gradientSigma << std::endl;

	double x,y,z, xv, yv, zv; //xvp, yvp, zvp, xvn, yvn, zvn;//, nxv, nyv,nzv;
	double nx,ny,nz;
		
	for (unsigned j=0;j<surf->nvertices;j++)
	{

		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
	
		std::vector<float> intensities(numberOfSteps*2+1,0);
		std::vector<float> magnitudes(numberOfSteps*2+1,0);
		nx =  surf->vertices[j].nx;
		ny =  surf->vertices[j].ny;
		nz =  surf->vertices[j].nz;
		float dist = sqrt(nx*nx + ny*ny + nz*nz);
		nx /= dist;
		ny /= dist;
		nz /= dist;

		double mag, val;
		for(int t=-this->numberOfSteps; t<=this->numberOfSteps;t++)
		{
			x=surf->vertices[j].x +nx*t*step;
			y=surf->vertices[j].y +ny*t*step;
			z=surf->vertices[j].z +nz*t*step;

			MRISsurfaceRASToVoxel(surf, images[0],x,y,z, &xv,&yv,&zv);
			MRIsampleVolume(images[0], xv, yv, zv, &val);

			double next_val, prev_val;
			float norm=0, valn=0, valp=0;
			for (int i = 1; i<3 ;i++) 
			{
				float w =   exp(-pow((i-1)*step,2) / (2 * pow(this->gradientSigma ,2)));
				norm += w;
				x= surf->vertices[j].x +nx*(t+i)*step;
				y=surf->vertices[j].y +ny*(t+i)*step;
				z=surf->vertices[j].z +nz*(t+i)*step;

				MRISsurfaceRASToVoxel(surf, images[0],x,y,z, &xv,&yv,&zv);
				MRIsampleVolume(images[0], xv, yv, zv, &next_val);
				valn += w * next_val;

				x= surf->vertices[j].x +nx*(t-i)*step;
				y=surf->vertices[j].y +ny*(t-i)*step;
				z=surf->vertices[j].z +nz*(t-i)*step;

				MRISsurfaceRASToVoxel(surf, images[0],x,y,z, &xv,&yv,&zv);
				MRIsampleVolume(images[0], xv, yv, zv, &prev_val);
				valp += w * prev_val;
				if(vertexDebug==j)
					std::cout << "prev  " << prev_val <<  " next " << next_val << std::endl;
			}
			valn /= norm;
			valp /= norm;
			mag = (valn - valp) ; 
			magnitudes[t+numberOfSteps]+= mag;
			intensities[t+numberOfSteps]+=val;
			if(vertexDebug==j)
				std::cout << "distance from ori " << t*step <<  " " << mag <<   " " << val << std::endl;

		}
		float opt_mag=0,opt_val;
		int opt_t;
		for(int i=0;i<intensities.size();i++)
		{
			//if(fabs(magnitudes[i])>=fabs(max_mag))
			
			if((magnitudes[i]>opt_mag && this->maximumGradient )||
				(magnitudes[i]<opt_mag && !this->maximumGradient ))
			{
				x=surf->vertices[j].x  +nx*(i -numberOfSteps)*step;
				y= surf->vertices[j].y +ny*(i-numberOfSteps)*step;
				z=surf->vertices[j].z +nz*(i-numberOfSteps)*step;

				//MRIsampleVolumeType(segmentation, xv, yv, zv, &val, SAMPLE_NEAREST);
				//MRISsurfaceRASToVoxel(surf, images[0],x,y,z, &xv,&yv,&zv);
				//MRIsampleVolume(images[0], xv, yv, zv, &val);
	
				//if(val <.5)
				{
					opt_mag = magnitudes[i];
					opt_val= intensities[i];
					opt_t =  i;
					//double xt,yt,zt;
					//MRIvoxelToSurfaceRAS(images[0], x,y,z, &xt,&yt,&zt);
					surf->vertices[j].targx = x ;
					surf->vertices[j].targy = y;
					surf->vertices[j].targz = z;

					if(vertexDebug==j)
					{
						std::cout << "TAKING: distance from ori " << (i-numberOfSteps)*step <<  " " << opt_mag <<   " " << opt_val << std::endl;

					}
				}
			}		
		}
	}
}
void MRIS_MultimodalRefinement::getTargetAdHoc(MRIS* surf )
{
	std::cout << " getTarget "<< std::endl;
	std::cout << " number of steps "<< this->numberOfSteps <<std::endl;
	std::cout << " step "<< step << std::endl;
	std::cout <<" gradient sigma " << this->gradientSigma << std::endl;
	std::cout <<" T2_max  " << this->T2_max<< std::endl;
	std::cout <<" T2_min " << this->T2_min << std::endl;
	std::cout <<" MAX_CSF " << this->MAX_CSF << std::endl;


	for (unsigned j=0;j<surf->nvertices;j++)
	{
		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
	}
	for (unsigned j=0;j<surf->nvertices;j++)
	{
		double nx,ny,nz;
		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].nx,surf->vertices[j].ny,surf->vertices[j].nz, &nx,&ny,&nz);
		float dist = sqrt(nx*nx + ny*ny + nz*nz);
		if( dist>0)
		{
			nx /= dist;
			ny /= dist;
			nz /= dist;
		}	

		double x,y,z, xv, yv, zv, xvp, yvp, zvp, xvn, yvn, zvn;
		double max_mag=0;
		double max_val=0;
		double first_val=0;
		double max_t=10000;

		bool csf_found=false;
		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &xv,&yv,&zv);
		MRIsampleVolume(images[0], xv,yv, zv,   &first_val);
		if(vertexDebug==j)
		{
			std::cout << "first intensity" << first_val << " " << surf->vertices[j].x <<" "<< surf->vertices[j].y << " "<< surf->vertices[j].z<< std::endl;

		}
		for (int d=-1; d<2;d+=2)
		{
			nx = -d*nx ;
			ny =-d*ny;
			nz =-d*nz;
			bool keep_going=true;
			//double mag,prev_mag, next_mag, val;
			for(int t=0; t<numberOfSteps&&keep_going;t++)
			{

				x=surf->vertices[j].x +nx*(t-1)*step;
				y=surf->vertices[j].y +ny*(t-1)*step;
				z=surf->vertices[j].z +nz*(t-1)*step;

				MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xvp,&yvp,&zvp);

				x=surf->vertices[j].x +nx*(t+1)*step;
				y=surf->vertices[j].y +ny*(t+1)*step;
				z=surf->vertices[j].z +nz*(t+1)*step;

				MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xvn,&yvn,&zvn);

				x=surf->vertices[j].x  +nx*t*step;
				y=surf->vertices[j].y  +ny*t*step;
				z=surf->vertices[j].z +nz*t*step;

				MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xv,&yv,&zv);


				for(unsigned int k=0;k<images.size();k++)
				{

					double mag, next_mag, prev_mag , val, prev_val, next_val;

					MRIsampleVolume(images[k], xv, yv, zv, &val);
					MRIsampleVolume(images[k], xvn,yvn,zvn, &next_val);
					MRIsampleVolume(images[k],xvp,yvp,zvp, &prev_val);


					MRIsampleVolumeDerivativeScale(images[k], xv, yv, zv, nx, ny, nz, &mag, this->gradientSigma);   // expensive
					MRIsampleVolumeDerivativeScale(images[k], xvn,yvn,zvn, nx, ny, nz, &next_mag, this->gradientSigma);   // expensive
					MRIsampleVolumeDerivativeScale(images[k], xvp,yvp,zvp, nx, ny, nz, &prev_mag, this->gradientSigma);   // expensive
					bool take=false;
					if ( continuousIncrease(prev_val, val, next_val) && score(val, max_val, mag, max_mag, t, max_t) >1 && isGrey(val) && isPial(next_val) &&( !csf_found)) 
					{
						take=true;
						if(vertexDebug==j)
						{
							std::cout << "Found bright " << t*step << std::endl;
						}
					}
					if ( continuousDecrease(prev_val ,val, next_val)  && score(val, max_val, mag, max_mag, t, max_t) >1 && isGrey(val)  && isPial(prev_val)  &&( !csf_found ) ) 
					{
						take=true;
						if(vertexDebug==j)
						{
							std::cout << "Found bright backwards " << t*step << std::endl;
						}
					}
					if (isGrey(val) && !isGrey(prev_val ) && !isGrey(next_val))
					{
						take=false;
					}

					if(( isPial(first_val) && isGrey(val) && (isPial(prev_val) ||isPial( next_val) )))
					{		
						if (t < max_t)
						{
							take=true;
						}
						else
						{
							take=false;
						}

						if(vertexDebug==j)
						{
							std::cout << "first val max intensity test  " << t*step << std::endl;
						}
					}

					if ( continuousDecrease(prev_val, val, next_val) && fabs(mag)>=max_mag  && isCSF(next_val)  && (isGrey(val)  || isGrey(prev_val ))) 
					{
						csf_found=true;
						take=true;
						if(vertexDebug==j)
						{
							std::cout << "Found CSF: " << t*step << std::endl;
						}
					}
					if(take ) //||( t==0 && d==-1))
					{
						max_mag = fabs(mag);
						max_val= val;
						max_t =  t;
						surf->vertices[j].targx =x ;
						surf->vertices[j].targy = y;
						surf->vertices[j].targz = z;

						if(vertexDebug==j)
						{
							std::cout << "TAKING: distance from ori " << t*step*d << std::endl;


						}	

						//keep_going=false;
					}
					if(isCSF(val) && continuousDecrease(prev_val ,val, next_val) )
					{
						t=numberOfSteps;
						if(vertexDebug==j)
						{
							std::cout << "Breaking due to reach CSF , this leaving sulcus" << t*step*d << std::endl;

						}
						keep_going=false;
					}
					if(isPial(val)  && isGrey(prev_val) )
					{
						if(vertexDebug==j)
						{
							std::cout << "Breaking due to reach max intensity , this leaving sulcus" << t*step*d << std::endl;

						}
						keep_going=false;
						t=numberOfSteps;
					}

					if(vertexDebug==j)
					{
						std::cout << "distance  " << t*step*d << std::endl;
						std::cout << "coordinates " << x << ", " << y << ", " <<z << std::endl;
						std::cout << "Intensity previous , val, next: " <<prev_val << " " << val << " " <<next_val << std::endl;
						std::cout << "Gradient previous, mag, next: " << prev_mag << " " << mag << next_mag << std::endl;

					}

				}
			}
		}
	}
}

void MRIS_MultimodalRefinement::refine (MRIS* surf, float lambda, int iterations)
{
	double max_thickness = 60;
	float step=0.4; //images[0].xsize/2.0; 
	//	double max_csf=55; //need?
	double T2_max=  160;
	//double T2_mid = 150;
	double T2_min=100;
	//	double T2_max=160;
	//	double T2_min=110;
	for (unsigned j=0;j<surf->nvertices;j++)
	{
		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
	}
	for (unsigned j=0;j<surf->nvertices;j++)
	{
		double nx,ny,nz;
		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].nx,surf->vertices[j].ny,surf->vertices[j].nz, &nx,&ny,&nz);
		float dist = sqrt(nx*nx + ny*ny + nz*nz);
		if( dist>0)
		{
			nx /= dist;
			ny /= dist;
			nz /= dist;
		}	

		double x,y,z, xv, yv, zv, xvp, yvp, zvp, xvn, yvn, zvn;
		for(int t=0; t<max_thickness*2;t++)
		{

			x=surf->vertices[j].x - nx*max_thickness*step +nx*(t-1)*step;
			y=surf->vertices[j].y - ny*max_thickness*step +ny*(t-1)*step;
			z=surf->vertices[j].z - nz*max_thickness*step +nz*(t-1)*step;

			MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xvp,&yvp,&zvp);

			x=surf->vertices[j].x - nx*max_thickness*step +nx*(t+1)*step;
			y=surf->vertices[j].y - ny*max_thickness*step +ny*(t+1)*step;
			z=surf->vertices[j].z - nz*max_thickness*step +nz*(t+1)*step;

			MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xvn,&yvn,&zvn);

			x=surf->vertices[j].x - nx*max_thickness*step +nx*t*step;
			y=surf->vertices[j].y - ny*max_thickness*step +ny*t*step;
			z=surf->vertices[j].z - nz*max_thickness*step+nz*t*step;

			MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xv,&yv,&zv);

			double max_mag=0;

			for(unsigned k=0;k<images.size();k++)
			{

				double mag, next_mag, prev_mag , val, prev_val, next_val;

				MRIsampleVolume(images[k], xv, yv, zv, &val);
				MRIsampleVolume(images[k], xvn,yvn,zvn, &next_val);
				MRIsampleVolume(images[k],xvp,yvp,zvp, &prev_val);

				if(val<prev_val)
				{
					nx *=-1;
					ny *=-1;
					nz *=-1;
				}


				MRIsampleVolumeDerivativeScale(images[k], xv, yv, zv, -nx, -ny, -nz, &mag, 1.0);   // expensive
				MRIsampleVolumeDerivativeScale(images[k], xvn,yvn,zvn, -nx, -ny, -nz, &next_mag, 1.0);   // expensive
				MRIsampleVolumeDerivativeScale(images[k], xvp,yvp,zvp, -nx, -ny, -nz, &prev_mag, 1.0);   // expensive
				//std::cout << mag << " " << prev_mag  << "  " << next_mag << " val "<< val <<std::endl;

				if ((fabs(mag) > fabs(prev_mag)) && (fabs(mag) > fabs(next_mag)) &&  (val <= T2_max) && (val >= T2_min) && fabs(mag)*1.2> max_mag && prev_val< val && val <next_val) 
				{
					max_mag = fabs(mag);
					std::cout << mag << std::endl;						
					/*surf->vertices[j].targx +=x ;
					  surf->vertices[j].targy += y;
					  surf->vertices[j].targz += z;
					  */

					surf->vertices[j].x = x;
					surf->vertices[j].y = surf->vertices[j].y -ny*max_thickness*step+t*step*ny ;
					surf->vertices[j].z = surf->vertices[j].z -nz*max_thickness*step+t*step*nz ;
				}
			}
		}
	}
}
