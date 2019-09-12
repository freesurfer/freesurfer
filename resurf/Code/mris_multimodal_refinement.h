

class MRIS_MultimodalRefinement {
	public:
		void refine(MRIS* surface, float lambda, int iterations);
		void  getTarget(MRIS* surface);
		void  getTargetAdHoc(MRIS* surface);
		std::vector<MRI*> images;
		MRI* segmentation;

		void addImage(MRI* image)
		{
			images.push_back(image);
		}
		void SetVertexDebug(int i){ this->vertexDebug =i;}
		void SetStep(float i){ this->step =i;}
		void SetNumberOfSteps(int i){ this->numberOfSteps =i;}
		void SetGradientSigma(float i){ this->gradientSigma =i;}
		void SetSegmentation(MRI* s) { this->segmentation=s;}

	private:
		int vertexDebug= -1;
		float step=0.4;
		int numberOfSteps =20;
		double T2_max=200;
		double T2_min_gray=120;
		double T2_min=120;
		double MAX_CSF=80;
		float gradientSigma=0.25;

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
		std::vector<float> intensities(numberOfSteps*2,0);
		std::vector<float> magnitudes(numberOfSteps*2,0);

		MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].nx,surf->vertices[j].ny,surf->vertices[j].nz, &nx,&ny,&nz);

		float dist = sqrt(nx*nx + ny*ny + nz*nz);
		if( dist>0)
		{
			nx /= dist;
			ny /= dist;
			nz /= dist;
		}	

		double x,y,z, xv, yv, zv;

		double mag, val;
		for(int t=-numberOfSteps; t<numberOfSteps;t++)
		{

			x=surf->vertices[j].x  +nx*t*step;
			y=surf->vertices[j].y  +ny*t*step;
			z=surf->vertices[j].z +nz*t*step;

			MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xv,&yv,&zv);

			for(int k=0;k<images.size();k++)
			{
				MRIsampleVolume(images[k], xv, yv, zv, &val);
				MRIsampleVolumeDerivativeScale(images[k], xv, yv, zv, nx, ny, nz, &mag, this->gradientSigma);   // expensive

				magnitudes[t+numberOfSteps]+=fabs(mag);
				intensities[t+numberOfSteps]+=val;
			}
							
		}
 		float max_mag=0,max_val;
		int max_t;
		for(int i=0;i<intensities.size();i++)
		{
			//if(fabs(magnitudes[i])>=fabs(max_mag))
			if(magnitudes[i]>=max_mag)
			{
				x=surf->vertices[j].x  +nx*(i -numberOfSteps)*step;
				y=surf->vertices[j].y  +ny*(i-numberOfSteps)*step;
				z=surf->vertices[j].z +nz*(i-numberOfSteps)*step;

				MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xv,&yv,&zv);
				MRIsampleVolumeType(segmentation, xv, yv, zv, &val, SAMPLE_NEAREST);
	
				if(val != 165 && val != 258 && val != 259)
				{
					max_mag = magnitudes[i];
					max_val= intensities[i];
					max_t =  i;
					surf->vertices[j].targx =x ;
					surf->vertices[j].targy = y;
					surf->vertices[j].targz = z;

					if(vertexDebug==j)
					{
						std::cout << "TAKING: distance from ori " << i*step << std::endl;

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
