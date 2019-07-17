

class MRIS_MultimodalRefinement {
	public:
		void refine(MRIS* surface, float lambda, int iterations);
		void  getTarget(MRIS* surface);
		std::vector<MRI*> images;
		void addImage(MRI* image)
		{
			images.push_back(image);
		}

		int vertexDebug= -1;
		float step=0.4;
		int numberOfSteps =20;
		double T2_max=180;
		double T2_min=130;
};
void MRIS_MultimodalRefinement::getTarget(MRIS* surf )
{
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
		for (int d=-1; d<2;d+=2)
		{
			//std::cout << d <<std::endl;
			nx = -d*nx ;
			ny =-d*ny;
			nz =-d*nz;
			//double mag,prev_mag, next_mag, val;
			for(int t=0; t<numberOfSteps;t++)
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


				for(unsigned k=0;k<images.size();k++)
				{

					double mag, next_mag, prev_mag , val, prev_val, next_val;

					MRIsampleVolume(images[k], xv, yv, zv, &val);
					MRIsampleVolume(images[k], xvn,yvn,zvn, &next_val);
					MRIsampleVolume(images[k],xvp,yvp,zvp, &prev_val);

					if(val<prev_val && prev_val > next_val)
					{
						t=numberOfSteps;
						if(vertexDebug==j)
						{
							std::cout << "Breaking due to intensity decreasing.  Max thickness: " << numberOfSteps*step << std::endl;
						}
						break;
					}


					MRIsampleVolumeDerivativeScale(images[k], xv, yv, zv, -nx, -ny, -nz, &mag, 1.0);   // expensive
					MRIsampleVolumeDerivativeScale(images[k], xvn,yvn,zvn, -nx, -ny, -nz, &next_mag, 1.0);   // expensive
					MRIsampleVolumeDerivativeScale(images[k], xvp,yvp,zvp, -nx, -ny, -nz, &prev_mag, 1.0);   // expensive
					//std::cout << mag << " " << prev_mag  << "  " << next_mag << " val "<< val <<std::endl;
					
					if ((fabs(mag) > fabs(prev_mag) || prev_val > T2_min) && (fabs(mag) > fabs(next_mag) || next_val < T2_max ) &&  (val <= T2_max) && (val >= T2_min) && (fabs(mag) > max_mag && val > max_val) ) // && prev_val< val && val <next_val) 
					{
						max_mag = fabs(mag);
						max_val= val;
						//std::cout << mag << std::endl;						
						surf->vertices[j].targx =x ;
						surf->vertices[j].targy = y;
						surf->vertices[j].targz = z;
					
						if(vertexDebug==j)
						{
							std::cout << "Found new maximum. coordinate: " << x << " " << y << " " << z << std::endl;
							std::cout << "Intensity val, next, and previous" << val << " " <<next_val << " " <<prev_val<< std::endl;
							std::cout << "Gradient mag, next, and previous" << mag << " " <<next_mag << " " <<prev_mag<< std::endl;

				
						}	
						
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
