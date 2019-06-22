class MRIS_MultimodalRefinement {

  public:
	void refine(MRIS* surface, float lambda, int iterations);
	std::vector<MRI*> images;
	void addImage(MRI* image)
	{
		images.push_back(image);
	}
};

void MRIS_MultimodalRefinement::refine (MRIS* surf, float lambda, int iterations)
{
	for (unsigned j=0;j<surf->nvertices;j++)
	{
		surf->vertices[j].targx = surf->vertices[j].x;
		surf->vertices[j].targy = surf->vertices[j].y;
		surf->vertices[j].targz = surf->vertices[j].z;
	}
	for(int i=0;i<iterations;i++)
	{ 
		for (unsigned j=0;j<surf->nvertices;j++)
		{

			for(unsigned k=0;k<images.size();k++)
			{
				double x,y,z;
				float pdx, pdy, pdz;

				if(surf->vertices[j].x >1)
				{
					MRISsurfaceRASToVoxel(surf, images[k], surf->vertices[j].targx,surf->vertices[j].targy,surf->vertices[j].targz, &x,&y,&z);
					float magnitud = MRIvoxelGradient(images[k], (float) x, (float) y,(float) z, &pdx,  &pdy, &pdz);	
					float norm= magnitud/images.size();
					if(norm>1)
					{
						//						std::cout << magnitud << " "<< pdx << " "<< pdy << " " << pdz << std::endl;
						surf->vertices[j].targx += lambda *  pdx /norm;
						surf->vertices[j].targy += lambda *  pdy/norm;
						surf->vertices[j].targz += lambda *  pdz/norm;
					}
				}


			}
		}

	}

}

