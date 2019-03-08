#ifndef _HierarchicalClusteringPruner_txx_
#define _HierarchicalClusteringPruner_txx_

#include "HierarchicalClusteringPruner.h"
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include "TrkVTKPolyDataFilter.txx"


template < class TOutputMesh, class TImage>
	HierarchicalClusteringPruner< TOutputMesh,TImage>
::HierarchicalClusteringPruner()
{
	this->SetNumberOfRequiredInputs (0);
}


template < class TOutputMesh, class TImage>
void
	HierarchicalClusteringPruner< TOutputMesh,TImage>
::GenerateData()
{
	this->m_outputBundles.clear();
	this->m_clustersIds.clear();
//	std::cout << " Hierarchy file : " <<  this->m_hierarchyFilename << std::endl;
	ifstream file (  this->m_hierarchyFilename); 
	std::string value;
	getline ( file, value, ',' ); 
	getline ( file, value, ',' ); 
	std::map<long long,long long> idTree_idCluster;
	std::map<long long,std::vector<long long>> idNode_idChilds;
	std::map<long long,long long>  idCluster_idClusterParent;
	std::set<long long> clustersIds;
	std::set<long long> leafIds;
	//int max = 0;
	std::map<long long,std::string> fiberFiles;
	while ( file.good() )
	{

		getline ( file, value, ',' );
		long long a  = atoll(value.c_str());
		//if(a>max)
		//	max = a;
	//}
	
	//while(max >-1)
	//{
		idNode_idChilds[a] = std::vector<long long>();
		//max--;
		//std::cout << " im looping here lalala"<< std::endl;
	}
	file.close();
	file.open ( this->m_hierarchyFilename); 
	getline ( file, value, ',' ); 
	getline ( file, value, ',' ); 
	//std::cout << "numbe rof clusters " <<this->m_numberOfClusters<<std::endl;
	while ( file.good() )
	{
		getline ( file, value, ',' );
		long long v1 = atoll(value.c_str());
		getline ( file, value, ',' ); 
		long long v2 = atoll(value.c_str());
		if(clustersIds.size() < this->m_numberOfClusters && v2 != 0)
		{
			//std::cout << v1 <<" "<< v2 << " " << clustersIds.size()<< " " << clustersIds.count(v1) <<std::endl;
			clustersIds.erase(v1);
			clustersIds.insert(v2);
		}
		idCluster_idClusterParent[v2] = v1;
		leafIds.insert(v2);
		leafIds.erase(v1);
//		std::cout << " v1 v2 " << v1  << " "<< v2 << std::endl;
	}  
//	std::cout << "num clusters found " << clustersIds.size() << std::endl;
//	std::cout << "num leaft found " << leafIds.size() << std::endl;
	for(std::set<long long>::iterator it= leafIds.begin();it!= leafIds.end();it++)
	{	 
		std::string a = (this->m_fiberFormat==FiberFormat::VTK)?"vtk":"trk";
//		std::cout << this->m_clustersPath << std::endl;
//		std::cout << a << std::endl;
		long long b = *it;
		std::string num = std::to_string(b);
//		std::cout << num<< std::endl;
		std::string clusterFile = this->m_clustersPath + "/" +num+"."+a;
		ifstream ifile(clusterFile);
		if (ifile)
		{	
			fiberFiles[*it] =clusterFile;
			long long lastNumber = *it;
			idNode_idChilds[lastNumber].push_back( lastNumber);
			while(idCluster_idClusterParent[lastNumber]!=0)
			{	
				lastNumber = idCluster_idClusterParent[lastNumber];
				idNode_idChilds[lastNumber].push_back( *it);
			}

			idNode_idChilds[0].push_back(*it);
		}
		else
		{
			std::cout << clusterFile << std::endl;		
		}	
	}
//	int w=0;
	for(std::set<long long>::iterator it= clustersIds.begin();it!= clustersIds.end();it++)
	{
		std::vector<vtkSmartPointer<vtkPolyData>> polydatas;
		//std::cout << *it << " hola" << std::endl;	
	
		for(int j=0;j<idNode_idChilds[*it].size();j++)
		{
			std::string file = fiberFiles[idNode_idChilds[*it][j]];			
			if(  std::string(file).find(std::string(".trk")) !=std::string::npos)
			{
				itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
				trkReader->SetTrkFileName(file.c_str());
				trkReader->SetReferenceImage(this->m_referenceImage);
				trkReader->TrkToVTK();
				polydatas.push_back( trkReader->GetOutputPolyData() );
			}
			else
			{
				vtkSmartPointer<vtkPolyDataReader> reader = vtkPolyDataReader::New();
				reader->SetFileName (file.c_str() );
#if VTK_MAJOR_VERSION > 5
				reader->Update();
#else
				reader->GetOutput()->Update();
#endif
				polydatas.push_back(reader->GetOutput());
			}
		}
//		std::cout << w << polydatas.size() << std::endl;
		AppendBundleFilter::Pointer appendBundles = AppendBundleFilter::New();
		appendBundles->SetNumberOfColours(polydatas.size());
		appendBundles->SetRepresentatives(false);
		appendBundles->SetInput(polydatas);
		appendBundles->Update();
		//std::cout << polydatas.size()<< std::endl;
		this->m_clustersIds.push_back(*it);
		this->m_outputBundles.push_back(appendBundles->GetOutput());
	}

}

#endif
