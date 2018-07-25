#include "AppendBundleFilter.h"
typedef struct {
	unsigned char r;
	unsigned char g;
	unsigned char b;
} color_triplet;

void AppendBundleFilter::Update()
{

	//std::cout << "AppendBundleFilter::Update()" << std::endl;
	int i;
	int red, green, blue;
	color_triplet table[256] = {{0,0,0}}; /* Initialize to all black */

	const color_triplet black = {0,0,0};
	const color_triplet white = {255,255,255};

	table[0]   = white;
	table[255] = black;		i = 20; /* first 20 and last 20 are reserved */
	for (red = 0; red <= 255; red+= 51) {/* the six values of red */
		for (green = 0; green <= 255; green += 51) {
			for (blue = 0; blue <= 255; blue+= 51) {
				table[i].r = red;
				table[i].g = green;
				table[i].b = blue;
				++i;
			}
		}
	}
	table[0]   = white; 
	table[255] = black;

	allBundles = vtkPolyData::New();
	vtkPoints   *allPoints  = vtkPoints::New();
	vtkIntArray *allLabels  = vtkIntArray::New();
	vtkUnsignedCharArray* allColors = vtkUnsignedCharArray::New();

	allColors->SetNumberOfComponents (3);
	allLabels->SetNumberOfComponents (1);
	allLabels->SetName ("Labels");

	allBundles->Initialize();
	allBundles->Allocate();

	int currentLabel = 0;
	//double val[3];

	vtkSmartPointer<vtkIntArray> intArrayRepresentativesWeights  = vtkIntArray::New();
	intArrayRepresentativesWeights->SetName("WeightsPerFiber");
	for (unsigned i=0; i<bundleList.size(); i++)
	{

		vtkSmartPointer<vtkPolyData> bundle = bundleList[i];

		//vtkSmartPointer<vtkPolyData> *bundle = bundleList[i];
		vtkCellArray *lines = bundle->GetLines();
		lines->InitTraversal();

		//  int index =((int)( 100./(bundle->GetNumberOfLines()))*((i%colorNumber)%(150)))%201;
		int index =((int)47.*((i%colorNumber)%(150)))%197+5;
		index = (int)(13*(i%colorNumber))%150+65;
		//  std::cout << "index " << index << std::endl;
		unsigned char color[3] = { table[index].r, table[index].g,table[index].b};

		vtkIdType pointCount=0, *pointBuf=0;
		if(rep)
		{
			vtkFieldData *fieldData = bundle->GetFieldData();
			//this is not fine, i knoww 0 will be the representatives index, or weight of fibers
			//vtkIntArray *arrayCellData = (vtkIntArray*) fieldData->GetArray("RepresentativeIndex");
			vtkIntArray *arrayCellData = (vtkIntArray*) fieldData->GetArray(0);

			if(arrayCellData != NULL)
			{
				if(   arrayCellData->GetNumberOfTuples() >1 ) 
				{
					std::cout << "WARNING: There are more than one representative index in field data " << std::endl;
				}
				//        for(int i=0; i< arrayCellData->GetNumberOfTuples();i++)
				//        {
				int cellId = arrayCellData->GetValue(0);
				for(int k=0;k<=cellId;k++)
				{
					lines->GetNextCell(pointCount, pointBuf);
				}
				//            std::cout <<" rep " << cellId << std::endl;
				for (vtkIdType k=0; k<pointCount; k++)
				{
					pointBuf[k] = allPoints->InsertNextPoint ( bundle->GetPoint ( pointBuf[k] ) );
				}
				intArrayRepresentativesWeights->InsertNextValue(bundle->GetNumberOfLines());
				allBundles->InsertNextCell (VTK_POLY_LINE, pointCount, pointBuf);
				allLabels->InsertNextTuple1 ( currentLabel );
				allColors->InsertNextTuple3 ( color[0], color[1], color[2] );

				//        }
			} 
		} else
		{
			/*if( bundle->GetFieldData())
			  {	
			  std::cout << bundle->GetFieldData() << std::endl ;
			  vtkFieldData *fieldData = bundle->GetFieldData();
			//this is not fine, i knoww 0 will be the representatives index, or weight of fibers
			vtkIntArray *arrayCellData = (vtkIntArray*) fieldData->GetArray(0);

			if(arrayCellData != NULL &&  arrayCellData->GetNumberOfTuples() ==  bundle->GetNumberOfLines() ) 
			{
			//        std::cout << " adding weights to file " << std::endl;
			for (int i=0;i<arrayCellData->GetNumberOfTuples();i++)
			{
			intArrayRepresentativesWeights->InsertNextValue(arrayCellData->GetValue(i));
			}
			}       
			}
			*/
			while ( lines->GetNextCell(pointCount, pointBuf) )
			{
				for (vtkIdType k=0; k<pointCount; k++)
				{
					pointBuf[k] = allPoints->InsertNextPoint ( bundle->GetPoint ( pointBuf[k] ) );
				}

				allBundles->InsertNextCell (VTK_POLY_LINE, pointCount, pointBuf);
				allLabels->InsertNextTuple1 ( currentLabel );
				allColors->InsertNextTuple3 ( color[0], color[1], color[2] );
			}
		}
		currentLabel++;
	}
	vtkSmartPointer<vtkFieldData> fieldData = vtkFieldData::New();
	fieldData->AddArray(intArrayRepresentativesWeights);

	allBundles->SetFieldData(fieldData);
	allBundles->SetPoints ( allPoints );
	allBundles->GetCellData()->SetScalars ( allColors );
	allBundles->GetCellData()->AddArray ( allLabels );
	//  vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();


	/*
	   vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	//  writer->SetFileTypeToBinary();
	writer->SetInput (allBundles);
	writer->SetFileName(output);
	writer->Update();

	allBundles->Delete();
	allPoints->Delete();
	allColors->Delete();
	allLabels->Delete();
	writer->Delete();
	*/

}
