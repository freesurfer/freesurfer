
import tkFileDialog

import vtk
import ImagePipeline
import Spline
import Region
import Globals

def ProfileExport():
    selector = ImagePipeline.SliceSelector(Globals.imagePipeline.volume)
    
    # open file for writing
    fileName = tkFileDialog.asksaveasfilename(defaultextension=".csv",
                                              filetypes=[ ( 'CSV files',
                                                            '*.csv' ) ] )
    if not fileName: return

    handle = open(fileName,'w')
    
    regions = Globals.objectSet.GetAllRegions()

    for r in regions:
        # initialize selector
        selector.SetOrientation(r.orientation)
        selector.sliceIndex = r.sliceIndex
        
        for i in range(len(r.profiles)):
            lstProps = [r.index,r.sliceIndex,r.orientation,i,
                        len(r.profiles[i].points),r.profiles[i].GetApproximateLength(),\
                        r.profiles[i].GetApproximateLength()/ float(len(r.profiles[i].points)) ]
            strLine = str(lstProps)[1:-1] + ','  \
                      + r.profiles[i].GetIntensities(selector) + '\n'
            handle.write(strLine)
        
    handle.close()
