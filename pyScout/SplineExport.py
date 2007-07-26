
import tkFileDialog

import vtk
import ImagePipeline
import SplineSet
from Spline import PlaneSpline
import Globals


def SplineExport():
    selector = ImagePipeline.SliceSelector(Globals.imagePipeline.volume)
    
    # open file for writing
    fileName = tkFileDialog.asksaveasfilename(defaultextension=".csv",
                                              filetypes=[ ( 'CSV files',
                                                            '*.csv' ) ] )
    if not fileName: return
    
    handle = open(fileName, 'w')

    splines = Globals.objectSet.GetAll()

    for s in splines:
        selector.SetOrientation(s.orientation)
        selector.sliceIndex = s.sliceIndex
        
        bufSpline = PlaneSpline(
            s.Resample( Globals.preferences["numberOfPoints"] ) )
        strLine = str( [s.globalIndex,s.label,s.sliceIndex,s.orientation, \
                        len(bufSpline.points),bufSpline.GetApproximateLength(),
                        bufSpline.GetApproximateLength()/ float(len(bufSpline.points)) ] )[1:-1] + ',' \
                        + bufSpline.GetIntensities(selector) + '\n'
        handle.write(strLine)

    handle.close()
