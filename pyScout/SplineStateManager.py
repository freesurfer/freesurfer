
import sys
import pickle

import tkFileDialog

import Globals
from SplineSet import IndexedSpline

def SaveState():
    """ saves all the concerned variables in the
    application """

    # open file for writing
    fileName = tkFileDialog.asksaveasfilename(defaultextension=".ssv",
                                              filetypes=[ ( 'Spline SAVE',
                                                            '*.ssv') ] )
    if not fileName: return

    handle = open(fileName, 'w')

    ### 1. save the rendered properties
    pickle.dump( Globals.renderProps, handle )

    # 2. save the slice index and the orientation
    pickle.dump( Globals.imagePipeline.GetSliceIndex(),
                 handle )
    pickle.dump( Globals.imagePipeline.GetOrientation(),
                 handle )

    # 3. save the splineSet

    # save the labels
    pickle.dump(Globals.objectSet.labels, handle)
    
    splines = Globals.objectSet.GetAll()
    pickle.dump( len(splines), handle )
    for s in splines:
        pickle.dump( s.globalIndex,  handle )
        pickle.dump( s.sliceIndex,   handle )
        pickle.dump( s.orientation,  handle )
        pickle.dump( s.label,        handle )

        pickle.dump( s.points, handle )

    handle.close()

def LoadState(sliceCallback,orientationCallback):
    """ loads all the concerned variables
    and restores the state of an application """

    # open file for reading
    handle = tkFileDialog.askopenfile( filetypes= [ ( 'Spline SAVE files',
                                                      '*.ssv' ) ] )
    if not handle: return

    ### 1. read the rendered properties
    Globals.renderProps = pickle.load(handle)

    ### 2. load the slice index and the orientation
    pipelineSliceIndex   = pickle.load( handle )
    pipelineOrientation  = pickle.load(handle)

    ### 3. load the splineSet
    Globals.objectSet.labels = pickle.load( handle )

    numberOfSplines = pickle.load( handle )
    for i in range(numberOfSplines):
        globalIndex = pickle.load(handle)
        sliceIndex  = pickle.load(handle)
        orientation = pickle.load(handle)
        label       = pickle.load(handle)

        points      = pickle.load(handle)

        spline = IndexedSpline(1)
        spline.sliceIndex  = sliceIndex
        spline.orientation = orientation
        spline.SetPoints( points )
        # hide all regions until updating
        spline.Hide()
        Globals.objectSet.AddSpline(spline, label,
                                    globalIndex)

    orientationCallback(pipelineOrientation)
    sliceCallback(pipelineSliceIndex)
    Globals.renWin.Render()
    
        
        
