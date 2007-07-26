
import sys
import pickle

import tkFileDialog

import Region
import Globals
import Spline

def SavePreferences():
    """ saves the preferences of the application """

    # open file for writing
    fileName = tkFileDialog.asksaveasfilename(defaultextension=".pref",
                                              filetypes=[ ( 'PREF files',
                                                            '*.pref') ] )
    if not fileName: return

    handle = open(fileName,'w')

    # write a version int for future reference
    pickle.dump(1,handle)
    
    pickle.dump(Globals.preferences,handle)

    handle.close()

def LoadPreferences():
    """ loads the preferences of the application """

    handle = tkFileDialog.askopenfile(filetypes=[ ( 'PREF files',
                                                    '*.pref') ] )
    if not handle: return

    # read the version tag
    n = pickle.load(handle)
    
    Globals.preferences = pickle.load(handle)

    handle.close()
    

def SaveState():
    """ saves all the concerned variables
    of the application """
    
    # open file for writing
    fileName = tkFileDialog.asksaveasfilename(defaultextension=".save",
                                              filetypes=[ ( 'SAVE files',
                                                            '*.save') ] )
    if not fileName: return
    
    handle = open(fileName,'w')    

    ### 1. save the rendered properties
    pickle.dump( Globals.renderProps, handle)
    
    # 2. save the slice index and the orientation
    pickle.dump( Globals.imagePipeline.GetSliceIndex(),
                 handle)
    pickle.dump( Globals.imagePipeline.GetOrientation(),
                 handle)
    
    # 3. save the regionSet, i.e. for each region - save zero and one
    regions = Globals.objectSet.GetAllRegions()
    pickle.dump( len(regions), handle )
    for r in regions:
        pickle.dump( r.index, handle)
        pickle.dump( r.sliceIndex, handle)
        pickle.dump( r.orientation, handle)
        pickle.dump( r.zero.points, handle)
        pickle.dump( r.one.points,  handle)
        pickle.dump( len(r.profiles), handle)
        for p in r.profiles:
            pickle.dump( p.points, handle)
            
            
    handle.close()
            

def LoadState(sliceCallback,orientationCallback):
    """ loads all the concerned variables
    and restores the state of an application """
    
    # open file for reading
    handle = tkFileDialog.askopenfile(filetypes=[ ( 'SAVE files',
                                                      '*.save') ] )
    if not handle: return

    ### 1. read the rendered properties
    Globals.renderProps = pickle.load(handle)
    
    # 2. load the slice index and the orientation
    pipelineSliceIndex = pickle.load(handle)
    pipelineOrientation = pickle.load(handle)
    
    # 3. load the region set
    Globals.objectSet.dic.clear()
    # load the number of regions
    numberOfRegions = pickle.load(handle)
    # load each region
    for i in range(numberOfRegions):
        regionIndex = pickle.load(handle)
        sliceIndex  = pickle.load(handle)
        orientation = pickle.load(handle)
        ptsZero     = pickle.load(handle)
        ptsOne      = pickle.load(handle)
        numberOfProfiles = pickle.load(handle)
        profiles = []
        for j in range(numberOfProfiles):
            p = pickle.load(handle)
            s = Spline.PlaneSpline(p)
            profiles.append(s)
        # commit each region
        zero = Spline.PlaneSpline(ptsZero)
        one  = Spline.PlaneSpline(ptsOne)
        region = Region.Region()
        region.SetZeroPoints(zero)
        region.SetOnePoints(one)
        region.sliceIndex = sliceIndex
        region.orientation = orientation
        region.profiles = profiles
        region.rendering.UpdateProfiles(region.profiles)
        Globals.objectSet.AddRegion( region, regionIndex )
        region.Hide()

    handle.close()
    orientationCallback(pipelineOrientation)
    sliceCallback(pipelineSliceIndex)
    Globals.renWin.Render()
