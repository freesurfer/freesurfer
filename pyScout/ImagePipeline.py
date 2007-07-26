
import vtk
import sys
import math

class LUTBuilder:

    def __init__(self, inMin, inMax, inSigma):
        self.rangeMin = inMin
        self.rangeMax = inMax
        self.sigma = inSigma

    def Build(self):
        lut = vtk.vtkLookupTable()
        lut.SetTableRange(self.rangeMin,
                          self.rangeMax)
        lut.SetSaturationRange(0,0)
        lut.SetHueRange(0,0)
        lut.SetValueRange(0,1)
        
        for i in range(256):
            v = math.pow( float(i)/256.0, self.sigma )
            lut.SetTableValue(i, v, v, v, 1.0)
        lut.Build()
        
        return lut
    

class SliceSelector:

    def __init__(self,vol):
        self.spacing = vol.GetSpacing()
        self.origin = vol.GetOrigin()

        extent = vol.GetWholeExtent()
        self.extent = []
        for a in extent: self.extent.append(a)

        self.InitDictionary()

        self.SetOrientation( "yz" )
        # if one of the dimensions has zero depth, display the view
        #    perpendicular to it
        i=-1
        for a in range(3):
            if self.extent[2*a]==self.extent[2*a+1]:
                i=a
        if i>-1:
            lstOri = [ "yz", "xz", "xy" ]
            self.SetOrientation( lstOri[i] )

        smin, smax = self.GetAxisExtent()
        self.sliceIndex = smin

    def InitDictionary(self):

        mxy = vtk.vtkMatrix4x4() # leave it as the origin

        myz = vtk.vtkMatrix4x4() 
        myz.Zero()
        myz.SetElement( 1,0, 1)
        myz.SetElement( 2,1, 1)
        myz.SetElement( 0,2, 1)

        mxz = vtk.vtkMatrix4x4()
        mxz.Zero()
        mxz.SetElement( 0,0, 1)
        mxz.SetElement( 1,2, 1)
        mxz.SetElement( 2,1, 1)

        self.dic = { "xy" : (2,mxy),
                     "xz" : (1,mxz),
                     "yz" : (0,myz)
                     }
        self.crtData = None
        self.invMatrix = vtk.vtkMatrix4x4()

    def SetOrientation(self,ori):
        self.crtData = self.dic[ori]
        self.orientation = ori
        
    def GetAxisExtent(self):
        i = self.crtData[0]
        return self.extent[2*i:2*(i+1)]
    def GetDirectionCosines(self):
        return self.crtData[1]
    def GetZCoordinate(self):
        """ return the active coordinate """
        i = self.crtData[0]
        return self.sliceIndex * self.spacing[i] + self.origin[i]
    def GetAxesOrigin(self):
        i = self.crtData[0]
        buf = []
        buf.extend(self.origin)
        buf[i] = self.GetZCoordinate()
        return buf

    def InferSliceFromDepth(self,d):
        i = self.crtData[0]
        return math.floor( (d-self.origin[i])
                           / (self.spacing[i]+.001) +.5)

class ImagePipeline:
    """ Manages volume data
    - encapsulates the pipeline ending in the actor being displayed
    - manages the selection of the orientation and the slice
    - manages the LUT for the contrast setup
    - is also a wrapper around the volume data (because of the display-world transformation)
    """

    def UpdateLUT(self):
        lut = self.lutBuilder.Build()
        self.colors.SetLookupTable( lut )

    def __init__(self, ren):
        self.volume = self.LoadVolume()
        self.ren = ren

        # next initialize the pipeline
        self.slicer = vtk.vtkImageReslice()
        self.slicer.SetInput( self.volume )
        self.slicer.SetOutputDimensionality(2)

        # the next filter provides a mechanism for slice selection
        self.selector = SliceSelector(self.volume)
        self.slicer.SetResliceAxes( self.selector.GetDirectionCosines() )
        self.slicer.SetResliceAxesOrigin( self.selector.GetAxesOrigin() )

        # setup link for adjusting the contrast of the image
        r = self.volume.GetScalarRange()
        self.lutBuilder = LUTBuilder(r[0],r[1],1)
        lut = self.lutBuilder.Build()

        self.colors = vtk.vtkImageMapToColors()
        self.colors.SetInputConnection( self.slicer.GetOutputPort() )
        self.colors.SetLookupTable( lut )

        self.actor = vtk.vtkImageActor()
        self.actor.SetInput( self.colors.GetOutput() )

    def SetSlice(self,si):
        """ bridge method to the slice selector """
        self.selector.sliceIndex = si
        self.slicer.SetResliceAxesOrigin( self.selector.GetAxesOrigin() )
    def GetAxisExtent(self):
        """  wrapper method """
        return self.selector.GetAxisExtent()
    def SetOrientation(self,ori):
        self.selector.SetOrientation(ori)
        
        self.slicer.SetResliceAxes( \
            self.selector.GetDirectionCosines() )
        self.SetSlice(0)
        self.ren.ResetCamera()
        self.ren.GetRenderWindow().Render()

    # the following are accessor methods
    def GetSliceIndex(self): return self.selector.sliceIndex
    def GetOrientation(self): return self.selector.orientation
    def GetWorldPoint(self,x,y):
        m = self.slicer.GetResliceAxes()
        return m.MultiplyPoint( [x,y,0,1] )

    def WorldPointToDisplayFrame(self,wp):
        m = vtk.vtkMatrix4x4()
        m.DeepCopy( self.slicer.GetResliceAxes() )
        m.Invert()
        framePoint = m.MultiplyPoint(wp)
        retOne = framePoint[0:2]
        retTwo = int( math.floor(
            framePoint[2] / self.selector.spacing[ 
            self.selector.crtData[0] ] ) )
        return (retOne, retTwo)
    def InferSliceFromDepth(self,d):
        return self.selector.InferSliceFromDepth(d)
    
    def GetSpacing(self): return self.volume.GetSpacing()

    def LoadVolume(self):
        """ loads the volume - this only happens once during the lifetime of the application """

        if len(sys.argv) < 2:
            print " Error - please provide a file name (.vtk)"
            exit

        fileName = sys.argv[1]

        extension = fileName.split('.')[-1]
        if extension == "vtk":
            reader = vtk.vtkStructuredPointsReader()
            reader.SetFileName( fileName )
            reader.Update()
        elif extension == "tif":
            reader = vtk.vtkTIFFReader()
            reader.SetFileName( fileName )
            reader.Update()
        else:
            print " Unknown file type ", extension

        return reader.GetOutput()
