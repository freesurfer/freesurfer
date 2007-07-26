from CrossHair import CrossHairWidget
from Spline import *

from ImagePipeline import ImagePipeline

# add text label to the spline
# hold it in this file

def GetTkColor(color):
    if color==None: return '#000000'
    
    strRed = hex(color[0])[2:]
    strGreen = hex(color[1])[2:]
    strBlue = hex(color[2])[2:]
    return '#' + strRed + strGreen + strBlue

class IndexedSpline(PlaneSpline):

    def __init__(self,showRendering):
        self.points = []
        
        self.globalIndex = None
        self.sliceIndex = None
        self.orientation = None
        self.label = None

        self.rendering = None
        if showRendering: self.AddRendering()

    def AddRendering(self):
        if self.rendering: return

        self.rendering = RenderSpline(1)
        self.rendering.SetInput( self )
        
    def AddPoint(self,wp):
        PlaneSpline.AddPoint(self,wp)
        if self.rendering: self.rendering.SetInput(self)

    def SetPoints(self,points):
        self.points = points
        if self.rendering:
            self.rendering.SetInput(self)
    
    def Hide(self):
        if self.rendering:
            self.rendering.Hide()
    def Show(self):
        if self.rendering:
            self.rendering.Show()
            
    def SetLabel(self,lbl,color=[1,1,1]):
        """ associates a color and a text label to the spline """

        self.rendering.profile.GetProperty().SetColor(color)
        self.label = lbl

def CommitSplineRenderingProperties():
    splines = Globals.objectSet.GetAll()
    for s in splines:
        s.rendering.UpdateProperties()

    Globals.renWin.Render()
        
class SplineSet:

    def __init__(self):
        self.dic = {}
        self.visibleSplines = []
        self.nextIndex = 0 # global counter for the splines added
        
        self.labels = { } # holds a dictionary with the available labels
        self.currentLabel = None
        # try to read the labels from a file when created
        self.ReadLabelsFile("labels.txt")

    def ReadLabelsFile(self,fileName):
        
        f = open(fileName, 'r')
        contents = f.read()
        f.close()

        lines = contents.split('\n')
        for line in lines:
            items = line.split()
            if len(items):
                color = []
                for c in items[1:]:
                    color.append( float(c))
                self.labels[items[0]] = color

    def SaveLabelsFile(self,fileName):
        f = open(fileName,'w')
        for k in self.labels:
            line = k + ' ' + str(self.labels[k]) + '\n'
            f.write(line)
        f.close()

    def DeleteLabel(self,lbl):
        if self.labels.has_key(lbl):
            del self.labels[lbl]
            splines = self.GetSplinesWithLabel(lbl)
            for s in splines:
                s.SetLabel(None)

    def DeleteVisibleSpline(self,splineId):
        spline = self.visibleSplines[splineId]
        if spline == None: return

        sliceIndex = spline.sliceIndex
        spline.Hide()
        self.dic[ sliceIndex ].remove(spline)
        if len(self.dic[sliceIndex])==0:
            del self.dic[sliceIndex]
        self.UpdateFromGUI()

    def SetLabel(self,lbl,color):
        if color==None:
            color = (1,1,1)
        else:
            buf = []
            for c in color:
                buf.append(float(c)/255.0)
            color = buf
            
        if self.labels.has_key(lbl):
            splines = self.GetSplinesWithLabel(lbl)
            for s in splines: s.SetLabel(color)
            Globals.renWin.Render()
        self.labels[lbl] = color

    def GetSplinesWithLabel(self,lbl):
        splines = []
        all = self.GetAll()
        for s in all:
            if s.label == lbl: splines.append(s)
        return splines

    def AddSpline(self, spline, lbl=None, index=None):
        """ adds a spline to the set
        the index is useful for loading a set """

        if index:
            spline.globalIndex = index
            self.nextIndex = max(index+1,self.nextIndex)
        else:
            spline.globalIndex = self.nextIndex
            self.nextIndex += 1

        if self.dic.has_key(spline.sliceIndex):
            self.dic[spline.sliceIndex].append(spline)
        else:
            self.dic[spline.sliceIndex] = [spline]

        self.visibleSplines.append(spline)
        if lbl:
            spline.SetLabel(lbl,
                            self.labels[lbl])
        elif self.currentLabel:
            spline.SetLabel(self.currentLabel,
                            self.labels[self.currentLabel])

    def UpdateFromGUI(self):
        """ only display splines matching the current position"""

        for s in self.visibleSplines: s.Hide()

        si  = Globals.imagePipeline.GetSliceIndex()
        ori = Globals.imagePipeline.GetOrientation()

        splines = self.GetSplines(si,ori)
        self.visibleSplines = splines
        
        for s in splines: s.Show()

    def GetSplines(self,sliceIndex,ori):
        splines = []
        if self.dic.has_key(sliceIndex):
            for s in self.dic[sliceIndex]:
                if s.orientation == ori:
                    splines.append(s)
        return splines

    def GetAll(self):
        splines = []
        for k in self.dic:
            lst = self.dic[k]
            splines.extend(lst)
        return splines
