
from ImagePipeline import ImagePipeline

renderProps = None

imagePipeline = None
crossHair = None
referenceSize = 1.0

objectSet = None

preferences = None

ren = None
renWin = None

root = None

def SetDefaultRenderProps():
    global renderProps

    renderProps = { "sphereSize" : 3,
                    "spherePhiResolution" : 10,
                    "sphereThetaResolution" : 10,
                    "profileNumberOfSides" : 10,
                    "tubeSize" : 1,
                    "opacity" : 1
                    }

def SetDefaultPreferences():
    global preferences

    preferences = { "laplaceResolution" : .5,
                    "laplaceConvergence" : 8,
                    "profileSpacing" : 2,
                    "profilePoints" : 100,
                    "highlightRed" : .5,
                    "highlightGreen" : .5,
                    "highlightBlue" : .5
                    }

def SetDefaultSplinePreferences():
    global preferences

    preferences = { "numberOfPoints" : 100 }

