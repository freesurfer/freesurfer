
from Tkinter import *

import Globals

smartCounter =0

class RenderingProperties(Toplevel):

    def __init__(self,parent,commitCallback):
        global smartCounter
        if smartCounter:return

        self.commitCallback = commitCallback
        smartCounter+=1
        
        Toplevel.__init__(self,parent)
        self.title("Add Region Status")

        self.parent = parent
        self.result = None

        self.build()

        self.protocol("WM_DELETE_WINDOW", self.close)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50 )
                      )

    def close(self):
        global smartCounter
        
        self.parent.focus_set()
        self.destroy()
        smartCounter = 0

    def body(self,master):
        pass

    def build(self):
        """ define the controls of the widget """

        dOpacity = DoubleVar()
        dOpacity.set( Globals.renderProps["opacity"] )
        Scale(self,
              from_=0, to=1,
              resolution=.05,
              variable=dOpacity,
              orient="horizontal",
              command=self.SetOpacity,
              label="Opacity").pack(side="top",
                                    fill="x",
                                    expand="false")

        dSphereSize = DoubleVar()
        dSphereSize.set( Globals.renderProps["sphereSize"] )
        Scale(self,
              from_=.1, to=5.0,
              resolution=.1,
              variable=dSphereSize,
              orient="horizontal",
              command=self.SetSphereSize,
              label="Sphere Size").pack(side="top",
                                        fill="x",
                                        expand="false")

        dTubeSize = DoubleVar()
        dTubeSize.set( Globals.renderProps["tubeSize"])
        Scale(self,
              from_=.1, to=5.0,
              resolution=.1,
              variable=dTubeSize,
              orient="horizontal",
              command=self.SetTubeSize,
              label="Tube Size").pack(side="top",
                                      fill="x",
                                      expand="false")

        frmButtons = Frame(self)
        frmButtons.pack(padx=3,pady=3,
                        side="top",
                        anchor="n",
                        fill="both",
                        expand="false"
                        )
        def InitProps():
            Globals.renderProps = { "sphereSize" : 3,
                                    "spherePhiResolution" : 10,
                                    "sphereThetaResolution" : 10,
                                    "profileNumberOfSides" : 10,
                                    "tubeSize" : 1,
                                    "opacity" : 1
                                    }
            self.commitCallback()
        btnReset = Button( frmButtons,
                           text="Reset",
                           command=InitProps)
        btnReset.pack()

    def SetOpacity(self,val):
        Globals.renderProps["opacity"] = float(val)
        self.commitCallback()
    def SetSphereSize(self,val):
        Globals.renderProps["sphereSize"] = float(val)
        self.commitCallback()
    def SetTubeSize(self,val):
        Globals.renderProps["tubeSize"] = float(val)
        self.commitCallback()
