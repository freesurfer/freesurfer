#! /usr/bin/python

# Gheorghe Postelnicu, 2006
#
#

import Tkinter
import bgDlg
import tkSimpleDialog
import tkFileDialog

import vtk
from vtk.tk.vtkTkRenderWindowInteractor import \
     vtkTkRenderWindowInteractor

import tkMessageBox

from ImagePipeline import ImagePipeline
import CrossHair
from RenderingProperties import RenderingProperties
from DlgImageContrast import DlgImageContrast

import SplineMode
from SplineSet import *
from SplineExport import SplineExport
from SplineStateManager import LoadState, SaveState
from SplineDlgPreferences import SplineDlgPreferences

import Globals

class Application:

    def __init__(self,root,
                 w=600,h=600):

        # VTK setup
        Globals.ren = vtk.vtkRenderer()
        Globals.renWin = vtk.vtkRenderWindow()
        Globals.renWin.AddRenderer( Globals.ren )

        # setup the main actor
        Globals.imagePipeline = ImagePipeline( Globals.ren )
        def SetReferenceSpacing():
            spacing = Globals.imagePipeline.GetSpacing()
            if spacing[2]>0: Globals.referenceSize = min( spacing )
            else: Globals.referenceSize = min(spacing[0:2])
        SetReferenceSpacing()
        Globals.ren.AddActor(Globals.imagePipeline.actor)

        Globals.SetDefaultRenderProps()
        Globals.SetDefaultSplinePreferences()
        Globals.crossHair = CrossHair.CrossHairWidget()
        Globals.objectSet = SplineSet()
        self.imageObservers = [ Globals.objectSet.UpdateFromGUI ]

        self.InitializeGUI(w,h)

        # initialize the interactor
        iren = self.renderWidget.GetRenderWindow().GetInteractor()
        Globals.ren.SetBackground(.5,.5,.5)

        Globals.ren.ResetCamera()

        iren.Initialize()
        self.mode = SplineMode.Navigation()
        self.mode.Start()
        Globals.renWin.Render()
        iren.Start()

        # set the initial position of the pointer in the center
        bounds = Globals.imagePipeline.slicer.GetOutput().GetBounds()
        centerX = .5 * ( bounds[0] + bounds[1] )
        centerY = .5 * ( bounds[2] + bounds[3] )
        Globals.crossHair.SetPosition(centerX, centerY)

    def InitializeGUI(self,w,h):

        self.top = Tkinter.Toplevel(Globals.root)

        # handle graceful exit
        def quit(obj=Globals.root): obj.quit()
        self.top.protocol("WM_DELETE_WINDOW",quit)

        self.frmDisplay = Tkinter.Frame(self.top)
        self.frmRender = Tkinter.Frame(self.frmDisplay)

        for f in (self.frmDisplay,self.frmRender):
            f.pack(padx=3,pady=3,
                   side="top",
                   anchor="n",
                   fill="both",
                   expand="false")

        rwi = vtkTkRenderWindowInteractor(self.frmRender,
                                          rw=Globals.renWin,
                                          width=w,height=h)
        rwi.pack()
        rwi.Render()
        self.renderWidget = rwi
        Globals.root.update()

        #############################
        #
        # setup modes
        frmModes = Tkinter.Frame(self.top)
        frmModes.pack()

        self.dicGuiModes = { 0:"Navigation",
                             1:"Add Spline"
                             }

        self.radioVar = Tkinter.IntVar()
        self.radioVar.set(0)
        for k in self.dicGuiModes.keys():
            Tkinter.Radiobutton(frmModes,
                                text=self.dicGuiModes[k],
                                variable=self.radioVar,
                                value=k,
                                command=self.ChangeMode).pack(
                side="left")
            
        #############################
        #
        # setup labels
        frmLabels = Tkinter.Frame(self.top)
        frmLabels.pack()

        self.dicLabels = {}
        counter = 0
        for lbl in Globals.objectSet.labels:
            self.dicLabels[counter] = lbl
            counter += 1
        self.lblVar = Tkinter.IntVar()
        self.lblVar.set(0)
        Globals.objectSet.currentLabel = self.dicLabels[0]
        for k in self.dicLabels.keys():
            Tkinter.Radiobutton(frmLabels,
                                text=self.dicLabels[k],
                                variable=self.lblVar,
                                value=k,
                                command=self.ChangeCurrentLabel).pack(
                side="left")

        ####################################
        # Add a slice scale to browse the current slice stack
        sliceNumber = Tkinter.IntVar()
        smin,smax = Globals.imagePipeline.GetAxisExtent()
        sliceNumber.set( Globals.imagePipeline.GetSliceIndex() )

        self.scaleSlice = Tkinter.Scale(self.top,
                                        from_=smin, to=smax,
                                        orient="horizontal",
                                        command=self.SetSlice,
                                        variable=sliceNumber,
                                        label="Slice")
        self.scaleSlice.pack(fill="x", expand="false")

        self.SetupMenubar()

    def SetupMenubar(self):

        self.menu = Tkinter.Menu(self.top)

        # File menu
        mnFile = Tkinter.Menu(self.menu, tearoff=0)
        # todo - SaveState
        #mnFile.add_command(label="Save",
        #                   command=SaveState)

        def ProxySlice(si):
            tksi = Tkinter.IntVar()
            tksi.set( si )
            self.scaleSlice.config( variable = tksi )
            self.SetSlice( si )
        def Load():
            LoadState(ProxySlice, self.SetOrientation)
        mnFile.add_command(label="Load",
                           command=Load)
        mnFile.add_command(label="Save State",
                           command=SaveState)

        def SaveCapture():
            fname = tkFileDialog.asksaveasfilename(defaultextension=".tif",
                                                   filetypes=[ ('TIF Files', '*.tif') ]
                                                   )

            if fname:
                w2i = vtk.vtkWindowToImageFilter()
                writer = vtk.vtkTIFFWriter()
                w2i.SetInput(Globals.renWin)
                writer.SetInputConnection( w2i.GetOutputPort() )
                writer.SetFileName(fname)
                Globals.renWin.Render()
                writer.Write()

        mnFile.add_command(label="Save TIF",
                           command=SaveCapture)
        mnFile.add_command(label="Export Intensities",
                           command=SplineExport)
        mnFile.add_command(label="Exit",
                           command=self.top.quit)
        self.menu.add_cascade(label="File",
                              menu=mnFile)

        # Edit menu
        mnEdit = Tkinter.Menu(self.menu,tearoff=0)
        def ShowPreferences(): SplineDlgPreferences(Globals.root,
                                                    "Preferences")
        mnEdit.add_command(label="Preferences",
                           command=ShowPreferences)
        self.menu.add_cascade(label="Edit",
                              menu=mnEdit)
        
        # View menu
        mnView = Tkinter.Menu(self.menu, tearoff=0)

        # orientation
        mnOrientation = Tkinter.Menu(mnView)
        def SetOrientationXY(): self.SetOrientation("xy")
        mnOrientation.add_command(label="XY",
                                  command=SetOrientationXY)
        def SetOrientationXZ(): self.SetOrientation("xz")
        mnOrientation.add_command(label="XZ",
                                  command=SetOrientationXZ)
        def SetOrientationYZ(): self.SetOrientation("yz")
        mnOrientation.add_command(label="YZ",
                                  command=SetOrientationYZ)
        mnView.add_cascade(label="Orientation",
                           menu=mnOrientation)
        mnView.add_separator()
        def ShowRenDialog(): RenderingProperties(Globals.root,
                                                 CommitSplineRenderingProperties)
        mnView.add_command(label="Rendering Properties",
                           command=ShowRenDialog)
        def ShowLevelsDlg(): DlgImageContrast(Globals.root,
                                              "Image Contrast")
        mnView.add_command(label="Image Contrast",
                           command=ShowLevelsDlg)
        
        self.menu.add_cascade(label="View",
                              menu=mnView)

        # Tools menu
        mnTools = Tkinter.Menu(self.menu, tearoff=0)
        def DeleteSpline():
            # 1. graft it as a new mode
            self.radioVar.set(0)
            tmpMode = SplineMode.CmdDeleteSpline(self)
            tmpMode.Start()
        mnTools.add_command(label="Delete Spline",
                            command=DeleteSpline)
        def GoToSlice():
            sliceIndex = tkSimpleDialog.askinteger("GoTo Slice"," Index : ")
            if not sliceIndex: return
            smin,smax = Globals.imagePipeline.GetAxisExtent()
            sliceIndex = min(smax, sliceIndex)
            sliceIndex = max(smin, sliceIndex)
            self.SetSlice(sliceIndex)
            self.scaleSlice.set(sliceIndex)
        mnTools.add_command(label="GoTo Slice",
                            command=GoToSlice)
        self.menu.add_cascade(label="Tools",
                              menu=mnTools)
        
        # Help menu
        mnHelp = Tkinter.Menu(self.menu,tearoff=0)

        def About():
            tkMessageBox.showinfo("About",
                                  "pyScout\n Gheorghe Postelnicu, 2006"
                                  )

        mnHelp.add_command(label="About",command=About)
        self.menu.add_cascade(label="Help",
                              menu=mnHelp)

        # set application menu
        self.top.config(menu=self.menu)

        
    def SetSlice(self,si):
        Globals.imagePipeline.SetSlice(int(si))
        self.NotifyObservers()
        Globals.renWin.Render()

    def SetOrientation(self,ori):
        refWP = Globals.imagePipeline.GetWorldPoint( Globals.crossHair.x,
                                                Globals.crossHair.y )
        Globals.imagePipeline.SetOrientation(ori)
        (modWP,si) = Globals.imagePipeline. \
                     WorldPointToDisplayFrame(refWP)
        Globals.crossHair.SetPosition( modWP[0], modWP[1] )

        self.SetSlice(si)
        sliceNo = Tkinter.IntVar()
        sliceNo.set( si )
        smin,smax = Globals.imagePipeline.GetAxisExtent()
        self.scaleSlice.config( variable=sliceNo,
                                from_=smin, to=smax )
        self.NotifyObservers()

        Globals.renWin.Render()

    def NotifyObservers(self):
        """ in charge of notifying the observers
        this is done in another thread """

        def update():
            for o in self.imageObservers: o()
        update() # multi-threading seems to cause update problems

    def ChangeMode(self):

        if self.mode:
            self.mode.Finish()

        mode = self.radioVar.get()
        if   mode == 0:  self.mode = SplineMode.Navigation()
        elif mode == 1:  self.mode = SplineMode.AddSpline(self)

        self.mode.Start()
    def ChangeCurrentLabel(self):

        lbl = int(self.lblVar.get())
        Globals.objectSet.currentLabel = self.dicLabels[lbl]


#################################################################

tkSimpleDialog.Dialog.__init__ = bgDlg.__init__

Globals.root = Tkinter.Tk()
Globals.root.withdraw()


app = Application(Globals.root,700,700)
Globals.root.mainloop()

