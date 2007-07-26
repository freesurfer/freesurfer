#! /usr/bin/python

# Gheorghe Postelnicu, 2006
#
#


import Tkinter
import bgDlg
import tkSimpleDialog
import tkFileDialog
import tkColorChooser

import vtk
from vtk.tk.vtkTkRenderWindowInteractor import \
     vtkTkRenderWindowInteractor

import sys
import tkMessageBox

from ImagePipeline import ImagePipeline
import Mode
import Globals

import CrossHair
import Region
import Laplace
import RenderingProperties
from StateManager import LoadState, SaveState, SavePreferences, LoadPreferences
from ProfileExport import ProfileExport
from DlgImageContrast import DlgImageContrast
from DlgPreferences import DlgPreferences

class Application:

    def __init__(self, root,
                 w=600,h=600):
        self.root = root

        # VTK setup
        Globals.ren = vtk.vtkRenderer()
        Globals.renWin = vtk.vtkRenderWindow()
        Globals.renWin.AddRenderer(Globals.ren)

        # setup main actor
        Globals.imagePipeline = ImagePipeline(Globals.ren)
        def SetReferenceSpacing():
            spacing = Globals.imagePipeline.GetSpacing()
            if spacing[2]>0: Globals.referenceSize = min( spacing )
            else: Globals.referenceSize = min(spacing[0:2])
        SetReferenceSpacing()
        Globals.ren.AddActor(Globals.imagePipeline.actor)

        Globals.SetDefaultRenderProps()
        Globals.SetDefaultPreferences()
        Globals.crossHair = CrossHair.CrossHairWidget()
        # setup observer chain for main events
        Globals.objectSet = Region.RegionSet()
        self.imageObservers = [ Globals.objectSet.UpdateFromGUI ]
        
        self.InitializeGUI(w,h)

        # initialize the interactor
        iren = self.renderWidget.GetRenderWindow().GetInteractor()
        Globals.ren.SetBackground(.5,.5,.5)

        Globals.ren.ResetCamera()

        iren.Initialize()
        self.mode = Mode.Navigation()
        self.mode.Start() # set the interactor style
        Globals.renWin.Render()
        iren.Start()

        # set the initial position of the pointer in the center
        bounds = Globals.imagePipeline.slicer.GetOutput().GetBounds()
        centerX = .5 * ( bounds[0] + bounds[1] )
        centerY = .5 * ( bounds[2] + bounds[3] )
        Globals.crossHair.SetPosition(centerX, centerY)


    def InitializeGUI(self,w,h):

        self.top = Tkinter.Toplevel(self.root)

        # make sure exit happens gracefully
        def quit(obj=self.root): obj.quit()
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
        self.root.update()

        frmModes = Tkinter.Frame(self.top)
        frmModes.pack()

        self.radioVar = Tkinter.IntVar()
        self.radioVar.set(0)

        self.dicGuiModes = { 0:"Navigation",
                             1:"Add Region",
                             2:"Laplace",
                             3:"Profile Inspector",
                             4:"Highlight",
                             5:"Delete Region"
                             }

        for k in self.dicGuiModes.keys():
            Tkinter.Radiobutton(frmModes,
                                text=self.dicGuiModes[k],
                                variable=self.radioVar,
                                value=k,
                                command=self.ChangeMode).pack(side="left")

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
        self.scaleSlice.pack(fill="x",expand="false")

        self.SetupMenubar()

    def SetupMenubar(self):

        self.menu = Tkinter.Menu(self.top)

        # File menu
        mnFile = Tkinter.Menu(self.menu,tearoff=0)
        mnFile.add_command(label="Save",
                           command=SaveState)
        def ProxySlice(si):
            tksi = Tkinter.IntVar()
            tksi.set( si )
            self.scaleSlice.config(variable=tksi)
            self.SetSlice(si)
        def Load():
            LoadState(ProxySlice,self.SetOrientation)
        mnFile.add_command(label="Load",
                           command=Load)

        mnFile.add_separator()

        mnFile.add_command(label="Save Prefs",
                           command=SavePreferences)
        mnFile.add_command(label="Load Prefs",
                           command=LoadPreferences)
        mnFile.add_separator()

        mnFile.add_command(label="Export Intensities",
                           command = ProfileExport)
        
        def SaveCapture():
            fname = tkFileDialog.asksaveasfilename(defaultextension=".tif",
                                                   filetypes=[ ('TIF Files','*.tif') ] )
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
        mnFile.add_command(label="Exit",
                           command=self.top.quit)
        self.menu.add_cascade(label="File",
                              menu=mnFile)

        mnEdit = Tkinter.Menu(self.menu, tearoff=0)
        def ShowPreferences(): DlgPreferences(Globals.root,
                                              "Preferences")
        mnEdit.add_command(label="Preferences",
                           command=ShowPreferences)
        def ChooseHighlightsColor():
            retVal = tkColorChooser.askcolor(parent=Globals.root,
                                             title="Highlights Color",
                                             initialcolor=( \
                int(Globals.preferences["highlightRed"]*255.0),
                int(Globals.preferences["highlightGreen"]*255.0),
                int(Globals.preferences["highlightBlue"]*255.0)
                ) )
            if retVal:
                Globals.preferences["highlightRed"] = float(retVal[0][0]) / 255.0
                Globals.preferences["highlightGreen"] = float(retVal[0][1]) / 255.0
                Globals.preferences["highlightBlue"] = float(retVal[0][2]) / 255.0

        mnEdit.add_command(label="Highlight Colors",
                           command=ChooseHighlightsColor)
                    
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
        def ShowRenDialog():
            dlg = RenderingProperties.RenderingProperties(Globals.root,
                                                          Region.CommitSplineRenderingProperties)
        mnView.add_command(label="Rendering Properties",
                           command=ShowRenDialog)
        def ShowLevelsDlg(): dlg=DlgImageContrast(self.root,
                                                  "Image Contrast")
        mnView.add_command(label="Image Contrast",
                           command=ShowLevelsDlg)
        self.menu.add_cascade(label="View",
                              menu=mnView)

        # Tools menu
        mnTools = Tkinter.Menu(self.menu, tearoff=0)
        def CompleteLaplace():
            regions = Globals.objectSet.GetAllRegions()
            for r in regions: r.ApplyLaplace()
            Globals.objectSet.UpdateAllProfiles()
        mnTools.add_command(label="Complete Laplace",
                            command=CompleteLaplace)
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
        smin,smax = Globals.imagePipeline.GetAxisExtent()
        self.NotifyObservers()
        Globals.renWin.Render()

    def SetOrientation(self,ori):
        refWP = Globals.imagePipeline.GetWorldPoint( Globals.crossHair.x,
                                                     Globals.crossHair.y )
        Globals.imagePipeline.SetOrientation(ori)
        (modWP,si) = Globals.imagePipeline. \
                     WorldPointToDisplayFrame(refWP)
        Globals.crossHair.SetPosition( modWP[0], modWP[1] )

        smin,smax = Globals.imagePipeline.GetAxisExtent()
        si = max( smin, min(smax, si))
        self.SetSlice(si)
        sliceNo = Tkinter.IntVar()
        sliceNo.set( si )
        self.scaleSlice.config(variable=sliceNo, \
                               from_ = smin,
                               to = smax)
        self.NotifyObservers()
        Globals.renWin.Render()

    def NotifyObservers(self):
        """ in charge of notifying the observers
        this is done in another thread """

        def update():
            for o in self.imageObservers:
                o()

        # multi-threading seems to cause some update problems
        update()

    def ChangeMode(self):
        # destructor for previous mode
        self.mode.Finish()
        
        mode = self.radioVar.get()
        if mode==0:
            self.mode = Mode.Navigation()
        elif mode==1:
            self.mode = Mode.AddRegion(self)
        elif mode==2:
            self.mode = Mode.Laplace(self)
        elif mode==3:
            self.mode = Mode.ProfileInspector(self)
        elif mode==4:
            self.mode = Mode.Highlight(self)
        else:
            self.mode = Mode.DeleteRegion(self)

        self.mode.Start()



####################################

Laplace.InitializePackage()
tkSimpleDialog.Dialog.__init__ = bgDlg.__init__

Globals.root = Tkinter.Tk()
Globals.root.withdraw()

app = Application(Globals.root,700,700)
Globals.root.mainloop()

Laplace.Finalize()

        
        
