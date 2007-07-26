
from Tkinter import *
import tkSimpleDialog

import Globals


class DlgPreferences(tkSimpleDialog.Dialog):

    def body(self, master):
        """ define the controls of the widget """

        frmLaplace = self.CreateFrame(master,
                    "Laplace Solver Options")

        self.laplaceResolution = StringVar()
        self.laplaceResolution.set( Globals.preferences["laplaceResolution"])
        self.CreateLabelAndEntry(frmLaplace,
                                 "Resolution",
                                 self.laplaceResolution)
        
        self.laplaceConvergence = StringVar()
        self.laplaceConvergence.set( \
            Globals.preferences["laplaceConvergence"] )
        self.CreateLabelAndEntry(frmLaplace,
                                 "Convergence",
                                 self.laplaceConvergence)
        
        frmProfiles = self.CreateFrame(master,
                 " Profile Options ")

        self.profileSpacing = StringVar()
        self.profileSpacing.set(
            Globals.preferences["profileSpacing"] )
        self.CreateLabelAndEntry(frmProfiles,
                                 "Spacing",
                                 self.profileSpacing)
        
        self.profilePoints = StringVar()
        self.profilePoints.set(
            Globals.preferences["profilePoints"] )
        self.CreateLabelAndEntry(frmProfiles,
                                 "Points",
                                 self.profilePoints)
        

    def CreateFrame(self,master,title):
        """ create and pack a new frame """

        frm = Frame(master)
        frm.pack( padx=3,pady=3,
                  side="top",
                  anchor="n",
                  fill="both",
                  expand="false"
                  )
        lbl = Label(frm,
                    text=title
                    )
        lbl.pack(side="top")
        return frm

    def CreateLabelAndEntry(self,
                            frmTop,
                            lblText,
                            var):
        frm = Frame(frmTop)
        frm.pack(side="top")
        
        lbl = Label(frm,
                    text=lblText)
        entry = Entry(
            frm,
            textvariable=var)

        for w in (lbl,entry):
            w.pack(side="left")


    def ok(self,event=None):

        needsLaplace = 0
        resolution = float(
            self.laplaceResolution.get() )
        convergence = float(
            self.laplaceConvergence.get() )
        if resolution != Globals.preferences["laplaceResolution"]:
            needsLaplace = 1
            Globals.preferences["laplaceResolution"] = resolution
        if convergence != Globals.preferences["laplaceConvergence"]:
            needsLaplace = 1
            Globals.preferences["laplaceConvergence"] = convergence

        needsProfiles = 0
        spacing = float(
            self.profileSpacing.get() )
        points = float(
            self.profilePoints.get() )
        if spacing != Globals.preferences["profileSpacing"]:
            needsProfiles = 1
            Globals.preferences["profileSpacing"] = spacing
        if points != Globals.preferences["profilePoints"]:
            needsProfiles = 1
            Globals.preferences["profilePoints"] = points

        # update all the necessary stuff
        if needsLaplace:
            for r in Globals.objectSet.GetAllRegions():
                r.ApplyLaplace()

        if needsLaplace or needsProfiles:
            Globals.objectSet.UpdateAllProfiles()
            Globals.renWin.Render()
            
        tkSimpleDialog.Dialog.cancel(self,event)
