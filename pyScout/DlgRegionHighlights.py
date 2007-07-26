
from Tkinter import *
from tkFont import *

class DlgRegionHighlights(Toplevel):

    def __init__(self,parent,
                 updateCallback,
                 gotoRegionCallback):

        self.updateCallback = updateCallback
        self.gotoRegionCallback = gotoRegionCallback

        Toplevel.__init__(self,parent)
        self.title("Highlight Region")

        self.parent = parent
        self.result = None

        self.build()

        def ignore(): pass
        self.protocol("WM_DELETE_WINDOW",ignore)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50 )
                      )

    def build(self):

        frmRegion = Frame( self )
        lblRegion = Label(frmRegion,
                          text = " Region " )
        self.region = StringVar()
        txtRegion = Label( frmRegion,
                           textvariable=self.region)
        for w in [lblRegion,txtRegion]: w.pack(side="left")

        frmProfiles = Frame( self )
        lblProfiles = Label(frmProfiles,
                            text = " Profiles ")
        self.profiles = StringVar()
        txtProfiles = Entry( frmProfiles,
                             textvariable=self.profiles)
        for w in [lblProfiles,txtProfiles]: w.pack(side="left")

        frmBtn = Frame( self )
        btnUpdate = Button( frmBtn,
                            text = "Update",
                            command=self.update )
        btnGoToRegion = Button( frmBtn,
                                text = "GoTo Region Index",
                                command=self.gotoRegion )
        for w in [btnUpdate]: w.pack(side="left")

        for f in [frmRegion,frmProfiles,frmBtn]:
            f.pack(padx=3,pady=3,
                   side="top",
                   anchor="n",
                   fill="both",
                   expand="false")

    def gotoRegion(self):

        # get the contents of the entry control
        val = int(self.region.get())
        self.gotoRegionCallback(val)

    def SetRegion(self,regionId):
        self.region.set(regionId)

    def update(self):

        # get the contents of the entry control and split it
        vals = self.profiles.get().split()
        self.updateCallback(vals)

    def SelectRegion(self,idx):

        self.region.set(idx)

    def cancel(self, event=None):
        self.parent.focus_set()
        self.destroy()

    def ClearData(self):
        self.region.set("")
        self.profiles.set("")
