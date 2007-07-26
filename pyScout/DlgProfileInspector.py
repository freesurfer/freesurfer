
from Tkinter import *

class DlgProfileInspector(Toplevel):

    def __init__(self,parent,
                 profileCallback,
                 pointCallback):

        # setup callbacks
        self.profileCallback = profileCallback
        self.pointCallback = pointCallback

        # setup state information
        self.regionIndex = None

        # setup widget 
        Toplevel.__init__(self,parent)
        self.title("Profile Inspector")

        self.parent = parent
        self.result = None

        self.build()

        def ignore(): pass
        self.protocol("WM_DELETE_WINDOW", ignore)

    def cancel(self, event=None):
        self.parent.focus_set()
        self.destroy()

    def body(self,master): pass

    def build(self):

        self.strRegion = StringVar()
        self.strRegion.set( " Region : " + str(self.regionIndex))
        Label(self,
              textvariable = self.strRegion).pack(
            padx=3,pady=3,
            side="top",
            anchor="n",
            fill="both",
            expand="false")

        self.scaleProfile = Scale(self,
                                  from_=0, to=0,
                                  orient="horizontal",
                                  command=self.SetProfile,
                                  label="Profile Index",
                                  state="disabled")
        self.scaleProfile.pack(
            padx=3,pady=3,
            side="top",
            anchor="n",
            fill="both",
            expand="false")

        self.scalePoint = Scale(self,
                                from_=0,to=0,
                                orient="horizontal",
                                command=self.SetPoint,
                                label="Point Index",
                                state="disabled")
        self.scalePoint.pack(
            padx=3,pady=3,
            side="top",
            anchor="n",
            fill="both",
            expand="false")

    def SetProfile(self,value): self.profileCallback( int(value) )
    def SetPoint(self,value): self.pointCallback( int(value) )


    

    def SetSelectedRegion(self,region=None,profile=None,point=None):

        if point==None:
            self.scaleProfile.config(state="disabled")
            self.scalePoint.config(state="disabled")
            self.strRegion.set(" Region : None")
            return

        pointState = self.scalePoint.cget("state")
        if pointState == "disabled":
            self.regionIndex = region.index
            self.strRegion.set(" Region : " + str(region.index))
            print " no profiles = ",  len(region.profiles)
            self.scaleProfile.config(state="active",
                                     from_=0,
                                     to=len(region.profiles)-1
                                     )
            
            self.scalePoint.config(state="active",
                                     from_=0,
                                     to=len(region.profiles[profile].points)-1
                                     )
            
        if self.regionIndex != region.index:
            self.strRegion.set(" Region : " + str(region.index))
            self.scaleProfile.config(from_=0,
                                     to=len(region.profiles)-1)
            self.scalePoint.config(from_=0,
                                   to=len(region.profiles[profile].points)-1)   

        # specify the values of the controls
        dprofile = DoubleVar()
        dprofile.set(profile)
        self.scaleProfile.config(variable=dprofile)
        
        dpoint = DoubleVar()
        dpoint.set(point)
        self.scalePoint.config(variable=dpoint)

            

    def ActivateScales(self):
        state = self.scaleProfile.cget("state")
        if state == "disabled":
            self.scaleProfile.config(state="active")
            self.scalePoint.config(state="active")
    def DeactivateScales(self):
        state = self.scaleProfile.cget("state")
        if state == "enabled":
            self.scaleProfile.config(state="disabled")
            self.scalePoint.config(state="disabled")
        
