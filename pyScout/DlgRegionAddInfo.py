
from Tkinter import *
from tkFont import *

class DlgRegionAddInfo(Toplevel):

    def __init__(self,parent,
                 nextCallback,
                 resetCallback,
                 delCallback):

        self.nextCallback = nextCallback
        self.resetCallback = resetCallback
        self.delCallback = delCallback

        self.sliceIndex = None
        self.orientation = None

        Toplevel.__init__(self,parent)
        self.title("Add Region Status")

        self.parent = parent
        self.result = None

        self.build()

        def ignore(): pass
        self.protocol("WM_DELETE_WINDOW", ignore)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50 )
                      )

    def SetRegionInformation(self,index,orientation):
        self.sliceIndex = index
        self.orientation = orientation

        self.lblSliceIndex.configure(text="Slice Index : " + str(self.sliceIndex))
        self.lblOrientation.configure(text="Orientation : " + str(self.orientation))

    def body(self,master):
        pass

    def build(self):
        # pack labels for the state
        frmState = Frame( self )
        frmState.pack(padx=3,pady=3,
                      side="top",
                      anchor="n",
                      fill="both",
                      expand="false"
                      )

        lblZero = Label(frmState,
                             text = "Add Zero Curve"
                             )
        lblOne = Label(frmState,
                            text = "Add One Curve"
                            )
        self.labels = ( lblZero,lblOne )
        self.stateIndex = 0
        for lbl in self.labels:
            lbl.pack(side="top")
        self.reset()
        # pack information
        frmInfo = Frame( self )
        frmInfo.pack(padx=3,pady=3,
                     side="top",
                     anchor="n",
                     fill="both",
                     expand="false"
                     )
        self.lblSliceIndex = Label(frmInfo,
                                   text = "Slice Index : " + str(self.sliceIndex)
                                   )
        self.lblOrientation = Label(frmInfo,
                                    text = "Orientation : " + str(self.orientation)
                                    )
        for lbl in (self.lblSliceIndex,self.lblOrientation):
            lbl.pack(side="left")

        # pack buttons
        frmBtns = Frame(self)
        frmBtns.pack(padx=3,pady=3,
                     side="top",
                     anchor="n",
                     fill="both",
                     expand="false")
        
        btnNext = Button(frmBtns,
                         text = "Next",
                         command=self.next)
        btnDel = Button(frmBtns,
                        text="Delete Last Point",
                        command=self.Del)
        btnReset = Button(frmBtns,
                          text = "Reset State",
                          command=self.reset)
        for btn in (btnNext, btnDel, btnReset):
            btn.pack(side="left")

    def reset(self):
        self.labels[self.stateIndex].config(relief="flat")
        self.stateIndex = 0
        self.labels[self.stateIndex].config(relief="raised")
        
        self.resetCallback()

    def next(self):
        """ will succeed if the underlying is valid, i.e. if points were added, for example """
        if self.nextCallback():
            self.labels[self.stateIndex].config(relief="flat")
            self.stateIndex = (self.stateIndex+1) % len(self.labels)
            self.labels[self.stateIndex].config(relief="raised")

    def Del(self):
        self.delCallback()

    def cancel(self, event=None):
        self.parent.focus_set()
        self.destroy()
