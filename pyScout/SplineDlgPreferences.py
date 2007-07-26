
from Tkinter import *
import tkSimpleDialog

import Globals

class SplineDlgPreferences(tkSimpleDialog.Dialog):

    def body(self,master):

        frm = Frame(master)
        frm.pack()

        lbl = Label(frm,
                    text="Number of Points")
        self.var = StringVar()
        self.var.set( Globals.preferences["numberOfPoints"] )
        entry = Entry(frm,
                      textvariable=self.var)

        for w in (lbl,entry):
            w.pack(side="left")

    def ok(self,event=None):
        Globals.preferences["numberOfPoints"] = int(self.var.get())
        tkSimpleDialog.Dialog.cancel(self,event)

