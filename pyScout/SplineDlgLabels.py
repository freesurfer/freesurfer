
from Tkinter import *
import tkColorChooser
import tkSimpleDialog

import Globals

from SplineSet import GetTkColor

class SplineDlgLabels(tkSimpleDialog.Dialog):

    def body(self,master):
        """ define the controls of the widget """

        frm = Frame(master)
        frm.pack(padx = 3, pady=3)

        self.text = StringVar()
        self.text.set("default")

        self.btn = Button(frm,
                          text="Choose Color",
                          command=self.SelectColor)
        self.color = None

    def SelectColor(self):
        colorTuple = tkColorChooser.askcolor(color=self.color)
        self.color = colorTuple[0]
        if self.color:
            self.btn.config(bg=GetTkColor(self.color))

    def ok(self,event=None):

        Globals.objectSet.SetLabel(self.text.get(),
                                   self.color)
        tkSimpleDialog.Dialog.cancel(self,event)

    
        
