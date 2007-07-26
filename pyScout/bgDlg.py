
import tkSimpleDialog
from Tkinter import *

def __init__(self, parent, title = None):
    '''Initialize a dialog.

    Arguments:

    parent -- a parent window (the application window)

    title -- the dialog title
    '''
    Toplevel.__init__(self, parent)

    if parent.winfo_viewable():  # XXX this condition is the only "fix".
        self.transient(parent)

    if title:
        self.title(title)

    self.parent = parent

    self.result = None

    body = Frame(self)
    self.initial_focus = self.body(body)
    body.pack(padx=5, pady=5)

    self.buttonbox()

    self.wait_visibility() # window needs to be visible for the grab
    self.grab_set()

    if not self.initial_focus:
        self.initial_focus = self

    self.protocol("WM_DELETE_WINDOW", self.cancel)

    if self.parent is not None:
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

    self.initial_focus.focus_set()

    self.wait_window(self)

