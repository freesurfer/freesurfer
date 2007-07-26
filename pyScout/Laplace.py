# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Laplace

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    def __repr__(self):
        return "<C std::vector<(double)> instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, DoubleVector, 'this', _Laplace.new_DoubleVector(*args))
        _swig_setattr(self, DoubleVector, 'thisown', 1)
    def __len__(*args): return _Laplace.DoubleVector___len__(*args)
    def __nonzero__(*args): return _Laplace.DoubleVector___nonzero__(*args)
    def clear(*args): return _Laplace.DoubleVector_clear(*args)
    def append(*args): return _Laplace.DoubleVector_append(*args)
    def pop(*args): return _Laplace.DoubleVector_pop(*args)
    def __getitem__(*args): return _Laplace.DoubleVector___getitem__(*args)
    def __getslice__(*args): return _Laplace.DoubleVector___getslice__(*args)
    def __setitem__(*args): return _Laplace.DoubleVector___setitem__(*args)
    def __setslice__(*args): return _Laplace.DoubleVector___setslice__(*args)
    def __delitem__(*args): return _Laplace.DoubleVector___delitem__(*args)
    def __delslice__(*args): return _Laplace.DoubleVector___delslice__(*args)
    def __del__(self, destroy=_Laplace.delete_DoubleVector):
        try:
            if self.thisown: destroy(self)
        except: pass

class DoubleVectorPtr(DoubleVector):
    def __init__(self, this):
        _swig_setattr(self, DoubleVector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DoubleVector, 'thisown', 0)
        _swig_setattr(self, DoubleVector,self.__class__,DoubleVector)
_Laplace.DoubleVector_swigregister(DoubleVectorPtr)

class WrapperLine(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, WrapperLine, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, WrapperLine, name)
    def __repr__(self):
        return "<C std::vector<(std::vector<(double)>)> instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, WrapperLine, 'this', _Laplace.new_WrapperLine(*args))
        _swig_setattr(self, WrapperLine, 'thisown', 1)
    def __len__(*args): return _Laplace.WrapperLine___len__(*args)
    def clear(*args): return _Laplace.WrapperLine_clear(*args)
    def append(*args): return _Laplace.WrapperLine_append(*args)
    def __nonzero__(*args): return _Laplace.WrapperLine___nonzero__(*args)
    def pop(*args): return _Laplace.WrapperLine_pop(*args)
    def __getitem__(*args): return _Laplace.WrapperLine___getitem__(*args)
    def __getslice__(*args): return _Laplace.WrapperLine___getslice__(*args)
    def __setitem__(*args): return _Laplace.WrapperLine___setitem__(*args)
    def __setslice__(*args): return _Laplace.WrapperLine___setslice__(*args)
    def __delitem__(*args): return _Laplace.WrapperLine___delitem__(*args)
    def __delslice__(*args): return _Laplace.WrapperLine___delslice__(*args)
    def __del__(self, destroy=_Laplace.delete_WrapperLine):
        try:
            if self.thisown: destroy(self)
        except: pass

class WrapperLinePtr(WrapperLine):
    def __init__(self, this):
        _swig_setattr(self, WrapperLine, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, WrapperLine, 'thisown', 0)
        _swig_setattr(self, WrapperLine,self.__class__,WrapperLine)
_Laplace.WrapperLine_swigregister(WrapperLinePtr)

class WrapperLineSet(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, WrapperLineSet, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, WrapperLineSet, name)
    def __repr__(self):
        return "<C std::vector<(std::vector<(std::vector<(double)>)>)> instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, WrapperLineSet, 'this', _Laplace.new_WrapperLineSet(*args))
        _swig_setattr(self, WrapperLineSet, 'thisown', 1)
    def __len__(*args): return _Laplace.WrapperLineSet___len__(*args)
    def clear(*args): return _Laplace.WrapperLineSet_clear(*args)
    def append(*args): return _Laplace.WrapperLineSet_append(*args)
    def __nonzero__(*args): return _Laplace.WrapperLineSet___nonzero__(*args)
    def pop(*args): return _Laplace.WrapperLineSet_pop(*args)
    def __getitem__(*args): return _Laplace.WrapperLineSet___getitem__(*args)
    def __getslice__(*args): return _Laplace.WrapperLineSet___getslice__(*args)
    def __setitem__(*args): return _Laplace.WrapperLineSet___setitem__(*args)
    def __setslice__(*args): return _Laplace.WrapperLineSet___setslice__(*args)
    def __delitem__(*args): return _Laplace.WrapperLineSet___delitem__(*args)
    def __delslice__(*args): return _Laplace.WrapperLineSet___delslice__(*args)
    def __del__(self, destroy=_Laplace.delete_WrapperLineSet):
        try:
            if self.thisown: destroy(self)
        except: pass

class WrapperLineSetPtr(WrapperLineSet):
    def __init__(self, this):
        _swig_setattr(self, WrapperLineSet, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, WrapperLineSet, 'thisown', 0)
        _swig_setattr(self, WrapperLineSet,self.__class__,WrapperLineSet)
_Laplace.WrapperLineSet_swigregister(WrapperLineSetPtr)


ComputeProfiles = _Laplace.ComputeProfiles

ComputeIsolines = _Laplace.ComputeIsolines

ComputeSolution = _Laplace.ComputeSolution

InitializePackage = _Laplace.InitializePackage

Finalize = _Laplace.Finalize
class WrapperPolyData(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, WrapperPolyData, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, WrapperPolyData, name)
    def __repr__(self):
        return "<C WrapperPolyData instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, WrapperPolyData, 'this', _Laplace.new_WrapperPolyData(*args))
        _swig_setattr(self, WrapperPolyData, 'thisown', 1)
    def __del__(self, destroy=_Laplace.delete_WrapperPolyData):
        try:
            if self.thisown: destroy(self)
        except: pass
    def InsertPoint(*args): return _Laplace.WrapperPolyData_InsertPoint(*args)
    def InsertNextCell(*args): return _Laplace.WrapperPolyData_InsertNextCell(*args)
    def InsertCellPoint(*args): return _Laplace.WrapperPolyData_InsertCellPoint(*args)
    def GetNumberOfCells(*args): return _Laplace.WrapperPolyData_GetNumberOfCells(*args)
    def InitTraversal(*args): return _Laplace.WrapperPolyData_InitTraversal(*args)
    def GoToNextCell(*args): return _Laplace.WrapperPolyData_GoToNextCell(*args)
    def GetCellNumberOfPoints(*args): return _Laplace.WrapperPolyData_GetCellNumberOfPoints(*args)
    def GetCellPoint(*args): return _Laplace.WrapperPolyData_GetCellPoint(*args)
    def GetPointX(*args): return _Laplace.WrapperPolyData_GetPointX(*args)
    def GetPointY(*args): return _Laplace.WrapperPolyData_GetPointY(*args)
    def GetPointZ(*args): return _Laplace.WrapperPolyData_GetPointZ(*args)
    def GetNumberOfPoints(*args): return _Laplace.WrapperPolyData_GetNumberOfPoints(*args)

class WrapperPolyDataPtr(WrapperPolyData):
    def __init__(self, this):
        _swig_setattr(self, WrapperPolyData, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, WrapperPolyData, 'thisown', 0)
        _swig_setattr(self, WrapperPolyData,self.__class__,WrapperPolyData)
_Laplace.WrapperPolyData_swigregister(WrapperPolyDataPtr)


GetPolyData = _Laplace.GetPolyData
class Line(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Line, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Line, name)
    def __repr__(self):
        return "<C Line instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Line, 'this', _Laplace.new_Line(*args))
        _swig_setattr(self, Line, 'thisown', 1)
    def __del__(self, destroy=_Laplace.delete_Line):
        try:
            if self.thisown: destroy(self)
        except: pass
    def AddPoint(*args): return _Laplace.Line_AddPoint(*args)
    def GetNumberOfPoints(*args): return _Laplace.Line_GetNumberOfPoints(*args)
    def GetPoint(*args): return _Laplace.Line_GetPoint(*args)

class LinePtr(Line):
    def __init__(self, this):
        _swig_setattr(self, Line, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Line, 'thisown', 0)
        _swig_setattr(self, Line,self.__class__,Line)
_Laplace.Line_swigregister(LinePtr)

class Locator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Locator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Locator, name)
    def __repr__(self):
        return "<C Locator instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Locator, 'this', _Laplace.new_Locator(*args))
        _swig_setattr(self, Locator, 'thisown', 1)
    def SetInputData(*args): return _Laplace.Locator_SetInputData(*args)
    def FindClosestPoint(*args): return _Laplace.Locator_FindClosestPoint(*args)
    def __del__(self, destroy=_Laplace.delete_Locator):
        try:
            if self.thisown: destroy(self)
        except: pass

class LocatorPtr(Locator):
    def __init__(self, this):
        _swig_setattr(self, Locator, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Locator, 'thisown', 0)
        _swig_setattr(self, Locator,self.__class__,Locator)
_Laplace.Locator_swigregister(LocatorPtr)


