import cx_Freeze
import sys
import os
import importlib
importlib.import_module('mpl_toolkits.mplot3d').Axes3D

base = None


executables = [cx_Freeze.Executable("Electric Field API.py", base=base, icon=os.getcwd()+"\\bin\\EFAPIicon.ico")]

cx_Freeze.setup(
    name = "Electric Field API",
    options = {"build_exe": {'includes': ['numpy.core._methods','numpy.lib.format','tkFileDialog','FileDialog'], 'packages': ['matplotlib.pyplot',"mpl_toolkits.mplot3d","matplotlib",'Tkinter','FileDialog','tkFileDialog'], "include_files":[os.getcwd()+"\\bin\\EFAPIicon.ico"]}},
    version = "1.3",
    description = "Electric Field Visualization",
    executables = executables
    )