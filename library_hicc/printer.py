#!/usr/bin/env python3
"""
This file exists because Slurm doesn't report printed values in logs when the memory limit is exceeded,
and I would like to know where these exceptions are occurring so I need to have the logs.

Thus this prints out strings into a file by writing the strings to it, which works fine
"""
import sys

class Printer:
    """
    Give it a path to create a log file to write to, use the write methods
    to write to it.
    """
    def __init__(self, path):
        self.f = open(path, 'w')
        return
    
    def write(self, instring):
        self.f.write(instring+'\n')
        return
    
    def writeTab(self, instring):
        self.f.write('\t')
        self.write(instring)
        return
    
    def writeMem(self, name, obj):
        self.f.write("%s has size %.3e MB \n"%(name, sys.getsizeof(obj)/1e6))
        return