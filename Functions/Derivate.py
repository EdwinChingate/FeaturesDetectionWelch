import numpy as np
def Derivate(RTVec,IntVec):
    #RTord=ChromData[:,0].argsort()
    #RTVec=ChromData[:,0][RTord]
    #IntVec=ChromData[:,9][RTord]
    dRT=RTVec[2:]-RTVec[:-2]
    dInt=IntVec[2:]-IntVec[:-2]
    dS=dInt/dRT
    return [RTVec[1:-1],dS]
