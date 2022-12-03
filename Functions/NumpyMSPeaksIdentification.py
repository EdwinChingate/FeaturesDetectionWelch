#This one is just beautifull! Next step will be to make it a function, and to integrate all data from different times
import numpy as np
from Derivate import *
from PondMZStats import *
from WelchTest import *
import matplotlib.pyplot as plt
def NumpyMSPeaksIdentification(RawSignals,NoiseTresInt=5000,MinTresRelDer=3e5,minMZbetweenPeaks=1e-3,MinInttobePeak=100,MinSignalstobePeak=4,MinPeaksSpectra=3,r2Filter=0.8,ConfidenceIntervalTolerance=30,MinMZ=0,MaxMZ=1500):    
    #Include the interval to analize, then I can run short analises for only one substance
    RawSignals=np.array(RawSignals)    
    dimen=np.shape(RawSignals)    
    if dimen[0]<dimen[1]:
        RawSignals=RawSignals.T     
    DenoisedLoc=np.where((RawSignals[:,1]>NoiseTresInt)&(RawSignals[:,0]>MinMZ)&(RawSignals[:,0]<MaxMZ))[0]
    DenoisedSignals=RawSignals[DenoisedLoc,:].copy()    
    
    
    
    dS=Derivate(DenoisedSignals[:,0],DenoisedSignals[:,1])

    SlocNeg=np.where(dS[1]<-MinTresRelDer)[0]       
    #DifMZNeg=DenoisedSignals[SlocNeg,0]
    DifMZNeg0=DenoisedSignals[SlocNeg,0]
    DifMZNeg0=np.append(DifMZNeg0,max(DenoisedSignals[:,0])+1)
    DifMZNeg=SlocNeg
    #DifMZNeg=np.append(DifMZNeg,max(DenoisedSignals[:,0])+1)
    DifMZNeg=np.append(DifMZNeg,len(DenoisedSignals[:,0])+3)    
    DifSNeg=DifMZNeg[1:]-DifMZNeg[:-1]       
    #DifSNegLoc=np.where((DifSNeg>minMZbetweenPeaks))[0]
    DifSNegLoc=np.where((DifSNeg>3))[0]
    MinMZ=min(DenoisedSignals[:,0])-0.1
    FirstPeak=True
    SpectrumPeaks=[]
    for mzp in DifSNegLoc:
        leftMZ=DifMZNeg0[:-1][mzp]
        rightMZ=DifMZNeg0[1:][mzp]
        ValleyMZLoc=np.where((DenoisedSignals[:,0]>=leftMZ)&(DenoisedSignals[:,0]<=rightMZ))[0]
        ValleyMZ=DenoisedSignals[ValleyMZLoc,:]
        if len(ValleyMZLoc)>1:                  
            minIntValley=min(ValleyMZ[:,1])
            minValleyLoc=np.where(ValleyMZ[:,1]==minIntValley)[0]
            MaxMZ=np.mean(ValleyMZ[minValleyLoc,0])
        else:
            MaxMZ=(leftMZ+rightMZ)/2
        print(MinMZ,MaxMZ)
        PeakLoc=np.where((DenoisedSignals[:,0]>MinMZ)&(DenoisedSignals[:,0]<MaxMZ))[0]
        PeakData=DenoisedSignals[PeakLoc,:]   
        if len(PeakData)>MinSignalstobePeak:
            PeakStats=PondMZStats(PeakData)
            SaveMinMZ=PeakStats[0]-PeakStats[3]
            SaveMaxMZ=PeakStats[0]+PeakStats[3]
            PeakStats.append(SaveMinMZ)
            PeakStats.append(SaveMaxMZ)
            if FirstPeak:     
                PeakStats.append(0)
                PeakStats.append(0)
                PeakStats.append(0)
                SpectrumPeaks.append(PeakStats)
                FirstPeak=False
            else:
                WelchVec=WelchTest(SpectrumPeaks[-1],PeakStats,alpha=0.01)
                PeakStats.append(WelchVec[1])
                PeakStats.append(WelchVec[2])
                PeakStats.append(WelchVec[3]) 
                if WelchVec[0]:
                    SpectrumPeaks.append(PeakStats)
                else:
                    MinMZ=SpectrumPeaks[-1][-5]
                    MaxMZ=SpectrumPeaks[-1][-4]
                    print(MinMZ,MaxMZ,WelchVec)
                    PeakLoc=np.where((DenoisedSignals[:,0]>MinMZ)&(DenoisedSignals[:,0]<MaxMZ))[0]
                    PrevDat=DenoisedSignals[PeakLoc,:]
                    PeakData=np.append(PrevDat,PeakData,axis=0)     
                    PeakStats=PondMZStats(PeakData)
                    PeakStats.append(MinMZ)
                    PeakStats.append(MaxMZ)
                    if len(SpectrumPeaks)>2:
                        WelchVec=WelchTest(SpectrumPeaks[-2],PeakStats,alpha=0.01)
                        PeakStats.append(WelchVec[1])
                        PeakStats.append(WelchVec[2])
                        PeakStats.append(WelchVec[3])  
                    else:
                        PeakStats.append(0)
                        PeakStats.append(0)
                        PeakStats.append(0)                        
                    SpectrumPeaks[-1]=PeakStats
        MinMZ=MaxMZ
    if len(SpectrumPeaks)<MinPeaksSpectra:
        return 0
    SpectrumPeaks=np.array(SpectrumPeaks)    
    LastFilterLoc=np.where((SpectrumPeaks[:,8]>MinInttobePeak)&(SpectrumPeaks[:,7]>r2Filter)&(SpectrumPeaks[:,5]<0)&(SpectrumPeaks[:,8]>0)&(SpectrumPeaks[:,4]<ConfidenceIntervalTolerance))[0]
    return SpectrumPeaks[LastFilterLoc,:]  
