
def GetHarnessRanges():
    HarnessRanges=[]
    #module +- 4,3,2,1
    ietamin=6
    ietamax=25
    iphimin=1
    iphimax=10
    for keta in range (0,4): #loop over eta
        for kphi in range (0,36):
            HarnessRange = [ietamin+keta*20, ietamax+keta*20, iphimin+kphi*10, iphimax+kphi*10]#EB+
            HarnessRanges.append(HarnessRange)
            HarnessRange = [-(ietamax+keta*20), -(ietamin+keta*20), iphimin+kphi*10, iphimax+kphi*10]#EB-
            HarnessRanges.append(HarnessRange) 

    #module +- 0
    ietamin=1
    ietamax=5
    iphimin=1
    iphimax=20
    for kphi in range (0,18):
        HarnessRange = [ietamin, ietamax, iphimin+kphi*20, iphimax+kphi*20]#EB+
        HarnessRanges.append(HarnessRange)
        HarnessRange = [-ietamax, -ietamin, iphimin+kphi*20, iphimax+kphi*20]#EB-
        HarnessRanges.append(HarnessRange)

    return HarnessRanges

def GetModuleRanges():
    ModuleRanges=[]
    #module +- 4,3,2,1
    ietamin=6
    ietamax=25
    iphimin=1
    iphimax=360
    for keta in range (0,4): #loop over eta
        ModuleRange = [ietamin+keta*20, ietamax+keta*20, iphimin, iphimax]#EB+
        ModuleRanges.append(ModuleRange)
        ModuleRange = [-(ietamax+keta*20), -(ietamin+keta*20), iphimin, iphimax]#EB-
        ModuleRanges.append(ModuleRange) 

    #module +- 0
    ModuleRange = [1, 5, iphimin, iphimax]#EB+
    ModuleRanges.append(ModuleRange)
    ModuleRange = [-5, -1, iphimin, iphimax]#EB-
    ModuleRanges.append(ModuleRange)

    return ModuleRanges


