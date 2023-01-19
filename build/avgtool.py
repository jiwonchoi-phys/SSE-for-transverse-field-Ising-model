import numpy as np
import glob

def AverageAll(filename_pattern):
    flist = glob.glob(filename_pattern)
    numcol = np.loadtxt(flist[0]).shape[-1]
    avglist = np.zeros([len(flist),numcol])
    for idx,f in enumerate(flist):
        avglist[idx,:] = np.mean(np.loadtxt(f),axis=0)
    return avglist

def AverageAndExport(filename_pattern,export_filename,fmt='%.6e'):
    avglist = AverageAll(filename_pattern)
    np.savetxt(export_filename,avglist,fmt=fmt)

if __name__=="__main__":
    AverageAndExport("sample0/*dS","sample0avg")
    #AverageAndExport("sample0/*dS")
