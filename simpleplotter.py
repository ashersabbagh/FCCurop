import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_eff(arr1, arr2, nbins=100):
    h1 = np.histogram(arr1, nbins)[0]
    h2, bins = np.histogram(arr2, nbins)
    mids = [(bins[i]+bins[i+1])/2 for i in range(0,nbins)]
    eff = []
    for i in range(nbins):
        if h2[i] != 0:
            eff.append(h1[i]/h2[i])
        else:
            eff.append(0)
    return eff, mids





def makePlots():
    df_detector = pd.read_hdf('/work/submit/amsabbag/anglefit/data/detector', 'data')
    df_nodetector = pd.read_hdf('/work/submit/amsabbag/anglefit/data/nodetector', 'data')
    
    for name in df_detector:
       
        plt.hist(df_detector[name], bins=500)
        plt.savefig('/work/submit/amsabbag/anglefit/plots/detector/' + name + '.png')
        plt.clf()

        
        plt.hist(df_nodetector[name], bins=500)
        plt.savefig('/work/submit/amsabbag/anglefit/plots/nodetector/' + name + '.png')
        plt.clf()

        eff, mids =  get_eff(df_nodetector[name], df_detector[name])
        plt.plot(mids, eff)
        plt.savefig('/work/submit/amsabbag/anglefit/plots/efficiencies/' + name + '_eff.png')
        plt.clf()



makePlots()
