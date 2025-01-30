import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
import math

path = '/work/submit/amsabbag/anglefit/data/'

def getEquation(params, fiterr):
    eqn = 'y = '
    n = len(params)
    for i in range(n):
        if  n - i > 2:
            eqn += f'({params[i]:.2f} ± {fiterr[i]:.2f})x^{n-1-i} + '
        elif n - i == 2:
            eqn += f'({params[i]:.2f} ± {fiterr[i]:.2f})x + '
        else:   
            eqn += f'{params[i]:.2f} ± {fiterr[i]:.2f}'
    return eqn
        



def lin_fit(x, a, b):
    return a*x + b

def quad_fit(x, a, b, c):
    return a*x**2 + b*x + c #use this to fit polarizations
    
def fit(x, y, yerr, model): 
    params, cov = curve_fit(model, x, y, sigma=yerr, absolute_sigma=True)
    fit_err = np.sqrt(np.diag(cov))
    x_fit = np.linspace(min(x), max(x), len(x))
    y_fit = model(x_fit, *params)
    


    chi2 = np.sum((y - y_fit)**2 / yerr**2)
    dof = len(y) - len(params)
    p = 1 - sp.stats.chi2.cdf(chi2, dof)

    return x_fit, y_fit, params, fit_err, chi2, dof, p


def getEff(names, arrs=None, nbins=50, detectorfile='clean_detector', nodetectorfile='Lb2Lmm_tree_nodetector'):
    df_detector = pd.read_hdf(path + detectorfile, 'data')
    df_nodetector = pd.read_hdf(path + nodetectorfile, 'data')

    if arrs == None:
        arr1 = df_nodetector[names]
        arr2 = df_detector[names]
    else:
        arr1 = arrs[0]
        arr2 = arrs[1]

    print(len(arr1))
    print(len(arr2))
    if len(names)==1:
        
        h1 = np.histogram(arr1, nbins, density=False)[0]
        h2, bins = np.histogram(arr2, nbins, density=False)
        mids = [(bins[i]+bins[i+1])/2 for i in range(0,nbins)]

    elif len(names)==2:
        h1 = np.histogram2d(arr1[names[0]], arr1[names[1]], nbins, density=False)[0]
        h2, xedges, yedges = np.histogram2d(arr2[names[0]], arr2[names[1]], nbins, density=False)

        mids = [[(xedges[i]+xedges[i+1])/2 for i in range(0,nbins)], [(yedges[i]+yedges[i+1])/2 for i in range(0,nbins)]]



    h1_nozeros = np.where(h1<1, np.ones_like(h1), h1)
    h2_nozeros = np.where(h2<1, np.ones_like(h2), h2)
    err1 = np.sqrt(h1)
    err2 = np.sqrt(h2)

    eff = h2_nozeros / h1_nozeros
    
    err = eff*(np.sqrt((err1/h1_nozeros)**2 + (err2/h2_nozeros)**2))

    return eff, mids, err

def get2DEff(name1, name2):
    eff, x, err = getEff([name1, name2])

    plt.figure(figsize=(8, 6))
    plt.imshow(eff, extent=[min(x[0]), max(x[0]), min(x[1]), max(x[1])], origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(label='total efficiency')
    plt.xlabel(name2)
    plt.ylabel(name1)
    plt.title(name2 + ' X ' + name1)
    plt.savefig('/work/submit/amsabbag/anglefit/plots/2D_efficiencies/' + name1 + '/' + name1 + 'X' + name2 + '_eff.png')
    plt.clf()



def makeAnglePlots(type, file, n=50):
    df = pd.read_hdf(path + file, 'data')


    for name in df:
        bin_counts, bin_edges = np.histogram(df[name], bins=n)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(bin_counts)  # Poisson error (sqrt of count)
        plt.bar(bin_centers, bin_counts, width=np.diff(bin_edges), edgecolor='black', color='skyblue', label='Histogram')
        plt.errorbar(bin_centers, bin_counts, yerr=errors, fmt='none', ecolor='red', capsize=5, label='Error Bars')

        plt.title(name + ' ' + type)
        if file == 'Lb2Lmm_phsp_tree_nodetector':
            plt.savefig('/work/submit/amsabbag/anglefit/plots/' + type + '/phasespace/' + name + '_' + file + '.png')
        else:
            plt.savefig('/work/submit/amsabbag/anglefit/plots/' + type + '/' + name + '_' + file + '.png')
        plt.clf()

def makeEffPlots(nbins=50, detectorfile='clean_detector', nodetectorfile='Lb2Lmm_tree_nodetector', q2Binning=False):
    df_detector = pd.read_hdf(path + detectorfile, 'data')

    if not q2Binning:
        for name in df_detector:
            if not (name == 'mass' or name == 'q2'):


                eff, mids, err =  getEff([name], nbins=nbins, detectorfile=detectorfile, nodetectorfile=nodetectorfile)
                x_fit, y_fit, params, fit_err, chi2, dof, p = fit(mids, eff, err, quad_fit)

                plt.errorbar(mids, eff, yerr=err, fmt='ko')
                plt.ylim(0, 1.5*max(y_fit))
                plt.plot(x_fit, y_fit, color='red')
                equation_text = getEquation(params, fit_err) + f'\n χ² = {chi2} \n df = {dof} \n p = {p:.5f}'
                plt.text(0.95, 0.95, equation_text, transform=plt.gca().transAxes,
                fontsize=12, color='black', backgroundcolor='white', 
                ha='right', va='top', bbox=dict(facecolor='white', edgecolor='blue', boxstyle='round,pad=0.5'))

                plt.savefig('/work/submit/amsabbag/anglefit/plots/efficiencies/' + name + '_eff.png')
                plt.clf()
    else:
        df_detector = pd.read_hdf(path + detectorfile, 'data')
        df_nodetector = pd.read_hdf(path + nodetectorfile, 'data')
        q2 = df_nodetector['q2']
        q2Bins = [[0,2], [2,4], [4,6], [6,8], [8,11], [11,12.5], [12.5,15], [15,20]]
        for name in df_detector:
            if not (name == 'mass' or name == 'q2'):
                for bins in q2Bins:
                    indices = np.where((q2 >= bins[0]) & (q2 <= bins[1]))
                    nodet_arr = df_nodetector[name].to_numpy()[indices]
                    det_arr = df_detector[name].to_numpy()[indices]
                    arrs = [nodet_arr, det_arr]
                    eff, mids, err =  getEff([name], arrs=arrs, nbins=n, detectorfile=detectorfile, nodetectorfile=nodetectorfile)
                    x_fit, y_fit, params, fit_err, chi2, dof, p = fit(mids, eff, err, quad_fit)

                    plt.errorbar(mids, eff, yerr=err, fmt='ko')
                    plt.ylim(0, 1.5*max(y_fit))
                    plt.plot(x_fit, y_fit, color='red')
                    equation_text = getEquation(params, fit_err) + f'\n χ² = {chi2} \n df = {dof} \n p = {p:.5f}'
                    plt.text(0.95, 0.95, equation_text, transform=plt.gca().transAxes,
                    fontsize=12, color='black', backgroundcolor='white', 
                    ha='right', va='top', bbox=dict(facecolor='white', edgecolor='blue', boxstyle='round,pad=0.5'))

                    plt.savefig('/work/submit/amsabbag/anglefit/plots/efficiencies/q2bins/' + str(bins) + '/' + name + '_eff.png')
                    plt.clf()

                    



def make2DEffPlots(detectorfile='clean_detector'):
    df_detector = pd.read_hdf(path + detectorfile, 'data')
    for name1 in df_detector:
        if not (name1 == 'mass'):
            for name2 in df_detector:
                if not (name2 == 'mass' or name1 == name2):
                    get2DEff(name1, name2)



def makeAllPlots(n=50):
    makeAnglePlots(n)
    makeEffPlots(n)
    make2DEffPlots()



#makeAllPlots()
#make2DEffPlots()
#makeEffPlots()
#makeAnglePlots(type='nodetector', file='Lb2Lmm_phsp_tree_nodetector')
#makeAnglePlots(type='nodetector', file='Lb2Lmm_tree_1M_nodetector')
#makeAnglePlots(type='detector', file='clean_detector')

makeEffPlots(q2Binning=False, nodetectorfile='Lb2Lmm_tree_1M_nodetector', nbins=20)
makeAnglePlots(type='nodetector', file='Lb2Lmm_tree_1M_nodetector')
