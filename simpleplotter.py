import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
import math

def lin_fit(x, a, b):
    return a*x + b
def quad_fit(x, a, b, c):
    return a*x**2 + b*x + c #use this to fit polarizations
    
def fit(x, y, model):
    params, cov = curve_fit(model, x, y)
    fit_err = np.sqrt(np.diag(cov))
    x_fit = np.linspace(min(x), max(x), len(x))
    y_fit = model(x_fit, *params)
    
    err = np.sqrt(y)
    zero_indices = np.where(err <= 1e-6)[0]
    filtered_err = np.delete(err, zero_indices)
    filtered_y = np.delete(y, zero_indices)
    filtered_y_fit = np.delete(y_fit, zero_indices)

    chi2 = np.sum((filtered_y - filtered_y_fit)**2 / filtered_err**2)
    dof = len(filtered_y) - len(params)
    p = 1 - sp.stats.chi2.cdf(chi2, dof)

    return x_fit, y_fit, params, fit_err, chi2, dof, p



def get_eff(arr1, arr2, nbins=100):
    h1 = np.histogram(arr1, nbins, density=True)[0]
    h2, bins = np.histogram(arr2, nbins, density=True)
    err1 = np.sqrt(h1)
    err2 = np.sqrt(h2)
    mids = [(bins[i]+bins[i+1])/2 for i in range(0,nbins)]
    eff = []
    for i in range(nbins):
        if not h2[i] <= 1e-6:
            eff.append(h1[i]/h2[i])
        else:
            eff.append(0)
    err = []

    for i in range(len(h1)):
        if not (h1[i] <= 1e-6 or h2[i] <= 1e-6):
            err.append(eff[i]*(np.sqrt(err1[i]/h1[i])**2 + (err2[i]/h2[i])**2))
        else:
            err.append(0)
    
    #err = eff*np.sqrt((err1/h1)**2 + (err2/h2)**2)
    err = np.array(err)
    zero_indices = np.where(err <= 1e-6)[0]
    filtered_err = np.delete(err, zero_indices)
    filtered_mids = np.delete(mids, zero_indices)
    filtered_eff = np.delete(eff, zero_indices)

    return filtered_eff, filtered_mids, filtered_err





def makePlots(n):
    df_detector = pd.read_hdf('/work/submit/amsabbag/anglefit/data/detector', 'data')
    df_nodetector = pd.read_hdf('/work/submit/amsabbag/anglefit/data/nodetector', 'data')
    
    for name in df_detector:
       
        
        bin_counts, bin_edges = np.histogram(df_detector[name], bins=n)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(bin_counts)  # Poisson error (sqrt of count)
        plt.bar(bin_centers, bin_counts, width=np.diff(bin_edges), edgecolor='black', color='skyblue', label='Histogram')
        plt.errorbar(bin_centers, bin_counts, yerr=errors, fmt='none', ecolor='red', capsize=5, label='Error Bars')

        plt.title(name + ' detector')
        plt.savefig('/work/submit/amsabbag/anglefit/plots/detector/' + name + '.png')
        plt.clf()

        bin_counts, bin_edges = np.histogram(df_nodetector[name], bins=n)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(bin_counts)  # Poisson error (sqrt of count)
        plt.bar(bin_centers, bin_counts, width=np.diff(bin_edges), edgecolor='black', color='skyblue', label='Histogram')
        plt.errorbar(bin_centers, bin_counts, yerr=errors, fmt='none', ecolor='red', capsize=5, label='Error Bars')

        plt.title(name + ' no detector')
        plt.savefig('/work/submit/amsabbag/anglefit/plots/nodetector/' + name + '.png')
        plt.clf()

        if not name == 'mass':
            eff, mids, err =  get_eff(df_nodetector[name], df_detector[name])
            x_fit, y_fit, params, fit_err, chi2, dof, p = fit(mids, eff, lin_fit)

            plt.errorbar(mids, eff, yerr=err)
            plt.plot(x_fit, y_fit, color='red')
            equation_text = f'Fitted line: y = ({params[0]:.2f} ± {fit_err[0]:.2f})x + {params[1]:.2f} ± {fit_err[1]:.2f} \n χ² = {chi2} \n df = {dof} \n p = {p:.4f}'
            plt.text(0.95, 0.95, equation_text, transform=plt.gca().transAxes,
            fontsize=12, color='black', backgroundcolor='white', 
            ha='right', va='top', bbox=dict(facecolor='white', edgecolor='blue', boxstyle='round,pad=0.5'))

            plt.savefig('/work/submit/amsabbag/anglefit/plots/efficiencies/' + name + '_eff.png')
            plt.clf()



makePlots(50)
