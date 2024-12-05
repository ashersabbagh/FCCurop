import math
import ROOT
from ROOT import TVector3, TFile, TCanvas, TLorentzVector
import uproot
import numpy as np
import pandas as pd






names_nodetector = ['Lambda0_TRUEP_E', 'Lambda0_TRUEP_X', 'Lambda0_TRUEP_Y', 'Lambda0_TRUEP_Z', 'Lambda_b0_TRUEP_E', 'Lambda_b0_TRUEP_X', 'Lambda_b0_TRUEP_Y', 
'Lambda_b0_TRUEP_Z', 'muplus_TRUEP_X', 'muplus_TRUEP_Y', 'muplus_TRUEP_Z', 'muplus_TRUEP_E', 'muminus_TRUEP_X', 
'muminus_TRUEP_Y', 'muminus_TRUEP_Z', 'muminus_TRUEP_E', 'pplus_TRUEP_X', 'pplus_TRUEP_Y', 'pplus_TRUEP_Z', 
'pplus_TRUEP_E', 'piminus_TRUEP_X', 'piminus_TRUEP_Y', 'piminus_TRUEP_Z', 'piminus_TRUEP_E', 'pplus_ID', 'piminus_ID', 'muplus_ID', 'muminus_ID']

names_detector = ['LCandidates_h1p', 'LCandidates_h1m', 'LCandidates_h1px', 'LCandidates_h1py', 'LCandidates_h1pz', 'LCandidates_h2p', 'LCandidates_h2m', 'LCandidates_h2px', 
'LCandidates_h2py', 'LCandidates_h2pz', 'LCandidates_h1q', 'LCandidates_mass', 'Muons_mu1p', 'Muons_mu1m', 'Muons_mu1px', 'Muons_mu1py', 
'Muons_mu1pz', 'Muons_mu1q', 'Muons_mu2p', 'Muons_mu2px', 'Muons_mu2py', 'Muons_mu2pz', 'Muons_mu2q', '(Muons_mu1p**2 + Muons_mu1m**2)**0.5', 
'(Muons_mu2p**2 + Muons_mu2m**2)**0.5', '(LCandidates_h1p**2 + LCandidates_h1m**2)**0.5', '(LCandidates_h2p**2 + LCandidates_h2m**2)**0.5']



def collect(filename):
    costhetap_array = []
    costhetamu_array = []
    phip_array = []
    phimu_array = []
    phitot_array = []
    costhetab_array = []
    costhetat_array = []


    
    with uproot.open('/work/submit/amsabbag/anglefit/data/' + filename) as f:
        names = ''
        if(filename == 'LbToLzMuMu.root'):
            type = 'nodetector'
            names = names_nodetector
            tree = f['MCDecayTreeTuple/MCDecayTree;1']
        elif(filename == 'clean.root'):
            type = 'detector'
            names = names_detector
            tree = f['events']
      
        
        dataDict = {}
        for n in names:
            dataDict[n] = tree.arrays([n])[n]
    
    for i in range(len(dataDict[names[0]])):
        if(filename == 'LbToLzMuMu.root'):
            #define particles and n direction

            pmass = np.sqrt(dataDict['pplus_TRUEP_E'][i]**2 - (dataDict['pplus_TRUEP_X'][i]**2 + dataDict['pplus_TRUEP_Y'][i]**2 + dataDict['pplus_TRUEP_Z'][i]**2))


            if(dataDict['muplus_ID'][i] > 0):
                muplus = TLorentzVector(dataDict['muplus_TRUEP_X'][i], dataDict['muplus_TRUEP_Y'][i], dataDict['muplus_TRUEP_Z'][i], dataDict['muplus_TRUEP_E'][i])
                muminus = TLorentzVector(dataDict['muminus_TRUEP_X'][i], dataDict['muminus_TRUEP_Y'][i], dataDict['muminus_TRUEP_Z'][i], dataDict['muminus_TRUEP_E'][i])
            else:
                muminus = TLorentzVector(dataDict['muplus_TRUEP_X'][i], dataDict['muplus_TRUEP_Y'][i], dataDict['muplus_TRUEP_Z'][i], dataDict['muplus_TRUEP_E'][i])
                muplus = TLorentzVector(dataDict['muminus_TRUEP_X'][i], dataDict['muminus_TRUEP_Y'][i], dataDict['muminus_TRUEP_Z'][i], dataDict['muminus_TRUEP_E'][i])

            
            if(pmass > 0.5):
                proton = TLorentzVector(dataDict['pplus_TRUEP_X'][i], dataDict['pplus_TRUEP_Y'][i], dataDict['pplus_TRUEP_Z'][i], dataDict['pplus_TRUEP_E'][i])
                Lzero = dataDict['pplus_ID'][i] > 0
                pion = TLorentzVector(dataDict['piminus_TRUEP_X'][i], dataDict['piminus_TRUEP_Y'][i], dataDict['piminus_TRUEP_Z'][i], dataDict['piminus_TRUEP_E'][i])
            else:
                print("pion")
                pion = TLorentzVector(dataDict['pplus_TRUEP_X'][i], dataDict['pplus_TRUEP_Y'][i], dataDict['pplus_TRUEP_Z'][i], dataDict['pplus_TRUEP_E'][i])
                Lzero = dataDict['pplus_ID'][i] < 0
                proton = TLorentzVector(dataDict['piminus_TRUEP_X'][i], dataDict['piminus_TRUEP_Y'][i], dataDict['piminus_TRUEP_Z'][i], dataDict['piminus_TRUEP_E'][i])
        
        
        elif(filename == 'clean.root'):
            #define particles and n direction




            if(dataDict['Muons_mu1q'][i][0] > 0):
                muplus = TLorentzVector(dataDict['Muons_mu1px'][i][0], dataDict['Muons_mu1py'][i][0], dataDict['Muons_mu1pz'][i][0], dataDict['(Muons_mu1p**2 + Muons_mu1m**2)**0.5'][i][0])
                muminus = TLorentzVector(dataDict['Muons_mu2px'][i][0], dataDict['Muons_mu2py'][i][0], dataDict['Muons_mu2pz'][i][0], dataDict['(Muons_mu2p**2 + Muons_mu2m**2)**0.5'][i][0])
            else:
                muminus = TLorentzVector(dataDict['Muons_mu1px'][i][0], dataDict['Muons_mu1py'][i][0], dataDict['Muons_mu1pz'][i][0], dataDict['(Muons_mu1p**2 + Muons_mu1m**2)**0.5'][i][0])
                muplus = TLorentzVector(dataDict['Muons_mu2px'][i][0], dataDict['Muons_mu2py'][i][0], dataDict['Muons_mu2pz'][i][0], dataDict['(Muons_mu2p**2 + Muons_mu2m**2)**0.5'][i][0])

            
            if(dataDict['LCandidates_h1m'][0] > 0.5):
                proton = TLorentzVector(dataDict['LCandidates_h1px'][i][0], dataDict['LCandidates_h1py'][i][0], dataDict['LCandidates_h1pz'][i][0], dataDict['(LCandidates_h1p**2 + LCandidates_h1m**2)**0.5'][i][0])
                Lzero = dataDict['LCandidates_h1q'][i][0] > 0
                pion = TLorentzVector(dataDict['LCandidates_h2px'][i][0], dataDict['LCandidates_h2py'][i][0], dataDict['LCandidates_h2pz'][i][0], dataDict['(LCandidates_h2p**2 + LCandidates_h2m**2)**0.5'][i][0])
            else:
                pion = TLorentzVector(dataDict['LCandidates_h1px'][i][0], dataDict['LCandidates_h1py'][i][0], dataDict['LCandidates_h1pz'][i][0], dataDict['(LCandidates_h1p**2 + LCandidates_h1m**2)**0.5'][i][0])
                Lzero = dataDict['LCandidates_h1q'][i][0] < 0
                proton = TLorentzVector(dataDict['LCandidates_h2px'][i][0], dataDict['LCandidates_h2py'][i][0], dataDict['LCandidates_h2pz'][i][0], dataDict['(LCandidates_h2p**2 + LCandidates_h2m**2)**0.5'][i][0])



        #get particle groups
        mumu = muplus + muminus
        lambda0 = proton + pion
        lambdab = mumu + lambda0
        lambdab3 = TVector3(lambda0[0], lambda0[1], lambda0[2])
        
        
        
        #define axes
        beam = TVector3(0, 0, 1)
        beamvec = TLorentzVector(0, beam[0], beam[1], beam[2])
        n = beam.Cross(lambdab3) / beam.Cross(lambdab3).Mag()
        nvec = TLorentzVector(0, n[0], n[1], n[2])



        #get boost functions
        mumuboost = TVector3(-mumu.BoostVector())
        lambda0boost = TVector3(-lambda0.BoostVector())
        lambdabboost = TVector3(-lambdab.BoostVector())

        #boost into appropriate frames
        muplus_mumuf = TLorentzVector(muplus)
        muplus_mumuf.Boost(mumuboost)
        muminus_mumuf = TLorentzVector(muminus)
        muminus_mumuf.Boost(mumuboost)

        proton_ppf = TLorentzVector(proton)
        proton_ppf.Boost(lambda0boost)

        lambda0_bf = TLorentzVector(lambda0)
        lambda0_bf.Boost(lambdabboost)

        lambdab_mumuf = TLorentzVector(lambdab)
        lambdab_mumuf.Boost(mumuboost)
        lambdab_ppf = TLorentzVector(lambdab)
        lambdab_ppf.Boost(lambda0boost)

        muplus_bf = TLorentzVector(muplus)
        muplus_bf.Boost(lambdabboost)
        muminus_bf = TLorentzVector(muminus)
        muminus_bf.Boost(lambdabboost)

        proton_bf = TLorentzVector(proton)
        proton_bf.Boost(lambdabboost)
        pion_bf = TLorentzVector(pion)
        pion_bf.Boost(lambdabboost)

        #get polarization angles
        costhetab_array.append(math.cos(nvec.Vect().Angle(lambda0_bf.Vect())))
        costhetat_array.append(math.cos(lambdab.Vect().Angle(lambda0_bf.Vect())))



        #getting costhetamu:
        if (Lzero):
            costhetamu_array.append(math.cos(muplus_mumuf.Vect().Angle(-lambdab_mumuf.Vect())))
        else:
            costhetamu_array.append(math.cos(muminus_mumuf.Vect().Angle(-lambdab_mumuf.Vect())))
        
        #getting costhetap:
        costhetap_array.append(math.cos(proton_ppf.Vect().Angle(-lambdab_ppf.Vect())))

        #now getting the phi angles:
        normalppi = proton_bf.Vect().Cross(pion_bf.Vect())
        normalmumu = muplus_bf.Vect().Cross(muminus_bf.Vect())
        
        phi_p = nvec.Vect().Angle(normalppi)
        if (normalppi.Cross(nvec.Vect()).Dot(lambda0_bf.Vect()) < 0.0):
                phi_p = -phi_p
    
        if Lzero:
            phi_mu = nvec.Vect().Angle(normalmumu)
            if (normalmumu.Cross(nvec.Vect()).Dot(lambda0_bf.Vect()) < 0.0):
                phi_mu = -phi_mu
        else:
            phi_mu = nvec.Vect().Angle(-normalmumu)
            if (normalmumu.Cross(nvec.Vect()).Dot(lambda0_bf.Vect()) < 0.0):
                phi_mu = -phi_mu
 
        phimu_array.append(phi_mu)
        phip_array.append(phi_p)
        phitot_array.append(phi_mu + phi_p)


    dict = {}
    dict['costhetamu'] = costhetamu_array
    dict['costhetap'] = costhetap_array
    dict['costhetab'] = costhetab_array
    dict['costhetat'] = costhetat_array
    dict['phimu'] = phimu_array
    dict['phip'] = phip_array
    dict['phitot'] = phitot_array
    
    df = pd.DataFrame.from_dict(dict)
    df.to_hdf('/work/submit/amsabbag/anglefit/data/' + type, key='data', mode='w')

    return (costhetamu_array, costhetap_array, phitot_array)



















