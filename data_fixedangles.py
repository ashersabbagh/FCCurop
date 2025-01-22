import math
import ROOT
from ROOT import TVector3, TFile, TCanvas, TLorentzVector
import uproot
import numpy as np
import pandas as pd





#NOTE: when using a new file, make sure to make the names array in this order!
names = ['muplus_px', 'muplus_py', 'muplus_pz', 'muplus_E', 'muplus_q', 'muminus_px', 'muminus_py', 'muminus_pz', 'muminus_E',
'proton_px', 'proton_py', 'proton_pz', 'proton_E', 'proton_q', 'pion_px', 'pion_py', 'pion_pz', 'pion_E']

names_detector = ['Muons_mu1px', 'Muons_mu1py', 'Muons_mu1pz', '(Muons_mu1p**2 + Muons_mu1m**2)**0.5', 'Muons_mu1q', 'Muons_mu2px', 
'Muons_mu2py', 'Muons_mu2pz', '(Muons_mu2p**2 + Muons_mu2m**2)**0.5', 'LCandidates_h1px', 'LCandidates_h1py', 'LCandidates_h1pz', 
'(LCandidates_h1p**2 + LCandidates_h1m**2)**0.5', 'LCandidates_h1q', 'LCandidates_h2px', 'LCandidates_h2py', 'LCandidates_h2pz', '(LCandidates_h2p**2 + LCandidates_h2m**2)**0.5']

names_nodetector = ['mup_0_PX_TRUE', 'mup_0_PY_TRUE', 'mup_0_PZ_TRUE', 'mup_0_E_TRUE', 'MISSING', 'mum_0_PX_TRUE', 'mum_0_PY_TRUE', 
'mum_0_PZ_TRUE', 'mum_0_E_TRUE', 'pp_0_PX_TRUE', 'pp_0_PY_TRUE', 'pp_0_PZ_TRUE', 'pp_0_E_TRUE', 'MISSING', 'pim_0_PX_TRUE', 'pim_0_PY_TRUE',
'pim_0_PZ_TRUE', 'pim_0_E_TRUE']

# names_detector = ['LCandidates_h1p', 'LCandidates_h1m', 'LCandidates_h1px', 'LCandidates_h1py', 'LCandidates_h1pz', 'LCandidates_h2p', 'LCandidates_h2m', 'LCandidates_h2px', 
# 'LCandidates_h2py', 'LCandidates_h2pz', 'LCandidates_h1q', 'LCandidates_mass', 'Muons_mu1p', 'Muons_mu1m', 'Muons_mu1px', 'Muons_mu1py', 
# 'Muons_mu1pz', 'Muons_mu1q', 'Muons_mu2p', 'Muons_mu2px', 'Muons_mu2py', 'Muons_mu2pz', 'Muons_mu2q', '(Muons_mu1p**2 + Muons_mu1m**2)**0.5', 
# '(Muons_mu2p**2 + Muons_mu2m**2)**0.5', '(LCandidates_h1p**2 + LCandidates_h1m**2)**0.5', '(LCandidates_h2p**2 + LCandidates_h2m**2)**0.5']


def collect(filename):
    costhetap_array = []
    costhetamu_array = []
    phi_array = []
    phimu_array = []
    phippi_array = []
    costhetab_array = []
    costhetat_array = []
    mass_array = []
    dataDict = {}



    
    with uproot.open('/work/submit/amsabbag/anglefit/data/' + filename) as f:
        if(filename == 'Lb2Lmm_tree.root'):
            type = 'nodetector'
            tree = f['DecayTree']
            for i in range(len(names)):
                if names_nodetector[i] == 'MISSING':
                    dataDict[names[i]] = np.full(len(dataDict[names[0]]), 1) 
                else:
                    dataDict[names[i]] = tree.arrays([names_nodetector[i]])[names_nodetector[i]]


        elif(filename == 'clean.root'):
            type = 'detector'
            tree = f['events']
            for i in range(len(names)):
                dataDict[names[i]] = np.array(tree.arrays([names_detector[i]])[names_detector[i]]).flatten()
      
        
        
    print('Getting data, this may take some time...')
    j = 0
    for i in range(len(dataDict[names[0]])):
        j += 1
        if j / len(dataDict[names[0]]) > 0.1:
            print(f'{100*round(i / len(dataDict[names[0]]), 2)}% complete...')
            j = 0

        #define particles and n direction:



        if(dataDict['muplus_q'][i] > 0):
            muplus = TLorentzVector(dataDict['muplus_px'][i], dataDict['muplus_py'][i], dataDict['muplus_pz'][i], dataDict['muplus_E'][i])
            muminus = TLorentzVector(dataDict['muminus_px'][i], dataDict['muminus_py'][i], dataDict['muminus_pz'][i], dataDict['muminus_E'][i])
        else:
            muminus = TLorentzVector(dataDict['muplus_px'][i], dataDict['muplus_py'][i], dataDict['muplus_pz'][i], dataDict['muplus_E'][i])
            muplus = TLorentzVector(dataDict['muminus_px'][i], dataDict['muminus_py'][i], dataDict['muminus_pz'][i], dataDict['muminus_E'][i])


        proton = TLorentzVector(dataDict['proton_px'][i], dataDict['proton_py'][i], dataDict['proton_pz'][i], dataDict['proton_E'][i])
        Lzero = dataDict['proton_q'][i] > 0
        pion = TLorentzVector(dataDict['pion_px'][i], dataDict['pion_py'][i], dataDict['pion_pz'][i], dataDict['pion_E'][i])

        
        
        #get particle groups
        mumu = muplus + muminus
        lambda0 = proton + pion
        lambdab = mumu + lambda0
        
        
        
        #define axes
        beam = TVector3(0, 0, 1)
        n = lambdab.Vect().Cross(beam)



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


        mumu_bf = TLorentzVector(mumu)
        mumu_bf.Boost(lambdabboost)

        muplus_bf = TLorentzVector(muplus)
        muplus_bf.Boost(lambdabboost)
        muminus_bf = TLorentzVector(muminus)
        muminus_bf.Boost(lambdabboost)

        proton_bf = TLorentzVector(proton)
        proton_bf.Boost(lambdabboost)
        pion_bf = TLorentzVector(pion)
        pion_bf.Boost(lambdabboost)

        #get polarization angles
        costhetab_array.append(math.cos(n.Angle(lambda0_bf.Vect())))
        costhetat_array.append(math.cos(lambdab.Vect().Angle(lambda0_bf.Vect())))



        #getting costhetamu:
        #NOTE: angle between dilepton in lb frame and muplus in dilepton frame 
        if (Lzero):
            costhetamu_array.append(math.cos(muplus_mumuf.Vect().Angle(mumu_bf.Vect())))
        else:
            costhetamu_array.append(math.cos(muminus_mumuf.Vect().Angle(mumu_bf.Vect())))
        
        #getting costhetap:
        costhetap_array.append(math.cos(proton_ppf.Vect().Angle(-mumu_bf.Vect())))

        #now getting the phi angles:
        # normalppi = proton_bf.Vect().Cross(pion_bf.Vect())
        # normalmumu = muplus_bf.Vect().Cross(muminus_bf.Vect())
        
        # if Lzero:
        #     phi = normalppi.Angle(normalmumu)
        #     if (normalmumu.Cross(normalppi).Dot(lambda0_bf.Vect()) < 0.0):
        #         phi = -phi
        # else:
        #     phi = normalppi.Angle(-normalmumu)
        #     if (normalmumu.Cross(normalppi).Dot(lambda0_bf.Vect()) < 0.0):
        #         phi = -phi
 
        # phi_array.append(phi)
        normalppi = proton_bf.Vect().Cross(pion_bf.Vect())
        normalmumu = muplus_bf.Vect().Cross(muminus_bf.Vect())
        
        phi_p = np.pi/2 - n.Angle(normalppi)

    
        if Lzero:
            phi_mu = np.pi/2 - n.Angle(normalmumu)
        else:
            phi_mu = np.pi/2 - n.Angle(-normalmumu)

 
        phimu_array.append(phi_mu)
        phippi_array.append(phi_p)
        phi_array.append(phi_mu - phi_p)
        
        mass_array.append(lambdab.Mag())

    dict = {}
    dict['costhetamu'] = costhetamu_array
    dict['costhetap'] = costhetap_array
    dict['costhetab'] = costhetab_array
    dict['costhetat'] = costhetat_array
    dict['phi'] = phi_array
    dict['phimu'] = phimu_array
    dict['phippi'] = phippi_array

    dict['mass'] = mass_array

    
    df = pd.DataFrame.from_dict(dict)
    df.to_hdf('/work/submit/amsabbag/anglefit/data/' + type, key='data', mode='w')
    print('Data complete.  ROOT file generated.')
    return (costhetamu_array, costhetap_array, phi_array, mass_array)



















