#include "ReadTree.hh"
#include "json/json.h"
#include "tools.hh"

#include <cstring>
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"


int main( int argc, char** argv ) {

  if ( argc != 2 ) {
    std::cout << "Usage: angles.exe file.json \n";
    return ( 1 );
  }

  Json::Value config;
  std::ifstream input( argv[1] );
  input >> config;
  input.close();

  std::string proton = config["proton"].asString();
  std::string pion   = config["pion"].asString();
  std::string muon1  = config["muon1"].asString();
  std::string muon2  = config["muon2"].asString();
  std::string lambda0 = config["lambda0"].asString();


  // bool isMC            = config["isMC"].asBool();
  bool rapidsim        = config["rapidsim"].asBool();
  std::string partName = config["particle"].asString();
  std::string inName   = config["inFile"].asString();
  std::string outName  = config["outFile"].asString();
  std::string tree     = config["tree"].asString();

  TFile* in = new TFile( inName.c_str(), "read" );
  TTree* inTree = (TTree*)in->Get( tree.c_str() );
 
  // Activate all branches
  inTree->SetBranchStatus("*", 1);

  TFile out(outName.c_str(), "recreate");
  TTree* newTree = new TTree( tree.c_str(), "" );
  newTree->SetDirectory( &out );
  newTree = inTree->CloneTree();

  int proton_q, muon1_q, h1_id;
  double proton_px, proton_py, proton_pz;
  double pion_px, pion_py, pion_pz;
  double muon1_px, muon1_py, muon1_pz;
  double muon2_px, muon2_py, muon2_pz;
  std::vector<int> *proton_q_vec=new std::vector<int> ();
  std::vector<int> *muon1_q_vec=new std::vector<int> ();
  std::vector<float> *proton_px_vec=new std::vector<float> ();
  std::vector<float> *proton_py_vec=new std::vector<float> ();
  std::vector<float> *proton_pz_vec=new std::vector<float> ();
  std::vector<float> *pion_px_vec=new std::vector<float> ();
  std::vector<float> *pion_py_vec=new std::vector<float> ();
  std::vector<float> *pion_pz_vec=new std::vector<float> ();
  std::vector<float> *muon1_px_vec=new std::vector<float> ();
  std::vector<float> *muon1_py_vec=new std::vector<float> ();
  std::vector<float> *muon1_pz_vec=new std::vector<float> ();
  std::vector<float> *muon2_px_vec=new std::vector<float> ();
  std::vector<float> *muon2_py_vec=new std::vector<float> ();
  std::vector<float> *muon2_pz_vec=new std::vector<float> ();
  std::vector<int> *h1_id_vec=new std::vector<int> ();
  if (rapidsim){
    inTree->SetBranchAddress( (proton+"PX").c_str(), &proton_px );
    inTree->SetBranchAddress( (proton+"PY").c_str(), &proton_py );
    inTree->SetBranchAddress( (proton+"PZ").c_str(), &proton_pz );
    inTree->SetBranchAddress( (pion+"PX").c_str(), &pion_px );
    inTree->SetBranchAddress( (pion+"PY").c_str(), &pion_py );
    inTree->SetBranchAddress( (pion+"PZ").c_str(), &pion_pz );
    inTree->SetBranchAddress( (muon1+"PX").c_str(), &muon1_px );
    inTree->SetBranchAddress( (muon1+"PY").c_str(), &muon1_py );
    inTree->SetBranchAddress( (muon1+"PZ").c_str(), &muon1_pz );
    inTree->SetBranchAddress( (muon2+"PX").c_str(), &muon2_px );
    inTree->SetBranchAddress( (muon2+"PY").c_str(), &muon2_py );
    inTree->SetBranchAddress( (muon2+"PZ").c_str(), &muon2_pz );
    proton_q = 1;
    muon1_q  = 1;
  } else {
    inTree->SetBranchAddress( (proton+"px").c_str(), &proton_px_vec );
    inTree->SetBranchAddress( (proton+"py").c_str(), &proton_py_vec );
    inTree->SetBranchAddress( (proton+"pz").c_str(), &proton_pz_vec );
    inTree->SetBranchAddress( (pion+"px").c_str(), &pion_px_vec );
    inTree->SetBranchAddress( (pion+"py").c_str(), &pion_py_vec );
    inTree->SetBranchAddress( (pion+"pz").c_str(), &pion_pz_vec );
    inTree->SetBranchAddress( (muon1+"px").c_str(), &muon1_px_vec );
    inTree->SetBranchAddress( (muon1+"py").c_str(), &muon1_py_vec );
    inTree->SetBranchAddress( (muon1+"pz").c_str(), &muon1_pz_vec );
    inTree->SetBranchAddress( (muon2+"px").c_str(), &muon2_px_vec );
    inTree->SetBranchAddress( (muon2+"py").c_str(), &muon2_py_vec );
    inTree->SetBranchAddress( (muon2+"pz").c_str(), &muon2_pz_vec );
    inTree->SetBranchAddress( (proton+"q").c_str(), &proton_q_vec );
    inTree->SetBranchAddress( (muon1+"q").c_str(), &muon1_q_vec );
    inTree->SetBranchAddress( (proton+"type").c_str(), &h1_id_vec );
  }

  double ctPerp = 0;
  double ctPara = 0;
  double ctL = 0;
  double ctH = 0;
  double phiL = 0;
  double phiH = 0;
  double dphi = 0;
  double q2 = 0;
  double mppi = 0;

  const double mproton = 0.938272081;
  const double mpion = 0.13957061;
  const double mmuon = 0.1056583745;

  newTree->Branch("cosThetaPerp", &ctPerp, "cosThetaPerp/D");
  newTree->Branch("cosThetaPara", &ctPara, "cosThetaPara/D");
  newTree->Branch("cosThetaL", &ctL, "cosThetaL/D");
  newTree->Branch("cosThetaH", &ctH, "cosThetaH/D");
  newTree->Branch("phiL", &phiL, "phiL/D");
  newTree->Branch("phiH", &phiH, "phiH/D");
  newTree->Branch("dphi", &dphi, "dphi/D");
  newTree->Branch( "q2", &q2, "q2/D" );
  newTree->Branch( "mppi", &mppi, "mppi/D" );

  TLorentzVector initialProton;
  initialProton.SetPxPyPzE(0,0,1.,0);

  int nevents = newTree->GetEntries();
  if ( !config["events"].empty() )
  {
    nevents = config["events"].asInt64();
  }

  for ( int i = 0; i < nevents; ++i ) {

    if (i%1000 == 0) {std::cout<<"Processing event = "<<nevents-i<<std::endl;}

    newTree->GetEntry( i );

    if (!rapidsim){
      proton_q = proton_q_vec->at(0);
      h1_id = h1_id_vec->at(0);
      // h1 may be the proton or the pion
      if (h1_id == 2212){
        proton_px = proton_px_vec->at(0);
        proton_py = proton_py_vec->at(0);
        proton_pz = proton_pz_vec->at(0);
        pion_px = pion_px_vec->at(0);
        pion_py = pion_py_vec->at(0);
        pion_pz = pion_pz_vec->at(0);
      } else if (h1_id == 211){
        pion_px = proton_px_vec->at(0);
        pion_py = proton_py_vec->at(0);
        pion_pz = proton_pz_vec->at(0);
        proton_px = pion_px_vec->at(0);
        proton_py = pion_py_vec->at(0);
        proton_pz = pion_pz_vec->at(0);
        proton_q = -proton_q;
      } else {
        std::cerr<<"h1_id = "<<h1_id<<std::endl;
        exit(1);
      }
      muon1_q = muon1_q_vec->at(0);
      muon1_px = muon1_px_vec->at(0);
      muon1_py = muon1_py_vec->at(0);
      muon1_pz = muon1_pz_vec->at(0);
      muon2_px = muon2_px_vec->at(0);
      muon2_py = muon2_py_vec->at(0);
      muon2_pz = muon2_pz_vec->at(0);
    }

    TLorentzVector LbP4, JpsiP4, LP4, pP4, piP4, mupP4, mumP4, q_diff_P4;

    tools::getFourMomentum( proton_px, proton_py, proton_pz, pP4, mproton);
    tools::getFourMomentum( pion_px, pion_py, pion_pz, piP4, mpion);

    if (muon1_q>0)
    {
      tools::getFourMomentum(muon1_px, muon1_py, muon1_pz, mupP4, mmuon);
      tools::getFourMomentum(muon2_px, muon2_py, muon2_pz, mumP4, mmuon);
    }
    else
    {
      tools::getFourMomentum(muon2_px, muon2_py, muon2_pz, mupP4, mmuon);
      tools::getFourMomentum(muon1_px, muon1_py, muon1_pz, mumP4, mmuon);
    }
    LbP4 = pP4 + piP4 + mupP4 + mumP4;
    JpsiP4 = mupP4 + mumP4;
    LP4 = pP4 + piP4;

    tools::LbPsiRAngles(initialProton, LbP4, JpsiP4, LP4, mupP4, mumP4, pP4,
                        piP4, proton_q, ctPerp, ctPara, ctL, ctH, phiL, phiH, dphi);

    q2 = JpsiP4.M2();
    mppi = LP4.M();

    newTree->Fill();
  }
  out.cd();
  newTree->Write();
  // out.Write();
  newTree->SetDirectory(nullptr);
  out.Close();
  in->Close();
  // this is to avoid the segfault
  // delete reader;
  return 0;
}
