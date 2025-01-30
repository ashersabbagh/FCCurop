#include "ReadTree.hh"
#include "json/json.h"
#include "tools.hh"

#include <cstring>
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"


int main( int argc, char** argv ) {

  // Arguments will have json file with list of info. Treat MC and data separately
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


  bool isMC            = config["isMC"].asBool();
  std::string partName = config["particle"].asString();
  // unsigned int year    = config["year"].asUInt();
  // Polarity: MgUp = -1, MgDn = 1; from tools/src/BeamCrossing.cc
  // int polarity         = config["polarity"].asInt();
  std::string inName   = config["inFile"].asString();
  std::string outName  = config["outFile"].asString();

  // std::string inputDir = std::string(std::getenv("LBANAMASTER")) + "/";
  // if (isMC) { inputDir = std::string(std::getenv("LBANAMC")) + "/"; }
  TreeReader* reader = new TreeReader( "DecayTree" );
  reader->AddFile(inName);
  reader->Initialize();

  // std::string outputDir = std::string(std::getenv("LBANAINTERMEDIATE")) + "/";
  TFile out(outName.c_str(), "recreate");
  TTree* newTree = new TTree( "DecayTree", "" );
  reader->BranchNewTree( newTree );

  bool isLb = true;
  if ( partName == "B0" ) {
    std::cout<<"PARTICLE is B0"<<std::endl;
    isLb = false;
  }

  // Anja did this
  // bool B03D=false;
  // if ( !config["B0_3D"].empty() )
  // {
  //   B03D = config["B0_3D"].asBool();
  //   std::cout << "B03D = " << (int)B03D << std::endl;
  //   if (B03D == 0) {std::cout << "Going to do 5D angles" << std::endl;}
  // }

  bool constrained = true;
  if ( !config["constrained"].empty() )
  {
    constrained = config["constrained"].asBool();
  }

  // Anja did this because FCC collisions are head-on, right?
  // bool BeamCrossing=true;
  // if ( !config["beamCrossing"].empty() ) {
  //   BeamCrossing = config["beamCrossing"].asBool();
  // }

  // Create new branches
  // reco
  double ctPerp;
  double ctPara;
  double ctL;
  double ctH;
  double phiL;
  double phiH;
  double dphi;

  //true
  double ctPerpT;
  double ctParaT;
  double ctLT;
  double ctHT;
  double phiLT;
  double phiHT;
  double dphiT;
  double true_q2;
  double true_mpk;

  //constrained momentum
  double cctPerp;
  double cctPara;
  double cctL;
  double cctH;
  double cphiL;
  double cphiH;
  double cdphi;

  double q2;
  // Anja did this
  // bool l0;
  // bool hlt1;
  // bool hlt2;
  // bool trigger;
  // double q2_diff;

  /////////////////////////////////////////////////////////////////////////////////////
  //this is normaly what is used for the B02JpsiKS with 3 angles
  // Anja did this
  // if(B03D)
  // {
  //   //////////////////////////////////////////////////////////////////////
  //   // reco
  //   newTree->Branch("cosThetaL", &ctL, "cosThetaL/D");
  //   newTree->Branch("cosThetaH", &ctH, "cosThetaH/D");
  //   newTree->Branch("dphi", &dphi, "dphi/D");
  //   if ( isMC )
  //   {
  //     ///////////////////////////////////////////////////////////////////
  //     //true
  //     newTree->Branch("true_cosThetaL", &ctLT, "true_cosThetaL/D");
  //     newTree->Branch("true_cosThetaH", &ctHT, "true_cosThetaH/D");
  //     newTree->Branch("true_dphi", &dphiT, "true_dphi/D");
  //   }
  // }
  // else
  {
  ////////////////////////////////////////////////////////////////////////////////////
  //for Lb2Lmumu and for Lb2JpsiL
    ///////////////////////////////////////////////////////////////////////////
    //reco
    newTree->Branch("cosThetaPerp", &ctPerp, "cosThetaPerp/D");
    newTree->Branch("cosThetaPara", &ctPara, "cosThetaPara/D");
    newTree->Branch("cosThetaL", &ctL, "cosThetaL/D");
    newTree->Branch("cosThetaH", &ctH, "cosThetaH/D");
    newTree->Branch("phiL", &phiL, "phiL/D");
    newTree->Branch("phiH", &phiH, "phiH/D");
    newTree->Branch("dphi", &dphi, "dphi/D");

    ////////////////////////////////////////////////////////////////////////////
    //const
    newTree->Branch("cosThetaPerp_c", &cctPerp, "cosThetaPerp_c/D");
    newTree->Branch("cosThetaPara_c", &cctPara, "cosThetaPara_c/D");
    newTree->Branch("cosThetaL_c", &cctL, "cosThetaL_c/D");
    newTree->Branch("cosThetaH_c", &cctH, "cosThetaH_c/D");
    newTree->Branch("phiL_c", &cphiL, "phiL_c/D");
    newTree->Branch("phiH_c", &cphiH, "phiH_c/D");
    newTree->Branch("dphi_c", &cdphi, "dphi_c/D");

    if ( isMC )
    {
      ////////////////////////////////////////////////////////////////////////////
      //true
      newTree->Branch("true_cosThetaPerp", &ctPerpT, "true_cosThetaPerp/D");
      newTree->Branch("true_cosThetaPara", &ctParaT, "true_cosThetaPara/D");
      newTree->Branch("true_cosThetaL", &ctLT, "true_cosThetaL/D");
      newTree->Branch("true_cosThetaH", &ctHT, "true_cosThetaH/D");
      newTree->Branch("true_phiL", &phiLT, "true_phiL/D");
      newTree->Branch("true_phiH", &phiHT, "true_phiH/D");
      newTree->Branch("true_dphi", &dphiT, "true_dphi/D");
      newTree->Branch("true_q2", &true_q2, "true_q2/D" );
      newTree->Branch("true_mpk", &true_mpk, "true_mpk/D" );
      ////////////////////////////////////////////////////////////////////////////
      //true q2 ranges
      // newTree->Branch("true_q2_jpsi",  &true_q2_jpsi,  "true_q2_jpsi/O");
      // newTree->Branch("true_q2_psi2s", &true_q2_psi2s, "true_q2_psi2s/O");
      // newTree->Branch("true_q2_rare",  &true_q2_rare,  "true_q2_rare/O");
    }
  }

  newTree->Branch( "q2", &q2, "q2/D" );
  // Anja did this
  // newTree->Branch( "q2_diff", &q2_diff, "q2_diff/D" );
  // newTree->Branch( "l0", &l0, "l0/O" );
  // newTree->Branch( "hlt1", &hlt1, "hlt1/O" );
  // newTree->Branch( "hlt2", &hlt2, "hlt2/O" );
  // newTree->Branch( "trigger", &trigger, "trigger/O" );

  ////////////////////////////////////////////////////////////////////////////
  // //q2 ranges
  // newTree->Branch("q2_jpsi",  &q2_jpsi,  "q2_jpsi/O");
  // newTree->Branch("q2_psi2s", &q2_psi2s, "q2_psi2s/O");
  // newTree->Branch("q2_rare",  &q2_rare,  "q2_rare/O");

  TLorentzVector initialProton;
  double charge=0;
  initialProton.SetPxPyPzE(0,0,1.,0);

  long nevents = reader->GetEntries();
  if ( !config["events"].empty() )
  {
    nevents = config["events"].asInt64();
  }

  for ( long i = 0; i < nevents; ++i ) {

    if (i%1000 == 0) {std::cout<<"Processing event = "<<nevents-i<<std::endl;}

    reader->GetEntry( i );

    TLorentzVector LbP4, JpsiP4, LP4, pP4, piP4, mupP4, mumP4, q_diff_P4;
    // tools::getFourMomentum( partName, reader, LbP4 );
    // tools::getFourMomentum( "Jpsi", reader, JpsiP4 );
    // tools::getFourMomentum( lambda0, reader, LP4 );
    tools::getFourMomentum( proton, reader, pP4 );
    tools::getFourMomentum( pion, reader, piP4 );

    if (muon1=="mum_0")//TODO!!!! (reader->GetValue(muon1+"_ID")<0)
    {
      tools::getFourMomentum(muon1,  reader, mupP4);
      tools::getFourMomentum(muon2, reader, mumP4);
    }
    else
    {
      tools::getFourMomentum(muon2, reader, mupP4);
      tools::getFourMomentum(muon1,  reader, mumP4);
    }
    LbP4 = pP4 + piP4 + mupP4 + mumP4;
    JpsiP4 = mupP4 + mumP4;
    LP4 = pP4 + piP4;


    ////////////////////////////////////////////////////////////////////////////////////////////
    // Anja did this
    // if (B03D)
    // {
    //   tools::B0PsiRAngles(LbP4, JpsiP4, LP4, mupP4, mumP4, pP4, piP4, ctL,
    //                       ctH, dphi);
    // }
    // else
    // {
      charge = 1;//TODO!!! reader->GetValue(proton+"_ID");
      // Anja did this because FCC collisions are head-on, right?
      // if (BeamCrossing)
      // {
      //   std::map<std::string, TLorentzVector> bc =
      //       tools::BeamCrossing(year, polarity);
      //   TVector3 boostPCOM = bc["pcom"].BoostVector();
      //   initialProton.SetPxPyPzE(bc["b1"][0], bc["b1"][1], bc["b1"][2],
      //                            bc["b1"][3]);
      //   initialProton.Boost(boostPCOM);
      //   mumP4.Boost(boostPCOM);
      //   mupP4.Boost(boostPCOM);
      //   piP4.Boost(boostPCOM);
      //   pP4.Boost(boostPCOM);

      //   JpsiP4 = mumP4 + mupP4;
      //   LP4 = piP4 + pP4;
      //   LbP4 = JpsiP4 + LP4;
      // }
      tools::LbPsiRAngles(initialProton, LbP4, JpsiP4, LP4, mupP4, mumP4, pP4,
                          piP4, charge, ctPerp, ctPara, ctL, ctH, phiL, phiH, dphi);
    // }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //true
    if (isMC)
    {
      //tools::getFourMomentum(partName,  reader, LbP4,false);
      //tools::getFourMomentum("Jpsi",    reader, JpsiP4,false);
      //tools::getFourMomentum("L0",      reader, LP4,false);
      tools::getFourMomentum(proton,   reader, pP4, false);
      tools::getFourMomentum(pion,  reader, piP4, false);

      if (muon1=="mum_0")//TODO!!!(reader->GetValue(muon1+"_TRUEID")<0)
      {
        tools::getFourMomentum(muon1,  reader, mupP4, false);
        tools::getFourMomentum(muon2, reader, mumP4, false);
      }
      else
      {
        tools::getFourMomentum(muon2, reader, mupP4, false);
        tools::getFourMomentum(muon1,  reader, mumP4, false);

      }
      JpsiP4=mumP4+mupP4;
      LP4=piP4+pP4;
      LbP4=JpsiP4+LP4;


      true_mpk=LP4.M() / 1e3;

      true_q2 = JpsiP4.M()*JpsiP4.M() / 1e6;

      // if (B03D)
      // {
      //   tools::B0PsiRAngles(LbP4, JpsiP4, LP4, mupP4, mumP4, pP4, piP4, ctLT,
      //                       ctHT, dphiT);
      // }
      // else
      // {
        charge = 1;//TODO!!!reader->GetValue(proton+"_TRUEID");
        tools::LbPsiRAngles(initialProton, LbP4, JpsiP4, LP4, mupP4, mumP4, pP4,
                            piP4, charge, ctPerpT, ctParaT, ctLT, ctHT, phiLT, phiHT, dphiT);
      // }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //constrained
    if ( constrained && (isLb) )
    {
      std::string head = "Lb";
      if (!isLb) {
        head = "B0";
      }
      // tools::getFourMomentum_C(head, "Jpsi",     reader, JpsiP4);
      if (isLb)
      {
        // tools::getFourMomentum_C(head, lambda0, reader, LP4);
        tools::getFourMomentum_C(head, proton, reader, pP4);
        tools::getFourMomentum_C(head, pion, reader, piP4);
      }
      else
      {
        // tools::getFourMomentum_C(head, lambda0, reader, LP4); // TODO: should be KS
        tools::getFourMomentum_C(head, proton, reader, pP4); // TODO: should be kaon
        tools::getFourMomentum_C(head, pion, reader, piP4);
      }

      if (muon1=="mum_0")//TODO!!!! (reader->GetValue(muon1+"_ID")<0)
      {
        tools::getFourMomentum_C(head, muon1,  reader, mupP4);
        tools::getFourMomentum_C(head, muon2, reader, mumP4);
      }
      else
      {
        tools::getFourMomentum_C(head, muon2, reader, mupP4);
        tools::getFourMomentum_C(head, muon1,  reader, mumP4);
      }

      LP4 = pP4 + piP4;

      if ( lambda0 == "L0" ) {
        double eL0 = TMath::Sqrt( LP4.P()*LP4.P() + 1.115683 * 1.115683 );
        LP4.SetE(eL0);
      } else if ( lambda0 == "Ks" ) {
        double eL0 = TMath::Sqrt( LP4.P()*LP4.P() + 0.497611 * 0.497611 );
        LP4.SetE(eL0);
      }

      JpsiP4 = mupP4 + mumP4;
      double eJpsi = TMath::Sqrt( JpsiP4.P()*JpsiP4.P() + 3.096900 * 3.096900 );
      JpsiP4.SetE(eJpsi);

      // Anja did this because FCC collisions are head-on, right?
      // if(BeamCrossing)
      // {
      //   std::map<std::string,TLorentzVector> bc = tools::BeamCrossing(year,polarity);
      //   TVector3 boostPCOM=bc["pcom"].BoostVector();
      //   initialProton.SetPxPyPzE(bc["b1"][0],bc["b1"][1],bc["b1"][2],bc["b1"][3]);
      //   initialProton.Boost(boostPCOM);
      //   mumP4.Boost(boostPCOM);
      //   mupP4.Boost(boostPCOM);
      //   piP4.Boost(boostPCOM);
      //   pP4.Boost(boostPCOM);

      // }
      // JpsiP4 = mumP4 + mupP4;
      // LP4 = piP4 + pP4;
      LbP4=JpsiP4+LP4;

      charge = 1;//TODO!!!reader->GetValue(proton+"_ID");
      tools::LbPsiRAngles(initialProton, LbP4, JpsiP4, LP4, mupP4, mumP4, pP4,
                          piP4, charge, cctPerp, cctPara, cctL, cctH, cphiL, cphiH, cdphi);

    }

    // Anja did this
    // // Trigger
    // trigger = tools::passTrigger( reader, "Jpsi", "Jpsi", year, //true,
    //                               l0, hlt1, hlt2 );

    // q2
    // q2 = reader->GetValue( "Jpsi_M" );
    // q2 = q2 / 1e3;  // Go to GeV^2 for simplicity
    // q2 = q2 * q2;

    // Dimuon from difference of Lb and L0
    // q_diff_P4 = LbP4 - LP4;
    // q2_diff = 1e-6*q_diff_P4.M2();

  
    // // q2 ranges
    // if ((  8   < q2 ) && ( q2 < 11 )) { q2_jpsi  = true; } else { q2_jpsi  = false; }
    // if (( 12.5 < q2 ) && ( q2 < 15 )) { q2_psi2s = true; } else { q2_psi2s = false; }
    // if (( q2 < 0.98 ) || (( 1.1 < q2 ) && ( q2 < 8 )) || (( 11 < q2 ) && (true_q2 < 12.5)) || ( 15 < true_q2 )) {
    //   q2_rare = true;
    // } else {
    //   q2_rare = false;
    // }

    // if (isMC) {
    //   if ((  8   < true_q2 ) && ( true_q2 < 11 )) { true_q2_jpsi  = true; } else { true_q2_jpsi  = false; }
    //   if (( 12.5 < true_q2 ) && ( true_q2 < 15 )) { true_q2_psi2s = true; } else { true_q2_psi2s = false; }
    //   if (( true_q2 < 0.98 ) || (( 1.1 < true_q2 ) && ( true_q2 < 8 )) || (( 11 < true_q2 ) && (true_q2 < 12.5)) || ( 15 < true_q2 )) {
    //     true_q2_rare = true;
    //   } else {
    //     true_q2_rare = false;
    //   }
    // }

    newTree->Fill();
  }
  newTree->Write();
  out.Close();
}
