#include "tools.hh"

namespace tools {

  void getFourMomentum( std::string particleName,
                        TreeReader* reader,
                        TLorentzVector& result,
                        bool reconstructed ) {
    double px, py, pz, pe;
    if ( reconstructed ) {
      px = reader->GetValue( ( particleName + "_PX" ).c_str() );
      py = reader->GetValue( ( particleName + "_PY" ).c_str() );
      pz = reader->GetValue( ( particleName + "_PZ" ).c_str() );
      pe = reader->GetValue( ( particleName + "_E" ).c_str() );
    } else {
      px = reader->GetValue( ( particleName + "_PX_TRUE" ).c_str() );
      py = reader->GetValue( ( particleName + "_PY_TRUE" ).c_str() );
      pz = reader->GetValue( ( particleName + "_PZ_TRUE" ).c_str() );
      pe = reader->GetValue( ( particleName + "_E_TRUE" ).c_str() );
    }
    result.SetPxPyPzE( px, py, pz, pe );

    return;
  }

  void getFourMomentum_C( std::string head,
                          std::string particleName,
                          TreeReader* reader,
                          TLorentzVector& result ) {
    double px, py, pz, pe;
    double mass = 0.;

    if ( particleName == "L0" ) {
      mass = 1.115683;
    }
    if ( particleName == "Ks" ) {
      mass = 0.497611;
    } else if ( particleName == "mup" || particleName == "mum" ) {
      mass = 0.1056583745;
    } else if ( head == "Lb" &&
                ( particleName == "pp_0" ) ) {
      mass = 0.938272081;
    } else if ( head == "B0" &&
                ( particleName == "pp_0" ) ) {
      mass = 0.13957061;
    } else if ( particleName == "pim" || particleName == "pim" ) {
      mass = 0.13957061;
    } else if ( particleName == "Jpsi" ) {
      mass = 3.096900;
    } else {
      mass = -10e-36;
    }
    //    std::cout << head << " " << particleName << " " << mass << std::endl;

    px = reader->GetValue(
        ( particleName + "_PX" ).c_str() );
    py = reader->GetValue(
        ( particleName + "_PY" ).c_str() );
    pz = reader->GetValue(
        ( particleName + "_PZ" ).c_str() );
    pe = TMath::Sqrt( mass * mass + px * px + py * py + pz * pz );
    //    pe=reader->GetValue( ("Lb_E_"+particleName+"_DauCon").c_str() );
    result.SetPxPyPzE( px, py, pz, pe );
    return;
  }

  void LbPsiRAngles( TLorentzVector initialProton,
                     TLorentzVector pB,
                     TLorentzVector pJpsi,
                     TLorentzVector pLambda,
                     TLorentzVector pmp,
                     TLorentzVector pmm,
                     TLorentzVector pp,
                     TLorentzVector ppi,
                     int pCharge,
                     double& cosThetaPerp,
                     double& cosThetaPara,
                     double& cosThetaL,
                     double& cosThetaH ,
                     double& phiL,
                     double& phiB,
                     double& dphi ) {
    bool decayLambda = true;
    if ( pp.M() <
         1e-5 )  // Should be zero exactly, but allow for numerical treatment
    {
      decayLambda = false;
    }
    TVector3 analyzer;
    analyzer = initialProton.Vect().Cross( pB.Vect() );
    analyzer *= 1.0 / analyzer.Mag();

    // Move everything to Lambda_b rest frame
    if ( pB.Vect().Mag() > 1e-6 ) {
      TVector3 boost = -1.0 * pB.BoostVector();
      pB.Boost( boost );
      pJpsi.Boost( boost );
      pLambda.Boost( boost );
      pmp.Boost( boost );
      pmm.Boost( boost );
      if ( decayLambda ) {
        pp.Boost( boost );
        ppi.Boost( boost );
      }
    }
    // First level
    cosThetaPerp =
        analyzer.Dot( pLambda.Vect() ) / analyzer.Mag() / pLambda.Vect().Mag();
    cosThetaPara = pB.Vect().Dot( pLambda.Vect() ) / pB.Vect().Mag() / pLambda.Vect().Mag();

    // Frame 1
    TVector3 z1 = pLambda.Vect() * ( 1. / pLambda.Vect().Mag() );
    TVector3 y1 = analyzer.Cross( z1 );
    y1          = y1 * ( 1. / y1.Mag() );
    TVector3 x1 = z1.Cross( y1 );
    x1          = x1 * ( 1. / x1.Mag() );
    // Proton in Lambda rest frame
    TLorentzVector ppLRest = pp;
    TVector3 h1;
    if ( decayLambda ) {
      ppLRest.Boost( -1.0 * pLambda.BoostVector() );
      cosThetaH  = z1.Dot( ppLRest.Vect() ) / ppLRest.Vect().Mag();
      TVector3 pperp =
          ppLRest.Vect() - 1.0 * ppLRest.Vect().Mag() * z1 * cosThetaH ;
      double cosPhi = pperp.Dot( x1 ) / pperp.Mag();
      double sinPhi = pperp.Dot( y1 ) / pperp.Mag();
      phiB          = TMath::ACos( cosPhi );
      if ( sinPhi < 0 )
        phiB = -phiB;

      h1 = ( ( ppLRest.Vect().Cross( pLambda.Vect() ) ) );
    }
    // Frame 2
    TVector3 z2 = pJpsi.Vect() * ( 1. / pJpsi.Vect().Mag() );
    TVector3 y2 = analyzer.Cross( z2 );
    y2          = y2 * ( 1. / y2.Mag() );
    TVector3 x2 = z2.Cross( y2 );
    x2          = x2 * ( 1. / x2.Mag() );
    // Proton in Lambda rest frame
    ppLRest = pmp;
    ppLRest.Boost( -1.0 * pJpsi.BoostVector() );
    cosThetaL = z2.Dot( ppLRest.Vect() ) / ppLRest.Vect().Mag();
    TVector3 pperp =
        ppLRest.Vect() - 1.0 * ppLRest.Vect().Mag() * z2 * cosThetaL;
    double cosPhi = pperp.Dot( x2 ) / pperp.Mag();
    double sinPhi = pperp.Dot( y2 ) / pperp.Mag();
    phiL          = TMath::ACos( cosPhi );
    if ( sinPhi < 0 )
      phiL = -phiL;

    if ( decayLambda ) {
      TVector3 h2( ( ppLRest.Vect().Cross( pJpsi.Vect() ) ) );

      // Calculate delta phi
      //  TVector3 h1((pmp.Vect()).Cross(pmm.Vect()));
      //  TVector3 h2((pp.Vect()).Cross(ppi.Vect()));
      TVector3 dir( h1.Cross( h2 ) );
      float a1 = dir.Dot( pJpsi.Vect() ) / dir.Mag() / pJpsi.Vect().Mag();
      a1       = TMath::ACos( a1 );
      dphi     = h1.Dot( h2 ) / h1.Mag() / h2.Mag();
      //  infoToStore.dphi=TMath::Pi()+TMath::ACos(infoToStore.dphi);
      dphi = TMath::ACos( dphi );
      if ( a1 > TMath::Pi() - 0.001 ) {
        dphi = -dphi;
      }
    }

    dphi = -999;
    if ( decayLambda ) {
      //
      TVector3 p3KPi    = pLambda.Vect();
      TVector3 p3K      = pp.Vect();
      TVector3 p3Psi    = pJpsi.Vect();
      TVector3 p3MuPlus = pmp.Vect();

      TVector3 aK = p3K - p3KPi * ( p3K.Dot( p3KPi ) / p3KPi.Mag2() );
      TVector3 aMuPlus =
          p3MuPlus - p3Psi * ( p3MuPlus.Dot( p3Psi ) / p3Psi.Mag2() );

      // angle between K* and Psi decay planes in B0 rest frame

      TVector3 p3Pi      = ppi.Vect();
      TVector3 p3MuMinus = pmm.Vect();

      TVector3 n_L = p3K.Cross( p3Pi ) * ( 1. / p3K.Cross( p3Pi ).Mag() );
      TVector3 n_J = p3MuMinus.Cross( p3MuPlus ) *
                     ( 1. / p3MuMinus.Cross( p3MuPlus ).Mag() );

      double sinT = ( n_L.Cross( n_J ) ).Dot( p3KPi ) / p3KPi.Mag();
      double cosT = n_L.Dot( n_J );
      dphi        = std::atan2( sinT, cosT );


      if ( std::isnan( dphi ) ) {
        // std::cout << "phi is nan" << std::endl;
        dphi = 0.;
      }
    }

    // Now we handle antiparticles (charge of the proton being negative)

    if ( pCharge < 0 ) {
      cosThetaL *= -1;
      phiB = -1 * phiB;
      if ( phiL > 0 ) {
        phiL = TMath::Pi() - phiL;
      } else {
        phiL = -TMath::Pi() - phiL;
      }
      if ( dphi > -999 ) {
        if ( dphi > 0 ) {
          dphi = TMath::Pi() - dphi;
        } else {
          dphi = -TMath::Pi() - dphi;
        }
      }
    }
  }

  // double coshel0_Tomasz( TLorentzVector particle,
  //                        TLorentzVector parent,
  //                        TLorentzVector grandparent ) {
  //   TVector3 boosttoparent = -( parent.BoostVector() );

  //   particle.Boost( boosttoparent );
  //   grandparent.Boost( boosttoparent );

  //   TVector3 particle3    = particle.Vect();
  //   TVector3 grandparent3 = -grandparent.Vect();

  //   Double_t numerator   = particle3.Dot( grandparent3 );
  //   Double_t denominator = ( particle3.Mag() ) * ( grandparent3.Mag() );
  //   Double_t temp        = numerator / denominator;

  //   return temp;
  // }

  // double planeAngle0_Tomasz( TLorentzVector particleFrame,
  //                            TLorentzVector particleA,
  //                            TLorentzVector particleB,
  //                            TLorentzVector particleC,
  //                            TLorentzVector particleD ) {
  //   TVector3 boostToFrame = -( particleFrame.BoostVector() );

  //   TLorentzVector boostedAxis = particleC + particleD;
  //   boostedAxis.Boost( boostToFrame );
  //   TLorentzVector boostedA = particleA;
  //   boostedA.Boost( boostToFrame );
  //   TLorentzVector boostedB = particleB;
  //   boostedB.Boost( boostToFrame );
  //   TLorentzVector boostedC = particleC;
  //   boostedC.Boost( boostToFrame );
  //   TLorentzVector boostedD = particleD;
  //   boostedD.Boost( boostToFrame );

  //   TVector3 vecA = ( boostedA.Vect() ).Unit();
  //   TVector3 vecB = ( boostedB.Vect() ).Unit();
  //   TVector3 vecC = ( boostedC.Vect() ).Unit();
  //   TVector3 vecD = ( boostedD.Vect() ).Unit();

  //   TVector3 el = ( vecA.Cross( vecB ) ).Unit();
  //   // TVector3 ek = ( vecC.Cross( vecD ) ).Unit() ;
  //   TVector3 ek = ( vecD.Cross( vecC ) ).Unit();

  //   TVector3 ez = ( boostedAxis.Vect() ).Unit();

  //   double cosPhi = ( ek.Dot( el ) );
  //   double sinPhi = ( el.Cross( ek ) ).Dot( ez );
  //   // cout <<"cosPhi  " << cosPhi <<endl;
  //   double phi = acos( cosPhi );

  //   return ( sinPhi > 0.0 ? phi : -phi );
  // }

}  // namespace tools
