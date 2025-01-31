#ifndef LHCB_MICHALS_TOOLS
#define LHCB_MICHALS_TOOLS

#include "TLorentzVector.h"
#include "ReadTree.hh"
#include <string>

namespace tools
{

  // When calling this function, user should choose right daughter to use.
  // Function takes care not to modify incoming vectors
  void armanteros( TLorentzVector& posDaughter,
                   TLorentzVector& negDaughter,
                   double* qt,
                   double* alpha );

  std::map< std::string, TLorentzVector > BeamCrossing( int year, int pol );

  void getFourMomentum( std::string partName,
                        TreeReader* reader,
                        TLorentzVector& result,
                        bool rapidsim,
                        bool reconstructed = true );

  void getFourMomentum_C( std::string head, std::string partName,
                        TreeReader* reader,
                        TLorentzVector& result,
                        bool rapidsim );

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
                     double& cosThetaH,
                     double& phiL,
                     double& phiH,
                     double& dphi );

  // void B0PsiRAngles( TLorentzVector pB,
  //                    TLorentzVector pJpsi,
  //                    TLorentzVector pR,
  //                    TLorentzVector pmp,
  //                    TLorentzVector pmm,
  //                    TLorentzVector ppos,
  //                    TLorentzVector pneg,
  //                    double& cosThetaL,
  //                    double& cosThetaH,
  //                    double& dphi );


}  // namespace tools

#endif
