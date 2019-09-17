//____________________________________________________________________________
/*!

\class    genie::ContactFormalismModel

\brief    Contact formalism approach. TODO: write this
          Implements the CorrelatedNuclearModelI interface.

\ref

\author   Edoardo Citi
          Steven Gardiner

\created  September 16, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CONTACT_FORMALISM_MODEL_H_
#define _CONTACT_FORMALISM_MODEL_H_

#include <map>

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class ContactFormalismModel : public NuclearModelI {

public:
  ContactFormalismModel();
  ContactFormalismModel(string config);
  virtual ~ContactFormalismModel();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;
  double         Prob            (double mom, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target &) const
  {
    return kNucmEffSpectralFunc;
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  // Implement the extra functions in the CorrelatedNuclearModelI interface
double SecondMomentum(void) const;
 // const TVector3& SecondMomentum3(void) const;
 double SecondRemovalEnergy(double p2, int nucleon_pdgc) const;
 bool HasSecondNucleon (double p1,const Target & target) const;

private:
  TH1D * ProbDistro (const Target & t) const;

  TH1D * MakeEffectiveSF(const Target & target) const;

  TH1D * MakeEffectiveSF(double bs, double bp, double alpha, double beta,
                         double c1, double c2, double c3,
                         const Target & target) const;

  double ReturnBindingEnergy(const Target & target) const;
  double GetTransEnh1p1hMod(const Target& target) const;
  void HighMomentumRecoilNucleon(double p1, double* P2, double* E1,double kf,  int type, const Target & target) const;
  void FillMomentumTail(TH1D* Histo_to_fill, double Kstart, double Kend, int nbins) const;

  double Returnf1p1h(const Target & target) const;
  void   LoadConfig (void);

  mutable map<string, TH1D *> fProbDistroMap;
  double fPMax;
  double fPCutOff;
  bool   fEjectSecondNucleon2p2h;
  bool   fUseElFFTransEnh;
  string fKFTable;

  // Map from PDG code to spectral function parameters
  map<int, double> fNucRmvE;
  map<int, double> f1p1hMap;
  map<int, std::vector<double> > fProbDistParams;
  map<int, double> fTransEnh1p1hMods;

  // Map from range of A (pair<lowA, highA> inclusive> to spectral
  // function parameters.
  map<pair<int, int>, double> fRangeNucRmvE;
  map<pair<int, int>, double> fRange1p1hMap;
  map<pair<int, int>, std::vector<double> > fRangeProbDistParams;
  map<pair<int, int>, double> fRangeTransEnh1p1hMods;
};

}         // genie namespace
#endif    // _EFFECTIVE_SF_H_
