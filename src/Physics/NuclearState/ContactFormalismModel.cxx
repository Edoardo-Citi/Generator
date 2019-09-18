//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Brian Coopersmith, University of Rochester

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>
#include <TRandom1.h>
#include <fstream>
#include <TH2.h>
#include <TParticlePDG.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/ContactFormalismModel.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/ConfigIsotopeMapUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using namespace genie::utils::config;
using namespace std;

//____________________________________________________________________________
ContactFormalismModel::ContactFormalismModel() :
NuclearModelI("genie::ContactFormalismModel")
{

}
//____________________________________________________________________________
ContactFormalismModel::ContactFormalismModel(string config) :
NuclearModelI("genie::CorrelatedNuclearModelI", config)
{

}
//____________________________________________________________________________
// Cleans up memory from probability distribution maps
//____________________________________________________________________________
ContactFormalismModel::~ContactFormalismModel()
{
  map<string, TH1D*>::iterator iter = fProbDistroMap.begin();
  for( ; iter != fProbDistroMap.begin(); ++iter) {
    TH1D * hst = iter->second;
    if(hst) {
      delete hst;
      hst=0;
    }
  }
  fProbDistroMap.clear();
}
//High momentum recoil nucleon
void ContactFormalismModel::HighMomentumRecoilNucleon (double p1, double* P2, double* E1,double kf,  int type, const Target & target) const{
 
  //Variable, must be checked and/or updated
  int       nbins = 100;
  int       nbinse = 1000;
  int       nbinsp = 1000;
  double    sigma_P_CM = 0.143;   //In GeV
  double    m=PDGLibrary::Instance()->Find(2112)->Mass();;      //Mass of the nucleon
  double    mp=PDGLibrary::Instance()->Find(2212)->Mass();;      //Mass of the nucleon 
  double DB = 0.033;      //Binding energy, must be UPDATED
  int    A = target.A();
  double P_max = fPCutOff;
  
    // I open the file with the Wave Function Squared and save it int arrays
    ifstream myfile;
    myfile.open("/genie/app/users/eciti/lamp_v3/Generator/src/Physics/NuclearState/all.txt");
    double k[100], phi2_nn[100], phi2_pn0[100], phi2_pn1[100];
    int i=0;
    
    for(i=0;i<100;i++){
      myfile >> k[i];
      k[i] = k[i] * (TMath::H()) * TMath::C() * 1e6/ (TMath::Qe()) / 2/ TMath::Pi(); // Conversion
      myfile >> phi2_nn[i];
      myfile >> phi2_pn0[i];
      myfile >> phi2_pn1[i];
    }
    myfile.close();
    
    //Open the histogram and start working for the wavefunction
    
    TH1D Histo_p_nn("histo_p_nn","histo_p_nn",nbins,k[0]/2,k[nbins-1]+k[0]/2);
    TH1D Histo_p_pn0("histo_p_pn0","histo_p_pn0",nbins,k[0]/2,k[nbins-1]+k[0]/2);
    TH1D Histo_p_pn1("histo_p_pn1","histo_p_pn1",nbins,k[0]/2,k[nbins-1]+k[0]/2);
   
    for(i=0;i<100;i++){
      Histo_p_nn.Fill(k[i],phi2_nn[i]);
      Histo_p_pn0.Fill(k[i],phi2_pn0[i]);
      Histo_p_pn1.Fill(k[i],phi2_pn1[i]);
    }
    
    //  Normalization and drawing
    Histo_p_nn.Scale(1/Histo_p_nn.Integral(),"width");
    Histo_p_pn0.Scale(1/Histo_p_pn0.Integral(),"width");
    Histo_p_pn1.Scale(1/Histo_p_pn1.Integral(),"width");
    
    //Here comes the hard work
    
    
    TH2D Histo_SF_nn("Spectral Function nn","Spectral Function nn",nbinsp,0,P_max,nbinse,(2 * m - sqrt(p1*p1 + m*m) -DB - (p1 +p1 )*(p1 + p1)/(2*m*(A-2))),m);
    TH2D Histo_SF_pn0("Spectral Function pn0","Spectral Function pn0",nbinsp,0,P_max,nbinse,(2 * m - sqrt(p1*p1 + m*m) -DB - (p1 +p1 )*(p1 + p1)/(2*m*(A-2))),m);
    TH2D Histo_SF_pn1("Spectral Function pn1","Spectral Function pn1",nbinsp,0,P_max,nbinse,(2 * m - sqrt(p1*p1 + m*m) -DB - (p1 +p1 )*(p1 + p1)/(2*m*(A-2))),m);
    
    
    
    //Definitions of variable and functions
    double dp =(P_max)/(double)nbinsp; //* (TMath::H()) * TMath::C() * 1e6/ (TMath::Qe()) / 2/ TMath::Pi();
    double de1 = m/(double)nbinse;
    double p2=0, e1 =0;
    double Weight_nn = 0, Weight_pn0 = 0, Weight_pn1 = 0;
    double Norm_SF_nn =0, Norm_SF_pn0=0, Norm_SF_pn1 =0;
    double P_R2 =0, P_CM2 =0, Cos_Theta_num =0, P_CM2_p=0;
    
    int p2_i=0, e1_start=0, e1_stop =0 ;
    
    for(p2_i=0; p2_i<nbinsp; p2_i++){
      p2 = dp * ((double)(p2_i + 0.5));
       
      //Kinematic bounds on energy, handle with care and be sure to UPDATE if you change the physics
        
      e1_start = (2 * m - sqrt(p2*p2 + m*m) -DB - (p1 +p2 )*(p1 + p2)/(2*m*(A-2)))/de1;
      if (e1_start<0){e1_start=0;}
        
      e1_stop =  (2 * m  - sqrt(p2*p2 + m*m) - DB) / de1;
      if (e1_stop > nbinse){e1_stop = nbinse;}
        
      for(i=e1_start; i<e1_stop; i++){
	e1 = de1 * ((double)(i + 0.5));
	P_CM2 = 2 * m * (A-2) * (2 * m - e1 - sqrt(p2*p2 + m*m) - DB);
	P_CM2_p = 2 * mp * (A-2) * (2 * mp - e1 - sqrt(p2*p2 + mp*mp) - DB);
	P_R2 = (p1*p1 +p2 *p2  - P_CM2/2)/2;
	Cos_Theta_num = (P_CM2 - p1*p1 -p2*p2);
               
	//Condizioni di non validitá
               
	if (P_CM2 <=0 || Cos_Theta_num < -(2*p1*p2) || Cos_Theta_num  > 2*p1*p2 || P_R2<=(kf*kf) ||P_R2>(k[99]*k[99])){
	  Weight_nn = 0;        //Not in domain
	  Weight_pn0 = 0;
	  Weight_pn1 = 0;
                            
	}
	else {
	  //Here lays the Physics, if you change the equations UPDATE here

	  Weight_nn = p2/p1 * mp  * (Histo_p_nn.GetBinContent((int)((sqrt(P_R2)/k[0]+0.5)))) / (  sigma_P_CM* sigma_P_CM* sigma_P_CM) * TMath::Exp((-1)*P_CM2_p/(2* sigma_P_CM*sigma_P_CM) );
	  Weight_pn0 = p2/p1 * m  * (Histo_p_pn0.GetBinContent((int)((sqrt(P_R2)/k[0]+0.5)))) / (  sigma_P_CM* sigma_P_CM* sigma_P_CM) * TMath::Exp((-1)*P_CM2/(2* sigma_P_CM*sigma_P_CM) );
                            
	  Weight_pn1 = p2/p1 * m  * (Histo_p_pn1.GetBinContent((int)((sqrt(P_R2)/k[0]+0.5)))) / (  sigma_P_CM* sigma_P_CM* sigma_P_CM) * TMath::Exp((-1)*P_CM2/(2* sigma_P_CM*sigma_P_CM) );

	}
	                    
	Histo_SF_nn.Fill (p2, e1, Weight_nn);
	Histo_SF_pn0.Fill (p2, e1, Weight_pn0);
	Histo_SF_pn1.Fill (p2, e1, Weight_pn1);
	  
	Norm_SF_nn += Weight_nn;
	Norm_SF_pn0+= Weight_pn0;
	Norm_SF_pn1+=Weight_pn1;
      }
    }
    
    //Normalizations
   
     Histo_SF_nn.Scale (1/Norm_SF_nn,"width");
     Histo_SF_pn0.Scale (1/Norm_SF_pn0,"width");
     Histo_SF_pn1.Scale (1/Norm_SF_pn1,"width");

    if (type == 1){Histo_SF_nn.GetRandom2(P2[0],E1[0]);}
    if (type == 2){Histo_SF_pn0.GetRandom2(P2[0],E1[0]);}
    if (type == 3){Histo_SF_pn1.GetRandom2(P2[0],E1[0]);}


}
//____________________________________________________________________________
// Set the removal energy, 3 momentum, and FermiMover interaction type
//____________________________________________________________________________
bool ContactFormalismModel::GenerateNucleon(const Target & target) const
{
  assert(target.HitNucIsSet());
  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);


  int target_pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  int nucleon_pdgc = target.HitNucPdg();
  //-- import Fermi momentum                                                                                                                                                         
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
  double  KF = kft->FindClosestKF(target_pdgc,nucleon_pdgc);

  //-- set fermi momentum vector
  //

  if ( target.A() > 1 ) {
    TH1D * prob = this->ProbDistro(target);
    if(!prob) {
      LOG("EffectiveSF", pNOTICE)
              << "Null nucleon momentum probability distribution";
      exit(1);
    }

    double p = prob->GetRandom();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("EffectiveSF", pDEBUG) << "|p,nucleon| = " << p;
#endif

    RandomGen * rnd = RandomGen::Instance();

    double costheta = -1. + 2. * rnd->RndGen().Rndm();
    double sintheta = TMath::Sqrt(1.-costheta*costheta);
    double fi       = 2 * kPi * rnd->RndGen().Rndm();
    double cosfi    = TMath::Cos(fi);
    double sinfi    = TMath::Sin(fi);

    double px = p*sintheta*cosfi;
    double py = p*sintheta*sinfi;
    double pz = p*costheta;

    fCurrMomentum.SetXYZ(px, py, pz);

  

  //Using Contact Formalism to predict SF, a recoil nucleon is generated

    if(HasSecondNucleon(p,target)==true){
    double p2 =0 ; //Momentum of the recoil nucleon
    double e1 =0 ; //Removal energy
    double coin = 0;//Random variable

    //Costants from CF, tells the program from which Isospin/Spin state extract te p2 distribution, must be UPDATED
    double Cpn0 = 0.;
    double Cpn1 = 0.25;
    double Cpp = 0.50;
    double m=0;
    TRandom1 * rx = new TRandom1();
    coin = rnd -> RndGen().Rndm();


    //The extract p2    
    if(Cpn0<= coin && coin < Cpn1){
      HighMomentumRecoilNucleon(p,&p2,&e1,KF,2,target);
      m=PDGLibrary::Instance()->Find(2112)->Mass();
      SCRPartnerPdgc[0] = 2112;
    }
    if(Cpn1<= coin && coin < Cpp){
      HighMomentumRecoilNucleon(p,&p2,&e1,KF,3,target);
      m=PDGLibrary::Instance()->Find(2112)->Mass();
      SCRPartnerPdgc[0] = 2112;
    }
    if(Cpp<= coin && coin <= 1){
      HighMomentumRecoilNucleon(p,&p2,&e1,KF,1,target); 
      m=PDGLibrary::Instance()->Find(2212)->Mass();
      SCRPartnerPdgc[0] = 2212;
    }
    
   SCRPartnerMomentum[0] = p2;

    //Now I compute the direction of p2

    //Some value for the calculations
    double DB = 0.33;  //Binding energy difference, mist be UPDATED
    double P_CM2 = 2 * m * (target.A()-2) * (2 * m - e1 - sqrt(p2*p2 + m*m) - DB);
    double Cos_Theta_2 = (P_CM2 - p*p -p2*p2)/(2*p*p2);
    double Sin_Theta_2 = TMath::Sqrt(1.-Cos_Theta_2);
    double fi2         = 2 * kPi * rnd->RndGen().Rndm();

    //now I compute p2 before the p1 rotation
    double p2zp1z = p2 * Cos_Theta_2;
    double p2xp1z = p2 * Sin_Theta_2 * TMath::Cos(fi2);
    double p2yp1z = p2 * Sin_Theta_2 * TMath::Sin(fi2);
 
    //And now I rotate p2
    double p2x = p2xp1z*costheta*cosfi - p2yp1z*sinfi + p2zp1z*sintheta*cosfi;
    double p2y = p2xp1z*sintheta*sinfi + p2yp1z*cosfi + p2zp1z*sintheta*sinfi;
    double p2z = -p2xp1z*sinfi + p2zp1z*costheta;

    SCRPartner3Momentum->SetXYZ(p2x,p2y,p2z);
  }

  }
  //
  //-- set removal energy
  //

  fCurrRemovalEnergy = this->ReturnBindingEnergy(target);
  double f1p1h = this->Returnf1p1h(target);
  // Since TE increases the QE peak via a 2p2h process, we decrease f1p1h
  // in order to increase the 2p2h interaction to account for this enhancement.
  f1p1h /= this->GetTransEnh1p1hMod(target);
  if ( RandomGen::Instance() -> RndGen().Rndm() < f1p1h) {
    fFermiMoverInteractionType = kFermiMoveEffectiveSF1p1h;
  } else if (fEjectSecondNucleon2p2h) {
    fFermiMoverInteractionType = kFermiMoveEffectiveSF2p2h_eject;
  } else {
    fFermiMoverInteractionType = kFermiMoveEffectiveSF2p2h_noeject;
  }

  return true;
}
//____________________________________________________________________________
// Returns the probability of the bin with given momentum. I don't know what w
// is supposed to be, but I copied its implementation from Bodek-Ritchie.
// Implements the interface.
//____________________________________________________________________________
double ContactFormalismModel::Prob(double mom, double w, const Target & target) const
{
  if(w < 0) {
     TH1D * prob_distr = this->ProbDistro(target);
     int bin = prob_distr->FindBin(mom);
     double y  = prob_distr->GetBinContent(bin);
     double dx = prob_distr->GetBinWidth(bin);
     double prob  = y * dx;
     return prob;
  }
  return 1;
}
//____________________________________________________________________________
// Check the map of nucleons to see if we have a probability distribution to
// compute with.  If not, make one.
//____________________________________________________________________________
TH1D * ContactFormalismModel::ProbDistro(const Target & target) const
{
  //-- return stored /if already computed/
  map<string, TH1D*>::iterator it = fProbDistroMap.find(target.AsString());
  if(it != fProbDistroMap.end()) return it->second;

  LOG("EffectiveSF", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString();
  LOG("EffectiveSF", pNOTICE)
               << "P(cut-off) = " << fPCutOff << ", P(max) = " << fPMax;

  //-- get information for the nuclear target
  int nucleon_pdgc = target.HitNucPdg();
  assert( pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc) );
  return this->MakeEffectiveSF(target);

}
//____________________________________________________________________________
// If transverse enhancement form factor modification is enabled, we must
// increase the 2p2h contribution to account for the QE peak enhancement.
// This gets that factor based on the target.
//____________________________________________________________________________
double ContactFormalismModel::GetTransEnh1p1hMod(const Target& target) const {
  double transEnhMod;
	if(GetValueFromNuclearMaps(target, fTransEnh1p1hMods,
	                           fRangeTransEnh1p1hMods, &transEnhMod) &&
	   transEnhMod > 0) {
	  return transEnhMod;
	}
  // If none specified, assume no enhancement to quasielastic peak
  return 1.0;
}
//____________________________________________________________________________
// Makes a momentum distribtuion for the given target using parameters
// from the config file.
//____________________________________________________________________________
TH1D * ContactFormalismModel::MakeEffectiveSF(const Target & target) const
{
  // First check for individually specified nuclei
  int pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  map<int,vector<double> >::const_iterator it = fProbDistParams.find(pdgc);
  if(it != fProbDistParams.end()) {
    vector<double> v = it->second;
    return this->MakeEffectiveSF(v[0], v[1], v[2], v[3],
                                 v[4], v[5], v[6], target);
  }

  // Then check in the ranges of A
  map<pair<int, int>, vector<double> >::const_iterator range_it = fRangeProbDistParams.begin();
  for(; range_it != fRangeProbDistParams.end(); ++range_it) {
    if (target.A() >= range_it->first.first && target.A() <= range_it->first.second) {
      vector<double> v = range_it->second;
      return this->MakeEffectiveSF(v[0], v[1], v[2], v[3],
                                   v[4], v[5], v[6], target);
    }
  }

  return NULL;
}
//____________________________________________________________________________
//Generates the tail of the momentum distribution
//____________________________________________________________________________
void ContactFormalismModel::FillMomentumTail(TH1D* Histo_to_fill, double Kstart, double Kend, int nbins)const {

    // I open the file with the Wave Function Squared and save it into arrays
  
    ifstream myfile;
    myfile.open("/genie/app/users/eciti/lamp_v3/Generator/src/Physics/NuclearState/all.txt");
    double k[100], phi2_nn[100], phi2_pn0[100], phi2_pn1[100];
     int i=0;
    
    for(i=0;i<100;i++){
      myfile >> k[i];
      k[i] = k[i] * (TMath::H()) * TMath::C() * 1e6/ (TMath::Qe()) / 2/ TMath::Pi(); // Conversion in GeV/C from fm-1
      myfile >> phi2_nn[i];
      myfile >> phi2_pn0[i];
      myfile >> phi2_pn1[i];
    }
    myfile.close();

    //Fill the Histogram for relative momentum distribution
    TH1D* histo_p_nn = new TH1D("histo_p_nn","histo_p_nn",101,k[1]/2,k[99]+k[1]/2);
    TH1D* histo_p_pn0 = new TH1D("histo_p_pn0","histo_p_pn0",101,k[1]/2,k[99]+k[1]/2);
     TH1D* histo_p_pn1 = new TH1D("histo_p_pn1","histo_p_pn1",101,k[1]/2,k[99]+k[1]/2);

    for(i=0;i<100;i++){
      histo_p_nn->Fill(k[i],phi2_nn[i]);
      histo_p_pn0->Fill(k[i],phi2_pn0[i]);
      histo_p_pn1->Fill(k[i],phi2_pn1[i]);
    } 

    //Some variables """"WARNING MUST BE UPDATED"""
    double    sigma_P_CM = 0.143;   //In Gev
    double    m=PDGLibrary::Instance()->Find(2112)->Mass();;      //Mass of the nucleon                                                                                                
    double    mp=PDGLibrary::Instance()->Find(2212)->Mass();;      //Mass of the nucleon  
    // int    A = target.A();
    double dp =(Kend-Kstart)/(double)nbins; 
    int nbinsp = 1e3, p1_i=0;
    double dPr = (k[99]-Kstart)/nbinsp, Pr=0;
    double p1=0;
    double DB = 0.033;      //Binding energy
    double Weight_pp = 0, Weight_pn0 = 0, Weight_pn1 = 0;
    double Cost_min = -1, Cost_min_p=-1, Cost_max = 1;

    TH1D* Histo_p1 = new TH1D ("P1","P1",nbins,Kstart,Kend);

    //Contact Formalism Constants """"WARNING MUST BE UPDATED""""
    double Cpn0 = 1;
    double Cpn1 =1;
    double Cpp = 1;
  
    for(p1_i=0; p1_i<nbins; p1_i++){
      p1 = dp * ((double)(p1_i + 0.5))+Kstart;
      Cost_min = ((2*m - DB)*(2*m - DB) - m*m - p1*p1 - Pr*Pr)/(2*p1*Pr); //Kinematic cut """CHECK IT IF YOU CHANGE THE PHYSICS"""
      Cost_min_p = ((2*mp - DB)*(2*mp - DB) - mp*mp - p1*p1 - Pr*Pr)/(2*p1*Pr); //Kinematic cut """CHECK IT IF YOU CHANGE THE PHYSICS""" 
      if (Cost_min<-1){Cost_min=-1;}
      if (Cost_min_p<-1){Cost_min_p=-1;}
      for(i=0; i<nbinsp; i++){
	Pr= dPr * ((double)(i + 0.5)) + Kstart;
  
	//Condizioni di non validitá
	if (Pr<Kstart||Pr>k[99]||Cost_min>Cost_max){
	  Weight_pp = 0;        //Not in domain
	  Weight_pn0 = 0;
	  Weight_pn1 = 0;
                            
	}
	else {
	  // Here lays the Physics, if you change the equation UPDATE here
	  if(Cost_min_p>Cost_max){Weight_pp=0;}
	  else{ Weight_pp =Cpp*2* Pr/p1  * (histo_p_nn -> Interpolate(Pr) )/ sigma_P_CM *(TMath::Exp((4*p1*Pr*Cost_max-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM))-TMath::Exp((4*p1*Pr*Cost_min_p-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM)));}
                        
	  Weight_pn0 =Cpn0* Pr/p1  * (histo_p_pn0 -> Interpolate(Pr) )/ sigma_P_CM *(TMath::Exp((4*p1*Pr*Cost_max-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM))-TMath::Exp((4*p1*Pr*Cost_min-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM)));
                            
	  Weight_pn1 =Cpn1* Pr/p1  * (histo_p_pn1 -> Interpolate(Pr) )/ sigma_P_CM *(TMath::Exp((4*p1*Pr*Cost_max-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM))-TMath::Exp((4*p1*Pr*Cost_min-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM)));

	}
	Histo_p1 ->Fill (p1, Weight_pp);
	  Histo_p1 ->Fill (p1, Weight_pn1);
	Histo_p1 ->Fill (p1, Weight_pn0);
                       
      } 
    }
    
    Histo_p1 ->Scale(0.2/Histo_p1->Integral(),"width");

		      //Now I Fill real Histogram
		      for(i=0;i<nbins;i++){
			Histo_to_fill->Fill(Histo_p1->GetBinCenter(i),Histo_p1->GetBinContent(i));
		      }
		      return;
  }
//____________________________________________________________________________
// Makes a momentum distribution using the factors below (see reference) and
// inserts it into the nucleus/momentum distribution map.
//____________________________________________________________________________
TH1D * ContactFormalismModel::MakeEffectiveSF(double bs, double bp, double alpha,
                                    double beta, double c1, double c2,
                                    double c3, const Target & target) const
{
  int target_pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  int nucleon_pdgc = target.HitNucPdg();

  //-- import Fermi momentum
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
  double  KF = kft->FindClosestKF(target_pdgc,nucleon_pdgc);
 



  //-- create the probability distribution
  int npbins = (int) (1000 * fPMax);

  TH1D * prob = new TH1D("", "", npbins, 0, fPMax);
  prob->SetDirectory(0);

  double dp = fPMax / (npbins-1);

  //Distribution below Fermi Momentum, usual Effective SF
  int CriticalIndex = KF/dp;
  for(int i = 0; i < CriticalIndex; i++) {
    double p  = i * dp;
    double y = p / 0.197;
    double as = c1 * exp(-pow(bs*y,2));
    double ap = c2 * pow(bp * y, 2) * exp(-pow(bp * y, 2));
    double at = c3 * pow(y, beta) * exp(-alpha * (y - 2));
    double rr = (3.14159265 / 4) * (as + ap + at) * pow(y, 2) / 0.197;
    double dP_dp = rr / 1.01691371;
    if(p>fPCutOff)
      dP_dp = 0;
    assert(dP_dp >= 0);
    // calculate probability density : dProbability/dp
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("EffectiveSF", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
#endif
    prob->Fill(p, dP_dp);
 }

  //-- normalize the probability distribution up to this point, !!!WARNING!!! for now using 0.8 normalization for first part, must put the correct value                                                                                                                                        
  prob->Scale( 0.8 / prob->Integral("width") );
 
 
  //Distribution behind Fermi Momentum, Contact Formalism
  int EndIndex = fPCutOff/dp;
  FillMomentumTail(prob , (CriticalIndex)*dp,(EndIndex)*dp , EndIndex-CriticalIndex);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("EffectiveSF", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
#endif

    //    prob->Scale( 1 / prob->Integral("width") );


  //-- store
  fProbDistroMap.insert(
      map<string, TH1D*>::value_type(target.AsString(),prob));
  return prob;
}
//---------------------------------------------------------------------------
//Check the emitted nucleon was in a pair or not
//--------------------------------------------------------------------------
bool ContactFormalismModel::HasSecondNucleon (double p1,const Target & target) const
{
  int target_pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  int nucleon_pdgc = target.HitNucPdg();

  //-- import Fermi momentum                                                                                                                                               
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
  double  KF = kft->FindClosestKF(target_pdgc,nucleon_pdgc);

  if(p1>=KF) return true;

  else return false;

}
//---------------------------------------------------------------------------
//Return removal evenrgy of the recoil nucleon
//---------------------------------------------------------------------------
double ContactFormalismModel::SecondRemovalEnergy(int nucleon_pdgc) const
{  
  double    m=PDGLibrary::Instance()->Find(nucleon_pdgc)->Mass();
  return sqrt(SCRPartnerMomentum[0] * SCRPartnerMomentum[0] +m*m);

}
//---------------------------------------------------------------------------
//Return Second momentum
//--------------------------------------------------------------------------
double ContactFormalismModel::SecondMomentum(void) const
{
  return SCRPartnerMomentum[0];
}
//--------------------------------------------------------------------------
//Return srecoil nucleon species
//--------------------------------------------------------------------------
int ContactFormalismModel::Secondpdgc(void) const
{
  return SCRPartnerPdgc[0];
}
//--------------------------------------------------------------------------
//Return 3 Vector for the recoil nucleon momentum
//-------------------------------------------------------------------------
const TVector3& ContactFormalismModel::SecondMomentum3(void) const
{
  return SCRPartner3Momentum[0];
}
//____________________________________________________________________________
// Returns the binding energy for a given nucleus.
//____________________________________________________________________________
double ContactFormalismModel::ReturnBindingEnergy(const Target & target) const
{
  double binding_en;
  if (GetValueFromNuclearMaps(target, fNucRmvE, fRangeNucRmvE, &binding_en) &&
      binding_en > 0) {
    return binding_en;
  }
  return 0;
}
//____________________________________________________________________________
// Returns the fraction of 1p1h events for a given nucleus.  All other events
// are 2p2h.
//____________________________________________________________________________
double ContactFormalismModel::Returnf1p1h(const Target & target) const
{
  double f1p1h;
  if (GetValueFromNuclearMaps(target, f1p1hMap, fRange1p1hMap, &f1p1h) &&
      f1p1h >= 0 && f1p1h <= 1) {
    return f1p1h;
  }
  return 1;
}
//____________________________________________________________________________
void ContactFormalismModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ContactFormalismModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
// Every parameter for this comes from the config files.
//____________________________________________________________________________
void ContactFormalismModel::LoadConfig(void)
{
  // Load Fermi Momentum Table 
  this->GetParam("FermiMomentumTable",fKFTable);




  this->GetParamDef("EjectSecondNucleon2p2h", fEjectSecondNucleon2p2h, false);

  this->GetParamDef("MomentumMax",    fPMax,    1.0);
  this->GetParamDef("MomentumCutOff", fPCutOff, 0.65);
  assert(fPMax > 0 && fPCutOff > 0 && fPCutOff <= fPMax);

  // Find out if Transverse enhancement is enabled to figure out whether to load
  // the 2p2h enhancement parameters.
  this->GetParam("UseElFFTransverseEnhancement", fUseElFFTransEnh );
  if (!fUseElFFTransEnh) {
    LOG("EffectiveSF", pINFO)
        << "Transverse enhancement not used; "
        << "Do not increase the 2p2h cross section.";
  }
  else {
    LoadAllIsotopesForKey("TransEnhf1p1hMod", "EffectiveSF",
                          GetOwnedConfig(), &fTransEnh1p1hMods);
    LoadAllNucARangesForKey("TransEnhf1p1hMod", "EffectiveSF",
                            GetOwnedConfig(), &fRangeTransEnh1p1hMods);
  }

  LoadAllIsotopesForKey("BindingEnergy", "EffectiveSF", GetOwnedConfig(), &fNucRmvE);
  LoadAllNucARangesForKey("BindingEnergy", "EffectiveSF",
                          GetOwnedConfig(), &fRangeNucRmvE);
  LoadAllIsotopesForKey("f1p1h", "EffectiveSF", GetOwnedConfig(), &f1p1hMap);
  LoadAllNucARangesForKey("f1p1h", "EffectiveSF", GetOwnedConfig(), &fRange1p1hMap);

  for (int Z = 1; Z < 140; Z++) {
    for (int A = Z; A < 3*Z; A++) {
      const int pdgc = pdg::IonPdgCode(A, Z);
      double bs, bp, alpha, beta, c1, c2, c3;
      if (GetDoubleKeyPDG("bs", pdgc, GetOwnedConfig(), &bs) &&
          GetDoubleKeyPDG("bp", pdgc, GetOwnedConfig(), &bp) &&
          GetDoubleKeyPDG("alpha", pdgc, GetOwnedConfig(), &alpha) &&
          GetDoubleKeyPDG("beta", pdgc, GetOwnedConfig(), &beta) &&
          GetDoubleKeyPDG("c1", pdgc, GetOwnedConfig(), &c1) &&
          GetDoubleKeyPDG("c2", pdgc, GetOwnedConfig(), &c2) &&
          GetDoubleKeyPDG("c3", pdgc, GetOwnedConfig(), &c3)) {
        vector<double> pars = vector<double>();
        pars.push_back(bs);
        pars.push_back(bp);
        pars.push_back(alpha);
        pars.push_back(beta);
        pars.push_back(c1);
        pars.push_back(c2);
        pars.push_back(c3);
        LOG("EffectiveSF", pINFO)
          << "Nucleus: " << pdgc << " -> using bs =  " << bs << "; bp = "<< bp
          << "; alpha = " << alpha << "; beta = "<<beta<<"; c1 = "<<c1
          <<"; c2 = "<<c2<< "; c3 = " << c3;
        fProbDistParams[pdgc] = pars;
      }
    }
  }
  for(int lowA = 1; lowA < 3 * 140; lowA++) {
    for(int highA = lowA; highA < 3 * 140; highA++) {
      double bs, bp, alpha, beta, c1, c2, c3;
      if (GetDoubleKeyRangeNucA("bs",   lowA, highA, GetOwnedConfig(), &bs)    &&
          GetDoubleKeyRangeNucA("bp",   lowA, highA, GetOwnedConfig(), &bp)    &&
          GetDoubleKeyRangeNucA("alpha",lowA, highA, GetOwnedConfig(), &alpha) &&
          GetDoubleKeyRangeNucA("beta", lowA, highA, GetOwnedConfig(), &beta) &&
          GetDoubleKeyRangeNucA("c1",   lowA, highA, GetOwnedConfig(), &c1) &&
          GetDoubleKeyRangeNucA("c2",   lowA, highA, GetOwnedConfig(), &c2) &&
          GetDoubleKeyRangeNucA("c3",   lowA, highA, GetOwnedConfig(), &c3)) {
        vector<double> pars = vector<double>();
        pars.push_back(bs);
        pars.push_back(bp);
        pars.push_back(alpha);
        pars.push_back(beta);
        pars.push_back(c1);
        pars.push_back(c2);
        pars.push_back(c3);
        LOG("EffectiveSF", pINFO) << "For "<< lowA - 1 <<" < A < " << highA + 1
          <<" -> using bs =  " << bs << "; bp = "<< bp
          << "; alpha = " << alpha << "; beta = "<<beta<<"; c1 = "<<c1
          <<"; c2 = "<<c2<< "; c3 = " << c3;
        fRangeProbDistParams[pair<int, int>(lowA, highA)] = pars;
      }
    }
  }
}
//____________________________________________________________________________
