//____________________________________________________________________________
/*!

\class    genie::GHepRecord

\brief    GENIE's GHEP MC event record.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 1, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GHEP_RECORD_H_
#define _GHEP_RECORD_H_

#include <ostream>
#include <vector>

#include <TClonesArray.h>
#include <TBits.h>

#include "Interaction/Interaction.h" 
#include "GHEP/GHepStatus.h"

class TLorentzVector;

using std::ostream;
using std::vector;

namespace genie {

class GHepParticle;

class GHepRecord : public TClonesArray {

public :

  GHepRecord();
  GHepRecord(int size);
  GHepRecord(const GHepRecord & record);
  virtual ~GHepRecord();

  //-- methods to attach / get summary information

  virtual Interaction * GetInteraction    (void) const;
  virtual void          AttachInteraction (Interaction * interaction);

  //-- common record operations

  virtual void Copy                    (const GHepRecord & record);
  virtual void ResetRecord             (void);
  virtual void CompactifyDaughterLists (void);
  virtual void Clear                   (Option_t * opt="");

  //-- provide a simplified wrapper of the 'new with placement'
  //   TClonesArray object insertion method

  //   ALWAYS use these methods to insert new particles as they check
  //   for the compactness of the daughter lists.
  //   Note that the record might be automatically re-arranged as the
  //   result of your GHepParticle insertion

  virtual void AddParticle (const GHepParticle & p);
  virtual void AddParticle (int pdg, GHepStatus_t ist,
                     int mom1, int mom2, int dau1, int dau2,
                        const TLorentzVector & p, const TLorentzVector & v);
  virtual void AddParticle (int pdg, GHepStatus_t ist,
                     int mom1, int mom2, int dau1, int dau2,
                           double px, double py, double pz, double E,
                                    double x, double y, double z, double t);

  //-- methods to search the GHEP (STDHEP-like) record

  virtual GHepParticle * Particle     (int position) const;
  virtual GHepParticle * FindParticle (int pdg, GHepStatus_t ist, int start) const;

  virtual int ParticlePosition (int pdg, GHepStatus_t i, int start=0) const;
  virtual int ParticlePosition (GHepParticle * particle, int start=0) const;

  virtual vector<int> * GetStableDescendants(int position) const;

  //-- easy access methods for the most frequently used GHEP entries

  virtual GHepParticle * Probe                            (void) const;
  virtual GHepParticle * TargetNucleus                    (void) const;
  virtual GHepParticle * RemnantNucleus                   (void) const;
  virtual GHepParticle * StruckNucleon                    (void) const;
  virtual GHepParticle * StruckElectron                   (void) const;
  virtual GHepParticle * FinalStatePrimaryLepton          (void) const;
  virtual GHepParticle * FinalStateHadronicSystem         (void) const;
  virtual int            ProbePosition                    (void) const;
  virtual int            TargetNucleusPosition            (void) const;
  virtual int            RemnantNucleusPosition           (void) const;
  virtual int            StruckNucleonPosition            (void) const;
  virtual int            StruckElectronPosition           (void) const;
  virtual int            FinalStatePrimaryLeptonPosition  (void) const;
  virtual int            FinalStateHadronicSystemPosition (void) const; 

  //-- number of GHepParticle occurences in GHEP

  virtual unsigned int NEntries (int pdg, GHepStatus_t ist, int start=0) const;
  virtual unsigned int NEntries (int pdg, int start=0) const;

  //-- methods to switch on/off and ask for event record flags

  virtual TBits * EventFlags   (void) const { return fEventFlags; }
  virtual bool    IsUnphysical (void) const { return (fEventFlags->CountBits()>0); }

  //-- methods to set/get the event weight and cross sections

  virtual double GetWeight   (void) const  { return fWeight;   }
  virtual double GetXSec     (void) const  { return fXSec;     }
  virtual double GetDiffXSec (void) const  { return fDiffXSec; }
  virtual void   SetWeight   (double wght) { fWeight   = (wght>0) ? wght : 0.; }
  virtual void   SetXSec     (double xsec) { fXSec     = (xsec>0) ? xsec : 0.; }
  virtual void   SetDiffXSec (double xsec) { fDiffXSec = (xsec>0) ? xsec : 0.; }

  //-- set/get event vertex in detector coordinate system

  virtual TLorentzVector * Vertex (void) const { return fVtx; }

  virtual void SetVertex (double x, double y, double z, double t);
  virtual void SetVertex (const TLorentzVector & vtx);

  //-- methods & operators to print the record

  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const GHepRecord & event);

protected:

  // Summary information for the Initial State, Process Type & Kinematics
  Interaction * fInteraction;

  // Vertex position (in the detector coordinate system)
  TLorentzVector * fVtx;

  // Flags for the generated event
  TBits * fEventFlags;    

  // Misc info associated with the generated event
  double fWeight;         ///< event weight
  double fXSec;           ///< cross section for selected event
  double fDiffXSec;       ///< differential cross section for selected event kinematics

  // Utility methods
  void InitRecord  (void);
  void CleanRecord (void);

  // Methods used by the daughter list compactifier
  virtual void UpdateDaughterLists    (void);
  virtual bool HasCompactDaughterList (int pos);
  virtual void SwapParticles          (int i, int j);
  virtual void FinalizeDaughterLists  (void);
  virtual int  FirstNonInitStateEntry (void);

private:

ClassDef(GHepRecord, 1)

};

}      // genie namespace

#endif // _GHEP_RECORD_H_
