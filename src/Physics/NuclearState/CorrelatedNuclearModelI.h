//____________________________________________________________________________
/*!

\class    genie::CorrelatedNuclearModelI

\brief    Pure abstract base class.
          Defines extensions to the NuclearModelI interface that allow
          information about a correlated second nucleon to be retrieved
          when appropriate.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  September 16, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef _CORRELATED_NUCLEAR_MODEL_I_H_
#define _CORRELATED_NUCLEAR_MODEL_I_H_

#include "Physics/NuclearState/NuclearModelI.h"

namespace genie {

class CorrelatedNuclearModelI : public NuclearModelI {

public:

  virtual ~CorrelatedNuclearModelI() {};

  virtual double SecondMomentum(void) const = 0;
  virtual const TVector3& SecondMomentum3(void) const = 0;
  virtual double SecondRemovalEnergy(void) const = 0;
  virtual bool HasSecondNucleon(void) const = 0;

protected:

  CorrelatedNuclearModelI() : NuclearModelI() {}
  CorrelatedNuclearModelI(std::string name) : NuclearModelI(name) {}
  CorrelatedNuclearModelI(std::string name, std::string config)
    : NuclearModelI(name, config) {}

};

} // genie namespace
#endif // _CORRELATED_NUCLEAR_MODEL_I_H_
