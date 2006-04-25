//____________________________________________________________________________
/*!

\class    genie::DISInteractionListGenerator

\brief    Concrete implementations of the InteractionListGeneratorI interface.
          Generate a list of all the Interaction (= event summary) objects that
          can be generated by the DIS EventGenerator.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 13, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _DIS_INTERACTION_LIST_GENERATOR_H_
#define _DIS_INTERACTION_LIST_GENERATOR_H_

#include "EVGCore/InteractionListGeneratorI.h"

namespace genie {

class DISInteractionListGenerator : public InteractionListGeneratorI {

public :

  DISInteractionListGenerator();
  DISInteractionListGenerator(string config);
  ~DISInteractionListGenerator();

  //-- implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfigData(void);

  bool fIsCC;
  bool fIsNC;
  bool fIsCharm;
};

}      // genie namespace

#endif // _DIS_INTERACTION_LIST_GENERATOR_H_
