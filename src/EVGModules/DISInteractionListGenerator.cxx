//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGModules/DISInteractionListGenerator.h"
#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
DISInteractionListGenerator::DISInteractionListGenerator() :
InteractionListGeneratorI("genie::DISInteractionListGenerator")
{

}
//___________________________________________________________________________
DISInteractionListGenerator::DISInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::DISInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DISInteractionListGenerator::~DISInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DISInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("InteractionList", pINFO)
                          << "InitialState = " << init_state.AsString();


  InteractionType_t inttype;
  if      (fIsCC) inttype = kIntWeakCC;
  else if (fIsNC) inttype = kIntWeakNC;
  else {
     LOG("InteractionList", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  int      nupdg  = init_state.GetProbePDGCode();
  Target * target = init_state.GetTargetPtr();

  if( !pdg::IsNeutrino(nupdg) && !pdg::IsAntiNeutrino(nupdg) ) {
     LOG("InteractionList", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  bool hasP = (target->Z() > 0);
  bool hasN = (target->N() > 0);

  InteractionList * intlist = new InteractionList;

  int nuclpdg[2] = { kPdgProton, kPdgNeutron };

  for(int inucl=0; inucl<2; inucl++) {

    int struck_nucleon = nuclpdg[inucl];

    if( (struck_nucleon == kPdgProton  && hasP) ||
        (struck_nucleon == kPdgNeutron && hasN) ) {

      ProcessInfo   proc_info(kScDeepInelastic, inttype);

      Interaction * interaction = new Interaction(init_state, proc_info);
      Target * target = interaction->GetInitialStatePtr()->GetTargetPtr();

      target->SetStruckNucleonPDGCode(struck_nucleon);

      if(fIsCharm) {
         XclsTag exclusive_tag;
         exclusive_tag.SetCharm();
         interaction->SetExclusiveTag(exclusive_tag);
      }
      intlist->push_back(interaction);
    }
  }

  if(intlist->size() == 0) {
     LOG("InteractionList", pERROR)
         << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete intlist;
     return 0;
  }
  return intlist;
}
//___________________________________________________________________________
void DISInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DISInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DISInteractionListGenerator::LoadConfigData(void)
{
  fIsCC    = fConfig->GetBoolDef("is-CC",    false);
  fIsNC    = fConfig->GetBoolDef("is-NC",    false);
  fIsCharm = fConfig->GetBoolDef("is-Charm", false);
}
//____________________________________________________________________________

