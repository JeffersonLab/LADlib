#ifndef __LADSimEvent_h
#define __LADSimEvent_h

//#include "gmn_dig_tree.h"
//#include "g4LAD_tree.h"

#include "lad_tree_digitized.h"

enum Exp_t    { kLAD}; //Functionality to change experiment type (and thus the type of tree used.)

class TTree;

class LADSimEvent {
 public:
  //LADSimEvent(){};                 // Default constructor, for ROOT I/O
  //Now we want to initialize the ROOT tree the same way 
  LADSimEvent(TTree* tree, Exp_t experiment=kLAD);//, std::vector<TString> det_list);
  
  virtual ~LADSimEvent(){};
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;

  ULong64_t RunID;
  ULong64_t EvtID;

  //AJRP: redesign for simplicity/portability/maintainability:

  //We probably don't actually need these getter and setter functions, but it doesn't hurt to have them.
  //void SetExperiment( TString expname ){ fExperiment = expname; }
  //TString GetExperiment() const { return fExperiment; }
  
  //TString fExperiment; //string that tells us which simulated experiment we are doing
  
  void SetExperiment( Exp_t exp ){ fExperiment = exp; }
  Exp_t GetExperiment() const { return fExperiment; }
  
  Exp_t fExperiment;
  
  //Auto-generated ROOT Tree classes for each experiment ROOT tree; generated using TTree::MakeClass()
  
  //Later on, any time we want to analyze a g4LAD root file whose format has changed, we can just run TTree::MakeClass on that root file with the appropriate
  //class name, and copy the source and header files into LAD-offline,
  //recompile, and voila: compatibility guaranteed:
  lad_tree_digitized *Tlad;
  
  
  ClassDef(LADSimEvent, 1) // Simulated data for one event
};

#endif
