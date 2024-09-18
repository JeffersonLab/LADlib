#include "THcLADGEM.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THcGlobals.h"
#include "THaCutList.h"
#include "THcParmList.h"
#include "VarDef.h"
#include "VarType.h"
#include "THaApparatus.h"
#include "THcHallCSpectrometer.h"

using namespace std;

//____________________________________________________________
THcLADGEM::THcLADGEM( const char* name, const char* description,
		      THaApparatus* apparatus) :
  THaTrackingDetector(name, description, apparatus)
{
  // constructor

  fModules.clear();
  fPlanes.clear();
  fModulesInitialized = false;

  fNModules = 0;
  fNPlanes = 0;

}

//____________________________________________________________
THcLADGEM::~THcLADGEM()
{
  // Destructor
  for(auto module : fModules ) {
    delete module;
  }
  fModules.clear();

  for(auto plane : fPlanes ) {
    delete plane;
  }
  fPlanes.clear();

}

//____________________________________________________________
inline
void THcLADGEM::ClearEvent()
{

  // Clear per-event data
  fNhits = 0;

}

//____________________________________________________________
THaAnalysisObject::EStatus THcLADGEM::Init( const TDatime& date)
{

  cout << "THcLADGEM::Init" << endl;

  EStatus status;
  if( (status = THaTrackingDetector::Init( date )) )
    return fStatus = status;

  fPresentP = nullptr;
  THaVar* vpresent = gHaVars->Find(Form("%s.present", GetApparatus()->GetName()));
  if(vpresent) {
    fPresentP = (Bool_t*) vpresent->GetValuePointer();
  }

  // Call Hall C style DetectorMap 
  // GEM channel maps to be defined in MAPS/LAD/detector.map 
  // e.g. gHcDetectorMap->FillMap(fDetMap, EngineID)
  // InitHitList(fDetMap, "THcLADGEMHit", fDetMap->GetTotNumChan()+1,0, RefTimeCut)

  // Create subdetectors
  Setup(GetName(), GetTitle());

  for(auto& module: fModules ) {
    cout << "Calling THcGEMModule->Init" << endl;
    status = module->Init(date);
    if( status != kOK )
      return fStatus = status;
  }

  return fStatus = kOK;

}

//____________________________________________________________
Int_t THcLADGEM::DefineVariables( EMode mode )
{

  // Initialize/register global variables
  /*

    RVarDef vars[] = {
    {"", "", ""},
    {"", "", ""},
    {nullptr}
    };

    return DefineVarsFromList( vars, mode );
   */

  return 0;

}

//____________________________________________________________________________________
void THcLADGEM::Setup(const char* name, const char* description)
{
  // initial values
  fNModules = 1;
  fNPlanes = 1;

  // FIXME: hard-coded for now, SBS GEM m0 - m33
  //  fNModules = 34;

  char prefix[2];
  
  prefix[0] = std::tolower(GetApparatus()->GetName()[0]);
  prefix[1] = '\0';

  DBRequest list[] = {
    {"gem_num_modules", &fNModules, kInt},
    {"gem_num_planes",  &fNPlanes,  kInt},
    {0}
  };
  gHcParms->LoadParmValues((DBRequest*)&list, prefix);

  for(int imod = 0; imod < fNModules; imod++) {
    cout << "Initialize GEM modules " << endl;
    THcLADGEMModule* new_module = new THcLADGEMModule(Form("m%d",imod), Form("m%d",imod), this);
    fModules.push_back(new_module);
  }
  
}


//____________________________________________________________
Int_t THcLADGEM::ReadDatabase( const TDatime& date )
{

  cout << "THcLADGEM::ReadDatabase" << endl;
  // Called by THaDetectorBase::Init()
  // Read parameters from THcParmList

  return kOK;

}

//____________________________________________________________
Int_t THcLADGEM::Decode( const THaEvData& evdata )
{

  ClearEvent();

  Int_t Nmpdfound = 0;

  //  auto *rdata = (UInt_t*) evdata.GetRawDataBuffer();
  //  AnalyzeSBSDataBuffer(rdata,  );

  // SBS GEM data decoding
  /*
  for( UInt_t i=0; i < fDetMap->GetSize(); i++) {
    
    THaDetMap::Module* d = fDetMap->GetModule(i);
    Decoder::MPDModule* impd = dynamic_cast<Decoder::MPDModule*>(evdata.GetModule(d->crate, d->slot));
    if(impd) {
      // found MPD
      Nmpdfound++;
    }      
  }
  */

  // Decode MPD data
  for(auto& module: fModules) {
    module->Decode(evdata);
  }

  Bool_t present = kTRUE;
  if(fPresentP) {
    present = *fPresentP;
  }

  return 0;
}

//____________________________________________________________
Int_t THcLADGEM::CoarseTrack( TClonesArray& tracks )
{
  // track finding called here and added them into the TClonesArray tracks
  // In hcana, tracks are defined in the detector coordinate system

  // Track finding algorithm for SBS GEM is defined in 
  // SBSGEMTrackerBase::find_tracks
  // which first needs to reconstruct 2D cluster hit by calling SBSGEMModule::find_2Dhits
  // We separate cluster hit finding to SBSGEMModule::CoarseProcess 
  // and will add track finding related functions here

  return 0;
}

//____________________________________________________________
Int_t THcLADGEM::FineTrack( TClonesArray& tracks )
{
  return 0;
}

ClassImp(THcLADGEM)

