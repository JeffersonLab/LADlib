#include "THcSignalHit.h"
#include "THcLADSpectrometer.h"
#include "THcHitList.h"
#include "THcLADHodoscope.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THaCutList.h"
#include "VarDef.h"
#include "VarType.h"
#include "THcGlobals.h"
#include "THcParmList.h"
#include "THcDetectorMap.h"

//_________________________________________________________________
THcLADHodoscope::THcLADHodoscope( const char* name, const char* description,
				  THaApparatus* apparatus) :
  THaNonTrackingDetector(name, description, apparatus)
{
  // Constructor
  fNPlanes = 1;

}

//_________________________________________________________________
THcLADHodoscope::~THcLADHodoscope()
{
  // Destructor

  for(int ip = 0; ip < fNPlanes; ip++)
    delete fPlanes[ip];
  delete [] fPlanes;

  delete [] fNPaddle;      fNPaddle      = NULL;
  delete [] fTdcOffset;    fTdcOffset    = NULL;
  delete [] fAdcTdcOffset; fAdcTdcOffset = NULL;
  delete [] fHodoSlop;     fHodoSlop     = NULL;

  delete [] fHodoNegAdcTimeWindowMin; fHodoNegAdcTimeWindowMin = NULL;
  delete [] fHodoNegAdcTimeWindowMax; fHodoNegAdcTimeWindowMax = NULL;
  delete [] fHodoPosAdcTimeWindowMin; fHodoPosAdcTimeWindowMin = NULL;
  delete [] fHodoPosAdcTimeWindowMax; fHodoPosAdcTimeWindowMax = NULL;

}

//_________________________________________________________________
void THcLADHodoscope::Clear( Option_t* opt )
{

}
//_________________________________________________________________

void THcLADHodoscope::Setup(const char* name, const char* description)
{
  /**
     Create the scintillator plane objects for the hodoscope.

     Uses the Xhodo_num_planes and Xhodo_plane_names to get the number of
     planes and their names.

  */
  if( IsZombie()) return;

  // fDebug = 1;  // Keep this at one while we're working on the code

  char prefix[2];

  prefix[0]=tolower(GetApparatus()->GetName()[0]);
  prefix[1]='\0';

  TString temp(prefix[0]);
  
  TString histname=temp+"_timehist";
  
  string planenamelist;
  DBRequest listextra[]={
    {"hodo_num_planes", &fNPlanes, kInt},
    {"hodo_plane_names",&planenamelist, kString},
    {"hodo_tdcrefcut", &fTDC_RefTimeCut, kInt, 0, 1},
    {"hodo_adcrefcut", &fADC_RefTimeCut, kInt, 0, 1},
    {0}
  };
  //fNPlanes = 4; 		// Default if not defined
  fTDC_RefTimeCut = 0;		// Minimum allowed reference times
  fADC_RefTimeCut = 0;
  gHcParms->LoadParmValues((DBRequest*)&listextra,prefix);

  cout << "Plane Name List : " << planenamelist << endl;

  vector<string> plane_names = Podd::vsplit(planenamelist);
  // Plane names
  if(plane_names.size() != (UInt_t) fNPlanes) {
    cout << "ERROR: Number of planes " << fNPlanes << " doesn't agree with number of plane names " << plane_names.size() << endl;
    // Should quit.  Is there an official way to quit?
  }

  fPlaneNames = new char* [fNPlanes];
  for(Int_t i=0;i<fNPlanes;i++) {
    fPlaneNames[i] = new char[plane_names[i].length()+1];
    strcpy(fPlaneNames[i], plane_names[i].c_str());
  }

  // Probably shouldn't assume that description is defined
  char* desc = new char[strlen(description)+100];
  fPlanes = new THcLADHodoPlane* [fNPlanes];
  for(Int_t i=0;i < fNPlanes;i++) {
    strcpy(desc, description);
    strcat(desc, " Plane ");
    strcat(desc, fPlaneNames[i]);
    fPlanes[i] = new THcLADHodoPlane(fPlaneNames[i], desc, i+1, this); // Number planes starting from zero!!
    cout << "Created Scintillator Plane " << fPlaneNames[i] << ", " << desc << endl;
  }
  // // Save the nominal particle mass
  // THcHallCSpectrometer *app = dynamic_cast<THcHallCSpectrometer*>(GetApparatus());
  // fPartMass = app->GetParticleMass();
  // fBetaNominal = app->GetBetaAtPcentral();


  delete [] desc;
}

//_________________________________________________________________
THaAnalysisObject::EStatus THcLADHodoscope::Init( const TDatime& date )
{
  Setup(GetName(), GetTitle());

  char EngineDID[] = "xSCIN";
  EngineDID[0] = toupper(GetApparatus()->GetName()[0]);
  if( gHcDetectorMap->FillMap(fDetMap, EngineDID) < 0 ) {
    static const char* const here = "Init()";
    Error( Here(here), "Error filling detectormap for %s.", EngineDID );
    return kInitError;
  }

  // Should probably put this in ReadDatabase as we will know the
  // maximum number of hits after setting up the detector map
  // But it needs to happen before the sub detectors are initialized
  // so that they can get the pointer to the hitlist.
  cout << " Hodo tdc ref time cut = " << fTDC_RefTimeCut << " " << fADC_RefTimeCut << endl;

  InitHitList(fDetMap, "THcRawHodoHit", fDetMap->GetTotNumChan()+1,
	      fTDC_RefTimeCut, fADC_RefTimeCut);

  EStatus status;
  if( (status = THaNonTrackingDetector::Init( date )) )
    return fStatus = status;
  
  for(Int_t ip=0;ip<fNPlanes;ip++) {
    if((status = fPlanes[ip]->Init( date ))) {
      return fStatus=status;
    }
  }


  // fNScinHits     = new Int_t [fNPlanes];
  // fGoodPlaneTime = new Bool_t [fNPlanes];
  // fNPlaneTime    = new Int_t [fNPlanes];
  // fSumPlaneTime  = new Double_t [fNPlanes];

  //  Double_t  fHitCnt4 = 0., fHitCnt3 = 0.;

  // Int_t m = 0;
  // fScinHit = new Double_t*[fNPlanes];
  // for ( m = 0; m < fNPlanes; m++ ){
  //   fScinHit[m] = new Double_t[fNPaddle[0]];
  // }

  // for (int ip=0; ip<fNPlanes; ++ip) {
  //   fScinHitPaddle.emplace_back(fNPaddle[ip], 0);
  // }

  fPresentP = 0;
  THaVar* vpresent = gHaVars->Find(Form("%s.present",GetApparatus()->GetName()));
  if(vpresent) {
    fPresentP = (Bool_t *) vpresent->GetValuePointer();
  }
  
  return kOK;
}

//_________________________________________________________________
Int_t THcLADHodoscope::End(THaRunBase* run)
{
  // Do we really need this function?

  return 0;
}

//_________________________________________________________________
Int_t THcLADHodoscope::DefineVariables( EMode mode )
{

  return 0;

}

//_________________________________________________________________
Int_t THcLADHodoscope::ReadDatabase( const TDatime& date )
{

  cout << "THcLADHodoscope::ReadDatabase()" << endl;
  char prefix[2];
  prefix[0] = tolower(GetApparatus()->GetName()[0]); // "lad"
  prefix[1] = '\0';

  // First, get # of planes
  DBRequest list[] = {
    {"hodo_num_planes", &fNPlanes, kInt},
    {0}
  };

  gHcParms->LoadParmValues((DBRequest*)&list, prefix);

  // since we define each hodoscope as a separate detector
  // we will need to use the detector name to load parameters
  // for each detector -- to be updated
  fNPaddle = new Int_t [fNPlanes];
  for(int ip = 0; ip < fNPlanes; ip++)
    {
      DBRequest list2[] = {
	{Form("hodo_%d_nr", ip),      &fNPaddle[ip], kInt},
	{0}
      };

      gHcParms->LoadParmValues((DBRequest*)&list2, prefix);
    }

  // for all planes
  fTdcOffset    = new Int_t[fNPlanes];
  fAdcTdcOffset = new Double_t[fNPlanes];
  fHodoSlop     = new Double_t[fNPlanes];
  for(int ip = 0; ip < fNPlanes; ip++)
    {
      fTdcOffset[ip] = 0;
      fAdcTdcOffset[ip] = 0;
      fHodoSlop[ip] = 0;
    }

  // for all elements
  fMaxHodoScin = fNPaddle[0] * fNPlanes;
  fHodoPosAdcTimeWindowMin = new Double_t [fMaxHodoScin];
  fHodoPosAdcTimeWindowMax = new Double_t [fMaxHodoScin];
  fHodoNegAdcTimeWindowMin = new Double_t [fMaxHodoScin];
  fHodoNegAdcTimeWindowMax = new Double_t [fMaxHodoScin];

  for(int ii=0; ii < fMaxHodoScin; ii++)
    {
      fHodoPosAdcTimeWindowMin[ii] = -1000.;
      fHodoPosAdcTimeWindowMax[ii] = 1000.;
      fHodoNegAdcTimeWindowMin[ii] = -1000.;
      fHodoNegAdcTimeWindowMax[ii] = 1000.;
    }      

  DBRequest list3[] = {
    {"cosmicflag",                       &fCosmicFlag,            kInt,            0,  1},
    {"hodo_tdc_offset",                  fTdcOffset,              kInt,     (UInt_t) fNPlanes, 1},
    {"hodo_adc_tdc_offset",              fAdcTdcOffset,           kDouble,  (UInt_t) fNPlanes, 1},
    {"hodo_PosAdcTimeWindowMin",         fHodoPosAdcTimeWindowMin, kDouble,  (UInt_t) fMaxHodoScin, 1},
    {"hodo_PosAdcTimeWindowMax",         fHodoPosAdcTimeWindowMax, kDouble,  (UInt_t) fMaxHodoScin, 1},
    {"hodo_NegAdcTimeWindowMin",         fHodoNegAdcTimeWindowMin, kDouble,  (UInt_t) fMaxHodoScin, 1},
    {"hodo_NegAdcTimeWindowMax",         fHodoNegAdcTimeWindowMax, kDouble,  (UInt_t) fMaxHodoScin, 1},
    {0}
  };

  fCosmicFlag = 0;
  fScinTdcMin = 0;
  fScinTdcMax = 0;
  fScinTdcToTime = 0;

  gHcParms->LoadParmValues((DBRequest*)&list3, prefix);

  return kOK;

}

//_________________________________________________________________
Int_t THcLADHodoscope::Decode( const THaEvData& evdata )
{

  // Decode raw data and pass it to hitlist
  // Read raw data -- THcHitList::DecodeToHitList
  // Processing hitlist for each plane -- THcLADHodoPlane::ProcessHits

  Bool_t present = kTRUE;
  if(fPresentP) {
    present = *fPresentP;
  }

  fNHits = DecodeToHitList(evdata, !present);

  // To analyze pedestal events -- Must define "Pedestal_event" cut in the cuts .def file
  // do we want to do this? or calculate pedestal for each event (using the first # of samples, e.g)
  // keeping it for now
  if(gHaCuts->Result("Pedestal_event")) {
    Int_t nexthit = 0;
    for(Int_t ip = 0; ip < fNPlanes; ip++) {
      nexthit = fPlanes[ip]->AccumulatePedestals(fRawHitList, nexthit);
    }
    fAnalyzePedestals = 1;    // Analyze pedestals first normal events
    return (0);
  }	
  if(fAnalyzePedestals) {
    for(Int_t ip =0; ip < fNPlanes; ip++) {
      fPlanes[ip]->CalculatePedestals();
    }
    fAnalyzePedestals = 0;    // Don't analyze pedestals next event
  }


  Int_t nexthit = 0;
  for(Int_t iplane = 0; iplane < fNPlanes; iplane++) {

    nexthit = fPlanes[iplane]->ProcessHits(fRawHitList, nexthit);

  }
    
  return fNHits;
}

//_________________________________________________________________
Int_t THcLADHodoscope::CoarseProcess( TClonesArray& tracks )
{

  // Loop over all tracks and get corrected time


  return 0;

}

//_________________________________________________________________
Int_t THcLADHodoscope::FineProcess( TClonesArray& tracks )
{

  return 0;
}

ClassImp(THcLADHodoscope)

