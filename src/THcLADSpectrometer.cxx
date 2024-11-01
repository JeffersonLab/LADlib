#include "THcLADSpectrometer.h"
#include "THaNonTrackingDetector.h"
#include "TDatime.h"
#include "TList.h"
#include "THaTrack.h"

using namespace std;

THcLADSpectrometer::THcLADSpectrometer( const char* name, const char* description ) :
  THaApparatus( name, description )
{
  // default constructor 
  fTracks = 0;
  fNonTrackingDetectors = new TList;

}

//____________________________________________________________________
THcLADSpectrometer::~THcLADSpectrometer()
{
  // Destructor

  DefineVariables( kDelete );
  
}

//____________________________________________________________________
void THcLADSpectrometer::ListInit()
{

  fNonTrackingDetectors->Clear();
  TIter next(fDetectors);
  while( THaDetector* theDetector =
	 static_cast<THaDetector*>( next() )) {

    // We don't really have a tracking detector 

    fNonTrackingDetectors->Add( theDetector );
  }

  fListInit = kTRUE;  
}

//____________________________________________________________________
Int_t THcLADSpectrometer::DefineVariables( EMode mode )
{
  if (mode == kDefine && fIsSetup) return kOK;
  fIsSetup = (mode == kDefine);

  return kOK;
}

//____________________________________________________________________
Int_t THcLADSpectrometer::ReadDatabase( const TDatime& date )
{

  return kOK;

}

//____________________________________________________________________
/*
Int_t THcLADSpectrometer::CoarseTrack()
{

  //  THaSpectrometer::CoarseTrack();
  TIter next( fTrackingDetectors );
  while( THaTrackingDetector* theTrackDetector =
	 static_cast<THaTrackingDetector*>( next() )) {
#ifdef WITH_DEBUG
    if( fDebug >1 ) cout << "Call CoarseProcess() for "
			 << theTrackDetector->GetName() << "... ";
#endif
    theTrackDetector->CoarseTrack( *fTracks );
#ifdef WITH_DEBUG
    if ( fDebug>1 ) cout << "done.\n";
#endif
  }

  return 0;

}
*/
//____________________________________________________________________

Int_t THcLADSpectrometer::Decode( const THaEvData& evdata )
{

  return THaApparatus::Decode(evdata);

}

//____________________________________________________________________
Int_t THcLADSpectrometer::Reconstruct()
{

  TIter next( fNonTrackingDetectors );
  while( THaNonTrackingDetector* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
#ifdef WITH_DEBUG
    if( fDebug > 1 ) cout << "Call FineProcess() for"
			  << theNonTrackDetector->GetName() << "... ";
#endif
    theNonTrackDetector->FineProcess( *fTracks );
#ifdef WITH_DEBUG
    if( fDebug > 1 ) cout << "done.\n";
#endif
  }

  return 0;
}

//____________________________________________________________________
Int_t THcLADSpectrometer::CoarseReconstruct()
{

  if( !fListInit )
    ListInit();

  TIter next( fNonTrackingDetectors );
  while( THaNonTrackingDetector* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
#ifdef WITH_DEBUG
    if( fDebug >1 ) cout << "Call CoarseProcess() for "
			 << theNonTrackDetector->GetName() << "... ";
#endif
    theNonTrackDetector->CoarseProcess( *fTracks );
#ifdef WITH_DEBUG
    if ( fDebug>1 ) cout << "done.\n";
#endif
  }

  return 0;

}

//____________________________________________________________________
THcLADSpectrometer::THcLADSpectrometer() {}

ClassImp(THcLADSpectrometer)

