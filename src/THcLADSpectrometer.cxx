#include "THcLADSpectrometer.h"
#include "TList.h"
#include "THaTrack.h"

using namespace std;

THcLADSpectrometer::THcLADSpectrometer( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // default constructor 

}

//____________________________________________________________________
THcLADSpectrometer::~THcLADSpectrometer()
{
  // Destructor

  DefineVariables( kDelete );
  
}

//____________________________________________________________________
void THcLADSpectrometer::Clear( Option_t* opt )
{

}

//____________________________________________________________________
Int_t THcLADSpectrometer::DefineVariables( EMode mode )
{

  return 0;
}

//____________________________________________________________________
Int_t THcLADSpectrometer::ReadDatabase( const TDatime& date )
{

  return kOK;

}

//____________________________________________________________________
Int_t THcLADSpectrometer::CoarseTrack()
{

  THaSpectrometer::CoarseTrack();
  return 0;

}


//____________________________________________________________________

Int_t THcLADSpectrometer::Decode( const THaEvData& evdata )
{

  return THaApparatus::Decode(evdata);

}

//____________________________________________________________________
Int_t THcLADSpectrometer::Track()
{

  THaSpectrometer::Track();
  return 0;

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
Int_t THcLADSpectrometer::FindVertices( TClonesArray& tracks )
{

  fNtracks = tracks.GetLast() + 1;
  
  // add more 
  // HallC spectrometer chooses between 3 different method to select best track
  // depends
  // (see THcHallCSpectrometer)
  
  return 0;
}

//____________________________________________________________________
Int_t THcLADSpectrometer::TrackCalc()
{

  for( Int_t t = 0; t < fTracks->GetLast()+1; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( fTracks->At(t) );
    //do nothing
  }
  return 0;

}

ClassImp(THcLADSpectrometer)

