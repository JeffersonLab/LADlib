#ifndef THcLADSpectrometer_h
#define THcLADSpectrometer_h

/////////////////////////////////////////////////
//
// Skeleton class for HallC LAD spectrometer
//
/////////////////////////////////////////////////

#include "THaAnalysisObject.h"
#include "THaApparatus.h"
#include "TClonesArray.h"

class THcLADSpectrometer : public THaApparatus {

 public:

  THcLADSpectrometer( const char* name, const char* description );
  virtual ~THcLADSpectrometer();
  
  virtual Int_t CoarseReconstruct();
  virtual Int_t Reconstruct();

  enum EStagesDone {
    kCoarseTrack = BIT(0),
    kCoarseRecon = BIT(1),
    kTracking    = BIT(2),
    kReconstruct = BIT(3)
  };

 protected:

  TClonesArray* fTracks;
  Int_t   fNtracks;

  UInt_t fStagesDone;

  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t ReadDatabase( const TDatime& date );

  TList* fNonTrackingDetectors;
  //  TList* fTrackingDetectors;

 private:
  Bool_t fListInit;
  void   ListInit();

  THcLADSpectrometer();
  ClassDef(THcLADSpectrometer,0)
    
};

#endif /* THcLADSpectrometer_h */
