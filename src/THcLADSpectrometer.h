#ifndef THcLADSpectrometer_h
#define THcLADSpectrometer_h

/////////////////////////////////////////////////
//
// Skeleton class for HallC LAD spectrometer
//
/////////////////////////////////////////////////

#include "THaApparatus.h"
#include "TClonesArray.h"

class THcLADSpectrometer : public THaApparatus {

 public:

  THcLADSpectrometer( const char* name, const char* description );
  virtual ~THcLADSpectrometer();
  
  virtual Int_t CoarseReconstruct();
  virtual Int_t Reconstruct();
  virtual Int_t Decode( const THaEvData& );

 protected:

  TClonesArray* fTracks;
  Int_t   fNtracks;

  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t ReadDatabase( const TDatime& date );

  TList* fNonTrackingDetectors;

 private:
  Bool_t fListInit;
  void   ListInit();

  THcLADSpectrometer();
  ClassDef(THcLADSpectrometer,0)
    
};

#endif /* THcLADSpectrometer_h */
