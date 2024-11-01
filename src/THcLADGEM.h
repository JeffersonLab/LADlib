#ifndef THcLADGEM_h
#define THcLADGEM_h

#include "THaNonTrackingDetector.h"
#include "THcHitList.h"
#include "THcLADGEMPlane.h"
#include "THcLADGEMModule.h"

class THcLADGEM : public THaNonTrackingDetector, public THcHitList {
 public:

  THcLADGEM( const char* name, const char* description = "",
	     THaApparatus* apparatus = nullptr );
  virtual ~THcLADGEM();

  virtual Int_t   Decode( const THaEvData& );
  virtual EStatus Init( const TDatime& date );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  //  virtual Int_t   End( THaRunBase* r=0 );

 protected:
  
  void            ClearEvent();
  void            Setup(const char* name, const char* description);
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   ReadDatabase( const TDatime& date );

  Int_t fNModules;  // total number of modules
  Int_t fNPlanes;   // total number of GEM planes
  Int_t fNhits;


  // pointer to global var indicatiing whether this spectrometer is triggered
  // for this event
  Bool_t* fPresentP;

  vector<THcLADGEMModule*> fModules;
  vector<THcLADGEMPlane*> fPlanes;

  bool fModulesInitialized;

  ClassDef(THcLADGEM,0)
};

#endif /* THcLADGEM_h */
