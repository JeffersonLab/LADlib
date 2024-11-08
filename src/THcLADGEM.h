#ifndef THcLADGEM_h
#define THcLADGEM_h

#include "THaNonTrackingDetector.h"
#include "THcHitList.h"
#include "THcLADGEMPlane.h"
#include "THcLADGEMModule.h"

struct ClusterOutputData {
  // handle 1d cluster output variables  
  std::vector<Int_t> layer;
  std::vector<Int_t> mpdid;
  std::vector<Int_t> axis;
  std::vector<Int_t> nstrip;
  std::vector<Int_t> maxstrip;
  std::vector<Double_t> adc;
  std::vector<Double_t> time;
  void clear(){
    layer.clear();
    mpdid.clear();
    axis.clear();
    nstrip.clear();
    maxstrip.clear();
    adc.clear();
    time.clear();
  }
};

class THcLADGEM : public THaNonTrackingDetector, public THcHitList {
 public:

  THcLADGEM( const char* name, const char* description = "",
	     THaApparatus* apparatus = nullptr );
  virtual ~THcLADGEM();

  virtual Int_t   Decode( const THaEvData& );
  virtual EStatus Init( const TDatime& date );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual void    Clear( Option_t* opt="" );
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

  // Cluster data output handling
  ClusterOutputData fClusOutData;

 class GEM2DHits {
  public:
   Int_t    layer;
   Double_t posX;
   Double_t posY;
   //    Double_t posZ;
   Double_t TimeMean; // average time
   Double_t TimeDiff;
   Double_t TimeCorr;
   Bool_t   IsGoodHit;
   //    Bool_t   Filtered;
   Double_t ADCMean; // average adc sum
   Double_t ADCasym;
  };
  std::vector<GEM2DHits> f2DHits;

 public:
  void Add2DHits(Int_t ilayer, Double_t x, Double_t y,
		 Double_t t, Double_t dt, Double_t tc,
		 Bool_t goodhit, Double_t adc, Double_t adcasy);
  std::vector<GEM2DHits> Get2DHits() { return f2DHits; }

  ClassDef(THcLADGEM,0)
};

#endif /* THcLADGEM_h */
