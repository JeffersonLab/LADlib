#ifndef THcLADHodoscope_h
#define THcLADHodoscope_h

#include "THaNonTrackingDetector.h"
#include "THcHitList.h"
#include "THcLADHodoHit.h"
#include "THcRawHodoHit.h"
#include "THcLADHodoPlane.h"
#include "TClonesArray.h"

class THcLADHodoscope : public THaNonTrackingDetector, public THcHitList {
 public:

  THcLADHodoscope( const char* name, const char* description = "",
		   THaApparatus* apparatus = nullptr );
  virtual ~THcLADHodoscope();

  virtual Int_t   Decode( const THaEvData& );
  virtual EStatus Init( const TDatime& date );
  virtual void    Clear( Option_t* opt="" );
  virtual Int_t   End(THaRunBase* run=0);

  virtual Int_t CoarseProcess( TClonesArray& tracks );
  virtual Int_t FineProcess( TClonesArray& tracks );


  Int_t GetScinIndex(Int_t nPlane, Int_t nPaddle) { return fNPlanes*nPaddle + nPlane; }
  Int_t GetScinIndex(Int_t nSide,  Int_t nPlane, Int_t nPaddle) { return nSide*fMaxHodoScin+fNPlanes*nPaddle+nPlane-1; }


  THcLADHodoPlane* GetPlane(Int_t iplane) { return fPlanes[iplane]; }

  Int_t     GetNPlanes() { return fNPlanes; }
  Int_t     GetTdcOffset(Int_t iplane)       const { return fTdcOffset[iplane]; }
  Double_t  GetHodoSlop(Int_t iplane)        const { return fHodoSlop[iplane]; }
  Double_t  GetAdcTdcOffset(Int_t iplane)    const { return fAdcTdcOffset[iplane]; }
  Double_t  GetTdcMin()                      const { return fScinTdcMin;}
  Double_t  GetTdcMax()                      const { return fScinTdcMax;}
  Double_t  GetTdcToTime()                   const { return fScinTdcToTime; }

  Double_t  GetHodoTopAdcTimeWindowMax(Int_t iii) const { return fHodoTopAdcTimeWindowMax[iii]; }
  Double_t  GetHodoTopAdcTimeWindowMin(Int_t iii) const { return fHodoTopAdcTimeWindowMin[iii]; }
  Double_t  GetHodoBtmAdcTimeWindowMax(Int_t iii) const { return fHodoBtmAdcTimeWindowMax[iii]; } 
  Double_t  GetHodoBtmAdcTimeWindowMin(Int_t iii) const { return fHodoBtmAdcTimeWindowMin[iii]; }

  Double_t GetHodoVelLight(Int_t iii) const {return fHodoVelLight[iii];}


  //Time walk
  Double_t GetHodoVelFit(Int_t iii) const {return fHodoVelFit[iii];}
  Double_t GetHodoCableFit(Int_t iii) const {return fHodoCableFit[iii];}
  Double_t GetHodoLCoeff(Int_t iii) const {return fHodo_LCoeff[iii];}
  Double_t GetHodoTop_c1(Int_t iii) const {return fHodoTop_c1[iii];}
  Double_t GetHodoBtm_c1(Int_t iii) const {return fHodoBtm_c1[iii];}
  Double_t GetHodoTop_c2(Int_t iii) const {return fHodoTop_c2[iii];}
  Double_t GetHodoBtm_c2(Int_t iii) const {return fHodoBtm_c2[iii];}

 protected:

  Int_t     fNPlanes;
  Int_t     fNHits;
  Int_t    *fNPaddle;
  Int_t     fMaxHodoScin;
  Int_t     fCosmicFlag;
  Double_t  fScinTdcMin, fScinTdcMax;
  Double_t  fScinTdcToTime;
  Int_t    *fTdcOffset;
  Double_t *fAdcTdcOffset;
  Double_t *fHodoSlop;

  Double_t *fHodoVelLight;

  //Time walk
  Double_t *fHodoVelFit;
  Double_t *fHodoCableFit;
  Double_t *fHodo_LCoeff;
  Double_t *fHodoTop_c1;
  Double_t *fHodoBtm_c1;
  Double_t *fHodoTop_c2;
  Double_t *fHodoBtm_c2;


  Int_t     fAnalyzePedestals;

  Double_t *fHodoBtmAdcTimeWindowMin;    // per element? per plane?
  Double_t *fHodoBtmAdcTimeWindowMax;
  Double_t *fHodoTopAdcTimeWindowMin;    
  Double_t *fHodoTopAdcTimeWindowMax;
  
  THcLADHodoPlane** fPlanes;
  char** fPlaneNames;

  Bool_t   *fPresentP;

  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t ReadDatabase( const TDatime& date );
  void    Setup(const char* name, const char* description);

  ClassDef(THcLADHodoscope,0)

};

#endif /* THcLADHodoscope_h */
