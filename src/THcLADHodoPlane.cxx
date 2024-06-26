#include "THcLADHodoPlane.h"
#include "THcGlobals.h"
#include "THcHitList.h"
#include "THcLADHodoHit.h"
#include "THcLADHodoscope.h"
#include "THcParmList.h"
#include "THcRawAdcHit.h"
#include "THcRawTdcHit.h"
#include "THcSignalHit.h"

// This is similar to THcScintillatorPlane from hcana

//_______________________________________________________________________________________
THcLADHodoPlane::THcLADHodoPlane(const char *name, const char *description, const Int_t planenum,
                                 THaDetectorBase *parent)
    : THaSubDetector(name, description, parent) {
  // constructor

  fPlaneNum = planenum;
  fNScinHits = 0;
  fNGoodHits = 0;
  fNScinGoodHits = 0;

  fHodoHits = new TClonesArray("THcLADHodoHit", 16);

  frPosAdcErrorFlag = new TClonesArray("THcSignalHit", 16);
  frNegAdcErrorFlag = new TClonesArray("THcSignalHit", 16);

  frPosTDCHits = new TClonesArray("THcSignalHit", 16);
  frNegTDCHits = new TClonesArray("THcSignalHit", 16);
  frPosADCHits = new TClonesArray("THcSignalHit", 16);
  frNegADCHits = new TClonesArray("THcSignalHit", 16);
  frPosADCSums = new TClonesArray("THcSignalHit", 16);
  frNegADCSums = new TClonesArray("THcSignalHit", 16);
  frPosADCPeds = new TClonesArray("THcSignalHit", 16);
  frNegADCPeds = new TClonesArray("THcSignalHit", 16);

  frPosTdcTimeRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcPedRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcPulseIntRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcPulseAmpRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcPulseTimeRaw = new TClonesArray("THcSignalHit", 16);

  frPosTdcTime = new TClonesArray("THcSignalHit", 16);
  frPosAdcPed = new TClonesArray("THcSignalHit", 16);
  frPosAdcPulseInt = new TClonesArray("THcSignalHit", 16);
  frPosAdcPulseAmp = new TClonesArray("THcSignalHit", 16);
  frPosAdcPulseTime = new TClonesArray("THcSignalHit", 16);

  frNegTdcTimeRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcPedRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcPulseIntRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcPulseAmpRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcPulseTimeRaw = new TClonesArray("THcSignalHit", 16);

  frNegTdcTime = new TClonesArray("THcSignalHit", 16);
  frNegAdcPed = new TClonesArray("THcSignalHit", 16);
  frNegAdcPulseInt = new TClonesArray("THcSignalHit", 16);
  frNegAdcPulseAmp = new TClonesArray("THcSignalHit", 16);
  frNegAdcPulseTime = new TClonesArray("THcSignalHit", 16);

  frPosAdcSampPedRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcSampPulseIntRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcSampPulseAmpRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcSampPulseTimeRaw = new TClonesArray("THcSignalHit", 16);
  frPosAdcSampPed = new TClonesArray("THcSignalHit", 16);
  frPosAdcSampPulseInt = new TClonesArray("THcSignalHit", 16);
  frPosAdcSampPulseAmp = new TClonesArray("THcSignalHit", 16);
  frPosAdcSampPulseTime = new TClonesArray("THcSignalHit", 16);

  frNegAdcSampPedRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcSampPulseIntRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcSampPulseAmpRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcSampPulseTimeRaw = new TClonesArray("THcSignalHit", 16);
  frNegAdcSampPed = new TClonesArray("THcSignalHit", 16);
  frNegAdcSampPulseInt = new TClonesArray("THcSignalHit", 16);
  frNegAdcSampPulseAmp = new TClonesArray("THcSignalHit", 16);
  frNegAdcSampPulseTime = new TClonesArray("THcSignalHit", 16);
}

//_______________________________________________________________________________________
THcLADHodoPlane::~THcLADHodoPlane() {
  // destructor
  if (fIsSetup)
    RemoveVariables();

  delete frPosAdcErrorFlag;
  frPosAdcErrorFlag = NULL;
  delete frNegAdcErrorFlag;
  frNegAdcErrorFlag = NULL;

  delete fHodoHits;
  delete frPosTDCHits;
  delete frNegTDCHits;
  delete frPosADCHits;
  delete frNegADCHits;
  delete frPosADCSums;
  delete frNegADCSums;
  delete frPosADCPeds;
  delete frNegADCPeds;

  delete frPosTdcTimeRaw;
  delete frPosAdcPedRaw;
  delete frPosAdcPulseIntRaw;
  delete frPosAdcPulseAmpRaw;
  delete frPosAdcPulseTimeRaw;

  delete frPosTdcTime;
  delete frPosAdcPed;
  delete frPosAdcPulseInt;
  delete frPosAdcPulseAmp;
  delete frPosAdcPulseTime;

  delete frNegTdcTimeRaw;
  delete frNegAdcPedRaw;
  delete frNegAdcPulseIntRaw;
  delete frNegAdcPulseAmpRaw;
  delete frNegAdcPulseTimeRaw;

  delete frNegTdcTime;
  delete frNegAdcPed;
  delete frNegAdcPulseInt;
  delete frNegAdcPulseAmp;
  delete frNegAdcPulseTime;

  delete frPosAdcSampPedRaw;
  frPosAdcSampPedRaw = NULL;
  delete frPosAdcSampPulseIntRaw;
  frPosAdcSampPulseIntRaw = NULL;
  delete frPosAdcSampPulseAmpRaw;
  frPosAdcSampPulseAmpRaw = NULL;
  delete frPosAdcSampPulseTimeRaw;
  frPosAdcSampPulseTimeRaw = NULL;
  delete frPosAdcSampPed;
  frPosAdcSampPed = NULL;
  delete frPosAdcSampPulseInt;
  frPosAdcSampPulseInt = NULL;
  delete frPosAdcSampPulseAmp;
  frPosAdcSampPulseAmp = NULL;
  delete frPosAdcSampPulseTime;
  frPosAdcSampPulseTime = NULL;

  delete frNegAdcSampPedRaw;
  frNegAdcSampPedRaw = NULL;
  delete frNegAdcSampPulseIntRaw;
  frNegAdcSampPulseIntRaw = NULL;
  delete frNegAdcSampPulseAmpRaw;
  frNegAdcSampPulseAmpRaw = NULL;
  delete frNegAdcSampPulseTimeRaw;
  frNegAdcSampPulseTimeRaw = NULL;
  delete frNegAdcSampPed;
  frNegAdcSampPed = NULL;
  delete frNegAdcSampPulseInt;
  frNegAdcSampPulseInt = NULL;
  delete frNegAdcSampPulseAmp;
  frNegAdcSampPulseAmp = NULL;
  delete frNegAdcSampPulseTime;
  frNegAdcSampPulseTime = NULL;

  delete[] fPosPedSum;
  fPosPedSum = 0;
  delete[] fPosPedSum2;
  fPosPedSum2 = 0;
  delete[] fPosPedLimit;
  fPosPedLimit = 0;
  delete[] fPosPedCount;
  fPosPedCount = 0;
  delete[] fNegPedSum;
  fNegPedSum = 0;
  delete[] fNegPedSum2;
  fNegPedSum2 = 0;
  delete[] fNegPedLimit;
  fNegPedLimit = 0;
  delete[] fNegPedCount;
  fNegPedCount = 0;
  delete[] fPosPed;
  fPosPed = 0;
  delete[] fNegPed;
  fNegPed = 0;
  delete[] fPosThresh;
  fPosThresh = 0;
  delete[] fNegThresh;
  fNegThresh = 0;

  /*

  delete [] fPosCenter; fPosCenter = 0;

  delete [] fHodoPosMinPh; fHodoPosMinPh = NULL;
  delete [] fHodoNegMinPh; fHodoNegMinPh = NULL;
  delete [] fHodoPosPhcCoeff; fHodoPosPhcCoeff = NULL;
  delete [] fHodoNegPhcCoeff; fHodoNegPhcCoeff = NULL;
  delete [] fHodoPosTimeOffset; fHodoPosTimeOffset = NULL;
  delete [] fHodoNegTimeOffset; fHodoNegTimeOffset = NULL;
  delete [] fHodoPosInvAdcOffset; fHodoPosInvAdcOffset = NULL;
  delete [] fHodoNegInvAdcOffset; fHodoNegInvAdcOffset = NULL;
  delete [] fHodoPosInvAdcLinear; fHodoPosInvAdcLinear = NULL;
  delete [] fHodoNegInvAdcLinear; fHodoNegInvAdcLinear = NULL;
  delete [] fHodoPosAdcTimeWindowMax; fHodoPosAdcTimeWindowMax = NULL;
  delete [] fHodoPosAdcTimeWindowMin; fHodoPosAdcTimeWindowMin = NULL;
  delete [] fHodoNegAdcTimeWindowMax; fHodoNegAdcTimeWindowMax = NULL;
  delete [] fHodoNegAdcTimeWindowMin; fHodoNegAdcTimeWindowMin = NULL;
  delete [] fHodoNegInvAdcLinear; fHodoNegInvAdcLinear = NULL;
  delete [] fHodoPosInvAdcAdc; fHodoPosInvAdcAdc = NULL;
  delete [] fHodoNegInvAdcAdc; fHodoNegInvAdcAdc = NULL;
  delete [] fHodoVelFit;                 fHodoVelFit = NULL;
  delete [] fHodoCableFit;               fHodoCableFit = NULL;
  delete [] fHodo_LCoeff;                fHodo_LCoeff = NULL;
  delete [] fHodoPos_c1;                 fHodoPos_c1 = NULL;
  delete [] fHodoNeg_c1;                 fHodoNeg_c1 = NULL;
  delete [] fHodoPos_c2;                 fHodoPos_c2 = NULL;
  delete [] fHodoNeg_c2;                 fHodoNeg_c2 = NULL;


  delete [] fHodoVelLight; fHodoVelLight = NULL;
  delete [] fHodoSigma; fHodoSigma = NULL;

  */
}

//_______________________________________________________________________________________
THaAnalysisObject::EStatus THcLADHodoPlane::Init(const TDatime &date) {

  if (IsZombie())
    return fStatus = kInitError;

  EStatus status;
  if ((status = THaSubDetector::Init(date)))
    return fStatus = status;
  return fStatus = kOK;
}

//_______________________________________________________________________________________
void THcLADHodoPlane::Clear(Option_t *opt) {

  fHodoHits->Clear();

  frPosAdcErrorFlag->Clear();
  frNegAdcErrorFlag->Clear();

  // Ref time
  fPosTdcRefTime = kBig;
  fPosAdcRefTime = kBig;
  fNegTdcRefTime = kBig;
  fNegAdcRefTime = kBig;
  fPosTdcRefDiffTime = kBig;
  fPosAdcRefDiffTime = kBig;
  fNegTdcRefDiffTime = kBig;
  fNegAdcRefDiffTime = kBig;

  // Several arrays, vectors reset here
  frPosTdcTimeRaw->Clear();
  frPosAdcPedRaw->Clear();
  frPosAdcPulseIntRaw->Clear();
  frPosAdcPulseAmpRaw->Clear();
  frPosAdcPulseTimeRaw->Clear();

  frPosTdcTime->Clear();
  frPosAdcPed->Clear();
  frPosAdcPulseInt->Clear();
  frPosAdcPulseAmp->Clear();
  frPosAdcPulseTime->Clear();

  frPosAdcSampPedRaw->Clear();
  frPosAdcSampPulseIntRaw->Clear();
  frPosAdcSampPulseAmpRaw->Clear();
  frPosAdcSampPulseTimeRaw->Clear();
  frPosAdcSampPed->Clear();
  frPosAdcSampPulseInt->Clear();
  frPosAdcSampPulseAmp->Clear();
  frPosAdcSampPulseTime->Clear();

  frNegTdcTimeRaw->Clear();
  frNegAdcPedRaw->Clear();
  frNegAdcPulseIntRaw->Clear();
  frNegAdcPulseAmpRaw->Clear();
  frNegAdcPulseTimeRaw->Clear();

  frNegTdcTime->Clear();
  frNegAdcPed->Clear();
  frNegAdcPulseInt->Clear();
  frNegAdcPulseAmp->Clear();
  frNegAdcPulseTime->Clear();

  frNegAdcSampPedRaw->Clear();
  frNegAdcSampPulseIntRaw->Clear();
  frNegAdcSampPulseAmpRaw->Clear();
  frNegAdcSampPulseTimeRaw->Clear();
  frNegAdcSampPed->Clear();
  frNegAdcSampPulseInt->Clear();
  frNegAdcSampPulseAmp->Clear();
  frNegAdcSampPulseTime->Clear();

  // Waveform vectors
  fPosAdcSampWaveform.clear();
  fNegAdcSampWaveform.clear();

  fCluster.clear();

  // Reset counters

  // Adc good hits
  fTotNumGoodPosAdcHits = 0;
  fTotNumGoodNegAdcHits = 0;
  fTotNumGoodAdcHits = 0;

  fTotNumPosAdcHits = 0;
  fTotNumNegAdcHits = 0;
  fTotNumAdcHits = 0;

  // Tdc good hits
  fTotNumGoodPosTdcHits = 0;
  fTotNumGoodNegTdcHits = 0;
  fTotNumGoodTdcHits = 0;

  fTotNumPosTdcHits = 0;
  fTotNumNegTdcHits = 0;
  fTotNumTdcHits = 0;

  fNScinHits = 0;
  fNGoodHits = 0;     // looks like not being used.. remove it?
  fNScinGoodHits = 0; // looks like not being used.. remove it?

  // Clear occupancies
  for (UInt_t ielem = 0; ielem < fNumGoodPosAdcHits.size(); ielem++)
    fNumGoodPosAdcHits.at(ielem) = 0;
  for (UInt_t ielem = 0; ielem < fNumGoodNegAdcHits.size(); ielem++)
    fNumGoodNegAdcHits.at(ielem) = 0;

  for (UInt_t ielem = 0; ielem < fNumGoodPosTdcHits.size(); ielem++)
    fNumGoodPosTdcHits.at(ielem) = 0;
  for (UInt_t ielem = 0; ielem < fNumGoodNegTdcHits.size(); ielem++)
    fNumGoodNegTdcHits.at(ielem) = 0;

  // Clear Ped/Amps/Int/Time
  for (UInt_t ielem = 0; ielem < fGoodPosAdcPed.size(); ielem++) {
    fGoodPosAdcPed.at(ielem) = 0.0;
    fGoodPosAdcMult.at(ielem) = 0.0;
    fGoodPosAdcHitUsed.at(ielem) = 0.0;
    fGoodPosAdcPulseInt.at(ielem) = 0.0;
    fGoodPosAdcPulseAmp.at(ielem) = 0.0;
    fGoodPosAdcPulseTime.at(ielem) = kBig;
    fGoodPosAdcTdcDiffTime.at(ielem) = kBig;
  }
  for (UInt_t ielem = 0; ielem < fGoodNegAdcPed.size(); ielem++) {
    fGoodNegAdcPed.at(ielem) = 0.0;
    fGoodNegAdcMult.at(ielem) = 0.0;
    fGoodNegAdcHitUsed.at(ielem) = 0.0;
    fGoodNegAdcPulseInt.at(ielem) = 0.0;
    fGoodNegAdcPulseAmp.at(ielem) = 0.0;
    fGoodNegAdcPulseTime.at(ielem) = kBig;
    fGoodNegAdcTdcDiffTime.at(ielem) = kBig;
  }

  /*
  //Clear Good TDC Variables
  for (UInt_t ielem = 0; ielem < fGoodPosTdcTimeUnCorr.size(); ielem++) {
    fGoodPosTdcTimeUnCorr.at(ielem)    = kBig;
    fGoodPosTdcTimeCorr.at(ielem)      = kBig;
    fGoodPosTdcTimeTOFCorr.at(ielem)   = kBig;
    fGoodPosTdcTimeWalkCorr.at(ielem)  = kBig;
  }

  for (UInt_t ielem = 0; ielem < fGoodNegTdcTimeUnCorr.size(); ielem++) {
    fGoodNegTdcTimeUnCorr.at(ielem)    = kBig;
    fGoodNegTdcTimeCorr.at(ielem)      = kBig;
    fGoodNegTdcTimeTOFCorr.at(ielem)   = kBig;
    fGoodNegTdcTimeWalkCorr.at(ielem)  = kBig;
  }

  for (UInt_t ielem = 0; ielem < fGoodDiffDistTrack.size(); ielem++) {
    fGoodDiffDistTrack.at(ielem) = kBig;
  }

  fpTime = -1.e4;
  fHitDistance = kBig;
  fScinYPos = kBig;
  fScinXPos = kBig;
  fTrackXPosition = kBig;
  fTrackYPosition = kBig;
  */
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::Decode(const THaEvData &evdata) {
  // Leave it to do nothing

  return 0;
}

//_______________________________________________________________________________________

Int_t THcLADHodoPlane::ReadDatabase(const TDatime &date) {

  // Read database files as needed here
  char prefix[2];
  prefix[0] = tolower(GetParent()->GetPrefix()[0]);
  prefix[1] = '\0';

  // Get # of element for each hodo detector
  string parname = "hodo_" + string(GetName()) + "_nr";
  DBRequest list_1[] = {{parname.c_str(), &fNelem, kInt}, {0}};
  gHcParms->LoadParmValues(list_1, prefix);

  DBRequest list[] = {{"hodo_adc_mode", &fADCMode, kInt, 0, 1},
                      {"hodo_adc_diag_cut", &fADCDiagCut, kInt, 0, 1},
                      {"cosmicflag", &fCosmicFlag, kInt, 0, 1},
                      {"hodo_debug_adc", &fDebugAdc, kInt, 0, 1},
                      {"hodo_SampThreshold", &fSampThreshold, kDouble, 0, 1},
                      {"hodo_SampNSA", &fSampNSA, kInt, 0, 1},
                      {"hodo_SampNSAT", &fSampNSAT, kInt, 0, 1},
                      {"hodo_SampNSB", &fSampNSB, kInt, 0, 1},
                      {"hodo_OutputSampWaveform", &fOutputSampWaveform, kInt, 0, 1},
                      {"hodo_UseSampWaveform", &fUseSampWaveform, kInt, 0, 1},
                      {0}};

  // Set Default values

  fDebugAdc = 0; // Set ADC debug parameter to false unless set in parameter file
  fADCMode = kADCDynamicPedestal;
  fADCDiagCut = 50.0;
  fCosmicFlag = 0;
  fSampThreshold = 5.;
  fSampNSA = 0;            // use value stored in event 125 info
  fSampNSB = 0;            // use value stored in event 125 info
  fSampNSAT = 2;           // default value in THcRawHit::SetF250Params
  fOutputSampWaveform = 0; // 0= no output , 1 = output Sample Waveform
  fUseSampWaveform = 0;    // 0= do not use , 1 = use Sample Waveform

  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  if (fCosmicFlag == 1)
    cout << " setup for cosmics in scint plane" << endl;

  // Retrieve parameters we need from parent class
  // Common for all planes

  THcLADHodoscope *parent = (THcLADHodoscope *)GetParent();

  fHodoSlop = parent->GetHodoSlop(fPlaneNum - 1);
  fTdcOffset = parent->GetTdcOffset(fPlaneNum - 1);
  fAdcTdcOffset = parent->GetAdcTdcOffset(fPlaneNum - 1);
  fScinTdcMin = parent->GetTdcMin();
  fScinTdcMax = parent->GetTdcMax();
  fScinTdcToTime = parent->GetTdcToTime();

  // fTofTolerance=parent->GetTofTolerance();
  // fBetaNominal=parent->GetBetaNominal();
  // fStartTimeCenter=parent->GetStartTimeCenter();
  // fStartTimeSlop=parent->GetStartTimeSlop();

  // Parameters for this plane
  fHodoPosAdcTimeWindowMin = new Double_t[fNelem];
  fHodoNegAdcTimeWindowMin = new Double_t[fNelem];
  fHodoPosAdcTimeWindowMax = new Double_t[fNelem];
  fHodoNegAdcTimeWindowMax = new Double_t[fNelem];

  /*
  fHodoPosMinPh = new Double_t[fNelem];
  fHodoNegMinPh = new Double_t[fNelem];
  fHodoPosPhcCoeff = new Double_t[fNelem];
  fHodoNegPhcCoeff = new Double_t[fNelem];
  fHodoPosTimeOffset = new Double_t[fNelem];
  fHodoNegTimeOffset = new Double_t[fNelem];
  fHodoVelLight = new Double_t[fNelem];
  fHodoPosInvAdcOffset = new Double_t[fNelem];
  fHodoNegInvAdcOffset = new Double_t[fNelem];
  fHodoPosInvAdcLinear = new Double_t[fNelem];
  fHodoNegInvAdcLinear = new Double_t[fNelem];
  fHodoPosInvAdcAdc = new Double_t[fNelem];
  fHodoNegInvAdcAdc = new Double_t[fNelem];
  fHodoSigma = new Double_t[fNelem];

  //New Time-Walk Calibration Parameters
  fHodoVelFit=new Double_t [fNelem];
  fHodoCableFit=new Double_t [fNelem];
  fHodo_LCoeff=new Double_t [fNelem];
  fHodoPos_c1=new Double_t [fNelem];
  fHodoNeg_c1=new Double_t [fNelem];
  fHodoPos_c2=new Double_t [fNelem];
  fHodoNeg_c2=new Double_t [fNelem];
  */

  for (Int_t j = 0; j < (Int_t)fNelem; j++) {
    Int_t index = parent->GetScinIndex(fPlaneNum - 1, j);
    fHodoPosAdcTimeWindowMin[j] = parent->GetHodoPosAdcTimeWindowMin(index);
    fHodoPosAdcTimeWindowMax[j] = parent->GetHodoPosAdcTimeWindowMax(index);
    fHodoNegAdcTimeWindowMin[j] = parent->GetHodoNegAdcTimeWindowMin(index);
    fHodoNegAdcTimeWindowMax[j] = parent->GetHodoNegAdcTimeWindowMax(index);
    /*
    fHodoPosMinPh[j] = parent->GetHodoPosMinPh(index);
    fHodoNegMinPh[j] = parent->GetHodoNegMinPh(index);
    fHodoPosPhcCoeff[j] = parent->GetHodoPosPhcCoeff(index);
    fHodoNegPhcCoeff[j] = parent->GetHodoNegPhcCoeff(index);
    fHodoPosTimeOffset[j] = parent->GetHodoPosTimeOffset(index);
    fHodoNegTimeOffset[j] = parent->GetHodoNegTimeOffset(index);
    fHodoPosInvAdcOffset[j] = parent->GetHodoPosInvAdcOffset(index);
    fHodoNegInvAdcOffset[j] = parent->GetHodoNegInvAdcOffset(index);
    fHodoPosInvAdcLinear[j] = parent->GetHodoPosInvAdcLinear(index);
    fHodoNegInvAdcLinear[j] = parent->GetHodoNegInvAdcLinear(index);
    fHodoPosInvAdcAdc[j] = parent->GetHodoPosInvAdcAdc(index);
    fHodoNegInvAdcAdc[j] = parent->GetHodoNegInvAdcAdc(index);
    fHodoVelLight[j] = parent->GetHodoVelLight(index);

    //Get Time-Walk correction param
    fHodoVelFit[j] = parent->GetHodoVelFit(index);
    fHodoCableFit[j] = parent->GetHodoCableFit(index);
    fHodo_LCoeff[j] =  parent->GetHodoLCoeff(index);
    fHodoPos_c1[j] = parent->GetHodoPos_c1(index);
    fHodoNeg_c1[j] = parent->GetHodoNeg_c1(index);
    fHodoPos_c2[j] = parent->GetHodoPos_c2(index);
    fHodoNeg_c2[j] = parent->GetHodoNeg_c2(index);
    Double_t possigma = parent->GetHodoPosSigma(index);
    Double_t negsigma = parent->GetHodoNegSigma(index);
    fHodoSigma[j] = TMath::Sqrt(possigma*possigma+negsigma*negsigma)/2.0;
    */
  }

  //  fTdc_Thrs = parent->GetTDCThrs();

  // Create arrays to hold results here
  InitializePedestals();

  // Initialize

  fNumGoodPosAdcHits = vector<Int_t>(fNelem, 0.0);
  fNumGoodNegAdcHits = vector<Int_t>(fNelem, 0.0);
  fNumGoodPosTdcHits = vector<Int_t>(fNelem, 0.0);
  fNumGoodNegTdcHits = vector<Int_t>(fNelem, 0.0);

  fGoodPosAdcPed = vector<Double_t>(fNelem, 0.0);
  fGoodNegAdcPed = vector<Double_t>(fNelem, 0.0);
  fGoodPosAdcMult = vector<Double_t>(fNelem, 0.0);
  fGoodNegAdcMult = vector<Double_t>(fNelem, 0.0);
  fGoodPosAdcHitUsed = vector<Double_t>(fNelem, 0.0);
  fGoodNegAdcHitUsed = vector<Double_t>(fNelem, 0.0);
  fGoodPosAdcPulseAmp = vector<Double_t>(fNelem, 0.0);
  fGoodNegAdcPulseAmp = vector<Double_t>(fNelem, 0.0);
  fGoodPosAdcPulseInt = vector<Double_t>(fNelem, 0.0);
  fGoodNegAdcPulseInt = vector<Double_t>(fNelem, 0.0);
  fGoodPosAdcPulseTime = vector<Double_t>(fNelem, 0.0);
  fGoodNegAdcPulseTime = vector<Double_t>(fNelem, 0.0);
  fGoodPosAdcTdcDiffTime = vector<Double_t>(fNelem, 0.0);
  fGoodNegAdcTdcDiffTime = vector<Double_t>(fNelem, 0.0);

  /*
  fGoodPosTdcTimeUnCorr  = vector<Double_t> (fNelem, 0.0);
  fGoodNegTdcTimeUnCorr  = vector<Double_t> (fNelem, 0.0);
  fGoodPosTdcTimeCorr    = vector<Double_t> (fNelem, 0.0);
  fGoodNegTdcTimeCorr    = vector<Double_t> (fNelem, 0.0);
  fGoodPosTdcTimeTOFCorr = vector<Double_t> (fNelem, 0.0);
  fGoodNegTdcTimeTOFCorr = vector<Double_t> (fNelem, 0.0);
  fGoodPosTdcTimeWalkCorr = vector<Double_t> (fNelem, 0.0);
  fGoodNegTdcTimeWalkCorr = vector<Double_t> (fNelem, 0.0);
  fGoodDiffDistTrack      = vector<Double_t> (fNelem, 0.0);
  */

  return 0;
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::DefineVariables(EMode mode) {

  /*
    Initialize global variables for histograms and Root tree
  */

  if (mode == kDefine && fIsSetup)
    return kOK;
  fIsSetup = (mode == kDefine);

  // Register variables in global list

  if (fDebugAdc) {
    RVarDef vars[] = {
        {"posAdcErrorFlag", "Error Flag for When FPGA Fails", "frPosAdcErrorFlag.THcSignalHit.GetData()"},
        {"negAdcErrorFlag", "Error Flag for When FPGA Fails", "frNegAdcErrorFlag.THcSignalHit.GetData()"},

        {"posTdcTimeRaw", "List of positive raw TDC values.", "frPosTdcTimeRaw.THcSignalHit.GetData()"},
        {"posAdcPedRaw", "List of positive raw ADC pedestals", "frPosAdcPedRaw.THcSignalHit.GetData()"},
        {"posAdcPulseIntRaw", "List of positive raw ADC pulse integrals.",
         "frPosAdcPulseIntRaw.THcSignalHit.GetData()"},
        {"posAdcPulseAmpRaw", "List of positive raw ADC pulse amplitudes.",
         "frPosAdcPulseAmpRaw.THcSignalHit.GetData()"},
        {"posAdcPulseTimeRaw", "List of positive raw ADC pulse times.", "frPosAdcPulseTimeRaw.THcSignalHit.GetData()"},

        {"posTdcTime", "List of positive TDC values.", "frPosTdcTime.THcSignalHit.GetData()"},
        {"posAdcPed", "List of positive ADC pedestals", "frPosAdcPed.THcSignalHit.GetData()"},
        {"posAdcPulseInt", "List of positive ADC pulse integrals.", "frPosAdcPulseInt.THcSignalHit.GetData()"},
        {"posAdcPulseAmp", "List of positive ADC pulse amplitudes.", "frPosAdcPulseAmp.THcSignalHit.GetData()"},
        {"posAdcPulseTime", "List of positive ADC pulse times.", "frPosAdcPulseTime.THcSignalHit.GetData()"},

        {"posAdcSampPedRaw", "Positive Raw Samp ADC pedestals", "frPosAdcSampPedRaw.THcSignalHit.GetData()"},
        {"posAdcSampPulseIntRaw", "Positive Raw Samp ADC pulse integrals",
         "frPosAdcSampPulseIntRaw.THcSignalHit.GetData()"},
        {"posAdcSampPulseAmpRaw", "Positive Raw Samp ADC pulse amplitudes",
         "frPosAdcSampPulseAmpRaw.THcSignalHit.GetData()"},
        {"posAdcSampPulseTimeRaw", "Positive Raw Samp ADC pulse times",
         "frPosAdcSampPulseTimeRaw.THcSignalHit.GetData()"},
        {"posAdcSampPed", "Positive Samp ADC pedestals", "frPosAdcSampPed.THcSignalHit.GetData()"},
        {"posAdcSampPulseInt", "Positive Samp ADC pulse integrals", "frPosAdcSampPulseInt.THcSignalHit.GetData()"},
        {"posAdcSampPulseAmp", "Positive Samp ADC pulse amplitudes", "frPosAdcSampPulseAmp.THcSignalHit.GetData()"},
        {"posAdcSampPulseTime", "Positive Samp ADC pulse times", "frPosAdcSampPulseTime.THcSignalHit.GetData()"},

        {"negTdcTimeRaw", "List of negative raw TDC values.", "frNegTdcTimeRaw.THcSignalHit.GetData()"},
        {"negAdcPedRaw", "List of negative raw ADC pedestals", "frNegAdcPedRaw.THcSignalHit.GetData()"},
        {"negAdcPulseIntRaw", "List of negative raw ADC pulse integrals.",
         "frNegAdcPulseIntRaw.THcSignalHit.GetData()"},
        {"negAdcPulseAmpRaw", "List of negative raw ADC pulse amplitudes.",
         "frNegAdcPulseAmpRaw.THcSignalHit.GetData()"},
        {"negAdcPulseTimeRaw", "List of negative raw ADC pulse times.", "frNegAdcPulseTimeRaw.THcSignalHit.GetData()"},

        {"negTdcTime", "List of negative TDC values.", "frNegTdcTime.THcSignalHit.GetData()"},
        {"negAdcPed", "List of negative ADC pedestals", "frNegAdcPed.THcSignalHit.GetData()"},
        {"negAdcPulseInt", "List of negative ADC pulse integrals.", "frNegAdcPulseInt.THcSignalHit.GetData()"},
        {"negAdcPulseAmp", "List of negative ADC pulse amplitudes.", "frNegAdcPulseAmp.THcSignalHit.GetData()"},
        {"negAdcPulseTime", "List of negative ADC pulse times.", "frNegAdcPulseTime.THcSignalHit.GetData()"},

        {"negAdcSampPedRaw", "Negative Raw Samp ADC pedestals", "frNegAdcSampPedRaw.THcSignalHit.GetData()"},
        {"negAdcSampPulseIntRaw", "Negative Raw Samp ADC pulse integrals",
         "frNegAdcSampPulseIntRaw.THcSignalHit.GetData()"},
        {"negAdcSampPulseAmpRaw", "Negative Raw Samp ADC pulse amplitudes",
         "frNegAdcSampPulseAmpRaw.THcSignalHit.GetData()"},
        {"negAdcSampPulseTimeRaw", "Negative Raw Samp ADC pulse times",
         "frNegAdcSampPulseTimeRaw.THcSignalHit.GetData()"},
        {"negAdcSampPed", "Negative Samp ADC pedestals", "frNegAdcSampPed.THcSignalHit.GetData()"},
        {"negAdcSampPulseInt", "Negative Samp ADC pulse integrals", "frNegAdcSampPulseInt.THcSignalHit.GetData()"},
        {"negAdcSampPulseAmp", "Negative Samp ADC pulse amplitudes", "frNegAdcSampPulseAmp.THcSignalHit.GetData()"},
        {"negAdcSampPulseTime", "Negative Samp ADC pulse times", "frNegAdcSampPulseTime.THcSignalHit.GetData()"},

        {"totNumPosAdcHits", "Total Number of Positive ADC Hits",
         "fTotNumPosAdcHits"}, // Hodo+ raw ADC multiplicity Int_t
        {"totNumNegAdcHits", "Total Number of Negative ADC Hits", "fTotNumNegAdcHits"}, // Hodo- raw ADC multiplicity ""
        {"totNumAdcHits", "Total Number of PMTs Hit (as measured by ADCs)",
         "fTotNumAdcHits"}, // Hodo raw ADC multiplicity  ""

        {"totNumPosTdcHits", "Total Number of Positive TDC Hits", "fTotNumPosTdcHits"}, // Hodo+ raw TDC multiplicity ""
        {"totNumNegTdcHits", "Total Number of Negative TDC Hits", "fTotNumNegTdcHits"}, // Hodo- raw TDC multiplicity ""
        {"totNumTdcHits", "Total Number of PMTs Hits (as measured by TDCs)",
         "fTotNumTdcHits"}, // Hodo raw TDC multiplicity  ""
        {0}};
    DefineVarsFromList(vars, mode);
  } // end debug statement

  if (fOutputSampWaveform == 1) {
    RVarDef vars[] = {{"adcNegSampWaveform", "FADC Neg ADCSample Waveform", "fNegAdcSampWaveform"},
                      {"adcPosSampWaveform", "FADC Pos ADCSample Waveform", "fPosAdcSampWaveform"},
                      {0}};
    DefineVarsFromList(vars, mode);
  }

  RVarDef vars[] = {
      {"nhits", "Number of paddle hits (passed TDC && ADC Min and Max cuts for either end)", "GetNScinHits() "},

      {"posTdcCounter", "List of positive TDC counter numbers.",
       "frPosTdcTimeRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo+ raw TDC occupancy
      {"posAdcCounter", "List of positive ADC counter numbers.",
       "frPosAdcPulseIntRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo+ raw ADC occupancy
      {"negTdcCounter", "List of negative TDC counter numbers.",
       "frNegTdcTimeRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo- raw TDC occupancy
      {"negAdcCounter", "List of negative ADC counter numbers.",
       "frNegAdcPulseIntRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo- raw ADC occupancy

      {"numGoodPosAdcHits", "Number of Good Positive ADC Hits Per PMT",
       "fNumGoodPosAdcHits"}, // Hodo+ good ADC occupancy - vector<Int_t>
      {"numGoodNegAdcHits", "Number of Good Negative ADC Hits Per PMT",
       "fNumGoodNegAdcHits"}, // Hodo- good ADC occupancy - vector <Int_t>

      {"numGoodPosTdcHits", "Number of Good Positive TDC Hits Per PMT",
       "fNumGoodPosTdcHits"}, // Hodo+ good TDC occupancy - vector<Int_t>
      {"numGoodNegTdcHits", "Number of Good Negative TDC Hits Per PMT",
       "fNumGoodNegTdcHits"}, // Hodo- good TDC occupancy - vector <Int_t>

      {"totNumGoodPosAdcHits", "Total Number of Good Positive ADC Hits",
       "fTotNumGoodPosAdcHits"}, // Hodo+ good ADC multiplicity - Int_t
      {"totNumGoodNegAdcHits", "Total Number of Good Negative ADC Hits",
       "fTotNumGoodNegAdcHits"}, // Hodo- good ADC multiplicity - Int_t
      {"totNumGoodAdcHits", "TotalNumber of Good ADC Hits Per PMT",
       "fTotNumGoodAdcHits"}, // Hodo good ADC multiplicity - Int_t

      {"totNumGoodPosTdcHits", "Total Number of Good Positive TDC Hits",
       "fTotNumGoodPosTdcHits"}, // Hodo+ good TDC multiplicity - Int_t
      {"totNumGoodNegTdcHits", "Total Number of Good Negative TDC Hits",
       "fTotNumGoodNegTdcHits"}, // Hodo- good TDC multiplicity - Int_t
      {"totNumGoodTdcHits", "TotalNumber of Good TDC Hits Per PMT",
       "fTotNumGoodTdcHits"}, // Hodo good TDC multiplicity - Int_t

      {"GoodPosAdcPed", "List of Positive ADC pedestals (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodPosAdcPed"}, // vector<Double_t>
      {"GoodNegAdcPed", "List of Negative ADC pedestals (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodNegAdcPed"}, // vector<Double_t>
      {"GoodPosAdcMult", "List of Positive ADC Mult (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodPosAdcMult"}, // vector<Double_t>
      {"GoodNegAdcMult", "List of Negative ADC Mult (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodNegAdcMult"}, // vector<Double_t>
      {"GoodPosAdcHitUsed", "List of Positive ADC Hit Used (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodPosAdcHitUsed"}, // vector<Double_t>
      {"GoodNegAdcHitUsed", "List of Negative ADC Hit Used (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodNegAdcHitUsed"}, // vector<Double_t>

      {"GoodPosAdcPulseInt", "List of positive ADC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodPosAdcPulseInt"},
      {"GoodPosAdcPulseAmp", "List of positive ADC peak amp (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodPosAdcPulseAmp"},
      {"GoodNegAdcPulseAmp", "List of Negative ADC peak amp (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodNegAdcPulseAmp"},
      {"GoodPosAdcPulseTime", "List of positive ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodPosAdcPulseTime"},
      {"GoodNegAdcPulseTime", "List of Negative ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodNegAdcPulseTime"},
      {"GoodPosAdcTdcDiffTime", "List of positive TDC - ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodPosAdcTdcDiffTime"},
      {"GoodNegAdcTdcDiffTime", "List of Negative TDC - ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodNegAdcTdcDiffTime"},

      /*
      // cluster variables
      {"NumClus",         "Number of clusters", "fNumberClusters"},
      {"Clus.Pos",        "Position of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterPosition()"},
      {"Clus.Size",       "Size of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterSize()"},
      {"Clus.Flag",       "Flag of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterFlag()"},
      {"Clus.UsedFlag",   "USed Flag of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterUsedFlag()"},
      */

      {"PosTdcRefTime", "Reference time of Pos TDC", "fPosTdcRefTime"},
      {"NegTdcRefTime", "Reference time of Neg TDC", "fNegTdcRefTime"},
      {"PosAdcRefTime", "Reference time of Pos ADC", "fPosAdcRefTime"},
      {"NegAdcRefTime", "Reference time of Neg aDC", "fNegAdcRefTime"},
      {"PosTdcRefDiffTime", "Reference Diff time of Pos TDC", "fPosTdcRefDiffTime"},
      {"NegTdcRefDiffTime", "Reference Diff time of Neg TDC", "fNegTdcRefDiffTime"},
      {"PosAdcRefDiffTime", "Reference Diff time of Pos ADC", "fPosAdcRefDiffTime"},
      {"NegAdcRefDiffTime", "Reference Diff time of Neg aDC", "fNegAdcRefDiffTime"},
      //{"ngoodhits", "Number of paddle hits (passed tof tolerance and used to determine the focal plane time )",
      //"GetNGoodHits() "},
      {0}};

  return DefineVarsFromList(vars, mode);
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::ProcessHits(TClonesArray *rawhits, Int_t nexthit) {

  /*! \brief Extract scintillator paddle hits from raw data starting at "nexthit"
   * - Called by THcHodoscope::Decode
   * - Loops through "rawhits" array  starting at index of "nexthit"
   * - Assumes that the hit list is sorted by plane and looping ends when plane number of hit doesn't match fPlaneNum
   * - Fills THcSignalHit objects frPosTDCHits and frNegTDCHits when TDC > 0
   * - Fills THcSignalHit objects frPosADCHits and frNegADCHit with pedestal subtracted ADC when value larger than 50
   * - For hits that have TDC value for either positive or negative PMT within  fScinTdcMin and fScinTdcMax
   *  + Creates new  fHodoHits[fNScinHits] =  THcHodoHit
   *  + Calculates pulse height correction to the positive and negative PMT times
   *  + Correct times for time traveled in paddle
   *  + Correct times for time of flight using beta from central spectrometer momentum and particle type
   *  + Calls  SetCorrectedTime method of THcHodoHit
   *  + Increments fNScinHits
   * - Returns value of nexthit + number of hits processed
   *
   */

  // Clear() is being called event by event

  // counters for Tdc/Adc hits
  UInt_t nrPosTdcHits = 0;
  UInt_t nrPosAdcHits = 0;
  UInt_t nrNegTdcHits = 0;
  UInt_t nrNegAdcHits = 0;
  UInt_t nrSampPosAdcHits = 0;
  UInt_t nrSampNegAdcHits = 0;
  UInt_t nrPosADCHits = 0;
  UInt_t nrNegADCHits = 0;

  Int_t nrawhits = rawhits->GetLast() + 1;
  Int_t ihit = nexthit;

  // A THcRawHodoHit contains all the information (tdc and adc for both
  // pmts) for a single paddle for a single trigger.  The tdc information
  // might include multiple hits if it uses a multihit tdc.
  // Use "ihit" as the index over THcRawHodoHit objects.  Use
  // "thit" to index over multiple tdc hits within an "ihit".

  Bool_t problem_flag = kFALSE; // check if fTdcRefTime is filled correctly or left initialized (kBig)

  while (ihit < nrawhits) {

    // I think we use THcRawHodoHit as it is
    THcRawHodoHit *hit = (THcRawHodoHit *)rawhits->At(ihit);

    if (hit->fPlane > fPlaneNum) {
      break;
    }

    Int_t padnum = hit->fCounter;
    Int_t index = padnum - 1;

    // Positive Tdc hits
    THcRawTdcHit &rawPosTdcHit = hit->GetRawTdcHitPos();
    if (rawPosTdcHit.GetNHits() > 0 && rawPosTdcHit.HasRefTime()) {

      // Assume fPosTdcRefTime is initialized
      if (fPosTdcRefTime == kBig) {
        fPosTdcRefTime = rawPosTdcHit.GetRefTime();
        fPosTdcRefDiffTime = rawPosTdcHit.GetRefDiffTime();
      }

      // Set problem_flag if the it's not set correctly
      if (fPosTdcRefTime != rawPosTdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    // Loop over multiple tdc hits within ihit
    for (UInt_t thit = 0; thit < rawPosTdcHit.GetNHits(); thit++) {
      ((THcSignalHit *)frPosTdcTimeRaw->ConstructedAt(nrPosTdcHits))->Set(padnum, rawPosTdcHit.GetTimeRaw(thit));
      ((THcSignalHit *)frPosTdcTime->ConstructedAt(nrPosTdcHits))->Set(padnum, rawPosTdcHit.GetTime(thit));

      nrPosTdcHits++; // FIXME: just use thit? or is it used somewhere else too?
      fTotNumTdcHits++;
      fTotNumPosTdcHits++;
    }

    // Now, repeat for the Neg end
    THcRawTdcHit &rawNegTdcHit = hit->GetRawTdcHitNeg();
    if (rawNegTdcHit.GetNHits() > 0 && rawNegTdcHit.HasRefTime()) {

      if (fNegTdcRefTime == kBig) {
        fNegTdcRefTime = rawNegTdcHit.GetRefTime();
        fNegTdcRefDiffTime = rawNegTdcHit.GetRefDiffTime();
      }

      if (fNegTdcRefTime != rawNegTdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    for (UInt_t thit = 0; thit < rawNegTdcHit.GetNHits(); thit++) {
      ((THcSignalHit *)frNegTdcTimeRaw->ConstructedAt(nrNegTdcHits))->Set(padnum, rawNegTdcHit.GetTimeRaw(thit));
      ((THcSignalHit *)frNegTdcTime->ConstructedAt(nrNegTdcHits))->Set(padnum, rawNegTdcHit.GetTime(thit));

      nrNegTdcHits++; // FIXME: just use thit? or is it used somewhere else too?
      fTotNumTdcHits++;
      fTotNumNegTdcHits++;
    } // thit loop

    // Positive ADC hits
    THcRawAdcHit &rawPosAdcHit = hit->GetRawAdcHitPos();

    if ((rawPosAdcHit.GetNPulses() > 0 || rawPosAdcHit.GetNSamples() > 0) && rawPosAdcHit.HasRefTime()) {

      if (fPosAdcRefTime == kBig) {
        fPosAdcRefTime = rawPosAdcHit.GetRefTime();
        fPosAdcRefDiffTime = rawPosAdcHit.GetRefDiffTime();
      }

      if (fPosAdcRefTime != rawPosAdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    if (fUseSampWaveform == 0) {

      for (UInt_t thit = 0; thit < rawPosAdcHit.GetNPulses(); thit++) {

        ((THcSignalHit *)frPosAdcPedRaw->ConstructedAt(nrPosAdcHits))->Set(padnum, rawPosAdcHit.GetPedRaw());
        ((THcSignalHit *)frPosAdcPed->ConstructedAt(nrPosAdcHits))->Set(padnum, rawPosAdcHit.GetPed());

        ((THcSignalHit *)frPosAdcPulseIntRaw->ConstructedAt(nrPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetPulseIntRaw(thit));
        ((THcSignalHit *)frPosAdcPulseInt->ConstructedAt(nrPosAdcHits))->Set(padnum, rawPosAdcHit.GetPulseInt(thit));

        ((THcSignalHit *)frPosAdcPulseAmpRaw->ConstructedAt(nrPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetPulseAmpRaw(thit));
        ((THcSignalHit *)frPosAdcPulseAmp->ConstructedAt(nrPosAdcHits))->Set(padnum, rawPosAdcHit.GetPulseAmp(thit));

        ((THcSignalHit *)frPosAdcPulseTimeRaw->ConstructedAt(nrPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetPulseTimeRaw(thit));
        ((THcSignalHit *)frPosAdcPulseTime->ConstructedAt(nrPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetPulseTime(thit) + fAdcTdcOffset);

        // Error flags 0-2
        if (rawPosAdcHit.GetPulseAmpRaw(thit) > 0)
          ((THcSignalHit *)frPosAdcErrorFlag->ConstructedAt(nrPosAdcHits))->Set(padnum, 0);
        if (rawPosAdcHit.GetPulseAmpRaw(thit) <= 0)
          ((THcSignalHit *)frPosAdcErrorFlag->ConstructedAt(nrPosAdcHits))->Set(padnum, 1);
        if (rawPosAdcHit.GetPulseAmpRaw(thit) <= 0 && rawPosAdcHit.GetNSamples() > 0)
          ((THcSignalHit *)frPosAdcErrorFlag->ConstructedAt(nrPosAdcHits))->Set(padnum, 2);

        nrPosAdcHits++;
        fTotNumAdcHits++;
        fTotNumPosAdcHits++;
      }
    }

    // Use Waveform -- Do we use it? Leave it just in case
    if (rawPosAdcHit.GetNSamples() > 0) {
      rawPosAdcHit.SetSampThreshold(fSampThreshold);
      if (fSampNSA == 0)
        fSampNSA = rawPosAdcHit.GetF250_NSA();
      if (fSampNSB == 0)
        fSampNSB = rawPosAdcHit.GetF250_NSB();

      rawPosAdcHit.SetF250Params(fSampNSA, fSampNSB, 4); // Set NPED =4
      if (fSampNSAT != 2)
        rawPosAdcHit.SetSampNSAT(fSampNSAT);
      rawPosAdcHit.SetSampIntTimePedestalPeak();
      fPosAdcSampWaveform.push_back(float(padnum));
      fPosAdcSampWaveform.push_back(float(rawPosAdcHit.GetNSamples()));

      for (UInt_t thit = 0; thit < rawPosAdcHit.GetNSamples(); thit++) {
        fPosAdcSampWaveform.push_back(rawPosAdcHit.GetSample(thit)); // ped subtracted sample (mV)
      }
      for (UInt_t thit = 0; thit < rawPosAdcHit.GetNSampPulses(); thit++) {
        ((THcSignalHit *)frPosAdcSampPedRaw->ConstructedAt(nrSampPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetSampPedRaw());
        ((THcSignalHit *)frPosAdcSampPed->ConstructedAt(nrSampPosAdcHits))->Set(padnum, rawPosAdcHit.GetSampPed());

        ((THcSignalHit *)frPosAdcSampPulseIntRaw->ConstructedAt(nrSampPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetSampPulseIntRaw(thit));
        ((THcSignalHit *)frPosAdcSampPulseInt->ConstructedAt(nrSampPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetSampPulseInt(thit));

        ((THcSignalHit *)frPosAdcSampPulseAmpRaw->ConstructedAt(nrSampPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetSampPulseAmpRaw(thit));
        ((THcSignalHit *)frPosAdcSampPulseAmp->ConstructedAt(nrSampPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetSampPulseAmp(thit));

        ((THcSignalHit *)frPosAdcSampPulseTimeRaw->ConstructedAt(nrSampPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetSampPulseTimeRaw(thit));
        ((THcSignalHit *)frPosAdcSampPulseTime->ConstructedAt(nrSampPosAdcHits))
            ->Set(padnum, rawPosAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);

        if (rawPosAdcHit.GetNPulses() == 0 || fUseSampWaveform == 1) {
          ((THcSignalHit *)frPosAdcPedRaw->ConstructedAt(nrPosAdcHits))->Set(padnum, rawPosAdcHit.GetSampPedRaw());
          ((THcSignalHit *)frPosAdcPed->ConstructedAt(nrPosAdcHits))->Set(padnum, rawPosAdcHit.GetSampPed());

          ((THcSignalHit *)frPosAdcPulseIntRaw->ConstructedAt(nrPosAdcHits))
              ->Set(padnum, rawPosAdcHit.GetSampPulseIntRaw(thit));
          ((THcSignalHit *)frPosAdcPulseInt->ConstructedAt(nrPosAdcHits))
              ->Set(padnum, rawPosAdcHit.GetSampPulseInt(thit));

          ((THcSignalHit *)frPosAdcPulseAmpRaw->ConstructedAt(nrPosAdcHits))
              ->Set(padnum, rawPosAdcHit.GetSampPulseAmpRaw(thit));
          ((THcSignalHit *)frPosAdcPulseAmp->ConstructedAt(nrPosAdcHits))
              ->Set(padnum, rawPosAdcHit.GetSampPulseAmp(thit));

          ((THcSignalHit *)frPosAdcPulseTimeRaw->ConstructedAt(nrPosAdcHits))
              ->Set(padnum, rawPosAdcHit.GetSampPulseTimeRaw(thit));
          ((THcSignalHit *)frPosAdcPulseTime->ConstructedAt(nrPosAdcHits))
              ->Set(padnum, rawPosAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);
          ((THcSignalHit *)frPosAdcErrorFlag->ConstructedAt(nrPosAdcHits))->Set(padnum, 3);
          if (fUseSampWaveform == 1)
            ((THcSignalHit *)frPosAdcErrorFlag->ConstructedAt(nrPosAdcHits))->Set(padnum, 0);

          ++nrPosAdcHits;
          fTotNumPosAdcHits++;
          fTotNumAdcHits++;
        }
        ++nrSampPosAdcHits;
      }
    }

    // Neg ADC hits
    THcRawAdcHit &rawNegAdcHit = hit->GetRawAdcHitNeg();
    if ((rawNegAdcHit.GetNPulses() > 0 || rawNegAdcHit.GetNSamples() > 0) && rawNegAdcHit.HasRefTime()) {

      if (fNegAdcRefTime == kBig) {
        fNegAdcRefTime = rawNegAdcHit.GetRefTime();
        fNegAdcRefDiffTime = rawNegAdcHit.GetRefDiffTime();
      }

      if (fNegAdcRefTime != rawNegAdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    if (fUseSampWaveform == 0) {

      for (UInt_t thit = 0; thit < rawNegAdcHit.GetNPulses(); thit++) {

        ((THcSignalHit *)frNegAdcPedRaw->ConstructedAt(nrNegAdcHits))->Set(padnum, rawNegAdcHit.GetPedRaw());
        ((THcSignalHit *)frNegAdcPed->ConstructedAt(nrNegAdcHits))->Set(padnum, rawNegAdcHit.GetPed());

        ((THcSignalHit *)frNegAdcPulseIntRaw->ConstructedAt(nrNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetPulseIntRaw(thit));
        ((THcSignalHit *)frNegAdcPulseInt->ConstructedAt(nrNegAdcHits))->Set(padnum, rawNegAdcHit.GetPulseInt(thit));

        ((THcSignalHit *)frNegAdcPulseAmpRaw->ConstructedAt(nrNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetPulseAmpRaw(thit));
        ((THcSignalHit *)frNegAdcPulseAmp->ConstructedAt(nrNegAdcHits))->Set(padnum, rawNegAdcHit.GetPulseAmp(thit));

        ((THcSignalHit *)frNegAdcPulseTimeRaw->ConstructedAt(nrNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetPulseTimeRaw(thit));
        ((THcSignalHit *)frNegAdcPulseTime->ConstructedAt(nrNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetPulseTime(thit) + fAdcTdcOffset);

        // Error flags 0-2
        if (rawNegAdcHit.GetPulseAmpRaw(thit) > 0)
          ((THcSignalHit *)frNegAdcErrorFlag->ConstructedAt(nrNegAdcHits))->Set(padnum, 0);
        if (rawNegAdcHit.GetPulseAmpRaw(thit) <= 0)
          ((THcSignalHit *)frNegAdcErrorFlag->ConstructedAt(nrNegAdcHits))->Set(padnum, 1);
        if (rawNegAdcHit.GetPulseAmpRaw(thit) <= 0 && rawNegAdcHit.GetNSamples() > 0)
          ((THcSignalHit *)frNegAdcErrorFlag->ConstructedAt(nrNegAdcHits))->Set(padnum, 2);

        nrNegAdcHits++;
        fTotNumAdcHits++;
        fTotNumNegAdcHits++;
      }
    }

    // Use Waveform -- Do we use it? Leave it just in case
    if (rawNegAdcHit.GetNSamples() > 0) {
      rawNegAdcHit.SetSampThreshold(fSampThreshold);
      if (fSampNSA == 0)
        fSampNSA = rawNegAdcHit.GetF250_NSA();
      if (fSampNSB == 0)
        fSampNSB = rawNegAdcHit.GetF250_NSB();

      rawNegAdcHit.SetF250Params(fSampNSA, fSampNSB, 4); // Set NPED =4
      if (fSampNSAT != 2)
        rawNegAdcHit.SetSampNSAT(fSampNSAT);
      rawNegAdcHit.SetSampIntTimePedestalPeak();
      fNegAdcSampWaveform.push_back(float(padnum));
      fNegAdcSampWaveform.push_back(float(rawNegAdcHit.GetNSamples()));

      for (UInt_t thit = 0; thit < rawNegAdcHit.GetNSamples(); thit++) {
        fNegAdcSampWaveform.push_back(rawNegAdcHit.GetSample(thit)); // ped subtracted sample (mV)
      }
      for (UInt_t thit = 0; thit < rawNegAdcHit.GetNSampPulses(); thit++) {
        ((THcSignalHit *)frNegAdcSampPedRaw->ConstructedAt(nrSampNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetSampPedRaw());
        ((THcSignalHit *)frNegAdcSampPed->ConstructedAt(nrSampNegAdcHits))->Set(padnum, rawNegAdcHit.GetSampPed());

        ((THcSignalHit *)frNegAdcSampPulseIntRaw->ConstructedAt(nrSampNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetSampPulseIntRaw(thit));
        ((THcSignalHit *)frNegAdcSampPulseInt->ConstructedAt(nrSampNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetSampPulseInt(thit));

        ((THcSignalHit *)frNegAdcSampPulseAmpRaw->ConstructedAt(nrSampNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetSampPulseAmpRaw(thit));
        ((THcSignalHit *)frNegAdcSampPulseAmp->ConstructedAt(nrSampNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetSampPulseAmp(thit));

        ((THcSignalHit *)frNegAdcSampPulseTimeRaw->ConstructedAt(nrSampNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetSampPulseTimeRaw(thit));
        ((THcSignalHit *)frNegAdcSampPulseTime->ConstructedAt(nrSampNegAdcHits))
            ->Set(padnum, rawNegAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);

        if (rawNegAdcHit.GetNPulses() == 0 || fUseSampWaveform == 1) {
          ((THcSignalHit *)frNegAdcPedRaw->ConstructedAt(nrNegAdcHits))->Set(padnum, rawNegAdcHit.GetSampPedRaw());
          ((THcSignalHit *)frNegAdcPed->ConstructedAt(nrNegAdcHits))->Set(padnum, rawNegAdcHit.GetSampPed());

          ((THcSignalHit *)frNegAdcPulseIntRaw->ConstructedAt(nrNegAdcHits))
              ->Set(padnum, rawNegAdcHit.GetSampPulseIntRaw(thit));
          ((THcSignalHit *)frNegAdcPulseInt->ConstructedAt(nrNegAdcHits))
              ->Set(padnum, rawNegAdcHit.GetSampPulseInt(thit));

          ((THcSignalHit *)frNegAdcPulseAmpRaw->ConstructedAt(nrNegAdcHits))
              ->Set(padnum, rawNegAdcHit.GetSampPulseAmpRaw(thit));
          ((THcSignalHit *)frNegAdcPulseAmp->ConstructedAt(nrNegAdcHits))
              ->Set(padnum, rawNegAdcHit.GetSampPulseAmp(thit));

          ((THcSignalHit *)frNegAdcPulseTimeRaw->ConstructedAt(nrNegAdcHits))
              ->Set(padnum, rawNegAdcHit.GetSampPulseTimeRaw(thit));
          ((THcSignalHit *)frNegAdcPulseTime->ConstructedAt(nrNegAdcHits))
              ->Set(padnum, rawNegAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);
          ((THcSignalHit *)frNegAdcErrorFlag->ConstructedAt(nrNegAdcHits))->Set(padnum, 3);
          if (fUseSampWaveform == 1)
            ((THcSignalHit *)frNegAdcErrorFlag->ConstructedAt(nrNegAdcHits))->Set(padnum, 0);

          ++nrNegAdcHits;
          fTotNumNegAdcHits++;
          fTotNumAdcHits++;
        }
        ++nrSampNegAdcHits;
      }
    }

    // Determine good tdc pos and neg times
    Bool_t btdcraw_pos = kFALSE;
    Bool_t btdcraw_neg = kFALSE;
    Int_t tdc_pos = -999;
    Int_t tdc_neg = -999;
    Double_t good_ielem_postdc = -1;
    Double_t good_ielem_negtdc = -1;

    // Good TDC Pos times
    for (UInt_t thit = 0; thit < hit->GetRawTdcHitPos().GetNHits(); thit++) {

      tdc_pos = hit->GetRawTdcHitPos().GetTime(thit) + fTdcOffset;

      if (tdc_pos >= fScinTdcMin && tdc_pos <= fScinTdcMax) {
        btdcraw_pos = kTRUE;
        good_ielem_postdc = thit;
        break;
      }
    }

    // Good TDC Neg times
    for (UInt_t thit = 0; thit < hit->GetRawTdcHitNeg().GetNHits(); thit++) {

      tdc_neg = hit->GetRawTdcHitNeg().GetTime(thit) + fTdcOffset;

      if (tdc_neg >= fScinTdcMin && tdc_neg <= fScinTdcMax) {
        btdcraw_neg = kTRUE;
        good_ielem_negtdc = thit;
        break;
      }
    }

    // Adc neg
    Bool_t badcraw_neg = kFALSE;
    Double_t adcped_neg = -999;
    Int_t adcmult_neg = 0;
    Int_t adchitused_neg = 0;
    Double_t adcint_neg = -999;
    Double_t adcamp_neg = -kBig;
    Double_t adctime_neg = kBig;
    Double_t adctdcdifftime_neg = kBig;
    Double_t good_ielem_negadc = -1;

    // Adc pos
    Bool_t badcraw_pos = kFALSE;
    Double_t adcped_pos = -999;
    Int_t adcmult_pos = 0;
    Int_t adchitused_pos = 0;
    Double_t adcint_pos = -999;
    Double_t adcamp_pos = -kBig;
    Double_t adctime_pos = kBig;
    Double_t adctdcdifftime_pos = kBig;
    Double_t good_ielem_posadc = -1;

    if (fADCMode == kADCDynamicPedestal) {

      // Loop Here over all hits per event for neg side of plane
      // Try to find the max pulseAmp and corresponding ielem
      // within the given time window for TdcAdcTimeDiff

      Int_t good_ielem_negadc_test2 = -1;

      if (good_ielem_negtdc != -1) {
        Double_t max_adcamp_test = -1000.;
        Double_t max_adctdcdiff_test = 1000.;

        for (UInt_t ielem = 0; ielem < rawNegAdcHit.GetNPulses(); ielem++) {

          Double_t pulseAmp = rawNegAdcHit.GetPulseAmp(ielem);
          Double_t pulseTime = rawNegAdcHit.GetPulseTime(ielem) + fAdcTdcOffset;
          Double_t TdcAdcTimeDiff = tdc_neg * fScinTdcToTime - pulseTime;

          if (rawNegAdcHit.GetPulseAmpRaw(ielem) <= 0)
            pulseAmp = 200.; // do we want to to this? or skip simply this element?

          Bool_t pulseTimeCut =
              (TdcAdcTimeDiff > fHodoNegAdcTimeWindowMin[index]) && (TdcAdcTimeDiff < fHodoNegAdcTimeWindowMax[index]);
          if (pulseTimeCut && pulseAmp > max_adcamp_test) {
            good_ielem_negadc = ielem;
            max_adcamp_test = pulseAmp;
          }
          if (abs(TdcAdcTimeDiff) < max_adctdcdiff_test) {
            good_ielem_negadc_test2 = ielem;
            max_adctdcdiff_test = abs(TdcAdcTimeDiff);
          }
        }
      }

      // good_ielem_negadc: select the pulse that has TdcAdcTimeDiff within the given range and with maximum pulseAmp
      // good_ielem_negadc_test2: selects the pulse with minimum value of TdcAdcTimeDiff
      if (good_ielem_negadc == -1 && good_ielem_negadc_test2 != -1)
        good_ielem_negadc = good_ielem_negadc_test2;
      if (good_ielem_negadc == -1 && good_ielem_negadc_test2 == -1 && rawNegAdcHit.GetNPulses() > 0)
        good_ielem_negadc = 0;

      if (good_ielem_negadc != -1 && good_ielem_negadc < rawNegAdcHit.GetNPulses()) {
        adcped_neg = rawNegAdcHit.GetPed();
        adcmult_neg = rawNegAdcHit.GetNPulses();
        adchitused_neg = good_ielem_negadc + 1;
        adcint_neg = rawNegAdcHit.GetPulseInt(good_ielem_negadc);
        adcamp_neg = rawNegAdcHit.GetPulseAmp(good_ielem_negadc);
        if (rawNegAdcHit.GetPulseAmpRaw(good_ielem_negadc) <= 0)
          adcamp_neg = 200.;
        adctime_neg = rawNegAdcHit.GetPulseTime(good_ielem_negadc) + fAdcTdcOffset;
        badcraw_neg = kTRUE;
        adctdcdifftime_neg = tdc_neg * fScinTdcToTime - adctime_neg;
      }

      // Loop Here over all hits per event for pos side of plane
      Int_t good_ielem_posadc_test2 = -1;

      if (good_ielem_postdc != -1) {
        Double_t max_adcamp_test = -1000.;
        Double_t max_adctdcdiff_test = 1000.;
        //
        for (UInt_t ielem = 0; ielem < rawPosAdcHit.GetNPulses(); ielem++) {
          Double_t pulseAmp = rawPosAdcHit.GetPulseAmp(ielem);
          Double_t pulseTime = rawPosAdcHit.GetPulseTime(ielem) + fAdcTdcOffset;
          Double_t TdcAdcTimeDiff = tdc_pos * fScinTdcToTime - pulseTime;
          Bool_t pulseTimeCut =
              (TdcAdcTimeDiff > fHodoPosAdcTimeWindowMin[index]) && (TdcAdcTimeDiff < fHodoPosAdcTimeWindowMax[index]);
          if (rawPosAdcHit.GetPulseAmpRaw(ielem) <= 0)
            pulseAmp = 200.;
          if (pulseTimeCut && pulseAmp > max_adcamp_test) {
            good_ielem_posadc = ielem;
            max_adcamp_test = pulseAmp;
          }
          if (abs(TdcAdcTimeDiff) < max_adctdcdiff_test) {
            good_ielem_posadc_test2 = ielem;
            max_adctdcdiff_test = abs(TdcAdcTimeDiff);
          }
        }
      }

      if (good_ielem_posadc == -1 && good_ielem_posadc_test2 != -1)
        good_ielem_posadc = good_ielem_posadc_test2;
      if (good_ielem_posadc == -1 && good_ielem_posadc_test2 == -1 && rawPosAdcHit.GetNPulses() > 0)
        good_ielem_posadc = 0;
      if (good_ielem_posadc != -1 && good_ielem_posadc < rawPosAdcHit.GetNPulses()) {
        adcped_pos = rawPosAdcHit.GetPed();
        adcmult_pos = rawPosAdcHit.GetNPulses();
        adchitused_pos = good_ielem_posadc + 1;
        adcint_pos = rawPosAdcHit.GetPulseInt(good_ielem_posadc);
        adcamp_pos = rawPosAdcHit.GetPulseAmp(good_ielem_posadc);
        if (rawPosAdcHit.GetPulseAmpRaw(good_ielem_posadc) <= 0)
          adcamp_pos = 200.;
        adctime_pos = rawPosAdcHit.GetPulseTime(good_ielem_posadc) + fAdcTdcOffset;
        badcraw_pos = kTRUE;
        adctdcdifftime_pos = tdc_pos * fScinTdcToTime - adctime_pos;
      }

    } // if kADCDynamicPedestal

    else if (fADCMode == kADCSampleIntegral) {
      adcint_pos = hit->GetRawAdcHitPos().GetSampleIntRaw() - fPosPed[index];
      adcint_neg = hit->GetRawAdcHitNeg().GetSampleIntRaw() - fNegPed[index];
      badcraw_pos = badcraw_neg = kTRUE;

    } else if (fADCMode == kADCSampIntDynPed) {
      adcint_pos = hit->GetRawAdcHitPos().GetSampleInt();
      adcint_neg = hit->GetRawAdcHitNeg().GetSampleInt();
      badcraw_pos = badcraw_neg = kTRUE;

    } else {
      // same as using kADCSampleIntegral
      adcint_pos = hit->GetRawAdcHitPos().GetPulseIntRaw() - fPosPed[index];
      adcint_neg = hit->GetRawAdcHitNeg().GetPulseIntRaw() - fNegPed[index];
      badcraw_pos = badcraw_neg = kTRUE;
    }

    // We don't really do anything with these arrays (frPosAdcHits, frNegAdcHits, ..)
    if (adcint_pos >= fADCDiagCut) {
      ((THcSignalHit *)frPosADCHits->ConstructedAt(nrPosADCHits))->Set(padnum, adcint_pos);
      Double_t samplesum = hit->GetRawAdcHitPos().GetSampleIntRaw();
      Double_t pedestal = hit->GetRawAdcHitPos().GetPedRaw();
      ((THcSignalHit *)frPosADCSums->ConstructedAt(nrPosADCHits))->Set(padnum, samplesum);
      ((THcSignalHit *)frPosADCPeds->ConstructedAt(nrPosADCHits++))->Set(padnum, pedestal);
    }
    if (adcint_neg >= fADCDiagCut) {
      ((THcSignalHit *)frNegADCHits->ConstructedAt(nrNegADCHits))->Set(padnum, adcint_neg);
      Double_t samplesum = hit->GetRawAdcHitNeg().GetSampleIntRaw();
      Double_t pedestal = hit->GetRawAdcHitNeg().GetPedRaw();
      ((THcSignalHit *)frNegADCSums->ConstructedAt(nrNegADCHits))->Set(padnum, samplesum);
      ((THcSignalHit *)frNegADCPeds->ConstructedAt(nrNegADCHits++))->Set(padnum, pedestal);
    }

    // Save counters
    if ((btdcraw_pos && badcraw_pos) || (btdcraw_neg && badcraw_neg)) {

      if (good_ielem_posadc != -1) {

        // good adc multiplicity
        fTotNumGoodPosAdcHits++;
        fTotNumGoodAdcHits++;

        // good adc occupancy
        fNumGoodPosAdcHits.at(padnum - 1) = padnum;
        fGoodPosAdcPed.at(padnum - 1) = adcped_pos;
        fGoodPosAdcMult.at(padnum - 1) = adcmult_pos;
        fGoodPosAdcHitUsed.at(padnum - 1) = adchitused_pos;
        fGoodPosAdcPulseInt.at(padnum - 1) = adcint_pos;
        fGoodPosAdcPulseAmp.at(padnum - 1) = adcamp_pos;
        fGoodPosAdcPulseTime.at(padnum - 1) = adctime_pos;
        fGoodPosAdcTdcDiffTime.at(padnum - 1) = adctdcdifftime_pos;
      }

      if (good_ielem_negadc != -1) {

        // good adc multiplicity
        fTotNumGoodNegAdcHits++;
        fTotNumGoodAdcHits++;

        // good adc occupancy
        fNumGoodNegAdcHits.at(padnum - 1) = padnum;
        fGoodNegAdcPed.at(padnum - 1) = adcped_neg;
        fGoodNegAdcMult.at(padnum - 1) = adcmult_neg;
        fGoodNegAdcHitUsed.at(padnum - 1) = adchitused_neg;
        fGoodNegAdcPulseInt.at(padnum - 1) = adcint_neg;
        fGoodNegAdcPulseAmp.at(padnum - 1) = adcamp_neg;
        fGoodNegAdcPulseTime.at(padnum - 1) = adctime_neg;
        fGoodNegAdcTdcDiffTime.at(padnum - 1) = adctdcdifftime_neg;
      }

      // DEFINE THE "GOOD +TDC Multiplicities and Occupancies"
      if (good_ielem_postdc != -1) {
        fTotNumGoodPosTdcHits++;
        fTotNumGoodTdcHits++;
        // good tdc occupancy
        fNumGoodPosTdcHits.at(padnum - 1) = padnum;
      }

      // DEFINE THE "GOOD -TDC  Multiplicities and Occupancies"
      if (good_ielem_negtdc != -1) {
        fTotNumGoodNegTdcHits++;
        fTotNumGoodTdcHits++;
        // good tdc occupancy
        fNumGoodNegTdcHits.at(padnum - 1) = padnum;
      }

      new ((*fHodoHits)[fNScinHits]) THcLADHodoHit(tdc_pos, tdc_neg, adcint_pos, adcint_neg, hit->fCounter, this);

      ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetPosADCpeak(adcamp_pos);
      ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetNegADCpeak(adcamp_neg);
      ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetPosADCtime(adctime_pos);
      ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetNegADCtime(adctime_neg);

      fNScinHits++;
    } // if valid tdc/adc at least in one end

    ihit++;
  } // while loop

  if (problem_flag) {
    cout << "THcScintillatorPlane::ProcessHits " << fPlaneNum << " " << nexthit << "/" << nrawhits << endl;
    cout << " Ref problem end *******" << endl;
  }

  return (ihit);
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::CoarseProcess(TClonesArray &tracks) {

  // Probably we won't do anything here

  return 0;
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::FineProcess(TClonesArray &tracks) {

  // Probably we won't do anything here

  return 0;
}

//_____________________________________________________________________________
Int_t THcLADHodoPlane::AccumulatePedestals(TClonesArray *rawhits, Int_t nexthit) {
  /*! \brief Extract the data for this plane from raw hit list THcRawHodoHit, accumulating into arrays for calculating
   * pedestals.
   *
   * - Loop through raw data for scintillator plane
   */
  Int_t nrawhits = rawhits->GetLast() + 1;
  // cout << "THcScintillatorPlane::AcculatePedestals " << fPlaneNum << " " << nexthit << "/" << nrawhits << endl;

  Int_t ihit = nexthit;
  while (ihit < nrawhits) {
    THcRawHodoHit *hit = (THcRawHodoHit *)rawhits->At(ihit);
    if (hit->fPlane > fPlaneNum) {
      break;
    }
    Int_t element = hit->fCounter - 1; // Should check if in range
    Int_t adcpos = hit->GetRawAdcHitPos().GetPulseIntRaw();
    Int_t adcneg = hit->GetRawAdcHitNeg().GetPulseIntRaw();

    if (adcpos <= fPosPedLimit[element]) {
      fPosPedSum[element] += adcpos;
      fPosPedSum2[element] += adcpos * adcpos;
      fPosPedCount[element]++;
      if (fPosPedCount[element] == fMinPeds / 5) {
        fPosPedLimit[element] = 100 + fPosPedSum[element] / fPosPedCount[element];
      }
    }
    if (adcneg <= fNegPedLimit[element]) {
      fNegPedSum[element] += adcneg;
      fNegPedSum2[element] += adcneg * adcneg;
      fNegPedCount[element]++;
      if (fNegPedCount[element] == fMinPeds / 5) {
        fNegPedLimit[element] = 100 + fNegPedSum[element] / fNegPedCount[element];
      }
    }
    ihit++;
  }

  fNPedestalEvents++;

  return (ihit);
}

//_______________________________________________________________________________________
void THcLADHodoPlane::CalculatePedestals() {
  /*! \brief   Calculate pedestals from arrays made in THcScintillatorPlane::AccumulatePedestals
   *
   * - Calculate pedestals from arrays made in THcScintillatorPlane::AccumulatePedestals
   * - In old fortran ENGINE code, a comparison was made between calculated pedestals and the pedestals read in by the
   * FASTBUS modules for zero supression. This is not implemented.
   */
  for (Int_t i = 0; i < fNelem; i++) {

    // Positive tubes
    fPosPed[i] = ((Double_t)fPosPedSum[i]) / TMath::Max(1, fPosPedCount[i]);
    fPosThresh[i] = fPosPed[i] + 15;

    // Negative tubes
    fNegPed[i] = ((Double_t)fNegPedSum[i]) / TMath::Max(1, fNegPedCount[i]);
    fNegThresh[i] = fNegPed[i] + 15;

    //    cout <<"Pedestals "<< i+1 << " " << fPosPed[i] << " " << fNegPed[i] << endl;
  }
  //  cout << " " << endl;
}

//_____________________________________________________________________________
void THcLADHodoPlane::InitializePedestals() {
  /*! \brief   called by THcHodoPlane::ReadDatabase
   *
   * - Initialize variables used in  THcScintillatorPlane::AccumulatePedestals and
   * THcScintillatorPlane::CalculatePedestals
   * - Minimum number of pedestal events needed for calculation, fMinPeds, hadrcoded to 500
   */
  fNPedestalEvents = 0;
  fMinPeds = 500; // In engine, this is set in parameter file
  fPosPedSum = new Int_t[fNelem];
  fPosPedSum2 = new Int_t[fNelem];
  fPosPedLimit = new Int_t[fNelem];
  fPosPedCount = new Int_t[fNelem];
  fNegPedSum = new Int_t[fNelem];
  fNegPedSum2 = new Int_t[fNelem];
  fNegPedLimit = new Int_t[fNelem];
  fNegPedCount = new Int_t[fNelem];

  fPosPed = new Double_t[fNelem];
  fNegPed = new Double_t[fNelem];
  fPosThresh = new Double_t[fNelem];
  fNegThresh = new Double_t[fNelem];
  for (Int_t i = 0; i < fNelem; i++) {
    fPosPedSum[i] = 0;
    fPosPedSum2[i] = 0;
    fPosPedLimit[i] = 1000; // In engine, this are set in parameter file
    fPosPedCount[i] = 0;
    fNegPedSum[i] = 0;
    fNegPedSum2[i] = 0;
    fNegPedLimit[i] = 1000; // In engine, this are set in parameter file
    fNegPedCount[i] = 0;
  }
}
//____________________________________________________________________________

ClassImp(THcLADHodoPlane)
